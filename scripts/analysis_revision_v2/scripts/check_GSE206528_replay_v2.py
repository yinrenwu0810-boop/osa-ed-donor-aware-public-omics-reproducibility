#!/usr/bin/env python3
"""Independent, fail-closed acceptance checker for GSE206528 Gate 03A."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import traceback
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
ATTEMPT_ROOT = (REVISION / "01_work/GSE206528_replay").resolve()
VALIDATION = REVISION / "04_validation"
RUNNER = REVISION / "scripts/run_GSE206528_repro_v2.py"
MANIFEST = REVISION / "00_protocol/gate_03_provenance_manifest.sha256.tsv"
MATRIX_RECORD = VALIDATION / "gate_03_matrix_to_raw_validation_attempt_02_PASS.json"
RAW = ROOT / "analysis/data/processed/GSE206528_raw.h5ad"
HISTORICAL = ROOT / "analysis/data/processed/GSE206528_clustered.h5ad"
HISTORICAL_SCORED = ROOT / "analysis/data/processed/GSE206528_scored.h5ad"
HISTORICAL_ANNOTATED = ROOT / "analysis/data/processed/GSE206528_annotated.h5ad"
LEGACY_COUNTS = ROOT / "analysis/results/GSE206528/pseudobulk_counts.csv"
LEGACY_SAMPLES = ROOT / "analysis/results/GSE206528/pseudobulk_samples.csv"
SAMPLE_DESIGN = ROOT / "analysis/data/metadata/core_sample_design.tsv"
SIGNATURES = ROOT / "analysis/config/GSE206528_celltype_signatures.json"
MAPPING = ROOT / "analysis/config/GSE206528_cluster_celltype_mapping.json"
CONFIG = ROOT / "analysis/config/GSE206528_scanpy.json"

EXPECTED_MANIFEST_SHA256 = "c832788e13681be8c2d0d73e3d46cbbaa8e7b1eb4f579a1faaff70b90671043b"
EXPECTED_MATRIX_SHA256 = "7883dee813a7549aa5ea2c4c94c9cb47e40be146b17269c291bd67832de65795"
EXPECTED_VERSIONS = {
    "scanpy": "1.12.1",
    "anndata": "0.13.1",
    "harmonypy": "0.0.10",
    "igraph": "1.0.0",
    "leidenalg": "0.12.0",
    "numpy": "2.4.6",
    "pandas": "3.0.3",
    "scipy": "1.18.0",
    "scikit-learn": "1.9.0",
}
EXPECTED_CONFIG = {
    "min_genes": 200,
    "max_genes": 8000,
    "min_cells": 3,
    "mt_threshold": 20,
    "scrublet": False,
    "target_sum": 10000,
    "n_top_genes": 3000,
    "hvg_flavor": "seurat",
    "scale": False,
    "regress_out": None,
    "n_pcs": 40,
    "n_neighbors": 15,
    "resolution": 0.6,
    "batch_key": "sample_accession",
    "batch_method": "harmony",
    "marker_method": "t-test",
    "skip_markers": False,
}
EXPECTED_HARMONY = {
    "random_state": 0,
    "sigma": 0.1,
    "tau": 0,
    "block_size": 0.05,
    "max_iter_harmony": 10,
    "max_iter_kmeans": 20,
    "epsilon_cluster": 1e-5,
    "epsilon_harmony": 1e-4,
}
EXPECTED_LEIDEN = {
    "flavor": "igraph",
    "n_iterations": 2,
    "directed": False,
    "random_state": 0,
}
EXPECTED_MATRIX_CHECKS = {
    "archive_has_eight_matching_members",
    "raw_shape_exact",
    "eight_accessions_unique",
    "cell_order_exact_all_samples",
    "gene_order_exact_all_samples",
    "matrix_counts_exact",
    "raw_X_and_counts_layer_exact",
    "counts_nonnegative_integers",
    "total_umi_conserved",
}
EXPECTED_MATRIX_TOTALS = {
    "cells": 64993,
    "genes_per_matrix": 33694,
    "source_vs_raw_mismatches": 0,
    "raw_X_vs_counts_layer_mismatches": 0,
    "max_absolute_difference": 0.0,
    "source_total_umi": 481449114,
    "raw_total_umi": 481449114,
    "counts_layer_total_umi": 481449114,
}
ARTIFACT_PATHS = {
    "clustered_h5ad": "GSE206528_clustered_v2.h5ad",
    "scored_h5ad": "GSE206528_scored_v2.h5ad",
    "annotated_h5ad": "GSE206528_annotated_v2.h5ad",
    "pseudobulk_counts": "pseudobulk_counts.csv",
    "pseudobulk_samples": "pseudobulk_samples.csv",
    "cluster_alignment": "cluster_alignment_v2.json",
    "signature_coverage": "signature_gene_coverage_v2.tsv",
    "donor_crosswalk": "legacy_to_canonical_donor_crosswalk.tsv",
    "per_donor_qc": "per_donor_QC_v2.tsv",
    "celltype_counts": "celltype_counts_by_donor_v2.tsv",
    "markers": "markers/markers_all_v2.csv",
}
REQUIRED_OUTPUTS = set(ARTIFACT_PATHS.values()) | {"replay_provenance_v2.json"}


def sha256(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def json_file(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_attempt_path(path: Path) -> Path:
    attempt = path.resolve()
    if attempt.parent != ATTEMPT_ROOT:
        raise ValueError(f"Attempt must be a direct child of {ATTEMPT_ROOT}")
    if re.fullmatch(r"attempt_[0-9]{2}", attempt.name) is None:
        raise ValueError(f"Invalid replay attempt name: {attempt.name}")
    if not attempt.is_dir():
        raise ValueError(f"Replay attempt directory is absent: {attempt}")
    return attempt


def write_numbered_attempt(prefix: str, status: str, payload: str) -> Path:
    VALIDATION.mkdir(parents=True, exist_ok=True)
    number = 1
    while True:
        if any(VALIDATION.glob(f"{prefix}_{number:02d}_*.json")):
            number += 1
            continue
        target = VALIDATION / f"{prefix}_{number:02d}_{status}.json"
        try:
            with target.open("x", encoding="utf-8", newline="\n") as handle:
                handle.write(payload)
            return target
        except FileExistsError:
            number += 1


def atomic_replace_json(path: Path, payload: str) -> None:
    counter = 0
    while True:
        temporary = path.with_name(f".{path.name}.{os.getpid()}.{counter}.tmp")
        try:
            handle = temporary.open("x", encoding="utf-8", newline="\n")
            break
        except FileExistsError:
            counter += 1
    try:
        with handle:
            handle.write(payload)
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def sparse_exact(left, right, chunk_size: int = 1000) -> tuple[int, float]:
    if left.shape != right.shape:
        return -1, float("inf")
    mismatches = 0
    max_difference = 0.0
    for start in range(0, left.shape[0], chunk_size):
        stop = min(start + chunk_size, left.shape[0])
        a = left[start:stop]
        b = right[start:stop]
        if not sparse.issparse(a):
            a = sparse.csr_matrix(a)
        if not sparse.issparse(b):
            b = sparse.csr_matrix(b)
        difference = a.astype(np.float64) - b.astype(np.float64)
        mismatches += int(difference.count_nonzero())
        if difference.nnz:
            max_difference = max(max_difference, float(np.max(np.abs(difference.data))))
    return mismatches, max_difference


def dense_columns_exact(left: pd.DataFrame, right: pd.DataFrame, columns: list[str]) -> bool:
    if not left.index.equals(right.index):
        return False
    for column in columns:
        if column not in left or column not in right:
            return False
        a = left[column].to_numpy(dtype=np.float64)
        b = right[column].to_numpy(dtype=np.float64)
        if not np.array_equal(a, b, equal_nan=True):
            return False
    return True


def array_is_finite(value) -> bool:
    array = np.asarray(value)
    return bool(np.isfinite(array).all())


def sparse_is_finite(value) -> bool:
    if sparse.issparse(value):
        return bool(np.isfinite(value.data).all())
    return array_is_finite(value)


def canonical_design() -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    design = pd.read_csv(SAMPLE_DESIGN, sep="\t", dtype=str)
    design = design.loc[design["series_accession"].eq("GSE206528")].copy()
    if len(design) != 8 or not design["sample_accession"].is_unique:
        raise RuntimeError("Canonical design must contain eight unique GSE206528 accessions")
    snapshot = {
        row["sample_accession"]: {"donor": row["subject"], "condition": row["disease"]}
        for _, row in design.iterrows()
    }
    return design, snapshot


def canonical_legacy_pseudobulk() -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = pd.read_csv(LEGACY_SAMPLES, index_col=0, dtype=str)
    counts = pd.read_csv(LEGACY_COUNTS, index_col=0)
    design, _ = canonical_design()
    donor = design.set_index("sample_accession")["subject"].to_dict()
    condition = design.set_index("sample_accession")["disease"].to_dict()
    samples["donor"] = samples["sample_accession"].map(donor)
    samples["condition"] = samples["sample_accession"].map(condition)
    key_columns = ["sample_accession", "donor", "condition", "cell_type"]
    samples["canonical_id"] = samples[key_columns].agg("|".join, axis=1)
    if samples["canonical_id"].duplicated().any():
        raise RuntimeError("Legacy pseudobulk canonical keys are not unique")
    counts.columns = samples["canonical_id"].tolist()
    samples = samples.set_index("canonical_id").sort_index()
    counts = counts.loc[:, samples.index]
    return counts, samples[key_columns]


def validate_provenance(attempt: Path, provenance: dict[str, object]) -> tuple[dict[str, bool], dict[str, object]]:
    checks: dict[str, bool] = {}
    metrics: dict[str, object] = {}
    _, design_snapshot = canonical_design()
    matrix_report = json_file(MATRIX_RECORD)
    current_runner_hash = sha256(RUNNER)
    current_matrix_hash = sha256(MATRIX_RECORD)
    current_manifest_hash = sha256(MANIFEST)

    checks["provenance_status_exact"] = (
        provenance.get("status") == "ARTIFACTS_READY_FOR_ACCEPTANCE_CHECK"
    )
    checks["provenance_attempt_exact"] = provenance.get("attempt") == attempt.name
    checks["provenance_attempt_root_exact"] = (
        Path(str(provenance.get("attempt_root", ""))).resolve() == ATTEMPT_ROOT
    )
    runner = provenance.get("runner", {})
    checks["runner_identity_exact"] = (
        isinstance(runner, dict)
        and Path(str(runner.get("path", ""))).resolve() == RUNNER.resolve()
        and runner.get("sha256") == current_runner_hash
    )
    frozen = provenance.get("frozen_inputs", {})
    checks["manifest_contract_exact"] = (
        isinstance(frozen, dict)
        and frozen.get("verified") is True
        and Path(str(frozen.get("manifest", ""))).resolve() == MANIFEST.resolve()
        and frozen.get("manifest_sha256") == EXPECTED_MANIFEST_SHA256
        and frozen.get("entries") == 31
        and frozen.get("workspace_entries") == 27
        and frozen.get("absolute_entries") == 4
        and current_manifest_hash == EXPECTED_MANIFEST_SHA256
    )
    lineage = provenance.get("matrix_to_raw_semantic_lineage", {})
    matrix_checks = matrix_report.get("checks", {}) if isinstance(matrix_report, dict) else {}
    matrix_totals = matrix_report.get("totals", {}) if isinstance(matrix_report, dict) else {}
    checks["matrix_record_contract_exact"] = (
        isinstance(lineage, dict)
        and Path(str(lineage.get("record", ""))).resolve() == MATRIX_RECORD.resolve()
        and lineage.get("record_sha256") == EXPECTED_MATRIX_SHA256
        and lineage.get("status") == "PASS"
        and set(lineage.get("checks", {})) == EXPECTED_MATRIX_CHECKS
        and all(lineage.get("checks", {}).values())
        and lineage.get("totals") == EXPECTED_MATRIX_TOTALS
        and current_matrix_hash == EXPECTED_MATRIX_SHA256
        and matrix_report.get("status") == "PASS"
        and set(matrix_checks) == EXPECTED_MATRIX_CHECKS
        and all(matrix_checks.values())
        and matrix_totals == EXPECTED_MATRIX_TOTALS
        and len(matrix_report.get("archive_members", [])) == 8
        and len(matrix_report.get("per_sample", [])) == 8
    )
    input_record = provenance.get("input", {})
    checks["raw_input_identity_exact"] = (
        isinstance(input_record, dict)
        and Path(str(input_record.get("path", ""))).resolve() == RAW.resolve()
        and input_record.get("sha256") == sha256(RAW)
    )
    sample_record = provenance.get("canonical_sample_design", {})
    checks["sample_design_identity_exact"] = (
        isinstance(sample_record, dict)
        and Path(str(sample_record.get("path", ""))).resolve() == SAMPLE_DESIGN.resolve()
        and sample_record.get("sha256") == sha256(SAMPLE_DESIGN)
        and sample_record.get("rows") == design_snapshot
    )
    checks["config_exact"] = (
        provenance.get("config") == EXPECTED_CONFIG
        and json_file(CONFIG) == EXPECTED_CONFIG
    )
    checks["harmony_parameters_exact"] = provenance.get("harmony_parameters") == EXPECTED_HARMONY
    checks["leiden_parameters_exact"] = provenance.get("leiden_parameters") == EXPECTED_LEIDEN
    checks["versions_exact"] = provenance.get("versions") == EXPECTED_VERSIONS
    protocol_minimums = provenance.get("protocol_resource_minimums_gb", {})
    requested = provenance.get("requested_resource_thresholds_gb", {})
    resources = provenance.get("resources_at_start", {})
    checks["resource_gate_exact"] = (
        isinstance(protocol_minimums, dict)
        and isinstance(requested, dict)
        and isinstance(resources, dict)
        and protocol_minimums == {"available_memory": 10.0, "available_disk": 15.0}
        and float(requested.get("available_memory", -1)) >= 10.0
        and float(requested.get("available_disk", -1)) >= 15.0
        and float(resources.get("available_memory_gb", -1)) >= 10.0
        and float(resources.get("available_disk_gb", -1)) >= 15.0
    )
    provenance_artifacts = provenance.get("output_artifacts", {})
    checks["artifact_key_set_exact"] = (
        isinstance(provenance_artifacts, dict)
        and set(provenance_artifacts) == set(ARTIFACT_PATHS)
    )
    artifact_details: dict[str, object] = {}
    artifact_checks = []
    if checks["artifact_key_set_exact"]:
        for key, relative_path in ARTIFACT_PATHS.items():
            record = provenance_artifacts[key]
            target = (attempt / relative_path).resolve()
            inside = False
            try:
                target.relative_to(attempt)
                inside = True
            except ValueError:
                pass
            current_bytes = target.stat().st_size if target.is_file() else None
            current_hash = sha256(target) if target.is_file() else None
            valid = (
                isinstance(record, dict)
                and record.get("relative_path") == relative_path
                and inside
                and target.is_file()
                and record.get("bytes") == current_bytes
                and record.get("sha256") == current_hash
            )
            artifact_checks.append(valid)
            artifact_details[key] = {
                "relative_path": relative_path,
                "identity_exact": valid,
                "bytes": current_bytes,
                "sha256": current_hash,
            }
    checks["all_artifact_hashes_exact"] = (
        len(artifact_checks) == len(ARTIFACT_PATHS) and all(artifact_checks)
    )
    alignment = json_file(attempt / ARTIFACT_PATHS["cluster_alignment"])
    checks["alignment_provenance_exact"] = provenance.get("cluster_alignment") == alignment
    metrics["current_runner_sha256"] = current_runner_hash
    metrics["artifact_details"] = artifact_details
    return checks, metrics


def validate_tables(attempt: Path, annotated) -> tuple[dict[str, bool], dict[str, object]]:
    checks: dict[str, bool] = {}
    metrics: dict[str, object] = {}
    design, _ = canonical_design()
    donor_map = design.set_index("sample_accession")["subject"].to_dict()
    condition_map = design.set_index("sample_accession")["disease"].to_dict()
    accessions = set(design["sample_accession"])
    signatures = json_file(SIGNATURES)

    coverage = pd.read_csv(attempt / ARTIFACT_PATHS["signature_coverage"], sep="\t")
    coverage_columns = {
        "signature", "declared", "present", "missing", "coverage",
        "present_genes", "missing_genes",
    }
    coverage_rows_ok = len(coverage) == len(signatures) == 11
    if coverage_rows_ok and set(coverage.columns) == coverage_columns:
        coverage = coverage.set_index("signature")
        coverage_rows_ok = set(coverage.index) == set(signatures)
        for name, genes in signatures.items():
            if name not in coverage.index:
                coverage_rows_ok = False
                continue
            row = coverage.loc[name]
            present_genes = str(row["present_genes"]).split(";") if str(row["present_genes"]) else []
            coverage_rows_ok &= (
                int(row["declared"]) == len(genes)
                and int(row["present"]) == len(genes)
                and int(row["missing"]) == 0
                and float(row["coverage"]) == 1.0
                and present_genes == genes
                and pd.isna(row["missing_genes"])
            )
    else:
        coverage_rows_ok = False
    checks["signature_coverage_exact"] = bool(coverage_rows_ok)

    crosswalk = pd.read_csv(attempt / ARTIFACT_PATHS["donor_crosswalk"], sep="\t", dtype=str)
    crosswalk_columns = {
        "sample_accession", "legacy_donor", "legacy_condition",
        "canonical_donor", "canonical_condition",
    }
    crosswalk_ok = (
        set(crosswalk.columns) == crosswalk_columns
        and len(crosswalk) == 8
        and crosswalk["sample_accession"].is_unique
        and set(crosswalk["sample_accession"]) == accessions
        and crosswalk["canonical_donor"].eq(crosswalk["sample_accession"].map(donor_map)).all()
        and crosswalk["canonical_condition"].eq(
            crosswalk["sample_accession"].map(condition_map)
        ).all()
    )
    checks["donor_crosswalk_exact"] = bool(crosswalk_ok)

    qc = pd.read_csv(attempt / ARTIFACT_PATHS["per_donor_qc"], sep="\t")
    qc_required = {
        "sample_accession", "donor", "condition", "cells", "n_genes_min",
        "n_genes_median", "n_genes_max", "total_counts_median",
        "pct_mt_median", "pct_mt_max",
    }
    qc_ok = (
        set(qc.columns) == qc_required
        and len(qc) == 8
        and qc["sample_accession"].is_unique
        and set(qc["sample_accession"]) == accessions
        and int(qc["cells"].sum()) == 64993
    )
    checks["per_donor_qc_structure_exact"] = bool(qc_ok)

    cell_counts = pd.read_csv(attempt / ARTIFACT_PATHS["celltype_counts"], sep="\t")
    count_required = {
        "sample_accession", "donor", "condition", "cell_type",
        "cell_count", "sample_total", "cell_fraction",
    }
    count_ok = (
        set(cell_counts.columns) == count_required
        and len(cell_counts) == 72
        and not cell_counts.duplicated(
            ["sample_accession", "donor", "condition", "cell_type"]
        ).any()
        and set(cell_counts["sample_accession"]) == accessions
        and int(cell_counts["cell_count"].sum()) == 64993
        and np.allclose(
            cell_counts["cell_fraction"],
            cell_counts["cell_count"] / cell_counts["sample_total"],
            rtol=0,
            atol=1e-15,
        )
    )
    checks["celltype_counts_structure_exact"] = bool(count_ok)

    markers = pd.read_csv(attempt / ARTIFACT_PATHS["markers"])
    marker_required = {
        "historical_aligned_cluster", "names", "scores",
        "logfoldchanges", "pvals", "pvals_adj",
    }
    marker_ok = (
        marker_required.issubset(markers.columns)
        and len(markers) == 450
        and markers["historical_aligned_cluster"].astype(str).nunique() == 18
        and markers["names"].notna().all()
    )
    checks["marker_table_structure_exact"] = bool(marker_ok)
    metrics.update({
        "signature_groups": len(coverage),
        "crosswalk_rows": len(crosswalk),
        "qc_rows": len(qc),
        "celltype_count_rows": len(cell_counts),
        "marker_rows": len(markers),
    })
    return checks, metrics


def validate_complete_attempt(attempt: Path) -> tuple[dict[str, bool], dict[str, object]]:
    provenance = json_file(attempt / "replay_provenance_v2.json")
    if not isinstance(provenance, dict):
        raise RuntimeError("Replay provenance must be a JSON object")
    checks, metrics = validate_provenance(attempt, provenance)
    alignment = json_file(attempt / ARTIFACT_PATHS["cluster_alignment"])
    cell_map = {str(key): value for key, value in json_file(MAPPING).items()}
    signature_names = list(json_file(SIGNATURES))
    score_columns = [f"{name}_score" for name in signature_names]

    objects = [
        ad.read_h5ad(attempt / ARTIFACT_PATHS["clustered_h5ad"], backed="r"),
        ad.read_h5ad(attempt / ARTIFACT_PATHS["scored_h5ad"], backed="r"),
        ad.read_h5ad(attempt / ARTIFACT_PATHS["annotated_h5ad"], backed="r"),
        ad.read_h5ad(HISTORICAL, backed="r"),
        ad.read_h5ad(HISTORICAL_SCORED, backed="r"),
        ad.read_h5ad(HISTORICAL_ANNOTATED, backed="r"),
    ]
    replay, scored, annotated, history, history_scored, history_annotated = objects
    try:
        checks["all_h5ad_shapes_exact"] = all(
            list(value.shape) == [64993, 24505] for value in objects
        )
        checks["all_h5ad_cell_gene_order_exact"] = all(
            value.obs_names.equals(replay.obs_names)
            and value.var_names.equals(replay.var_names)
            for value in objects[1:]
        )
        counts_history, counts_history_max = sparse_exact(
            replay.layers["counts"], history.layers["counts"]
        )
        counts_scored, counts_scored_max = sparse_exact(
            replay.layers["counts"], scored.layers["counts"]
        )
        counts_annotated, counts_annotated_max = sparse_exact(
            replay.layers["counts"], annotated.layers["counts"]
        )
        normalized_history, normalized_history_max = sparse_exact(replay.X, history.X)
        normalized_scored, normalized_scored_max = sparse_exact(replay.X, scored.X)
        normalized_annotated, normalized_annotated_max = sparse_exact(replay.X, annotated.X)
        checks["counts_chain_exact"] = (
            counts_history == counts_scored == counts_annotated == 0
            and counts_history_max == counts_scored_max == counts_annotated_max == 0
        )
        checks["normalized_X_chain_exact"] = (
            normalized_history == normalized_scored == normalized_annotated == 0
            and normalized_history_max == normalized_scored_max == normalized_annotated_max == 0
        )
        history_hvg = set(history.var_names[np.asarray(history.var["highly_variable"]).astype(bool)])
        replay_hvg = set(replay.var_names[np.asarray(replay.var["highly_variable"]).astype(bool)])
        hvg_intersection = len(history_hvg & replay_hvg)
        hvg_union = len(history_hvg | replay_hvg)
        checks["HVG_exact"] = len(replay_hvg) == 3000 and replay_hvg == history_hvg
        checks["embeddings_shape_and_finite"] = (
            list(replay.obsm["X_pca"].shape) == [64993, 50]
            and list(replay.obsm["X_pca_harmony"].shape) == [64993, 50]
            and list(replay.obsm["X_umap"].shape) == [64993, 2]
            and array_is_finite(replay.obsm["X_pca"])
            and array_is_finite(replay.obsm["X_pca_harmony"])
            and array_is_finite(replay.obsm["X_umap"])
        )
        expected_neighbors = {
            "method": "umap",
            "metric": "euclidean",
            "n_neighbors": 15,
            "n_pcs": 40,
            "random_state": 0,
            "use_rep": "X_pca_harmony",
        }
        neighbor_params = dict(replay.uns["neighbors"]["params"])
        checks["neighbor_graph_contract_exact"] = (
            neighbor_params == expected_neighbors
            and replay.uns["neighbors"]["connectivities_key"] == "connectivities"
            and replay.uns["neighbors"]["distances_key"] == "distances"
            and list(replay.obsp["connectivities"].shape) == [64993, 64993]
            and list(replay.obsp["distances"].shape) == [64993, 64993]
            and sparse_is_finite(replay.obsp["connectivities"])
            and sparse_is_finite(replay.obsp["distances"])
        )
        history_labels = history.obs["leiden"].astype(str)
        replay_raw = replay.obs["leiden_replay"].astype(str)
        replay_aligned = replay.obs["leiden"].astype(str)
        ari = float(adjusted_rand_score(history_labels, replay_raw))
        nmi = float(normalized_mutual_info_score(history_labels, replay_raw))
        aligned_rate = float(replay_aligned.eq(history_labels).mean())
        checks["cluster_acceptance_pass"] = (
            replay_raw.nunique() == 18
            and ari >= 0.95
            and aligned_rate >= 0.95
            and alignment.get("replayed_clusters") == 18
            and alignment.get("historical_clusters") == 18
            and np.isclose(float(alignment.get("ari", -1)), ari)
            and np.isclose(float(alignment.get("nmi", -1)), nmi)
            and np.isclose(
                float(alignment.get("hungarian_cell_match_rate", -1)), aligned_rate
            )
            and not alignment.get("unmapped_replay_clusters")
        )
        checks["aligned_leiden_chain_exact"] = (
            scored.obs["leiden"].astype(str).equals(replay_aligned)
            and annotated.obs["leiden"].astype(str).equals(replay_aligned)
        )
        checks["score_column_set_exact"] = (
            sorted(column for column in scored.obs if column.endswith("_score"))
            == sorted(score_columns)
            and sorted(column for column in annotated.obs if column.endswith("_score"))
            == sorted(score_columns)
        )
        checks["score_values_scored_to_annotated_exact"] = dense_columns_exact(
            scored.obs, annotated.obs, score_columns
        )
        checks["score_values_historical_exact"] = dense_columns_exact(
            scored.obs, history_scored.obs, score_columns
        )
        expected_cell_type = annotated.obs["leiden"].astype(str).map(cell_map)
        cell_type_match = float(
            annotated.obs["cell_type"].astype(str).eq(
                history_annotated.obs["cell_type"].astype(str)
            ).mean()
        )
        checks["cell_type_map_exact"] = (
            annotated.obs["cell_type"].astype(str).eq(expected_cell_type).all()
            and cell_type_match >= 0.95
        )
        table_checks, table_metrics = validate_tables(attempt, annotated)
        checks.update(table_checks)
        metrics.update(table_metrics)

        new_samples = pd.read_csv(
            attempt / ARTIFACT_PATHS["pseudobulk_samples"], index_col=0, dtype=str
        ).sort_index()
        new_counts = pd.read_csv(
            attempt / ARTIFACT_PATHS["pseudobulk_counts"], index_col=0
        ).loc[:, new_samples.index]
        legacy_counts, legacy_samples = canonical_legacy_pseudobulk()
        key_columns = ["sample_accession", "donor", "condition", "cell_type"]
        sample_metadata_exact = (
            list(new_samples.columns) == key_columns
            and new_samples.index.tolist() == legacy_samples.index.tolist()
            and new_samples[key_columns].equals(legacy_samples[key_columns])
        )
        genes_exact = (
            new_counts.index.astype(str).tolist()
            == legacy_counts.index.astype(str).tolist()
        )
        if sample_metadata_exact and genes_exact:
            difference = (
                new_counts.to_numpy(dtype=np.float64)
                - legacy_counts.to_numpy(dtype=np.float64)
            )
            pb_mismatches = int(np.count_nonzero(difference))
            pb_max_difference = (
                float(np.max(np.abs(difference))) if pb_mismatches else 0.0
            )
        else:
            pb_mismatches = -1
            pb_max_difference = float("inf")
        pb_total = float(new_counts.to_numpy(dtype=np.float64).sum())
        checks["pseudobulk_structure_exact"] = (
            list(new_counts.shape) == [24505, 72]
            and sample_metadata_exact
            and genes_exact
        )
        checks["pseudobulk_values_exact"] = (
            pb_mismatches == 0
            and pb_max_difference == 0
            and pb_total == 481445434.0
        )
        metrics.update({
            "shape": list(replay.shape),
            "hvg_replay": len(replay_hvg),
            "hvg_historical": len(history_hvg),
            "hvg_intersection": hvg_intersection,
            "hvg_jaccard": hvg_intersection / hvg_union,
            "clusters_replay": int(replay_raw.nunique()),
            "ari": ari,
            "nmi": nmi,
            "hungarian_cell_match_rate": aligned_rate,
            "cell_type_match_rate": cell_type_match,
            "counts_mismatches": {
                "historical": counts_history,
                "scored": counts_scored,
                "annotated": counts_annotated,
            },
            "normalized_X_mismatches": {
                "historical": normalized_history,
                "scored": normalized_scored,
                "annotated": normalized_annotated,
            },
            "pseudobulk_shape": list(new_counts.shape),
            "pseudobulk_total_umi": pb_total,
            "pseudobulk_value_mismatches": pb_mismatches,
            "pseudobulk_max_absolute_difference": pb_max_difference,
        })
    finally:
        for value in objects:
            value.file.close()
    return checks, metrics


def incomplete_status(attempt_status: object, attempt_name: str) -> str:
    if not isinstance(attempt_status, dict):
        return "FAIL"
    if attempt_status.get("attempt") != attempt_name:
        return "FAIL"
    return "HOLD" if attempt_status.get("status") == "HOLD" else "FAIL"


def evaluate(path: Path) -> dict[str, object]:
    attempt = validate_attempt_path(path)
    missing = sorted(name for name in REQUIRED_OUTPUTS if not (attempt / name).is_file())
    status_path = attempt / "attempt_status.json"
    if missing:
        try:
            upstream = json_file(status_path) if status_path.is_file() else None
        except Exception as exc:
            upstream = {"parse_error": repr(exc)}
        status = incomplete_status(upstream, attempt.name)
        return {
            "status": status,
            "gate": "03A_historical_technical_replay",
            "attempt": attempt.name,
            "evaluated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "missing_outputs": missing,
            "upstream_attempt_status": upstream,
            "checks": {
                "complete_output_set": False,
                "incomplete_attempt_is_explicit_hold": status == "HOLD",
                "incomplete_attempt_can_pass": False,
            },
            "gate_03_overall": "HOLD_PENDING_GATE_03A_AND_GATE_03B",
            "evidence_boundary": "Missing replay outputs can never receive PASS.",
        }
    if status_path.exists():
        return {
            "status": "FAIL",
            "gate": "03A_historical_technical_replay",
            "attempt": attempt.name,
            "evaluated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "checks": {
                "complete_output_set": True,
                "no_conflicting_attempt_status": False,
            },
            "error": "A completed artifact set conflicts with attempt_status.json",
            "gate_03_overall": "HOLD_PENDING_GATE_03A_AND_GATE_03B",
        }
    checks, metrics = validate_complete_attempt(attempt)
    status = "PASS" if checks and all(checks.values()) else "FAIL"
    return {
        "status": status,
        "gate": "03A_historical_technical_replay",
        "attempt": attempt.name,
        "evaluated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks": {"complete_output_set": True, **checks},
        "metrics": metrics,
        "gate_03_overall": (
            "HOLD_PENDING_GATE_03B_DOUBLET_SENSITIVITY"
            if status == "PASS"
            else "HOLD_PENDING_GATE_03A_AND_GATE_03B"
        ),
        "evidence_boundary": (
            "A PASS is technical equivalence to a result-known historical workflow, "
            "not independent biological validation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("attempt_dir", type=Path)
    args = parser.parse_args()
    try:
        report = evaluate(args.attempt_dir)
    except Exception as exc:
        report = {
            "status": "FAIL",
            "gate": "03A_historical_technical_replay",
            "attempt": args.attempt_dir.name,
            "evaluated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "gate_03_overall": "HOLD_PENDING_GATE_03A_AND_GATE_03B",
        }
    report["checker"] = {
        "path": str(Path(__file__).resolve()),
        "sha256": sha256(Path(__file__).resolve()),
    }
    payload = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    record = write_numbered_attempt("gate_03A_replay_attempt", report["status"], payload)
    atomic_replace_json(VALIDATION / "gate_03A_replay_latest.json", payload)
    report["validation_record"] = str(record)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["status"] == "PASS" else 2)


if __name__ == "__main__":
    main()
