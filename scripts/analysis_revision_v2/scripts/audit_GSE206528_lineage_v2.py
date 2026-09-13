#!/usr/bin/env python3
"""Memory-bounded lineage audit for the historical GSE206528 AnnData products.

This audit verifies that the stored raw, clustered, scored, annotated, and
pseudobulk products form a coherent chain. It does not replay Harmony, UMAP,
Leiden, marker testing, or annotation and therefore cannot by itself pass the
full reproducibility gate.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
RESULTS = REVISION / "02_results"
VALIDATION = REVISION / "04_validation"
DATA = ROOT / "analysis/data/processed"
GSE_RESULTS = ROOT / "analysis/results/GSE206528"


def write_numbered_attempt(
    validation_dir: Path, prefix: str, status: str, payload: str
) -> Path:
    """Create the first unused numbered attempt without overwriting history."""
    attempt_number = 1
    while True:
        number = f"{attempt_number:02d}"
        if any(validation_dir.glob(f"{prefix}_{number}_*.json")):
            attempt_number += 1
            continue
        attempt = validation_dir / f"{prefix}_{number}_{status}.json"
        try:
            with attempt.open("x", encoding="utf-8", newline="\n") as handle:
                handle.write(payload)
            return attempt
        except FileExistsError:
            attempt_number += 1


def atomic_replace_json(path: Path, payload: str) -> None:
    """Replace a latest pointer atomically via a PID-tagged sibling file."""
    counter = 0
    while True:
        suffix = f".{os.getpid()}.{counter}.tmp"
        temporary = path.with_name(f".{path.name}{suffix}")
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


def text_hash(values) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def qc_from_backed_raw(path: Path, chunk_size: int = 1000) -> dict[str, object]:
    adata = ad.read_h5ad(path, backed="r")
    mt_mask = np.asarray(adata.var_names.str.startswith(("MT-", "mt-", "Mt-")))
    gene_ncells = np.zeros(adata.n_vars, dtype=np.int64)
    cell_ngenes = np.zeros(adata.n_obs, dtype=np.int64)
    cell_totals = np.zeros(adata.n_obs, dtype=np.float64)
    cell_mt = np.zeros(adata.n_obs, dtype=np.float64)
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        matrix = adata.X[start:stop]
        if not sparse.issparse(matrix):
            matrix = sparse.csr_matrix(matrix)
        matrix = matrix.tocsr()
        gene_ncells += np.asarray(matrix.getnnz(axis=0)).ravel()
        cell_ngenes[start:stop] = np.asarray(matrix.getnnz(axis=1)).ravel()
        cell_totals[start:stop] = np.asarray(matrix.sum(axis=1)).ravel()
        if mt_mask.any():
            cell_mt[start:stop] = np.asarray(matrix[:, mt_mask].sum(axis=1)).ravel()
    pct_mt = np.divide(
        cell_mt * 100.0,
        cell_totals,
        out=np.zeros_like(cell_mt),
        where=cell_totals > 0,
    )
    cell_pass = (cell_ngenes >= 200) & (cell_ngenes < 8000) & (pct_mt < 20)
    gene_pass = gene_ncells >= 3
    result = {
        "shape": [adata.n_obs, adata.n_vars],
        "obs_names_hash": text_hash(adata.obs_names),
        "var_names_hash": text_hash(adata.var_names),
        "cell_pass_count": int(cell_pass.sum()),
        "cell_fail_count": int((~cell_pass).sum()),
        "gene_pass_count": int(gene_pass.sum()),
        "filtered_gene_names": adata.var_names[gene_pass].astype(str).tolist(),
        "cell_n_genes_min": int(cell_ngenes.min()),
        "cell_n_genes_max": int(cell_ngenes.max()),
        "cell_pct_mt_max": float(pct_mt.max()),
    }
    adata.file.close()
    return result


def inspect_processed(path: Path) -> dict[str, object]:
    adata = ad.read_h5ad(path, backed="r")
    result = {
        "shape": [adata.n_obs, adata.n_vars],
        "obs_names_hash": text_hash(adata.obs_names),
        "var_names_hash": text_hash(adata.var_names),
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "obsp_keys": list(adata.obsp.keys()),
        "layers": [str(key) for key in adata.layers.keys()],
        "raw_shape": None if adata.raw is None else list(adata.raw.shape),
    }
    if "leiden" in adata.obs:
        result["leiden_clusters"] = int(adata.obs["leiden"].nunique())
        result["leiden_counts"] = {
            str(key): int(value)
            for key, value in adata.obs["leiden"].value_counts().sort_index().items()
        }
    if "highly_variable" in adata.var:
        result["highly_variable_genes"] = int(adata.var["highly_variable"].sum())
    if "X_pca" in adata.obsm:
        result["X_pca_shape"] = list(adata.obsm["X_pca"].shape)
    if "X_pca_harmony" in adata.obsm:
        result["X_pca_harmony_shape"] = list(adata.obsm["X_pca_harmony"].shape)
    if "X_umap" in adata.obsm:
        result["X_umap_shape"] = list(adata.obsm["X_umap"].shape)
    result["obs"] = adata.obs.copy()
    result["var_names"] = adata.var_names.astype(str).tolist()
    adata.file.close()
    return result


def audit_pseudobulk(annotated_path: Path) -> dict[str, object]:
    sample_table = pd.read_csv(GSE_RESULTS / "pseudobulk_samples.csv", index_col=0)
    legacy_counts = pd.read_csv(GSE_RESULTS / "pseudobulk_counts.csv", index_col=0)
    annotated = ad.read_h5ad(annotated_path, backed="r")

    grouping = ["sample_accession", "donor", "condition", "cell_type"]
    sample_keys = [tuple(str(row[column]) for column in grouping) for _, row in sample_table.iterrows()]
    key_to_code = {key: index for index, key in enumerate(sample_keys)}
    cell_keys = [tuple(str(row[column]) for column in grouping) for _, row in annotated.obs.iterrows()]
    missing_keys = sorted(set(cell_keys) - set(key_to_code))
    if missing_keys:
        raise RuntimeError(f"Annotated cells contain pseudobulk groups absent from sample table: {missing_keys}")
    codes = np.fromiter((key_to_code[key] for key in cell_keys), dtype=np.int64, count=len(cell_keys))

    design = sparse.csr_matrix(
        (np.ones(annotated.n_obs, dtype=np.int64), (codes, np.arange(annotated.n_obs))),
        shape=(len(sample_keys), annotated.n_obs),
    )
    counts_layer = annotated.layers["counts"]
    if not sparse.issparse(counts_layer):
        counts_layer = sparse.csr_matrix(counts_layer)
    aggregated = (design @ counts_layer).tocsr()

    gene_order_match = legacy_counts.index.astype(str).tolist() == annotated.var_names.astype(str).tolist()
    sample_order_match = legacy_counts.columns.astype(str).tolist() == sample_table.index.astype(str).tolist()
    expected = legacy_counts.to_numpy(dtype=np.float64).T
    difference = aggregated.astype(np.float64) - sparse.csr_matrix(expected)
    mismatch_count = int(difference.count_nonzero())
    max_absolute_difference = float(np.max(np.abs(difference.data))) if mismatch_count else 0.0

    observed_group_counts = np.bincount(codes, minlength=len(sample_keys))
    historical_cell_counts = pd.read_csv(GSE_RESULTS / "celltype_counts_by_donor.tsv", sep="\t")
    count_map = {
        tuple(str(row[column]) for column in grouping): int(row["cell_count"])
        for _, row in historical_cell_counts.iterrows()
    }
    expected_group_counts = np.array([count_map[key] for key in sample_keys], dtype=np.int64)
    cell_count_mismatch = int(np.sum(observed_group_counts != expected_group_counts))

    annotated.file.close()
    return {
        "genes": int(legacy_counts.shape[0]),
        "pseudobulk_samples": int(legacy_counts.shape[1]),
        "gene_order_match": bool(gene_order_match),
        "sample_order_match": bool(sample_order_match),
        "count_mismatch_cells": mismatch_count,
        "max_absolute_count_difference": max_absolute_difference,
        "celltype_group_count_mismatches": cell_count_mismatch,
        "total_counts_recomputed": float(aggregated.sum()),
        "total_counts_legacy": float(expected.sum()),
    }


def main() -> None:
    raw_path = DATA / "GSE206528_raw.h5ad"
    clustered_path = DATA / "GSE206528_clustered.h5ad"
    scored_path = DATA / "GSE206528_scored.h5ad"
    annotated_path = DATA / "GSE206528_annotated.h5ad"

    raw_qc = qc_from_backed_raw(raw_path)
    clustered = inspect_processed(clustered_path)
    scored = inspect_processed(scored_path)
    annotated = inspect_processed(annotated_path)

    mapping = json.loads(
        (ROOT / "analysis/config/GSE206528_cluster_celltype_mapping.json").read_text(encoding="utf-8")
    )
    mapped = annotated["obs"]["leiden"].astype(str).map(mapping)
    annotation_match = mapped.astype(str).equals(annotated["obs"]["cell_type"].astype(str))
    signature_names = json.loads(
        (ROOT / "analysis/config/GSE206528_celltype_signatures.json").read_text(encoding="utf-8")
    )
    required_scores = {f"{name}_score" for name in signature_names}
    scored_columns_present = required_scores.issubset(set(scored["obs_columns"]))
    pseudobulk = audit_pseudobulk(annotated_path)

    filtered_gene_hash = text_hash(raw_qc.pop("filtered_gene_names"))
    checks = {
        "raw_shape_64993_by_33694": raw_qc["shape"] == [64993, 33694],
        "all_64993_cells_pass_recorded_QC_thresholds": raw_qc["cell_pass_count"] == 64993,
        "gene_filter_yields_24505_genes": raw_qc["gene_pass_count"] == 24505,
        "clustered_shape_64993_by_24505": clustered["shape"] == [64993, 24505],
        "filtered_gene_order_matches_clustered": filtered_gene_hash == clustered["var_names_hash"],
        "obs_order_preserved_raw_to_clustered": raw_qc["obs_names_hash"] == clustered["obs_names_hash"],
        "obs_order_preserved_all_processed_objects": len(
            {clustered["obs_names_hash"], scored["obs_names_hash"], annotated["obs_names_hash"]}
        ) == 1,
        "var_order_preserved_all_processed_objects": len(
            {clustered["var_names_hash"], scored["var_names_hash"], annotated["var_names_hash"]}
        ) == 1,
        "HVG_count_3000": clustered.get("highly_variable_genes") == 3000,
        "PCA_50_and_Harmony_50": clustered.get("X_pca_shape") == [64993, 50]
        and clustered.get("X_pca_harmony_shape") == [64993, 50],
        "UMAP_two_dimensions": clustered.get("X_umap_shape") == [64993, 2],
        "Leiden_18_clusters": clustered.get("leiden_clusters") == 18,
        "all_signature_score_columns_present": scored_columns_present,
        "annotation_matches_frozen_cluster_map": annotation_match,
        "pseudobulk_gene_and_sample_order_match": pseudobulk["gene_order_match"]
        and pseudobulk["sample_order_match"],
        "pseudobulk_counts_exact": pseudobulk["count_mismatch_cells"] == 0
        and pseudobulk["max_absolute_count_difference"] == 0,
        "celltype_counts_exact": pseudobulk["celltype_group_count_mismatches"] == 0,
    }

    serializable_objects = {}
    for name, value in {"clustered": clustered, "scored": scored, "annotated": annotated}.items():
        serializable_objects[name] = {
            key: item
            for key, item in value.items()
            if key not in {"obs", "var_names"}
        }

    audit = {
        "status": "HOLD",
        "historical_lineage_checks_pass": all(checks.values()),
        "checks": checks,
        "raw_qc": raw_qc,
        "objects": serializable_objects,
        "pseudobulk": pseudobulk,
        "remaining_requirements": [
            "Replay QC, normalization, HVG, PCA, Harmony, neighbors, UMAP, Leiden, signature scoring, annotation, and pseudobulk in the isolated v2 namespace.",
            "Compare replayed cell, gene, cluster, annotation, and pseudobulk outputs with frozen acceptance criteria.",
            "Run and report a sample-aware doublet-detection sensitivity analysis; the historical configuration explicitly set scrublet=false.",
            "Record per-donor QC distributions and annotation evidence in shareable tables.",
        ],
        "evidence_boundary": (
            "Stored products are internally coherent and pseudobulk counts are traceable to the annotated counts layer, "
            "but product coherence is not an end-to-end replay. Gate 03 remains HOLD."
        ),
    }
    (RESULTS / "GSE206528_historical_lineage_audit_v2.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    VALIDATION.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(audit, indent=2, ensure_ascii=False) + "\n"
    write_numbered_attempt(
        VALIDATION,
        "gate_03_lineage_audit_attempt",
        audit["status"],
        payload,
    )
    atomic_replace_json(
        VALIDATION / "gate_03_validation_latest.json",
        payload,
    )
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
