#!/usr/bin/env python3
"""Guarded historical replay of the GSE206528 Scanpy chain.

The replay is deliberately isolated from v1.  It starts from the frozen raw
AnnData only after the Gate 03 input manifest has been verified, preserves the
historical Scanpy/Leiden semantics, aligns replayed cluster labels to the
historical labels before applying the result-known manual cell-type map, and
writes to a new immutable attempt directory.

This script produces replay artifacts.  Gate status is assigned separately by
``check_GSE206528_replay_v2.py``.  A sample-aware doublet sensitivity is Gate
03B and is not part of the historical scrublet=false replay.
"""

from __future__ import annotations

import argparse
import copy
import csv
import ctypes
import gc
import hashlib
import json
import os
import platform
import shutil
import sys
import traceback
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
ATTEMPT_ROOT = REVISION / "01_work/GSE206528_replay"
DEFAULT_INPUT = ROOT / "analysis/data/processed/GSE206528_raw.h5ad"
HISTORICAL_CLUSTERED = ROOT / "analysis/data/processed/GSE206528_clustered.h5ad"
CONFIG = ROOT / "analysis/config/GSE206528_scanpy.json"
SIGNATURES = ROOT / "analysis/config/GSE206528_celltype_signatures.json"
MAPPING = ROOT / "analysis/config/GSE206528_cluster_celltype_mapping.json"
SAMPLE_DESIGN = ROOT / "analysis/data/metadata/core_sample_design.tsv"
GATE03_MANIFEST = REVISION / "00_protocol/gate_03_provenance_manifest.sha256.tsv"
MATRIX_TO_RAW_VALIDATION = (
    REVISION / "04_validation/gate_03_matrix_to_raw_validation_attempt_02_PASS.json"
)
EXPECTED_MANIFEST_SHA256 = "c832788e13681be8c2d0d73e3d46cbbaa8e7b1eb4f579a1faaff70b90671043b"
EXPECTED_MANIFEST_ENTRIES = 31
EXPECTED_MATRIX_VALIDATION_SHA256 = "7883dee813a7549aa5ea2c4c94c9cb47e40be146b17269c291bd67832de65795"
MINIMUM_FREE_MEMORY_GB = 10.0
MINIMUM_FREE_DISK_GB = 15.0
REQUIRED_MANIFEST_COLUMNS = {"path_kind", "path", "role", "bytes", "sha256"}
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
EXTERNAL_SNAPSHOT_ENV = "GSE206528_SCANPY_SNAPSHOT_DIR"

HARMONY_PARAMETERS = {
    "random_state": 0,
    "sigma": 0.1,
    "tau": 0,
    "block_size": 0.05,
    "max_iter_harmony": 10,
    "max_iter_kmeans": 20,
    "epsilon_cluster": 1e-5,
    "epsilon_harmony": 1e-4,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    temporary.replace(path)


def resolve_manifest_target(path_kind: str, recorded_path: str) -> Path:
    if path_kind == "workspace":
        target = (ROOT / recorded_path).resolve()
        try:
            target.relative_to(ROOT.resolve())
        except ValueError as exc:
            raise RuntimeError(f"Workspace manifest path escapes the workspace: {recorded_path}") from exc
        return target
    if path_kind == "absolute":
        target = Path(recorded_path)
        if not target.is_absolute():
            raise RuntimeError(f"Manifest absolute path is not absolute: {recorded_path}")
        if not target.is_file():
            snapshot_value = os.environ.get(EXTERNAL_SNAPSHOT_ENV)
            if snapshot_value:
                snapshot_root = Path(snapshot_value).resolve()
                candidate = (snapshot_root / target.name).resolve()
                try:
                    candidate.relative_to(snapshot_root)
                except ValueError as exc:
                    raise RuntimeError(
                        f"External-script snapshot escapes its root: {candidate}"
                    ) from exc
                return candidate
        return target.resolve()
    raise RuntimeError(f"Unsupported manifest path_kind: {path_kind}")


def verify_manifest(path: Path) -> dict[str, object]:
    if not path.exists():
        raise RuntimeError(f"Gate 03 manifest is absent: {path}")
    manifest_hash = sha256(path)
    if manifest_hash != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError(
            f"Gate 03 manifest identity drift: {manifest_hash} != {EXPECTED_MANIFEST_SHA256}"
        )
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != REQUIRED_MANIFEST_COLUMNS:
            raise RuntimeError(
                f"Gate 03 manifest columns must be exactly {sorted(REQUIRED_MANIFEST_COLUMNS)}"
            )
        rows = list(reader)
    if len(rows) != EXPECTED_MANIFEST_ENTRIES:
        raise RuntimeError(
            f"Gate 03 manifest must contain {EXPECTED_MANIFEST_ENTRIES} entries, observed {len(rows)}"
        )
    failures: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    workspace_entries = 0
    absolute_entries = 0
    for row in rows:
        identity = (row["path_kind"], row["path"])
        if identity in seen:
            failures.append({"path": row["path"], "reason": "duplicate"})
            continue
        seen.add(identity)
        target = resolve_manifest_target(*identity)
        workspace_entries += row["path_kind"] == "workspace"
        absolute_entries += row["path_kind"] == "absolute"
        if not target.is_file():
            failures.append({"path": row["path"], "reason": "missing"})
            continue
        expected_bytes = int(row["bytes"])
        observed_bytes = target.stat().st_size
        if observed_bytes != expected_bytes:
            failures.append(
                {
                    "path": row["path"],
                    "reason": "byte_mismatch",
                    "expected": expected_bytes,
                    "observed": observed_bytes,
                }
            )
        expected_hash = row["sha256"].lower()
        observed_hash = sha256(target)
        if observed_hash != expected_hash:
            failures.append(
                {
                    "path": row["path"],
                    "reason": "hash_mismatch",
                    "expected": expected_hash,
                    "observed": observed_hash,
                }
            )
    if failures:
        raise RuntimeError(f"Gate 03 frozen-input verification failed: {json.dumps(failures)}")
    return {
        "manifest": str(path.resolve()),
        "manifest_sha256": manifest_hash,
        "entries": len(rows),
        "workspace_entries": workspace_entries,
        "absolute_entries": absolute_entries,
        "verified": True,
    }


def verify_matrix_to_raw_record(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"Frozen matrix-to-raw validation record is absent: {path}")
    record_hash = sha256(path)
    if record_hash != EXPECTED_MATRIX_VALIDATION_SHA256:
        raise RuntimeError(
            f"Matrix-to-raw validation identity drift: {record_hash} != {EXPECTED_MATRIX_VALIDATION_SHA256}"
        )
    report = json.loads(path.read_text(encoding="utf-8"))
    required_checks = {
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
    checks = report.get("checks", {})
    totals = report.get("totals", {})
    expected_totals = {
        "cells": 64993,
        "genes_per_matrix": 33694,
        "source_vs_raw_mismatches": 0,
        "raw_X_vs_counts_layer_mismatches": 0,
        "max_absolute_difference": 0.0,
        "source_total_umi": 481449114,
        "raw_total_umi": 481449114,
        "counts_layer_total_umi": 481449114,
    }
    if (
        report.get("status") != "PASS"
        or set(checks) != required_checks
        or not all(checks.values())
        or totals != expected_totals
        or len(report.get("archive_members", [])) != 8
        or len(report.get("per_sample", [])) != 8
    ):
        raise RuntimeError("Frozen matrix-to-raw PASS record does not meet the complete acceptance contract")
    return {
        "record": str(path.resolve()),
        "record_sha256": record_hash,
        "status": report["status"],
        "checks": checks,
        "totals": totals,
    }


def artifact_record(path: Path, attempt: Path) -> dict[str, object]:
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(attempt.resolve())
    except ValueError as exc:
        raise RuntimeError(f"Output artifact escapes replay attempt: {path}") from exc
    if not resolved.is_file():
        raise RuntimeError(f"Required replay artifact is absent: {path}")
    return {
        "relative_path": relative.as_posix(),
        "bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def available_memory_gb() -> float | None:
    try:
        import psutil

        return psutil.virtual_memory().available / 1024**3
    except Exception:
        pass
    if platform.system() == "Windows":
        class MemoryStatusEx(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        status = MemoryStatusEx()
        status.dwLength = ctypes.sizeof(status)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            return status.ullAvailPhys / 1024**3
    if platform.system() == "Linux":
        try:
            fields = {}
            for line in Path("/proc/meminfo").read_text().splitlines():
                key, value = line.split(":", 1)
                fields[key] = float(value.strip().split()[0])
            return fields["MemAvailable"] / 1024**2
        except Exception:
            return None
    return None


def next_attempt_dir(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    numbers = []
    for path in root.glob("attempt_*"):
        try:
            numbers.append(int(path.name.split("_")[1]))
        except (IndexError, ValueError):
            continue
    number = max(numbers, default=0) + 1
    while True:
        attempt = root / f"attempt_{number:02d}"
        try:
            attempt.mkdir(parents=False, exist_ok=False)
            return attempt
        except FileExistsError:
            number += 1


def atomic_write_h5ad(adata: ad.AnnData, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    adata.strings_to_categoricals()
    adata.write_h5ad(temporary, compression="gzip")
    check = ad.read_h5ad(temporary, backed="r")
    shape = check.shape
    check.file.close()
    if shape != adata.shape:
        raise RuntimeError(f"Written object shape mismatch for {path}: {shape} != {adata.shape}")
    temporary.replace(path)


def canonicalize_metadata(adata: ad.AnnData, out: Path) -> None:
    design = pd.read_csv(SAMPLE_DESIGN, sep="\t", dtype=str)
    design = design.loc[design["series_accession"].eq("GSE206528")].copy()
    if len(design) != 8 or not design["sample_accession"].is_unique:
        raise RuntimeError("Expected eight unique GSE206528 accessions in core_sample_design.tsv")
    expected = set(design["sample_accession"])
    observed = set(adata.obs["sample_accession"].astype(str))
    if observed != expected:
        raise RuntimeError(f"Raw AnnData accessions differ from sample design: {observed ^ expected}")
    donor_map = design.set_index("sample_accession")["subject"].to_dict()
    condition_map = design.set_index("sample_accession")["disease"].to_dict()
    legacy = adata.obs[["sample_accession", "donor", "condition"]].astype(str).drop_duplicates()
    legacy = legacy.rename(columns={"donor": "legacy_donor", "condition": "legacy_condition"})
    legacy["canonical_donor"] = legacy["sample_accession"].map(donor_map)
    legacy["canonical_condition"] = legacy["sample_accession"].map(condition_map)
    legacy.to_csv(out / "legacy_to_canonical_donor_crosswalk.tsv", sep="\t", index=False)
    accession = adata.obs["sample_accession"].astype(str)
    adata.obs["donor_legacy"] = adata.obs["donor"].astype(str).values
    adata.obs["donor"] = accession.map(donor_map).values
    adata.obs["condition"] = accession.map(condition_map).values


def sample_design_snapshot() -> dict[str, dict[str, str]]:
    """Return the canonical accession map used by v2.

    The eight rows are frozen in the 31-item Gate 03 manifest and are also
    captured explicitly in the replay provenance.
    """
    design = pd.read_csv(SAMPLE_DESIGN, sep="\t", dtype=str)
    design = design.loc[design["series_accession"].eq("GSE206528")].copy()
    if len(design) != 8 or not design["sample_accession"].is_unique:
        raise RuntimeError("Expected eight unique GSE206528 sample-design rows")
    return {
        row["sample_accession"]: {"donor": row["subject"], "condition": row["disease"]}
        for _, row in design.iterrows()
    }


def align_clusters(replayed: pd.Series, historical: pd.Series) -> tuple[pd.Series, dict[str, object]]:
    if not replayed.index.equals(historical.index):
        raise RuntimeError("Replayed and historical cell indices differ before cluster alignment")
    replay_labels = sorted(replayed.astype(str).unique())
    history_labels = sorted(historical.astype(str).unique())
    contingency = pd.crosstab(replayed.astype(str), historical.astype(str)).reindex(
        index=replay_labels, columns=history_labels, fill_value=0
    )
    row, column = linear_sum_assignment(-contingency.to_numpy())
    mapping = {replay_labels[i]: history_labels[j] for i, j in zip(row, column, strict=True)}
    aligned = replayed.astype(str).map(mapping)
    matched = aligned.eq(historical.astype(str))
    metrics = {
        "replayed_clusters": len(replay_labels),
        "historical_clusters": len(history_labels),
        "ari": float(adjusted_rand_score(historical.astype(str), replayed.astype(str))),
        "nmi": float(normalized_mutual_info_score(historical.astype(str), replayed.astype(str))),
        "hungarian_cell_match_rate": float(matched.mean()),
        "replay_to_historical_cluster": mapping,
        "unmapped_replay_clusters": sorted(set(replay_labels) - set(mapping)),
        "contingency": contingency.to_dict(),
    }
    if metrics["unmapped_replay_clusters"]:
        raise RuntimeError(f"Unmapped replay clusters: {metrics['unmapped_replay_clusters']}")
    return aligned, metrics


def export_pseudobulk(adata: ad.AnnData, out: Path) -> None:
    grouping = ["sample_accession", "donor", "condition", "cell_type"]
    metadata = adata.obs[grouping].astype(str).drop_duplicates().sort_values(grouping)
    if len(metadata) != 72 or metadata.duplicated(grouping).any():
        raise RuntimeError(f"Expected 72 unique pseudobulk groups, observed {len(metadata)}")
    keys = [tuple(row[column] for column in grouping) for _, row in metadata.iterrows()]
    key_to_code = {key: index for index, key in enumerate(keys)}
    cell_keys = [tuple(row[column] for column in grouping) for _, row in adata.obs[grouping].astype(str).iterrows()]
    codes = np.fromiter((key_to_code[key] for key in cell_keys), dtype=np.int64, count=adata.n_obs)
    design = sparse.csr_matrix(
        (np.ones(adata.n_obs, dtype=np.int64), (codes, np.arange(adata.n_obs))),
        shape=(len(keys), adata.n_obs),
    )
    counts = adata.layers["counts"]
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    aggregated = (design @ counts).tocsr()
    sample_ids = ["|".join(key) for key in keys]
    if len(sample_ids) != len(set(sample_ids)):
        raise RuntimeError("Canonical pseudobulk IDs are not unique")
    pd.DataFrame(keys, columns=grouping, index=sample_ids).to_csv(out / "pseudobulk_samples.csv")
    table = pd.DataFrame.sparse.from_spmatrix(aggregated.T, index=adata.var_names, columns=sample_ids)
    table.index.name = "gene_symbol"
    table.to_csv(out / "pseudobulk_counts.csv")


def export_evidence(adata: ad.AnnData, out: Path) -> None:
    qc = (
        adata.obs.groupby(["sample_accession", "donor", "condition"], observed=True)
        .agg(
            cells=("sample_accession", "size"),
            n_genes_min=("n_genes_by_counts", "min"),
            n_genes_median=("n_genes_by_counts", "median"),
            n_genes_max=("n_genes_by_counts", "max"),
            total_counts_median=("total_counts", "median"),
            pct_mt_median=("pct_counts_mt", "median"),
            pct_mt_max=("pct_counts_mt", "max"),
        )
        .reset_index()
    )
    qc.to_csv(out / "per_donor_QC_v2.tsv", sep="\t", index=False)
    counts = (
        adata.obs.groupby(["sample_accession", "donor", "condition", "cell_type"], observed=True)
        .size().rename("cell_count").reset_index()
    )
    counts["sample_total"] = counts.groupby("sample_accession")["cell_count"].transform("sum")
    counts["cell_fraction"] = counts["cell_count"] / counts["sample_total"]
    counts.to_csv(out / "celltype_counts_by_donor_v2.tsv", sep="\t", index=False)


def package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--attempt-root", type=Path, default=ATTEMPT_ROOT)
    parser.add_argument("--minimum-free-memory-gb", type=float, default=10.0)
    parser.add_argument("--minimum-free-disk-gb", type=float, default=15.0)
    args = parser.parse_args()

    if args.input.resolve() != DEFAULT_INPUT.resolve():
        parser.error(f"Gate 03 replay input is frozen to {DEFAULT_INPUT}")
    if args.attempt_root.resolve() != ATTEMPT_ROOT.resolve():
        parser.error(f"Gate 03 replay output is restricted to {ATTEMPT_ROOT}")
    attempt = next_attempt_dir(args.attempt_root)
    started = utc_now()
    try:
        if args.minimum_free_memory_gb < MINIMUM_FREE_MEMORY_GB:
            raise RuntimeError(
                f"Protocol violation: memory gate cannot be below {MINIMUM_FREE_MEMORY_GB:.2f} GB"
            )
        if args.minimum_free_disk_gb < MINIMUM_FREE_DISK_GB:
            raise RuntimeError(
                f"Protocol violation: disk gate cannot be below {MINIMUM_FREE_DISK_GB:.2f} GB"
            )
        manifest = verify_manifest(GATE03_MANIFEST)
        matrix_to_raw = verify_matrix_to_raw_record(MATRIX_TO_RAW_VALIDATION)
        cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
        expected_config = {
            "scrublet": False,
            "scale": False,
            "regress_out": None,
            "batch_method": "harmony",
            "skip_markers": False,
            "n_pcs": 40,
            "n_neighbors": 15,
            "resolution": 0.6,
            "batch_key": "sample_accession",
            "hvg_flavor": "seurat",
            "n_top_genes": 3000,
        }
        for key, expected in expected_config.items():
            if cfg.get(key) != expected:
                raise RuntimeError(f"Historical configuration drift for {key}: {cfg.get(key)!r} != {expected!r}")
        observed_versions = {name: package_version(name) for name in EXPECTED_VERSIONS}
        if observed_versions != EXPECTED_VERSIONS:
            raise RuntimeError(
                f"Frozen replay environment drift: {observed_versions!r} != {EXPECTED_VERSIONS!r}"
            )
        free_memory = available_memory_gb()
        free_disk = shutil.disk_usage(attempt).free / 1024**3
        if free_memory is None or free_memory < args.minimum_free_memory_gb:
            raise RuntimeError(
                f"HOLD: available RAM {free_memory!r} GB; at least {args.minimum_free_memory_gb:.2f} GB required"
            )
        if free_disk < args.minimum_free_disk_gb:
            raise RuntimeError(
                f"HOLD: available disk {free_disk:.2f} GB; at least {args.minimum_free_disk_gb:.2f} GB required"
            )

        signatures = json.loads(SIGNATURES.read_text(encoding="utf-8"))
        cell_map = {str(key): value for key, value in json.loads(MAPPING.read_text(encoding="utf-8")).items()}
        figures = attempt / "figures"
        markers = attempt / "markers"
        figures.mkdir()
        markers.mkdir()
        np.random.seed(0)
        sc.settings.verbosity = 3
        sc.settings.figdir = figures
        sc.settings.set_figure_params(dpi=100, facecolor="white")

        design_snapshot = sample_design_snapshot()
        adata = ad.read_h5ad(args.input)
        canonicalize_metadata(adata, attempt)
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-", "Mt-"))
        qc_vars = ["mt"] if adata.var["mt"].any() else []
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True)
        sc.pp.filter_cells(adata, min_genes=int(cfg["min_genes"]))
        sc.pp.filter_genes(adata, min_cells=int(cfg["min_cells"]))
        adata = adata[adata.obs["n_genes_by_counts"] < int(cfg["max_genes"]), :].copy()
        if qc_vars:
            adata = adata[adata.obs["pct_counts_mt"] < float(cfg["mt_threshold"]), :].copy()

        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=float(cfg["target_sum"]))
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata, n_top_genes=int(cfg["n_top_genes"]), flavor=cfg["hvg_flavor"], batch_key=cfg["batch_key"]
        )
        adata.raw = adata

        hvg_mask = np.asarray(adata.var["highly_variable"]).astype(bool)
        work = ad.AnnData(
            X=adata.X[:, hvg_mask].copy(),
            obs=adata.obs.copy(),
            var=adata.var.loc[hvg_mask].copy(),
        )
        sc.tl.pca(work, svd_solver="arpack", random_state=0)  # historical call; yields 50 components here
        sc.external.pp.harmony_integrate(work, cfg["batch_key"], **HARMONY_PARAMETERS)
        sc.pp.neighbors(
            work, n_neighbors=int(cfg["n_neighbors"]), n_pcs=int(cfg["n_pcs"]),
            use_rep="X_pca_harmony", random_state=0,
        )
        sc.tl.umap(work, random_state=0)
        # Exact semantics recovered from the 2026-07-12 historical session.
        sc.tl.leiden(
            work, resolution=float(cfg["resolution"]), flavor="igraph",
            n_iterations=2, directed=False, random_state=0,
        )

        historical = ad.read_h5ad(HISTORICAL_CLUSTERED, backed="r")
        historical_labels = historical.obs["leiden"].copy()
        if not work.obs_names.equals(historical.obs_names):
            historical.file.close()
            raise RuntimeError("Cell order differs from the historical clustered object")
        aligned, alignment = align_clusters(work.obs["leiden"], historical_labels)
        historical.file.close()
        atomic_write_json(attempt / "cluster_alignment_v2.json", alignment)

        adata.obs["leiden_replay"] = work.obs["leiden"].astype(str).values
        adata.obs["leiden"] = pd.Categorical(aligned.values, categories=sorted(cell_map))
        adata.obsm["X_pca"] = np.asarray(work.obsm["X_pca"]).copy()
        adata.obsm["X_pca_harmony"] = np.asarray(work.obsm["X_pca_harmony"]).copy()
        adata.obsm["X_umap"] = np.asarray(work.obsm["X_umap"]).copy()
        adata.uns["neighbors"] = copy.deepcopy(work.uns["neighbors"])
        for key, value in work.obsp.items():
            adata.obsp[key] = value.copy()
        del work
        gc.collect()

        sc.tl.rank_genes_groups(adata, "leiden", method=cfg["marker_method"], use_raw=True)
        marker_frames = []
        for cluster in adata.obs["leiden"].cat.categories:
            table = sc.get.rank_genes_groups_df(adata, group=cluster).head(25)
            table.insert(0, "historical_aligned_cluster", cluster)
            marker_frames.append(table)
        pd.concat(marker_frames, ignore_index=True).to_csv(markers / "markers_all_v2.csv", index=False)
        atomic_write_h5ad(adata, attempt / "GSE206528_clustered_v2.h5ad")

        coverage = []
        for name, genes in signatures.items():
            present = [gene for gene in genes if gene in adata.raw.var_names]
            missing = [gene for gene in genes if gene not in adata.raw.var_names]
            coverage.append({
                "signature": name, "declared": len(genes), "present": len(present),
                "missing": len(missing), "coverage": len(present) / len(genes),
                "present_genes": ";".join(present), "missing_genes": ";".join(missing),
            })
            if missing:
                raise RuntimeError(f"Signature {name} has missing genes: {missing}")
            sc.tl.score_genes(adata, present, score_name=f"{name}_score", use_raw=True, random_state=0)
        pd.DataFrame(coverage).to_csv(attempt / "signature_gene_coverage_v2.tsv", sep="\t", index=False)
        atomic_write_h5ad(adata, attempt / "GSE206528_scored_v2.h5ad")

        adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(cell_map).astype("category")
        if adata.obs["cell_type"].isna().any():
            raise RuntimeError("Frozen cell-type map does not cover all aligned clusters")
        atomic_write_h5ad(adata, attempt / "GSE206528_annotated_v2.h5ad")
        export_pseudobulk(adata, attempt)
        export_evidence(adata, attempt)

        artifact_paths = {
            "clustered_h5ad": attempt / "GSE206528_clustered_v2.h5ad",
            "scored_h5ad": attempt / "GSE206528_scored_v2.h5ad",
            "annotated_h5ad": attempt / "GSE206528_annotated_v2.h5ad",
            "pseudobulk_counts": attempt / "pseudobulk_counts.csv",
            "pseudobulk_samples": attempt / "pseudobulk_samples.csv",
            "cluster_alignment": attempt / "cluster_alignment_v2.json",
            "signature_coverage": attempt / "signature_gene_coverage_v2.tsv",
            "donor_crosswalk": attempt / "legacy_to_canonical_donor_crosswalk.tsv",
            "per_donor_qc": attempt / "per_donor_QC_v2.tsv",
            "celltype_counts": attempt / "celltype_counts_by_donor_v2.tsv",
            "markers": markers / "markers_all_v2.csv",
        }
        output_artifacts = {
            name: artifact_record(path, attempt) for name, path in artifact_paths.items()
        }
        provenance = {
            "status": "ARTIFACTS_READY_FOR_ACCEPTANCE_CHECK",
            "started_at": started,
            "completed_at": utc_now(),
            "attempt": attempt.name,
            "attempt_root": str(ATTEMPT_ROOT.resolve()),
            "command": sys.argv,
            "runner": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256(Path(__file__).resolve()),
            },
            "frozen_inputs": manifest,
            "matrix_to_raw_semantic_lineage": matrix_to_raw,
            "input": {"path": str(args.input.resolve()), "sha256": sha256(args.input)},
            "canonical_sample_design": {
                "path": str(SAMPLE_DESIGN.resolve()), "sha256": sha256(SAMPLE_DESIGN), "rows": design_snapshot
            },
            "config": cfg,
            "harmony_parameters": HARMONY_PARAMETERS,
            "leiden_parameters": {"flavor": "igraph", "n_iterations": 2, "directed": False, "random_state": 0},
            "cluster_alignment": alignment,
            "protocol_resource_minimums_gb": {
                "available_memory": MINIMUM_FREE_MEMORY_GB,
                "available_disk": MINIMUM_FREE_DISK_GB,
            },
            "requested_resource_thresholds_gb": {
                "available_memory": args.minimum_free_memory_gb,
                "available_disk": args.minimum_free_disk_gb,
            },
            "resources_at_start": {
                "available_memory_gb": free_memory,
                "available_disk_gb": free_disk,
            },
            "versions": observed_versions,
            "python": sys.version,
            "python_executable": sys.executable,
            "thread_environment": {key: os.environ.get(key) for key in [
                "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"
            ]},
            "output_artifacts": output_artifacts,
            "evidence_boundary": (
                "Result-known technical replay using a frozen manual annotation map; not independent biological validation."
            ),
        }
        atomic_write_json(attempt / "replay_provenance_v2.json", provenance)
        print(json.dumps(provenance, indent=2, ensure_ascii=False))
    except Exception as exc:
        failure = {
            "status": "HOLD" if str(exc).startswith("HOLD:") else "FAIL",
            "started_at": started,
            "stopped_at": utc_now(),
            "attempt": attempt.name,
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "command": sys.argv,
            "runner": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256(Path(__file__).resolve()),
            },
            "prechecks_completed": {
                "frozen_inputs": locals().get("manifest"),
                "matrix_to_raw_semantic_lineage": locals().get("matrix_to_raw"),
                "versions": locals().get("observed_versions"),
            },
            "resource_gate": {
                "protocol_minimum_memory_gb": MINIMUM_FREE_MEMORY_GB,
                "protocol_minimum_disk_gb": MINIMUM_FREE_DISK_GB,
                "requested_minimum_memory_gb": args.minimum_free_memory_gb,
                "requested_minimum_disk_gb": args.minimum_free_disk_gb,
                "available_memory_gb": locals().get("free_memory"),
                "available_disk_gb": locals().get("free_disk"),
            },
        }
        atomic_write_json(attempt / "attempt_status.json", failure)
        print(json.dumps(failure, indent=2, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
