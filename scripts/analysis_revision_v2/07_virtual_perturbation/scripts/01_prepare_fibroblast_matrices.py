#!/usr/bin/env python3
"""Prepare frozen, doublet-excluded fibroblast matrices for VP-G02.

The source H5AD is opened only in backed read-only mode.  This script never
reclusters, reannotates, alters source files, or reads virtual-KO results.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from scipy.io import mmwrite


MAIN_DONORS = (
    ("normal", "Normal_1", 2377),
    ("normal", "Normal_2", 3254),
    ("normal", "Normal_3", 7149),
    ("organic_ED_nonDM", "non-DM_1", 2042),
    ("organic_ED_nonDM", "non-DM_2", 1586),
    ("organic_ED_nonDM", "non-DM_3", 4854),
)
TARGETS = ("TYMS", "EFNB2", "LRRC17")
HALLMARK_HYPOXIA = "HALLMARK_HYPOXIA"
TOTAL_CELLS_EXPECTED = 64993
DOUBLET_CELLS_EXPECTED = 2653
RETAINED_CELLS_EXPECTED = 62340


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_mtx_gz(path: Path, matrix: sparse.spmatrix) -> None:
    temporary_mtx = path.with_name(f".{path.name}.matrix-{os.getpid()}")
    temporary_gz = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary_mtx.open("wb") as matrix_handle:
            mmwrite(matrix_handle, matrix.tocoo(), field="real", precision=8)
        with temporary_mtx.open("rb") as source, temporary_gz.open("wb") as raw_target:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_target, mtime=0) as compressed:
                shutil.copyfileobj(source, compressed)
        temporary_gz.replace(path)
    finally:
        temporary_mtx.unlink(missing_ok=True)
        temporary_gz.unlink(missing_ok=True)


def read_hallmark_hypoxia(gmt_path: Path) -> set[str]:
    with gmt_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == HALLMARK_HYPOXIA:
                return {gene.strip() for gene in fields[2:] if gene.strip()}
    raise ValueError(f"{HALLMARK_HYPOXIA} is absent from {gmt_path}")


def require_g01_pass(vp_root: Path) -> None:
    validation_path = vp_root / "validation" / "VP_G01_freeze_PASS.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if validation.get("status") != "VP_G01_freeze_PASS":
        raise RuntimeError("VP-G01 is not PASS; VP-G02 cannot start.")


def load_doublets(path: Path, obs: pd.DataFrame) -> pd.Series:
    calls = pd.read_csv(path, sep="\t", dtype={"cell_id": "string", "predicted_doublet": "string"}, keep_default_na=False)
    required = {"cell_id", "predicted_doublet", "donor", "condition", "cell_type"}
    missing = required.difference(calls.columns)
    if missing:
        raise ValueError(f"Doublet table missing required columns: {sorted(missing)}")
    if calls["cell_id"].isna().any() or calls["cell_id"].duplicated().any():
        raise ValueError("Doublet table cell_id must be present and unique.")
    if not calls["predicted_doublet"].isin(["True", "False"]).all():
        raise ValueError("predicted_doublet must contain only literal True or False values.")
    calls = calls.set_index("cell_id", verify_integrity=True)
    if set(calls.index) != set(obs.index):
        only_calls = len(set(calls.index).difference(obs.index))
        only_obs = len(set(obs.index).difference(calls.index))
        raise ValueError(f"Cell-barcode alignment failed: {only_calls} table-only and {only_obs} H5AD-only barcodes.")
    aligned = calls.reindex(obs.index)
    for column in ("donor", "condition", "cell_type"):
        if not (aligned[column].astype(str).to_numpy() == obs[column].astype(str).to_numpy()).all():
            raise ValueError(f"Doublet metadata differs from frozen H5AD for {column}.")
    return aligned["predicted_doublet"].eq("True")


def backed_to_csr(view, layer: str) -> sparse.csr_matrix:
    matrix = view.layers[layer]
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    return matrix.tocsr()


def normalized_log1p(counts: sparse.csr_matrix) -> sparse.csr_matrix:
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    if np.any(totals <= 0):
        raise ValueError("A retained fibroblast cell has non-positive total raw counts.")
    normalized = counts.astype(np.float64).multiply((10000.0 / totals)[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    return normalized


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    vp_root = project_root / "revision_v2" / "07_virtual_perturbation"
    prepared_root = vp_root / "04_prepared"
    validation_root = vp_root / "validation"
    h5ad_path = project_root / "revision_v2/01_work/GSE206528_replay/attempt_07/GSE206528_annotated_v2.h5ad"
    doublet_path = project_root / "revision_v2/01_work/GSE206528_doublet_sensitivity/attempt_01/doublet_calls_by_cell.tsv"
    hallmark_path = project_root / "analysis/data/processed/gene_sets/Hallmark.gmt"
    started = now()
    started_monotonic = time.monotonic()
    process = psutil.Process()
    initial_rss = process.memory_info().rss
    adata = None
    try:
        require_g01_pass(vp_root)
        if (prepared_root / "fibroblast").exists() or (prepared_root / "gene_universe.tsv").exists():
            raise FileExistsError("VP-G02 final artifacts already exist; refusing to overwrite a frozen attempt.")
        prepared_root.mkdir(parents=True, exist_ok=True)
        validation_root.mkdir(parents=True, exist_ok=True)
        staging = prepared_root / f".VP_G02_staging_{os.getpid()}"
        staging.mkdir(parents=False, exist_ok=False)

        adata = ad.read_h5ad(h5ad_path, backed="r")
        try:
            if adata.shape[0] != TOTAL_CELLS_EXPECTED:
                raise ValueError(f"H5AD cell count drift: {adata.shape[0]} != {TOTAL_CELLS_EXPECTED}")
            for required in ("donor", "condition", "cell_type"):
                if required not in adata.obs.columns:
                    raise ValueError(f"Frozen H5AD lacks obs.{required}")
            if "counts" not in adata.layers:
                raise ValueError("Frozen H5AD lacks layers['counts']")
            if "highly_variable" not in adata.var.columns:
                raise ValueError("Frozen H5AD lacks var.highly_variable")
            if not adata.obs_names.is_unique or not adata.var_names.is_unique:
                raise ValueError("Frozen H5AD requires unique cell and gene identifiers.")
            doublet_mask = load_doublets(doublet_path, adata.obs)
            if int(doublet_mask.sum()) != DOUBLET_CELLS_EXPECTED:
                raise ValueError(f"Doublet count drift: {int(doublet_mask.sum())} != {DOUBLET_CELLS_EXPECTED}")
            retained_mask = ~doublet_mask.to_numpy()
            if int(retained_mask.sum()) != RETAINED_CELLS_EXPECTED:
                raise ValueError(f"Retained-cell count drift: {int(retained_mask.sum())} != {RETAINED_CELLS_EXPECTED}")

            obs = adata.obs.copy()
            var_names = np.asarray(adata.var_names.astype(str))
            hvg_mask = adata.var["highly_variable"].fillna(False).astype(bool).to_numpy()
            if int(hvg_mask.sum()) != 3000:
                raise ValueError(f"HVG count drift: {int(hvg_mask.sum())} != 3000")
            hvg_indices = np.flatnonzero(hvg_mask)
            donor_positions: dict[tuple[str, str], np.ndarray] = {}
            detection_by_donor = np.zeros((len(MAIN_DONORS), len(hvg_indices)), dtype=np.int32)
            observed_donor_counts: dict[str, int] = {}
            for donor_i, (condition, donor, expected_count) in enumerate(MAIN_DONORS):
                donor_mask = retained_mask & (obs["condition"].astype(str).to_numpy() == condition) & (obs["donor"].astype(str).to_numpy() == donor) & (obs["cell_type"].astype(str).to_numpy() == "fibroblast")
                positions = np.flatnonzero(donor_mask)
                donor_positions[(condition, donor)] = positions
                observed_donor_counts[donor] = len(positions)
                if len(positions) != expected_count:
                    raise ValueError(f"Fibroblast count drift for {donor}: {len(positions)} != {expected_count}")
                hvg_counts = backed_to_csr(adata[positions, hvg_indices], "counts")
                detection_by_donor[donor_i, :] = np.asarray(hvg_counts.getnnz(axis=0)).ravel()
            selected_hvg_indices = hvg_indices[(detection_by_donor >= 25).sum(axis=0) >= 4]

            target_indices = []
            for target in TARGETS:
                matches = np.flatnonzero(var_names == target)
                if len(matches) != 1:
                    raise ValueError(f"Target {target} must occur exactly once in frozen H5AD genes.")
                target_indices.append(int(matches[0]))
            hypoxia_genes = read_hallmark_hypoxia(hallmark_path)
            all_six_positions = np.concatenate(list(donor_positions.values()))
            hypoxia_indices = np.flatnonzero(np.isin(var_names, list(hypoxia_genes)))
            hypoxia_detectable: set[int] = set()
            if len(hypoxia_indices):
                hypoxia_counts = backed_to_csr(adata[all_six_positions, hypoxia_indices], "counts")
                hypoxia_detectable = set(hypoxia_indices[np.asarray(hypoxia_counts.getnnz(axis=0)).ravel() > 0].tolist())
            universe_indices = np.array(sorted(set(selected_hvg_indices.tolist()).union(target_indices).union(hypoxia_detectable)), dtype=np.int64)
            if not set(target_indices).issubset(set(universe_indices.tolist())):
                raise ValueError("At least one target was lost from the frozen gene universe.")

            gene_rows: list[dict[str, object]] = []
            selected_hvg_set = set(selected_hvg_indices.tolist())
            target_set = set(target_indices)
            for order, index in enumerate(universe_indices, start=1):
                reasons = []
                if int(index) in selected_hvg_set:
                    reasons.append("HVG_DETECTED_GE25_IN_GE4_OF_6_DONORS")
                if int(index) in target_set:
                    reasons.append("FORCED_TARGET")
                if int(index) in hypoxia_detectable:
                    reasons.append("HALLMARK_HYPOXIA_DETECTABLE_IN_MAIN_FIBROBLASTS")
                gene_rows.append({"gene_order": order, "gene": var_names[index], "source_var_index_0based": int(index), "selection_reason": ";".join(reasons)})
            gene_path_staging = staging / "gene_universe.tsv"
            write_tsv(gene_path_staging, list(gene_rows[0]), gene_rows)
            gene_sha = sha256(gene_path_staging)
            source_h5ad_sha = sha256(h5ad_path)
            source_doublet_sha = sha256(doublet_path)

            fibroblast_staging = staging / "fibroblast"
            fibroblast_staging.mkdir()
            target_detection_all: list[dict[str, object]] = []
            for condition, donor, expected_count in MAIN_DONORS:
                donor_dir = fibroblast_staging / condition / donor
                donor_dir.mkdir(parents=True)
                positions = donor_positions[(condition, donor)]
                raw_counts = backed_to_csr(adata[positions, universe_indices], "counts")
                normalized = normalized_log1p(raw_counts)
                write_mtx_gz(donor_dir / "matrix.mtx.gz", normalized.T)
                write_tsv(
                    donor_dir / "genes.tsv",
                    ["gene_order", "gene"],
                    [{"gene_order": row["gene_order"], "gene": row["gene"]} for row in gene_rows],
                )
                cells = obs.iloc[positions]
                cell_rows = [
                    {"cell_order": cell_order, "cell_id": cell_id, "donor": donor, "condition": condition, "cell_type": "fibroblast", "predicted_doublet": "False"}
                    for cell_order, cell_id in enumerate(cells.index.astype(str), start=1)
                ]
                write_tsv(donor_dir / "cells.tsv", list(cell_rows[0]), cell_rows)
                for target in TARGETS:
                    target_column = int(np.where(var_names[universe_indices] == target)[0][0])
                    detected_cells = int(raw_counts[:, target_column].getnnz())
                    primary = (target == "TYMS" and condition == "normal") or (target in {"EFNB2", "LRRC17"} and condition == "organic_ED_nonDM")
                    status = "PASS" if (not primary or detected_cells >= 25) else "TARGET_LOW_DETECTION"
                    row = {
                        "gene": target,
                        "condition": condition,
                        "donor": donor,
                        "fibroblast_cells": expected_count,
                        "detected_cells_raw_counts": detected_cells,
                        "detection_rate": f"{detected_cells / expected_count:.10f}",
                        "is_primary_gene_donor_combination": str(primary),
                        "minimum_detected_cells_required": 25 if primary else "NOT_PRIMARY",
                        "status": status,
                    }
                    target_detection_all.append(row)
                donor_target_rows = [row for row in target_detection_all if row["condition"] == condition and row["donor"] == donor]
                write_tsv(donor_dir / "target_detection.tsv", list(donor_target_rows[0]), donor_target_rows)
                provenance = {
                    "gate": "VP-G02",
                    "source_h5ad": str(h5ad_path),
                    "source_h5ad_sha256": source_h5ad_sha,
                    "source_doublet_calls": str(doublet_path),
                    "source_doublet_calls_sha256": source_doublet_sha,
                    "cell_type_label_source": "frozen H5AD obs.cell_type; no reclustering or reannotation",
                    "donor": donor,
                    "condition": condition,
                    "cell_count": expected_count,
                    "gene_universe_file": "../../../../gene_universe.tsv",
                    "gene_universe_sha256": gene_sha,
                    "matrix_orientation": "gene_by_cell",
                    "matrix_layer": "counts",
                    "normalization": "per-cell total-count normalization to 10000 followed by log1p",
                    "cell_id_order_sha256": hashlib.sha256("\n".join(cells.index.astype(str)).encode("utf-8")).hexdigest(),
                    "created_at": now(),
                }
                atomic_json(donor_dir / "provenance.json", provenance)

            primary_detection_rows = [row for row in target_detection_all if row["is_primary_gene_donor_combination"] == "True"]
            if len(primary_detection_rows) != 9 or any(row["status"] != "PASS" for row in primary_detection_rows):
                raise ValueError("One or more frozen primary target×donor combinations failed the >=25-cell detection gate.")
            write_tsv(staging / "target_detection_all.tsv", list(target_detection_all[0]), target_detection_all)
            (staging / "gene_universe.sha256").write_text(f"{gene_sha}  gene_universe.tsv\n", encoding="utf-8", newline="\n")

            final_fibroblast = prepared_root / "fibroblast"
            staging_fibroblast = staging / "fibroblast"
            staging_fibroblast.replace(final_fibroblast)
            for name in ("gene_universe.tsv", "gene_universe.sha256", "target_detection_all.tsv"):
                (staging / name).replace(prepared_root / name)
            staging.rmdir()
            elapsed_seconds = round(time.monotonic() - started_monotonic, 3)
            validation = {
                "gate": "VP-G02",
                "status": "VP_G02_matrix_PASS",
                "validated_at": now(),
                "checks": {
                    "source_h5ad_backed_read_only": True,
                    "source_h5ad_total_cells": int(adata.shape[0]),
                    "doublet_cells_excluded": int(doublet_mask.sum()),
                    "retained_cells": int(retained_mask.sum()),
                    "frozen_cell_type_label_used_without_reannotation": True,
                    "hvg_count": int(hvg_mask.sum()),
                    "fixed_gene_universe_count": int(len(universe_indices)),
                    "six_donor_fibroblast_counts": observed_donor_counts,
                    "nine_primary_target_detection_combinations_pass": True,
                    "all_six_matrices_share_gene_order": True,
                    "normalization": "counts_to_10000_then_log1p",
                },
                "artifacts": {
                    "gene_universe": "04_prepared/gene_universe.tsv",
                    "gene_universe_sha256": gene_sha,
                    "target_detection_all": "04_prepared/target_detection_all.tsv",
                    "fibroblast_matrix_root": "04_prepared/fibroblast",
                },
                "resource_usage": {
                    "process_rss_start_bytes": initial_rss,
                    "process_rss_end_bytes": process.memory_info().rss,
                    "elapsed_seconds_since_process_start": elapsed_seconds,
                },
            }
            atomic_json(validation_root / "VP_G02_matrix_PASS.json", validation)
            print(json.dumps({"gate": "VP-G02", "status": "VP_G02_matrix_PASS", "gene_universe_count": len(universe_indices)}, ensure_ascii=False))
        finally:
            if adata is not None:
                adata.file.close()
    except Exception as error:
        validation_root.mkdir(parents=True, exist_ok=True)
        atomic_json(
            validation_root / "VP_G02_matrix_FAIL.json",
            {"gate": "VP-G02", "status": "FAIL_RETAINED", "failed_at": now(), "error_type": type(error).__name__, "error": str(error)},
        )
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
