"""Aggregate fixed Gate 03A cell types after frozen Gate 03B doublet exclusion."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import psutil
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
REPLAY = REVISION / "01_work" / "GSE206528_replay" / "attempt_07"
CALLING = REVISION / "01_work" / "GSE206528_doublet_sensitivity" / "attempt_01"
INPUT = REPLAY / "GSE206528_annotated_v2.h5ad"
CALLS = CALLING / "doublet_calls_by_cell.tsv"
CALLING_PROVENANCE = CALLING / "doublet_calling_provenance.json"
ORIGINAL_COUNTS = REPLAY / "pseudobulk_counts.csv"
ORIGINAL_SAMPLES = REPLAY / "pseudobulk_samples.csv"
PROTOCOL = REVISION / "00_protocol" / "GATE_03B_DOUBLET_SENSITIVITY_PROTOCOL.md"
ATTEMPT_ROOT = REVISION / "01_work" / "GSE206528_doublet_sensitivity"
MIN_FREE_MEMORY_GB = 8.0
MIN_FREE_DISK_GB = 10.0
GROUPING = ["sample_accession", "donor", "condition", "cell_type"]


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    temporary.replace(path)


def next_attempt() -> Path:
    numbers = []
    for candidate in ATTEMPT_ROOT.glob("attempt_*"):
        try:
            numbers.append(int(candidate.name.split("_")[1]))
        except (IndexError, ValueError):
            continue
    number = max(numbers, default=0) + 1
    while True:
        candidate = ATTEMPT_ROOT / f"attempt_{number:02d}"
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            number += 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--attempt-root", type=Path, default=ATTEMPT_ROOT)
    args = parser.parse_args()
    if args.attempt_root.resolve() != ATTEMPT_ROOT.resolve():
        parser.error("Gate 03B output is restricted to its v2 attempt root")
    attempt = next_attempt()
    started = now()
    try:
        for path in [INPUT, CALLS, CALLING_PROVENANCE, ORIGINAL_COUNTS, ORIGINAL_SAMPLES, PROTOCOL]:
            if not path.is_file():
                raise RuntimeError(f"Required frozen input is absent: {path}")
        memory = psutil.virtual_memory().available / 1024**3
        disk = shutil.disk_usage(attempt).free / 1024**3
        if memory < MIN_FREE_MEMORY_GB:
            raise RuntimeError(f"HOLD: available RAM {memory:.2f} GB is below {MIN_FREE_MEMORY_GB:.2f} GB")
        if disk < MIN_FREE_DISK_GB:
            raise RuntimeError(f"HOLD: available disk {disk:.2f} GB is below {MIN_FREE_DISK_GB:.2f} GB")
        calling_provenance = json.loads(CALLING_PROVENANCE.read_text(encoding="utf-8"))
        if calling_provenance.get("status") != "CALLING_COMPLETE_AWAITING_DOWNSTREAM_SENSITIVITY":
            raise RuntimeError("Gate 03B calling provenance is not a completed frozen calling stage")
        adata = ad.read_h5ad(INPUT)
        if "counts" not in adata.layers:
            raise RuntimeError("Gate 03A annotated object has no counts layer")
        calls = pd.read_csv(CALLS, sep="\t", dtype={"cell_id": str})
        if calls["cell_id"].duplicated().any() or set(calls["cell_id"]) != set(adata.obs_names.astype(str)):
            raise RuntimeError("Doublet calls do not cover Gate 03A cells exactly once")
        calls = calls.set_index("cell_id").reindex(adata.obs_names.astype(str))
        if calls["predicted_doublet"].isna().any():
            raise RuntimeError("Doublet calls are missing after exact cell alignment")
        raw_prediction = calls["predicted_doublet"]
        if raw_prediction.dtype == bool:
            predicted = raw_prediction.to_numpy()
        else:
            parsed = raw_prediction.astype(str).str.strip().str.lower().map({"true": True, "false": False})
            if parsed.isna().any():
                raise RuntimeError("Predicted-doublet column is not strict boolean data")
            predicted = parsed.to_numpy(dtype=bool)
        keep = ~predicted
        retained = adata[keep].copy()
        original_metadata = adata.obs[GROUPING].astype(str).drop_duplicates().sort_values(GROUPING)
        keys = [tuple(row[column] for column in GROUPING) for _, row in original_metadata.iterrows()]
        if len(keys) != 72 or len(set(keys)) != 72:
            raise RuntimeError("The 03A reference does not contain exactly 72 pseudobulk keys")
        key_to_code = {key: index for index, key in enumerate(keys)}
        cell_keys = [tuple(row[column] for column in GROUPING) for _, row in retained.obs[GROUPING].astype(str).iterrows()]
        codes = np.fromiter((key_to_code[key] for key in cell_keys), dtype=np.int64, count=retained.n_obs)
        design = sparse.csr_matrix((np.ones(retained.n_obs, dtype=np.int64), (codes, np.arange(retained.n_obs))), shape=(72, retained.n_obs))
        counts = retained.layers["counts"]
        if not sparse.issparse(counts):
            counts = sparse.csr_matrix(counts)
        aggregated = (design @ counts).tocsr()
        sample_ids = ["|".join(key) for key in keys]
        samples = pd.DataFrame(keys, columns=GROUPING, index=sample_ids)
        samples.to_csv(attempt / "pseudobulk_samples.csv")
        table = pd.DataFrame.sparse.from_spmatrix(aggregated.T, index=retained.var_names, columns=sample_ids)
        table.index.name = "gene_symbol"
        table.to_csv(attempt / "pseudobulk_counts.csv")
        cell_counts = (retained.obs.groupby(GROUPING, observed=True).size().rename("cell_count").reset_index())
        cell_counts = original_metadata.merge(cell_counts, on=GROUPING, how="left")
        cell_counts["cell_count"] = cell_counts["cell_count"].fillna(0).astype(int)
        cell_counts["sample_total"] = cell_counts.groupby("sample_accession")["cell_count"].transform("sum")
        cell_counts["cell_fraction"] = np.where(cell_counts["sample_total"] > 0, cell_counts["cell_count"] / cell_counts["sample_total"], np.nan)
        cell_counts.to_csv(attempt / "celltype_counts_by_donor.tsv", sep="\t", index=False)
        original = pd.read_csv(ORIGINAL_COUNTS, index_col=0).reindex(index=table.index, columns=table.columns)
        dense = table.sparse.to_dense().astype(np.int64)
        difference = dense.to_numpy(dtype=np.int64) - original.to_numpy(dtype=np.int64)
        comparison = {
            "shape": [int(dense.shape[0]), int(dense.shape[1])],
            "value_mismatches": int(np.count_nonzero(difference)),
            "max_absolute_difference": int(np.abs(difference).max()),
            "original_total_umi": int(original.to_numpy(dtype=np.int64).sum()),
            "doublet_excluded_total_umi": int(dense.to_numpy(dtype=np.int64).sum()),
            "retained_cells": int(retained.n_obs),
            "excluded_predicted_doublets": int(predicted.sum()),
            "unavailable_pseudobulk_groups": int((cell_counts["cell_count"] == 0).sum()),
        }
        write_json(attempt / "pseudobulk_comparison_to_03A.json", comparison)
        provenance = {
            "status": "PSEUDOBULK_REAGGREGATED_AWAITING_EDGER_GSEA",
            "started_at": started,
            "completed_at": now(),
            "inputs": {str(path.resolve()): sha256(path) for path in [INPUT, CALLS, CALLING_PROVENANCE, ORIGINAL_COUNTS, ORIGINAL_SAMPLES, PROTOCOL]},
            "resources_at_start_gb": {"available_memory": memory, "available_disk": disk},
            "grouping": GROUPING,
            "fixed_annotation": "cell_type and leiden labels from Gate 03A attempt_07; no reclustering or relabelling",
            "comparison": comparison,
            "evidence_boundary": "This reaggregation isolates the effect of excluding result-known predicted doublets. It is a sensitivity analysis, not independent validation.",
        }
        write_json(attempt / "doublet_exclusion_pseudobulk_provenance.json", provenance)
        print(json.dumps(provenance, ensure_ascii=False, indent=2))
    except Exception as exc:
        status = "HOLD" if str(exc).startswith("HOLD:") else "FAIL"
        write_json(attempt / "attempt_status.json", {"status": status, "started_at": started, "stopped_at": now(), "error": repr(exc), "traceback": traceback.format_exc()})
        raise SystemExit(2)


if __name__ == "__main__":
    main()
