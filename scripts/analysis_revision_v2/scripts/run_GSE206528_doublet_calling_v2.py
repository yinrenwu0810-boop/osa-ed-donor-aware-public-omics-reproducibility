"""Run the frozen, sample-wise Gate 03B Scrublet calling stage.

This does not replace Gate 03A and does not itself declare Gate 03 complete.
It creates an immutable attempt containing per-cell calls and required summaries.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
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
import scanpy as sc


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
INPUT = REVISION / "01_work" / "GSE206528_replay" / "attempt_07" / "GSE206528_annotated_v2.h5ad"
PROTOCOL = REVISION / "00_protocol" / "GATE_03B_DOUBLET_SENSITIVITY_PROTOCOL.md"
ATTEMPT_ROOT = REVISION / "01_work" / "GSE206528_doublet_sensitivity"
MIN_FREE_MEMORY_GB = 8.0
MIN_FREE_DISK_GB = 10.0
PARAMETERS = {
    "expected_doublet_rate": 0.05,
    "stdev_doublet_rate": 0.02,
    "sim_doublet_ratio": 2.0,
    "synthetic_doublet_umi_subsampling": 1.0,
    "n_prin_comps": 30,
    "random_state": 0,
    "threshold": None,
    "use_approx_neighbors": False,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    temporary.replace(path)


def next_attempt() -> Path:
    ATTEMPT_ROOT.mkdir(parents=True, exist_ok=True)
    existing = []
    for candidate in ATTEMPT_ROOT.glob("attempt_*"):
        try:
            existing.append(int(candidate.name.split("_")[1]))
        except (IndexError, ValueError):
            continue
    number = max(existing, default=0) + 1
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
    started = utc_now()
    try:
        if not INPUT.is_file() or not PROTOCOL.is_file():
            raise RuntimeError("Required Gate 03A input or frozen Gate 03B protocol is absent")
        memory = psutil.virtual_memory().available / 1024**3
        disk = shutil.disk_usage(attempt).free / 1024**3
        if memory < MIN_FREE_MEMORY_GB:
            raise RuntimeError(f"HOLD: available RAM {memory:.2f} GB is below {MIN_FREE_MEMORY_GB:.2f} GB")
        if disk < MIN_FREE_DISK_GB:
            raise RuntimeError(f"HOLD: available disk {disk:.2f} GB is below {MIN_FREE_DISK_GB:.2f} GB")
        adata = ad.read_h5ad(INPUT)
        if "counts" not in adata.layers:
            raise RuntimeError("Gate 03A annotated input has no counts layer")
        required_obs = {"sample_accession", "donor", "condition", "leiden", "cell_type"}
        if not required_obs.issubset(adata.obs.columns):
            raise RuntimeError(f"Required observation columns absent: {sorted(required_obs - set(adata.obs.columns))}")
        accessions = sorted(adata.obs["sample_accession"].astype(str).unique())
        if len(accessions) != 8:
            raise RuntimeError(f"Expected eight sample accessions, observed {accessions}")
        calls: list[pd.DataFrame] = []
        summaries: list[dict[str, object]] = []
        for accession in accessions:
            mask = adata.obs["sample_accession"].astype(str).eq(accession).to_numpy()
            observed = ad.AnnData(
                X=adata.layers["counts"][mask].copy(),
                obs=adata.obs.loc[mask, ["sample_accession", "donor", "condition", "leiden", "cell_type"]].copy(),
                var=adata.var.copy(),
            )
            observed.obs_names = adata.obs_names[mask].copy()
            sc.pp.scrublet(observed, verbose=False, **PARAMETERS)
            details = observed.uns.get("scrublet", {})
            predicted = observed.obs["predicted_doublet"].astype(bool)
            summary = {
                "sample_accession": accession,
                "donor": str(observed.obs["donor"].iloc[0]),
                "condition": str(observed.obs["condition"].iloc[0]),
                "input_cells": int(observed.n_obs),
                "predicted_doublets": int(predicted.sum()),
                "predicted_doublet_fraction": float(predicted.mean()),
                "expected_doublet_rate": PARAMETERS["expected_doublet_rate"],
                "automatic_threshold": float(details["threshold"]),
                "simulated_doublets": int(len(details["doublet_scores_sim"])),
            }
            summaries.append(summary)
            frame = observed.obs.copy()
            frame.insert(0, "cell_id", observed.obs_names.astype(str))
            frame.insert(1, "automatic_threshold", summary["automatic_threshold"])
            calls.append(frame.reset_index(drop=True))
            del observed
        calls_table = pd.concat(calls, ignore_index=True)
        if len(calls_table) != adata.n_obs or calls_table["cell_id"].nunique() != adata.n_obs:
            raise RuntimeError("Per-sample calling does not cover each Gate 03A cell exactly once")
        calls_table.to_csv(attempt / "doublet_calls_by_cell.tsv", sep="\t", index=False)
        sample_table = pd.DataFrame(summaries).sort_values("sample_accession")
        sample_table.to_csv(attempt / "doublet_summary_by_sample.tsv", sep="\t", index=False)
        for column, filename in [("leiden", "doublet_summary_by_historical_cluster.tsv"), ("cell_type", "doublet_summary_by_cell_type.tsv")]:
            table = (calls_table.groupby(column, observed=True)["predicted_doublet"].agg([("input_cells", "size"), ("predicted_doublets", "sum"), ("predicted_doublet_fraction", "mean")]).reset_index())
            table.to_csv(attempt / filename, sep="\t", index=False)
        provenance = {
            "status": "CALLING_COMPLETE_AWAITING_DOWNSTREAM_SENSITIVITY",
            "stage": "Gate_03B_samplewise_doublet_calling",
            "started_at": started,
            "completed_at": utc_now(),
            "input": {"path": str(INPUT.resolve()), "sha256": sha256(INPUT), "shape": [int(adata.n_obs), int(adata.n_vars)]},
            "protocol": {"path": str(PROTOCOL.resolve()), "sha256": sha256(PROTOCOL)},
            "parameters": PARAMETERS,
            "resources_at_start_gb": {"available_memory": memory, "available_disk": disk},
            "software": {name: importlib.metadata.version(name) for name in ["scanpy", "anndata", "scrublet", "scikit-image", "scikit-learn", "numpy", "scipy", "pandas"]},
            "annoy_status": "not_installed; protocol-fixed use_approx_neighbors=false uses exact scikit-learn nearest neighbors",
            "samples": int(len(sample_table)),
            "cells": int(adata.n_obs),
            "predicted_doublets": int(calls_table["predicted_doublet"].sum()),
            "evidence_boundary": "Result-known technical sensitivity analysis only. It does not establish independent biological validation or complete Gate 03 until downstream sensitivity is completed.",
        }
        atomic_json(attempt / "doublet_calling_provenance.json", provenance)
        print(json.dumps(provenance, ensure_ascii=False, indent=2))
    except Exception as exc:
        status = "HOLD" if str(exc).startswith("HOLD:") else "FAIL"
        atomic_json(attempt / "attempt_status.json", {"status": status, "started_at": started, "stopped_at": utc_now(), "error": repr(exc), "traceback": traceback.format_exc()})
        raise SystemExit(2)


if __name__ == "__main__":
    main()
