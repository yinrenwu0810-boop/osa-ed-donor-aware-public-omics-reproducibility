#!/usr/bin/env python3
"""Publish a completed G04 run after a post-model empty-log hashing failure.

No scTenifoldKnk call is made by this recovery tool.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


EMPTY_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--temporary-attempt", required=True, type=Path)
    parser.add_argument("--final-attempt", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    args = parser.parse_args()
    temporary = args.temporary_attempt.resolve()
    final = args.final_attempt.resolve()
    registry = args.registry.resolve()
    if not temporary.is_dir() or final.exists() or registry.exists():
        raise RuntimeError("Recovery requires a retained temporary attempt, absent final attempt, and a fresh registry path.")
    required = {"differential_regulation.tsv", "run_parameters.json", "session_info.txt", "stdout.log", "stderr.log", "resource_usage.tsv"}
    if {p.name for p in temporary.iterdir() if p.is_file()} != required:
        raise RuntimeError("Unexpected temporary artifact set.")
    d = pd.read_csv(temporary / "differential_regulation.tsv", sep="\t")
    expected = ["gene", "distance", "Z", "FC", "p.value", "p.adj"]
    if list(d.columns) != expected or d.shape != (1388, 6) or not d["gene"].is_unique:
        raise RuntimeError("Differential-regulation structural check failed.")
    if not np.isfinite(d.iloc[:, 1:].to_numpy(dtype=float)).all() or not d["p.value"].between(0, 1).all() or not d["p.adj"].between(0, 1).all():
        raise RuntimeError("Differential-regulation numeric check failed.")
    p = json.loads((temporary / "run_parameters.json").read_text(encoding="utf-8"))
    if not (p["target"] == "TYMS" and p["seed"] == 2026082701 and p["n_cores"] == 4 and p["package_version"] == "1.0.3" and p["frozen_parameters"]["qc"] is False):
        raise RuntimeError("Frozen run-parameter check failed.")
    stderr = (temporary / "stderr.log").read_text(encoding="utf-8")
    if "Finished scTenifoldKnk for TYMS" not in stderr or "Error" in stderr:
        raise RuntimeError("Completion-log check failed.")
    if sha256(temporary / "stdout.log") != EMPTY_SHA256:
        raise RuntimeError("Expected empty stdout log hash is not standard SHA-256.")
    recovery = {
        "gate": "VP-G04",
        "run_id": "VP_MAIN_01",
        "status": "RECOVERED_POSTMODEL_HASHING_FAILURE",
        "finalized_at": now(),
        "model_reexecuted": False,
        "failure_cause": "certutil returns ERROR_FILE_INVALID for zero-byte stdout.log on this Windows host.",
        "repair": "The runner now assigns the standard SHA-256 empty-byte-stream digest before invoking certutil.",
        "technical_checks": {"parameters": True, "output_shape": True, "unique_genes": True, "finite_values": True, "p_value_ranges": True, "completion_log": True},
    }
    write_json(temporary / "recovery_provenance.json", recovery)
    files = sorted(p for p in temporary.iterdir() if p.is_file())
    pd.DataFrame([{"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)} for p in files]).to_csv(temporary / "artifact_manifest.sha256.tsv", sep="\t", index=False, lineterminator="\n")
    final.parent.mkdir(parents=True, exist_ok=True)
    temporary.replace(final)
    registry.parent.mkdir(parents=True, exist_ok=True)
    fields = ["run_order", "run_id", "gene", "condition", "donor", "seed", "started_at_utc", "finished_at_utc", "status", "exit_code", "output_attempt", "note"]
    with registry.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow({"run_order": 1, "run_id": "VP_MAIN_01", "gene": "TYMS", "condition": "normal", "donor": "Normal_1", "seed": 2026082701, "started_at_utc": "2026-08-27T21:31:07.9429193Z", "finished_at_utc": now(), "status": "PASS_TECHNICAL", "exit_code": 0, "output_attempt": str(final), "note": "Model ran once; independently recovered after a zero-byte stdout hashing failure. No biological interpretation performed."})
    print(json.dumps({"run_id": "VP_MAIN_01", "status": "PASS_TECHNICAL", "model_reexecuted": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
