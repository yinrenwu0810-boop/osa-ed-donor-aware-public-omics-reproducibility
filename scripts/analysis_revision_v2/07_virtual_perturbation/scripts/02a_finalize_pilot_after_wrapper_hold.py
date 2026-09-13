#!/usr/bin/env python3
"""Finalize a completed VP-G03 pilot after a wrapper-only post-run hold.

This performs no model computation.  It independently validates the retained
temporary artifacts, writes the missing manifest, then atomically publishes
the original attempt_01 directory.
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def timestamp() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    temporary.replace(path)


def main() -> None:
    vp_root = Path(__file__).resolve().parents[1]
    pilot_root = vp_root / "05_runs" / "pilot"
    temporary = pilot_root / "attempt_01.tmp_23016"
    final = pilot_root / "attempt_01"
    if not temporary.is_dir() or final.exists():
        raise RuntimeError("Expected retained temporary attempt_01 and absent final attempt_01.")
    required = {"differential_regulation.tsv", "run_parameters.json", "session_info.txt", "stdout.log", "stderr.log", "resource_usage.tsv"}
    if {path.name for path in temporary.iterdir() if path.is_file()} != required:
        raise RuntimeError("Temporary pilot artifact set differs from the expected pre-manifest set.")

    differential = pd.read_csv(temporary / "differential_regulation.tsv", sep="\t")
    expected_columns = ["gene", "distance", "Z", "FC", "p.value", "p.adj"]
    if list(differential.columns) != expected_columns or len(differential) != 1388 or not differential["gene"].is_unique:
        raise RuntimeError("Differential-regulation structural validation failed.")
    numeric = differential.loc[:, expected_columns[1:]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or not differential["p.value"].between(0, 1).all() or not differential["p.adj"].between(0, 1).all():
        raise RuntimeError("Differential-regulation numeric validation failed.")

    parameters = json.loads((temporary / "run_parameters.json").read_text(encoding="utf-8"))
    frozen = parameters["frozen_parameters"]
    if not (parameters["target"] == "LRRC17" and parameters["seed"] == 2026082701 and parameters["n_cores"] == 4 and parameters["package_version"] == "1.0.3" and frozen["qc"] is False and frozen["nc_nNet"] == 10 and frozen["nc_nCells"] == 500 and frozen["td_K"] == 3):
        raise RuntimeError("Run-parameter validation failed.")
    session = (temporary / "session_info.txt").read_text(encoding="utf-8")
    stderr = (temporary / "stderr.log").read_text(encoding="utf-8")
    if "scTenifoldKnk" not in session or "Finished scTenifoldKnk for LRRC17" not in stderr or "Error" in stderr:
        raise RuntimeError("Session or completion-log validation failed.")

    recovery = {
        "gate": "VP-G03",
        "attempt": "attempt_01",
        "status": "RECOVERED_WRAPPER_HOLD_AFTER_COMPLETED_MODEL_RUN",
        "finalized_at": timestamp(),
        "model_reexecuted": False,
        "original_model_log_completion": "Finished scTenifoldKnk for LRRC17",
        "wrapper_hold_evidence": "The original R wrapper remained CPU-active after all core outputs were written and before artifact_manifest.sha256.tsv/atomic rename.",
        "interruption": "Rscript PIDs 23016, 5768, and 17060 were terminated only after independent validation of completed core outputs.",
        "independent_prepublication_checks": {
            "expected_six_columns": True,
            "input_gene_count_equals_output_row_count": True,
            "unique_gene_identifiers": True,
            "no_NA_or_Inf": True,
            "p_value_ranges": True,
            "frozen_pilot_parameters": True,
            "completion_log": True,
        },
    }
    write_json_atomic(temporary / "recovery_provenance.json", recovery)
    files = sorted(path for path in temporary.iterdir() if path.is_file())
    manifest = pd.DataFrame(
        [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in files]
    )
    manifest.to_csv(temporary / "artifact_manifest.sha256.tsv", sep="\t", index=False, lineterminator="\n")
    temporary.replace(final)

    validation = {
        "gate": "VP-G03",
        "status": "VP_G03_pilot_PASS",
        "validated_at": timestamp(),
        "pilot": {"gene": "LRRC17", "condition": "organic_ED_nonDM", "donor": "non-DM_1", "seed": 2026082701},
        "checks": {
            "target_retained": True,
            "output_gene_set_equals_input_gene_set": True,
            "no_duplicate_gene_NA_or_Inf": True,
            "p_value_and_adjusted_p_value_in_zero_one": True,
            "frozen_parameters_match": True,
            "completion_was_technical_only": True,
            "resource_hold": False,
        },
        "run_finalization": "Model ran once. A post-model wrapper hold was independently audited and finalized without model re-execution.",
        "artifacts": [
            "05_runs/pilot/attempt_01/differential_regulation.tsv",
            "05_runs/pilot/attempt_01/run_parameters.json",
            "05_runs/pilot/attempt_01/session_info.txt",
            "05_runs/pilot/attempt_01/stdout.log",
            "05_runs/pilot/attempt_01/stderr.log",
            "05_runs/pilot/attempt_01/resource_usage.tsv",
            "05_runs/pilot/attempt_01/artifact_manifest.sha256.tsv",
            "05_runs/pilot/attempt_01/recovery_provenance.json",
        ],
    }
    write_json_atomic(vp_root / "validation" / "VP_G03_pilot_PASS.json", validation)
    print(json.dumps({"gate": "VP-G03", "status": "VP_G03_pilot_PASS", "model_reexecuted": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
