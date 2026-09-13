#!/usr/bin/env python3
"""Append-only corrected VP-ES01-00 preflight for attempt_20260910_02."""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


PROJECT_ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP_ROOT = PROJECT_ROOT / "revision_v2" / "07_virtual_perturbation"
ATTEMPT = VP_ROOT / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_02"
ENGINE = VP_ROOT / "scripts" / "09a_preflight_exploratory_hypoxia_sensitivity.py"


def load_engine():
    spec = importlib.util.spec_from_file_location("es01_preflight_attempt01_engine", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load retained preflight engine: {ENGINE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rel(path: Path) -> str:
    return str(path.relative_to(VP_ROOT)).replace("\\", "/")


def main() -> int:
    if Path.cwd().resolve() != PROJECT_ROOT.resolve():
        raise SystemExit(f"Refusing to run outside frozen project root: {Path.cwd()}")
    if ATTEMPT.exists():
        raise SystemExit(f"Refusing to overwrite existing corrected attempt directory: {ATTEMPT}")
    engine = load_engine()
    registry = engine.read_registry(VP_ROOT / "05_runs" / "run_registry.tsv")
    main_outputs = [rel(Path(row["output_attempt"]) / "differential_regulation.tsv") for row in registry]
    control_root = VP_ROOT / "05_runs" / "control_batches_v2"
    control_outputs = [rel(path) for path in sorted(control_root.glob("*/*/*/attempt_01/controls/*/attempt_01/differential_regulation.tsv"))]
    calibration_files = [rel(path) for path in sorted((VP_ROOT / "06_consensus" / "negative_control_calibration_v2").rglob("*")) if path.is_file()]
    added_declared_inputs = [
        "06_consensus/all_virtual_KO_rankings.tsv.gz",
        "06_consensus/candidate_donor_consensus.tsv",
        "06_consensus/donor_stability.tsv",
    ]
    manifest_files = list(dict.fromkeys(engine.MANIFEST_FILES + added_declared_inputs + calibration_files + main_outputs + control_outputs))
    if len(main_outputs) != 45 or len(control_outputs) != 90:
        raise SystemExit(f"expected 45 main and 90 control outputs before correction; found {len(main_outputs)} and {len(control_outputs)}")
    engine.ATTEMPT_DIR = ATTEMPT
    engine.MANIFEST_FILES = manifest_files
    result = engine.main()
    print({"gate": "VP-ES01-00", "attempt": "attempt_20260910_02", "manifest_entries": len(manifest_files), "status_code": result})
    return result


if __name__ == "__main__":
    sys.exit(main())
