#!/usr/bin/env python3
"""Append-only VP-ES01-02 retry wrapper for attempt_20260910_02."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


PROJECT_ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP_ROOT = PROJECT_ROOT / "revision_v2" / "07_virtual_perturbation"
ATTEMPT = VP_ROOT / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_02"
ENGINE = VP_ROOT / "scripts" / "09b_rebuild_seed_rank_qc.py"


def main() -> int:
    spec = importlib.util.spec_from_file_location("es01_rank_attempt01_engine", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load retained ES01-02 engine: {ENGINE}")
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    engine.ATTEMPT = ATTEMPT
    engine.PREFLIGHT = ATTEMPT / "00_preflight.json"
    engine.SCOPE = ATTEMPT / "01_scope_freeze.json"
    engine.OUTPUTS = [
        ATTEMPT / "02_seed_rank_qc.tsv",
        ATTEMPT / "02_seed_pair_stability.tsv",
        ATTEMPT / "02_seed_aggregate_reconciliation.tsv",
    ]
    engine.FAILURE = ATTEMPT / "02_seed_rank_reconstruction_FAILURE_RETAINED.json"
    return engine.main()


if __name__ == "__main__":
    sys.exit(main())
