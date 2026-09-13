#!/usr/bin/env python3
"""Freeze corrected ES01 attempt before the retry reads seed-result tables."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP_ROOT = PROJECT_ROOT / "revision_v2" / "07_virtual_perturbation"
OLD_ATTEMPT = VP_ROOT / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_01"
ATTEMPT = VP_ROOT / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_02"
PREFLIGHT = ATTEMPT / "00_preflight.json"
OUTPUT = ATTEMPT / "01_scope_freeze.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    if Path.cwd().resolve() != PROJECT_ROOT.resolve():
        raise SystemExit(f"Refusing to run outside frozen project root: {Path.cwd()}")
    if not PREFLIGHT.is_file() or OUTPUT.exists():
        raise SystemExit("Corrected preflight is missing or corrected scope freeze already exists")
    preflight = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    if preflight.get("status") != "PASS_INPUTS_READ_ONLY" or preflight.get("errors"):
        raise SystemExit("Corrected preflight is not a clean PASS_INPUTS_READ_ONLY record")
    template_path = OLD_ATTEMPT / "01_scope_freeze.json"
    retained_hold = OLD_ATTEMPT / "02_seed_rank_reconstruction_FAILURE_RETAINED.json"
    if not template_path.is_file() or not retained_hold.is_file():
        raise SystemExit("Retained attempt_01 scope or HOLD record is absent")
    scope = json.loads(template_path.read_text(encoding="utf-8"))
    run_script = VP_ROOT / "scripts" / "09d_rebuild_seed_rank_qc_es01_attempt02.py"
    preflight_script = VP_ROOT / "scripts" / "09c_preflight_es01_attempt02.py"
    engine_preflight = VP_ROOT / "scripts" / "09a_preflight_exploratory_hypoxia_sensitivity.py"
    engine_rebuild = VP_ROOT / "scripts" / "09b_rebuild_seed_rank_qc.py"
    freeze_script = Path(__file__).resolve()
    scope.update({
        "gate": "VP-ES01-01",
        "stage": "corrected_exploratory_analysis_scope_freeze",
        "status": "FROZEN_EXPLORATORY_SPECIFICATION",
        "created_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "attempt": "attempt_20260910_02",
        "preflight_record": {"path": str(PREFLIGHT.relative_to(VP_ROOT)).replace("\\", "/"), "sha256": sha256(PREFLIGHT), "status": preflight["status"]},
        "correction_lineage": {
            "retained_attempt": "attempt_20260910_01",
            "retained_hold_path": str(retained_hold.relative_to(VP_ROOT)).replace("\\", "/"),
            "retained_hold_sha256": sha256(retained_hold),
            "correction": "Attempt 01 input_manifest omitted all_virtual_KO_rankings.tsv.gz. Attempt 02 records every declared high-level input plus all 45 main and 90 frozen control differential-regulation files before retrying ES01-02.",
            "scope_template_sha256": sha256(template_path),
        },
        "implementation": {
            "scope_freeze_script": str(freeze_script.relative_to(VP_ROOT)).replace("\\", "/"),
            "scope_freeze_script_sha256": sha256(freeze_script),
            "preflight_script": "scripts/09c_preflight_es01_attempt02.py",
            "preflight_script_sha256": sha256(preflight_script),
            "preflight_engine_script": "scripts/09a_preflight_exploratory_hypoxia_sensitivity.py",
            "preflight_engine_script_sha256": sha256(engine_preflight),
            "es01_02_runner_script": "scripts/09d_rebuild_seed_rank_qc_es01_attempt02.py",
            "es01_02_runner_script_sha256": sha256(run_script),
            "es01_02_engine_script": "scripts/09b_rebuild_seed_rank_qc.py",
            "es01_02_engine_script_sha256": sha256(engine_rebuild),
            "result_tables_read_before_freeze": False,
            "scope_change_rule": "This corrected scope is immutable. Any subsequent correction requires a distinct append-only attempt.",
        },
    })
    OUTPUT.write_text(json.dumps(scope, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
