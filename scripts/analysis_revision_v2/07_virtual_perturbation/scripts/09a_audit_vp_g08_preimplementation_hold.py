#!/usr/bin/env python3
"""Audit the VP-G08 preimplementation HOLD without converting it into completion."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
FREEZE = VP / "01_protocol" / "VP_G08_experimental_handoff_scope_freeze_v1.json"
G07 = VP / "validation" / "VP_G07_external_scope_COMPLETE.json"
HOLD = VP / "validation" / "VP_G08_preimplementation_HOLD.json"
AUDIT = VP / "validation" / "VP_G08_preimplementation_AUDIT01.json"
EXPECTED_G07 = "630a638be38c01ff36cd6448d5b4fb5e3697422e51f9cabd75aa193f96698255"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if AUDIT.exists():
        raise SystemExit("Refusing to overwrite VP-G08 preimplementation audit")
    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    hold = json.loads(HOLD.read_text(encoding="utf-8"))
    checks = [
        {"name":"g07_completion_hash_matches", "pass":G07.is_file() and sha256(G07) == EXPECTED_G07},
        {"name":"freeze_is_preimplementation_hold", "pass":freeze.get("status") == "FROZEN_PREIMPLEMENTATION_HOLD"},
        {"name":"two_by_two_design_and_three_targets_recorded", "pass":"2x2" in freeze.get("fixed_design", {}).get("core_factorial_structure", "") and len(freeze.get("target_direction_mapping", [])) == 3},
        {"name":"independent_unit_and_blocking_recorded", "pass":"donor" in freeze.get("fixed_design", {}).get("unit_of_biological_inference", "").lower() and len(freeze.get("fixed_design", {}).get("blocked_factors", [])) >= 4},
        {"name":"operational_holds_are_explicit", "pass":hold.get("status") == "HOLD_OPERATIONAL_PARAMETERS_NOT_FROZEN" and len(hold.get("blocking_requirements", [])) >= 5},
        {"name":"no_wet_lab_or_completion_claim", "pass":"wet-lab execution" in hold.get("not_authorized", []) and not (VP / "validation" / "VP_G08_experimental_handoff_COMPLETE.json").exists()}
    ]
    passed = all(item["pass"] for item in checks)
    payload = {"gate":"VP-G08", "audit_id":"PREIMPLEMENTATION_AUDIT01", "status":"PASS_HOLD_STATE_VALIDATED" if passed else "HOLD_AUDIT_FAILED", "generated_at_utc":datetime.now(timezone.utc).isoformat(), "checks":checks, "freeze_sha256":sha256(FREEZE), "hold_sha256":sha256(HOLD), "interpretation":"A passing audit validates that the HOLD is explicit and scientifically bounded; it is not VP-G08 completion and does not authorize experiments."}
    AUDIT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status":payload["status"], "audit":str(AUDIT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
