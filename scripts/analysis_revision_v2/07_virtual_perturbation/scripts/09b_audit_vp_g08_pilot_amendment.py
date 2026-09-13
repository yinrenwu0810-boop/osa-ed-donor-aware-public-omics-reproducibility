#!/usr/bin/env python3
"""Audit the append-only VP-G08 feasibility-pilot amendment without claiming completion."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
V1 = VP / "01_protocol" / "VP_G08_experimental_handoff_scope_freeze_v1.json"
V2 = VP / "01_protocol" / "VP_G08_preimplementation_amendment_v2.json"
G07 = VP / "validation" / "VP_G07_external_scope_COMPLETE.json"
AUDIT = VP / "validation" / "VP_G08_pilot_amendment_AUDIT01.json"
EXPECTED_V1 = "bb77df1b70176b50e0902f32b97d7021b6203b98af8255e226e684d5c1637257"
EXPECTED_G07 = "630a638be38c01ff36cd6448d5b4fb5e3697422e51f9cabd75aa193f96698255"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if AUDIT.exists():
        raise SystemExit("Refusing to overwrite VP-G08 pilot-amendment audit")
    amendment = json.loads(V2.read_text(encoding="utf-8"))
    ih = amendment["newly_frozen_parameters"]["IH_cycle_programme_requested"]
    pilot = amendment["pilot_design"]
    checks = [
        {"name": "retained_v1_freeze_hash_matches", "pass": V1.is_file() and sha256(V1) == EXPECTED_V1 and amendment.get("v1_freeze_sha256") == EXPECTED_V1},
        {"name": "g07_completion_hash_matches", "pass": G07.is_file() and sha256(G07) == EXPECTED_G07},
        {"name": "ih_programme_recorded", "pass": amendment["newly_frozen_parameters"]["normoxia"].get("oxygen_percent") == 21 and ih.get("descent_to_1_percent_O2_seconds_approx") == 30 and ih.get("hypoxia_1_percent_O2_hold_minutes") == 5 and ih.get("reoxygenation_to_21_percent_O2_seconds_approx") == 30 and ih.get("normoxia_21_percent_O2_hold_minutes") == 5 and ih.get("observation_timepoints_hours") == [0, 8, 16, 24]},
        {"name": "twelve_well_capacity_not_misrepresented", "pass": "4 conditions x 3 technical wells" in pilot.get("plate_unit", "") and amendment["newly_frozen_parameters"]["pilot_capacity"].get("available_wells_per_batch") == 12},
        {"name": "cell_source_boundary_is_explicit", "pass": "No primary human corpus-cavernosum fibroblasts" in amendment["cell_source_decision_boundary"].get("current_state", "") and "non-T1" in amendment["cell_source_decision_boundary"].get("non_CCFB_commercial_fibroblast_path", "")},
        {"name": "pilot_not_completion_or_wet_lab_authorization", "pass": amendment.get("status") == "PILOT_SCOPE_PARTIALLY_FROZEN_REMAINING_HOLDS" and "wet-lab execution" in amendment.get("authorization_boundary", "") and not (VP / "validation" / "VP_G08_experimental_handoff_COMPLETE.json").exists()},
        {"name": "unresolved_endpoints_and_operational_holds_retained", "pass": len(pilot["pilot_endpoint_hierarchy"].get("unresolved_gene_specific_endpoints", [])) == 3 and len(amendment.get("unresolved_holds", [])) >= 6},
    ]
    passed = all(item["pass"] for item in checks)
    payload = {
        "gate": "VP-G08",
        "audit_id": "PILOT_AMENDMENT_AUDIT01",
        "status": "PASS_AMENDMENT_BOUNDARY_VALIDATED" if passed else "HOLD_AUDIT_FAILED",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "v1_freeze_sha256": sha256(V1),
        "v2_amendment_sha256": sha256(V2),
        "interpretation": "PASS validates only that supplied feasibility parameters were appended without relaxing the operational HOLD. It is not VP-G08 completion and does not authorize experiments or purchase."
    }
    AUDIT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "audit": str(AUDIT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
