#!/usr/bin/env python3
"""Audit VP-G08 v3 endpoint/coculture amendment without lifting operational holds."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
V2 = VP / "01_protocol" / "VP_G08_preimplementation_amendment_v2.json"
V3 = VP / "01_protocol" / "VP_G08_endpoint_coculture_amendment_v3.json"
AUDIT = VP / "validation" / "VP_G08_endpoint_coculture_AUDIT01.json"
EXPECTED_V2 = "af44f9f6a00216a94a6eb6f58938185afbff1bbba65765819ea5788c8ebc6efd"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if AUDIT.exists():
        raise SystemExit("Refusing to overwrite VP-G08 endpoint/coculture audit")
    amendment = json.loads(V3.read_text(encoding="utf-8"))
    endpoints = {item["gene"]: item for item in amendment["single_functional_primary_endpoints"]}
    coculture = amendment["EFNB2_coculture_specification"]
    checks = [
        {"name": "retained_v2_hash_matches", "pass": V2.is_file() and sha256(V2) == EXPECTED_V2 and amendment["retains_prior_records"].get("v2_amendment_sha256") == EXPECTED_V2},
        {"name": "candidate_cell_and_TYMS_sirna_recorded", "pass": amendment["cell_source_candidate"]["fibroblast"].get("catalogue_number_as_listed") == "QuiCell-Pri-8089" and amendment["perturbation_choice"].get("TYMS", "").startswith("siRNA")},
        {"name": "three_single_functional_endpoints_have_scales", "pass": set(endpoints) == {"TYMS", "EFNB2", "LRRC17"} and "percent" in endpoints["TYMS"]["operational_measure"] and "micrometres per well" in endpoints["EFNB2"]["operational_measure"] and "final gel area" in endpoints["LRRC17"]["operational_measure"]},
        {"name": "EFNB2_modes_share_one_endpoint_and_define_contact_boundary", "pass": coculture.get("same_functional_endpoint_in_both_modes") == "CD31-positive endothelial network total length per well" and "0.4-micrometre" in coculture.get("Transwell_mode", "") and "interaction" in coculture.get("interpretation_rule", "")},
        {"name": "twelve_well_rule_preserves_technical_replication_boundary", "pass": "2 oxygen conditions x 2 intervention conditions x 3 technical wells" in coculture.get("plate_capacity_rule", "")},
        {"name": "operational_holds_and_no_execution_authorization_retained", "pass": amendment.get("status") == "PILOT_ENDPOINTS_FROZEN_REMAINING_OPERATIONAL_HOLDS" and "does not authorize purchase" in amendment.get("authorization_boundary", "") and len(amendment.get("remaining_holds_before_any_execution", [])) >= 6 and not (VP / "validation" / "VP_G08_experimental_handoff_COMPLETE.json").exists()}
    ]
    passed = all(item["pass"] for item in checks)
    payload = {
        "gate": "VP-G08",
        "audit_id": "ENDPOINT_COCULTURE_AUDIT01",
        "status": "PASS_AMENDMENT_BOUNDARY_VALIDATED" if passed else "HOLD_AUDIT_FAILED",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "v2_amendment_sha256": sha256(V2),
        "v3_amendment_sha256": sha256(V3),
        "interpretation": "PASS verifies append-only endpoint/coculture documentation and retained holds. It is not approval to purchase or execute, and not VP-G08 completion."
    }
    AUDIT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "audit": str(AUDIT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
