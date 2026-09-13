"""Freeze the user-approved VP-G05 relaxed control-matching amendment."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
OUT = VP / "01_protocol" / "VP_G05_negative_control_rule_freeze_v2.json"
V1 = VP / "01_protocol" / "VP_G05_negative_control_rule_freeze_v1.json"
V1_HOLD = VP / "validation" / "VP_G05_negative_control_selection_HOLD_INSUFFICIENT_MATCHED_CONTROLS_ATTEMPT01.json"
DEGREE_AUDIT = VP / "validation" / "VP_G05_degree_reference_AUDIT01.json"
IMPLEMENTATIONS = [
    VP / "scripts" / "05b_select_matched_negative_controls_v2.py",
    VP / "scripts" / "05c_validate_negative_control_selection_v2.py",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite frozen v2 rule: {OUT}")
    required = [V1, V1_HOLD, DEGREE_AUDIT, *IMPLEMENTATIONS]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"Required retained evidence missing: {missing}")
    if json.loads(DEGREE_AUDIT.read_text(encoding="utf-8")).get("status") != "PASS":
        raise RuntimeError("Degree-reference independent audit is not PASS.")
    rule = {
        "gate": "VP-G05",
        "status": "FROZEN_AMENDMENT_BEFORE_V2_SELECTION",
        "version": "v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User explicitly confirmed the recommended v2 rule in this conversation.",
        "amendment_reason": "The retained v1 rule produced 1/0/0 eligible controls for TYMS/EFNB2/LRRC17 and entered HOLD. V2 is a transparent feasibility amendment, not the original prespecified rule.",
        "retained_v1": {"rule_sha256": sha256(V1), "hold_sha256": sha256(V1_HOLD)},
        "degree_reference_audit_sha256": sha256(DEGREE_AUDIT),
        "target_definitions": [
            {"target_gene": "TYMS", "condition": "normal", "donors": ["Normal_1", "Normal_2", "Normal_3"]},
            {"target_gene": "EFNB2", "condition": "organic_ED_nonDM", "donors": ["non-DM_1", "non-DM_2", "non-DM_3"]},
            {"target_gene": "LRRC17", "condition": "organic_ED_nonDM", "donors": ["non-DM_1", "non-DM_2", "non-DM_3"]},
        ],
        "matching_rule": "In every one of the three paired donors, abs(control detection decile - target detection decile) <= 1 and abs(control outdegree decile - target outdegree decile) <= 1.",
        "all_donors_required": True,
        "detection_decile_tolerance": 1,
        "outdegree_decile_tolerance": 1,
        "degree_and_detection_definitions": "Unchanged from v1; use the independently audited frozen degree_reference.tsv.",
        "candidate_pool": "Genes present in all paired donor matrices before exclusions.",
        "exclusions": ["TYMS", "EFNB2", "LRRC17", "All genes marked HALLMARK_HYPOXIA_DETECTABLE_IN_MAIN_FIBROBLASTS in the frozen gene universe"],
        "controls_per_target": 10,
        "sampling": {"seed": 2026082799, "algorithm": "Python random.Random initialized once; target definitions in listed order; sample without replacement within each gene-symbol-sorted eligible pool."},
        "cross_target_overlap": "Allowed. A shared control remains target-specific analytically; identical condition/donor/seed KO tasks are run once and cross-referenced.",
        "initial_control_KO_seed": 2026082701,
        "interpretation": "V2 controls are relaxed matched reference genes, not biologically inert negative genes. V1 HOLD must be reported with V2 calibration.",
        "prohibited": ["additional result-driven exclusions", "pathway inspection before selection is sealed", "rule widening after selection", "overwriting v1 or v2 evidence"],
        "implementation_sha256": {str(path.relative_to(VP)).replace("\\", "/"): sha256(path) for path in IMPLEMENTATIONS},
    }
    OUT.write_text(json.dumps(rule, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": rule["status"], "path": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
