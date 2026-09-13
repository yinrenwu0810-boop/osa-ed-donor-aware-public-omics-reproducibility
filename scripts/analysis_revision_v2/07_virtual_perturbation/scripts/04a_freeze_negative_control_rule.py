"""Freeze VP-G05 negative-control matching before computing eligible controls.

This script deliberately contains no virtual-KO result handling. It records
the user-approved, donor-specific detection and outdegree matching definition
and binds that definition to the implementation scripts by SHA-256.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
OUT = VP / "01_protocol" / "VP_G05_negative_control_rule_freeze_v1.json"
PRECONTROL = VP / "validation" / "VP_G05_stability_PRECONTROL.json"
IMPLEMENTATIONS = [
    VP / "scripts" / "04b_build_negative_control_degree_reference.R",
    VP / "scripts" / "04c_select_matched_negative_controls.py",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite frozen rule: {OUT}")
    report = json.loads(PRECONTROL.read_text(encoding="utf-8"))
    if report.get("status") != "PASS_STABILITY_PRECONTROL":
        raise RuntimeError("VP-G05 pre-control stability audit is not PASS.")
    missing = [str(path) for path in IMPLEMENTATIONS if not path.is_file()]
    if missing:
        raise RuntimeError(f"Implementation missing before freeze: {missing}")
    rule = {
        "gate": "VP-G05",
        "status": "FROZEN_BEFORE_NEGATIVE_CONTROL_SELECTION",
        "version": "v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User confirmed the proposed outdegree definition in this conversation.",
        "target_definitions": [
            {"target_gene": "TYMS", "condition": "normal", "donors": ["Normal_1", "Normal_2", "Normal_3"]},
            {"target_gene": "EFNB2", "condition": "organic_ED_nonDM", "donors": ["non-DM_1", "non-DM_2", "non-DM_3"]},
            {"target_gene": "LRRC17", "condition": "organic_ED_nonDM", "donors": ["non-DM_1", "non-DM_2", "non-DM_3"]},
        ],
        "outdegree_definition": {
            "network_source": "Per donor, 10 raw scTenifoldNet::makeNetworks PCNet GRNs built from its frozen VP-G02 matrix.",
            "random_seed_per_donor": 2026082701,
            "model_parameters": {"qc": False, "nNet": 10, "nCells": 500, "nComp": 3, "scaleScores": True, "symmetric": False, "q": 0.9, "nCores": 4},
            "edge_orientation": "pcNet rows are targets and columns are regulators; therefore outgoing degree is the number of non-zero entries in each regulator column.",
            "formula": "mean_{network=1..10}(count_nonzero(network[, gene]))",
            "zero_threshold": "exact numeric zero; an edge counts iff its stored PCNet value is not equal to 0.",
            "aggregation": "No cross-donor degree aggregation. Matching is evaluated separately in every donor paired to the target.",
        },
        "detection_definition": "Per donor, raw-count detection rate = number of cells with count > 0 divided by total fibroblast cells in the frozen matrix.",
        "decile_definition": "Within each donor and its full frozen matrix gene universe, rank ascending by metric then gene symbol ascending; decile = floor(10*(rank-1)/n_genes)+1. This deterministically resolves ties.",
        "matching_rule": "A candidate must have the same detection-rate decile and an outdegree decile differing by at most 1 from the target in every paired donor.",
        "candidate_pool": "Genes present in all paired donor matrices before exclusions. A pool with fewer than 10 eligible genes yields HOLD_INSUFFICIENT_MATCHED_CONTROLS; no relaxation is permitted.",
        "exclusions": ["TYMS", "EFNB2", "LRRC17", "Every gene whose frozen gene-universe selection_reason contains HALLMARK_HYPOXIA_DETECTABLE_IN_MAIN_FIBROBLASTS"],
        "controls_per_target": 10,
        "sampling": {"seed": 2026082799, "algorithm": "Python random.Random initialized once; target definitions in listed order; each target samples 10 without replacement from its gene-symbol-sorted eligible pool."},
        "prohibited": ["pathway inspection", "differential-regulation result inspection", "result-driven rule widening", "overwriting this freeze"],
        "implementation_sha256": {str(path.relative_to(VP)).replace("\\", "/"): sha256(path) for path in IMPLEMENTATIONS},
    }
    OUT.write_text(json.dumps(rule, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": rule["status"], "path": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
