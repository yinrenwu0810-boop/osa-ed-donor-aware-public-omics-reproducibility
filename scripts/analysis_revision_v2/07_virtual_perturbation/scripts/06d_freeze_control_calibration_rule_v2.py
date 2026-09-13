#!/usr/bin/env python3
"""Freeze VP-G05 v2 non-pathway control calibration before calculation."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "01_protocol" / "VP_G05_control_calibration_rule_freeze_v2.json"
CONTROL_AUDIT = ROOT / "validation" / "VP_G05_control_runs_v2_AUDIT01.json"
PRECONTROL = ROOT / "validation" / "VP_G05_stability_PRECONTROL.json"
SELECTION_AUDIT = ROOT / "validation" / "VP_G05_negative_control_selection_v2_AUDIT01.json"
V1_HOLD = ROOT / "validation" / "VP_G05_negative_control_selection_HOLD_INSUFFICIENT_MATCHED_CONTROLS_ATTEMPT01.json"
IMPLEMENTATIONS = [
    "scripts/06b_calibrate_controls_v2.py",
    "scripts/06c_validate_control_calibration_v2.py",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite sealed freeze: {OUT}")
    required = {
        CONTROL_AUDIT: "PASS",
        PRECONTROL: "PASS_STABILITY_PRECONTROL",
        SELECTION_AUDIT: "PASS",
        V1_HOLD: "HOLD_INSUFFICIENT_MATCHED_CONTROLS",
    }
    for path, expected in required.items():
        status = json.loads(path.read_text(encoding="utf-8")).get("status")
        if status != expected:
            raise SystemExit(f"Unexpected prerequisite status {path.name}: {status}")
    freeze = {
        "gate": "VP-G05",
        "status": "FROZEN_BEFORE_CALIBRATION",
        "version": "v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User confirmed the v2 relaxed matched-control rule and instructed continuation through the unfinished VP-G05 work.",
        "primary_seed": 2026082701,
        "independent_unit": "donor",
        "target_run_selection": "For each target and donor, use only the PASS_TECHNICAL direct main KO at primary seed 2026082701.",
        "control_comparison": "For each target and paired donor, compare with its 10 v2 selected matched-reference control KOs at the same condition, donor, and seed.",
        "distance_universe": "All 1388 genes in each differential_regulation.tsv; include the perturbed gene; require finite nonnegative distance.",
        "top_fraction": 0.05,
        "top_n_rule": "ceil(0.05 * n_genes), giving 70 of 1388 genes; rank distance descending and break exact ties by gene symbol ascending.",
        "primary_metric": "top5pct_distance_sum",
        "secondary_metrics": ["top5pct_distance_mean", "distance_sum_all", "distance_mean_all", "distance_median", "distance_q95", "distance_q99"],
        "aggregation": "Compute each run metric first. For every target or control gene, aggregate its three paired donors by the median; donor-level target/control-median ratios are retained descriptively.",
        "empirical_percentile": "count(control donor-medians <= target donor-median) / 10",
        "conservative_upper_tail_empirical_value": "(1 + count(control donor-medians >= target donor-median)) / 11; minimum attainable value is 1/11.",
        "ties": "Inclusive comparisons are used in both tails; no randomized tie breaking.",
        "missing_values": "No missing/nonfinite distance is allowed; a violation fails the calculation rather than being imputed.",
        "inference_boundary": "All calibration quantities are descriptive prioritization aids, not causal estimates, biologically inert-control proofs, or confirmatory frequentist p-values. No multiple-testing significance claim is permitted.",
        "pathway_calibration": "DEFERRED_TO_VP_G06; do not inspect gene-set or pathway results while closing VP-G05.",
        "v1_history": "The original strict rule and its design-feasibility HOLD remain authoritative retained history and must be reported alongside v2.",
        "control_run_audit_sha256": sha256(CONTROL_AUDIT),
        "precontrol_sha256": sha256(PRECONTROL),
        "selection_v2_audit_sha256": sha256(SELECTION_AUDIT),
        "v1_hold_sha256": sha256(V1_HOLD),
        "implementation_sha256": {script: sha256(ROOT / script) for script in IMPLEMENTATIONS},
        "expected_outputs": [
            "06_consensus/negative_control_calibration_v2/run_metric_summaries.tsv",
            "06_consensus/negative_control_calibration_v2/donor_calibration.tsv",
            "06_consensus/negative_control_calibration_v2/negative_control_calibration.tsv",
            "06_consensus/negative_control_calibration_v2/calibration_report.json",
            "06_consensus/negative_control_calibration_v2/artifact_manifest.sha256.tsv",
            "validation/VP_G05_control_calibration_v2_AUDIT01.json",
        ],
    }
    OUT.write_text(json.dumps(freeze, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(freeze, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
