#!/usr/bin/env python3
"""Final independent acceptance audit for VP-G05."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
COMPLETE = VALIDATION / "VP_G05_stability_COMPLETE.json"
FAIL = VALIDATION / "VP_G05_stability_COMPLETE_ATTEMPT01_FAIL_RETAINED.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def main() -> None:
    if COMPLETE.exists() or FAIL.exists():
        raise SystemExit("Refusing to overwrite an existing VP-G05 final-attempt record")
    errors: list[str] = []
    prerequisites = {
        "VP_G04_main_runs_COMPLETE.json": "PASS",
        "VP_G05_stability_PRECONTROL.json": "PASS_STABILITY_PRECONTROL",
        "VP_G05_degree_reference_AUDIT01.json": "PASS",
        "VP_G05_negative_control_selection_HOLD_INSUFFICIENT_MATCHED_CONTROLS_ATTEMPT01.json": "HOLD_INSUFFICIENT_MATCHED_CONTROLS",
        "VP_G05_negative_control_selection_v2_AUDIT01.json": "PASS",
        "VP_G05_control_execution_freeze_v2_AUDIT01.json": "PASS",
        "VP_G05_control_launcher_compatibility_erratum_AUDIT01.json": "PASS",
        "VP_G05_control_runs_v2_AUDIT01.json": "PASS",
        "VP_G05_control_calibration_v2_AUDIT01.json": "PASS",
    }
    evidence: dict[str, dict[str, str]] = {}
    for name, expected in prerequisites.items():
        path = VALIDATION / name
        if not path.is_file():
            errors.append(f"Missing prerequisite: {name}")
            continue
        actual = load(path).get("status")
        evidence[name] = {"status": actual, "sha256": sha256(path)}
        if actual != expected:
            errors.append(f"Unexpected status for {name}: {actual}, expected {expected}")

    ranking_freeze_path = ROOT / "06_consensus" / "VP_G05_ranking_rule_freeze.json"
    ranking_freeze = load(ranking_freeze_path)
    if ranking_freeze.get("donor_consensus_rule") != "gene is consensus when it enters top 5% in at least 2/3 donors":
        errors.append("Frozen 2/3 donor consensus rule mismatch")
    if "pathway inspection" not in ranking_freeze.get("prohibited", []):
        errors.append("Ranking freeze does not prohibit pathway inspection")

    precontrol = load(VALIDATION / "VP_G05_stability_PRECONTROL.json")
    target_status = precontrol.get("target_status", [])
    if len(target_status) != 3 or {r.get("target_gene") for r in target_status} != {"TYMS", "EFNB2", "LRRC17"}:
        errors.append("Precontrol target membership mismatch")
    if any(r.get("donors") != 3 or r.get("stable_donors") != 3 or r.get("status") != "PASS_STABILITY" for r in target_status):
        errors.append("Not all targets passed stability in all three donors")

    rankings_path = ROOT / "06_consensus" / "all_virtual_KO_rankings.tsv.gz"
    consensus_path = ROOT / "06_consensus" / "candidate_donor_consensus.tsv"
    donor_stability_path = ROOT / "06_consensus" / "donor_stability.tsv"
    ranking_rows = rows(rankings_path)
    consensus_rows = rows(consensus_path)
    stability_rows = rows(donor_stability_path)
    if len(ranking_rows) != 3 * 3 * 5 * 1388:
        errors.append(f"Full ranking row count mismatch: {len(ranking_rows)}")
    if len(consensus_rows) != 3 * 1388:
        errors.append(f"Donor consensus row count mismatch: {len(consensus_rows)}")
    if len(stability_rows) != 9 or any(r.get("seed_stability_status") != "PASS" for r in stability_rows):
        errors.append("Donor stability table count/status mismatch")

    control_audit = load(VALIDATION / "VP_G05_control_runs_v2_AUDIT01.json")
    if control_audit.get("audited_unique_control_tasks") != 90 or control_audit.get("numerical_equivalence_benchmarks_passed") != 6:
        errors.append("Control audit task/benchmark count mismatch")
    calibration_audit = load(VALIDATION / "VP_G05_control_calibration_v2_AUDIT01.json")
    if calibration_audit.get("recomputed_run_metric_rows") != 99 or calibration_audit.get("recomputed_donor_calibration_rows") != 63 or calibration_audit.get("recomputed_aggregate_rows") != 21:
        errors.append("Calibration recomputation row counts mismatch")

    calibration_path = ROOT / "06_consensus" / "negative_control_calibration_v2" / "negative_control_calibration.tsv"
    calibration_rows = rows(calibration_path)
    primary = [r for r in calibration_rows if r["metric"] == "top5pct_distance_sum"]
    if len(primary) != 3:
        errors.append("Primary calibration target count mismatch")
    primary_summary = [
        {
            "target_gene": r["target_gene"],
            "condition": r["condition"],
            "target_to_control_median_ratio": float(r["ratio_target_to_control_median"]),
            "empirical_percentile": float(r["empirical_percentile"]),
            "conservative_upper_tail_empirical_value": float(r["conservative_upper_tail_empirical_value"]),
            "matched_controls": int(r["control_count"]),
        }
        for r in sorted(primary, key=lambda x: x["target_gene"])
    ]

    retained_failures = [
        "VP_G05_degree_reference_ATTEMPT01_INTERRUPTED_RETAINED.json",
        "VP_G05_negative_control_selection_HOLD_INSUFFICIENT_MATCHED_CONTROLS_ATTEMPT01.json",
        "VP_G05_control_launcher_attempt_01_COMPATIBILITY_FAIL_RETAINED.json",
    ]
    for name in retained_failures:
        if not (VALIDATION / name).is_file():
            errors.append(f"Retained history missing: {name}")

    result = {
        "gate": "VP-G05",
        "stage": "stability_donor_consensus_and_matched_control_final_acceptance",
        "status": "PASS" if not errors else "FAIL_RETAINED",
        "gate_state": "COMPLETE_WITH_VERSIONED_CONTROL_AMENDMENT" if not errors else "INCOMPLETE",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "acceptance": {
            "seed_and_donor_stability": "PASS_3_OF_3_DONORS_FOR_ALL_TARGETS" if not errors else "SEE_ERRORS",
            "full_rankings_retained": len(ranking_rows) == 3 * 3 * 5 * 1388,
            "donor_consensus_rule": "top 5% in at least 2/3 donors",
            "original_v1_control_design": "HOLD_INSUFFICIENT_MATCHED_CONTROLS_RETAINED",
            "v2_control_selection": "PASS_VERSIONED_FEASIBILITY_AMENDMENT",
            "v2_control_runs": "PASS_90_OF_90_TECHNICAL",
            "shared_wt_direct_equivalence": "PASS_6_OF_6_AT_ABSOLUTE_TOLERANCE_1E-7",
            "non_pathway_empirical_calibration": "PASS_INDEPENDENT_RECOMPUTATION",
        },
        "primary_metric": "top5pct_distance_sum over ceil(5% of 1388)=70 genes per run",
        "primary_metric_summary": primary_summary,
        "evidence": evidence,
        "ranking_artifacts": {
            "all_virtual_KO_rankings.tsv.gz": {"rows": len(ranking_rows), "sha256": sha256(rankings_path)},
            "candidate_donor_consensus.tsv": {"rows": len(consensus_rows), "sha256": sha256(consensus_path)},
            "donor_stability.tsv": {"rows": len(stability_rows), "sha256": sha256(donor_stability_path)},
            "VP_G05_ranking_rule_freeze.json": {"sha256": sha256(ranking_freeze_path)},
        },
        "calibration_artifact": {"rows": len(calibration_rows), "sha256": sha256(calibration_path)},
        "retained_nonpassing_history": retained_failures,
        "scope_amendment": "The workbook's gene-set enrichment-strength calibration is deferred to VP-G06 because VP-G05 was frozen to close before pathway inspection. VP-G05 closes on non-pathway network-distance calibration only.",
        "scientific_boundary": "PASS establishes technical integrity, seed/donor stability, and descriptive comparison with relaxed v2 matched-reference genes. It does not establish biologically inert controls, pathway direction, expression reversal, therapeutic benefit, or causality.",
        "downstream_gate_status": "VP-G06_NOT_STARTED_REQUIRES_EXPLICIT_AUTHORIZATION",
        "errors": errors,
        "implementation_sha256": {"scripts/06e_finalize_vp_g05.py": sha256(Path(__file__))},
    }
    destination = COMPLETE if not errors else FAIL
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
