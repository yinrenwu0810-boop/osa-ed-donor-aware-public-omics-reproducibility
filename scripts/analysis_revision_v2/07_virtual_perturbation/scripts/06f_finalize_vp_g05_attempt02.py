#!/usr/bin/env python3
"""Corrected append-only VP-G05 final audit after the retained attempt-01 count error."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "validation"
OUT = VAL / "VP_G05_stability_COMPLETE.json"
ATTEMPT1 = VAL / "VP_G05_stability_COMPLETE_ATTEMPT01_FAIL_RETAINED.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def table(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite sealed completion record: {OUT}")
    errors: list[str] = []
    failed = load(ATTEMPT1)
    expected_error = ["Full ranking row count mismatch: 12492"]
    if failed.get("status") != "FAIL_RETAINED" or failed.get("errors") != expected_error:
        errors.append("Attempt-01 retained failure is not the expected single count-model error")
    old_script = ROOT / "scripts" / "06e_finalize_vp_g05.py"
    if failed.get("implementation_sha256", {}).get("scripts/06e_finalize_vp_g05.py") != sha256(old_script):
        errors.append("Attempt-01 auditor no longer matches its recorded SHA-256")

    required = {
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
    evidence = {}
    for name, expected in required.items():
        path = VAL / name
        obj = load(path)
        evidence[name] = {"status": obj.get("status"), "sha256": sha256(path)}
        if obj.get("status") != expected:
            errors.append(f"Prerequisite status mismatch: {name}")

    rankings_path = ROOT / "06_consensus" / "all_virtual_KO_rankings.tsv.gz"
    consensus_path = ROOT / "06_consensus" / "candidate_donor_consensus.tsv"
    stability_path = ROOT / "06_consensus" / "donor_stability.tsv"
    rankings = table(rankings_path)
    consensus = table(consensus_path)
    stability = table(stability_path)
    expected_targets = {"TYMS", "EFNB2", "LRRC17"}
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rankings:
        groups[(row["target_gene"], row["donor"])].append(row)
    if len(rankings) != 3 * 3 * 1388 or len(groups) != 9:
        errors.append("Seed-aggregated full ranking dimensions are not 3 targets x 3 donors x 1388 genes")
    if {key[0] for key in groups} != expected_targets:
        errors.append("Full ranking target membership mismatch")
    for key, group in groups.items():
        genes = {r["gene"] for r in group}
        ranks = {int(r["consensus_rank"]) for r in group}
        top_count = sum(r["in_donor_top5pct"] == "True" for r in group)
        if len(group) != 1388 or len(genes) != 1388 or ranks != set(range(1, 1389)) or top_count != 70:
            errors.append(f"Full ranking group integrity mismatch: {key}")

    consensus_groups = Counter(r["target_gene"] for r in consensus)
    if len(consensus) != 3 * 1388 or consensus_groups != Counter({gene: 1388 for gene in expected_targets}):
        errors.append("Across-donor consensus dimensions/membership mismatch")
    if any(int(r["donors_total"]) != 3 or int(r["donors_top5pct"]) not in range(4) for r in consensus):
        errors.append("Across-donor consensus donor counts invalid")
    if any((r["is_2of3_top5pct_consensus"] == "True") != (int(r["donors_top5pct"]) >= 2) for r in consensus):
        errors.append("Frozen 2-of-3 donor consensus flag mismatch")
    if len(stability) != 9 or any(r["seed_stability_status"] != "PASS" for r in stability):
        errors.append("Seed stability is not PASS for all nine target-donor units")

    control_audit = load(VAL / "VP_G05_control_runs_v2_AUDIT01.json")
    calibration_audit = load(VAL / "VP_G05_control_calibration_v2_AUDIT01.json")
    if (control_audit.get("audited_unique_control_tasks"), control_audit.get("numerical_equivalence_benchmarks_passed")) != (90, 6):
        errors.append("Control run audit counts mismatch")
    if (calibration_audit.get("recomputed_run_metric_rows"), calibration_audit.get("recomputed_donor_calibration_rows"), calibration_audit.get("recomputed_aggregate_rows")) != (99, 63, 21):
        errors.append("Calibration audit recomputation counts mismatch")

    cal_path = ROOT / "06_consensus" / "negative_control_calibration_v2" / "negative_control_calibration.tsv"
    calibration = table(cal_path)
    primary = sorted((r for r in calibration if r["metric"] == "top5pct_distance_sum"), key=lambda r: r["target_gene"])
    if len(primary) != 3:
        errors.append("Primary calibration target count mismatch")
    summary = [{
        "target_gene": r["target_gene"],
        "condition": r["condition"],
        "target_to_control_median_ratio": float(r["ratio_target_to_control_median"]),
        "empirical_percentile": float(r["empirical_percentile"]),
        "conservative_upper_tail_empirical_value": float(r["conservative_upper_tail_empirical_value"]),
        "matched_controls": int(r["control_count"]),
    } for r in primary]

    retained = [
        "VP_G05_degree_reference_ATTEMPT01_INTERRUPTED_RETAINED.json",
        "VP_G05_negative_control_selection_HOLD_INSUFFICIENT_MATCHED_CONTROLS_ATTEMPT01.json",
        "VP_G05_control_launcher_attempt_01_COMPATIBILITY_FAIL_RETAINED.json",
        "VP_G05_stability_COMPLETE_ATTEMPT01_FAIL_RETAINED.json",
    ]
    if any(not (VAL / name).is_file() for name in retained):
        errors.append("One or more required nonpassing history records are missing")

    result = {
        "gate": "VP-G05",
        "stage": "stability_donor_consensus_and_matched_control_final_acceptance",
        "attempt": 2,
        "status": "PASS" if not errors else "FAIL_RETAINED",
        "gate_state": "COMPLETE_WITH_VERSIONED_CONTROL_AMENDMENT" if not errors else "INCOMPLETE",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "attempt_01_correction": {
            "retained_record_sha256": sha256(ATTEMPT1),
            "reason": "Attempt 01 incorrectly expected per-seed rather than frozen seed-aggregated ranking rows. The correct expected count is 3 targets x 3 donors x 1388 genes = 12492.",
            "data_artifacts_changed": False,
        },
        "acceptance": {
            "seed_and_donor_stability": "PASS_3_OF_3_DONORS_FOR_ALL_TARGETS",
            "full_seed_aggregated_rankings": "PASS_12492_ROWS",
            "donor_consensus_rule": "PASS_TOP_5_PERCENT_IN_AT_LEAST_2_OF_3_DONORS",
            "original_v1_control_design": "HOLD_INSUFFICIENT_MATCHED_CONTROLS_RETAINED",
            "v2_control_selection": "PASS_VERSIONED_FEASIBILITY_AMENDMENT",
            "v2_control_runs": "PASS_90_OF_90_TECHNICAL",
            "shared_wt_direct_equivalence": "PASS_6_OF_6_AT_ABSOLUTE_TOLERANCE_1E-7",
            "non_pathway_empirical_calibration": "PASS_INDEPENDENT_RECOMPUTATION",
        },
        "primary_metric": "top5pct_distance_sum over ceil(5% of 1388)=70 genes per run",
        "primary_metric_summary": summary,
        "evidence": evidence,
        "ranking_artifacts": {
            "all_virtual_KO_rankings.tsv.gz": {"rows": len(rankings), "sha256": sha256(rankings_path)},
            "candidate_donor_consensus.tsv": {"rows": len(consensus), "sha256": sha256(consensus_path)},
            "donor_stability.tsv": {"rows": len(stability), "sha256": sha256(stability_path)},
        },
        "calibration_artifact": {"rows": len(calibration), "sha256": sha256(cal_path)},
        "retained_nonpassing_history": retained,
        "scope_amendment": "Gene-set enrichment-strength calibration is deferred to VP-G06 because the frozen gate boundary closes VP-G05 before pathway inspection. VP-G05 closes on non-pathway network-distance calibration only.",
        "scientific_boundary": "PASS establishes technical integrity, seed/donor stability, and descriptive comparison with relaxed v2 matched-reference genes. It does not establish biologically inert controls, pathway direction, expression reversal, therapeutic benefit, or causality.",
        "downstream_gate_status": "VP-G06_NOT_STARTED_REQUIRES_EXPLICIT_AUTHORIZATION",
        "errors": errors,
        "implementation_sha256": {"scripts/06f_finalize_vp_g05_attempt02.py": sha256(Path(__file__))},
    }
    if errors:
        failure = VAL / "VP_G05_stability_COMPLETE_ATTEMPT02_FAIL_RETAINED.json"
        result["status"] = "FAIL_RETAINED"
        result["gate_state"] = "INCOMPLETE"
        failure.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, ensure_ascii=False, indent=2))
        raise SystemExit(1)
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
