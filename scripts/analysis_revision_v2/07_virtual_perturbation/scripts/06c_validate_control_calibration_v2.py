#!/usr/bin/env python3
"""Independent recomputation audit for VP-G05 v2 matched-control calibration."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "01_protocol" / "VP_G05_control_calibration_rule_freeze_v2.json"
CONTROL_AUDIT = ROOT / "validation" / "VP_G05_control_runs_v2_AUDIT01.json"
SELECTION = ROOT / "06_consensus" / "negative_control_selection_v2" / "negative_control_selection.tsv"
TASK_MATRIX = ROOT / "01_protocol" / "control_execution_v2" / "control_task_matrix.tsv"
BATCH_MATRIX = ROOT / "01_protocol" / "control_execution_v2" / "control_batch_matrix.tsv"
REGISTRY = ROOT / "05_runs" / "run_registry.tsv"
CAL = ROOT / "06_consensus" / "negative_control_calibration_v2"
OUT = ROOT / "validation" / "VP_G05_control_calibration_v2_AUDIT01.json"
METRICS = ["top5pct_distance_sum", "top5pct_distance_mean", "distance_sum_all", "distance_mean_all", "distance_median", "distance_q95", "distance_q99"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def quantile(values: list[float], p: float) -> float:
    values = sorted(values)
    pos = (len(values) - 1) * p
    lo, hi = math.floor(pos), math.ceil(pos)
    return values[lo] if lo == hi else values[lo] + (values[hi] - values[lo]) * (pos - lo)


def metrics(path: Path) -> tuple[int, int, dict[str, float]]:
    rows = table(path)
    ordered = sorted(((float(r["distance"]), r["gene"]) for r in rows), key=lambda x: (-x[0], x[1]))
    values = [x[0] for x in ordered]
    top_n = math.ceil(len(values) * 0.05)
    top = values[:top_n]
    return len(values), top_n, {
        "top5pct_distance_sum": sum(top),
        "top5pct_distance_mean": statistics.fmean(top),
        "distance_sum_all": sum(values),
        "distance_mean_all": statistics.fmean(values),
        "distance_median": statistics.median(values),
        "distance_q95": quantile(values, 0.95),
        "distance_q99": quantile(values, 0.99),
    }


def close(a: str, b: float) -> bool:
    return math.isclose(float(a), b, rel_tol=1e-12, abs_tol=1e-15)


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite sealed audit: {OUT}")
    errors: list[str] = []
    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    report = json.loads((CAL / "calibration_report.json").read_text(encoding="utf-8"))
    manifest = table(CAL / "artifact_manifest.sha256.tsv")
    expected_manifest = {"run_metric_summaries.tsv", "donor_calibration.tsv", "negative_control_calibration.tsv", "calibration_report.json"}
    if {r["file"] for r in manifest} != expected_manifest:
        errors.append("Calibration manifest membership mismatch")
    for row in manifest:
        path = CAL / row["file"]
        if not path.is_file() or int(row["bytes"]) != path.stat().st_size or row["sha256"] != sha256(path):
            errors.append(f"Calibration manifest mismatch: {row['file']}")
    if report.get("status") != "PASS_CALIBRATION_COMPUTED" or report.get("rule_freeze_sha256") != sha256(FREEZE):
        errors.append("Calibration report status/freeze provenance mismatch")
    if report.get("control_run_audit_sha256") != sha256(CONTROL_AUDIT):
        errors.append("Calibration report control-audit provenance mismatch")
    for script, expected in freeze.get("implementation_sha256", {}).items():
        if sha256(ROOT / script) != expected:
            errors.append(f"Frozen implementation hash mismatch: {script}")

    selection = table(SELECTION)
    tasks = table(TASK_MATRIX)
    batches = {int(r["batch_order"]): r for r in table(BATCH_MATRIX)}
    registry = table(REGISTRY, delimiter=",")
    primary_seed = int(freeze["primary_seed"])
    targets: dict[tuple[str, str], dict] = {}
    for row in selection:
        key = (row["target_gene"], row["condition"])
        obj = targets.setdefault(key, {"donors": set(), "controls": {}})
        obj["donors"].add(row["donor"])
        obj["controls"][int(row["selection_order"])] = row["control_gene"]

    produced_runs = table(CAL / "run_metric_summaries.tsv")
    produced_index = {(r["target_gene"], r["condition"], r["donor"], r["role"], r["perturbation_gene"]): r for r in produced_runs}
    recomputed: dict[tuple[str, str, str, str, str], dict[str, float]] = {}
    expected_run_count = 0
    for (target, condition), spec in sorted(targets.items()):
        donors = sorted(spec["donors"])
        controls = [spec["controls"][i] for i in sorted(spec["controls"])]
        for donor in donors:
            matches = [r for r in registry if r["gene"] == target and r["condition"] == condition and r["donor"] == donor and int(r["seed"]) == primary_seed and r["status"] == "PASS_TECHNICAL"]
            if len(matches) != 1:
                errors.append(f"Target registry membership mismatch: {target}/{donor}")
                continue
            sources = [("target", target, Path(matches[0]["output_attempt"]) / "differential_regulation.tsv")]
            for control in controls:
                found = [r for r in tasks if r["control_gene"] == control and r["condition"] == condition and r["donor"] == donor and int(r["seed"]) == primary_seed and target in r["matched_targets"].split(";")]
                if len(found) != 1:
                    errors.append(f"Control task membership mismatch: {target}/{control}/{donor}")
                    continue
                batch = batches[int(found[0]["batch_order"])]
                path = ROOT / batch["output_batch"] / "controls" / control / "attempt_01" / "differential_regulation.tsv"
                sources.append(("control", control, path))
            for role, gene, path in sources:
                expected_run_count += 1
                key = (target, condition, donor, role, gene)
                if key not in produced_index:
                    errors.append(f"Missing run metric row: {key}")
                    continue
                row = produced_index[key]
                n, top_n, values = metrics(path)
                recomputed[key] = values
                if int(row["n_genes"]) != n or int(row["top_n"]) != top_n or row["source_sha256"] != sha256(path):
                    errors.append(f"Run provenance/count mismatch: {key}")
                for metric in METRICS:
                    if not close(row[metric], values[metric]):
                        errors.append(f"Run metric mismatch: {key}/{metric}")
    if len(produced_runs) != expected_run_count or len(produced_index) != expected_run_count:
        errors.append(f"Run metric row count/uniqueness mismatch: expected {expected_run_count}, found {len(produced_runs)}/{len(produced_index)}")

    donor_rows = table(CAL / "donor_calibration.tsv")
    donor_index = {(r["target_gene"], r["condition"], r["donor"], r["metric"]): r for r in donor_rows}
    aggregate_rows = table(CAL / "negative_control_calibration.tsv")
    aggregate_index = {(r["target_gene"], r["condition"], r["metric"]): r for r in aggregate_rows}
    expected_donor = 0
    expected_aggregate = 0
    for (target, condition), spec in sorted(targets.items()):
        donors = sorted(spec["donors"])
        controls = [spec["controls"][i] for i in sorted(spec["controls"])]
        for donor in donors:
            for metric in METRICS:
                expected_donor += 1
                row = donor_index.get((target, condition, donor, metric))
                if row is None:
                    errors.append(f"Missing donor calibration row: {target}/{donor}/{metric}")
                    continue
                target_value = recomputed[(target, condition, donor, "target", target)][metric]
                control_values = [recomputed[(target, condition, donor, "control", c)][metric] for c in controls]
                center = statistics.median(control_values)
                if not close(row["target_value"], target_value) or not close(row["control_median"], center):
                    errors.append(f"Donor value mismatch: {target}/{donor}/{metric}")
                expected_ratio = target_value / center if center else None
                if expected_ratio is None:
                    if row["ratio_target_to_control_median"] != "NA": errors.append(f"Donor zero ratio mismatch: {target}/{donor}/{metric}")
                elif not close(row["ratio_target_to_control_median"], expected_ratio):
                    errors.append(f"Donor ratio mismatch: {target}/{donor}/{metric}")
                if int(row["controls_le_target"]) != sum(x <= target_value for x in control_values) or int(row["controls_ge_target"]) != sum(x >= target_value for x in control_values) or int(row["control_count"]) != 10:
                    errors.append(f"Donor empirical count mismatch: {target}/{donor}/{metric}")
        for metric in METRICS:
            expected_aggregate += 1
            row = aggregate_index.get((target, condition, metric))
            if row is None:
                errors.append(f"Missing aggregate row: {target}/{metric}")
                continue
            target_median = statistics.median([recomputed[(target, condition, d, "target", target)][metric] for d in donors])
            control_medians = [statistics.median([recomputed[(target, condition, d, "control", c)][metric] for d in donors]) for c in controls]
            center = statistics.median(control_medians)
            le, ge = sum(x <= target_median for x in control_medians), sum(x >= target_median for x in control_medians)
            numeric_checks = {
                "target_median_across_donors": target_median,
                "control_median_of_donor_medians": center,
                "ratio_target_to_control_median": target_median / center,
                "empirical_percentile": le / 10,
                "conservative_upper_tail_empirical_value": (1 + ge) / 11,
            }
            if any(not close(row[field], value) for field, value in numeric_checks.items()):
                errors.append(f"Aggregate numeric mismatch: {target}/{metric}")
            if int(row["controls_le_target"]) != le or int(row["controls_ge_target"]) != ge or int(row["control_count"]) != 10:
                errors.append(f"Aggregate empirical count mismatch: {target}/{metric}")
    if len(donor_rows) != expected_donor or len(donor_index) != expected_donor:
        errors.append("Donor calibration row count/uniqueness mismatch")
    if len(aggregate_rows) != expected_aggregate or len(aggregate_index) != expected_aggregate:
        errors.append("Aggregate calibration row count/uniqueness mismatch")

    audit = {
        "gate": "VP-G05",
        "stage": "non_pathway_matched_control_calibration_v2_independent_audit",
        "status": "PASS" if not errors else "FAIL_RETAINED",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rule_freeze_sha256": sha256(FREEZE),
        "calibration_manifest_sha256": sha256(CAL / "artifact_manifest.sha256.tsv"),
        "recomputed_run_metric_rows": expected_run_count,
        "recomputed_donor_calibration_rows": expected_donor,
        "recomputed_aggregate_rows": expected_aggregate,
        "errors": errors,
        "boundary": "Independent numerical recomputation of non-pathway matched-reference calibration; no pathway result, causal effect, or confirmatory significance claim was assessed.",
    }
    OUT.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
