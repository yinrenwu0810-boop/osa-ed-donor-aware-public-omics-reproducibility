#!/usr/bin/env python3
"""Apply the frozen non-pathway VP-G05 v2 matched-control calibration."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
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
OUT = ROOT / "06_consensus" / "negative_control_calibration_v2"
METRICS = [
    "top5pct_distance_sum",
    "top5pct_distance_mean",
    "distance_sum_all",
    "distance_mean_all",
    "distance_median",
    "distance_q95",
    "distance_q99",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def summarize(path: Path) -> tuple[int, int, dict[str, float]]:
    rows = read_table(path)
    pairs = sorted(((float(row["distance"]), row["gene"]) for row in rows), key=lambda x: (-x[0], x[1]))
    distances = [item[0] for item in pairs]
    if not distances or any((not math.isfinite(x) or x < 0) for x in distances):
        raise ValueError(f"Invalid distance values in {path}")
    top_n = math.ceil(0.05 * len(distances))
    top = distances[:top_n]
    metrics = {
        "top5pct_distance_sum": sum(top),
        "top5pct_distance_mean": statistics.fmean(top),
        "distance_sum_all": sum(distances),
        "distance_mean_all": statistics.fmean(distances),
        "distance_median": statistics.median(distances),
        "distance_q95": quantile(distances, 0.95),
        "distance_q99": quantile(distances, 0.99),
    }
    return len(distances), top_n, metrics


def fmt(value: float) -> str:
    return format(value, ".17g")


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite sealed calibration directory: {OUT}")
    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    control_audit = json.loads(CONTROL_AUDIT.read_text(encoding="utf-8"))
    if freeze.get("status") != "FROZEN_BEFORE_CALIBRATION":
        raise SystemExit("Calibration rule is not frozen")
    if control_audit.get("status") != "PASS":
        raise SystemExit("Control-run audit is not PASS")
    if freeze.get("control_run_audit_sha256") != sha256(CONTROL_AUDIT):
        raise SystemExit("Control-run audit hash differs from frozen input")
    for script, expected in freeze.get("implementation_sha256", {}).items():
        if sha256(ROOT / script) != expected:
            raise SystemExit(f"Implementation hash mismatch: {script}")

    primary_seed = int(freeze["primary_seed"])
    selection_rows = read_table(SELECTION)
    task_rows = read_table(TASK_MATRIX)
    batch_rows = read_table(BATCH_MATRIX)
    registry = read_table(REGISTRY, delimiter=",")
    batches = {int(row["batch_order"]): row for row in batch_rows}

    targets: dict[tuple[str, str], dict] = {}
    for row in selection_rows:
        key = (row["target_gene"], row["condition"])
        entry = targets.setdefault(key, {"donors": [], "controls": {}})
        if row["donor"] not in entry["donors"]:
            entry["donors"].append(row["donor"])
        entry["controls"][int(row["selection_order"])] = row["control_gene"]
    for entry in targets.values():
        entry["donors"].sort()

    task_lookup: dict[tuple[str, str, str, int], Path] = {}
    for task in task_rows:
        batch = batches[int(task["batch_order"])]
        output = ROOT / batch["output_batch"] / "controls" / task["control_gene"] / "attempt_01" / "differential_regulation.tsv"
        for target in task["matched_targets"].split(";"):
            task_lookup[(target, task["condition"], task["donor"], int(task["seed"]))] = output

    tmp = OUT.with_name(OUT.name + f".tmp_{os.getpid()}")
    if tmp.exists():
        raise SystemExit(f"Temporary output already exists: {tmp}")
    tmp.mkdir(parents=True)

    run_rows: list[dict] = []
    numeric: dict[tuple[str, str, str, str, str], dict[str, float]] = {}
    source_hashes: dict[str, str] = {}
    for (target, condition), entry in sorted(targets.items()):
        controls = [entry["controls"][i] for i in sorted(entry["controls"])]
        if len(entry["donors"]) != 3 or len(controls) != 10:
            raise ValueError(f"Unexpected donor/control count for {target}")
        for donor in entry["donors"]:
            matches = [
                row for row in registry
                if row["gene"] == target and row["condition"] == condition and row["donor"] == donor
                and int(row["seed"]) == primary_seed and row["status"] == "PASS_TECHNICAL"
            ]
            if len(matches) != 1:
                raise ValueError(f"Expected one primary target run for {target}/{donor}, found {len(matches)}")
            target_path = Path(matches[0]["output_attempt"]) / "differential_regulation.tsv"
            sources = [("target", target, target_path)]
            for control in controls:
                control_path = task_lookup.get((target, condition, donor, primary_seed))
                # The physical task lookup is target-specific; resolve the named control directly from its batch.
                task_match = [
                    row for row in task_rows
                    if row["control_gene"] == control and row["condition"] == condition
                    and row["donor"] == donor and int(row["seed"]) == primary_seed
                    and target in row["matched_targets"].split(";")
                ]
                if len(task_match) != 1:
                    raise ValueError(f"Expected one control task for {target}/{control}/{donor}")
                batch = batches[int(task_match[0]["batch_order"])]
                control_path = ROOT / batch["output_batch"] / "controls" / control / "attempt_01" / "differential_regulation.tsv"
                sources.append(("control", control, control_path))

            for role, perturbation, path in sources:
                if not path.is_file():
                    raise FileNotFoundError(path)
                n_genes, top_n, metrics = summarize(path)
                relative = path.relative_to(ROOT).as_posix()
                source_hashes[relative] = sha256(path)
                numeric[(target, condition, donor, role, perturbation)] = metrics
                row = {
                    "role": role,
                    "target_gene": target,
                    "condition": condition,
                    "donor": donor,
                    "perturbation_gene": perturbation,
                    "source_relative_path": relative,
                    "source_sha256": source_hashes[relative],
                    "n_genes": n_genes,
                    "top_n": top_n,
                }
                row.update({metric: fmt(metrics[metric]) for metric in METRICS})
                run_rows.append(row)

    donor_rows: list[dict] = []
    aggregate_rows: list[dict] = []
    report_targets: list[dict] = []
    for (target, condition), entry in sorted(targets.items()):
        controls = [entry["controls"][i] for i in sorted(entry["controls"])]
        for donor in entry["donors"]:
            for metric in METRICS:
                target_value = numeric[(target, condition, donor, "target", target)][metric]
                control_values = [numeric[(target, condition, donor, "control", control)][metric] for control in controls]
                control_median = statistics.median(control_values)
                donor_rows.append({
                    "target_gene": target,
                    "condition": condition,
                    "donor": donor,
                    "metric": metric,
                    "target_value": fmt(target_value),
                    "control_median": fmt(control_median),
                    "ratio_target_to_control_median": fmt(target_value / control_median) if control_median != 0 else "NA",
                    "controls_le_target": sum(value <= target_value for value in control_values),
                    "controls_ge_target": sum(value >= target_value for value in control_values),
                    "control_count": len(control_values),
                })

        primary_summary = None
        for metric in METRICS:
            target_values = [numeric[(target, condition, donor, "target", target)][metric] for donor in entry["donors"]]
            target_median = statistics.median(target_values)
            control_medians = [
                statistics.median([numeric[(target, condition, donor, "control", control)][metric] for donor in entry["donors"]])
                for control in controls
            ]
            le = sum(value <= target_median for value in control_medians)
            ge = sum(value >= target_median for value in control_medians)
            control_center = statistics.median(control_medians)
            aggregate = {
                "target_gene": target,
                "condition": condition,
                "metric": metric,
                "target_median_across_donors": fmt(target_median),
                "control_median_of_donor_medians": fmt(control_center),
                "ratio_target_to_control_median": fmt(target_median / control_center) if control_center != 0 else "NA",
                "controls_le_target": le,
                "controls_ge_target": ge,
                "control_count": len(control_medians),
                "empirical_percentile": fmt(le / len(control_medians)),
                "conservative_upper_tail_empirical_value": fmt((1 + ge) / (len(control_medians) + 1)),
            }
            aggregate_rows.append(aggregate)
            if metric == freeze["primary_metric"]:
                primary_summary = aggregate.copy()
        report_targets.append(primary_summary)

    run_fields = [
        "role", "target_gene", "condition", "donor", "perturbation_gene", "source_relative_path",
        "source_sha256", "n_genes", "top_n", *METRICS,
    ]
    donor_fields = [
        "target_gene", "condition", "donor", "metric", "target_value", "control_median",
        "ratio_target_to_control_median", "controls_le_target", "controls_ge_target", "control_count",
    ]
    aggregate_fields = [
        "target_gene", "condition", "metric", "target_median_across_donors",
        "control_median_of_donor_medians", "ratio_target_to_control_median", "controls_le_target",
        "controls_ge_target", "control_count", "empirical_percentile",
        "conservative_upper_tail_empirical_value",
    ]
    write_table(tmp / "run_metric_summaries.tsv", run_rows, run_fields)
    write_table(tmp / "donor_calibration.tsv", donor_rows, donor_fields)
    write_table(tmp / "negative_control_calibration.tsv", aggregate_rows, aggregate_fields)

    report = {
        "gate": "VP-G05",
        "stage": "non_pathway_matched_control_calibration_v2",
        "status": "PASS_CALIBRATION_COMPUTED",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rule_freeze_sha256": sha256(FREEZE),
        "control_run_audit_sha256": sha256(CONTROL_AUDIT),
        "primary_seed": primary_seed,
        "donor_is_independent_unit": True,
        "targets": 3,
        "target_donor_runs": 9,
        "analytical_control_donor_runs": 90,
        "metrics": METRICS,
        "primary_metric": freeze["primary_metric"],
        "primary_metric_summary": report_targets,
        "source_artifact_count": len(source_hashes),
        "pathway_calibration": "DEFERRED_TO_VP_G06",
        "interpretation": "Descriptive matched-reference calibration only. Empirical percentiles and upper-tail values are not causal estimates or confirmatory frequentist p-values.",
    }
    (tmp / "calibration_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    manifest_rows = []
    for name in ["run_metric_summaries.tsv", "donor_calibration.tsv", "negative_control_calibration.tsv", "calibration_report.json"]:
        path = tmp / name
        manifest_rows.append({"file": name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_table(tmp / "artifact_manifest.sha256.tsv", manifest_rows, ["file", "bytes", "sha256"])
    tmp.rename(OUT)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
