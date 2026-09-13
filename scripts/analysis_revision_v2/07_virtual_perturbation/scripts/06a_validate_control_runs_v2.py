#!/usr/bin/env python3
"""Independent structural and numerical audit of all VP-G05 v2 control KO runs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "01_protocol" / "control_execution_v2"
RUN_ROOT = ROOT / "05_runs" / "control_batches_v2"
OUT = ROOT / "validation" / "VP_G05_control_runs_v2_AUDIT01.json"
EXPECTED_COLUMNS = ["gene", "distance", "Z", "FC", "p.value", "p.adj"]
REQUIRED_ARTIFACTS = {
    "differential_regulation.tsv",
    "run_parameters.json",
    "session_info.txt",
    "stderr.log",
    "stdout.log",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite sealed audit: {OUT}")

    errors: list[str] = []
    batch_rows = read_tsv(PROTOCOL / "control_batch_matrix.tsv")
    task_rows = read_tsv(PROTOCOL / "control_task_matrix.tsv")
    freeze_path = PROTOCOL / "VP_G05_control_execution_freeze_v2.json"
    freeze_hash = sha256(freeze_path)

    if len(batch_rows) != 6:
        errors.append(f"Expected 6 batches, found {len(batch_rows)}")
    if len(task_rows) != 90:
        errors.append(f"Expected 90 control tasks, found {len(task_rows)}")

    tasks_by_batch: dict[int, list[dict[str, str]]] = {}
    for row in task_rows:
        tasks_by_batch.setdefault(int(row["batch_order"]), []).append(row)

    audited_tasks = 0
    audited_batches = 0
    benchmark_passes = 0
    row_counts: set[int] = set()
    seen_keys: set[tuple[str, str, str, int]] = set()

    for batch in batch_rows:
        order = int(batch["batch_order"])
        condition = batch["condition"]
        donor = batch["donor"]
        seed = int(batch["seed"])
        expected_tasks = tasks_by_batch.get(order, [])
        batch_dir = ROOT / Path(batch["output_batch"])
        audited_batches += 1

        if not batch_dir.is_dir():
            errors.append(f"Missing batch directory: {rel(batch_dir)}")
            continue

        complete_path = batch_dir / "batch_complete.json"
        params_path = batch_dir / "batch_parameters.json"
        benchmark_path = batch_dir / "direct_equivalence_benchmark.json"
        progress_path = batch_dir / "progress.tsv"
        shared_wt_path = batch_dir / "shared_wt_reference.rds"
        for required in [complete_path, params_path, benchmark_path, progress_path, shared_wt_path]:
            if not required.is_file():
                errors.append(f"Missing batch artifact: {rel(required)}")
        if any(not p.is_file() for p in [complete_path, params_path, benchmark_path, progress_path, shared_wt_path]):
            continue

        complete = load_json(complete_path)
        params = load_json(params_path)
        benchmark = load_json(benchmark_path)
        progress = read_tsv(progress_path)
        shared_wt_hash = sha256(shared_wt_path)
        benchmark_hash = sha256(benchmark_path)

        if complete.get("gate") != "VP-G05" or complete.get("status") != "PASS_TECHNICAL":
            errors.append(f"Batch {order}: invalid completion gate/status")
        for field, expected in [("condition", condition), ("donor", donor), ("seed", seed)]:
            if complete.get(field) != expected:
                errors.append(f"Batch {order}: completion {field} mismatch")
            if params.get(field) != expected:
                errors.append(f"Batch {order}: parameters {field} mismatch")
        if complete.get("controls_expected") != len(expected_tasks):
            errors.append(f"Batch {order}: controls_expected mismatch")
        if complete.get("controls_completed") != len(expected_tasks):
            errors.append(f"Batch {order}: controls_completed mismatch")
        if complete.get("shared_wt_sha256") != shared_wt_hash:
            errors.append(f"Batch {order}: shared WT hash mismatch")
        if complete.get("benchmark_sha256") != benchmark_hash:
            errors.append(f"Batch {order}: benchmark hash mismatch")
        if params.get("freeze_sha256") != freeze_hash:
            errors.append(f"Batch {order}: execution freeze hash mismatch")

        if benchmark.get("gate") != "VP-G05" or benchmark.get("status") != "PASS_NUMERICAL_EQUIVALENCE":
            errors.append(f"Batch {order}: numerical-equivalence benchmark not PASS")
        else:
            benchmark_passes += 1
        tolerance = float(benchmark.get("tolerance_absolute", float("nan")))
        if not math.isfinite(tolerance) or tolerance != 1e-7:
            errors.append(f"Batch {order}: unexpected equivalence tolerance {tolerance}")
        for metric, value in benchmark.get("maximum_absolute_difference", {}).items():
            numeric = float(value)
            if not math.isfinite(numeric) or numeric > tolerance:
                errors.append(f"Batch {order}: benchmark {metric} exceeds tolerance")
        if benchmark.get("shared_wt_sha256") != shared_wt_hash:
            errors.append(f"Batch {order}: benchmark shared WT hash mismatch")
        direct_path = Path(str(benchmark.get("direct_result", "")))
        if not direct_path.is_file() or benchmark.get("direct_result_sha256") != sha256(direct_path):
            errors.append(f"Batch {order}: direct benchmark source/hash mismatch")

        expected_progress = {(r["control_gene"], r["matched_targets"]) for r in expected_tasks}
        actual_progress = {(r["control_gene"], r["matched_targets"]) for r in progress}
        if len(progress) != len(expected_tasks) or actual_progress != expected_progress:
            errors.append(f"Batch {order}: progress membership/count mismatch")
        if any(r.get("status") != "PASS_TECHNICAL" for r in progress):
            errors.append(f"Batch {order}: progress contains non-PASS task")

        for task in expected_tasks:
            gene = task["control_gene"]
            matched_targets = task["matched_targets"]
            key = (gene, condition, donor, seed)
            if key in seen_keys:
                errors.append(f"Duplicate task key: {key}")
            seen_keys.add(key)
            attempt = batch_dir / "controls" / gene / "attempt_01"
            manifest_path = attempt / "artifact_manifest.sha256.tsv"
            if not manifest_path.is_file():
                errors.append(f"Missing manifest: {rel(manifest_path)}")
                continue

            manifest = read_tsv(manifest_path)
            manifest_names = {r["file"] for r in manifest}
            if manifest_names != REQUIRED_ARTIFACTS:
                errors.append(f"{gene}/{donor}: manifest membership mismatch")
            for item in manifest:
                artifact = attempt / item["file"]
                if not artifact.is_file():
                    errors.append(f"Missing listed artifact: {rel(artifact)}")
                    continue
                if int(item["bytes"]) != artifact.stat().st_size:
                    errors.append(f"{gene}/{donor}: byte size mismatch for {item['file']}")
                if item["sha256"] != sha256(artifact):
                    errors.append(f"{gene}/{donor}: SHA-256 mismatch for {item['file']}")

            run_params_path = attempt / "run_parameters.json"
            diff_path = attempt / "differential_regulation.tsv"
            if not run_params_path.is_file() or not diff_path.is_file():
                continue
            run_params = load_json(run_params_path)
            expected_fields = {
                "gate": "VP-G05",
                "condition": condition,
                "donor": donor,
                "control_gene": gene,
                "matched_targets": matched_targets,
                "seed": seed,
                "shared_wt_sha256": shared_wt_hash,
                "direct_equivalence_benchmark_sha256": benchmark_hash,
                "freeze_sha256": freeze_hash,
            }
            for field, expected in expected_fields.items():
                if run_params.get(field) != expected:
                    errors.append(f"{gene}/{donor}: run parameter {field} mismatch")

            with diff_path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != EXPECTED_COLUMNS:
                    errors.append(f"{gene}/{donor}: differential-regulation schema mismatch")
                rows = list(reader)
            row_counts.add(len(rows))
            if len(rows) != 1388:
                errors.append(f"{gene}/{donor}: expected 1388 differential rows, found {len(rows)}")
            genes = [r["gene"] for r in rows]
            if len(set(genes)) != len(genes) or any(not g for g in genes):
                errors.append(f"{gene}/{donor}: blank or duplicated gene identifiers")
            for row_number, row in enumerate(rows, start=2):
                for column in EXPECTED_COLUMNS[1:]:
                    try:
                        value = float(row[column])
                    except (TypeError, ValueError):
                        errors.append(f"{gene}/{donor}:{row_number}: nonnumeric {column}")
                        continue
                    if not math.isfinite(value):
                        errors.append(f"{gene}/{donor}:{row_number}: nonfinite {column}")
                    if column in {"distance", "p.value", "p.adj"} and value < 0:
                        errors.append(f"{gene}/{donor}:{row_number}: negative {column}")
                    if column in {"p.value", "p.adj"} and value > 1:
                        errors.append(f"{gene}/{donor}:{row_number}: {column} exceeds 1")
            audited_tasks += 1

    report = {
        "gate": "VP-G05",
        "stage": "matched_control_KO_v2_independent_audit",
        "status": "PASS" if not errors else "FAIL_RETAINED",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_freeze_sha256": freeze_hash,
        "expected_batches": 6,
        "audited_batches": audited_batches,
        "expected_unique_control_tasks": 90,
        "audited_unique_control_tasks": audited_tasks,
        "numerical_equivalence_benchmarks_passed": benchmark_passes,
        "differential_row_counts_observed": sorted(row_counts),
        "errors": errors,
        "boundary": "Artifact, provenance, schema, and numerical-equivalence audit only; no pathway result or biological effect was interpreted.",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
