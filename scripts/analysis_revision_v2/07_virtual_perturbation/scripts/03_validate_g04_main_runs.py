"""Independent technical audit for frozen VP-G04 main virtual-KO runs.

This validator never changes run outputs.  It checks the frozen run matrix,
registry, required artifacts, artifact manifests, and per-run parameters.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
MATRIX = VP / "01_protocol" / "run_matrix.tsv"
REGISTRY = VP / "05_runs" / "run_registry.tsv"
VALIDATION = VP / "validation"
ERRATUM = VALIDATION / "VP_G04_metadata_erratum_v1.json"
RUNNER_SOURCE = VP / "scripts" / "02_run_sctenifoldknk.R"
REQUIRED = (
    "differential_regulation.tsv",
    "run_parameters.json",
    "session_info.txt",
    "stdout.log",
    "stderr.log",
    "resource_usage.tsv",
    "artifact_manifest.sha256.tsv",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def norm(path_text: str) -> Path:
    return Path(path_text.replace("\\", "/")).resolve()


def main() -> int:
    matrix = read_table(MATRIX, "\t")
    registry = read_table(REGISTRY, ",")
    errors: list[dict[str, object]] = []
    warnings: list[dict[str, object]] = []
    metadata_holds: list[dict[str, object]] = []
    metadata_erratum_verified: list[dict[str, object]] = []
    erratum_entries: dict[str, dict[str, object]] = {}
    if ERRATUM.is_file():
        try:
            erratum = json.loads(ERRATUM.read_text(encoding="utf-8"))
            erratum_entries = {entry["run_id"]: entry for entry in erratum.get("entries", [])}
            if erratum.get("registry_sha256") != sha256(REGISTRY):
                errors.append({"check": "metadata_erratum_registry_binding"})
            if erratum.get("runner_source_sha256_after_correction") != sha256(RUNNER_SOURCE):
                errors.append({"check": "metadata_erratum_corrected_runner_binding"})
        except (json.JSONDecodeError, OSError, KeyError):
            errors.append({"check": "metadata_erratum_readable"})

    expected_orders = list(range(1, 46))
    matrix_orders = [int(row["run_order"]) for row in matrix]
    registry_orders = [int(row["run_order"]) for row in registry]
    if matrix_orders != expected_orders:
        errors.append({"check": "frozen_matrix_orders", "observed": matrix_orders})
    if sorted(registry_orders) != expected_orders or len(registry_orders) != len(set(registry_orders)):
        errors.append({"check": "registry_unique_final_status", "observed_orders": registry_orders})
    if any(row["status"] != "PASS_TECHNICAL" for row in registry):
        errors.append({"check": "registry_all_technical_pass", "counts": Counter(row["status"] for row in registry)})

    matrix_by_order = {int(row["run_order"]): row for row in matrix}
    matrix_by_run_id = {row["run_id"]: row for row in matrix}
    paused_run_ids = set()
    for pause_event in (VP / "05_runs" / "pause_events").glob("*_PAUSED_*.json"):
        try:
            paused_run_ids.add(json.loads(pause_event.read_text(encoding="utf-8")).get("run_id"))
        except (json.JSONDecodeError, OSError):
            warnings.append({"check": "pause_event_readable", "path": str(pause_event)})
    paused_bases = {
        norm(str(VP / matrix_by_run_id[run_id]["output_attempt_directory"]))
        for run_id in paused_run_ids
        if run_id in matrix_by_run_id
    }
    observed_attempts: set[Path] = set()
    checked_manifest_entries = 0
    input_hash_sets: dict[str, set[str]] = {"matrix_mtx_gz": set(), "genes_tsv": set(), "cells_tsv": set(), "provenance_json": set()}

    for row in registry:
        order = int(row["run_order"])
        frozen = matrix_by_order.get(order)
        if frozen is None:
            errors.append({"check": "registry_order_in_matrix", "run_order": order})
            continue
        for field, expected in (("run_id", frozen["run_id"]), ("gene", frozen["gene"]), ("condition", frozen["condition"]), ("donor", frozen["donor"]), ("seed", frozen["seed"])):
            if row[field] != expected:
                errors.append({"check": "registry_matches_frozen_matrix", "run_order": order, "field": field, "expected": expected, "observed": row[field]})
        attempt = norm(row["output_attempt"])
        expected_attempt = norm(str(VP / frozen["output_attempt_directory"]))
        if attempt != expected_attempt:
            errors.append({"check": "output_attempt_matches_frozen_matrix", "run_order": order, "expected": str(expected_attempt), "observed": str(attempt)})
        observed_attempts.add(attempt)
        missing = [name for name in REQUIRED if not (attempt / name).is_file()]
        if missing:
            errors.append({"check": "required_artifacts", "run_order": order, "missing": missing})
            continue

        params = json.loads((attempt / "run_parameters.json").read_text(encoding="utf-8"))
        expected_matrix = norm(str(VP / frozen["matrix_directory"]))
        parameter_checks = {
            "target": frozen["gene"],
            "seed": int(frozen["seed"]),
            "n_cores": 4,
            "package": "scTenifoldKnk",
            "package_version": "1.0.3",
            "matrix_dir": str(expected_matrix).replace("\\", "/"),
        }
        for key, expected in parameter_checks.items():
            observed = params.get(key)
            if key == "matrix_dir":
                observed = str(norm(str(observed))).replace("\\", "/")
            if observed != expected:
                errors.append({"check": "frozen_parameter", "run_order": order, "field": key, "expected": expected, "observed": observed})
        if params.get("frozen_parameters", {}).get("nCores") != 4:
            errors.append({"check": "frozen_parameter", "run_order": order, "field": "frozen_parameters.nCores", "expected": 4, "observed": params.get("frozen_parameters", {}).get("nCores")})
        raw_gate = params.get("gate")
        raw_purpose = params.get("purpose")
        if raw_gate != "VP-G04" or raw_purpose != "Main frozen virtual-KO; biological interpretation prohibited.":
            erratum_entry = erratum_entries.get(frozen["run_id"])
            corrected = (
                raw_gate == "VP-G03"
                and raw_purpose == "Technical pilot only; biological interpretation prohibited."
                and erratum_entry is not None
                and erratum_entry.get("historical_run_parameters_sha256") == sha256(attempt / "run_parameters.json")
                and erratum_entry.get("corrected_metadata") == {"gate": "VP-G04", "purpose": "Main frozen virtual-KO; biological interpretation prohibited."}
            )
            if corrected:
                metadata_erratum_verified.append({"run_order": order, "run_id": frozen["run_id"]})
            else:
                metadata_holds.append({"run_order": order, "run_id": frozen["run_id"], "gate": raw_gate, "purpose": raw_purpose})

        for key, filename in (("matrix_mtx_gz", "matrix.mtx.gz"), ("genes_tsv", "genes.tsv"), ("cells_tsv", "cells.tsv"), ("provenance_json", "provenance.json")):
            observed_hash = params.get("input_sha256", {}).get(key)
            actual_hash = sha256(expected_matrix / filename)
            input_hash_sets[key].add(actual_hash)
            if observed_hash != actual_hash:
                errors.append({"check": "input_hash", "run_order": order, "field": key, "expected": actual_hash, "observed": observed_hash})

        for manifest_row in read_table(attempt / "artifact_manifest.sha256.tsv", "\t"):
            artifact = attempt / manifest_row["file"]
            if not artifact.is_file() or int(manifest_row["bytes"]) != artifact.stat().st_size or manifest_row["sha256"] != sha256(artifact):
                errors.append({"check": "artifact_manifest", "run_order": order, "file": manifest_row["file"]})
            checked_manifest_entries += 1

    actual_attempts = {path.resolve() for path in (VP / "05_runs" / "main").glob("**/attempt_*") if path.is_dir()}
    retained_paused_temporary_attempts = []
    extras = []
    for path in sorted(actual_attempts - observed_attempts):
        if any(path.parent == base.parent and path.name.startswith(base.name + ".tmp_") for base in paused_bases):
            retained_paused_temporary_attempts.append(str(path))
        else:
            extras.append(str(path))
    if extras:
        errors.append({"check": "no_unregistered_biological_attempts", "extra_attempt_paths": extras})

    status = "PASS" if not errors and not metadata_holds else "HOLD_METADATA_PROVENANCE" if not errors else "FAIL_TECHNICAL"
    report = {
        "gate": "VP-G04",
        "audit_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "technical_execution": "PASS" if not errors else "FAIL",
        "metadata_provenance": "HOLD" if metadata_holds else "PASS",
        "frozen_tasks_expected": 45,
        "registry_rows": len(registry),
        "manifest_entries_verified": checked_manifest_entries,
        "unique_input_hash_counts": {key: len(value) for key, value in input_hash_sets.items()},
        "metadata_holds": metadata_holds,
        "metadata_erratum_verified": metadata_erratum_verified,
        "retained_paused_temporary_attempts": retained_paused_temporary_attempts,
        "errors": errors,
        "warnings": warnings,
        "boundary": "Technical execution and artifact integrity only; no biological interpretation or causal conclusion was performed.",
    }
    VALIDATION.mkdir(parents=True, exist_ok=True)
    output = VALIDATION / ("VP_G04_main_runs_COMPLETE.json" if status == "PASS" else "VP_G04_main_runs_HOLD_METADATA_PROVENANCE.json")
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "output": str(output), "errors": len(errors), "metadata_holds": len(metadata_holds)}, ensure_ascii=False))
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
