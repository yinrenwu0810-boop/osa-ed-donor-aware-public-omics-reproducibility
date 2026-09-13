#!/usr/bin/env python3
"""Read-only VP-ES01-00 integrity gate.

Creates one append-only preflight record.  It never modifies sealed VP-G05--G08
artifacts, virtual-KO outputs, or the formal v1 enrichment directory.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP_ROOT = PROJECT_ROOT / "revision_v2" / "07_virtual_perturbation"
ATTEMPT_DIR = VP_ROOT / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_01"

ANCHORS = {
    "validation/VP_G05_stability_COMPLETE.json": "8af733e965b5486994fc0829febc14d57f5abcb28295ad1ba851e47787574fa7",
    "validation/VP_G06_enrichment_COMPLETE.json": "67740aa1ab08bbd067ca9c85aead06cbbccce84db2e3e0c3308f16e6dfe2a75e",
    "validation/VP_G06_enrichment_AUDIT01.json": "e88929bbb3ea6abb79a8a19f19055921d7060b432c3461d3c8d719b97c069177",
    "01_protocol/VP_G06_enrichment_rule_freeze_v1.json": "38f1f81b4215002ee0f74154b24813e3ee5f704dbda21e6bfb3d7af1b5c47c12",
    "07_enrichment/v1/artifact_manifest.sha256.tsv": "575782ce9e89d54306184ef0affd1ddec0d5ebd706468151e1648c0faaa21450",
}

TAKEOVER_FILES = [
    "00_manual/VIRTUAL_PERTURBATION_WORKBOOK_zh.md",
    "01_protocol/VP_G06_enrichment_rule_freeze_v1.json",
    "validation/VP_G05_stability_COMPLETE.json",
    "validation/VP_G06_enrichment_COMPLETE.json",
    "validation/VP_G06_enrichment_AUDIT01.json",
    "HANDOFF_20260909_VP_G07_COMPLETE.md",
    "HANDOFF_20260909_VP_G08_PREIMPLEMENTATION_HOLD.md",
    "PLAN_20260910_EXPLORATORY_HYPOXIA_SENSITIVITY.md",
]

MANIFEST_FILES = TAKEOVER_FILES + [
    "05_runs/run_registry.tsv",
    "06_consensus/negative_control_selection_v2/negative_control_selection.tsv",
    "04_prepared/gene_universe.tsv",
    "07_enrichment/v1/hypoxia_primary_summary.tsv",
    "07_enrichment/v1/hypoxia_donor_enrichment.tsv",
    "07_enrichment/v1/hypoxia_donor_consistency.tsv",
    "07_enrichment/v1/hypoxia_matched_control_enrichment.tsv",
    "07_enrichment/v1/hypoxia_control_calibration.tsv",
    "07_enrichment/v1/target_consensus_rank_inputs.tsv",
    "07_enrichment/v1/matched_profile_rank_inputs.tsv",
    "07_enrichment/v1/gene_set_eligibility.tsv",
]

EXPECTED_SEEDS = {"2026082701", "2026082702", "2026082703", "2026082704", "2026082705"}
EXPECTED_GROUPS = {
    ("TYMS", "normal", donor) for donor in ("Normal_1", "Normal_2", "Normal_3")
} | {
    (gene, "organic_ED_nonDM", donor)
    for gene in ("EFNB2", "LRRC17")
    for donor in ("non-DM_1", "non-DM_2", "non-DM_3")
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_entry(relative_path: str) -> dict[str, object]:
    path = VP_ROOT / relative_path
    entry: dict[str, object] = {"path": relative_path, "exists": path.is_file()}
    if path.is_file():
        entry.update({"bytes": path.stat().st_size, "sha256": sha256(path)})
    return entry


def read_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return list(csv.DictReader(handle, dialect=dialect))


def validate_main_runs(errors: list[str]) -> dict[str, object]:
    registry_path = VP_ROOT / "05_runs" / "run_registry.tsv"
    if not registry_path.is_file():
        errors.append("missing run_registry.tsv")
        return {"row_count": None}
    try:
        rows = read_registry(registry_path)
    except (csv.Error, UnicodeError) as exc:
        errors.append(f"run_registry.tsv unreadable: {exc}")
        return {"row_count": None}
    required = {"run_id", "gene", "condition", "donor", "seed", "status", "output_attempt"}
    fields = set(rows[0]) if rows else set()
    missing_fields = sorted(required - fields)
    if missing_fields:
        errors.append(f"run_registry.tsv missing columns: {missing_fields}")
    if len(rows) != 45:
        errors.append(f"main registry row count is {len(rows)}, expected 45")
    if len({row.get("run_id") for row in rows}) != len(rows):
        errors.append("duplicate run_id in main registry")
    compound = [(row.get("gene"), row.get("condition"), row.get("donor"), row.get("seed")) for row in rows]
    if len(set(compound)) != len(compound):
        errors.append("duplicate gene/condition/donor/seed combination in main registry")
    nonpass = [row.get("run_id") for row in rows if row.get("status") != "PASS_TECHNICAL"]
    if nonpass:
        errors.append(f"non-PASS_TECHNICAL main runs: {nonpass}")
    observed_groups = {(row.get("gene"), row.get("condition"), row.get("donor")) for row in rows}
    if observed_groups != EXPECTED_GROUPS:
        errors.append("main registry gene/condition/donor groups differ from frozen 3x3 design")
    group_seed_counts = {}
    for group in EXPECTED_GROUPS:
        seeds = {row.get("seed") for row in rows if (row.get("gene"), row.get("condition"), row.get("donor")) == group}
        group_seed_counts["|".join(group)] = sorted(seeds)
        if seeds != EXPECTED_SEEDS:
            errors.append(f"main run seeds differ from freeze for {'|'.join(group)}: {sorted(seeds)}")
    missing_output, empty_output = [], []
    for row in rows:
        output = Path(row.get("output_attempt", "")) / "differential_regulation.tsv"
        if not output.is_file():
            missing_output.append(row.get("run_id"))
        elif output.stat().st_size == 0:
            empty_output.append(row.get("run_id"))
    if missing_output:
        errors.append(f"missing main differential_regulation.tsv: {missing_output}")
    if empty_output:
        errors.append(f"empty main differential_regulation.tsv: {empty_output}")
    return {
        "row_count": len(rows),
        "pass_technical_count": len(rows) - len(nonpass),
        "unique_run_ids": len({row.get("run_id") for row in rows}),
        "unique_gene_condition_donor_seed": len(set(compound)),
        "expected_group_count": len(EXPECTED_GROUPS),
        "group_seed_sets": group_seed_counts,
        "missing_output_run_ids": missing_output,
        "empty_output_run_ids": empty_output,
    }


def validate_control_runs(errors: list[str]) -> dict[str, object]:
    selection_path = VP_ROOT / "06_consensus" / "negative_control_selection_v2" / "negative_control_selection.tsv"
    if not selection_path.is_file():
        errors.append("missing frozen v2 negative-control selection")
        return {"selection_row_count": None}
    with selection_path.open("r", encoding="utf-8-sig", newline="") as handle:
        selection = list(csv.DictReader(handle, delimiter="\t"))
    needed = {"target_gene", "condition", "control_gene", "donor"}
    fields = set(selection[0]) if selection else set()
    if needed - fields:
        errors.append(f"negative-control selection missing columns: {sorted(needed - fields)}")
    expected = {(r["target_gene"], r["condition"], r["donor"], r["control_gene"]) for r in selection}
    if len(selection) != 90 or len(expected) != 90:
        errors.append(f"frozen negative-control selection has {len(selection)} rows and {len(expected)} unique rows; expected 90")
    actual = set()
    missing, empty = [], []
    control_root = VP_ROOT / "05_runs" / "control_batches_v2"
    for path in control_root.glob("*/*/*/attempt_01/controls/*/attempt_01/differential_regulation.tsv"):
        parts = path.relative_to(control_root).parts
        condition, donor, seed, attempt, controls, control_gene, control_attempt, filename = parts
        key = (control_gene, condition, donor, control_gene)
        if seed != "2026082701" or attempt != "attempt_01" or controls != "controls" or control_attempt != "attempt_01":
            errors.append(f"unexpected frozen control path: {path}")
        if path.stat().st_size == 0:
            empty.append(str(path))
        for target_gene, sel_condition, sel_donor, sel_control in expected:
            if (sel_condition, sel_donor, sel_control) == (condition, donor, control_gene):
                actual.add((target_gene, sel_condition, sel_donor, sel_control))
    missing = sorted(expected - actual)
    extra_paths = 0
    for path in control_root.glob("*/*/*/attempt_01/controls/*/attempt_01/differential_regulation.tsv"):
        rel = path.relative_to(control_root).parts
        if not any((r["condition"], r["donor"], r["control_gene"]) == (rel[0], rel[1], rel[5]) for r in selection):
            extra_paths += 1
    if len(actual) != 90:
        errors.append(f"control differential_regulation coverage is {len(actual)}, expected 90")
    if missing:
        errors.append(f"missing frozen control outputs: {missing}")
    if empty:
        errors.append(f"empty frozen control outputs: {empty}")
    if extra_paths:
        errors.append(f"unselected control output paths found: {extra_paths}")
    group_counts = Counter((r["target_gene"], r["condition"], r["donor"]) for r in selection)
    if set(group_counts) != EXPECTED_GROUPS or any(count != 10 for count in group_counts.values()):
        errors.append("negative-control selection does not contain 10 controls per frozen target/donor group")
    return {
        "selection_row_count": len(selection),
        "unique_frozen_control_rows": len(expected),
        "resolved_control_output_count": len(actual),
        "expected_control_output_count": 90,
        "missing_control_outputs": missing,
        "empty_control_outputs": empty,
        "unselected_control_output_path_count": extra_paths,
        "controls_per_target_donor": {"|".join(key): value for key, value in sorted(group_counts.items())},
    }


def main() -> int:
    if Path.cwd().resolve() != PROJECT_ROOT.resolve():
        raise SystemExit(f"Refusing to run outside frozen project root: {Path.cwd()}")
    if ATTEMPT_DIR.exists():
        raise SystemExit(f"Refusing to overwrite existing attempt directory: {ATTEMPT_DIR}")

    errors: list[str] = []
    ATTEMPT_DIR.mkdir(parents=True, exist_ok=False)
    timestamp = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    try:
        takeover = [file_entry(path) for path in TAKEOVER_FILES]
        for entry in takeover:
            if not entry["exists"]:
                errors.append(f"missing takeover file: {entry['path']}")
        anchor_results = []
        for relative_path, expected_sha in ANCHORS.items():
            entry = file_entry(relative_path)
            entry["expected_sha256"] = expected_sha
            entry["matches_expected_sha256"] = entry.get("sha256") == expected_sha
            anchor_results.append(entry)
            if not entry["matches_expected_sha256"]:
                errors.append(f"input hash drift: {relative_path}")
        v1_dir = VP_ROOT / "07_enrichment" / "v1"
        v1_check = {"path": "07_enrichment/v1", "exists": v1_dir.is_dir(), "tmp_directories_read": False}
        if not v1_dir.is_dir():
            errors.append("formal 07_enrichment/v1 directory missing")
        main_runs = validate_main_runs(errors)
        control_runs = validate_control_runs(errors)
        manifest_entries = [file_entry(path) for path in MANIFEST_FILES]
        for entry in manifest_entries:
            if not entry["exists"]:
                errors.append(f"missing declared read-only input: {entry['path']}")
        status = "PASS_INPUTS_READ_ONLY" if not errors else (
            "HOLD_INPUT_DRIFT" if any("hash drift" in item or "missing declared read-only input" in item for item in errors)
            else "HOLD_RUN_REGISTRY_MISMATCH"
        )
        record = {
            "gate": "VP-ES01-00",
            "stage": "recovery_and_input_integrity",
            "status": status,
            "created_at_utc": timestamp,
            "project_root": str(PROJECT_ROOT),
            "working_directory": str(Path.cwd().resolve()),
            "read_only_assertion": "No sealed VP-G05--G08 artifact, virtual-KO run, or formal v1 enrichment file was modified.",
            "takeover_files": takeover,
            "anchor_hash_checks": anchor_results,
            "formal_v1_enrichment_directory": v1_check,
            "main_run_registry": main_runs,
            "matched_reference_runs": control_runs,
            "errors": errors,
        }
        (ATTEMPT_DIR / "00_preflight.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        with (ATTEMPT_DIR / "input_manifest.sha256.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path", "exists", "bytes", "sha256"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest_entries)
        return 0 if status == "PASS_INPUTS_READ_ONLY" else 2
    except Exception as exc:
        hold = {
            "gate": "VP-ES01-00",
            "stage": "recovery_and_input_integrity",
            "status": "HOLD_RUN_REGISTRY_MISMATCH",
            "created_at_utc": timestamp,
            "project_root": str(PROJECT_ROOT),
            "errors": errors + [f"preflight implementation exception: {type(exc).__name__}: {exc}"],
        }
        (ATTEMPT_DIR / "00_preflight.json").write_text(json.dumps(hold, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        return 2


if __name__ == "__main__":
    sys.exit(main())
