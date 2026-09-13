"""Build the sealed, de-duplicated VP-G05 v2 control execution batches."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
RULE = VP / "01_protocol" / "VP_G05_negative_control_rule_freeze_v2.json"
AUDIT = VP / "validation" / "VP_G05_negative_control_selection_v2_AUDIT01.json"
SELECTION = VP / "06_consensus" / "negative_control_selection_v2" / "negative_control_selection.tsv"
REGISTRY = VP / "05_runs" / "run_registry.tsv"
RUNNER = VP / "scripts" / "05d_run_control_batch_shared_wt.R"
LAUNCHER = VP / "scripts" / "05f_run_control_batches_v2.ps1"
OUT = VP / "01_protocol" / "control_execution_v2"


def sha256(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def read_tsv(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle, delimiter=delimiter))
def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def main() -> None:
    if OUT.exists(): raise RuntimeError(f"Refusing to overwrite execution freeze: {OUT}")
    if json.loads(AUDIT.read_text(encoding="utf-8")).get("status") != "PASS": raise RuntimeError("V2 selection audit is not PASS.")
    if not RUNNER.is_file() or not LAUNCHER.is_file(): raise RuntimeError("Control runner or launcher missing before freeze.")
    selected, main_registry = read_tsv(SELECTION), read_tsv(REGISTRY, ",")
    memberships = {(r["control_gene"], r["condition"], r["donor"], r["target_gene"]) for r in selected}
    grouped: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for control, condition, donor, target in memberships: grouped[(condition, donor)][control].add(target)
    temp = OUT.with_name(OUT.name + ".tmp"); temp.mkdir(parents=True, exist_ok=False)
    (temp / "batch_controls").mkdir()
    batch_rows, task_rows = [], []
    definitions = [("normal", f"Normal_{i}", "TYMS") for i in range(1, 4)] + [("organic_ED_nonDM", f"non-DM_{i}", "EFNB2") for i in range(1, 4)]
    task_order = 1
    for batch_order, (condition, donor, benchmark_gene) in enumerate(definitions, 1):
        controls = grouped[(condition, donor)]
        control_rows = [{"control_gene": gene, "matched_targets": ";".join(sorted(targets))} for gene, targets in sorted(controls.items())]
        control_file = temp / "batch_controls" / f"batch_{batch_order:02d}_{condition}_{donor}.tsv"
        write_tsv(control_file, control_rows, ["control_gene", "matched_targets"])
        benchmark = [r for r in main_registry if r["gene"] == benchmark_gene and r["condition"] == condition and r["donor"] == donor and r["seed"] == "2026082701" and r["status"] == "PASS_TECHNICAL"]
        if len(benchmark) != 1: raise RuntimeError(f"Expected one direct benchmark for {condition}/{donor}/{benchmark_gene}")
        batch_rows.append({"batch_order": batch_order, "condition": condition, "donor": donor, "seed": 2026082701, "control_count": len(control_rows), "controls_file": str(control_file.relative_to(temp)).replace("\\", "/"), "matrix_directory": f"04_prepared/fibroblast/{condition}/{donor}", "benchmark_gene": benchmark_gene, "benchmark_file": str((Path(benchmark[0]["output_attempt"]) / "differential_regulation.tsv").relative_to(VP)).replace("\\", "/"), "output_batch": f"05_runs/control_batches_v2/{condition}/{donor}/2026082701/attempt_01"})
        for row in control_rows:
            task_rows.append({"task_order": task_order, "control_gene": row["control_gene"], "matched_targets": row["matched_targets"], "condition": condition, "donor": donor, "seed": 2026082701, "batch_order": batch_order}); task_order += 1
    write_tsv(temp / "control_batch_matrix.tsv", batch_rows, list(batch_rows[0]))
    write_tsv(temp / "control_task_matrix.tsv", task_rows, list(task_rows[0]))
    freeze = {"gate": "VP-G05", "status": "FROZEN_BEFORE_CONTROL_KO", "version": "v2", "created_at_utc": datetime.now(timezone.utc).isoformat(), "rule_sha256": sha256(RULE), "selection_audit_sha256": sha256(AUDIT), "selection_sha256": sha256(SELECTION), "primary_seed": 2026082701, "batch_count": len(batch_rows), "unique_control_tasks": len(task_rows), "execution_design": "Build one deterministic WT network per condition/donor/seed, preserve RNG state immediately before manifold alignment, validate against a retained direct main-KO result, then reuse it for all controls in that batch.", "direct_equivalence_tolerance_absolute": 1e-7, "runner_sha256": sha256(RUNNER), "launcher_sha256": sha256(LAUNCHER), "boundary": "Shared-WT reuse is a computational optimization only; model parameters, donor context, seed, KO operation, manifold alignment, and differential-regulation calculation are unchanged."}
    (temp / "VP_G05_control_execution_freeze_v2.json").write_text(json.dumps(freeze, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    manifest = []
    for path in sorted(p for p in temp.rglob("*") if p.is_file()): manifest.append({"file": str(path.relative_to(temp)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_tsv(temp / "artifact_manifest.sha256.tsv", manifest, ["file", "bytes", "sha256"])
    temp.rename(OUT)
    print(json.dumps({"status": freeze["status"], "batches": len(batch_rows), "tasks": len(task_rows)}, ensure_ascii=False))


if __name__ == "__main__": main()
