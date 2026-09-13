"""Audit the sealed VP-G05 v2 control task and batch matrices before launch."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
PROTOCOL = VP / "01_protocol" / "control_execution_v2"
FREEZE = PROTOCOL / "VP_G05_control_execution_freeze_v2.json"
SELECTION = VP / "06_consensus" / "negative_control_selection_v2" / "negative_control_selection.tsv"
OUT = VP / "validation" / "VP_G05_control_execution_freeze_v2_AUDIT01.json"


def sha256(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if OUT.exists(): raise RuntimeError(f"Refusing to overwrite audit: {OUT}")
    freeze = json.loads(FREEZE.read_text(encoding="utf-8")); errors = []
    manifest = read_tsv(PROTOCOL / "artifact_manifest.sha256.tsv")
    for row in manifest:
        path = PROTOCOL / row["file"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]: errors.append(f"manifest mismatch: {row['file']}")
    if sha256(VP / "scripts" / "05d_run_control_batch_shared_wt.R") != freeze["runner_sha256"]: errors.append("runner hash mismatch")
    if sha256(VP / "scripts" / "05f_run_control_batches_v2.ps1") != freeze["launcher_sha256"]: errors.append("launcher hash mismatch")
    batches, tasks, selected = read_tsv(PROTOCOL / "control_batch_matrix.tsv"), read_tsv(PROTOCOL / "control_task_matrix.tsv"), read_tsv(SELECTION)
    expected_memberships = {(r["control_gene"], r["condition"], r["donor"], r["target_gene"]) for r in selected}
    actual_memberships = set()
    for row in tasks:
        for target in row["matched_targets"].split(";"): actual_memberships.add((row["control_gene"], row["condition"], row["donor"], target))
    if actual_memberships != expected_memberships: errors.append("task matrix does not reproduce sealed selection memberships")
    if len(batches) != 6 or len(tasks) != 90: errors.append(f"expected 6 batches/90 tasks; found {len(batches)}/{len(tasks)}")
    if [int(r["control_count"]) for r in batches] != [10, 10, 10, 20, 20, 20]: errors.append("unexpected batch control counts")
    for row in batches:
        if not (PROTOCOL / row["controls_file"]).is_file(): errors.append(f"missing controls file: {row['controls_file']}")
        if not (VP / row["benchmark_file"]).is_file(): errors.append(f"missing benchmark: {row['benchmark_file']}")
        if not (VP / row["matrix_directory"]).is_dir(): errors.append(f"missing matrix: {row['matrix_directory']}")
    report = {"gate": "VP-G05", "stage": "control_execution_freeze_v2_independent_audit", "status": "PASS" if not errors else "FAIL_RETAINED", "created_at_utc": datetime.now(timezone.utc).isoformat(), "freeze_sha256": sha256(FREEZE), "protocol_manifest_sha256": sha256(PROTOCOL / "artifact_manifest.sha256.tsv"), "batch_count": len(batches), "unique_control_tasks": len(tasks), "selection_memberships": len(expected_memberships), "errors": errors, "boundary": "Pre-execution structural and provenance audit; no control KO was launched or interpreted."}
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "errors": len(errors)}, ensure_ascii=False))
    if errors: raise SystemExit(2)


if __name__ == "__main__": main()
