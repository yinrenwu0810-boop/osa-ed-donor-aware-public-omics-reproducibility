"""Independently audit the VP-G05 outdegree reference before control selection."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
RULE = VP / "01_protocol" / "VP_G05_negative_control_rule_freeze_v1.json"
REF = VP / "06_consensus" / "negative_control_degree_reference"
OUT = VP / "validation" / "VP_G05_degree_reference_AUDIT01.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite prior audit: {OUT}")
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    metadata = json.loads((REF / "metadata.json").read_text(encoding="utf-8"))
    manifest = rows(REF / "artifact_manifest.sha256.tsv")
    degree_rows, qc_rows = rows(REF / "degree_reference.tsv"), rows(REF / "network_qc.tsv")
    errors: list[str] = []
    for entry in manifest:
        path = REF / entry["file"]
        if not path.is_file() or path.stat().st_size != int(entry["bytes"]) or sha256(path) != entry["sha256"]:
            errors.append(f"manifest mismatch: {entry['file']}")
    donor_keys = {(r["condition"], r["donor"]) for r in degree_rows}
    per_donor = {key: sum((r["condition"], r["donor"]) == key for r in degree_rows) for key in donor_keys}
    if len(donor_keys) != 6 or len(set(per_donor.values())) != 1:
        errors.append("degree reference does not contain six donor matrices with one common gene count")
    if len(qc_rows) != 60 or any(sum((r["condition"], r["donor"]) == key for r in qc_rows) != 10 for key in donor_keys):
        errors.append("network QC does not contain exactly 10 networks per donor")
    if any(int(r["detection_decile"]) not in range(1, 11) or int(r["outdegree_decile"]) not in range(1, 11) for r in degree_rows):
        errors.append("out-of-range matching decile")
    if len({(r["condition"], r["donor"], r["gene"]) for r in degree_rows}) != len(degree_rows):
        errors.append("duplicate degree-reference gene key")
    if metadata.get("status") != "PASS_TECHNICAL":
        errors.append("degree-reference metadata is not technical PASS")
    if metadata.get("rule_sha256") != sha256(RULE):
        errors.append("degree-reference rule hash does not match frozen rule")
    report = {
        "gate": "VP-G05",
        "stage": "negative_control_outdegree_reference_independent_audit",
        "status": "PASS" if not errors else "FAIL_RETAINED",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rule_sha256": sha256(RULE),
        "reference_manifest_sha256": sha256(REF / "artifact_manifest.sha256.tsv"),
        "donor_count": len(donor_keys),
        "genes_per_donor": sorted(set(per_donor.values())),
        "network_qc_row_count": len(qc_rows),
        "errors": errors,
        "boundary": "Technical audit of control-matching inputs only; no perturbation ranking or pathway result is inspected.",
    }
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "errors": len(errors)}, ensure_ascii=False))
    if errors:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
