"""Independent audit of VP-G05 frozen negative-control selection."""

from __future__ import annotations

import csv
import hashlib
import json
import random
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
RULE = VP / "01_protocol" / "VP_G05_negative_control_rule_freeze_v1.json"
REFERENCE = VP / "06_consensus" / "negative_control_degree_reference" / "degree_reference.tsv"
UNIVERSE = VP / "04_prepared" / "gene_universe.tsv"
SELECTION = VP / "06_consensus" / "negative_control_selection_v1"
OUT = VP / "validation" / "VP_G05_negative_control_selection_AUDIT01.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Audit output already exists: {OUT}")
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    hash_match = {rel: sha256(VP / rel) == digest for rel, digest in rule["implementation_sha256"].items()}
    reference, universe = read_tsv(REFERENCE), read_tsv(UNIVERSE)
    selected = read_tsv(SELECTION / "negative_control_selection.tsv")
    excluded = {"TYMS", "EFNB2", "LRRC17"} | {r["gene"] for r in universe if "HALLMARK_HYPOXIA_DETECTABLE_IN_MAIN_FIBROBLASTS" in r["selection_reason"]}
    by_key = {(r["condition"], r["donor"], r["gene"]): r for r in reference}
    genes, rng, errors, target_summary = sorted({r["gene"] for r in reference}), random.Random(rule["sampling"]["seed"]), [], []
    for spec in rule["target_definitions"]:
        target, condition, donors = spec["target_gene"], spec["condition"], spec["donors"]
        eligible = []
        for gene in genes:
            if gene in excluded:
                continue
            try:
                valid = all(int(by_key[(condition, donor, target)]["detection_decile"]) == int(by_key[(condition, donor, gene)]["detection_decile"]) and abs(int(by_key[(condition, donor, target)]["outdegree_decile"]) - int(by_key[(condition, donor, gene)]["outdegree_decile"])) <= 1 for donor in donors)
            except KeyError:
                valid = False
            if valid:
                eligible.append(gene)
        expected = rng.sample(eligible, rule["controls_per_target"])
        rows = [r for r in selected if r["target_gene"] == target and r["condition"] == condition]
        actual = [r["control_gene"] for r in sorted(rows, key=lambda r: (int(r["selection_order"]), r["donor"]))[::len(donors)]]
        if len(rows) != rule["controls_per_target"] * len(donors):
            errors.append(f"{target}: unexpected selected row count {len(rows)}")
        if actual != expected:
            errors.append(f"{target}: selected controls do not reproduce fixed-seed draw")
        target_summary.append({"target_gene": target, "eligible_pool_size": len(eligible), "selected_control_count": len(actual), "selected_controls_sha256": hashlib.sha256("\n".join(actual).encode()).hexdigest()})
    report = {"gate": "VP-G05", "stage": "negative_control_selection_independent_audit", "status": "PASS" if not errors and all(hash_match.values()) else "FAIL_RETAINED", "created_at_utc": datetime.now(timezone.utc).isoformat(), "rule_sha256": sha256(RULE), "degree_reference_sha256": sha256(REFERENCE), "selection_manifest_sha256": sha256(SELECTION / "artifact_manifest.sha256.tsv"), "implementation_hash_match": hash_match, "target_summary": target_summary, "errors": errors, "boundary": "This audit validates rule-concordant selection only; it does not inspect KO perturbation or pathway results."}
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "errors": len(errors)}, ensure_ascii=False))
    if report["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
