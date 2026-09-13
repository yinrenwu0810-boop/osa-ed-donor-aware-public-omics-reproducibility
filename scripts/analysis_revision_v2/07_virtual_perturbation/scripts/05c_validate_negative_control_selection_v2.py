"""Independent reproduction audit of the sealed VP-G05 v2 control selection."""

from __future__ import annotations

import csv
import hashlib
import json
import random
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
RULE = VP / "01_protocol" / "VP_G05_negative_control_rule_freeze_v2.json"
REFERENCE = VP / "06_consensus" / "negative_control_degree_reference" / "degree_reference.tsv"
UNIVERSE = VP / "04_prepared" / "gene_universe.tsv"
SELECTION = VP / "06_consensus" / "negative_control_selection_v2"
OUT = VP / "validation" / "VP_G05_negative_control_selection_v2_AUDIT01.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite v2 audit: {OUT}")
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    reference, universe = read_tsv(REFERENCE), read_tsv(UNIVERSE)
    selected = read_tsv(SELECTION / "negative_control_selection.tsv")
    excluded = {"TYMS", "EFNB2", "LRRC17"} | {r["gene"] for r in universe if "HALLMARK_HYPOXIA_DETECTABLE_IN_MAIN_FIBROBLASTS" in r["selection_reason"]}
    by_key = {(r["condition"], r["donor"], r["gene"]): r for r in reference}
    genes, rng, errors, summaries = sorted({r["gene"] for r in reference}), random.Random(rule["sampling"]["seed"]), [], []
    for target in rule["target_definitions"]:
        focal, condition, donors = target["target_gene"], target["condition"], target["donors"]
        eligible = []
        for gene in genes:
            if gene in excluded: continue
            try:
                valid = all(abs(int(by_key[(condition, donor, focal)]["detection_decile"]) - int(by_key[(condition, donor, gene)]["detection_decile"])) <= 1 and abs(int(by_key[(condition, donor, focal)]["outdegree_decile"]) - int(by_key[(condition, donor, gene)]["outdegree_decile"])) <= 1 for donor in donors)
            except KeyError:
                valid = False
            if valid: eligible.append(gene)
        expected = rng.sample(eligible, rule["controls_per_target"])
        observed_rows = [r for r in selected if r["target_gene"] == focal and r["condition"] == condition]
        observed = [r["control_gene"] for r in sorted(observed_rows, key=lambda r: (int(r["selection_order"]), r["donor"]))[::len(donors)]]
        if observed != expected: errors.append(f"{focal}: fixed-seed selection mismatch")
        if len(observed_rows) != 30: errors.append(f"{focal}: expected 30 donor-expanded rows, found {len(observed_rows)}")
        summaries.append({"target_gene": focal, "eligible_pool_size": len(eligible), "selected_controls": observed})
    hash_match = {rel: sha256(VP / rel) == digest for rel, digest in rule["implementation_sha256"].items()}
    if not all(hash_match.values()): errors.append("frozen implementation hash mismatch")
    report = {"gate": "VP-G05", "stage": "matched_negative_control_selection_v2_independent_audit", "status": "PASS" if not errors else "FAIL_RETAINED", "created_at_utc": datetime.now(timezone.utc).isoformat(), "rule_sha256": sha256(RULE), "selection_manifest_sha256": sha256(SELECTION / "artifact_manifest.sha256.tsv"), "implementation_hash_match": hash_match, "target_summary": summaries, "errors": errors, "boundary": "Selection reproducibility audit only; no control KO or pathway result inspected."}
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "errors": len(errors)}, ensure_ascii=False))
    if errors: raise SystemExit(2)


if __name__ == "__main__":
    main()
