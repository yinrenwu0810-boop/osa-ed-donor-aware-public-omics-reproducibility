"""Apply the frozen VP-G05 v2 matched-control rule without KO-result access."""

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
OUT = VP / "06_consensus" / "negative_control_selection_v2"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, records: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite sealed v2 selection: {OUT}")
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    if rule.get("status") != "FROZEN_AMENDMENT_BEFORE_V2_SELECTION":
        raise RuntimeError("V2 control rule is not frozen.")
    reference, universe = read_tsv(REFERENCE), read_tsv(UNIVERSE)
    excluded = {"TYMS", "EFNB2", "LRRC17"} | {r["gene"] for r in universe if "HALLMARK_HYPOXIA_DETECTABLE_IN_MAIN_FIBROBLASTS" in r["selection_reason"]}
    by_key = {(r["condition"], r["donor"], r["gene"]): r for r in reference}
    genes = sorted({r["gene"] for r in reference})
    rng = random.Random(rule["sampling"]["seed"])
    selected_rows, eligibility_rows, summaries = [], [], []
    for target in rule["target_definitions"]:
        focal, condition, donors = target["target_gene"], target["condition"], target["donors"]
        eligible = []
        for gene in genes:
            reasons = []
            if gene in excluded:
                reasons.append("EXCLUDED_PREDECLARED")
            for donor in donors:
                target_row, candidate_row = by_key.get((condition, donor, focal)), by_key.get((condition, donor, gene))
                if target_row is None or candidate_row is None:
                    reasons.append(f"MISSING_{donor}"); continue
                if abs(int(target_row["detection_decile"]) - int(candidate_row["detection_decile"])) > 1:
                    reasons.append(f"DETECTION_DECILE_{donor}")
                if abs(int(target_row["outdegree_decile"]) - int(candidate_row["outdegree_decile"])) > 1:
                    reasons.append(f"OUTDEGREE_DECILE_{donor}")
            is_eligible = not reasons
            eligibility_rows.append({"target_gene": focal, "condition": condition, "candidate_gene": gene, "eligible": is_eligible, "reason": ";".join(sorted(set(reasons))) if reasons else "ELIGIBLE"})
            if is_eligible:
                eligible.append(gene)
        if len(eligible) < rule["controls_per_target"]:
            raise RuntimeError(f"V2 HOLD: {focal} has {len(eligible)} eligible controls.")
        chosen = rng.sample(eligible, rule["controls_per_target"])
        summaries.append({"target_gene": focal, "condition": condition, "eligible_pool_size": len(eligible), "selected_controls": len(chosen), "status": "PASS_SELECTION"})
        for order, control in enumerate(chosen, 1):
            for donor in donors:
                target_row, candidate_row = by_key[(condition, donor, focal)], by_key[(condition, donor, control)]
                selected_rows.append({"target_gene": focal, "condition": condition, "control_gene": control, "selection_order": order, "donor": donor, "target_detection_decile": target_row["detection_decile"], "control_detection_decile": candidate_row["detection_decile"], "target_outdegree_decile": target_row["outdegree_decile"], "control_outdegree_decile": candidate_row["outdegree_decile"], "target_detection_rate": target_row["detection_rate"], "control_detection_rate": candidate_row["detection_rate"], "target_mean_outdegree": target_row["mean_outdegree_raw_10net"], "control_mean_outdegree": candidate_row["mean_outdegree_raw_10net"]})
    temp = OUT.with_name(OUT.name + ".tmp")
    temp.mkdir(parents=True, exist_ok=False)
    try:
        write_tsv(temp / "negative_control_selection.tsv", selected_rows, list(selected_rows[0]))
        write_tsv(temp / "eligibility_audit.tsv", eligibility_rows, list(eligibility_rows[0]))
        report = {"gate": "VP-G05", "stage": "matched_negative_control_selection_v2", "status": "PASS_SELECTION", "created_at_utc": datetime.now(timezone.utc).isoformat(), "rule_sha256": sha256(RULE), "degree_reference_sha256": sha256(REFERENCE), "selection_report": summaries, "boundary": "V2 feasibility amendment applied without perturbation-ranking or pathway-result access."}
        (temp / "selection_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        manifest = [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(temp.iterdir())]
        write_tsv(temp / "artifact_manifest.sha256.tsv", manifest, ["file", "bytes", "sha256"])
        temp.rename(OUT)
        print(json.dumps({"status": report["status"], "targets": len(summaries)}, ensure_ascii=False))
    except Exception:
        failed = OUT.with_name(OUT.name + ".attempt_01_FAIL")
        if not failed.exists(): temp.rename(failed)
        raise


if __name__ == "__main__":
    main()
