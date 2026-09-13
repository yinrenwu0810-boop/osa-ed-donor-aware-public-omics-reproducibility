"""Select VP-G05 matched negative controls from the frozen degree reference."""

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
OUT = VP / "06_consensus" / "negative_control_selection_v1"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite retained selection: {OUT}")
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    if rule.get("status") != "FROZEN_BEFORE_NEGATIVE_CONTROL_SELECTION":
        raise RuntimeError("Negative-control matching rule is not frozen.")
    records, universe = read_tsv(REFERENCE), read_tsv(UNIVERSE)
    hypoxia = {r["gene"] for r in universe if "HALLMARK_HYPOXIA_DETECTABLE_IN_MAIN_FIBROBLASTS" in r["selection_reason"]}
    excluded = hypoxia | {"TYMS", "EFNB2", "LRRC17"}
    by_donor = {(r["condition"], r["donor"], r["gene"]): r for r in records}
    all_genes = sorted({r["gene"] for r in records})
    rng = random.Random(rule["sampling"]["seed"])
    selected_rows: list[dict[str, object]] = []
    eligibility_rows: list[dict[str, object]] = []
    reports: list[dict[str, object]] = []
    for target in rule["target_definitions"]:
        target_gene, condition, donors = target["target_gene"], target["condition"], target["donors"]
        eligible: list[str] = []
        for gene in all_genes:
            reasons: list[str] = []
            if gene in excluded:
                reasons.append("EXCLUDED_PREDECLARED")
            for donor in donors:
                focal = by_donor.get((condition, donor, target_gene))
                candidate = by_donor.get((condition, donor, gene))
                if focal is None or candidate is None:
                    reasons.append(f"MISSING_{donor}")
                    continue
                if int(focal["detection_decile"]) != int(candidate["detection_decile"]):
                    reasons.append(f"DETECTION_DECILE_{donor}")
                if abs(int(focal["outdegree_decile"]) - int(candidate["outdegree_decile"])) > 1:
                    reasons.append(f"OUTDEGREE_DECILE_{donor}")
            is_eligible = not reasons
            eligibility_rows.append({"target_gene": target_gene, "condition": condition, "candidate_gene": gene, "eligible": is_eligible, "reason": ";".join(sorted(set(reasons))) if reasons else "ELIGIBLE"})
            if is_eligible:
                eligible.append(gene)
        if len(eligible) < rule["controls_per_target"]:
            raise RuntimeError(f"HOLD_INSUFFICIENT_MATCHED_CONTROLS: {target_gene} has {len(eligible)} eligible genes.")
        chosen = rng.sample(eligible, rule["controls_per_target"])
        reports.append({"target_gene": target_gene, "condition": condition, "eligible_pool_size": len(eligible), "selected_controls": len(chosen), "status": "PASS_SELECTION"})
        for selection_order, control in enumerate(chosen, start=1):
            for donor in donors:
                focal = by_donor[(condition, donor, target_gene)]
                candidate = by_donor[(condition, donor, control)]
                selected_rows.append({"target_gene": target_gene, "condition": condition, "control_gene": control, "selection_order": selection_order, "donor": donor, "target_detection_rate": focal["detection_rate"], "control_detection_rate": candidate["detection_rate"], "target_detection_decile": focal["detection_decile"], "control_detection_decile": candidate["detection_decile"], "target_mean_outdegree": focal["mean_outdegree_raw_10net"], "control_mean_outdegree": candidate["mean_outdegree_raw_10net"], "target_outdegree_decile": focal["outdegree_decile"], "control_outdegree_decile": candidate["outdegree_decile"]})
    temp = OUT.with_name(OUT.name + ".tmp")
    temp.mkdir(parents=True, exist_ok=False)
    try:
        write_tsv(temp / "negative_control_selection.tsv", selected_rows, list(selected_rows[0]))
        write_tsv(temp / "eligibility_audit.tsv", eligibility_rows, list(eligibility_rows[0]))
        report = {"gate": "VP-G05", "stage": "matched_negative_control_selection", "status": "PASS_SELECTION", "created_at_utc": datetime.now(timezone.utc).isoformat(), "rule_sha256": sha256(RULE), "degree_reference_sha256": sha256(REFERENCE), "selection_report": reports, "boundary": "Selection uses only frozen detection and degree strata; no perturbation ranking or pathway result was inspected."}
        (temp / "selection_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        manifest_rows = [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(temp.iterdir())]
        write_tsv(temp / "artifact_manifest.sha256.tsv", manifest_rows, ["file", "bytes", "sha256"])
        temp.rename(OUT)
        print(json.dumps({"status": report["status"], "output": str(OUT), "targets": len(reports)}, ensure_ascii=False))
    except Exception:
        failed = OUT.with_name(OUT.name + ".attempt_01_FAIL")
        if not failed.exists():
            temp.rename(failed)
        raise


if __name__ == "__main__":
    main()
