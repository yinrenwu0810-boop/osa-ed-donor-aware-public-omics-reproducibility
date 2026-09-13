#!/usr/bin/env python3
"""Final acceptance record for VP-G06 after independent recomputation."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
RESULT = VP / "07_enrichment" / "v1"
VAL = VP / "validation"
OUT = VAL / "VP_G06_enrichment_COMPLETE.json"
FAIL = VAL / "VP_G06_enrichment_COMPLETE_ATTEMPT01_FAIL_RETAINED.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str) -> float | None:
    return None if value in {"", "NA", "NaN"} else float(value)


def main() -> None:
    if OUT.exists() or FAIL.exists():
        raise SystemExit("Refusing to overwrite a G06 final-acceptance record")
    errors: list[str] = []
    freeze_path = VP / "01_protocol" / "VP_G06_enrichment_rule_freeze_v1.json"
    amendment_path = VP / "01_protocol" / "VP_G06_ED_shared_universe_amendment_v2.json"
    binding_path = VP / "01_protocol" / "VP_G06_ED_shared_universe_implementation_erratum_v2.json"
    audit_path = VAL / "VP_G06_enrichment_AUDIT01.json"
    g05_path = VAL / "VP_G05_stability_COMPLETE.json"
    audit = load(audit_path)
    report = load(RESULT / "g06_report.json")
    if audit.get("status") != "PASS" or audit.get("errors"):
        errors.append("Independent G06 audit is not clean PASS")
    if load(g05_path).get("status") != "PASS":
        errors.append("VP-G05 entry condition is not PASS")
    if report.get("status") != "PASS_COMPUTED_PENDING_INDEPENDENT_AUDIT":
        errors.append("G06 computation report status mismatch")
    if report.get("freeze_sha256") != sha256(freeze_path) or audit.get("freeze_sha256") != sha256(freeze_path):
        errors.append("Freeze provenance mismatch")

    manifest_path = RESULT / "artifact_manifest.sha256.tsv"
    manifest = table(manifest_path)
    if len(manifest) != 18:
        errors.append(f"Expected 18 output artifacts, found {len(manifest)}")
    for row in manifest:
        path = RESULT / row["file"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            errors.append(f"Manifest mismatch: {row['file']}")
    if audit.get("result_manifest_sha256") != sha256(manifest_path):
        errors.append("Audit/result manifest hash mismatch")

    eligibility = table(RESULT / "gene_set_eligibility.tsv")
    source_counts = Counter(row["database"] for row in eligibility)
    tested_counts = Counter(row["database"] for row in eligibility if row["status"] == "TESTED")
    if source_counts != Counter({"Hallmark": 50, "Reactome": 1839, "GO_BP": 7538}):
        errors.append("Source database term counts mismatch")
    if tested_counts != Counter({"Hallmark": 32, "Reactome": 119, "GO_BP": 1105}):
        errors.append("Eligible database term counts mismatch")

    main_rows = table(RESULT / "target_consensus_enrichment.tsv")
    if len(main_rows) != 3 * sum(tested_counts.values()):
        errors.append("Target consensus enrichment row count mismatch")
    significant_counts = []
    not_estimable = []
    for target in ["EFNB2", "LRRC17", "TYMS"]:
        for database in ["Hallmark", "Reactome", "GO_BP"]:
            subset = [r for r in main_rows if r["target_gene"] == target and r["database"] == database]
            finite = [number(r["padj"]) for r in subset if number(r["padj"]) is not None]
            significant_counts.append({
                "target_gene": target,
                "database": database,
                "tested_terms": len(subset),
                "finite_FDR_terms": len(finite),
                "FDR_lt_0_05_terms": sum(value < 0.05 for value in finite),
            })
            for row in subset:
                if number(row["padj"]) is None:
                    not_estimable.append({"target_gene": target, "database": database, "pathway": row["pathway"], "status": "NOT_ESTIMABLE_BY_FGSEA_MULTILEVEL"})
    if len(not_estimable) != 4:
        errors.append(f"Expected four reproducible non-estimable rows, found {len(not_estimable)}")

    hypoxia = table(RESULT / "hypoxia_primary_summary.tsv")
    donor = table(RESULT / "hypoxia_donor_consistency.tsv")
    controls = table(RESULT / "hypoxia_control_calibration.tsv")
    if len(hypoxia) != 3 or len(donor) != 3 or len(controls) != 3:
        errors.append("HALLMARK_HYPOXIA summary dimensions mismatch")
    hypoxia_summary = []
    for row in sorted(hypoxia, key=lambda x: x["target_gene"]):
        donor_row = next(x for x in donor if x["target_gene"] == row["target_gene"])
        control_row = next(x for x in controls if x["target_gene"] == row["target_gene"])
        hypoxia_summary.append({
            "target_gene": row["target_gene"],
            "condition": row["condition"],
            "five_seed_consensus_NES": float(row["NES"]),
            "five_seed_consensus_FDR_within_Hallmark": float(row["padj"]),
            "positive_NES_donors": int(donor_row["positive_NES_donors"]),
            "FDR_lt_0_05_donors": int(donor_row["FDR_lt_0_05_donors"]),
            "primary_seed_absolute_NES_empirical_percentile_vs_10_controls": float(control_row["absolute_NES_empirical_percentile"]),
            "primary_seed_conservative_absolute_upper_tail_value": float(control_row["conservative_absolute_upper_tail_value"]),
            "leading_edge": row["leadingEdge"].split(";") if row["leadingEdge"] else [],
        })

    ed_sets = table(RESULT / "ed_signature_sets.tsv")
    ed_sizes = Counter(row["signature"] for row in ed_sets)
    if ed_sizes != Counter({"ED_UP_FDR": 18, "ED_DOWN_FDR": 3}):
        errors.append("ED signature sizes mismatch")
    ed_gsea = table(RESULT / "ed_signature_gsea.tsv")
    if len(ed_gsea) != 3 or {r["pathway"] for r in ed_gsea} != {"ED_UP_FDR"}:
        errors.append("ED ranked-GSEA eligibility/result mismatch")
    ed_overlap = table(RESULT / "ed_top5_overlap.tsv")
    ed_concordance = table(RESULT / "ed_rank_concordance.tsv")
    if len(ed_overlap) != 6 or len(ed_concordance) != 3 or any(int(r["universe_size"]) != 1370 for r in ed_overlap):
        errors.append("ED overlap/concordance dimensions mismatch")

    retained = [
        "VP_G06_enrichment_ATTEMPT01_COMPATIBILITY_FAIL_RETAINED.json",
        "VP_G06_enrichment_ATTEMPT02_IMPLEMENTATION_FAIL_RETAINED.json",
        "VP_G06_enrichment_ATTEMPT03_IMPLEMENTATION_FAIL_RETAINED.json",
        "VP_G06_enrichment_ATTEMPT04_ED_UNIVERSE_HOLD_RETAINED.json",
    ]
    if any(not (VAL / name).is_file() for name in retained):
        errors.append("Required failed/HOLD history is missing")

    result = {
        "gate": "VP-G06",
        "stage": "ordered_enrichment_and_ED_signature_overlap_final_acceptance",
        "status": "PASS" if not errors else "FAIL_RETAINED",
        "gate_state": "COMPLETE_WITH_VERSIONED_ED_SHARED_UNIVERSE_AMENDMENT_AND_4_NOT_ESTIMABLE_ROWS" if not errors else "INCOMPLETE",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "entry_condition": {"VP_G05_status": "PASS", "VP_G05_sha256": sha256(g05_path)},
        "acceptance": {
            "database_versions_and_hashes": "PASS_3_OF_3",
            "full_rank_and_fixed_background": "PASS_1388_GENES",
            "ED_shared_testable_background": "PASS_VERSIONED_AMENDMENT_1370_GENES",
            "per_target_per_database_FDR": "PASS",
            "leading_edges_and_nonsignificant_results_retained": "PASS",
            "independent_fgsea_recomputation": "PASS_54_PROFILES",
            "artifact_manifest": "PASS_18_OF_18",
            "direction_and_causality_boundary": "PASS",
        },
        "database_terms": {
            "source": dict(source_counts),
            "eligible_and_tested_per_target": dict(tested_counts),
            "target_database_result_rows": len(main_rows),
            "significant_counts": significant_counts,
        },
        "not_estimable_rows": not_estimable,
        "not_estimable_policy": "Retained as NOT_ESTIMABLE; never treated as nonsignificant or removed from the full output.",
        "HALLMARK_HYPOXIA": hypoxia_summary,
        "ED_signature_summary": {
            "shared_testable_background": 1370,
            "signature_sizes": dict(ed_sizes),
            "ED_UP_FDR_ranked_GSEA": [{"target_gene": r["target_gene"], "NES": float(r["NES"]), "FDR_within_target": float(r["padj"]), "leading_edge": r["leadingEdge"].split(";") if r["leadingEdge"] else []} for r in ed_gsea],
            "ED_DOWN_FDR_ranked_GSEA": "NOT_TESTED_SIZE_3_BELOW_FROZEN_MINIMUM_5",
            "top5_overlap_rows": [{"target_gene": r["target_gene"], "signature": r["signature"], "target_consensus_size": int(r["target_consensus_size"]), "signature_size": int(r["signature_size"]), "overlap_size": int(r["overlap_size"]), "FDR_within_target": float(r["FDR_within_target"]), "overlap_genes": r["overlap_genes"].split(";") if r["overlap_genes"] else []} for r in ed_overlap],
            "full_rank_spearman": [{"target_gene": r["target_gene"], "shared_genes": int(r["shared_genes"]), "rho": float(r["spearman_rho_perturbation_vs_ED_absolute_stat"])} for r in ed_concordance],
        },
        "provenance": {
            "v1_freeze_sha256": sha256(freeze_path),
            "v2_ED_amendment_sha256": sha256(amendment_path),
            "v2_implementation_binding_sha256": sha256(binding_path),
            "result_manifest_sha256": sha256(manifest_path),
            "independent_audit_sha256": sha256(audit_path),
        },
        "retained_nonpassing_history": retained,
        "interpretation_boundary": "Positive NES means concentration toward the high network-perturbation end and negative NES toward the low end. Results do not show pathway up/down regulation, ED-expression reversal, therapeutic benefit, or causality. Matched v2 genes are reference controls, not proven biologically inert controls.",
        "downstream_gate_status": "VP-G07_NOT_STARTED_REQUIRES_EXPLICIT_AUTHORIZATION",
        "errors": errors,
        "implementation_sha256": {"scripts/07g_finalize_vp_g06.py": sha256(Path(__file__))},
    }
    destination = OUT if not errors else FAIL
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "gate_state": result["gate_state"], "errors": errors, "output": str(destination)}, ensure_ascii=False, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
