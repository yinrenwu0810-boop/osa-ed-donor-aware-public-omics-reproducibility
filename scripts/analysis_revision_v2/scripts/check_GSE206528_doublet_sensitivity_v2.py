"""Fail-closed completion checker for the frozen Gate 03B sensitivity chain."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
BASE = REVISION / "01_work" / "GSE206528_doublet_sensitivity"
CALL = BASE / "attempt_01"
DOWNSTREAM = BASE / "attempt_02"
VALIDATION = REVISION / "04_validation"


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_new(path: Path, value: object) -> None:
    payload = json.dumps(value, ensure_ascii=False, indent=2) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(payload)


def next_record() -> Path:
    VALIDATION.mkdir(parents=True, exist_ok=True)
    number = 1
    while True:
        candidate = VALIDATION / f"gate_03B_doublet_sensitivity_attempt_{number:02d}_PASS.json"
        if not candidate.exists():
            return candidate
        number += 1


def main() -> None:
    required = [
        CALL / "doublet_calling_provenance.json", CALL / "doublet_calls_by_cell.tsv", CALL / "doublet_summary_by_sample.tsv",
        CALL / "doublet_summary_by_historical_cluster.tsv", CALL / "doublet_summary_by_cell_type.tsv",
        DOWNSTREAM / "doublet_exclusion_pseudobulk_provenance.json", DOWNSTREAM / "pseudobulk_comparison_to_03A.json",
        DOWNSTREAM / "pseudobulk_counts.csv", DOWNSTREAM / "pseudobulk_samples.csv", DOWNSTREAM / "celltype_counts_by_donor.tsv",
        DOWNSTREAM / "edgeR_Hallmark" / "analysis_summary.tsv", DOWNSTREAM / "edgeR_Hallmark" / "TYMS_EFNB2_LRRC17_sensitivity.tsv",
        DOWNSTREAM / "edgeR_Hallmark" / "Hallmark_sensitivity_comparison.tsv", DOWNSTREAM / "integration" / "integration_summary.tsv",
        DOWNSTREAM / "integration" / "TYMS_EFNB2_LRRC17_integrated_sensitivity.tsv",
        DOWNSTREAM / "integration" / "IH_ED_Hallmark_celltype_evidence_doublet_sensitivity.tsv",
    ]
    checks = {f"exists::{path.relative_to(REVISION)}": path.is_file() for path in required}
    if not all(checks.values()):
        raise SystemExit("Required Gate 03B output is absent")
    calling = json.loads((CALL / "doublet_calling_provenance.json").read_text(encoding="utf-8"))
    pseudo = json.loads((DOWNSTREAM / "pseudobulk_comparison_to_03A.json").read_text(encoding="utf-8"))
    calls = pd.read_csv(CALL / "doublet_calls_by_cell.tsv", sep="\t")
    samples = pd.read_csv(CALL / "doublet_summary_by_sample.tsv", sep="\t")
    pb_samples = pd.read_csv(DOWNSTREAM / "pseudobulk_samples.csv")
    summary = pd.read_csv(DOWNSTREAM / "edgeR_Hallmark" / "analysis_summary.tsv", sep="\t")
    focus = pd.read_csv(DOWNSTREAM / "integration" / "TYMS_EFNB2_LRRC17_integrated_sensitivity.tsv", sep="\t")
    hallmark = pd.read_csv(DOWNSTREAM / "edgeR_Hallmark" / "Hallmark_sensitivity_comparison.tsv", sep="\t")
    checks.update({
        "calling_status_exact": calling.get("status") == "CALLING_COMPLETE_AWAITING_DOWNSTREAM_SENSITIVITY",
        "eight_samples_exact": len(samples) == 8 and samples["sample_accession"].nunique() == 8,
        "all_64993_cells_called_once": len(calls) == 64993 and calls["cell_id"].nunique() == 64993,
        "pseudobulk_shape_exact": pseudo.get("shape") == [24505, 72],
        "pseudobulk_72_groups_retained": len(pb_samples) == 72 and pseudo.get("unavailable_pseudobulk_groups") == 0,
        "doublet_exclusion_changed_matrix": pseudo.get("value_mismatches", 0) > 0,
        "edger_hallmark_contrasts_complete": len(summary) == 14 and summary[["cell_type", "comparison"]].duplicated().sum() == 0,
        "three_candidate_integration_rows": set(focus["gene_symbol"]) == {"TYMS", "EFNB2", "LRRC17"} and len(focus) == 3,
        "hallmark_comparison_nonempty": len(hallmark) > 0,
    })
    if not all(checks.values()):
        raise SystemExit("Gate 03B completion checks failed")
    checks = {key: bool(value) for key, value in checks.items()}
    old_sig = hallmark["old_padj"].lt(0.05)
    new_sig = hallmark["new_padj"].lt(0.05)
    result = {
        "status": "PASS",
        "gate": "03B_samplewise_doublet_sensitivity",
        "evaluated_at": now(),
        "checks": checks,
        "metrics": {
            "called_cells": int(len(calls)),
            "predicted_doublets": int(calls["predicted_doublet"].sum()),
            "retained_cells": int(pseudo["retained_cells"]),
            "pseudobulk_value_mismatches_vs_03A": int(pseudo["value_mismatches"]),
            "candidate_tiers": focus.set_index("gene_symbol")[["historical_best_tier", "sensitivity_best_tier"]].to_dict(orient="index"),
            "hallmark_significant_old": int(old_sig.sum()),
            "hallmark_significant_doublet_excluded": int(new_sig.sum()),
            "hallmark_significance_status_changed": int((old_sig != new_sig).sum()),
        },
        "gate_03_overall": "PASS",
        "evidence_boundary": "PASS means the frozen result-known doublet sensitivity was completed and fully reported. It does not establish independent biological validation, remove GSE206528 source confounding, or establish causality.",
    }
    record = next_record()
    write_new(record, result)
    latest = VALIDATION / "gate_03B_doublet_sensitivity_latest.json"
    temporary = latest.with_name(f".{latest.name}.{os.getpid()}.partial")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2); handle.write("\n")
    temporary.replace(latest)
    print(json.dumps({**result, "validation_record": str(record.resolve())}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
