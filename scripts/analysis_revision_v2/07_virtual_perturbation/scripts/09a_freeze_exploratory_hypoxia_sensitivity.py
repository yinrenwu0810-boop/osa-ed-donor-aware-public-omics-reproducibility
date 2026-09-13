#!/usr/bin/env python3
"""VP-ES01-01 append-only scope freeze; no result tables are read here."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP_ROOT = PROJECT_ROOT / "revision_v2" / "07_virtual_perturbation"
ATTEMPT_DIR = VP_ROOT / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_01"
PREFLIGHT = ATTEMPT_DIR / "00_preflight.json"
SCOPE_FREEZE = ATTEMPT_DIR / "01_scope_freeze.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    if Path.cwd().resolve() != PROJECT_ROOT.resolve():
        raise SystemExit(f"Refusing to run outside frozen project root: {Path.cwd()}")
    if not PREFLIGHT.is_file():
        raise SystemExit(f"Required VP-ES01-00 preflight record is absent: {PREFLIGHT}")
    if SCOPE_FREEZE.exists():
        raise SystemExit(f"Refusing to overwrite existing scope freeze: {SCOPE_FREEZE}")

    preflight = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    if preflight.get("status") != "PASS_INPUTS_READ_ONLY" or preflight.get("errors"):
        raise SystemExit("VP-ES01-00 is not a clean PASS_INPUTS_READ_ONLY record; scope cannot be frozen.")

    script_path = Path(__file__).resolve()
    timestamp = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    scope = {
        "gate": "VP-ES01-01",
        "stage": "exploratory_analysis_scope_freeze",
        "status": "FROZEN_EXPLORATORY_SPECIFICATION",
        "created_at_utc": timestamp,
        "project_root": str(PROJECT_ROOT),
        "preflight_record": {
            "path": str(PREFLIGHT.relative_to(VP_ROOT)).replace("\\", "/"),
            "sha256": sha256(PREFLIGHT),
            "status": preflight["status"],
        },
        "analysis_label": "EXPLORATORY_SENSITIVITY_NOT_CONFIRMATORY",
        "formal_result_protection": {
            "formal_baseline": "VP-G06 HALLMARK_HYPOXIA remains nonsignificant for EFNB2, LRRC17, and TYMS.",
            "rule": "No VP-ES01 result may replace, rename, delete, downgrade, or otherwise alter the VP-G06 formal result.",
        },
        "targets_and_primary_backgrounds": [
            {"target_gene": "EFNB2", "primary_background": "organic_ED_nonDM"},
            {"target_gene": "LRRC17", "primary_background": "organic_ED_nonDM"},
            {"target_gene": "TYMS", "primary_background": "normal"},
        ],
        "units": {
            "biological_independent_unit": "donor",
            "technical_random_unit": "network random seed",
            "technical_unit_boundary": "Seeds estimate network stability and must never be analyzed or described as biological replication.",
        },
        "primary_ranking": {
            "metric": "distance",
            "sort_order": "descending",
            "tie_break": "gene symbol ascending",
            "gene_universe": "fixed 1,388-gene network universe",
        },
        "primary_gene_set": {
            "name": "HALLMARK_HYPOXIA",
            "source": "frozen VP-G06 Hallmark GMT",
            "formal_multiple_testing_family": "retain the original VP-G06 within-Hallmark FDR family; do not shrink it",
            "diagnostic_single_set_p_rule": "If calculated solely for diagnosis, record separately as diagnostic_single_set_p and never mix it with the formal Hallmark FDR.",
        },
        "allowed_specifications_in_required_order": [
            {
                "id": "S0",
                "ranking_and_aggregation": "distance descending; five-seed median standardized rank within donor; three-donor consensus",
                "purpose": "read-only formal-baseline reproduction",
                "evidence_label": "FORMAL_BASELINE",
            },
            {
                "id": "S1",
                "ranking_and_aggregation": "distance descending; five-seed mean standardized rank within donor; three-donor mean rank",
                "purpose": "test dependence on median aggregation",
                "evidence_label": "EXPLORATORY",
            },
            {
                "id": "S2",
                "ranking_and_aggregation": "Z descending; otherwise the S0 aggregation rule",
                "purpose": "test dependence on scTenifoldKnk output metric",
                "evidence_label": "EXPLORATORY",
            },
            {
                "id": "S3",
                "ranking_and_aggregation": "enrich each donor before descriptive aggregation of donor NES",
                "purpose": "test dependence on aggregating rankings before enrichment",
                "evidence_label": "EXPLORATORY",
            },
        ],
        "disallowed_analytic_changes": [
            "No additional specification, donor exclusion, seed selection, direction reversal, pathway-family reduction, control reselection, virtual-KO rerun, external-data retrieval, wet experiment, paid run, or causal endpoint may be added in attempt_20260910_01.",
            "A specification defect requires retaining this freeze and creating a separate attempt_02 freeze before any new calculation.",
        ],
        "planned_output_contract": {
            "02_seed_rank_qc.tsv": ["target_gene", "condition", "donor", "seed", "gene_universe_size", "target_present", "finite_distance", "status"],
            "02_seed_pair_stability.tsv": ["target_gene", "condition", "donor", "seed_a", "seed_b", "spearman_rho", "top200_jaccard"],
            "02_seed_aggregate_reconciliation.tsv": ["target_gene", "condition", "donor", "gene", "reconstructed_rank", "frozen_rank", "comparison_status"],
            "03_target_donor_seed_hypoxia.tsv": ["target_gene", "condition", "donor", "seed", "ES", "NES", "p_value", "FDR_within_Hallmark", "leading_edge", "leading_edge_size", "hallmark_rank"],
            "03_target_donor_seed_hallmark_full.tsv.gz": ["target_gene", "condition", "donor", "seed", "pathway", "ES", "NES", "p_value", "FDR_within_Hallmark", "leading_edge"],
            "03_seed_hypoxia_summary.tsv": ["target_gene", "condition", "donor", "NES_median", "NES_range", "positive_NES_seed_count", "FDR_lt_0_05_seed_count"],
            "04_leave_one_seed_out.tsv": ["target_gene", "condition", "donor", "excluded_seed", "ES", "NES", "p_value", "FDR_within_Hallmark", "leading_edge_jaccard_vs_formal"],
            "04_leave_one_donor_out.tsv": ["target_gene", "condition", "excluded_donor", "ES", "NES", "p_value", "FDR_within_Hallmark", "leading_edge_jaccard_vs_formal"],
            "04_donor_heterogeneity.tsv": ["target_gene", "condition", "donor_NES_median", "donor_NES_min", "donor_NES_max", "donor_NES_MAD", "positive_NES_donor_count", "heterogeneity_description"],
            "05_target_vs_matched_controls.tsv": ["target_gene", "condition", "control_gene", "signed_NES_percentile", "absolute_NES_percentile", "conservative_absolute_upper_tail", "donor_positive_NES_proportion", "seed_NES_stability", "leading_edge_size"],
            "05_leading_edge_stability.tsv": ["comparison_type", "target_gene", "comparison_id_a", "comparison_id_b", "leading_edge_jaccard", "leading_edge_size_a", "leading_edge_size_b"],
            "05_leading_edge_overlap_with_ED.tsv": ["target_gene", "reference_set", "reference_status", "overlap_size", "overlap_genes", "not_available_reason"],
            "06_hypoxia_mapping_qc.tsv": ["gene_set", "original_member_count", "network_universe_member_count", "mapping_fraction", "missing_members", "mapping_issue_status"],
            "06_analysis_multiverse.tsv": ["specification", "evidence_label", "target_gene", "condition", "ES", "NES", "p_value", "FDR_within_Hallmark", "result_status"],
            "06_specification_concordance.tsv": ["target_gene", "specification_a", "specification_b", "NES_difference", "direction_concordant", "leading_edge_jaccard"],
            "Figure_ES01_hypoxia_robustness.png": "Regenerable from machine-readable donor/seed and multiverse tables.",
            "Figure_ES02_leading_edge_stability.png": "Regenerable from machine-readable leading-edge stability tables.",
            "VP_ES01_exploratory_sensitivity_REPORT.md": "Chinese report presenting effect size, direction, raw P, Hallmark FDR, empirical references, uncertainty, and scope boundary together.",
        },
        "predeclared_failure_and_final_statuses": [
            "HOLD_INPUT_DRIFT",
            "HOLD_RUN_REGISTRY_MISMATCH",
            "FAIL_RECONCILIATION",
            "HOLD_POTENTIAL_IMPLEMENTATION_DEFECT",
            "SEED_OR_NETWORK_UNSTABLE",
            "DONOR_SENSITIVE_SIGNAL",
            "ROBUST_POSITIVE_DIRECTION_BUT_NONSPECIFIC",
            "TARGET_SPECIFIC_EXPLORATORY_SIGNAL",
            "NO_CONSISTENT_HYPOXIA_SIGNAL",
        ],
        "prohibited_interpretations": [
            "Sensitivity analysis makes FDR significant.",
            "The three genes regulate the hypoxia pathway.",
            "OSA causes ED through the three genes.",
            "CPAP reverses ED through the three genes.",
            "Virtual KO shows pathway upregulation or downregulation.",
            "More seeds are more biological samples.",
            "Exploratory significance is confirmatory validation.",
        ],
        "implementation": {
            "script": str(script_path.relative_to(VP_ROOT)).replace("\\", "/"),
            "script_sha256": sha256(script_path),
            "result_tables_read_before_freeze": False,
            "scope_change_rule": "This file is immutable after creation. Any correction requires a new append-only attempt directory.",
        },
    }
    SCOPE_FREEZE.write_text(json.dumps(scope, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
