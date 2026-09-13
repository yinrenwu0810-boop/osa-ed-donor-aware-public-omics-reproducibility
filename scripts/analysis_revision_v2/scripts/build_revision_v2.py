#!/usr/bin/env python3
"""Build evidence-bounded scientific revision v2 tables without altering v1.

This script consumes only existing v1 analysis products. It does not re-fit the
underlying differential-expression models and it does not query live services.
Its purpose is to correct classification, denominators, and missingness labels
identified during the 2026-08-12 scientific audit.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
PROTOCOL = REVISION / "00_protocol"
RESULTS = REVISION / "02_results"
MANUSCRIPT = REVISION / "03_manuscript"


def write_numbered_attempt(
    validation_dir: Path, prefix: str, status: str, payload: str
) -> Path:
    number = 1
    while True:
        if any(validation_dir.glob(f"{prefix}_{number:02d}_*.json")):
            number += 1
            continue
        target = validation_dir / f"{prefix}_{number:02d}_{status}.json"
        try:
            with target.open("x", encoding="utf-8", newline="\n") as handle:
                handle.write(payload)
            return target
        except FileExistsError:
            number += 1


def atomic_replace_json(path: Path, payload: str) -> None:
    counter = 0
    while True:
        temporary = path.with_name(f".{path.name}.{os.getpid()}.{counter}.tmp")
        try:
            handle = temporary.open("x", encoding="utf-8", newline="\n")
            break
        except FileExistsError:
            counter += 1
    try:
        with handle:
            handle.write(payload)
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


INPUTS = [
    "HANDOFF.md",
    "analysis/R/50_integrate_bulk_evidence.R",
    "analysis/R/51_integrate_IH_ED_celltypes.R",
    "analysis/R/52_prepare_external_candidates.R",
    "analysis/R/53_finalize_candidate_ranking.R",
    "analysis/scripts/71_endothelial_robustness.py",
    "analysis/scripts/query_opentargets_candidates.py",
    "analysis/scripts/query_gwas_catalog_traits.py",
    "analysis/scripts/query_reactome_candidates.py",
    "analysis/results/phase3/endothelial_robustness/robustness_with_lodo.tsv",
    "analysis/results/phase3/endothelial_robustness/leave_one_study_out_meta.tsv",
    "analysis/results/integration/IH_ED_celltype_gene_evidence.tsv",
    "analysis/results/integration/IH_ED_prioritized_genes.tsv",
    "analysis/results/integration/final_candidate_ranking.tsv",
    "analysis/results/external_databases/candidate_genes_for_lookup.tsv",
    "analysis/results/external_databases/opentargets_candidate_evidence.tsv",
    "analysis/results/external_databases/gwas_candidate_overlaps.tsv",
    "analysis/results/external_databases/reactome_candidate_pathways.tsv",
    "analysis/results/external_databases/opentargets_candidate_provenance.json",
    "analysis/results/external_databases/gwas_catalog_provenance.json",
    "analysis/results/external_databases/reactome_candidate_provenance.json",
    "paper_rewriting_output/manuscript_en.src.md",
    "paper_rewriting_output/manuscript_zh.src.md",
    "paper_rewriting_output/results_validation.md",
    "paper_rewriting_output/reviewer_audit.md",
    "paper_rewriting_output/integrity_audit.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(relative: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / relative, sep="\t", low_memory=False)


def bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().map({"true": True, "false": False}).fillna(False)


def freeze_inputs() -> pd.DataFrame:
    rows = []
    for relative in INPUTS:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"Frozen input missing: {relative}")
        rows.append(
            {
                "relative_path": relative.replace("\\", "/"),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest = pd.DataFrame(rows).sort_values("relative_path")
    manifest.to_csv(PROTOCOL / "input_manifest.sha256.tsv", sep="\t", index=False)
    return manifest


def build_lodo_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    full = read_tsv(
        "analysis/results/phase3/endothelial_robustness/robustness_with_lodo.tsv"
    )
    lodo = read_tsv(
        "analysis/results/phase3/endothelial_robustness/leave_one_study_out_meta.tsv"
    )

    fdr_wide = lodo.pivot(
        index="gene_symbol", columns="omitted_study", values="lodo_meta_FDR"
    )
    fdr_wide.columns = [f"lodo_FDR_excluding_{column}" for column in fdr_wide.columns]
    fdr_wide = fdr_wide.reset_index()

    out = full.merge(fdr_wide, on="gene_symbol", how="left", validate="one_to_one")
    fdr_columns = [column for column in out if column.startswith("lodo_FDR_excluding_")]
    if len(fdr_columns) != 3:
        raise RuntimeError(f"Expected three LODO FDR columns, found {len(fdr_columns)}")

    out["full_meta_FDR_pass"] = out["meta_FDR"].lt(0.05)
    out["three_study_direction_concordant"] = bool_series(out["all_same_direction"])
    out["directional_stability_screen"] = (
        out["full_meta_FDR_pass"] & out["three_study_direction_concordant"]
    )
    out["all_three_lodo_FDR_pass"] = out[fdr_columns].lt(0.05).all(axis=1)
    out["strict_lodo_FDR_candidate"] = (
        out["directional_stability_screen"] & out["all_three_lodo_FDR_pass"]
    )
    out["revised_evidence_class"] = np.select(
        [
            out["strict_lodo_FDR_candidate"],
            out["directional_stability_screen"],
        ],
        ["STRICT_LODO_FDR", "DIRECTIONAL_STABILITY_ONLY"],
        default="OTHER",
    )

    first = [
        "gene_symbol",
        "revised_evidence_class",
        "strict_lodo_FDR_candidate",
        "directional_stability_screen",
        "full_meta_FDR_pass",
        "three_study_direction_concordant",
        "all_three_lodo_FDR_pass",
        "meta_z",
        "meta_p",
        "meta_FDR",
    ]
    remaining = [column for column in out.columns if column not in first]
    out = out[first + remaining].sort_values(
        ["strict_lodo_FDR_candidate", "directional_stability_screen", "meta_FDR", "gene_symbol"],
        ascending=[False, False, True, True],
    )
    out.to_csv(RESULTS / "endothelial_three_study_evidence_v2.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            ("common_genes", len(out)),
            ("full_meta_FDR_lt_0_05", int(out["full_meta_FDR_pass"].sum())),
            (
                "three_study_all_same_direction",
                int(out["three_study_direction_concordant"].sum()),
            ),
            (
                "directional_stability_screen",
                int(out["directional_stability_screen"].sum()),
            ),
            (
                "strict_all_three_lodo_FDR_lt_0_05",
                int(out["strict_lodo_FDR_candidate"].sum()),
            ),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(RESULTS / "endothelial_three_study_summary_v2.tsv", sep="\t", index=False)
    return out, summary


def build_candidate_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = read_tsv("analysis/results/integration/IH_ED_prioritized_genes.tsv")
    cell = read_tsv("analysis/results/integration/IH_ED_celltype_gene_evidence.tsv")
    legacy_shortlist = read_tsv(
        "analysis/results/external_databases/candidate_genes_for_lookup.tsv"
    )

    candidate_rows = cell[cell["evidence_tier"].isin(["A", "B", "C"])].copy()
    gene_stats = (
        candidate_rows.groupby("gene_symbol", as_index=False)
        .agg(
            minimum_ED_p=("ED_p", "min"),
            minimum_ED_FDR=("ED_FDR", "min"),
            best_legacy_tier=("evidence_tier", "min"),
            supporting_cell_types=("cell_type", lambda x: ";".join(sorted(set(x)))),
            nested_celltype_readout_count=("cell_type", "nunique"),
        )
    )
    out = genes.merge(gene_stats, on="gene_symbol", how="left", validate="one_to_one")
    out["revised_candidate_class"] = np.where(
        out["best_tier"].eq("A"),
        "ED_FDR_SUPPORTED",
        "NOMINAL_P_EXPLORATORY",
    )
    out["celltype_count_is_independent_evidence"] = False
    out["legacy_top40_shortlist"] = out["gene_symbol"].isin(legacy_shortlist["gene_symbol"])
    out["legacy_shortlist_rank"] = out["gene_symbol"].map(
        legacy_shortlist.set_index("gene_symbol")["candidate_rank"]
    )

    class_order = {"ED_FDR_SUPPORTED": 0, "NOMINAL_P_EXPLORATORY": 1}
    out["_class_order"] = out["revised_candidate_class"].map(class_order)
    out = out.sort_values(
        ["_class_order", "endothelial_meta_FDR", "minimum_ED_p", "gene_symbol"]
    ).drop(columns="_class_order")
    out.insert(0, "display_order_v2", np.arange(1, len(out) + 1))

    preferred = [
        "display_order_v2",
        "gene_symbol",
        "revised_candidate_class",
        "best_legacy_tier",
        "endothelial_meta_z",
        "endothelial_meta_FDR",
        "minimum_ED_p",
        "minimum_ED_FDR",
        "supporting_cell_types",
        "nested_celltype_readout_count",
        "celltype_count_is_independent_evidence",
        "legacy_top40_shortlist",
        "legacy_shortlist_rank",
    ]
    remaining = [column for column in out.columns if column not in preferred]
    out = out[preferred + remaining]
    out.to_csv(RESULTS / "IH_ED_candidate_classification_v2.tsv", sep="\t", index=False)

    summary = (
        out.groupby("revised_candidate_class", as_index=False)
        .agg(candidate_genes=("gene_symbol", "nunique"))
        .sort_values("revised_candidate_class")
    )
    summary.to_csv(RESULTS / "IH_ED_candidate_classification_summary_v2.tsv", sep="\t", index=False)
    return out, summary


def build_external_scope(candidate_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    shortlist = read_tsv("analysis/results/external_databases/candidate_genes_for_lookup.tsv")
    ot = read_tsv("analysis/results/external_databases/opentargets_candidate_evidence.tsv")
    gwas = read_tsv("analysis/results/external_databases/gwas_candidate_overlaps.tsv")
    reactome = read_tsv("analysis/results/external_databases/reactome_candidate_pathways.tsv")
    legacy_final = read_tsv("analysis/results/integration/final_candidate_ranking.tsv")

    ot_genes = set(ot["gene_symbol"].dropna())
    ot_ed = set(ot.loc[(ot["disease_label"] == "erectile_dysfunction") & (ot["evidence_count"] > 0), "gene_symbol"])
    ot_osa = set(ot.loc[(ot["disease_label"] == "obstructive_sleep_apnea") & (ot["evidence_count"] > 0), "gene_symbol"])
    gwas_genes = set(gwas["candidate_gene"].dropna())
    reactome_genes = set(reactome["gene_symbol"].dropna())
    reactome_queried = set(
        shortlist.loc[shortlist["ensembl_id"].notna()].head(20)["gene_symbol"]
    )
    legacy_final_top20 = legacy_final.nsmallest(20, "final_rank")

    out = candidate_table[["gene_symbol", "revised_candidate_class", "legacy_top40_shortlist", "legacy_shortlist_rank"]].copy()
    out["opentargets_status"] = np.where(out["gene_symbol"].isin(ot_genes), "QUERIED", "NOT_QUERIED")
    out["opentargets_ED_support"] = np.where(
        out["opentargets_status"].eq("QUERIED"),
        np.where(out["gene_symbol"].isin(ot_ed), "PRESENT", "ABSENT_IN_QUERY"),
        "NOT_QUERIED",
    )
    out["opentargets_OSA_support"] = np.where(
        out["opentargets_status"].eq("QUERIED"),
        np.where(out["gene_symbol"].isin(ot_osa), "PRESENT", "ABSENT_IN_QUERY"),
        "NOT_QUERIED",
    )
    out["gwas_exact_trait_filter_status"] = np.where(
        out["legacy_top40_shortlist"], "EVALUATED_IN_LEGACY_TOP40", "NOT_EVALUATED"
    )
    out["gwas_mapped_support"] = np.where(
        out["legacy_top40_shortlist"],
        np.where(out["gene_symbol"].isin(gwas_genes), "PRESENT", "ABSENT_IN_FILTER"),
        "NOT_EVALUATED",
    )
    out["reactome_status"] = np.where(
        out["gene_symbol"].isin(reactome_genes),
        "MAPPED_PATHWAY_PRESENT",
        np.where(
            out["gene_symbol"].isin(reactome_queried),
            "QUERIED_NO_PATHWAY_RETURNED",
            "NOT_QUERIED",
        ),
    )
    out.to_csv(RESULTS / "candidate_external_evidence_scope_v2.tsv", sep="\t", index=False)

    summary_rows = [
        ("expression_candidates_total", len(candidate_table), "all expression candidates"),
        ("legacy_shortlist_total", len(shortlist), "result-informed legacy shortlist"),
        ("opentargets_candidates_queried", len(ot_genes), "mapped members of legacy shortlist"),
        ("opentargets_ED_supported_among_queried", len(ot_ed), "denominator is queried candidates only"),
        ("opentargets_OSA_supported_among_queried", len(ot_osa), "denominator is queried candidates only"),
        ("gwas_candidates_evaluated", len(shortlist), "local filter against legacy shortlist"),
        ("gwas_supported_among_evaluated", len(gwas_genes), "unique mapped candidate genes"),
        ("reactome_candidates_queried", len(reactome_queried), "first 20 mapped legacy candidates"),
        ("reactome_mapped_among_queried", len(reactome_genes), "unique genes with returned pathways"),
        (
            "legacy_final_top20_reactome_queried",
            int(legacy_final_top20["gene_symbol"].isin(reactome_queried).sum()),
            "only this subset of the re-ranked top 20 was actually queried",
        ),
        (
            "legacy_final_top20_recorded_mapping",
            int(legacy_final_top20["reactome_pathway_count"].gt(0).sum()),
            "recorded mappings, not a complete 20-gene coverage rate",
        ),
    ]
    summary = pd.DataFrame(summary_rows, columns=["metric", "value", "denominator_note"])
    summary.to_csv(RESULTS / "external_evidence_scope_summary_v2.tsv", sep="\t", index=False)
    return out, summary


def copy_manuscript_sources() -> None:
    copies = {
        ROOT / "paper_rewriting_output/manuscript_en.src.md": MANUSCRIPT / "manuscript_en.v2.src.md",
        ROOT / "paper_rewriting_output/manuscript_zh.src.md": MANUSCRIPT / "manuscript_zh.v2.src.md",
    }
    for source, destination in copies.items():
        if not destination.exists():
            shutil.copyfile(source, destination)


def validate(
    manifest: pd.DataFrame,
    lodo_summary: pd.DataFrame,
    candidate_summary: pd.DataFrame,
    external_summary: pd.DataFrame,
) -> dict[str, object]:
    metrics = dict(zip(lodo_summary["metric"], lodo_summary["value"]))
    classes = dict(
        zip(candidate_summary["revised_candidate_class"], candidate_summary["candidate_genes"])
    )
    external = dict(zip(external_summary["metric"], external_summary["value"]))
    checks = {
        "manifest_entries_26": len(manifest) == 26,
        "common_genes_5139": metrics.get("common_genes") == 5139,
        "directional_screen_1035": metrics.get("directional_stability_screen") == 1035,
        "strict_lodo_fdr_258": metrics.get("strict_all_three_lodo_FDR_lt_0_05") == 258,
        "ED_FDR_supported_3": classes.get("ED_FDR_SUPPORTED") == 3,
        "nominal_exploratory_563": classes.get("NOMINAL_P_EXPLORATORY") == 563,
        "opentargets_queried_38": external.get("opentargets_candidates_queried") == 38,
        "gwas_evaluated_40": external.get("gwas_candidates_evaluated") == 40,
        "reactome_queried_20": external.get("reactome_candidates_queried") == 20,
        "reactome_mapped_17_of_queried20": external.get("reactome_mapped_among_queried") == 17,
        "legacy_final_top20_only18_queried": external.get("legacy_final_top20_reactome_queried") == 18,
        "legacy_final_top20_recorded_mapping15": external.get("legacy_final_top20_recorded_mapping") == 15,
    }
    passed = all(checks.values())
    report = {
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "evidence_boundary": (
            "This gate corrects labels and denominators from existing v1 outputs. "
            "It does not constitute independent validation, causal evidence, or a fresh database query."
        ),
    }
    validation_dir = REVISION / "04_validation"
    payload = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    write_numbered_attempt(
        validation_dir, "gate_01_validation_attempt", report["status"], payload
    )
    atomic_replace_json(
        validation_dir / "gate_01_validation_latest.json", payload
    )
    if not passed:
        raise RuntimeError(json.dumps(report, ensure_ascii=False))
    return report


def main() -> None:
    for directory in (PROTOCOL, RESULTS, MANUSCRIPT, REVISION / "04_validation"):
        directory.mkdir(parents=True, exist_ok=True)

    manifest = freeze_inputs()
    _, lodo_summary = build_lodo_table()
    candidate_table, candidate_summary = build_candidate_table()
    _, external_summary = build_external_scope(candidate_table)
    copy_manuscript_sources()
    report = validate(manifest, lodo_summary, candidate_summary, external_summary)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
