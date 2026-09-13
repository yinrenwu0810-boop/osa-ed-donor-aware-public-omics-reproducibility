#!/usr/bin/env python3
"""Freeze the complete VP-G06 analysis design before reading G06 results."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
PROJECT = VP.parents[1]
OUT = VP / "01_protocol" / "VP_G06_enrichment_rule_freeze_v1.json"
RSCRIPT = Path(r"D:\R\R-4.6.1\bin\Rscript.exe")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def status(path: Path) -> str:
    return json.loads(path.read_text(encoding="utf-8"))["status"]


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite sealed freeze: {OUT}")
    complete = VP / "validation" / "VP_G05_stability_COMPLETE.json"
    control_audit = VP / "validation" / "VP_G05_control_runs_v2_AUDIT01.json"
    if status(complete) != "PASS" or status(control_audit) != "PASS":
        raise SystemExit("VP-G05 prerequisite is not PASS")
    result_dir = VP / "07_enrichment" / "v1"
    if result_dir.exists():
        raise SystemExit("G06 result directory already exists before freeze")

    inputs = [
        VP / "04_prepared" / "gene_universe.tsv",
        VP / "06_consensus" / "candidate_donor_consensus.tsv",
        VP / "06_consensus" / "all_virtual_KO_rankings.tsv.gz",
        VP / "06_consensus" / "negative_control_selection_v2" / "negative_control_selection.tsv",
        VP / "01_protocol" / "control_execution_v2" / "control_task_matrix.tsv",
        VP / "01_protocol" / "control_execution_v2" / "control_batch_matrix.tsv",
        VP / "05_runs" / "run_registry.tsv",
        PROJECT / "analysis" / "data" / "processed" / "gene_sets" / "Hallmark.gmt",
        PROJECT / "analysis" / "data" / "processed" / "gene_sets" / "Reactome.gmt",
        PROJECT / "analysis" / "data" / "processed" / "gene_sets" / "GO_BP.gmt",
        PROJECT / "revision_v2" / "01_work" / "GSE206528_doublet_sensitivity" / "attempt_02" / "edgeR_Hallmark" / "differential_fibroblast_nonDM_ED_vs_normal.tsv",
        PROJECT / "revision_v2" / "01_work" / "GSE206528_doublet_sensitivity" / "attempt_02" / "edgeR_Hallmark" / "analysis_manifest.tsv",
    ]
    missing = [str(path) for path in inputs if not path.is_file()]
    if missing:
        raise SystemExit(f"Missing frozen inputs: {missing}")

    vp_g01_manifest = VP / "01_protocol" / "input_manifest.sha256.tsv"
    with vp_g01_manifest.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    prior_gmt = {row["relative_path"]: row for row in rows if row["relative_path"].endswith(".gmt")}
    for name in ["Hallmark.gmt", "Reactome.gmt", "GO_BP.gmt"]:
        key = f"analysis/data/processed/gene_sets/{name}"
        path = PROJECT / key
        if prior_gmt.get(key, {}).get("status") != "PASS" or prior_gmt[key]["sha256_observed"] != sha256(path):
            raise SystemExit(f"VP-G01 GMT provenance mismatch: {name}")

    r_version = subprocess.run([str(RSCRIPT), "--version"], capture_output=True, text=True, encoding="utf-8", errors="replace").stderr.strip()
    package_cmd = (
        f".libPaths(c('{VP.as_posix()}/02_env/Rlib',.libPaths()));"
        "cat(as.character(packageVersion('fgsea')),as.character(packageVersion('data.table')),"
        "as.character(packageVersion('jsonlite')),as.character(packageVersion('digest')),sep='|')"
    )
    packages = subprocess.run([str(RSCRIPT), "-e", package_cmd], capture_output=True, text=True, encoding="utf-8", errors="replace", check=True).stdout.strip().split("|")
    implementations = ["scripts/07_run_g06_enrichment.R", "scripts/07a_validate_g06_enrichment.R"]
    resource_path = VP / "03_preflight" / "resource_detection_VP_G06_20260908.json"
    compatibility_path = VP / "03_preflight" / "resource_detection_VP_G06_20260908_COMPATIBILITY_NOTE.json"
    freeze = {
        "gate": "VP-G06",
        "status": "FROZEN_BEFORE_ENRICHMENT",
        "version": "v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User explicitly authorized starting the next gate after VP-G05 completion.",
        "entry_condition": {"VP_G05_status": "PASS", "VP_G05_sha256": sha256(complete)},
        "independent_unit": "donor; genes and pathways are not biological replicates",
        "fixed_gene_universe": "Exactly the 1388 genes in 04_prepared/gene_universe.tsv; no whole-genome background and no significance-only subset.",
        "target_rank_input": "For each target, rank all 1388 genes by the VP-G05 across-donor median standardized perturbation rank descending; gene symbol ascending breaks ties; convert to a unique ordinal score (n-rank+1)/n.",
        "donor_rank_input": "For HALLMARK_HYPOXIA donor consistency, use each target-donor five-seed consensus full ranking with the same deterministic ordinal conversion.",
        "databases": {
            "Hallmark": {"path": "analysis/data/processed/gene_sets/Hallmark.gmt", "expected_sets": 50},
            "Reactome": {"path": "analysis/data/processed/gene_sets/Reactome.gmt", "expected_sets": 1839},
            "GO_BP": {"path": "analysis/data/processed/gene_sets/GO_BP.gmt", "expected_sets": 7538},
        },
        "gene_set_eligibility": "Intersect each set with the fixed universe; test overlap sizes 15 through 500 inclusive. Preserve all tested non-significant results and an eligibility row for every source set.",
        "method": {
            "engine": "fgsea::fgseaMultilevel",
            "minSize": 15,
            "maxSize": 500,
            "eps": 0,
            "scoreType": "std",
            "gseaParam": 1,
            "nproc": 1,
            "random_seed": 2026090806,
            "seed_assignment": "Increment by one for every call in sorted target, database, donor, and matched-control order; record call_seed per result row.",
            "multiple_testing": "BH FDR computed independently by fgsea within each target-profile and database call.",
        },
        "hypoxia_primary": "Report HALLMARK_HYPOXIA ES, NES, nominal p, within-Hallmark FDR, set size, leading edge, and 3-donor NES consistency.",
        "matched_control_hypoxia": "For fair calibration, aggregate primary-seed 2026082701 ranks across three donors for each target and each of its 10 v2 controls; run the complete Hallmark collection before extracting HALLMARK_HYPOXIA. Report signed and absolute NES percentiles plus the conservative absolute upper-tail value.",
        "primary_seed": 2026082701,
        "ED_signatures": {
            "source": "doublet-excluded fibroblast nonDM-ED versus normal edgeR table",
            "ED_UP_FDR": "FDR < 0.05 and logFC > 0 within the fixed 1388-gene universe",
            "ED_DOWN_FDR": "FDR < 0.05 and logFC < 0 within the fixed 1388-gene universe",
            "ranked_GSEA": "Test ED_UP_FDR and ED_DOWN_FDR against each target full perturbation rank with minSize 5, maxSize 500, and BH within target.",
            "top5_overlap": "One-sided hypergeometric overlap with the frozen 2-of-3-donor top-5% consensus; background N=1388; BH across the two ED signatures within each target.",
            "full_rank_concordance": "Descriptive Spearman correlation between perturbation score and absolute signed sqrt(F) ED statistic across all 1388 genes.",
        },
        "direction_boundary": "Positive NES means pathway members concentrate at the high network-perturbation end; negative NES means concentration at the low end. Neither is pathway up/down regulation. ED overlap is not expression reversal.",
        "prohibited": ["thresholding the target list before GSEA", "whole-genome background", "cross-database pooled FDR", "calling negative NES down-regulation", "calling ED overlap expression reversal", "causal or therapeutic claims"],
        "environment": {"Rscript": r_version, "fgsea": packages[0], "data.table": packages[1], "jsonlite": packages[2], "digest": packages[3]},
        "resource_preflight": {"sha256": sha256(resource_path), "compatibility_note_sha256": sha256(compatibility_path), "execution": "Sequential fgsea calls with nproc=1 for deterministic auditing; detected 10.29 GB available RAM and 184.39 GB free disk."},
        "input_sha256": {str(path.relative_to(PROJECT)).replace("\\", "/"): sha256(path) for path in inputs},
        "control_run_audit_sha256": sha256(control_audit),
        "implementation_sha256": {script: sha256(VP / script) for script in implementations},
        "expected_result_directory": "07_enrichment/v1",
        "next_gate": "VP-G07 requires explicit authorization after VP_G06_enrichment_COMPLETE.json passes.",
    }
    OUT.write_text(json.dumps(freeze, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(freeze, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
