#!/usr/bin/env python3
"""Freeze VP-G07 public-perturbation coverage rules before any result search."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
OUT = VP / "01_protocol" / "VP_G07_external_scope_rule_freeze_v1.json"
G06 = VP / "validation" / "VP_G06_enrichment_COMPLETE.json"
RESULT_DIR = VP / "08_external_reference"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite sealed freeze: {OUT}")
    if RESULT_DIR.exists():
        raise SystemExit(f"Result directory already exists before rule freeze: {RESULT_DIR}")
    g06 = json.loads(G06.read_text(encoding="utf-8"))
    if g06.get("status") != "PASS":
        raise SystemExit("VP-G06 prerequisite is not PASS")

    freeze = {
        "gate": "VP-G07",
        "status": "FROZEN_BEFORE_PUBLIC_SEARCH",
        "version": "v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User explicitly authorized starting VP-G07 on 2026-09-08.",
        "entry_condition": {
            "VP_G06_status": "PASS",
            "VP_G06_sha256": sha256(G06),
        },
        "objective": "Audit public real genetic-perturbation coverage for TYMS, EFNB2, and LRRC17; adjudicate GEARS and UniPert-G2CP applicability without treating cross-cell predictions as external validation.",
        "targets": ["EFNB2", "LRRC17", "TYMS"],
        "scope": {
            "organisms": "All organisms are recorded; human is required for any human-context applicability claim.",
            "assays": "Single-cell transcriptomic genetic perturbation is primary. Bulk transcriptomic genetic perturbation may be retained as secondary evidence but cannot satisfy Perturb-seq or GEARS training coverage.",
            "genetic_modalities": ["CRISPR knockout", "CRISPRi", "CRISPRa", "RNAi/siRNA/shRNA", "overexpression", "other direct genetic perturbation"],
            "excluded_as_target_coverage": ["gene merely measured or mentioned", "compound exposure without direct genetic perturbation", "in-silico-only perturbation", "association-only data", "non-transcriptomic perturbation without an expression readout"],
            "search_cutoff": "Searches are current through the execution date recorded in the search log.",
        },
        "source_priority": [
            "scPerturb official catalogue and source records",
            "PerturbLab official portal and source records",
            "GEARS official repository, processed datasets, and linked primary studies",
            "NCBI GEO/GDS and SRA-linked primary accessions",
            "PubMed primary research articles and linked repositories",
            "bioRxiv primary preprints when no peer-reviewed version is identified",
            "official UniPert repository and the primary UniPert-G2CP article",
        ],
        "minimum_search_families": {
            "catalogue": ["scPerturb", "PerturbLab", "GEARS processed data"],
            "repository": ["NCBI GEO/GDS", "linked SRA or project repository"],
            "literature": ["PubMed", "bioRxiv or publisher primary article"],
        },
        "frozen_query_concepts": {
            "per_target": "<TARGET> AND (Perturb-seq OR CROP-seq OR single-cell CRISPR OR CRISPRi OR CRISPRa OR knockout OR knockdown OR overexpression)",
            "context": "(fibroblast OR endothelial OR vascular OR smooth muscle OR pericyte OR corpus cavernosum OR penile)",
            "catalogue_membership": "Exact canonical target symbol membership in dataset perturbation annotations; aliases require an authoritative canonical-symbol cross-check.",
        },
        "record_unit": "One row per canonical target x dataset/accession x cell context x perturbation modality x time x dose/intensity; split distinct contexts and retain combinations as separate rows.",
        "required_columns": [
            "target_gene", "canonical_symbol_verified", "perturbation_type", "single_or_combination", "cell_line_or_type", "cell_context_tier", "species", "time", "dose_or_intensity", "dataset", "accession", "has_controls", "downloadable", "download_url", "matches_corpus_cavernosum_fibroblast", "evidence_status", "source_url", "primary_citation", "notes"
        ],
        "missing_metadata_policy": "Use NOT_REPORTED or NOT_VERIFIED; never infer time, dose, controls, target membership, or cell identity from silence.",
        "deduplication": "Group the same experiment by stable accession/DOI/repository dataset identifier; preserve multiple source URLs in the evidence table and emit one coverage row per frozen record unit.",
        "downloadability_rule": "YES only when a public landing page or file endpoint for raw or processed expression data is reachable without private credentials at audit time; registration-only or paper-only evidence is recorded separately.",
        "control_rule": "YES only when explicit non-targeting, untreated, vehicle, mock, wild-type, or other stated comparator cells/samples are documented for the expression experiment.",
        "cell_context_tiers": {
            "T1_EXACT": "Human corpus cavernosum fibroblast with direct target perturbation.",
            "T2_RELATED_HUMAN": "Human non-cavernosal fibroblast, endothelial, vascular smooth-muscle, or pericyte context.",
            "T3_OTHER_PRIMARY_HUMAN": "Other primary or stem-cell-derived human context.",
            "T4_HUMAN_CELL_LINE": "Immortalized or cancer human cell line, including K562, RPE1, A549, HEK293, or similar.",
            "T5_NONHUMAN": "Non-human context.",
            "TX_UNKNOWN": "Cell context or species not verified.",
        },
        "coverage_labels": {
            "VERIFIED_DATASET_MEMBER": "Exact target perturbation is verified in machine-readable annotations or the primary experiment record.",
            "PAPER_ONLY": "Primary article supports direct perturbation but downloadable expression data or target membership could not be verified.",
            "SEARCH_HIT_NOT_COVERAGE": "Target appears only in text, measurement, or a nonqualifying context.",
            "NO_VERIFIED_HIT": "No qualifying hit found in the frozen searched sources as of the audit date; this is not proof of global absence.",
        },
        "gears_adjudication": {
            "DIRECT_EXTERNAL_REFERENCE": "Requires T1_EXACT, real target perturbation, single-cell transcriptome, controls, and downloadable data.",
            "RELATED_CONTEXT_REFERENCE": "T2_RELATED_HUMAN with the same evidence requirements; biological context mismatch must remain explicit.",
            "ALGORITHMIC_OR_GENERAL_FUNCTION_REFERENCE": "T3-T5 or unrelated cell lines; may support reproduction or general molecular-function comparison only.",
            "NOT_VALIDATION": "Target absent from actual training perturbations, zero-shot-only prediction, unavailable controls/data, or in-silico-only evidence.",
            "combination_rule": "Combination prediction is externally assessable only when matched combination perturbation observations exist; single-gene data do not validate combinations.",
        },
        "unipert_scope": {
            "allowed": ["256-dimensional gene/protein representation", "nearest-neighbour ranking in the frozen representation space", "version or initialization stability when reproducible assets exist"],
            "not_required_for_gate": "Encoding execution is optional; software/checkpoint availability and scope are audited. Absence of a reproducible checkpoint is NOT_APPLICABLE_YET, not failure of target biology.",
            "prohibited": ["claiming expression reversal", "claiming efficacy in corpus cavernosum cells", "calling representation proximity an experimental validation", "drug recommendation without matched effect data"],
        },
        "g2cp_start_conditions": [
            "same human target-cell background",
            "real TYMS/EFNB2/LRRC17 genetic-perturbation transcriptomes",
            "small matched-background candidate-drug transcriptome set",
            "fixed time, dose, oxygen condition, and negative controls",
            "held-out unseen-drug and unseen-donor evaluation",
        ],
        "g2cp_failure_label": "If any start condition is unmet, label NOT_APPLICABLE_YET rather than FAIL.",
        "negative_result_policy": "No exact or related public dataset is a valid coverage result and does not block technical completion of VP-G07.",
        "completeness_boundary": "Targeted multi-source coverage audit, not a claim of exhaustive global nonexistence. Search counts, returned records, screened records, exclusions, endpoints, parameters, and access dates must be logged.",
        "expected_outputs": [
            "08_external_reference/public_perturbation_target_coverage.tsv",
            "08_external_reference/search_log.tsv",
            "08_external_reference/source_evidence.tsv",
            "08_external_reference/gears_applicability.tsv",
            "08_external_reference/unipert_g2cp_applicability.json",
            "08_external_reference/VP_G07_external_scope_report.md",
            "08_external_reference/artifact_manifest.sha256.tsv",
            "validation/VP_G07_external_scope_AUDIT01.json",
            "validation/VP_G07_external_scope_COMPLETE.json",
        ],
        "acceptance": [
            "All three targets have explicit coverage rows, including NO_VERIFIED_HIT rows where applicable.",
            "Cell context and exact corpus-cavernosum fibroblast match are separately recorded.",
            "GEARS is not used as cross-cell-type external validation or as zero-shot evidence.",
            "Any UniPert output remains representation/ranking only.",
            "All five G2CP start conditions are adjudicated.",
            "Search provenance, source URLs, counts, missing metadata, retained exclusions, and artifact hashes are preserved.",
        ],
        "next_gate": "VP-G08 remains NOT_STARTED and requires separate explicit authorization after VP-G07 completion.",
        "implementation_sha256": {
            "scripts/08a_freeze_vp_g07.py": sha256(Path(__file__)),
        },
    }
    OUT.write_text(json.dumps(freeze, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": freeze["status"], "path": str(OUT), "sha256": sha256(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
