#!/usr/bin/env python
"""Fail-closed acceptance checks for Gate 04 v2 figures."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "02_results"
FIGURES = RESULTS / "figures"
SOURCE = RESULTS / "figure_source_data"
VALIDATION = ROOT / "04_validation"


def next_attempt_path(status: str) -> Path:
    existing = list(VALIDATION.glob("gate_04_figure_rebuild_attempt_*_*.json"))
    numbers = []
    for path in existing:
        try:
            numbers.append(int(path.name.split("_attempt_")[1].split("_")[0]))
        except (IndexError, ValueError):
            continue
    return VALIDATION / f"gate_04_figure_rebuild_attempt_{max(numbers, default=0) + 1:02d}_{status}.json"


def require(condition: bool, name: str, checks: dict, detail: str = "") -> None:
    checks[name] = {"pass": bool(condition), "detail": detail}
    if not condition:
        raise AssertionError(f"{name}: {detail}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manual-render-review", action="store_true",
                        help="Record only after reviewer has visually inspected PNG renders.")
    args = parser.parse_args()
    checks: dict = {}
    try:
        for stem in ("Figure1_parallel_triangulation_v2", "Figure3_candidate_sensitivity_v2"):
            for suffix in (".svg", ".pdf", ".png", ".tiff"):
                path = FIGURES / f"{stem}{suffix}"
                require(path.is_file() and path.stat().st_size > 10_000,
                        f"export_{stem}{suffix}", checks, str(path))
            svg = (FIGURES / f"{stem}.svg").read_text(encoding="utf-8")
            require("<text" in svg and "<path" in svg,
                    f"editable_svg_text_{stem}", checks,
                    "SVG must retain text elements and vector paths.")

        heat = pd.read_csv(SOURCE / "Figure1_Hallmark_bridge_source_data_v2.tsv", sep="\t")
        require(len(heat) == 48, "figure1_heatmap_rows", checks, "12 pathways x 4 contexts")
        require(heat["pathway"].nunique() == 12 and heat["context"].nunique() == 4,
                "figure1_complete_directional_grid", checks, "Expected 12x4 matrix")
        require(heat["strict_clinical_bridge_pass"].sum() == 4 and
                set(heat.loc[heat["strict_clinical_bridge_pass"], "pathway"]) == {"HALLMARK_HYPOXIA"},
                "figure1_unique_strict_bridge", checks, "Only HALLMARK_HYPOXIA may be strict.")
        require(heat["FDR_lt_0_05_star"].equals(heat["FDR"].lt(0.05)),
                "figure1_fdr_star_rule", checks, "Every star flag is computed from FDR < 0.05.")

        counts = pd.read_csv(SOURCE / "Figure3_candidate_class_count_source_data_v2.tsv", sep="\t")
        count_map = dict(zip(counts["revised_candidate_class"], counts["candidate_genes"]))
        require(count_map == {"ED_FDR_SUPPORTED": 3, "NOMINAL_P_EXPLORATORY": 563},
                "figure3_candidate_strata", checks, "Frozen candidate class counts")
        effect = pd.read_csv(SOURCE / "Figure3_fibroblast_effect_source_data_v2.tsv", sep="\t")
        require(set(effect["gene_symbol"]) == {"TYMS", "EFNB2", "LRRC17"} and len(effect) == 3,
                "figure3_named_candidates", checks, "Three FDR-supported candidates only")
        require((effect["old_FDR"] < 0.05).all() and (effect["new_FDR"] < 0.05).all(),
                "figure3_doublet_fdr_retention", checks, "All three retain FDR < 0.05")
        script_text = (ROOT / "scripts" / "build_gate04_figures_v2.py").read_text(encoding="utf-8")
        banned = ("final_candidate_ranking.tsv", "final_rank", "priority_index")
        require(not any(token in script_text for token in banned),
                "figure3_no_legacy_ranking", checks, "No legacy ranking input or field")
        require(args.manual_render_review, "manual_png_layout_review", checks,
                "Reviewer must inspect both generated PNG renders before PASS.")
        status = "PASS"
    except Exception as exc:  # Preserve evidence of the failed attempt.
        status = "FAIL"
        checks.setdefault("exception", {"pass": False, "detail": str(exc)})
    result = {
        "gate": "04",
        "status": status,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "backend": "Python / Matplotlib only",
        "manual_png_layout_review": bool(args.manual_render_review),
        "checks": checks,
        "evidence_boundary": (
            "This validates v2 figure provenance, exports and displayed evidence labels. "
            "It does not establish causal effects, independent replication, or target-journal acceptance."
        ),
    }
    path = next_attempt_path(status)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    latest = VALIDATION / "gate_04_figure_rebuild_latest.json"
    latest.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
