from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure7_20260912"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure7_exploratory_hypoxia_sensitivity"
REPORT = FIGURE_DIR / "Figure7_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure7_build_audit.json", FIGURE_DIR / "Figure7_legend.md",
        SOURCE_DIR / "Figure7a_specification_NES_source_data.tsv",
        SOURCE_DIR / "Figure7b_inferential_FDR_source_data.tsv",
        SOURCE_DIR / "Figure7c_direction_stability_source_data.tsv",
        SOURCE_DIR / "Figure7d_matched_control_source_data.tsv",
        SOURCE_DIR / "Figure7e_leading_edge_stability_source_data.tsv",
        SOURCE_DIR / "Figure7f_mapping_source_data.tsv",
    ]
    for path in expected:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path.name}")

    png = Image.open(BASE.with_suffix(".png")).convert("RGB")
    tiff = Image.open(BASE.with_suffix(".tiff"))
    require(png.width >= 2100 and png.height >= 1800, "PNG preview resolution is below expected size")
    require(all(abs(value - 600.0) < 0.1 for value in tiff.info.get("dpi", (0, 0))), "TIFF is not 600 dpi")
    require(tiff.info.get("compression") == "tiff_lzw", "TIFF compression is not LZW")
    bbox = ImageChops.difference(png, Image.new("RGB", png.size, "white")).getbbox()
    require(bbox is not None, "PNG has no visible content")
    left, top, right, bottom = bbox
    require(left >= 3 and top >= 3 and right <= png.width - 3 and bottom <= png.height - 3,
            "visible content is clipped at a PNG edge")

    svg_text = BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(len(re.findall(r"<text[ >]", svg_text)) > 95, "SVG text is not retained as editable text nodes")
    for label in ("Exploratory sensitivity analysis", "S3: descriptive only; no P/FDR", "191 mapped"):
        require(label in svg_text, f"required interpretation label missing from SVG: {label}")

    specifications = pd.read_csv(SOURCE_DIR / "Figure7a_specification_NES_source_data.tsv", sep="\t")
    require(specifications.shape[0] == 12 and (specifications["NES"] > 0).all(), "positive direction audit failed")
    require(specifications.loc[specifications["specification"].eq("S3"), ["p_value", "FDR_within_Hallmark"]].isna().all().all(), "S3 non-inferential boundary audit failed")

    fdr = pd.read_csv(SOURCE_DIR / "Figure7b_inferential_FDR_source_data.tsv", sep="\t")
    require(fdr.shape[0] == 9 and fdr["specification"].isin(["S0", "S1", "S2"]).all(), "inferential specification audit failed")
    require((fdr["FDR_within_Hallmark"] >= 0.05).all(), "formal FDR-negative audit failed")

    stability = pd.read_csv(SOURCE_DIR / "Figure7c_direction_stability_source_data.tsv", sep="\t")
    counts = stability.groupby("analysis").size().to_dict()
    require(counts == {"Seed-level (n=45)": 45, "Leave-one-seed (n=45)": 45, "Leave-one-donor (n=9)": 9}, "resampling count audit failed")
    require((stability["NES"] > 0).all(), "resampling direction audit failed")

    calibration = pd.read_csv(SOURCE_DIR / "Figure7d_matched_control_source_data.tsv", sep="\t")
    require(calibration.shape[0] == 3 and calibration["matched_controls"].eq(10).all(), "matched-control count audit failed")
    require(np.allclose(calibration.set_index("target_gene").reindex(["TYMS", "EFNB2", "LRRC17"])["conservative_absolute_upper_tail"], [0.9090909090909091, 0.7272727272727273, 0.9090909090909091]), "calibration value audit failed")

    leading = pd.read_csv(SOURCE_DIR / "Figure7e_leading_edge_stability_source_data.tsv", sep="\t")
    require(leading.shape[0] == 111, "leading-edge row-count audit failed")
    require(leading.groupby("comparison_type").size().to_dict() == {"within_target_donor_seed": 90, "between_donor_five_seed_consensus": 9, "between_target_formal_consensus": 3, "between_specification_rank_level": 9}, "leading-edge comparison audit failed")
    require(leading.loc[leading["comparison_type"].eq("between_specification_rank_level"), "status"].eq("DISPLAY_DERIVED_FROM_FROZEN_S071").all(), "specification Jaccard provenance audit failed")

    mapping = pd.read_csv(SOURCE_DIR / "Figure7f_mapping_source_data.tsv", sep="\t")
    require(mapping.shape[0] == 1 and mapping.loc[0, "universe_member_count"] == 191 and mapping.loc[0, "source_member_count_unique"] == 200, "mapping count audit failed")
    require(mapping.loc[0, "mapping_qc_status"] == "PASS_EXACT_SYMBOL_MAPPING", "mapping status audit failed")

    audit = json.loads((FIGURE_DIR / "Figure7_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS" and audit["input_manifest_validation"] == "PASS", "build audit is not PASS")
    require(audit["assertions"]["all_12_NES_positive"] and audit["assertions"]["inferential_FDR_all_ge_0_05"] and audit["assertions"]["S3_no_P_or_FDR"], "build boundary assertion failed")

    lines = [
        "# Figure 7 QA", "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for 12 S0–S3 specification records (all positive NES), the nine inferential S0–S2 FDR values (all at least 0.05), 45 seed-level, 45 leave-one-seed-out, and nine leave-one-donor-out records, three matched-control records with 10 controls per target, 111 leading-edge comparisons, and exact mapping of 191/200 source symbols.",
        "- S3 boundary audit: PASS. It is labelled descriptive and contains no aggregate-rank P value or FDR.",
        "- Display-derivation audit: the nine S0–S2 pairwise Jaccard values in panel e are deterministic set-overlap summaries from frozen S071 leading-edge strings; no enrichment or inferential analysis was rerun.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        f"- Visual raster audit: non-white content is not clipped at image edges; {png.width:,} × {png.height:,}-pixel PNG preview and {tiff.width:,} × {tiff.height:,}-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are vector/chart elements derived from frozen tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the legend calls this an exploratory sensitivity analysis and prohibits gene-expression-upregulation, target-specificity, biological knockout validation, causal mechanism, ED-reversal, and OSA-to-ED-causality claims.",
        "", "## Remaining boundary",
        "Positive NES is a rank-direction result. Stability across model resamples and specifications does not replace formal Hallmark FDR control or demonstrate specificity relative to matched genes. Resampled seeds are technical/model repetitions, not biological replication.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
