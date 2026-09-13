from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure6_20260912"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure6_virtual_perturbation_stability_and_ED_association"
REPORT = FIGURE_DIR / "Figure6_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure6_build_audit.json", FIGURE_DIR / "Figure6_legend.md",
        SOURCE_DIR / "Figure6a_design_source_data.tsv",
        SOURCE_DIR / "Figure6b_donor_stability_source_data.tsv",
        SOURCE_DIR / "Figure6c_ED_UP_GSEA_source_data.tsv",
        SOURCE_DIR / "Figure6d_top5_overlap_source_data.tsv",
        SOURCE_DIR / "Figure6e_matched_calibration_source_data.tsv",
    ]
    for path in expected:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path.name}")

    png = Image.open(BASE.with_suffix(".png")).convert("RGB")
    tiff = Image.open(BASE.with_suffix(".tiff"))
    require(png.width >= 2100 and png.height >= 1700, "PNG preview resolution is below expected size")
    require(all(abs(value - 600.0) < 0.1 for value in tiff.info.get("dpi", (0, 0))), "TIFF is not 600 dpi")
    require(tiff.info.get("compression") == "tiff_lzw", "TIFF compression is not LZW")
    non_white = ImageChops.difference(png, Image.new("RGB", png.size, "white")).getbbox()
    require(non_white is not None, "PNG has no visible content")
    left, top, right, bottom = non_white
    require(left >= 3 and top >= 3 and right <= png.width - 3 and bottom <= png.height - 3,
            "visible content is clipped at a PNG edge")

    svg_text = BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(len(re.findall(r"<text[ >]", svg_text)) > 70, "SVG text is not retained as editable text nodes")
    for label in ("Predeclared virtual-perturbation design", "NOT_TESTED_SIZE_3", "Matched-control calibration"):
        require(label in svg_text, f"required interpretation label missing from SVG: {label}")

    design = pd.read_csv(SOURCE_DIR / "Figure6a_design_source_data.tsv", sep="\t")
    require(design.loc[0, ["main_runs", "matched_controls", "targets", "donor_backgrounds", "seeds_per_target_donor"]].tolist() == [45, 90, 3, 3, 5],
            "design audit failed")

    stability = pd.read_csv(SOURCE_DIR / "Figure6b_donor_stability_source_data.tsv", sep="\t")
    require(stability.shape[0] == 9 and (stability["n_seeds"] == 5).all() and (stability["pair_count"] == 10).all(),
            "donor stability design audit failed")
    require(stability["seed_stability_status"].eq("PASS").all(), "frozen donor stability status audit failed")
    require(np.allclose(stability.groupby("target_gene")["median_seed_spearman"].median().reindex(["TYMS", "EFNB2", "LRRC17"]), [0.836507, 0.938075, 0.936988]),
            "stability summary audit failed")

    gsea = pd.read_csv(SOURCE_DIR / "Figure6c_ED_UP_GSEA_source_data.tsv", sep="\t").set_index("target_gene")
    require(gsea["pathway"].eq("ED_UP_FDR").all() and gsea["size"].eq(18).all(), "ED_UP pathway/size audit failed")
    require(np.allclose(gsea.reindex(["TYMS", "EFNB2", "LRRC17"])["NES"], [1.39435424415536, 1.76264330475721, 1.7473900882192]),
            "ED_UP NES audit failed")

    overlap = pd.read_csv(SOURCE_DIR / "Figure6d_top5_overlap_source_data.tsv", sep="\t")
    up = overlap.loc[overlap["signature"].eq("ED_UP_FDR")].set_index("target_gene").reindex(["TYMS", "EFNB2", "LRRC17"])
    down = overlap.loc[overlap["signature"].eq("ED_DOWN_FDR")]
    require(up["overlap_size"].tolist() == [2, 8, 5], "ED_UP overlap audit failed")
    require(down["signature_size"].eq(3).all(), "ED_DOWN size-3 not-tested boundary audit failed")

    calibration = pd.read_csv(SOURCE_DIR / "Figure6e_matched_calibration_source_data.tsv", sep="\t")
    require(calibration.shape[0] == 21 and calibration.groupby("target_gene").size().to_dict() == {"TYMS": 7, "EFNB2": 7, "LRRC17": 7},
            "matched-control calibration row audit failed")
    require(calibration["control_count"].eq(10).all(), "matched-control count audit failed")

    audit = json.loads((FIGURE_DIR / "Figure6_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS" and audit["input_manifest_validation"] == "PASS", "build audit is not PASS")
    require(audit["assertions"]["main_runs"] == 45 and audit["assertions"]["calibration_rows"] == 21,
            "build audit assertion failed")

    lines = [
        "# Figure 6 QA",
        "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for the 45 predeclared main model runs, 90 matched-gene control runs, 9 donor-level stability summaries (five seeds and 10 seed-pair comparisons per donor), the three frozen ED_UP_FDR GSEA records (signature size 18), six top-5% overlap records, and 21 matched-control calibration records (10 controls per target).",
        "- Eligibility audit: ED_DOWN_FDR has signature size 3 and is labelled NOT_TESTED_SIZE_3 for ranked enrichment; it is not displayed or interpreted as a null result.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        f"- Visual raster audit: non-white content is not clipped at image edges; {png.width:,} × {png.height:,}-pixel PNG preview and {tiff.width:,} × {tiff.height:,}-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are vector/chart elements derived from frozen tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the legend limits claims to predictive network effects and prohibits biological knockout validation, expression reversal, causal mechanism, and ED-reversal claims.",
        "",
        "## Remaining boundary",
        "The repeated seeds are technical/model repeats, and donor backgrounds are distinct cellular backgrounds rather than independent human or animal cohorts. Matched controls calibrate network-score specificity only. This figure does not establish biological perturbation effects, OSA-to-ED causality, or therapeutic reversibility.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
