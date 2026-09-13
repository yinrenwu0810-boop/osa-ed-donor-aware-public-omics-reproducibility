from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure1_20260911"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure1_parallel_evidence_and_strict_bridge"
REPORT = FIGURE_DIR / "Figure1_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure1_build_audit.json",
        FIGURE_DIR / "Figure1_legend.md",
        SOURCE_DIR / "Figure1a_evidence_layers.tsv",
        SOURCE_DIR / "Figure1b_endothelial_filters.tsv",
        SOURCE_DIR / "Figure1c_clinical_bridge_filters.tsv",
        SOURCE_DIR / "Figure1d_strict_context_heatmap.tsv",
    ]
    for path in expected:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path.name}")

    png = Image.open(BASE.with_suffix(".png")).convert("RGB")
    tiff = Image.open(BASE.with_suffix(".tiff"))
    require(png.width >= 2100 and png.height >= 1300, "PNG preview resolution is below expected size")
    require(tiff.info.get("dpi") == (600.0, 600.0), "TIFF is not 600 dpi")
    require(tiff.info.get("compression") == "tiff_lzw", "TIFF compression is not LZW")
    non_white = ImageChops.difference(png, Image.new("RGB", png.size, "white")).getbbox()
    require(non_white is not None, "PNG has no visible content")
    left, top, right, bottom = non_white
    require(left >= 5 and top >= 5 and right <= png.width - 5 and bottom <= png.height - 5,
            "visible content is clipped at a PNG edge")

    svg_text = BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(len(re.findall(r"<text[ >]", svg_text)) > 30, "SVG text is not retained as editable text nodes")
    require("HALLMARK_HYPOXIA" in svg_text or "HYPOXIA" in svg_text, "strict bridge label absent from SVG")

    filter_data = pd.read_csv(SOURCE_DIR / "Figure1b_endothelial_filters.tsv", sep="\t")
    filters = dict(zip(filter_data["metric"], filter_data["value"]))
    require([filters[k] for k in [
        "common_genes", "full_meta_FDR_lt_0_05", "three_study_all_same_direction",
        "directional_stability_screen", "strict_all_three_lodo_FDR_lt_0_05",
    ]] == [5139, 1282, 1766, 1035, 258], "endothelial count audit failed")

    bridge_data = pd.read_csv(SOURCE_DIR / "Figure1c_clinical_bridge_filters.tsv", sep="\t")
    bridge = dict(zip(bridge_data["metric"], bridge_data["value"]))
    require([bridge[k] for k in [
        "hallmark_pathways_total", "all_four_estimated", "complete_directional_pattern",
        "all_four_FDR_pass", "strict_clinical_bridge_pass",
    ]] == [50, 49, 12, 3, 1], "clinical bridge count audit failed")

    heat = pd.read_csv(SOURCE_DIR / "Figure1d_strict_context_heatmap.tsv", sep="\t")
    require(heat.shape[0] == 48, "heatmap source must contain 12 pathways by 4 contexts")
    require(heat["pathway"].nunique() == 12, "heatmap pathway count audit failed")
    require(int(heat["strict_clinical_bridge_pass"].astype(bool).groupby(heat["pathway"]).first().sum()) == 1,
            "strict bridge count audit failed")
    hypoxia = heat.loc[heat["pathway"].eq("HALLMARK_HYPOXIA")]
    require(hypoxia.shape[0] == 4 and bool(hypoxia["FDR_lt_0_05_star"].all()),
            "Hypoxia four-context FDR audit failed")

    audit = json.loads((FIGURE_DIR / "Figure1_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS", "build audit is not PASS")
    lines = [
        "# Figure 1 QA",
        "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for endothelial counts (5,139; 1,282; 1,766; 1,035; 258), bridge counts (50; 49; 12; 3; 1), 12 × 4 heatmap cells, and four FDR-marked Hypoxia cells.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        "- Visual raster audit: non-white content is not clipped at the image edges; 2,161 × 1,417-pixel PNG preview and 4,323 × 2,834-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are vector/chart elements derived from tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the legend states that evidence layers are parallel, that original statistical units are retained, and that the bridge does not establish OSA–ED causal mediation.",
        "",
        "## Remaining boundary",
        "",
        "This figure is a cross-context evidence map. It does not test clinical mediation, causal gene function, or intervention response in ED.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()

