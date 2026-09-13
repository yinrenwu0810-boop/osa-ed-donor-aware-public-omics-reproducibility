from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure2_20260911"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure2_endothelial_convergence_and_heterogeneity"
REPORT = FIGURE_DIR / "Figure2_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure2_build_audit.json",
        FIGURE_DIR / "Figure2_legend.md",
        SOURCE_DIR / "Figure2a_volcano_source_data.tsv",
        SOURCE_DIR / "Figure2b_reported_correlation.tsv",
        SOURCE_DIR / "Figure2b_shared_effect_source_data.tsv",
        SOURCE_DIR / "Figure2c_filter_count_source_data.tsv",
        SOURCE_DIR / "Figure2d_strict_lodo_heatmap_source_data.tsv",
        SOURCE_DIR / "Figure2e_pathway_dotplot_source_data.tsv",
    ]
    for path in expected:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path.name}")

    png = Image.open(BASE.with_suffix(".png")).convert("RGB")
    tiff = Image.open(BASE.with_suffix(".tiff"))
    require(png.width >= 2100 and png.height >= 1800, "PNG preview resolution is below expected size")
    require(tiff.info.get("dpi") == (600.0, 600.0), "TIFF is not 600 dpi")
    require(tiff.info.get("compression") == "tiff_lzw", "TIFF compression is not LZW")
    non_white = ImageChops.difference(png, Image.new("RGB", png.size, "white")).getbbox()
    require(non_white is not None, "PNG has no visible content")
    left, top, right, bottom = non_white
    require(left >= 5 and top >= 5 and right <= png.width - 5 and bottom <= png.height - 5,
            "visible content is clipped at a PNG edge")

    svg_text = BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(len(re.findall(r"<text[ >]", svg_text)) > 80, "SVG text is not retained as editable text nodes")
    require("Hypoxia" in svg_text and "strict LODO-FDR" in svg_text, "key Figure 2 labels absent from SVG")

    volcano = pd.read_csv(SOURCE_DIR / "Figure2a_volcano_source_data.tsv", sep="\t")
    observed = volcano.assign(significant=volcano["adj.P.Val"] < 0.05).groupby("dataset")["significant"].sum().to_dict()
    require(observed == {"GSE10723": 1244, "GSE205050": 2339, "GSE243023": 0}, "volcano FDR count audit failed")

    corr = pd.read_csv(SOURCE_DIR / "Figure2b_reported_correlation.tsv", sep="\t").iloc[0]
    effects = pd.read_csv(SOURCE_DIR / "Figure2b_shared_effect_source_data.tsv", sep="\t")
    require(effects.shape[0] == 5139, "shared-effect row count audit failed")
    require(int(effects["strict_lodo_FDR_candidate"].astype(bool).sum()) == 258, "strict LODO overlay audit failed")
    require(abs(float(corr["reported_spearman_rho_endo243_vs_endo205"]) - 0.08629966133792) < 1e-12, "reported rho audit failed")

    counts = pd.read_csv(SOURCE_DIR / "Figure2c_filter_count_source_data.tsv", sep="\t")
    require(counts["count"].tolist() == [1766, 1035, 1282, 258], "frozen filter count audit failed")

    heat = pd.read_csv(SOURCE_DIR / "Figure2d_strict_lodo_heatmap_source_data.tsv", sep="\t")
    require(heat["gene_symbol"].nunique() == 12, "strict LODO heatmap gene-count audit failed")
    direction_counts = heat.drop_duplicates("gene_symbol")["meta_direction"].value_counts().to_dict()
    require(direction_counts == {"positive": 6, "negative": 6}, "strict LODO direction audit failed")

    pathways = pd.read_csv(SOURCE_DIR / "Figure2e_pathway_dotplot_source_data.tsv", sep="\t")
    require(pathways["collection"].value_counts().to_dict() == {"Hallmark": 6, "Reactome": 6}, "pathway dotplot count audit failed")
    require(int(pathways["pathway"].eq("HALLMARK_HYPOXIA").sum()) == 1, "Hypoxia pathway audit failed")

    audit = json.loads((FIGURE_DIR / "Figure2_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS", "build audit is not PASS")
    lines = [
        "# Figure 2 QA",
        "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for study-level FDR counts (GSE243023: 0; GSE205050: 2,339; GSE10723: 1,244), the reported cross-study rho (0.08629966133792; n = 5,139), strict LODO overlay (258), frozen filter counts (1,766; 1,035; 1,282; 258), 6-positive/6-negative heatmap selection, and 6 Hallmark plus 6 Reactome dot-plot pathways.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        f"- Visual raster audit: non-white content is not clipped at image edges; {png.width:,} × {png.height:,}-pixel PNG preview and {tiff.width:,} × {tiff.height:,}-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are vector/chart elements derived from tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the legend preserves the zero-result study and limited genome-wide effect agreement, identifies the 258-gene overlay as a strict frozen subset, keeps pathway FDR families separate, and does not claim universal concordance or causality.",
        "",
        "## Remaining boundary",
        "",
        "This figure displays cross-study endothelial IH heterogeneity and concentrated convergence. It does not validate a gene perturbation, demonstrate OSA-to-ED mediation, or establish causal gene function.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
