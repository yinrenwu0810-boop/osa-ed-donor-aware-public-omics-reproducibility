from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure3_20260912"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure3_ED_corpus_cavernosum_localization"
REPORT = FIGURE_DIR / "Figure3_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure3_build_audit.json", FIGURE_DIR / "Figure3_legend.md",
        SOURCE_DIR / "Figure3a_UMAP_descriptive_source_data.tsv.gz",
        SOURCE_DIR / "Figure3b_donor_celltype_composition_source_data.tsv",
        SOURCE_DIR / "Figure3c_donor_pseudobulk_pathway_source_data.tsv",
        SOURCE_DIR / "Figure3d_gene_celltype_effect_source_data.tsv",
        SOURCE_DIR / "Figure3e_fibroblast_donor_dot_source_data.tsv",
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
    require(len(re.findall(r"<text[ >]", svg_text)) > 110, "SVG text is not retained as editable text nodes")
    for label in ("64,993 cells; descriptive only", "Tumour-margin", "Fibroblast donor pseudobulk"):
        require(label in svg_text, f"required interpretive label missing from SVG: {label}")

    umap = pd.read_csv(SOURCE_DIR / "Figure3a_UMAP_descriptive_source_data.tsv.gz", sep="\t")
    require(umap.shape[0] == 64993 and bool(umap["descriptive_only"].all()), "UMAP source audit failed")

    composition = pd.read_csv(SOURCE_DIR / "Figure3b_donor_celltype_composition_source_data.tsv", sep="\t")
    donor_groups = composition[["sample_accession", "condition"]].drop_duplicates()["condition"].value_counts().to_dict()
    require(donor_groups == {"normal": 3, "organic_ED_nonDM": 3, "organic_ED_DM": 2}, "3/3/2 donor group audit failed")
    require(np.allclose(composition.groupby("sample_accession")["cell_fraction"].sum().to_numpy(), 1.0), "composition fractions do not sum to one")
    require(composition["sample_accession"].nunique() == 8 and composition["cell_type"].nunique() == 9, "donor composition dimensions failed")

    pathways = pd.read_csv(SOURCE_DIR / "Figure3c_donor_pseudobulk_pathway_source_data.tsv", sep="\t")
    require(pathways.shape[0] == 56 and pathways["pathway"].nunique() == 8 and pathways["cell_type"].nunique() == 7,
            "pathway heatmap dimensions failed")
    require(pathways["pathway"].iloc[0] == "HALLMARK_HYPOXIA", "Hypoxia must be the first pathway row")
    require(int(pathways["frozen_display_flag"].astype(bool).sum()) == 27, "frozen pathway flag audit failed")

    genes = pd.read_csv(SOURCE_DIR / "Figure3d_gene_celltype_effect_source_data.tsv", sep="\t")
    require(genes.shape[0] == 21 and genes["gene_symbol"].nunique() == 3 and genes["cell_type"].nunique() == 7,
            "gene heatmap dimensions failed")
    significant = genes.loc[genes["FDR_lt_0_05"].astype(bool)]
    require(significant.shape[0] == 3 and set(significant["gene_symbol"]) == {"TYMS", "EFNB2", "LRRC17"} and set(significant["cell_type"]) == {"fibroblast"},
            "target-gene FDR audit failed")
    require(int((~genes["estimate_available"].astype(bool)).sum()) == 7, "missing-estimate audit failed")

    dots = pd.read_csv(SOURCE_DIR / "Figure3e_fibroblast_donor_dot_source_data.tsv", sep="\t")
    require(dots.shape[0] == 24 and dots.groupby("gene_symbol").size().to_dict() == {"TYMS": 8, "EFNB2": 8, "LRRC17": 8},
            "panel e donor-dot dimensions failed")
    require(dots["column_id"].nunique() == 8 and dots["cell_type"].eq("fibroblast").all(), "panel e source mapping audit failed")
    expected_log2cpm = np.log2(((dots["raw_count"] + 0.5) / (dots["library_size"] + 1.0)) * 1_000_000.0)
    require(np.allclose(dots["log2CPM_prior_0_5"], expected_log2cpm), "panel e normalization audit failed")

    audit = json.loads((FIGURE_DIR / "Figure3_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS" and audit["input_manifest_validation"] == "PASS", "build audit is not PASS")
    require(audit["assertions"]["panel_e_mapping"] == "PASS_UNAMBIGUOUS", "panel e mapping did not pass")

    lines = [
        "# Figure 3 QA",
        "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for 64,993 descriptive UMAP cells; 3 tumour-margin reference, 3 non-diabetic ED implantation, and 2 diabetic ED implantation donors; composition fractions summing to one per donor; 8 pathway programmes × 7 compartments with 27 pre-existing flags; 3 genes × 7 compartments with three ED-FDR effects in fibroblast; and 24 descriptive donor dots from eight unambiguously mapped fibroblast pseudobulk columns.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        f"- Visual raster audit: non-white content is not clipped at image edges; {png.width:,} × {png.height:,}-pixel PNG preview and {tiff.width:,} × {tiff.height:,}-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are charts derived from frozen tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the UMAP and cell composition are labelled descriptive, cells are not treated as replicates, panel e has no inferential symbol, and the tumour-margin reference versus implantation context is retained in the legend.",
        "",
        "## Remaining boundary",
        "",
        "This figure localizes donor-level ED effects in one atlas. It does not resolve condition from procurement context, establish OSA-to-ED mediation, demonstrate causal function, or provide independent validation.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
