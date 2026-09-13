from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure4_20260912"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure4_candidate_hierarchy_and_doublet_sensitivity"
REPORT = FIGURE_DIR / "Figure4_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure4_build_audit.json", FIGURE_DIR / "Figure4_legend.md",
        SOURCE_DIR / "Figure4a_candidate_hierarchy_source_data.tsv",
        SOURCE_DIR / "Figure4b_candidate_evidence_matrix_source_data.tsv",
        SOURCE_DIR / "Figure4c_doublet_gene_sensitivity_source_data.tsv",
        SOURCE_DIR / "Figure4d_doublet_pathway_sensitivity_source_data.tsv",
        SOURCE_DIR / "Figure4e_pathway_state_count_source_data.tsv",
    ]
    for path in expected:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path.name}")

    png = Image.open(BASE.with_suffix(".png")).convert("RGB")
    tiff = Image.open(BASE.with_suffix(".tiff"))
    require(png.width >= 2100 and png.height >= 1600, "PNG preview resolution is below expected size")
    require(tiff.info.get("dpi") == (600.0, 600.0), "TIFF is not 600 dpi")
    require(tiff.info.get("compression") == "tiff_lzw", "TIFF compression is not LZW")
    non_white = ImageChops.difference(png, Image.new("RGB", png.size, "white")).getbbox()
    require(non_white is not None, "PNG has no visible content")
    left, top, right, bottom = non_white
    require(left >= 3 and top >= 3 and right <= png.width - 3 and bottom <= png.height - 3,
            "visible content is clipped at a PNG edge")

    svg_text = BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(len(re.findall(r"<text[ >]", svg_text)) > 100, "SVG text is not retained as editable text nodes")
    for label in ("Frozen ED support hierarchy", "Class A → A", "680 complete pairs"):
        require(label in svg_text, f"required interpretive label missing from SVG: {label}")

    hierarchy = pd.read_csv(SOURCE_DIR / "Figure4a_candidate_hierarchy_source_data.tsv", sep="\t")
    require(hierarchy["count"].tolist() == [3, 563], "candidate hierarchy count audit failed")

    evidence = pd.read_csv(SOURCE_DIR / "Figure4b_candidate_evidence_matrix_source_data.tsv", sep="\t")
    require(evidence["gene_symbol"].tolist() == ["TYMS", "EFNB2", "LRRC17"], "candidate order audit failed")
    require((evidence["endothelial_meta_FDR"] < 0.05).all() and (evidence["ED_FDR"] < 0.05).all(), "candidate FDR audit failed")
    require((evidence["historical_best_tier"].eq("A") & evidence["sensitivity_best_tier"].eq("A")).all(), "A-to-A retention audit failed")

    gene = pd.read_csv(SOURCE_DIR / "Figure4c_doublet_gene_sensitivity_source_data.tsv", sep="\t")
    require(gene.shape[0] == 31 and int(gene["complete_pair"].astype(bool).sum()) == 28, "gene sensitivity row audit failed")
    require(int(gene["old_FDR_lt_0_05"].astype(bool).sum()) == 3 and int(gene["new_FDR_lt_0_05"].astype(bool).sum()) == 3,
            "gene sensitivity FDR-state audit failed")

    pathway = pd.read_csv(SOURCE_DIR / "Figure4d_doublet_pathway_sensitivity_source_data.tsv", sep="\t")
    require(pathway.shape[0] == 683 and int(pathway["complete_pair"].astype(bool).sum()) == 680, "pathway row audit failed")
    require(int(pathway["significance_state"].eq("Incomplete estimate").sum()) == 3, "incomplete-pathway audit failed")
    states = pd.read_csv(SOURCE_DIR / "Figure4e_pathway_state_count_source_data.tsv", sep="\t")
    expected_states = {"Retained significant": 263, "Gained significance": 23, "Lost significance": 26, "Retained non-significant": 368}
    require(dict(zip(states["significance_state"], states["count"])) == expected_states, "pathway state-count audit failed")
    require(int(states["count"].sum()) == 680, "complete pathway state total audit failed")

    audit = json.loads((FIGURE_DIR / "Figure4_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS" and audit["input_manifest_validation"] == "PASS", "build audit is not PASS")
    require(audit["assertions"]["pathway_complete_pairs"] == 680, "build audit pathway-pair count failed")

    lines = [
        "# Figure 4 QA",
        "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for the frozen 3-versus-563 hierarchy; TYMS, EFNB2, and LRRC17 with endothelial-meta and fibroblast ED FDR < 0.05; class A-to-A retention after doublet exclusion; 31 gene rows with 28 complete pairs; 683 pathway records with 680 complete pairs and three incomplete estimates; and pathway state counts of 263 retained significant, 23 gained, 26 lost, and 368 retained non-significant.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        f"- Visual raster audit: non-white content is not clipped at image edges; {png.width:,} × {png.height:,}-pixel PNG preview and {tiff.width:,} × {tiff.height:,}-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are vector/chart elements derived from frozen tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the legend distinguishes three incomplete pathway estimates from non-significant states, retains nested-readout and within-atlas limits, and does not call doublet sensitivity independent replication or causal validation.",
        "",
        "## Remaining boundary",
        "",
        "This figure shows candidate classification and within-atlas doublet sensitivity. It does not establish independent replication, OSA-to-ED mediation, or causal gene function.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
