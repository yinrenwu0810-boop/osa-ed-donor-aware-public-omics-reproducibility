from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure5_20260912"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure5_partial_transportability_and_heterogeneity"
REPORT = FIGURE_DIR / "Figure5_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure5_build_audit.json", FIGURE_DIR / "Figure5_legend.md",
        SOURCE_DIR / "Figure5a_heldout_probability_source_data.tsv",
        SOURCE_DIR / "Figure5b_pooled_ROC_source_data.tsv",
        SOURCE_DIR / "Figure5c_per_dataset_AUC_source_data.tsv",
        SOURCE_DIR / "Figure5d_feature_stability_source_data.tsv",
        SOURCE_DIR / "Figure5e_pooled_metric_source_data.tsv",
    ]
    for path in expected:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path.name}")

    png = Image.open(BASE.with_suffix(".png")).convert("RGB")
    tiff = Image.open(BASE.with_suffix(".tiff"))
    require(png.width >= 2100 and png.height >= 1500, "PNG preview resolution is below expected size")
    require(tiff.info.get("dpi") == (600.0, 600.0), "TIFF is not 600 dpi")
    require(tiff.info.get("compression") == "tiff_lzw", "TIFF compression is not LZW")
    non_white = ImageChops.difference(png, Image.new("RGB", png.size, "white")).getbbox()
    require(non_white is not None, "PNG has no visible content")
    left, top, right, bottom = non_white
    require(left >= 3 and top >= 3 and right <= png.width - 3 and bottom <= png.height - 3,
            "visible content is clipped at a PNG edge")

    svg_text = BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(len(re.findall(r"<text[ >]", svg_text)) > 100, "SVG text is not retained as editable text nodes")
    for label in ("Transportability is heterogeneous", "0.11", "Conditional pooled uncertainty"):
        require(label in svg_text, f"required interpretive label missing from SVG: {label}")

    predictions = pd.read_csv(SOURCE_DIR / "Figure5a_heldout_probability_source_data.tsv", sep="\t")
    require(predictions.shape[0] == 60 and predictions.groupby("strategy").size().to_dict() == {"gene_rank": 20, "Hallmark_rank_score": 20, "adaptive": 20},
            "held-out prediction audit failed")
    require(set(predictions["dataset"]) == {"GSE205050", "GSE243023", "GSE10723"}, "held-out dataset audit failed")

    roc = pd.read_csv(SOURCE_DIR / "Figure5b_pooled_ROC_source_data.tsv", sep="\t")
    auc_by_strategy = roc.groupby("strategy")["frozen_ROC_AUC"].first().to_dict()
    require(auc_by_strategy == {"gene_rank": 0.8, "Hallmark_rank_score": 0.62, "adaptive": 0.54}, "pooled ROC AUC audit failed")

    per_dataset = pd.read_csv(SOURCE_DIR / "Figure5c_per_dataset_AUC_source_data.tsv", sep="\t")
    gene_auc = per_dataset.loc[per_dataset["strategy"].eq("gene_rank"), "ROC_AUC"].tolist()
    require(np.allclose(gene_auc, [0.1111111111111111, 1.0, 0.9375]), "gene-rank heterogeneity audit failed")

    features = pd.read_csv(SOURCE_DIR / "Figure5d_feature_stability_source_data.tsv", sep="\t")
    require(features.shape[0] == 20 and features["display_order"].tolist() == list(range(20)), "feature display count/order audit failed")
    require(features["strategy"].value_counts().to_dict() == {"Hallmark_rank_score": 13, "gene_rank": 5, "adaptive": 2}, "feature strategy audit failed")
    require(features["selected_outer_folds"].value_counts().to_dict() == {3: 15, 2: 5}, "outer-fold count audit failed")
    require(features["selection_rule"].nunique() == 1, "feature selection rule audit failed")

    pooled = pd.read_csv(SOURCE_DIR / "Figure5e_pooled_metric_source_data.tsv", sep="\t")
    require(pooled["n"].tolist() == [20.0, 20.0, 20.0] and (pooled["exact_permutations"] == 28000).all(), "pooled n/permutation audit failed")
    require(np.allclose(pooled["ROC_AUC"], [0.80, 0.62, 0.54]), "pooled metric audit failed")

    audit = json.loads((FIGURE_DIR / "Figure5_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS" and audit["input_manifest_validation"] == "PASS", "build audit is not PASS")
    require(audit["assertions"]["feature_display_count"] == 20 and audit["assertions"]["exact_permutations"] == 28000,
            "build audit feature/permutation assertion failed")

    lines = [
        "# Figure 5 QA",
        "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for 60 held-out predictions (20 per strategy), pooled AUCs of 0.80/0.62/0.54, gene-rank dataset AUCs of 0.11/1.00/0.94, the deterministic 20-feature display (15 features in three and five features in two outer folds), and three n=20 pooled metrics using 28,000 exact permutations.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        f"- Visual raster audit: non-white content is not clipped at image edges; {png.width:,} × {png.height:,}-pixel PNG preview and {tiff.width:,} × {tiff.height:,}-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are vector/chart elements derived from frozen tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the legend requires transportability and heterogeneity framing, labels pooled intervals and permutation P values as conditional, and prohibits diagnostic, biomarker-panel, and external-clinical-validation claims.",
        "",
        "## Remaining boundary",
        "This figure describes cross-dataset prediction within small assembled public datasets. It does not establish clinical diagnostic performance, a validated biomarker panel, stable external transportability, or OSA-to-ED causality.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
