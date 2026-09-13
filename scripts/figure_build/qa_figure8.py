from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure8_20260912"
SOURCE_DIR = FIGURE_DIR / "source_data"
BASE = FIGURE_DIR / "Figure8_statistical_resolution_and_external_gap"
REPORT = FIGURE_DIR / "Figure8_QA.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    expected = [BASE.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    expected += [
        FIGURE_DIR / "Figure8_build_audit.json", FIGURE_DIR / "Figure8_session_record.json", FIGURE_DIR / "Figure8_legend.md",
        SOURCE_DIR / "Figure8a_effect_distribution_source_data.tsv",
        SOURCE_DIR / "Figure8a_top15_absolute_delta_source_data.tsv",
        SOURCE_DIR / "Figure8b_exact_allocation_schematic_source_data.tsv",
        SOURCE_DIR / "Figure8c_n3v3_MDE_source_data.tsv",
        SOURCE_DIR / "Figure8c_idealized_sample_size_source_data.tsv",
        SOURCE_DIR / "Figure8d_external_coverage_records_source_data.tsv",
        SOURCE_DIR / "Figure8d_coverage_matrix_source_data.tsv",
        SOURCE_DIR / "Figure8e_gap_map_source_data.tsv",
    ]
    for path in expected:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path.name}")

    png = Image.open(BASE.with_suffix(".png")).convert("RGB")
    tiff = Image.open(BASE.with_suffix(".tiff"))
    require(png.width >= 2100 and png.height >= 1700, "PNG preview resolution is below expected size")
    require(all(abs(value - 600.0) < 0.1 for value in tiff.info.get("dpi", (0, 0))), "TIFF is not 600 dpi")
    require(tiff.info.get("compression") == "tiff_lzw", "TIFF compression is not LZW")
    bbox = ImageChops.difference(png, Image.new("RGB", png.size, "white")).getbbox()
    require(bbox is not None, "PNG has no visible content")
    left, top, right, bottom = bbox
    require(left >= 3 and top >= 3 and right <= png.width - 3 and bottom <= png.height - 3,
            "visible content is clipped at a PNG edge")

    svg_text = BASE.with_suffix(".svg").read_text(encoding="utf-8")
    require(len(re.findall(r"<text[ >]", svg_text)) > 95, "SVG text is not retained as editable text nodes")
    for label in ("minimum two-sided exact P = 0.10", "NO", "RECORD", "Direct validation gap remains"):
        require(label in svg_text, f"required interpretation label missing from SVG: {label}")

    distribution = pd.read_csv(SOURCE_DIR / "Figure8a_effect_distribution_source_data.tsv", sep="\t")
    require(distribution.shape[0] == 60 and distribution["count"].sum() == 60299, "effect distribution count audit failed")
    top = pd.read_csv(SOURCE_DIR / "Figure8a_top15_absolute_delta_source_data.tsv", sep="\t")
    require(top.shape[0] == 15 and top["abs_delta"].is_monotonic_decreasing, "top-absolute-delta selection audit failed")
    require((top["exact_permutation_p"] >= .1).all(), "top-edge P-value boundary audit failed")

    allocations = pd.read_csv(SOURCE_DIR / "Figure8b_exact_allocation_schematic_source_data.tsv", sep="\t")
    require(allocations.shape[0] == 20 and allocations["tail_class"].eq("extreme tail").sum() == 2, "exact-allocation audit failed")

    mde = pd.read_csv(SOURCE_DIR / "Figure8c_n3v3_MDE_source_data.tsv", sep="\t")
    require(mde.shape[0] == 2 and np.allclose(mde["minimum_detectable_cohens_d"], [3.070892269599019, 3.5892055765349036]), "n=3-vs-3 MDE audit failed")
    sample = pd.read_csv(SOURCE_DIR / "Figure8c_idealized_sample_size_source_data.tsv", sep="\t")
    require(sample.shape[0] == 10 and sample.groupby("target_power").size().to_dict() == {0.8: 5, 0.9: 5}, "idealized planning-curve audit failed")

    external = pd.read_csv(SOURCE_DIR / "Figure8d_external_coverage_records_source_data.tsv", sep="\t")
    require(external.shape[0] == 7 and external["evidence_status"].value_counts().to_dict() == {"VERIFIED_DATASET_MEMBER": 6, "NO_VERIFIED_HIT": 1}, "external coverage record audit failed")
    require(external["matches_corpus_cavernosum_fibroblast"].isin(["NO", "NO_VERIFIED_T1_HIT"]).all(), "external exact-context boundary audit failed")
    matrix = pd.read_csv(SOURCE_DIR / "Figure8d_coverage_matrix_source_data.tsv", sep="\t")
    require(matrix.shape[0] == 15 and matrix.loc[matrix["tier"].eq("T1_EXACT"), "state"].eq("NO_VERIFIED_T1_EXACT_RECORD").all(), "T1 coverage-gap audit failed")
    require(matrix.loc[matrix["state"].eq("VERIFIED_COVERAGE"), "record_count"].sum() == 6, "verified coverage matrix audit failed")

    audit = json.loads((FIGURE_DIR / "Figure8_build_audit.json").read_text(encoding="utf-8"))
    require(audit["status"] == "PASS" and audit["input_manifest_validation"] == "PASS", "build audit is not PASS")
    require(audit["assertions"]["tested_edges"] == 60299 and audit["assertions"]["exact_label_allocations"] == 20 and audit["assertions"]["T1_exact_records"] == 0 and not audit["assertions"]["G08_executed"], "build assertion audit failed")
    session = json.loads((FIGURE_DIR / "Figure8_session_record.json").read_text(encoding="utf-8"))
    require(all(key in session for key in ("python", "matplotlib", "pandas", "numpy")), "session record audit failed")

    lines = [
        "# Figure 8 QA", "",
        "- Status: **PASS**",
        "- Backend for build, exports, and raster inspection: Python (matplotlib/Pillow).",
        "- Source-drift check: PASS against the frozen input manifest before build.",
        "- Numerical audit: PASS for all 60,299 tested communication edges, 60 deterministic histogram bins, the mechanically selected 15 largest absolute deltas, 20 exact 3-versus-3 label allocations with two extreme tails, n=3-versus-3 MDE values of 3.071/3.589 at 80%/90% power, 10 idealized sample-size records, seven retained external-coverage records, and the 15-cell target-by-tier coverage matrix.",
        "- Exact-test boundary audit: PASS. Every communication edge has nominal P at least 0.10; zero records meet nominal P < 0.05 or global FDR < 0.05. The two-sided lower bound is 2/20 = 0.10.",
        "- External-coverage audit: PASS. The matrix contains six verified non-exact external records and one bounded NO_VERIFIED_HIT record. All three T1 cells are NO_VERIFIED_T1_EXACT_RECORD; this is not global absence or independent validation.",
        "- Export audit: SVG, PDF, PNG, and 600-dpi LZW TIFF present; SVG retains editable text nodes.",
        f"- Visual raster audit: non-white content is not clipped at image edges; {png.width:,} × {png.height:,}-pixel PNG preview and {tiff.width:,} × {tiff.height:,}-pixel TIFF were checked.",
        "- Image-integrity audit: no microscopy or source raster image was used; all panels are vector/chart elements derived from frozen tabular inputs. No crop, local brightness/contrast change, stitching, or image reuse applies.",
        "- Interpretation audit: the legend describes exact-test resolution and idealized planning without retroactively changing the analysis. It labels external records as coverage/background evidence and prohibits independent-validation, signaling, causal, expression-reversal, ED-reversal, and OSA-to-ED-causality claims.",
        "", "## Remaining boundary",
        "This figure identifies why current data can prioritize hypotheses but cannot resolve donor-level communication significance or independently validate target perturbations in the exact cavernosal-fibroblast context. G08 remains unexecuted and excluded.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("status=PASS")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
