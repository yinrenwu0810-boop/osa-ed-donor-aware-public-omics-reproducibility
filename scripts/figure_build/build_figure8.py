from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 7, "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": .8, "legend.frameon": False,
})
OUTPUT_ROOT = Path(__file__).resolve().parents[1]
ROOT = OUTPUT_ROOT.parent
OUT = OUTPUT_ROOT / "figures" / "Figure8_20260912"
SOURCE_DIR = OUT / "source_data"
MANIFEST = OUTPUT_ROOT / "figure_contracts" / "input_manifest.sha256.tsv"
SOURCES = {
    "S080": "analysis/results/phase3/cell_communication/all_tested_ligand_receptor_edges.tsv",
    "S081": "analysis/results/phase3/cell_communication/analysis_summary.json",
    "S082": "analysis/results/phase3/power_sensitivity/minimum_detectable_effects.tsv",
    "S083": "analysis/results/phase3/power_sensitivity/sample_size_for_standardized_effects.tsv",
    "S084": "revision_v2/07_virtual_perturbation/08_external_reference/public_perturbation_target_coverage.tsv",
    "S085": "revision_v2/07_virtual_perturbation/08_external_reference/VP_G07_external_scope_report.md",
}
TARGETS = ["TYMS", "EFNB2", "LRRC17"]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_inputs() -> dict[str, str]:
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str).set_index("source_id")
    out = {}
    for source_id, relative_path in SOURCES.items():
        path = ROOT / relative_path
        actual = digest(path)
        if manifest.loc[source_id, "relative_path"] != relative_path or actual.lower() != manifest.loc[source_id, "sha256"].lower():
            raise ValueError(f"input drift: {source_id}")
        out[source_id] = actual
    return out


def panel_label(ax, label: str, x: float = -0.14) -> None:
    ax.text(x, 1.05, label, transform=ax.transAxes, fontsize=8.2, fontweight="bold")


def add_box(ax, xy, width, height, text, face, edge="#555555", fontsize=4.1):
    box = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.012,rounding_size=0.02", facecolor=face, edgecolor=edge, linewidth=.7)
    ax.add_patch(box)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=fontsize, wrap=True)


def main() -> None:
    hashes = verify_inputs()
    edges = pd.read_csv(ROOT / SOURCES["S080"], sep="\t", usecols=["ligand", "receptor", "source_cell", "target_cell", "delta_ED_minus_normal", "exact_permutation_p"])
    summary = json.loads((ROOT / SOURCES["S081"]).read_text(encoding="utf-8"))
    mde = pd.read_csv(ROOT / SOURCES["S082"], sep="\t")
    sample = pd.read_csv(ROOT / SOURCES["S083"], sep="\t")
    coverage = pd.read_csv(ROOT / SOURCES["S084"], sep="\t")
    scope = (ROOT / SOURCES["S085"]).read_text(encoding="utf-8")
    if edges.shape[0] != int(summary["tested_cell_specific_edges"]) or edges.shape[0] != 60299:
        raise ValueError("communication edge count failed")
    if not (edges["exact_permutation_p"] >= .1).all() or summary["minimum_exact_two_sided_p_with_3_vs_3"] != .1:
        raise ValueError("exact-test resolution boundary failed")
    if summary["nominal_p_lt_0_05"] != 0 or summary["FDR_lt_0_05"] != 0:
        raise ValueError("communication significance boundary failed")
    design = mde.loc[mde["analysis"].eq("GSE206528_nonDM_ED_vs_normal")].sort_values("target_power")
    if not np.allclose(design["minimum_detectable_cohens_d"], [3.070892269599019, 3.5892055765349036]):
        raise ValueError("n=3-vs-3 MDE audit failed")
    if not {"NO_VERIFIED_HIT", "VERIFIED_DATASET_MEMBER"}.issubset(set(coverage["evidence_status"])) or "没有任何记录为人海绵体成纤维细胞" not in scope:
        raise ValueError("external coverage boundary failed")

    edges = edges.assign(abs_delta=edges["delta_ED_minus_normal"].abs())
    top = edges.sort_values(["abs_delta", "ligand", "receptor", "source_cell", "target_cell"], ascending=[False, True, True, True, True]).head(15).copy()
    top["edge_label"] = [f"{i + 1}. {ligand}–{receptor}" for i, (ligand, receptor) in enumerate(zip(top["ligand"], top["receptor"]))]
    hist_counts, hist_bins = np.histogram(edges["delta_ED_minus_normal"], bins=60)
    hist = pd.DataFrame({"bin_left": hist_bins[:-1], "bin_right": hist_bins[1:], "count": hist_counts})
    allocations = pd.DataFrame({"allocation_rank": range(1, 21), "tail_class": ["extreme tail"] + ["non-extreme"] * 18 + ["extreme tail"]})
    external = coverage.loc[coverage["evidence_status"].isin(["VERIFIED_DATASET_MEMBER", "NO_VERIFIED_HIT"])].copy()
    external = external.drop_duplicates(["accession", "target_gene", "cell_context_tier"])
    tiers = ["T1_EXACT", "T2_RELATED_HUMAN", "T3_OTHER_PRIMARY_HUMAN", "T4_HUMAN_CELL_LINE", "T5"]
    matrix = []
    for target in TARGETS:
        for tier in tiers:
            q = external.loc[(external["target_gene"].eq(target)) & (external["cell_context_tier"].eq(tier))]
            if tier == "T1_EXACT":
                state, label = "NO_VERIFIED_T1_EXACT_RECORD", "NO\nRECORD"
            elif q.empty:
                state, label = "NO_VERIFIED_RECORD", "—"
            else:
                state = "VERIFIED_COVERAGE"
                modalities = q["perturbation_type"].str.replace(" knockdown", " KD", regex=False).str.replace("adenoviral ", "Ad-", regex=False).tolist()
                label = "\n".join(modalities)
            matrix.append({"target_gene": target, "tier": tier, "state": state, "display_label": label, "record_count": int(q.shape[0])})
    matrix = pd.DataFrame(matrix)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    for frame, name in [
        (hist, "Figure8a_effect_distribution_source_data.tsv"), (top, "Figure8a_top15_absolute_delta_source_data.tsv"),
        (allocations, "Figure8b_exact_allocation_schematic_source_data.tsv"), (design, "Figure8c_n3v3_MDE_source_data.tsv"),
        (sample, "Figure8c_idealized_sample_size_source_data.tsv"), (external, "Figure8d_external_coverage_records_source_data.tsv"),
        (matrix, "Figure8d_coverage_matrix_source_data.tsv"), (matrix, "Figure8e_gap_map_source_data.tsv"),
    ]:
        frame.to_csv(SOURCE_DIR / name, sep="\t", index=False)

    fig = plt.figure(figsize=(7.205, 5.709))
    grid = fig.add_gridspec(2, 6, height_ratios=[1, 1.12])
    fig.subplots_adjust(left=.075, right=.98, top=.94, bottom=.10, wspace=.85, hspace=.86)
    a = fig.add_subplot(grid[0, :3]); b = fig.add_subplot(grid[0, 3:5]); c = fig.add_subplot(grid[0, 5])
    d = fig.add_subplot(grid[1, :4]); e = fig.add_subplot(grid[1, 4:])

    # a: full distribution + mechanically selected top 15 inset
    a.hist(edges["delta_ED_minus_normal"], bins=60, color="#8FB3D2", edgecolor="white", linewidth=.15)
    a.axvline(0, color="#555555", lw=.7)
    a.set_xlabel("ED minus normal compatibility-score delta", fontsize=4.7); a.set_ylabel("Tested edges", fontsize=4.7); a.tick_params(labelsize=4.2)
    a.set_title("All communication edges: effect sizes without significant P values", loc="left", fontsize=6.3, fontweight="bold")
    a.text(.02, .94, f"{len(edges):,} tested edges\nall nominal P ≥ 0.10; FDR ≥ 0.05", transform=a.transAxes, fontsize=4.0, va="top", color="#555555")
    inset = a.inset_axes([.57, .18, .41, .66])
    order = top.iloc[::-1]
    bar_colors = np.where(order["delta_ED_minus_normal"] >= 0, "#C76E62", "#5C90B8")
    inset.barh(np.arange(15), order["delta_ED_minus_normal"], color=bar_colors, height=.65)
    inset.axvline(0, color="#555555", lw=.55); inset.set_yticks(np.arange(15)); inset.set_yticklabels(order["edge_label"], fontsize=3.25)
    inset.tick_params(axis="x", labelsize=3.1); inset.set_title("Top 15 |Δ| (P/direction not used for selection)", fontsize=3.55, loc="left")
    inset.spines[["right", "top"]].set_visible(False)
    panel_label(a, "a", x=-.09)

    # b: exact combinatorial resolution
    x = np.arange(1, 21)
    colors = np.array(["#D65F5F"] + ["#D9D9D9"] * 18 + ["#D65F5F"])
    b.hlines(0, 1, 20, color="#9E9E9E", lw=.8); b.scatter(x, np.zeros(20), s=22, color=colors, edgecolor="#666666", linewidth=.25, zorder=3)
    b.annotate("extreme\nallocation", xy=(1, 0), xytext=(3.2, .42), arrowprops={"arrowstyle": "-", "lw": .5}, fontsize=3.8, ha="center")
    b.annotate("extreme\nallocation", xy=(20, 0), xytext=(17.8, .42), arrowprops={"arrowstyle": "-", "lw": .5}, fontsize=3.8, ha="center")
    b.text(10.5, -.42, "2 extreme tails / 20 label allocations\nminimum two-sided exact P = 0.10", ha="center", fontsize=4.5, fontweight="bold")
    b.text(10.5, -.74, "3 normal vs 3 organic-ED donors; no simulated P values", ha="center", fontsize=3.6, color="#666666")
    b.set_xlim(.2, 20.8); b.set_ylim(-.95, .82); b.axis("off"); b.set_title("Why P = 0.10 is the attainable two-sided minimum", loc="left", fontsize=6.1, fontweight="bold")
    panel_label(b, "b", x=-.13)

    # c: n=3-vs-3 MDE points plus idealized equal-group planning curves
    for power, color in [(0.8, "#739E82"), (0.9, "#406D5A")]:
        q = sample.loc[sample["target_power"].eq(power)].sort_values("cohens_d")
        c.plot(q["cohens_d"], q["n_per_group_equal_design"], marker="o", ms=2.5, lw=.85, color=color, label=f"{int(power*100)}%")
    c.axhline(3, color="#888888", lw=.55, ls="--")
    for power, value, color in zip(design["target_power"], design["minimum_detectable_cohens_d"], ["#739E82", "#406D5A"]):
        c.scatter(value, 3, s=20, color=color, zorder=3)
        c.text(value, 7, f"d={value:.2f}", fontsize=3.3, ha="center", rotation=55)
    c.set_xlim(.4, 3.9); c.set_ylim(0, 92); c.set_xlabel("Cohen's d", fontsize=4.1); c.set_ylabel("Idealized n/group", fontsize=4.1); c.tick_params(labelsize=3.5)
    c.set_title("n=3 vs 3\nresolution + planning", loc="left", fontsize=5.3, fontweight="bold")
    c.legend(fontsize=3.4, loc="upper right", title="power", title_fontsize=3.3, handlelength=.9)
    c.text(.02, .04, "MDE points at n=3; α=0.05", transform=c.transAxes, fontsize=3.1, color="#666666")
    panel_label(c, "c", x=-.34)

    # d: T1–T5 coverage matrix
    state_colors = {"VERIFIED_COVERAGE": "#8FB3D2", "NO_VERIFIED_T1_EXACT_RECORD": "#F2C7C7", "NO_VERIFIED_RECORD": "#EEEEEE"}
    short_tier = {"T1_EXACT": "T1\nexact", "T2_RELATED_HUMAN": "T2\nrelated", "T3_OTHER_PRIMARY_HUMAN": "T3\nprimary", "T4_HUMAN_CELL_LINE": "T4\ncell line", "T5": "T5"}
    for yi, target in enumerate(TARGETS):
        for xi, tier in enumerate(tiers):
            row = matrix.loc[(matrix["target_gene"].eq(target)) & (matrix["tier"].eq(tier))].iloc[0]
            d.add_patch(plt.Rectangle((xi, yi), 1, 1, facecolor=state_colors[row["state"]], edgecolor="white", linewidth=1.2))
            if row["state"] == "VERIFIED_COVERAGE":
                d.text(xi+.5, yi+.5, row["display_label"], ha="center", va="center", fontsize=3.1, wrap=True)
            elif tier == "T1_EXACT":
                d.text(xi+.5, yi+.5, "NO\nRECORD", ha="center", va="center", fontsize=3.7, color="#8B1A1A", fontweight="bold")
            else:
                d.text(xi+.5, yi+.5, "—", ha="center", va="center", fontsize=5.5, color="#888888")
    d.set_xlim(0, 5); d.set_ylim(3, 0); d.set_xticks(np.arange(5)+.5); d.set_xticklabels([short_tier[t] for t in tiers], fontsize=4.6); d.set_yticks(np.arange(3)+.5); d.set_yticklabels(["TYMS", "EFNB2", "LRRC17\n(no verified hit)"], fontsize=4.5)
    d.tick_params(length=0); d.set_title("Frozen external perturbation coverage: records are not validation", loc="left", fontsize=6.2, fontweight="bold")
    d.text(0, 1.09, "Blue: verified external record; red: no verified T1 exact record in frozen search scope; gray: no tiered record.", transform=d.transAxes, fontsize=3.7, color="#555555")
    panel_label(d, "d", x=-.06)

    # e: concise gap map
    e.axis("off"); e.set_title("Direct validation gap remains", loc="left", fontsize=6.2, fontweight="bold")
    add_box(e, (.03, .54), .94, .27, "Required exact context (T1): human cavernosal fibroblast genetic perturbation\nNO VERIFIED RECORD in the frozen search scope", "#F4CCCC", edge="#B24C4C", fontsize=4.25)
    add_box(e, (.03, .18), .45, .23, "Verified but non-exact coverage\nEFNB2: T2, T3, T4\nTYMS: T3, T4", "#D9E8F3", edge="#5A86A8", fontsize=4.0)
    add_box(e, (.52, .18), .45, .23, "Unresolved target\nLRRC17: NO_VERIFIED_HIT\n(bounded scope; not global absence)", "#EEEEEE", edge="#777777", fontsize=3.85)
    e.text(.5, .05, "Coverage/background evidence only; G08 was unexecuted and is excluded.", transform=e.transAxes, ha="center", fontsize=3.6, color="#555555")
    panel_label(e, "e", x=-.10)

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "Figure8_statistical_resolution_and_external_gap"
    for ext, kwargs in [(".svg", {}), (".pdf", {}), (".png", {"dpi": 300}), (".tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
        fig.savefig(base.with_suffix(ext), facecolor="white", **kwargs)
    plt.close(fig)
    session = {"python": sys.version, "platform": platform.platform(), "matplotlib": matplotlib.__version__, "pandas": pd.__version__, "numpy": np.__version__}
    (OUT / "Figure8_session_record.json").write_text(json.dumps(session, indent=2) + "\n", encoding="utf-8")
    audit = {
        "status": "PASS", "figure": "Figure 8", "backend": "Python/matplotlib", "input_manifest_validation": "PASS", "input_sha256": hashes,
        "assertions": {"tested_edges": int(edges.shape[0]), "all_nominal_P_ge_0_1": True, "nominal_P_lt_0_05": 0, "FDR_lt_0_05": 0,
                       "exact_label_allocations": 20, "two_extreme_tails": 2, "n3v3_MDE_at_80pct": float(design.iloc[0]["minimum_detectable_cohens_d"]),
                       "n3v3_MDE_at_90pct": float(design.iloc[1]["minimum_detectable_cohens_d"]), "verified_external_records": int(external.loc[external["evidence_status"].eq("VERIFIED_DATASET_MEMBER")].shape[0]),
                       "T1_exact_records": 0, "G08_executed": False},
        "output_dimensions_mm": [183, 145],
        "interpretation_boundary": "Exact-test resolution, idealized planning calculations, and external coverage prioritize future work only; external records are not independent biological validation and NO_VERIFIED_T1_EXACT_RECORD is bounded to the frozen search scope.",
    }
    (OUT / "Figure8_build_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}))


if __name__ == "__main__":
    main()
