#!/usr/bin/env python
"""Build Gate 04 v2 figures solely from frozen revision_v2 result tables."""
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "02_results"
WORK = ROOT / "01_work" / "GSE206528_doublet_sensitivity" / "attempt_02"
OUT = RESULTS / "figures"
SOURCE = RESULTS / "figure_source_data"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.titleweight"] = "bold"

COLORS = {
    "navy": "#1F4E79",
    "blue": "#4C78A8",
    "sky": "#86BBD8",
    "teal": "#2A9D8F",
    "green": "#4E9F6E",
    "orange": "#E9A03B",
    "red": "#C8553D",
    "purple": "#7158A3",
    "grey": "#777777",
    "lightgrey": "#E8E8E8",
    "dark": "#202124",
}


def add_panel_label(ax, label: str) -> None:
    ax.text(-0.075, 1.14, label, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left")


def clean_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_all(fig: plt.Figure, stem: str) -> None:
    for suffix, dpi in ((".svg", None), (".pdf", None), (".png", 600), (".tiff", 600)):
        path = OUT / f"{stem}{suffix}"
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if dpi is not None:
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)
    plt.close(fig)


def draw_box(ax, xy, width, height, title, body, color) -> None:
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0.012,rounding_size=0.025",
                           linewidth=1.3, edgecolor=color, facecolor="#FFFFFF")
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height * 0.68, title, ha="center", va="center",
            fontsize=6.8, fontweight="bold", color=color, linespacing=1.02)
    ax.text(x + width / 2, y + height * 0.34, body, ha="center", va="center",
            fontsize=5.8, color=COLORS["dark"], linespacing=1.06)


def make_figure1() -> None:
    hallmark = pd.read_csv(RESULTS / "Hallmark_clinical_bridge_v2.tsv", sep="\t")
    complete = hallmark.loc[hallmark["complete_directional_pattern"].astype(bool)].copy()
    complete["strict_clinical_bridge_pass"] = complete["strict_clinical_bridge_pass"].astype(bool)
    complete = complete.sort_values(
        ["strict_clinical_bridge_pass", "pathway"], ascending=[False, True]
    ).reset_index(drop=True)
    if len(complete) != 12 or complete["strict_clinical_bridge_pass"].sum() != 1:
        raise ValueError("Expected 12 complete-direction pathways and exactly one strict bridge.")

    contexts = [
        ("endo243", "Endothelial\nIH dataset 1"),
        ("endo205", "Endothelial\nIH dataset 2"),
        ("osaAM", "OSA vs control\n(AM)"),
        ("cpapAM", "CPAP vs OSA\n(AM)"),
    ]
    source_rows = []
    for _, row in complete.iterrows():
        for key, label in contexts:
            source_rows.append({
                "pathway": row["pathway"],
                "revised_bridge_class": row["revised_bridge_class"],
                "strict_clinical_bridge_pass": bool(row["strict_clinical_bridge_pass"]),
                "context": key,
                "context_label": label.replace("\n", " "),
                "NES": row[f"{key}_NES"],
                "FDR": row[f"{key}_FDR"],
                "FDR_lt_0_05_star": bool(row[f"{key}_FDR"] < 0.05),
            })
    heat_source = pd.DataFrame(source_rows)
    heat_source.to_csv(SOURCE / "Figure1_Hallmark_bridge_source_data_v2.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(7.20, 6.70))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.35])
    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.10, hspace=0.55, wspace=0.42)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    ax_d = fig.add_subplot(grid[1, :])

    # a. Parallel evidence, intentionally without arrows.
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis("off")
    draw_box(ax_a, (0.03, 0.66), 0.42, 0.20, "In-vitro IH",
             "controlled hypoxia", COLORS["blue"])
    draw_box(ax_a, (0.55, 0.66), 0.42, 0.20, "OSA/CPAP",
             "supportive clinical context", COLORS["orange"])
    draw_box(ax_a, (0.29, 0.39), 0.42, 0.20, "ED target tissue",
             "donor-pseudobulk", COLORS["green"])
    draw_box(ax_a, (0.22, 0.10), 0.56, 0.18, "Triangulation",
             "hypothesis generating", COLORS["purple"])
    for x0, y0, x1, y1 in [(0.24, 0.66, 0.48, 0.28), (0.76, 0.66, 0.52, 0.28), (0.50, 0.39, 0.50, 0.28)]:
        ax_a.plot([x0, x1], [y0, y1], color=COLORS["grey"], linewidth=1.1, zorder=0)
    ax_a.text(0.5, 0.96, "Parallel evidence layers", ha="center", fontsize=8.2, fontweight="bold")
    ax_a.text(0.5, -0.08, "Connections indicate triangulation, not causality.",
              ha="center", va="top", fontsize=7.0, style="italic", transform=ax_a.transAxes)
    add_panel_label(ax_a, "a")

    # b. Frozen evidence numbers, separated by context rather than a rank ladder.
    ax_b.axis("off")
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.text(0.02, 0.96, "Frozen v2 evidence labels", fontsize=8.5, fontweight="bold", va="top")
    cards = [
        (0.03, 0.55, "Endothelial IH", "1,035", "direction-stable\nscreen", COLORS["blue"]),
        (0.52, 0.55, "Endothelial IH", "258", "strict\nLODO-FDR", COLORS["navy"]),
        (0.03, 0.14, "IH–ED candidates", "3", "ED-FDR\nsupported", COLORS["green"]),
        (0.52, 0.14, "IH–ED candidates", "563", "nominal-P\nexploratory", COLORS["orange"]),
    ]
    for x, y, group, count, label, color in cards:
        ax_b.add_patch(FancyBboxPatch((x, y), 0.43, 0.25, boxstyle="round,pad=0.012",
                                      linewidth=0.8, edgecolor=color, facecolor="#FFFFFF"))
        ax_b.text(x + 0.03, y + 0.195, group, fontsize=5.8, color=COLORS["grey"])
        ax_b.text(x + 0.03, y + 0.095, count, fontsize=14, fontweight="bold", color=color)
        ax_b.text(x + 0.29, y + 0.048, label, fontsize=5.0, va="center", ha="center", linespacing=0.95)
    ax_b.text(0.02, 0.02, "Counts label distinct evidence filters; they are not a causal sequence.",
              fontsize=6.8, style="italic")
    add_panel_label(ax_b, "b")

    # c. Predeclared strict bridge rule.
    ax_c.axis("off")
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.text(0.00, 0.96, "Predeclared clinical-bridge screen", fontsize=8.3,
              fontweight="bold", va="bottom")
    rules = [
        ("50", "Hallmark pathways evaluated", COLORS["lightgrey"]),
        ("12", "complete directional pattern", COLORS["sky"]),
        ("1", "all four FDR < 0.05", COLORS["green"]),
    ]
    ys = [0.73, 0.46, 0.19]
    for (count, label, color), y in zip(rules, ys):
        ax_c.add_patch(FancyBboxPatch((0.04, y), 0.18, 0.17, boxstyle="round,pad=0.01",
                                      facecolor=color, edgecolor=COLORS["dark"], linewidth=0.6))
        ax_c.text(0.13, y + 0.085, count, ha="center", va="center", fontsize=11, fontweight="bold")
        ax_c.text(0.29, y + 0.085, label, va="center", fontsize=6.7)
        if y != ys[-1]:
            ax_c.plot([0.13, 0.13], [y - 0.06, y], color=COLORS["grey"], linewidth=1.0)
    ax_c.text(0.04, 0.03, "1 strict pathway: HALLMARK_HYPOXIA", fontsize=5.5,
              color=COLORS["green"], fontweight="bold")
    add_panel_label(ax_c, "c")

    # d. NES heatmap with vector FDR stars and strict bridge outline.
    mat = np.array([[row[f"{key}_NES"] for key, _ in contexts] for _, row in complete.iterrows()])
    max_abs = max(2.6, float(np.abs(mat).max()))
    im = ax_d.imshow(mat, aspect="auto", cmap="RdBu_r",
                     norm=TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs))
    ax_d.set_xticks(range(len(contexts)))
    ax_d.set_xticklabels([label for _, label in contexts], fontsize=7.0)
    ax_d.set_yticks(range(len(complete)))
    ax_d.set_yticklabels(complete["pathway"].str.replace("HALLMARK_", "", regex=False)
                         .str.replace("_", " ", regex=False), fontsize=6.8)
    ax_d.tick_params(length=0)
    for i, (_, row) in enumerate(complete.iterrows()):
        for j, (key, _) in enumerate(contexts):
            if row[f"{key}_FDR"] < 0.05:
                ax_d.text(j, i, "*", ha="center", va="center", fontsize=10.5, color="black", fontweight="bold")
    strict_i = int(complete.index[complete["strict_clinical_bridge_pass"]][0])
    ax_d.add_patch(Rectangle((-0.5, strict_i - 0.5), len(contexts), 1,
                             fill=False, edgecolor=COLORS["green"], linewidth=2.0, clip_on=False))
    cbar = fig.colorbar(im, ax=ax_d, fraction=0.026, pad=0.012)
    cbar.set_label("NES", fontsize=7.2)
    cbar.ax.tick_params(labelsize=6.5)
    ax_d.set_title("Complete directional-pattern pathways", fontsize=9.0, pad=7)
    ax_d.text(0.01, -0.16, "* FDR < 0.05; green outline = only strict bridge", transform=ax_d.transAxes,
              fontsize=6.5, va="top")
    add_panel_label(ax_d, "d")
    save_all(fig, "Figure1_parallel_triangulation_v2")


def make_figure3() -> None:
    candidates = pd.read_csv(RESULTS / "IH_ED_candidate_classification_v2.tsv", sep="\t")
    summary = pd.read_csv(RESULTS / "IH_ED_candidate_classification_summary_v2.tsv", sep="\t")
    sensitivity = pd.read_csv(WORK / "integration" / "TYMS_EFNB2_LRRC17_integrated_sensitivity.tsv", sep="\t")
    effects = pd.read_csv(WORK / "edgeR_Hallmark" / "TYMS_EFNB2_LRRC17_sensitivity.tsv", sep="\t")
    supported = candidates.loc[candidates["revised_candidate_class"] == "ED_FDR_SUPPORTED"].copy()
    supported = supported.sort_values("gene_symbol").reset_index(drop=True)
    if len(supported) != 3 or set(supported["gene_symbol"]) != {"TYMS", "EFNB2", "LRRC17"}:
        raise ValueError("Expected three named ED-FDR-supported candidates.")
    if dict(zip(summary["revised_candidate_class"], summary["candidate_genes"])) != {
        "ED_FDR_SUPPORTED": 3, "NOMINAL_P_EXPLORATORY": 563
    }:
        raise ValueError("Candidate class summary differs from frozen v2 definition.")
    selected = effects.loc[effects["file"] == "differential_fibroblast_nonDM_ED_vs_normal.tsv"].copy()
    selected = selected.set_index("gene_symbol").loc[["TYMS", "EFNB2", "LRRC17"]].reset_index()
    if not (selected["new_FDR"] < 0.05).all():
        raise ValueError("All three candidate sensitivity effects must retain FDR support in fibroblast.")

    sensitivity.to_csv(SOURCE / "Figure3_candidate_sensitivity_class_source_data_v2.tsv", sep="\t", index=False)
    selected.to_csv(SOURCE / "Figure3_fibroblast_effect_source_data_v2.tsv", sep="\t", index=False)
    supported[["gene_symbol", "revised_candidate_class", "supporting_cell_types", "celltype_count_is_independent_evidence"]].to_csv(
        SOURCE / "Figure3_descriptive_celltype_map_source_data_v2.tsv", sep="\t", index=False
    )
    summary.to_csv(SOURCE / "Figure3_candidate_class_count_source_data_v2.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(7.20, 3.25), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[0.78, 1.30, 1.22])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    # a. Evidence strata, deliberately no ordering within either stratum.
    count_map = dict(zip(summary["revised_candidate_class"], summary["candidate_genes"]))
    labels = ["ED-FDR\nsupported", "Nominal-P\nexploratory"]
    values = [count_map["ED_FDR_SUPPORTED"], count_map["NOMINAL_P_EXPLORATORY"]]
    bars = ax_a.bar(labels, values, color=[COLORS["green"], COLORS["orange"]], edgecolor=COLORS["dark"], linewidth=0.7)
    ax_a.set_yscale("log")
    ax_a.set_ylim(1, 1000)
    ax_a.set_ylabel("Candidate genes (log scale)", fontsize=8.5)
    ax_a.tick_params(axis="x", labelsize=7.7)
    ax_a.tick_params(axis="y", labelsize=7.2)
    for bar, value in zip(bars, values):
        ax_a.text(bar.get_x() + bar.get_width() / 2, value * 1.25, str(value), ha="center", va="bottom",
                  fontsize=9.8, fontweight="bold")
    ax_a.set_title("Evidence strata", fontsize=9.8, pad=7)
    ax_a.text(0.5, -0.28, "No within-stratum ranking", transform=ax_a.transAxes, ha="center",
              fontsize=6.8, style="italic")
    clean_axis(ax_a)
    add_panel_label(ax_a, "a")

    # b. Descriptive candidate-to-cell-type map; edges have no direction/weight.
    ax_b.axis("off")
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    gene_y = {"TYMS": 0.76, "EFNB2": 0.50, "LRRC17": 0.24}
    cell_y = {"fibroblast": 0.82, "pericyte": 0.62, "smooth_muscle": 0.39, "vascular_endothelial": 0.17}
    for gene, y in gene_y.items():
        ax_b.add_patch(FancyBboxPatch((0.04, y - 0.065), 0.28, 0.13, boxstyle="round,pad=0.015",
                                      facecolor="#EAF4EC", edgecolor=COLORS["green"], linewidth=1.0))
        ax_b.text(0.18, y, gene, ha="center", va="center", fontsize=9.2, fontweight="bold")
    for cell, y in cell_y.items():
        ax_b.add_patch(FancyBboxPatch((0.68, y - 0.055), 0.28, 0.11, boxstyle="round,pad=0.012",
                                      facecolor="#F4F4F4", edgecolor=COLORS["grey"], linewidth=0.8))
        ax_b.text(0.82, y, cell.replace("_", " "), ha="center", va="center", fontsize=7.0)
    for _, row in supported.iterrows():
        gene = row["gene_symbol"]
        for cell in str(row["supporting_cell_types"]).split(";"):
            ax_b.plot([0.32, 0.68], [gene_y[gene], cell_y[cell]], color=COLORS["grey"], linewidth=1.1, zorder=0)
    ax_b.text(0.50, 1.01, "Descriptive concordant cell-type readouts", ha="center", fontsize=9.4, fontweight="bold")
    ax_b.text(0.50, -0.10, "Lines are non-causal; cell types are nested donor-pseudobulk readouts.",
              ha="center", fontsize=6.5, style="italic", transform=ax_b.transAxes)
    add_panel_label(ax_b, "b")

    # c. Donor-level effects pre/post doublet exclusion; FDR shown explicitly.
    ordered = selected.set_index("gene_symbol").loc[["TYMS", "EFNB2", "LRRC17"]].reset_index()
    ypos = np.arange(len(ordered))[::-1]
    for y, (_, row) in zip(ypos, ordered.iterrows()):
        color = COLORS["blue"] if row["new_logFC"] < 0 else COLORS["red"]
        ax_c.plot([row["old_logFC"], row["new_logFC"]], [y, y], color=COLORS["grey"], linewidth=1.4, zorder=1)
        ax_c.scatter(row["old_logFC"], y, s=38, color="white", edgecolor=COLORS["dark"], linewidth=0.8, zorder=2, label="Historical" if y == ypos[0] else None)
        ax_c.scatter(row["new_logFC"], y, s=42, color=color, edgecolor="white", linewidth=0.5, zorder=3, label="Doublet excluded" if y == ypos[0] else None)
        ax_c.text(2.52, y, f"FDR {row['old_FDR']:.3f} → {row['new_FDR']:.3f}", va="center", fontsize=6.8)
    ax_c.axvline(0, color=COLORS["grey"], linestyle="--", linewidth=0.8, zorder=0)
    ax_c.set_yticks(ypos)
    ax_c.set_yticklabels(ordered["gene_symbol"], fontsize=8.4, fontweight="bold")
    ax_c.set_xlim(-2.8, 3.9)
    ax_c.set_xlabel("Donor-pseudobulk log2FC\n(non-DM ED vs normal, fibroblast)", fontsize=7.7)
    ax_c.set_title("Doublet sensitivity", fontsize=9.8, pad=7)
    ax_c.legend(loc="lower left", fontsize=6.7, frameon=False, handletextpad=0.4)
    ax_c.tick_params(axis="x", labelsize=7.0)
    ax_c.tick_params(axis="y", length=0)
    clean_axis(ax_c)
    ax_c.text(0.0, -0.29, "All three retain FDR < 0.05 after exclusion.", transform=ax_c.transAxes,
              fontsize=6.5, style="italic")
    add_panel_label(ax_c, "c")
    save_all(fig, "Figure3_candidate_sensitivity_v2")


if __name__ == "__main__":
    make_figure1()
    make_figure3()
    print(f"Wrote figures to {OUT}")
