from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = OUTPUT_ROOT.parent
OUT = OUTPUT_ROOT / "figures" / "Figure4_20260912"
SOURCE_OUT = OUT / "source_data"
INPUT_MANIFEST = OUTPUT_ROOT / "figure_contracts" / "input_manifest.sha256.tsv"
SOURCE_IDS = {
    "S034": "analysis/results/integration/IH_ED_celltype_gene_evidence.tsv",
    "S040": "revision_v2/02_results/IH_ED_candidate_classification_summary_v2.tsv",
    "S041": "revision_v2/02_results/IH_ED_candidate_classification_v2.tsv",
    "S042": "revision_v2/01_work/GSE206528_doublet_sensitivity/attempt_02/edgeR_Hallmark/TYMS_EFNB2_LRRC17_sensitivity.tsv",
    "S043": "revision_v2/01_work/GSE206528_doublet_sensitivity/attempt_02/edgeR_Hallmark/Hallmark_sensitivity_comparison.tsv",
    "S044": "revision_v2/02_results/figure_source_data/Figure3_candidate_sensitivity_class_source_data_v2.tsv",
}
GENE_ORDER = ["TYMS", "EFNB2", "LRRC17"]
GENE_COLORS = {"TYMS": "#2A9D8F", "EFNB2": "#E69F00", "LRRC17": "#9C6FB6"}
STATE_COLORS = {
    "Retained significant": "#3D76B7", "Gained significance": "#009E73",
    "Lost significance": "#C65B50", "Retained non-significant": "#BEBEBE",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inputs() -> dict[str, str]:
    manifest = pd.read_csv(INPUT_MANIFEST, sep="\t", dtype=str, keep_default_na=False).set_index("source_id")
    verified: dict[str, str] = {}
    for source_id, relative_path in SOURCE_IDS.items():
        if source_id not in manifest.index or manifest.loc[source_id, "relative_path"] != relative_path:
            raise ValueError(f"{source_id} path does not match frozen manifest")
        observed = sha256(PROJECT_ROOT / relative_path)
        if observed.lower() != manifest.loc[source_id, "sha256"].lower():
            raise ValueError(f"input drift detected for {source_id}")
        verified[source_id] = observed
    return verified


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8.2, fontweight="bold", ha="left", va="bottom", color="#272727")


def parse_gene_file(value: str) -> tuple[str, str]:
    stem = value.replace("differential_", "").replace(".tsv", "")
    for suffix, display in [("_DMED_vs_normal", "Diabetic ED vs reference"), ("_nonDM_ED_vs_normal", "Non-diabetic ED vs reference")]:
        if stem.endswith(suffix):
            cell_type = stem[: -len(suffix)]
            return cell_type, display
    raise ValueError(f"unrecognized gene comparison file: {value}")


def short_cell_type(value: str) -> str:
    return {
        "ACKR1_endothelial": "ACKR1 EC", "vascular_endothelial": "Vascular EC", "macrophage_myeloid": "Macrophage",
        "smooth_muscle": "Smooth muscle", "fibroblast": "Fibroblast", "pericyte": "Pericyte", "T_cell": "T cell",
    }.get(value, value.replace("_", " "))


def select_inputs() -> dict[str, pd.DataFrame]:
    summary = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S040"], sep="\t")
    count_map = dict(zip(summary["revised_candidate_class"], summary["candidate_genes"]))
    if count_map != {"ED_FDR_SUPPORTED": 3, "NOMINAL_P_EXPLORATORY": 563}:
        raise ValueError("candidate hierarchy does not match frozen 3/563 counts")
    count_display = pd.DataFrame([
        {"candidate_class": "ED-FDR supported", "count": 3, "display_order": 0},
        {"candidate_class": "Nominal-P exploratory", "count": 563, "display_order": 1},
    ])

    candidate = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S041"], sep="\t")
    target = candidate.loc[candidate["revised_candidate_class"].eq("ED_FDR_SUPPORTED")].copy()
    target = target.set_index("gene_symbol").reindex(GENE_ORDER).reset_index()
    if target["gene_symbol"].tolist() != GENE_ORDER or target.shape[0] != 3:
        raise ValueError("ED-FDR candidate hierarchy does not retain the three prespecified genes")
    if target["minimum_ED_FDR"].ge(0.05).any() or target["endothelial_meta_FDR"].ge(0.05).any():
        raise ValueError("target candidate FDR state differs from frozen hierarchy")

    celltype_gene = pd.read_csv(PROJECT_ROOT / "analysis/results/integration/IH_ED_celltype_gene_evidence.tsv", sep="\t", low_memory=False)
    fibro = celltype_gene.loc[(celltype_gene["gene_symbol"].isin(GENE_ORDER)) & (celltype_gene["cell_type"].eq("fibroblast")),
                              ["gene_symbol", "ED_logFC", "ED_FDR"]].set_index("gene_symbol").reindex(GENE_ORDER).reset_index()
    if fibro.shape[0] != 3 or fibro["ED_FDR"].ge(0.05).any():
        raise ValueError("fibroblast ED evidence does not retain all three FDR-supported candidates")

    doublet_class = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S044"], sep="\t").set_index("gene_symbol").reindex(GENE_ORDER).reset_index()
    if not (doublet_class["historical_best_tier"].eq("A") & doublet_class["sensitivity_best_tier"].eq("A")).all():
        raise ValueError("doublet sensitivity class does not retain A-to-A status")
    evidence_matrix = target.loc[:, ["gene_symbol", "endothelial_meta_z", "endothelial_meta_FDR"]].merge(fibro, on="gene_symbol").merge(
        doublet_class.loc[:, ["gene_symbol", "historical_best_tier", "sensitivity_best_tier"]], on="gene_symbol"
    )

    gene_sensitivity = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S042"], sep="\t")
    parsed = gene_sensitivity["file"].map(parse_gene_file)
    gene_sensitivity[["cell_type", "contrast"]] = pd.DataFrame(parsed.tolist(), index=gene_sensitivity.index)
    gene_sensitivity["cell_type_label"] = gene_sensitivity["cell_type"].map(short_cell_type)
    gene_sensitivity["old_FDR_lt_0_05"] = gene_sensitivity["old_FDR"] < 0.05
    gene_sensitivity["new_FDR_lt_0_05"] = gene_sensitivity["new_FDR"] < 0.05
    gene_sensitivity["complete_pair"] = gene_sensitivity[["old_logFC", "new_logFC", "old_FDR", "new_FDR"]].notna().all(axis=1)
    if gene_sensitivity.shape[0] != 31 or int(gene_sensitivity["complete_pair"].sum()) != 28:
        raise ValueError("doublet gene table does not retain 31 rows and 28 complete pairs")
    if int(gene_sensitivity["old_FDR_lt_0_05"].sum()) != 3 or int(gene_sensitivity["new_FDR_lt_0_05"].sum()) != 3:
        raise ValueError("doublet gene FDR states differ from frozen tables")

    pathway = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S043"], sep="\t")
    pathway["complete_pair"] = pathway[["old_NES", "new_NES", "old_padj", "new_padj"]].notna().all(axis=1)
    complete = pathway.loc[pathway["complete_pair"]].copy()
    old_sig = complete["old_padj"] < 0.05
    new_sig = complete["new_padj"] < 0.05
    pathway.loc[pathway["complete_pair"], "significance_state"] = np.select(
        [old_sig & new_sig, ~old_sig & new_sig, old_sig & ~new_sig],
        ["Retained significant", "Gained significance", "Lost significance"],
        default="Retained non-significant",
    )
    pathway.loc[~pathway["complete_pair"], "significance_state"] = "Incomplete estimate"
    state_order = ["Retained significant", "Gained significance", "Lost significance", "Retained non-significant"]
    state_counts = pathway.loc[pathway["complete_pair"], "significance_state"].value_counts().reindex(state_order, fill_value=0).reset_index()
    state_counts.columns = ["significance_state", "count"]
    state_counts["display_order"] = range(len(state_order))
    if state_counts["count"].tolist() != [263, 23, 26, 368] or pathway.shape[0] != 683 or complete.shape[0] != 680:
        raise ValueError("doublet pathway contingency states differ from frozen sensitivity table")
    return {
        "count_display": count_display, "evidence_matrix": evidence_matrix, "gene_sensitivity": gene_sensitivity,
        "pathway": pathway, "state_counts": state_counts,
    }


def write_source_data(data: dict[str, pd.DataFrame]) -> None:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    files = {
        "count_display": "Figure4a_candidate_hierarchy_source_data.tsv",
        "evidence_matrix": "Figure4b_candidate_evidence_matrix_source_data.tsv",
        "gene_sensitivity": "Figure4c_doublet_gene_sensitivity_source_data.tsv",
        "pathway": "Figure4d_doublet_pathway_sensitivity_source_data.tsv",
        "state_counts": "Figure4e_pathway_state_count_source_data.tsv",
    }
    for key, filename in files.items():
        data[key].to_csv(SOURCE_OUT / filename, sep="\t", index=False)


def plot_candidate_matrix(ax: plt.Axes, data: pd.DataFrame) -> None:
    n = data.shape[0]
    for row, item in data.reset_index(drop=True).iterrows():
        endo_color = "#B94B4A" if item["endothelial_meta_z"] > 0 else "#3C78B5"
        ed_color = "#B94B4A" if item["ED_logFC"] > 0 else "#3C78B5"
        cells = [
            (0, endo_color, f"z = {item['endothelial_meta_z']:.2f}\nFDR = {item['endothelial_meta_FDR']:.3g}"),
            (1, ed_color, f"Fibroblast\nlog2FC = {item['ED_logFC']:.2f}\nFDR = {item['ED_FDR']:.3g}"),
            (2, "#8264A5", "Class A → A\nretained after\ndoublet exclusion"),
        ]
        for col, color, text in cells:
            ax.add_patch(Rectangle((col - 0.48, row - 0.46), 0.96, 0.92, facecolor=color, edgecolor="white", linewidth=0.8))
            ax.text(col, row, text, ha="center", va="center", fontsize=4.0, color="white", linespacing=1.05)
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Endothelial\nmeta", "ED cell-type\neffect", "Doublet\nsensitivity"], fontsize=4.6)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_yticks(range(n))
    ax.set_yticklabels(data["gene_symbol"], fontsize=5.1, fontweight="bold")
    ax.tick_params(axis="y", length=0, pad=1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Frozen evidence matrix", fontsize=6.1, fontweight="bold", loc="left", pad=2)


def make_figure(data: dict[str, pd.DataFrame]) -> None:
    count_display = data["count_display"]
    evidence = data["evidence_matrix"]
    gene = data["gene_sensitivity"]
    pathway = data["pathway"]
    states = data["state_counts"]
    fig = plt.figure(figsize=(7.205, 5.512), constrained_layout=False)
    grid = fig.add_gridspec(2, 6, height_ratios=[0.93, 1.42], width_ratios=[1, 1, 1, 1, 1, 1])
    fig.subplots_adjust(left=0.112, right=0.958, top=0.950, bottom=0.105, wspace=0.90, hspace=0.78)
    ax_a = fig.add_subplot(grid[0, :2])
    ax_b = fig.add_subplot(grid[0, 2:4])
    ax_e = fig.add_subplot(grid[0, 4:])
    c_grid = grid[1, :3].subgridspec(1, 2, wspace=0.33)
    ax_c = [fig.add_subplot(c_grid[0, i]) for i in range(2)]
    ax_d = fig.add_subplot(grid[1, 3:])

    # a — the frozen two-level candidate hierarchy.
    colors = ["#704B87", "#A6A6A6"]
    y = np.arange(2)
    ax_a.barh(y, count_display["count"], color=colors, height=0.54)
    for ypos, value in zip(y, count_display["count"]):
        ax_a.text(value * 1.09, ypos, f"{int(value):,}", va="center", fontsize=7.0, fontweight="bold")
    ax_a.set_xscale("log")
    ax_a.set_xlim(1, 1300)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(count_display["candidate_class"], fontsize=5.0)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Candidate genes (log10 scale)", fontsize=4.6, labelpad=1)
    ax_a.tick_params(axis="x", labelsize=4.3, length=2, pad=1)
    ax_a.tick_params(axis="y", length=0, pad=1)
    ax_a.set_title("Frozen ED support hierarchy", fontsize=6.2, fontweight="bold", loc="left", pad=2)
    ax_a.text(0.0, -0.31, "566 genes satisfy the upstream candidate screen; class area is not used as an encoding.",
              transform=ax_a.transAxes, fontsize=3.8, color="#5F5F5F", va="top")
    panel_label(ax_a, "a")

    # b — compact, non-composite candidate evidence matrix.
    plot_candidate_matrix(ax_b, evidence)
    panel_label(ax_b, "b", x=-0.20)

    # e — deterministic state contingency; no test is introduced.
    state_short = ["Retained\nsignificant", "Gained", "Lost", "Retained\nnon-significant"]
    positions = np.arange(states.shape[0])
    bars = ax_e.bar(positions, states["count"], color=[STATE_COLORS[x] for x in states["significance_state"]], width=0.70)
    ax_e.set_ylim(0, 430)
    for bar, value in zip(bars, states["count"]):
        ax_e.text(bar.get_x() + bar.get_width() / 2, value + 11, f"{int(value)}", ha="center", va="bottom", fontsize=5.4, fontweight="bold")
    ax_e.set_xticks(positions)
    ax_e.set_xticklabels(state_short, fontsize=4.0)
    ax_e.set_ylabel("Hallmark states", fontsize=4.8, labelpad=1)
    ax_e.tick_params(axis="y", labelsize=4.2, length=2, pad=1)
    ax_e.set_title("FDR-state contingency", fontsize=6.0, fontweight="bold", loc="left", pad=2)
    ax_e.text(0.0, -0.31, "680 complete old/new pairs; 3 incomplete estimates are not classified.", transform=ax_e.transAxes,
              fontsize=3.7, color="#5F5F5F", va="top")
    panel_label(ax_e, "e", x=-0.18)

    # c — paired old/new effects, faceted by contrast; FDR ring is copied from sealed analyses.
    contrasts = ["Diabetic ED vs reference", "Non-diabetic ED vs reference"]
    cell_order = ["ACKR1_endothelial", "fibroblast", "macrophage_myeloid", "pericyte", "smooth_muscle", "T_cell", "vascular_endothelial"]
    for index, (ax, contrast) in enumerate(zip(ax_c, contrasts)):
        d = gene.loc[(gene["contrast"].eq(contrast)) & gene["complete_pair"]].copy()
        for gi, gene_symbol in enumerate(GENE_ORDER):
            g = d.loc[d["gene_symbol"].eq(gene_symbol)]
            for row in g.itertuples(index=False):
                xpos = cell_order.index(row.cell_type) + (gi - 1) * 0.18
                ax.plot([xpos - 0.045, xpos + 0.045], [row.old_logFC, row.new_logFC], color=GENE_COLORS[gene_symbol], linewidth=0.65, alpha=0.75, zorder=1)
                old_edge = "#202020" if row.old_FDR_lt_0_05 else "white"
                new_edge = "#202020" if row.new_FDR_lt_0_05 else GENE_COLORS[gene_symbol]
                ax.scatter(xpos - 0.045, row.old_logFC, s=16, color=GENE_COLORS[gene_symbol], edgecolor=old_edge, linewidth=0.7, zorder=2)
                ax.scatter(xpos + 0.045, row.new_logFC, s=16, facecolor="white", edgecolor=new_edge, linewidth=0.9, zorder=2)
        ax.axhline(0, color="#7C7C7C", linewidth=0.55, zorder=0)
        ax.set_xlim(-0.55, len(cell_order) - 0.45)
        ax.set_xticks(range(len(cell_order)))
        ax.set_xticklabels([short_cell_type(x) for x in cell_order], fontsize=3.55, rotation=38, ha="right")
        ax.tick_params(axis="x", length=0, pad=1)
        ax.tick_params(axis="y", labelsize=4.2, length=2, pad=1)
        ax.set_title(contrast, fontsize=4.9, fontweight="bold", pad=2)
        if index == 0:
            ax.set_ylabel("ED log2FC", fontsize=4.7, labelpad=1)
        else:
            ax.set_yticklabels([])
    ax_c[0].text(0.0, 1.17, "Doublet-exclusion sensitivity: sealed gene effects", transform=ax_c[0].transAxes,
                 fontsize=5.7, fontweight="bold", va="bottom")
    fig.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GENE_COLORS["TYMS"], markeredgecolor="white", markersize=4, label="TYMS"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GENE_COLORS["EFNB2"], markeredgecolor="white", markersize=4, label="EFNB2"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GENE_COLORS["LRRC17"], markeredgecolor="white", markersize=4, label="LRRC17"),
        Line2D([0], [0], marker="o", color="#666666", markerfacecolor="#666666", markersize=3.5, label="before"),
        Line2D([0], [0], marker="o", color="#666666", markerfacecolor="white", markersize=3.5, label="after"),
        Line2D([0], [0], marker="o", color="#202020", markerfacecolor="white", markersize=4, label="black ring: FDR < 0.05"),
    ], loc="lower left", bbox_to_anchor=(0.108, 0.006), ncol=3, fontsize=3.3, columnspacing=0.6, handletextpad=0.25, borderpad=0.1)
    panel_label(ax_c[0], "c", x=-0.26, y=1.18)

    # d — complete pathway NES pairs only, with FDR state transitions copied from S043.
    complete_pathway = pathway.loc[pathway["complete_pair"]].copy()
    limit = float(np.ceil(max(complete_pathway["old_NES"].abs().max(), complete_pathway["new_NES"].abs().max())))
    for state in ["Retained non-significant", "Retained significant", "Gained significance", "Lost significance"]:
        d = complete_pathway.loc[complete_pathway["significance_state"].eq(state)]
        ax_d.scatter(d["old_NES"], d["new_NES"], s=5.0 if state != "Retained non-significant" else 3.4,
                     color=STATE_COLORS[state], alpha=0.74 if state != "Retained non-significant" else 0.38,
                     linewidths=0, rasterized=True, label=state)
    ax_d.plot([-limit, limit], [-limit, limit], color="#555555", linewidth=0.65, linestyle="--", zorder=0)
    ax_d.axhline(0, color="#A0A0A0", linewidth=0.45, zorder=0)
    ax_d.axvline(0, color="#A0A0A0", linewidth=0.45, zorder=0)
    ax_d.set_xlim(-limit, limit)
    ax_d.set_ylim(-limit, limit)
    ax_d.set_xlabel("Before exclusion NES", fontsize=4.9, labelpad=1)
    ax_d.set_ylabel("After exclusion NES", fontsize=4.9, labelpad=1)
    ax_d.tick_params(labelsize=4.3, length=2, pad=1)
    ax_d.set_title("Hallmark pathway sensitivity", fontsize=6.0, fontweight="bold", loc="left", pad=2)
    ax_d.text(0.03, 0.96, "n = 680 complete pairs\nidentity line", transform=ax_d.transAxes, ha="left", va="top", fontsize=4.2,
              bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.0})
    ax_d.legend(loc="lower right", fontsize=3.5, markerscale=1.5, labelspacing=0.25, borderpad=0.2, handletextpad=0.3)
    panel_label(ax_d, "d")

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "Figure4_candidate_hierarchy_and_doublet_sensitivity"
    fig.savefig(base.with_suffix(".svg"), facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white")
    fig.savefig(base.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main() -> None:
    verified = verify_inputs()
    data = select_inputs()
    write_source_data(data)
    make_figure(data)
    audit = {
        "status": "PASS", "figure": "Figure 4", "backend": "Python/matplotlib",
        "input_manifest_validation": "PASS", "input_sha256": verified,
        "assertions": {
            "candidate_counts": {"ED_FDR_SUPPORTED": 3, "NOMINAL_P_EXPLORATORY": 563},
            "candidate_genes": GENE_ORDER,
            "doublet_class_retention": {gene: "A_to_A" for gene in GENE_ORDER},
            "doublet_gene_rows": 31, "doublet_gene_complete_pairs": 28,
            "pathway_rows": 683, "pathway_complete_pairs": 680,
            "pathway_state_counts": dict(zip(data["state_counts"]["significance_state"], data["state_counts"]["count"])),
        },
        "output_dimensions_mm": [183, 140],
        "interpretation_boundary": "Doublet exclusion is a within-atlas sensitivity analysis; nested cell-type readouts are not independent confirmations, and retained class does not establish replication or causality.",
    }
    (OUT / "Figure4_build_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
