from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = OUTPUT_ROOT.parent
OUT = OUTPUT_ROOT / "figures" / "Figure2_20260911"
SOURCE_OUT = OUT / "source_data"
INPUT_MANIFEST = OUTPUT_ROOT / "figure_contracts" / "input_manifest.sha256.tsv"

SOURCE_IDS = {
    "S010": "revision_v2/02_results/endothelial_three_study_summary_v2.tsv",
    "S011": "revision_v2/02_results/endothelial_three_study_evidence_v2.tsv",
    "S020": "analysis/results/GSE243023/differential_IH_CTRL_vs_normoxia_CTRL.tsv",
    "S021": "analysis/results/GSE205050/differential_IH_vehicle_vs_normoxia_vehicle.tsv",
    "S022": "analysis/results/GSE10723/differential_IH_vs_control.tsv",
    "S024": "analysis/results/integration/bulk_effect_spearman_correlations.tsv",
    "S025": "analysis/results/integration/Hallmark_cross_dataset_evidence.tsv",
    "S026": "analysis/results/phase3/reactome_multilevel/Reactome_multilevel_summary.tsv",
    "S027": "analysis/results/phase3/reactome_multilevel/Endothelial_IH_meta_Reactome_fgseaMultilevel.tsv",
}

COLORS = {
    "ink": "#272727",
    "muted": "#767676",
    "grey": "#CFCFCF",
    "grey_dark": "#656565",
    "blue": "#3775BA",
    "red": "#B64342",
    "purple": "#9A4D8E",
    "soft_blue": "#D9EAF7",
    "soft_purple": "#E9DDF1",
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
        if source_id not in manifest.index:
            raise ValueError(f"{source_id} missing from input manifest")
        if manifest.loc[source_id, "relative_path"] != relative_path:
            raise ValueError(f"{source_id} path does not match the frozen contract")
        observed = sha256(PROJECT_ROOT / relative_path)
        if observed.lower() != manifest.loc[source_id, "sha256"].lower():
            raise ValueError(f"input drift detected for {source_id}")
        verified[source_id] = observed
    return verified


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.07, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8.2, fontweight="bold", color=COLORS["ink"],
            ha="left", va="bottom")


def short_label(pathway: str) -> str:
    labels = {
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNFα signaling via NF-κB",
        "HALLMARK_UV_RESPONSE_DN": "UV response DN",
        "HALLMARK_HYPOXIA": "Hypoxia",
        "HALLMARK_DNA_REPAIR": "DNA repair",
        "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint",
        "HALLMARK_E2F_TARGETS": "E2F targets",
        "REACTOME_CELL_CYCLE_MITOTIC": "Cell cycle, mitotic",
        "REACTOME_MITOTIC_METAPHASE_AND_ANAPHASE": "Mitotic metaphase and anaphase",
        "REACTOME_M_PHASE": "M phase",
        "REACTOME_CELL_CYCLE_CHECKPOINTS": "Cell cycle checkpoints",
        "REACTOME_SYNTHESIS_OF_DNA": "Synthesis of DNA",
        "REACTOME_MITOTIC_PROMETAPHASE": "Mitotic prometaphase",
    }
    return labels.get(pathway, pathway.replace("HALLMARK_", "").replace("REACTOME_", "").replace("_", " ").title())


def select_inputs() -> dict[str, pd.DataFrame | float | int]:
    endothelial_summary = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S010"], sep="\t")
    evidence = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S011"], sep="\t", low_memory=False)
    volcano_specs = [
        ("GSE243023", SOURCE_IDS["S020"]),
        ("GSE205050", SOURCE_IDS["S021"]),
        ("GSE10723", SOURCE_IDS["S022"]),
    ]
    volcano_parts = []
    for dataset, relative_path in volcano_specs:
        frame = pd.read_csv(PROJECT_ROOT / relative_path, sep="\t")
        required = {"gene_symbol", "logFC", "P.Value", "adj.P.Val"}
        if not required.issubset(frame.columns):
            raise ValueError(f"volcano source lacks required columns: {dataset}")
        frame = frame.loc[:, ["gene_symbol", "logFC", "P.Value", "adj.P.Val"]].copy()
        frame["dataset"] = dataset
        frame["FDR_lt_0_05"] = frame["adj.P.Val"] < 0.05
        volcano_parts.append(frame)
    volcano = pd.concat(volcano_parts, ignore_index=True)

    b_data = evidence.loc[:, [
        "gene_symbol", "strict_lodo_FDR_candidate", "meta_FDR", "GSE243023_logFC", "GSE205050_logFC", "GSE10723_logFC",
    ]].copy()
    b_data["strict_lodo_FDR_candidate"] = b_data["strict_lodo_FDR_candidate"].astype(bool)
    b_data = b_data.dropna(subset=["GSE243023_logFC", "GSE205050_logFC", "GSE10723_logFC"])
    if b_data.shape[0] != 5139:
        raise ValueError("shared three-study effect table must contain 5,139 genes")
    correlation = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S024"], sep="\t")
    rho_row = correlation.loc[
        (correlation["dataset_x"].eq("endo243_logFC")) & (correlation["dataset_y"].eq("endo205_logFC")),
        "spearman_rho",
    ]
    if rho_row.shape[0] != 1:
        raise ValueError("frozen endothelial effect correlation is unavailable")
    rho = float(rho_row.iloc[0])

    summary_values = dict(zip(endothelial_summary["metric"], endothelial_summary["value"]))
    count_data = pd.DataFrame([
        {"filter": "All three\nsame direction", "count": summary_values["three_study_all_same_direction"], "order": 1},
        {"filter": "Directional\nstability screen", "count": summary_values["directional_stability_screen"], "order": 2},
        {"filter": "Full meta\nFDR < 0.05", "count": summary_values["full_meta_FDR_lt_0_05"], "order": 3},
        {"filter": "All-three LODO\nFDR < 0.05", "count": summary_values["strict_all_three_lodo_FDR_lt_0_05"], "order": 4},
    ])
    if count_data["count"].tolist() != [1766, 1035, 1282, 258]:
        raise ValueError("frozen filter counts differ from Figure 2 contract")

    strict = evidence.loc[evidence["strict_lodo_FDR_candidate"].astype(bool)].copy()
    positive = strict.loc[strict["meta_z"] > 0].sort_values(["meta_FDR", "gene_symbol"], ascending=[True, True]).head(6)
    negative = strict.loc[strict["meta_z"] < 0].sort_values(["meta_FDR", "gene_symbol"], ascending=[True, True]).head(6)
    if positive.shape[0] != 6 or negative.shape[0] != 6:
        raise ValueError("strict LODO subset does not provide six positive and six negative genes")
    strict_heat = pd.concat([positive.assign(meta_direction="positive"), negative.assign(meta_direction="negative")], ignore_index=True)
    strict_heat = strict_heat.loc[:, [
        "gene_symbol", "meta_direction", "meta_FDR", "GSE243023_logFC", "GSE205050_logFC", "GSE10723_logFC",
    ]]

    hallmark = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S025"], sep="\t")
    hallmark["endothelial_max_FDR"] = pd.to_numeric(hallmark["endothelial_max_FDR"], errors="coerce")
    hypoxia = hallmark.loc[hallmark["pathway"].eq("HALLMARK_HYPOXIA")].copy()
    if hypoxia.shape[0] != 1:
        raise ValueError("HALLMARK_HYPOXIA must occur exactly once")
    selected_hallmark = hallmark.loc[~hallmark["pathway"].eq("HALLMARK_HYPOXIA")].sort_values(
        ["endothelial_max_FDR", "pathway"], ascending=[True, True]
    ).head(5)
    hallmark_plot = pd.concat([hypoxia, selected_hallmark], ignore_index=True).loc[:, [
        "pathway", "endothelial_NES_mean", "endothelial_max_FDR",
    ]].copy()
    hallmark_plot = hallmark_plot.rename(columns={"endothelial_NES_mean": "NES", "endothelial_max_FDR": "FDR"})
    hallmark_plot["collection"] = "Hallmark"
    hallmark_plot["selection_rule"] = np.where(
        hallmark_plot["pathway"].eq("HALLMARK_HYPOXIA"), "prespecified strict bridge", "lowest endothelial_max_FDR excluding Hypoxia"
    )

    reactome = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S027"], sep="\t")
    reactome["padj"] = pd.to_numeric(reactome["padj"], errors="coerce")
    reactome["NES"] = pd.to_numeric(reactome["NES"], errors="coerce")
    reactome_plot = reactome.loc[reactome["padj"] < 0.05].copy()
    reactome_plot["abs_NES"] = reactome_plot["NES"].abs()
    reactome_plot = reactome_plot.sort_values(["abs_NES", "pathway"], ascending=[False, True]).head(6)
    if reactome_plot.shape[0] != 6:
        raise ValueError("Reactome selection must contain six FDR-significant pathways")
    reactome_plot = reactome_plot.loc[:, ["pathway", "NES", "padj"]].rename(columns={"padj": "FDR"})
    reactome_plot["collection"] = "Reactome"
    reactome_plot["selection_rule"] = "FDR < 0.05; largest absolute NES; alphabetical tie-break"

    pathway_plot = pd.concat([hallmark_plot, reactome_plot], ignore_index=True)
    pathway_plot["minus_log10_FDR"] = -np.log10(pathway_plot["FDR"].clip(lower=np.finfo(float).tiny))
    return {
        "volcano": volcano,
        "b_data": b_data,
        "rho": rho,
        "count_data": count_data,
        "strict_heat": strict_heat,
        "pathway_plot": pathway_plot,
    }


def write_source_data(data: dict[str, pd.DataFrame | float | int]) -> None:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    cast = lambda key: data[key]  # noqa: E731
    assert isinstance(cast("volcano"), pd.DataFrame)
    cast("volcano").to_csv(SOURCE_OUT / "Figure2a_volcano_source_data.tsv", sep="\t", index=False)
    cast("b_data").to_csv(SOURCE_OUT / "Figure2b_shared_effect_source_data.tsv", sep="\t", index=False)
    pd.DataFrame([{"reported_spearman_rho_endo243_vs_endo205": float(cast("rho")), "n_shared_genes": 5139}]).to_csv(
        SOURCE_OUT / "Figure2b_reported_correlation.tsv", sep="\t", index=False
    )
    cast("count_data").to_csv(SOURCE_OUT / "Figure2c_filter_count_source_data.tsv", sep="\t", index=False)
    cast("strict_heat").to_csv(SOURCE_OUT / "Figure2d_strict_lodo_heatmap_source_data.tsv", sep="\t", index=False)
    cast("pathway_plot").to_csv(SOURCE_OUT / "Figure2e_pathway_dotplot_source_data.tsv", sep="\t", index=False)


def plot_volcano(ax: plt.Axes, data: pd.DataFrame, dataset: str, xlim: float, ylim: float) -> None:
    d = data.loc[data["dataset"].eq(dataset)]
    y = -np.log10(d["P.Value"].clip(lower=np.finfo(float).tiny))
    sig = d["FDR_lt_0_05"].to_numpy(dtype=bool)
    ax.scatter(d.loc[~sig, "logFC"], y.loc[~sig], s=2.2, color=COLORS["grey"], alpha=0.45, linewidths=0, rasterized=True)
    pos = sig & (d["logFC"].to_numpy() > 0)
    neg = sig & (d["logFC"].to_numpy() < 0)
    ax.scatter(d.loc[pos, "logFC"], y.loc[pos], s=3.2, color=COLORS["red"], alpha=0.72, linewidths=0, rasterized=True)
    ax.scatter(d.loc[neg, "logFC"], y.loc[neg], s=3.2, color=COLORS["blue"], alpha=0.72, linewidths=0, rasterized=True)
    ax.axvline(0, color=COLORS["muted"], linewidth=0.55)
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(0, ylim)
    ax.set_title(dataset, fontsize=6.3, fontweight="bold", pad=2)
    ax.text(0.98, 0.96, f"FDR < 0.05: {int(sig.sum()):,}", transform=ax.transAxes, ha="right", va="top", fontsize=4.7)
    ax.tick_params(labelsize=4.8, length=2, pad=1)
    ax.set_xlabel("log2FC", fontsize=5.0, labelpad=1)


def dot_panel(ax: plt.Axes, data: pd.DataFrame, collection: str, norm: TwoSlopeNorm) -> plt.Collection:
    d = data.loc[data["collection"].eq(collection)].copy().sort_values(["NES", "pathway"], ascending=[True, True])
    y = np.arange(d.shape[0])
    size = 18 + 16 * np.clip(d["minus_log10_FDR"].to_numpy(dtype=float), 0, 10)
    scatter = ax.scatter(d["NES"], y, s=size, c=d["NES"], cmap="RdBu_r", norm=norm, edgecolor="white", linewidth=0.45)
    ax.axvline(0, color=COLORS["muted"], linewidth=0.6, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([short_label(p) for p in d["pathway"]], fontsize=4.9)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=4.8, length=2, pad=1)
    ax.set_xlabel("Endothelial meta NES", fontsize=5.0, labelpad=1)
    ax.set_title(collection, fontsize=6.2, fontweight="bold", loc="left", pad=2)
    return scatter


def make_figure(data: dict[str, pd.DataFrame | float | int]) -> None:
    volcano = data["volcano"]
    b_data = data["b_data"]
    count_data = data["count_data"]
    strict_heat = data["strict_heat"]
    pathway_plot = data["pathway_plot"]
    assert isinstance(volcano, pd.DataFrame) and isinstance(b_data, pd.DataFrame)
    assert isinstance(count_data, pd.DataFrame) and isinstance(strict_heat, pd.DataFrame) and isinstance(pathway_plot, pd.DataFrame)
    rho = float(data["rho"])

    fig = plt.figure(figsize=(7.205, 6.102), constrained_layout=False)
    grid = fig.add_gridspec(3, 6, height_ratios=[1.03, 1.24, 1.17], width_ratios=[1, 1, 1, 1, 1, 1])
    fig.subplots_adjust(left=0.116, right=0.940, top=0.955, bottom=0.075, wspace=0.85, hspace=0.80)
    a_grid = grid[0, :4].subgridspec(1, 3, wspace=0.38)
    ax_a = [fig.add_subplot(a_grid[0, i]) for i in range(3)]
    ax_c = fig.add_subplot(grid[0, 4:])
    ax_b = fig.add_subplot(grid[1, :4])
    ax_d = fig.add_subplot(grid[1, 4:])
    e_grid = grid[2, :].subgridspec(1, 2, wspace=0.58)
    ax_eh = fig.add_subplot(e_grid[0, 0])
    ax_er = fig.add_subplot(e_grid[0, 1])

    # a — aligned volcanoes, retaining null FDR findings.
    xlim = float(np.ceil(volcano["logFC"].abs().max()))
    ylim = float(np.ceil((-np.log10(volcano["P.Value"].clip(lower=np.finfo(float).tiny))).max()))
    for i, dataset in enumerate(["GSE243023", "GSE205050", "GSE10723"]):
        plot_volcano(ax_a[i], volcano, dataset, xlim, ylim)
        if i == 0:
            ax_a[i].set_ylabel("−log10 P", fontsize=5.0, labelpad=1)
        else:
            ax_a[i].set_yticklabels([])
    add_panel_label(ax_a[0], "a", x=-0.24, y=1.10)

    # b — hero effect density. The rho is copied from its frozen summary table.
    x = b_data["GSE243023_logFC"].to_numpy(dtype=float)
    y = b_data["GSE205050_logFC"].to_numpy(dtype=float)
    limit = max(1.8, float(np.ceil(np.max(np.abs(np.concatenate([x, y]))))))
    hb = ax_b.hexbin(x, y, gridsize=58, mincnt=1, cmap="Greys", linewidths=0, bins="log", rasterized=True)
    strict = b_data["strict_lodo_FDR_candidate"].to_numpy(dtype=bool)
    ax_b.scatter(x[strict], y[strict], s=4.5, facecolors="none", edgecolors=COLORS["purple"], linewidths=0.45, alpha=0.85, label="258 strict LODO-FDR genes")
    ax_b.axhline(0, color=COLORS["muted"], linewidth=0.6)
    ax_b.axvline(0, color=COLORS["muted"], linewidth=0.6)
    ax_b.set_xlim(-limit, limit)
    ax_b.set_ylim(-limit, limit)
    ax_b.set_xlabel("GSE243023 log2FC", fontsize=5.4, labelpad=1)
    ax_b.set_ylabel("GSE205050 log2FC", fontsize=5.4, labelpad=1)
    ax_b.tick_params(labelsize=5.0, length=2, pad=1)
    ax_b.text(0.03, 0.96, f"reported Spearman ρ = {rho:.3f}\nn = 5,139 shared genes", transform=ax_b.transAxes,
              ha="left", va="top", fontsize=5.2, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 1.2})
    ax_b.legend(loc="lower right", fontsize=4.8, handletextpad=0.3, borderpad=0.2)
    add_panel_label(ax_b, "b")

    # c — non-equivalent frozen filter counts, deliberately not drawn as a set intersection.
    count_colors = [COLORS["grey_dark"], "#88A9C9", COLORS["blue"], COLORS["purple"]]
    positions = np.arange(count_data.shape[0])
    bars = ax_c.bar(positions, count_data["count"], color=count_colors, width=0.72, edgecolor="white", linewidth=0.45)
    ymax = float(count_data["count"].max()) * 1.25
    ax_c.set_ylim(0, ymax)
    for bar, value in zip(bars, count_data["count"]):
        ax_c.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + ymax * 0.025, f"{int(value):,}", ha="center", va="bottom", fontsize=5.5, fontweight="bold")
    ax_c.set_xticks(positions)
    ax_c.set_xticklabels(count_data["filter"], fontsize=4.8)
    ax_c.set_ylabel("Genes", fontsize=5.2, labelpad=1)
    ax_c.tick_params(axis="y", labelsize=4.7, length=2, pad=1)
    ax_c.text(0.50, -0.36, "Frozen filters are not interchangeable set intersections.", transform=ax_c.transAxes,
              ha="center", va="top", fontsize=4.5, style="italic", color=COLORS["muted"])
    add_panel_label(ax_c, "c", x=-0.18, y=1.10)

    # d — deterministic 6-up / 6-down strict LODO heatmap.
    heat_cols = ["GSE243023_logFC", "GSE205050_logFC", "GSE10723_logFC"]
    matrix = strict_heat[heat_cols].to_numpy(dtype=float)
    heat_limit = max(1.0, float(np.ceil(np.abs(matrix).max() * 2) / 2))
    im = ax_d.imshow(matrix, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-heat_limit, vcenter=0, vmax=heat_limit))
    ax_d.set_xticks(range(3))
    ax_d.set_xticklabels(["IH 1", "IH 2", "IH 3"], fontsize=5.0)
    ax_d.xaxis.tick_top()
    ax_d.tick_params(axis="x", length=0, pad=2)
    ax_d.set_yticks(range(strict_heat.shape[0]))
    ax_d.set_yticklabels(strict_heat["gene_symbol"], fontsize=4.7)
    ax_d.tick_params(axis="y", length=0, pad=1)
    for i, direction in enumerate(strict_heat["meta_direction"]):
        if direction == "positive":
            ax_d.add_patch(plt.Rectangle((-0.56, i - 0.5), 0.05, 1, color=COLORS["red"], clip_on=False))
        else:
            ax_d.add_patch(plt.Rectangle((-0.56, i - 0.5), 0.05, 1, color=COLORS["blue"], clip_on=False))
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax_d, fraction=0.050, pad=0.06)
    cb.set_label("log2FC", fontsize=5.0)
    cb.ax.tick_params(labelsize=4.7, length=2)
    ax_d.text(0.0, -0.15, "Top 6 positive + top 6 negative strict genes\nby meta FDR", transform=ax_d.transAxes,
              fontsize=4.4, color=COLORS["muted"], va="top")
    add_panel_label(ax_d, "d", x=-0.28, y=1.10)

    # e — separate Hallmark and Reactome dot panels, retaining their distinct FDR families.
    nes_limit = max(2.0, float(np.ceil(pathway_plot["NES"].abs().max())))
    norm = TwoSlopeNorm(vmin=-nes_limit, vcenter=0, vmax=nes_limit)
    sc_h = dot_panel(ax_eh, pathway_plot, "Hallmark", norm)
    dot_panel(ax_er, pathway_plot, "Reactome", norm)
    cb2 = fig.colorbar(sc_h, ax=[ax_eh, ax_er], fraction=0.026, pad=0.02)
    cb2.set_label("NES", fontsize=5.0)
    cb2.ax.tick_params(labelsize=4.7, length=2)
    ax_eh.text(0.0, -0.28, "Dot area = −log10 FDR within the original collection-specific family.",
               transform=ax_eh.transAxes, fontsize=4.5, color=COLORS["muted"], va="top")
    add_panel_label(ax_eh, "e", x=-0.13, y=1.10)

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "Figure2_endothelial_convergence_and_heterogeneity"
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
        "status": "PASS",
        "figure": "Figure 2",
        "backend": "Python/matplotlib",
        "input_manifest_validation": "PASS",
        "input_sha256": verified,
        "assertions": {
            "shared_three_study_genes": 5139,
            "reported_endo243_vs_endo205_spearman_rho": float(data["rho"]),
            "frozen_filter_counts": data["count_data"]["count"].tolist(),
            "strict_lodo_heatmap_genes": 12,
            "strict_lodo_heatmap_positive_negative": [6, 6],
            "pathway_dotplot_counts": {"Hallmark": 6, "Reactome": 6},
        },
        "output_dimensions_mm": [183, 155],
        "interpretation_boundary": "Cross-study convergence is concentrated in a strict subset and pathway-level results; it does not indicate universal gene-level concordance.",
    }
    (OUT / "Figure2_build_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
