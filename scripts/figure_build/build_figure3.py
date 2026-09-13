from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
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
OUT = OUTPUT_ROOT / "figures" / "Figure3_20260912"
SOURCE_OUT = OUT / "source_data"
INPUT_MANIFEST = OUTPUT_ROOT / "figure_contracts" / "input_manifest.sha256.tsv"

SOURCE_IDS = {
    "S030": "analysis/figures/phase2_manuscript/source_data/Figure2a_UMAP_coordinates.tsv.gz",
    "S031": "analysis/results/GSE206528/celltype_counts_by_donor.tsv",
    "S033": "analysis/results/integration/IH_ED_Hallmark_celltype_evidence.tsv",
    "S034": "analysis/results/integration/IH_ED_celltype_gene_evidence.tsv",
    "S035": "analysis/results/GSE206528/pseudobulk_counts.csv",
    "S036": "analysis/results/GSE206528/pseudobulk_samples.csv",
}

CELLTYPE_ORDER = [
    "ACKR1_endothelial", "vascular_endothelial", "pericyte", "fibroblast", "smooth_muscle",
    "macrophage_myeloid", "T_cell",
]
ALL_CELLTYPE_ORDER = CELLTYPE_ORDER[:1] + ["Schwann_glial", "T_fibroblast_mixed"] + CELLTYPE_ORDER[1:]
CELLTYPE_LABEL = {
    "ACKR1_endothelial": "ACKR1 EC", "vascular_endothelial": "Vascular EC", "pericyte": "Pericyte",
    "fibroblast": "Fibroblast", "smooth_muscle": "Smooth muscle", "macrophage_myeloid": "Macrophage",
    "T_cell": "T cell", "T_fibroblast_mixed": "T-fibroblast", "Schwann_glial": "Schwann/glial",
}
CELL_COLORS = {
    "ACKR1_endothelial": "#E69F00", "vascular_endothelial": "#0072B2", "pericyte": "#56B4E9",
    "fibroblast": "#009E73", "smooth_muscle": "#CC79A7", "macrophage_myeloid": "#D55E00",
    "T_cell": "#F0E442", "T_fibroblast_mixed": "#8A7CC8", "Schwann_glial": "#8C8C8C",
}
GROUP_ORDER = ["normal", "organic_ED_nonDM", "organic_ED_DM"]
GROUP_LABEL = {
    "normal": "Reference\ntumour margin", "organic_ED_nonDM": "Non-diabetic ED\nimplantation", "organic_ED_DM": "Diabetic ED\nimplantation",
}
GROUP_COLOR = {"normal": "#7A7A7A", "organic_ED_nonDM": "#0072B2", "organic_ED_DM": "#CC79A7"}
GENE_ORDER = ["TYMS", "EFNB2", "LRRC17"]
PATHWAY_LABEL = {
    "HALLMARK_HYPOXIA": "Hypoxia", "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNFα signaling via NF-κB",
    "HALLMARK_DNA_REPAIR": "DNA repair", "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
    "HALLMARK_KRAS_SIGNALING_UP": "KRAS signaling up", "HALLMARK_MYC_TARGETS_V1": "MYC targets v1",
    "HALLMARK_MYC_TARGETS_V2": "MYC targets v2", "HALLMARK_UV_RESPONSE_DN": "UV response DN",
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
            raise ValueError(f"{source_id} path does not match the frozen contract")
        observed = sha256(PROJECT_ROOT / relative_path)
        if observed.lower() != manifest.loc[source_id, "sha256"].lower():
            raise ValueError(f"input drift detected for {source_id}")
        verified[source_id] = observed
    return verified


def panel_label(ax: plt.Axes, label: str, x: float = -0.10, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8.2, fontweight="bold", ha="left", va="bottom", color="#272727")


def compact_donor_label(sample_accession: str, condition: str) -> str:
    number = str(sample_accession).replace("GSM62559", "")
    if condition == "normal":
        return f"R{int(number) - 6}"
    if condition == "organic_ED_nonDM":
        return f"ED{int(number) - 9}"
    return f"D-ED{int(number) - 12}"


def select_inputs() -> dict[str, pd.DataFrame | int]:
    umap = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S030"], sep="\t", compression="infer")
    if umap.shape[0] != 64993 or not {"UMAP1", "UMAP2", "cell_type"}.issubset(umap.columns):
        raise ValueError("S030 does not match the Figure 3 UMAP contract")
    umap = umap.assign(descriptive_only=True)

    composition = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S031"], sep="\t")
    required_composition = {"sample_accession", "donor", "condition", "cell_type", "cell_count", "sample_total", "cell_fraction"}
    if not required_composition.issubset(composition.columns):
        raise ValueError("S031 lacks donor-composition columns")
    donor_meta = composition.loc[:, ["sample_accession", "condition"]].drop_duplicates()
    if donor_meta.shape[0] != 8 or donor_meta["condition"].value_counts().to_dict() != {"normal": 3, "organic_ED_nonDM": 3, "organic_ED_DM": 2}:
        raise ValueError("donor composition does not retain the 3/3/2 atlas structure")
    donor_order = donor_meta.assign(group_order=donor_meta["condition"].map({g: i for i, g in enumerate(GROUP_ORDER)})).sort_values(["group_order", "sample_accession"])
    donor_order["donor_label"] = [compact_donor_label(s, c) for s, c in zip(donor_order["sample_accession"], donor_order["condition"])]
    composition = composition.merge(donor_order.loc[:, ["sample_accession", "donor_label"]], on="sample_accession", how="left")
    composition["cell_type"] = pd.Categorical(composition["cell_type"], categories=ALL_CELLTYPE_ORDER, ordered=True)
    composition = composition.sort_values(["sample_accession", "cell_type"]).reset_index(drop=True)
    if not np.allclose(composition.groupby("sample_accession")["cell_fraction"].sum().to_numpy(), 1.0):
        raise ValueError("cell fractions do not sum to one per donor")

    pathways = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S033"], sep="\t")
    required_pathway = {"pathway", "cell_type", "NES", "padj", "IH_ED_direction_concordant", "cross_context_significant"}
    if not required_pathway.issubset(pathways.columns):
        raise ValueError("S033 lacks pathway-evidence columns")
    pathways = pathways.loc[pathways["cell_type"].isin(CELLTYPE_ORDER)].copy()
    # Existing `cross_context_significant` is used only as a frozen display flag; no new test is performed.
    selected_pathways = pathways.loc[pathways["cross_context_significant"].astype(bool), "pathway"].drop_duplicates().tolist()
    if "HALLMARK_HYPOXIA" not in selected_pathways or len(selected_pathways) != 8:
        raise ValueError("expected eight pre-existing cross-context-significant pathway programmes")
    counts = pathways.loc[pathways["cross_context_significant"].astype(bool)].groupby("pathway")["cell_type"].nunique()
    pathway_order = ["HALLMARK_HYPOXIA"] + sorted(
        [p for p in selected_pathways if p != "HALLMARK_HYPOXIA"], key=lambda p: (-int(counts[p]), p)
    )
    pathway_display = pathways.loc[pathways["pathway"].isin(pathway_order)].copy()
    grid = pd.MultiIndex.from_product([pathway_order, CELLTYPE_ORDER], names=["pathway", "cell_type"]).to_frame(index=False)
    pathway_display = grid.merge(pathway_display, on=["pathway", "cell_type"], how="left", validate="one_to_one")
    if pathway_display["NES"].isna().any():
        raise ValueError("pathway display contains missing cell-type estimates")
    pathway_display["pathway_label"] = pathway_display["pathway"].map(PATHWAY_LABEL)
    pathway_display["cell_type_label"] = pathway_display["cell_type"].map(CELLTYPE_LABEL)
    pathway_display["display_rule"] = "All pathways with >=1 pre-existing cross_context_significant row; Hypoxia first, then descending number of flagged compartments and alphabetical ties."
    pathway_display["frozen_display_flag"] = pathway_display["cross_context_significant"].astype(bool)

    genes = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S034"], sep="\t", low_memory=False)
    required_gene = {"gene_symbol", "cell_type", "ED_logFC", "ED_FDR", "IH_ED_direction_concordant"}
    if not required_gene.issubset(genes.columns):
        raise ValueError("S034 lacks gene-evidence columns")
    gene_grid = pd.MultiIndex.from_product([GENE_ORDER, CELLTYPE_ORDER], names=["gene_symbol", "cell_type"]).to_frame(index=False)
    gene_display = gene_grid.merge(
        genes.loc[genes["gene_symbol"].isin(GENE_ORDER) & genes["cell_type"].isin(CELLTYPE_ORDER),
                  ["gene_symbol", "cell_type", "ED_logFC", "ED_FDR", "IH_ED_direction_concordant", "evidence_tier"]],
        on=["gene_symbol", "cell_type"], how="left", validate="one_to_one"
    )
    gene_display["estimate_available"] = gene_display["ED_logFC"].notna()
    gene_display["FDR_lt_0_05"] = gene_display["ED_FDR"].lt(0.05).fillna(False)
    gene_display["cell_type_label"] = gene_display["cell_type"].map(CELLTYPE_LABEL)
    if int(gene_display["FDR_lt_0_05"].sum()) != 3 or set(gene_display.loc[gene_display["FDR_lt_0_05"], "cell_type"]) != {"fibroblast"}:
        raise ValueError("three FDR-supported gene effects must occur in fibroblast")

    pseudobulk = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S035"])
    sample_map = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S036"])
    sample_map = sample_map.rename(columns={"Unnamed: 0": "column_id"})
    if sample_map["column_id"].duplicated().any() or not set(sample_map["column_id"]).issubset(pseudobulk.columns):
        raise ValueError("pseudobulk columns cannot be mapped unambiguously to donor metadata")
    fibro_map = sample_map.loc[sample_map["cell_type"].eq("fibroblast")].copy()
    if fibro_map.shape[0] != 8 or fibro_map["condition"].value_counts().to_dict() != {"normal": 3, "organic_ED_nonDM": 3, "organic_ED_DM": 2}:
        raise ValueError("panel e fibroblast donor mapping is incomplete")
    count_index = pseudobulk.set_index("gene_symbol")
    if not set(GENE_ORDER).issubset(count_index.index):
        raise ValueError("panel e target genes are absent from pseudobulk count table")
    library_sizes = pseudobulk.drop(columns="gene_symbol").sum(axis=0)
    e_rows = []
    for gene in GENE_ORDER:
        for row in fibro_map.itertuples(index=False):
            raw_count = float(count_index.loc[gene, row.column_id])
            library_size = float(library_sizes[row.column_id])
            log2cpm = float(np.log2(((raw_count + 0.5) / (library_size + 1.0)) * 1_000_000.0))
            e_rows.append({
                "gene_symbol": gene, "cell_type": "fibroblast", "column_id": row.column_id,
                "sample_accession": row.sample_accession, "donor": row.donor, "condition": row.condition,
                "donor_label": compact_donor_label(row.sample_accession, row.condition), "raw_count": raw_count,
                "library_size": library_size, "log2CPM_prior_0_5": log2cpm,
                "display_rule": "Descriptive fibroblast donor dots only; log2(((count + 0.5)/(library size + 1))*1e6); no hypothesis test.",
            })
    donor_dots = pd.DataFrame(e_rows)
    if donor_dots.shape[0] != 24 or donor_dots.groupby("gene_symbol").size().to_dict() != {g: 8 for g in GENE_ORDER}:
        raise ValueError("panel e does not contain three genes by eight donors")
    return {
        "umap": umap, "composition": composition, "donor_order": donor_order,
        "pathway_display": pathway_display, "gene_display": gene_display, "donor_dots": donor_dots,
    }


def write_source_data(data: dict[str, pd.DataFrame | int]) -> None:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    for key, filename in [
        ("umap", "Figure3a_UMAP_descriptive_source_data.tsv.gz"),
        ("composition", "Figure3b_donor_celltype_composition_source_data.tsv"),
        ("pathway_display", "Figure3c_donor_pseudobulk_pathway_source_data.tsv"),
        ("gene_display", "Figure3d_gene_celltype_effect_source_data.tsv"),
        ("donor_dots", "Figure3e_fibroblast_donor_dot_source_data.tsv"),
    ]:
        frame = data[key]
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(SOURCE_OUT / filename, sep="\t", index=False, compression="infer")


def make_figure(data: dict[str, pd.DataFrame | int]) -> None:
    umap = data["umap"]
    composition = data["composition"]
    donor_order = data["donor_order"]
    pathways = data["pathway_display"]
    genes = data["gene_display"]
    dots = data["donor_dots"]
    assert all(isinstance(x, pd.DataFrame) for x in [umap, composition, donor_order, pathways, genes, dots])

    fig = plt.figure(figsize=(7.205, 6.102), constrained_layout=False)
    grid = fig.add_gridspec(3, 6, height_ratios=[1.05, 1.05, 1.15], width_ratios=[1, 1, 1, 1, 1, 1])
    fig.subplots_adjust(left=0.130, right=0.955, top=0.955, bottom=0.090, wspace=1.10, hspace=0.86)
    ax_a = fig.add_subplot(grid[0, :3])
    ax_b = fig.add_subplot(grid[0, 3:])
    ax_c = fig.add_subplot(grid[1:, :4])
    ax_d = fig.add_subplot(grid[1, 4:])
    e_grid = grid[2, 4:].subgridspec(1, 3, wspace=0.12)
    ax_e = [fig.add_subplot(e_grid[0, i]) for i in range(3)]

    # a — descriptive, non-inferential UMAP of all cells.
    for cell_type in ALL_CELLTYPE_ORDER:
        d = umap.loc[umap["cell_type"].eq(cell_type)]
        if d.empty:
            continue
        ax_a.scatter(d["UMAP1"], d["UMAP2"], s=0.34, color=CELL_COLORS[cell_type], linewidths=0, alpha=0.68,
                     rasterized=True, label=CELLTYPE_LABEL[cell_type])
    ax_a.set_xlabel("UMAP1", fontsize=5.1, labelpad=1)
    ax_a.set_ylabel("UMAP2", fontsize=5.1, labelpad=1)
    ax_a.tick_params(labelsize=4.6, length=2, pad=1)
    ax_a.set_title("ED corpus-cavernosum atlas", fontsize=6.2, fontweight="bold", loc="left", pad=2)
    ax_a.text(0.02, 0.03, "64,993 cells; descriptive only", transform=ax_a.transAxes, fontsize=4.2,
              bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.0})
    ax_a.legend(loc="upper left", ncol=2, fontsize=3.25, markerscale=3.5, handletextpad=0.2,
                columnspacing=0.55, borderpad=0.2, labelspacing=0.25)
    panel_label(ax_a, "a")

    # b — donor-composition plot; no cell-level statistical test.
    x = np.arange(donor_order.shape[0])
    bottom = np.zeros(donor_order.shape[0])
    for cell_type in ALL_CELLTYPE_ORDER:
        values = []
        for accession in donor_order["sample_accession"]:
            value = composition.loc[(composition["sample_accession"].eq(accession)) & (composition["cell_type"].eq(cell_type)), "cell_fraction"]
            values.append(float(value.iloc[0]) if value.shape[0] else 0.0)
        ax_b.bar(x, values, bottom=bottom, width=0.78, color=CELL_COLORS[cell_type], edgecolor="white", linewidth=0.20)
        bottom += np.asarray(values)
    ax_b.set_ylim(0, 1.0)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(donor_order["donor_label"], fontsize=4.3)
    ax_b.set_ylabel("Cell fraction", fontsize=5.1, labelpad=1)
    ax_b.tick_params(axis="y", labelsize=4.5, length=2, pad=1)
    ax_b.axvline(2.5, color="#9A9A9A", linewidth=0.55, linestyle="--")
    ax_b.axvline(5.5, color="#9A9A9A", linewidth=0.55, linestyle="--")
    for xpos, text in [(1, "Tumour-margin\nreference"), (4, "Non-diabetic ED\nimplantation"), (6.5, "Diabetic ED\nimplantation")]:
        ax_b.text(xpos, -0.26, text, ha="center", va="top", transform=ax_b.get_xaxis_transform(), fontsize=3.6, color="#5F5F5F")
    ax_b.set_title("Donor-level cell composition (n = 3/3/2)", fontsize=5.7, fontweight="bold", loc="left", pad=2)
    panel_label(ax_b, "b")

    # c — all pathways with pre-existing cross-context-significant rows, no new test.
    pathway_order = pathways["pathway"].drop_duplicates().tolist()
    matrix = pathways.pivot(index="pathway", columns="cell_type", values="NES").reindex(index=pathway_order, columns=CELLTYPE_ORDER)
    flag = pathways.pivot(index="pathway", columns="cell_type", values="frozen_display_flag").reindex(index=pathway_order, columns=CELLTYPE_ORDER).astype(bool)
    limit = max(2.0, float(np.ceil(np.nanmax(np.abs(matrix.to_numpy(dtype=float))) * 2) / 2))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#F4F4F4")
    im = ax_c.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
    yy, xx = np.where(flag.to_numpy())
    ax_c.scatter(xx, yy, s=11, facecolors="none", edgecolors="white", linewidths=0.70, zorder=3)
    ax_c.set_xticks(range(len(CELLTYPE_ORDER)))
    ax_c.set_xticklabels([CELLTYPE_LABEL[c] for c in CELLTYPE_ORDER], fontsize=4.35, rotation=35, ha="right")
    ax_c.set_yticks(range(len(pathway_order)))
    ax_c.set_yticklabels([PATHWAY_LABEL[p] for p in pathway_order], fontsize=4.5)
    ax_c.tick_params(length=0, pad=1)
    for spine in ax_c.spines.values():
        spine.set_visible(False)
    cb_c = fig.colorbar(im, ax=ax_c, fraction=0.032, pad=0.035)
    cb_c.set_label("ED NES", fontsize=4.7)
    cb_c.ax.tick_params(labelsize=4.1, length=2)
    ax_c.set_title("Donor-pseudobulk IH-linked programme localization", fontsize=6.4, fontweight="bold", loc="left", pad=4)
    panel_label(ax_c, "c")

    # d — all target-gene compartment effects; unavailable estimates remain blank.
    gene_matrix = genes.pivot(index="gene_symbol", columns="cell_type", values="ED_logFC").reindex(index=GENE_ORDER, columns=CELLTYPE_ORDER)
    gene_fdr = genes.pivot(index="gene_symbol", columns="cell_type", values="FDR_lt_0_05").reindex(index=GENE_ORDER, columns=CELLTYPE_ORDER).fillna(False).astype(bool)
    gene_available = genes.pivot(index="gene_symbol", columns="cell_type", values="estimate_available").reindex(index=GENE_ORDER, columns=CELLTYPE_ORDER).fillna(False).astype(bool)
    g_limit = max(1.0, float(np.ceil(np.nanmax(np.abs(gene_matrix.to_numpy(dtype=float))) * 2) / 2))
    im_g = ax_d.imshow(gene_matrix.to_numpy(dtype=float), aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-g_limit, vcenter=0, vmax=g_limit))
    yy, xx = np.where(gene_fdr.to_numpy())
    ax_d.scatter(xx, yy, s=22, facecolors="none", edgecolors="white", linewidths=0.9, zorder=3)
    ax_d.set_xticks(range(len(CELLTYPE_ORDER)))
    ax_d.set_xticklabels([CELLTYPE_LABEL[c] for c in CELLTYPE_ORDER], fontsize=3.9, rotation=38, ha="right")
    ax_d.set_yticks(range(len(GENE_ORDER)))
    ax_d.set_yticklabels(GENE_ORDER, fontsize=5.0, fontweight="bold")
    ax_d.tick_params(length=0, pad=1)
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    cb_g = fig.colorbar(im_g, ax=ax_d, fraction=0.050, pad=0.08)
    cb_g.set_label("ED log2FC", fontsize=4.5)
    cb_g.ax.tick_params(labelsize=4.0, length=2)
    ax_d.set_title("Prioritized-gene compartment effects", fontsize=5.8, fontweight="bold", loc="left", pad=4)
    ax_d.text(0.0, -0.53, "White ring: ED FDR < 0.05; blank: no estimate.", transform=ax_d.transAxes,
              fontsize=3.8, color="#5F5F5F", va="top")
    panel_label(ax_d, "d", x=-0.22)

    # e — ruled descriptive donor dots from unambiguously mapped fibroblast pseudobulk columns.
    y_limits = (float(np.floor(dots["log2CPM_prior_0_5"].min() - 0.4)), float(np.ceil(dots["log2CPM_prior_0_5"].max() + 0.4)))
    for ax, gene in zip(ax_e, GENE_ORDER):
        d = dots.loc[dots["gene_symbol"].eq(gene)].copy()
        for xpos, condition in enumerate(GROUP_ORDER):
            group = d.loc[d["condition"].eq(condition)].sort_values("sample_accession")
            offsets = np.linspace(-0.11, 0.11, group.shape[0]) if group.shape[0] > 1 else np.array([0.0])
            ax.scatter(np.full(group.shape[0], xpos) + offsets, group["log2CPM_prior_0_5"], s=15,
                       color=GROUP_COLOR[condition], edgecolor="white", linewidth=0.35, zorder=3)
        ax.set_title(gene, fontsize=5.0, fontweight="bold", pad=2)
        ax.set_ylim(*y_limits)
        ax.set_xticks(range(3))
        ax.set_xticklabels(["Ref\n(n=3)", "non-DM\n(n=3)", "DM\n(n=2)"], fontsize=3.55)
        ax.tick_params(axis="x", length=0, pad=1)
        ax.tick_params(axis="y", labelsize=3.9, length=2, pad=1)
        ax.grid(axis="y", color="#E3E3E3", linewidth=0.45, zorder=0)
        if gene == GENE_ORDER[0]:
            ax.set_ylabel("log2CPM", fontsize=4.4, labelpad=1)
        else:
            ax.set_yticklabels([])
    ax_e[0].text(0.0, 1.18, "Fibroblast donor pseudobulk", transform=ax_e[0].transAxes, fontsize=5.2, fontweight="bold", va="bottom")
    ax_e[0].text(0.0, -0.35, "0.5 prior count; descriptive only; no P values.", transform=ax_e[0].transAxes, fontsize=3.65, color="#5F5F5F", va="top")
    panel_label(ax_e[0], "e", x=-0.40, y=1.23)

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "Figure3_ED_corpus_cavernosum_localization"
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
        "figure": "Figure 3",
        "backend": "Python/matplotlib",
        "input_manifest_validation": "PASS",
        "input_sha256": verified,
        "assertions": {
            "umap_cells": 64993,
            "donors_by_group": {"tumour_margin_reference": 3, "non_diabetic_ED_implantation": 3, "diabetic_ED_implantation": 2},
            "pathway_programmes_shown": 8,
            "pathway_compartments_shown": 7,
            "target_gene_cells": 21,
            "ED_FDR_lt_0_05_target_effects": 3,
            "panel_e_mapping": "PASS_UNAMBIGUOUS",
            "panel_e_donor_dots": 24,
            "panel_e_normalization": "log2(((count + 0.5)/(library size + 1))*1e6)",
        },
        "output_dimensions_mm": [183, 155],
        "interpretation_boundary": "Cells are not independent replicates; UMAP and composition are descriptive, donor group is inseparable from procurement context, and this atlas does not establish causality or independent validation.",
    }
    (OUT / "Figure3_build_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
