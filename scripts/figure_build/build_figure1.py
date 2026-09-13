from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


# Required publication settings: editable SVG text and unified journal typography.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = OUTPUT_ROOT.parent
OUT = OUTPUT_ROOT / "figures" / "Figure1_20260911"
SOURCE_OUT = OUT / "source_data"
INPUT_MANIFEST = OUTPUT_ROOT / "figure_contracts" / "input_manifest.sha256.tsv"

SOURCE_IDS = {
    "S010": "revision_v2/02_results/endothelial_three_study_summary_v2.tsv",
    "S012": "revision_v2/02_results/Hallmark_clinical_bridge_summary_v2.tsv",
    "S013": "revision_v2/02_results/Hallmark_clinical_bridge_v2.tsv",
    "S014": "revision_v2/02_results/figure_source_data/Figure1_Hallmark_bridge_source_data_v2.tsv",
}

COLORS = {
    "ink": "#272727",
    "muted": "#767676",
    "border": "#4D4D4D",
    "neutral": "#E8E8E8",
    "ih": "#D9EAF7",
    "ih_dark": "#0F4D92",
    "clinical": "#E9DDF1",
    "clinical_dark": "#7A4C96",
    "ed": "#F3E3D5",
    "ed_dark": "#A45D2C",
    "model": "#DDEFF0",
    "model_dark": "#287F88",
    "strict": "#9A4D8E",
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
            raise ValueError(f"{source_id} missing from source manifest")
        declared_path = manifest.loc[source_id, "relative_path"]
        if declared_path != relative_path:
            raise ValueError(f"{source_id} manifest path differs from figure contract")
        path = PROJECT_ROOT / relative_path
        observed = sha256(path)
        expected = manifest.loc[source_id, "sha256"]
        if observed.lower() != expected.lower():
            raise ValueError(f"input drift for {source_id}")
        verified[source_id] = observed
    return verified


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.03, label, transform=ax.transAxes, fontsize=8.2, fontweight="bold",
            ha="left", va="bottom", color=COLORS["ink"])


def rounded_box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float,
                title: str, body: str, fill: str, edge: str, title_size: float = 6.6,
                body_size: float = 5.4) -> None:
    x, y = xy
    ax.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.018,rounding_size=0.018",
                                facecolor=fill, edgecolor=edge, linewidth=0.8, zorder=2))
    ax.text(x + width / 2, y + height * 0.66, title, ha="center", va="center",
            fontsize=title_size, fontweight="bold", color=COLORS["ink"], zorder=3)
    ax.text(x + width / 2, y + height * 0.32, body, ha="center", va="center",
            fontsize=body_size, linespacing=1.03, color=COLORS["ink"], zorder=3)


def write_source_data(endothelial_summary: pd.DataFrame, bridge_summary: pd.DataFrame,
                      heat_source: pd.DataFrame) -> None:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    evidence_layers = pd.DataFrame([
        {
            "layer": "IH endothelial discovery",
            "datasets_or_context": "GSE243023; GSE205050; GSE10723",
            "independent_unit": "study-specific bulk samples",
            "role": "exposure-associated discovery",
            "causal_interpretation": "not a causal link",
        },
        {
            "layer": "OSA/CPAP clinical context",
            "datasets_or_context": "OSA versus control (AM); CPAP versus OSA (AM)",
            "independent_unit": "participants in original clinical comparisons",
            "role": "strict pathway bridge",
            "causal_interpretation": "not an ED cohort",
        },
        {
            "layer": "ED corpus-cavernosum target tissue",
            "datasets_or_context": "GSE206528 donor pseudobulk",
            "independent_unit": "donor",
            "role": "target-tissue localization",
            "causal_interpretation": "procurement context is group-confounded",
        },
        {
            "layer": "Model-based prioritization",
            "datasets_or_context": "predeclared virtual perturbation",
            "independent_unit": "model output; seeds are not biological replicates",
            "role": "hypothesis prioritization",
            "causal_interpretation": "not experimental validation",
        },
    ])
    evidence_layers.to_csv(SOURCE_OUT / "Figure1a_evidence_layers.tsv", sep="\t", index=False)
    endothelial_summary.to_csv(SOURCE_OUT / "Figure1b_endothelial_filters.tsv", sep="\t", index=False)
    bridge_summary.to_csv(SOURCE_OUT / "Figure1c_clinical_bridge_filters.tsv", sep="\t", index=False)
    heat_source.to_csv(SOURCE_OUT / "Figure1d_strict_context_heatmap.tsv", sep="\t", index=False)


def make_figure(endothelial_summary: pd.DataFrame, bridge_summary: pd.DataFrame,
                complete: pd.DataFrame, heat_source: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(7.205, 4.724), constrained_layout=False)
    grid = fig.add_gridspec(2, 3, width_ratios=[0.92, 0.92, 1.45], height_ratios=[0.86, 1.14])
    fig.subplots_adjust(left=0.075, right=0.945, top=0.92, bottom=0.135, wspace=0.48, hspace=0.55)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :2])
    ax_d = fig.add_subplot(grid[:, 2])

    # a — parallel evidence layers. Connectors have no arrowheads by design.
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis("off")
    rounded_box(ax_a, (0.02, 0.65), 0.30, 0.20, "IH endothelium", "3 datasets\nexposure", COLORS["ih"], COLORS["ih_dark"], 5.3, 4.6)
    rounded_box(ax_a, (0.35, 0.65), 0.30, 0.20, "OSA / CPAP", "clinical\nbridge", COLORS["clinical"], COLORS["clinical_dark"], 5.3, 4.6)
    rounded_box(ax_a, (0.68, 0.65), 0.30, 0.20, "ED tissue", "donor-pseudobulk\nlocalization", COLORS["ed"], COLORS["ed_dark"], 5.3, 4.4)
    rounded_box(ax_a, (0.24, 0.26), 0.52, 0.18, "Evidence triangulation", "parallel support; no causal sequence", COLORS["model"], COLORS["model_dark"], 5.8, 4.7)
    for x in (0.17, 0.50, 0.83):
        ax_a.plot([x, 0.50], [0.65, 0.44], color=COLORS["muted"], linewidth=0.8, zorder=1)
    ax_a.text(0.50, 0.02, "ED atlas reference: tumour-margin context", ha="center", va="bottom", fontsize=4.7,
              color=COLORS["muted"])
    add_panel_label(ax_a, "a")

    # b — two explicitly separate filtering lanes.
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.axis("off")
    endothelial = dict(zip(endothelial_summary["metric"], endothelial_summary["value"]))
    upper = [
        ("5,139", "common\ngenes", COLORS["neutral"]),
        ("1,282", "full-meta\nFDR < 0.05", COLORS["ih"]),
    ]
    lower = [
        ("1,766", "all three\nsame direction", "#EFF4F8"),
        ("1,035", "directional\nstability", "#E4EFF8"),
        ("258", "all-three\nLODO FDR", COLORS["ih"]),
    ]
    expected_values = [
        endothelial["common_genes"], endothelial["full_meta_FDR_lt_0_05"],
        endothelial["three_study_all_same_direction"], endothelial["directional_stability_screen"],
        endothelial["strict_all_three_lodo_FDR_lt_0_05"],
    ]
    if expected_values != [5139, 1282, 1766, 1035, 258]:
        raise ValueError("frozen endothelial summary differs from Figure 1 contract")
    for i, (number, label, fill) in enumerate(upper):
        x = 0.05 + i * 0.48
        rounded_box(ax_b, (x, 0.60), 0.40, 0.27, number, label, fill, COLORS["ih_dark"], 8.0, 5.0)
        if i == 0:
            ax_b.plot([0.45, 0.53], [0.735, 0.735], color=COLORS["muted"], linewidth=0.8)
    for i, (number, label, fill) in enumerate(lower):
        x = 0.02 + i * 0.325
        rounded_box(ax_b, (x, 0.21), 0.29, 0.25, number, label, fill, COLORS["ih_dark"], 7.0, 4.6)
        if i < 2:
            ax_b.plot([x + 0.29, x + 0.325], [0.335, 0.335], color=COLORS["muted"], linewidth=0.8)
    ax_b.text(0.50, 0.03, "Separate evidence filters—not a causal sequence", ha="center", va="bottom",
              fontsize=4.8, style="italic", color=COLORS["muted"])
    add_panel_label(ax_b, "b")

    # c — the predeclared clinical-bridge screen.
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.axis("off")
    bridge = dict(zip(bridge_summary["metric"], bridge_summary["value"]))
    bridge_steps = [
        (str(bridge["hallmark_pathways_total"]), "Hallmark pathways\nevaluated", COLORS["neutral"]),
        (str(bridge["all_four_estimated"]), "all four contexts\nestimated", "#EEF2F6"),
        (str(bridge["complete_directional_pattern"]), "complete directional\npattern", COLORS["ih"]),
        (str(bridge["all_four_FDR_pass"]), "all-four\nFDR < 0.05", COLORS["clinical"]),
        (str(bridge["strict_clinical_bridge_pass"]), "strict bridge\nHALLMARK_HYPOXIA", "#EAD5E6"),
    ]
    expected_bridge = [50, 49, 12, 3, 1]
    observed_bridge = [
        bridge["hallmark_pathways_total"], bridge["all_four_estimated"], bridge["complete_directional_pattern"],
        bridge["all_four_FDR_pass"], bridge["strict_clinical_bridge_pass"],
    ]
    if observed_bridge != expected_bridge:
        raise ValueError("frozen bridge summary differs from Figure 1 contract")
    for i, (number, label, fill) in enumerate(bridge_steps):
        x = 0.025 + i * 0.197
        rounded_box(ax_c, (x, 0.34), 0.15, 0.34, number, label, fill,
                    COLORS["strict"] if i == 4 else COLORS["border"], 8.0, 4.9)
        if i < len(bridge_steps) - 1:
            ax_c.plot([x + 0.15, x + 0.197], [0.51, 0.51], color=COLORS["muted"], linewidth=0.8)
    ax_c.text(0.50, 0.08, "Prespecified four-context screen; the final rule is more restrictive than FDR alone.",
              ha="center", va="bottom", fontsize=5.2, color=COLORS["muted"], style="italic")
    add_panel_label(ax_c, "c")

    # d — hero heatmap, with FDR stars and strict-pathway outline.
    contexts = [
        ("endo243", "Endothelial\nIH 1"),
        ("endo205", "Endothelial\nIH 2"),
        ("osaAM", "OSA vs\ncontrol (AM)"),
        ("cpapAM", "CPAP vs\nOSA (AM)"),
    ]
    matrix = np.array([[row[f"{key}_NES"] for key, _ in contexts] for _, row in complete.iterrows()], dtype=float)
    max_abs = max(2.5, float(np.abs(matrix).max()))
    image = ax_d.imshow(matrix, cmap="RdBu_r", aspect="auto",
                        norm=TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs))
    ax_d.set_xticks(range(len(contexts)))
    ax_d.set_xticklabels([label for _, label in contexts], fontsize=5.0)
    ax_d.xaxis.tick_top()
    ax_d.tick_params(axis="x", length=0, pad=2)
    labels = [str(pathway).replace("HALLMARK_", "").replace("_", " ") for pathway in complete["pathway"]]
    ax_d.set_yticks(range(len(labels)))
    ax_d.set_yticklabels(labels, fontsize=4.7)
    ax_d.tick_params(axis="y", length=0, pad=2)
    for i, (_, row) in enumerate(complete.iterrows()):
        for j, (key, _) in enumerate(contexts):
            if float(row[f"{key}_FDR"]) < 0.05:
                ax_d.text(j, i, "*", ha="center", va="center", fontsize=7.2, fontweight="bold", color=COLORS["ink"])
    strict_index = int(np.flatnonzero(complete["strict_clinical_bridge_pass"].to_numpy(dtype=bool))[0])
    ax_d.add_patch(Rectangle((-0.5, strict_index - 0.5), len(contexts), 1,
                             fill=False, edgecolor=COLORS["strict"], linewidth=1.5, clip_on=False))
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax_d, fraction=0.036, pad=0.04)
    colorbar.set_label("NES", fontsize=6.2)
    colorbar.ax.tick_params(labelsize=5.5, length=2)
    ax_d.text(0.0, -0.11, "* FDR < 0.05; outline = sole strict bridge", transform=ax_d.transAxes,
              fontsize=4.8, color=COLORS["muted"], va="top")
    add_panel_label(ax_d, "d")

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "Figure1_parallel_evidence_and_strict_bridge"
    fig.savefig(base.with_suffix(".svg"), facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white")
    fig.savefig(base.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main() -> None:
    verified = verify_inputs()
    endothelial_summary = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S010"], sep="\t")
    bridge_summary = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S012"], sep="\t")
    bridge_full = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S013"], sep="\t")
    existing_long = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S014"], sep="\t")

    complete = bridge_full.loc[bridge_full["complete_directional_pattern"].astype(bool)].copy()
    complete["strict_clinical_bridge_pass"] = complete["strict_clinical_bridge_pass"].astype(bool)
    complete = complete.sort_values(["strict_clinical_bridge_pass", "pathway"], ascending=[False, True]).reset_index(drop=True)
    if len(complete) != 12 or int(complete["strict_clinical_bridge_pass"].sum()) != 1:
        raise ValueError("expected 12 complete-direction pathways and one strict bridge")

    contexts = [
        ("endo243", "Endothelial IH dataset 1"),
        ("endo205", "Endothelial IH dataset 2"),
        ("osaAM", "OSA vs control (AM)"),
        ("cpapAM", "CPAP vs OSA (AM)"),
    ]
    source_rows = []
    for _, row in complete.iterrows():
        for key, label in contexts:
            source_rows.append({
                "pathway": row["pathway"],
                "revised_bridge_class": row["revised_bridge_class"],
                "strict_clinical_bridge_pass": bool(row["strict_clinical_bridge_pass"]),
                "context": key,
                "context_label": label,
                "NES": float(row[f"{key}_NES"]),
                "FDR": float(row[f"{key}_FDR"]),
                "FDR_lt_0_05_star": bool(row[f"{key}_FDR"] < 0.05),
            })
    heat_source = pd.DataFrame(source_rows)
    canonical = heat_source.sort_values(["pathway", "context"]).reset_index(drop=True)
    prior = existing_long.sort_values(["pathway", "context"]).reset_index(drop=True)
    compare_columns = ["pathway", "revised_bridge_class", "strict_clinical_bridge_pass", "context", "NES", "FDR", "FDR_lt_0_05_star"]
    pd.testing.assert_frame_equal(
        canonical[compare_columns], prior[compare_columns], check_dtype=False, check_exact=False, atol=1e-12, rtol=0
    )

    write_source_data(endothelial_summary, bridge_summary, heat_source)
    make_figure(endothelial_summary, bridge_summary, complete, heat_source)
    audit = {
        "status": "PASS",
        "figure": "Figure 1",
        "backend": "Python/matplotlib",
        "input_manifest_validation": "PASS",
        "input_sha256": verified,
        "assertions": {
            "endothelial_summary_counts": [5139, 1282, 1766, 1035, 258],
            "bridge_summary_counts": [50, 49, 12, 3, 1],
            "complete_directional_pathways": 12,
            "strict_bridge_pathways": 1,
            "long_source_equivalence_to_prior_sealed_source": "PASS",
        },
        "output_dimensions_mm": [183, 120],
        "interpretation_boundary": "Parallel evidence layers and a strict clinical bridge do not establish OSA–ED causal mediation.",
    }
    (OUT / "Figure1_build_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
