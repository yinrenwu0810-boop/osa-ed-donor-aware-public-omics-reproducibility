from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
OUT = OUTPUT_ROOT / "figures" / "Figure5_20260912"
SOURCE_OUT = OUT / "source_data"
INPUT_MANIFEST = OUTPUT_ROOT / "figure_contracts" / "input_manifest.sha256.tsv"
SOURCE_IDS = {
    "S050": "analysis/results/phase2/machine_learning/nested_LODO_predictions.tsv",
    "S051": "analysis/results/phase2/machine_learning/nested_LODO_pooled_metrics.tsv",
    "S052": "analysis/results/phase2/machine_learning/nested_LODO_per_dataset_metrics.tsv",
    "S053": "analysis/results/phase2/machine_learning/feature_stability.tsv",
    "S054": "analysis/results/phase2/machine_learning/selected_features_by_outer_fold.tsv",
}
STRATEGY_ORDER = ["gene_rank", "Hallmark_rank_score", "adaptive"]
STRATEGY_LABEL = {"gene_rank": "Gene rank", "Hallmark_rank_score": "Hallmark rank", "adaptive": "Adaptive"}
STRATEGY_COLOR = {"gene_rank": "#3775BA", "Hallmark_rank_score": "#E69F00", "adaptive": "#8A6BB8"}
DATASET_ORDER = ["GSE205050", "GSE243023", "GSE10723"]
DATASET_COLOR = {"GSE205050": "#C65B50", "GSE243023": "#009E73", "GSE10723": "#6E8FC7"}


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


def panel_label(ax: plt.Axes, label: str, x: float = -0.13, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8.2, fontweight="bold", ha="left", va="bottom", color="#272727")


def roc_curve_from_predictions(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    ordered = frame.sort_values("probability_IH", ascending=False, kind="mergesort")
    labels = ordered["true_label"].to_numpy(dtype=int)
    positives = labels.sum()
    negatives = labels.shape[0] - positives
    if positives == 0 or negatives == 0:
        raise ValueError("ROC requires both labels")
    tp = np.cumsum(labels == 1)
    fp = np.cumsum(labels == 0)
    return np.r_[0.0, fp / negatives, 1.0], np.r_[0.0, tp / positives, 1.0]


def feature_label(value: str) -> str:
    labels = {
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNFα signaling via NF-κB",
        "HALLMARK_MYC_TARGETS_V2": "MYC targets v2",
        "HALLMARK_E2F_TARGETS": "E2F targets",
        "HALLMARK_UV_RESPONSE_DN": "UV response DN",
        "HALLMARK_HYPOXIA": "Hypoxia",
        "HALLMARK_IL2_STAT5_SIGNALING": "IL2–STAT5 signaling",
        "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "Epithelial–mesenchymal transition",
        "HALLMARK_PROTEIN_SECRETION": "Protein secretion",
        "HALLMARK_COMPLEMENT": "Complement",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "WNT/β-catenin signaling",
        "HALLMARK_TGF_BETA_SIGNALING": "TGF-β signaling",
        "HALLMARK_FATTY_ACID_METABOLISM": "Fatty acid metabolism",
    }
    return labels.get(value, value.replace("HALLMARK_", "").replace("_", " "))


def select_inputs() -> dict[str, pd.DataFrame]:
    predictions = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S050"], sep="\t")
    required_predictions = {"strategy", "sample_accession", "dataset", "true_label", "probability_IH"}
    if not required_predictions.issubset(predictions.columns):
        raise ValueError("prediction table lacks required columns")
    if predictions.shape[0] != 60 or predictions.groupby("strategy").size().to_dict() != {s: 20 for s in STRATEGY_ORDER}:
        raise ValueError("prediction table does not retain 20 held-out samples per strategy")
    predictions["strategy_label"] = predictions["strategy"].map(STRATEGY_LABEL)
    predictions["true_label_display"] = predictions["true_label"].map({0: "Control", 1: "IH"})
    predictions["display_group"] = predictions["dataset"] + "\n" + predictions["true_label_display"]

    pooled = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S051"], sep="\t").set_index("strategy").reindex(STRATEGY_ORDER).reset_index()
    if pooled["strategy"].tolist() != STRATEGY_ORDER or pooled["n"].tolist() != [20.0, 20.0, 20.0]:
        raise ValueError("pooled metric table does not retain three strategies and n=20")
    expected_auc = [0.80, 0.62, 0.54]
    if not np.allclose(pooled["ROC_AUC"], expected_auc):
        raise ValueError("pooled AUC values differ from frozen metrics")
    if not (pooled["exact_permutations"] == 28000).all():
        raise ValueError("exact permutation count differs from frozen metrics")
    pooled["strategy_label"] = pooled["strategy"].map(STRATEGY_LABEL)

    roc_rows = []
    for strategy in STRATEGY_ORDER:
        fpr, tpr = roc_curve_from_predictions(predictions.loc[predictions["strategy"].eq(strategy)])
        calculated_auc = float(np.trapezoid(tpr, fpr))
        frozen_auc = float(pooled.loc[pooled["strategy"].eq(strategy), "ROC_AUC"].iloc[0])
        if abs(calculated_auc - frozen_auc) > 1e-12:
            raise ValueError(f"ROC derivation does not reproduce frozen AUC for {strategy}")
        roc_rows.extend({"strategy": strategy, "strategy_label": STRATEGY_LABEL[strategy], "FPR": x, "TPR": y, "frozen_ROC_AUC": frozen_auc} for x, y in zip(fpr, tpr))
    roc = pd.DataFrame(roc_rows)

    per_dataset = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S052"], sep="\t")
    per_dataset = per_dataset.set_index(["strategy", "dataset"]).reindex(pd.MultiIndex.from_product([STRATEGY_ORDER, DATASET_ORDER], names=["strategy", "dataset"])).reset_index()
    if per_dataset.shape[0] != 9 or per_dataset["ROC_AUC"].isna().any():
        raise ValueError("per-dataset AUC table is incomplete")
    gene_auc = per_dataset.loc[per_dataset["strategy"].eq("gene_rank"), "ROC_AUC"].tolist()
    if not np.allclose(gene_auc, [0.1111111111111111, 1.0, 0.9375]):
        raise ValueError("gene-rank heterogeneity values differ from frozen metrics")
    per_dataset["strategy_label"] = per_dataset["strategy"].map(STRATEGY_LABEL)

    stability = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S053"], sep="\t")
    selected_features = pd.read_csv(PROJECT_ROOT / SOURCE_IDS["S054"], sep="\t")
    if stability.shape[0] != 5240 or selected_features.shape[0] != 10452:
        raise ValueError("feature source table dimensions differ from frozen files")
    top_features = stability.loc[stability["selected_outer_folds"].ge(2)].copy()
    top_features = top_features.sort_values(["selected_outer_folds", "median_abs_coefficient", "feature"], ascending=[False, False, True], kind="mergesort").head(20).reset_index(drop=True)
    if top_features.shape[0] != 20:
        raise ValueError("deterministic feature selection did not produce 20 rows")
    top_features["display_order"] = np.arange(top_features.shape[0])
    top_features["feature_label"] = top_features["feature"].map(feature_label)
    duplicate_labels = top_features["feature_label"].value_counts()
    top_features["feature_display_label"] = [
        label if duplicate_labels[label] == 1 else f"{label} [{STRATEGY_LABEL[strategy]}]"
        for label, strategy in zip(top_features["feature_label"], top_features["strategy"])
    ]
    top_features["selection_rule"] = "selected_outer_folds >= 2; then selected_outer_folds descending, median_abs_coefficient descending, feature alphabetical; capped at 20"
    if top_features["strategy"].value_counts().to_dict() != {"Hallmark_rank_score": 13, "gene_rank": 5, "adaptive": 2}:
        raise ValueError("top feature strategy composition differs from frozen deterministic ordering")
    return {"predictions": predictions, "pooled": pooled, "roc": roc, "per_dataset": per_dataset, "top_features": top_features}


def write_source_data(data: dict[str, pd.DataFrame]) -> None:
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    files = {
        "predictions": "Figure5a_heldout_probability_source_data.tsv",
        "roc": "Figure5b_pooled_ROC_source_data.tsv",
        "per_dataset": "Figure5c_per_dataset_AUC_source_data.tsv",
        "top_features": "Figure5d_feature_stability_source_data.tsv",
        "pooled": "Figure5e_pooled_metric_source_data.tsv",
    }
    for key, filename in files.items():
        data[key].to_csv(SOURCE_OUT / filename, sep="\t", index=False)


def make_figure(data: dict[str, pd.DataFrame]) -> None:
    predictions = data["predictions"]
    pooled = data["pooled"]
    roc = data["roc"]
    per_dataset = data["per_dataset"]
    features = data["top_features"]
    fig = plt.figure(figsize=(7.205, 5.315), constrained_layout=False)
    grid = fig.add_gridspec(2, 6, height_ratios=[0.93, 1.48], width_ratios=[1, 1, 1, 1, 1, 1])
    fig.subplots_adjust(left=0.112, right=0.958, top=0.950, bottom=0.105, wspace=0.90, hspace=0.82)
    ax_a = fig.add_subplot(grid[0, :2])
    ax_b = fig.add_subplot(grid[0, 2:4])
    ax_c = fig.add_subplot(grid[0, 4:])
    ax_d = fig.add_subplot(grid[1, :4])
    ax_e = fig.add_subplot(grid[1, 4:])

    # a — all frozen held-out probabilities across strategies, datasets, and true labels.
    group_order = [(dataset, label) for dataset in DATASET_ORDER for label in [0, 1]]
    for gi, (dataset, label) in enumerate(group_order):
        for si, strategy in enumerate(STRATEGY_ORDER):
            d = predictions.loc[(predictions["dataset"].eq(dataset)) & (predictions["true_label"].eq(label)) & (predictions["strategy"].eq(strategy))].sort_values("sample_accession")
            offsets = np.linspace(-0.24, 0.24, d.shape[0]) if d.shape[0] > 1 else np.array([0.0])
            ax_a.scatter(np.full(d.shape[0], gi) + offsets + (si - 1) * 0.045, d["probability_IH"], s=9, color=STRATEGY_COLOR[strategy], alpha=0.80, linewidths=0, zorder=3)
    ax_a.axhline(0.5, color="#7C7C7C", linewidth=0.55, linestyle="--")
    ax_a.set_xlim(-0.55, 5.55)
    ax_a.set_ylim(-0.03, 1.03)
    ax_a.set_xticks(range(6))
    ax_a.set_xticklabels(["Ctrl", "IH", "Ctrl", "IH", "Ctrl", "IH"], fontsize=4.1)
    ax_a.set_ylabel("Held-out IH probability", fontsize=4.7, labelpad=1)
    ax_a.tick_params(axis="y", labelsize=4.2, length=2, pad=1)
    ax_a.tick_params(axis="x", length=0, pad=1)
    for boundary in [1.5, 3.5]:
        ax_a.axvline(boundary, color="#C5C5C5", linewidth=0.55)
    for xpos, dataset in zip([0.5, 2.5, 4.5], DATASET_ORDER):
        ax_a.text(xpos, -0.29, dataset, transform=ax_a.get_xaxis_transform(), ha="center", va="top", fontsize=3.8, color="#5F5F5F")
    ax_a.set_title("Held-out probabilities", fontsize=6.0, fontweight="bold", loc="left", pad=2)
    ax_a.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=STRATEGY_COLOR[s], markersize=3.7, label=STRATEGY_LABEL[s]) for s in STRATEGY_ORDER],
                loc="upper left", fontsize=3.35, borderpad=0.15, handletextpad=0.25, labelspacing=0.25)
    panel_label(ax_a, "a")

    # b — pooled ROC from frozen held-out predictions; AUCs copied from S051.
    ax_b.plot([0, 1], [0, 1], color="#8A8A8A", linewidth=0.55, linestyle="--", zorder=0)
    for strategy in STRATEGY_ORDER:
        d = roc.loc[roc["strategy"].eq(strategy)]
        auc = float(d["frozen_ROC_AUC"].iloc[0])
        ax_b.step(d["FPR"], d["TPR"], where="post", color=STRATEGY_COLOR[strategy], linewidth=1.0,
                  label=f"{STRATEGY_LABEL[strategy]} ({auc:.2f})")
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1.02)
    ax_b.set_xlabel("False-positive rate", fontsize=4.6, labelpad=1)
    ax_b.set_ylabel("True-positive rate", fontsize=4.6, labelpad=1)
    ax_b.tick_params(labelsize=4.0, length=2, pad=1)
    ax_b.set_title("Pooled held-out ROC", fontsize=6.0, fontweight="bold", loc="left", pad=2)
    ax_b.legend(loc="lower right", fontsize=3.3, borderpad=0.18, handlelength=1.2, labelspacing=0.25, handletextpad=0.3)
    panel_label(ax_b, "b", x=-0.20)

    # c — per-dataset AUC range to make heterogeneity prominent.
    y_positions = np.arange(len(STRATEGY_ORDER))
    for yi, strategy in enumerate(STRATEGY_ORDER):
        d = per_dataset.loc[per_dataset["strategy"].eq(strategy)].copy()
        ax_c.hlines(yi, d["ROC_AUC"].min(), d["ROC_AUC"].max(), color="#8A8A8A", linewidth=1.0, zorder=1)
        for row in d.itertuples(index=False):
            ax_c.scatter(row.ROC_AUC, yi, s=20, color=DATASET_COLOR[row.dataset], edgecolor="white", linewidth=0.35, zorder=2)
            if strategy == "gene_rank":
                if row.dataset == "GSE205050":
                    x_offset, y_offset, ha = 0.025, 0.18, "left"
                elif row.dataset == "GSE243023":
                    x_offset, y_offset, ha = -0.025, -0.22, "right"
                else:
                    x_offset, y_offset, ha = -0.025, 0.22, "right"
                ax_c.text(row.ROC_AUC + x_offset, yi + y_offset, f"{row.ROC_AUC:.2f}", ha=ha, va="center", fontsize=4.4, fontweight="bold")
    ax_c.set_xlim(0, 1.05)
    ax_c.set_ylim(-0.55, 2.55)
    ax_c.set_yticks(y_positions)
    ax_c.set_yticklabels([STRATEGY_LABEL[s] for s in STRATEGY_ORDER], fontsize=4.15)
    ax_c.invert_yaxis()
    ax_c.set_xlabel("Held-out dataset AUC", fontsize=4.6, labelpad=1)
    ax_c.tick_params(axis="x", labelsize=4.0, length=2, pad=1)
    ax_c.tick_params(axis="y", length=0, pad=1)
    ax_c.set_title("Transportability is heterogeneous", fontsize=6.0, fontweight="bold", loc="left", pad=2)
    ax_c.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=DATASET_COLOR[d], markersize=3.7, label=d) for d in DATASET_ORDER],
                loc="lower right", fontsize=3.15, borderpad=0.15, handletextpad=0.25, labelspacing=0.2)
    panel_label(ax_c, "c", x=-0.22)

    # d — capped deterministic stable feature set; no refitting.
    display = features.sort_values("display_order", ascending=False).copy()
    y = np.arange(display.shape[0])
    ax_d.hlines(y, 0, display["median_abs_coefficient"], color="#A8A8A8", linewidth=0.55, zorder=1)
    for strategy in STRATEGY_ORDER:
        d = display.loc[display["strategy"].eq(strategy)]
        ax_d.scatter(d["median_abs_coefficient"], d.index.map(lambda i: list(display.index).index(i)), s=15 + 7 * d["selected_outer_folds"],
                     color=STRATEGY_COLOR[strategy], edgecolor="white", linewidth=0.35, zorder=2, label=STRATEGY_LABEL[strategy])
    ax_d.set_yticks(y)
    ax_d.set_yticklabels(display["feature_display_label"], fontsize=3.45)
    ax_d.set_xlabel("Median absolute coefficient", fontsize=4.7, labelpad=1)
    ax_d.tick_params(axis="x", labelsize=4.0, length=2, pad=1)
    ax_d.tick_params(axis="y", length=0, pad=1)
    ax_d.set_title("Deterministic outer-fold feature display", fontsize=6.0, fontweight="bold", loc="left", pad=2)
    ax_d.text(0.0, -0.19, "All features selected in ≥2 outer folds; capped at 20 by fold count, coefficient, then feature name.",
              transform=ax_d.transAxes, fontsize=3.75, color="#5F5F5F", va="top")
    ax_d.legend(loc="lower right", fontsize=3.45, borderpad=0.15, handletextpad=0.25, labelspacing=0.25)
    panel_label(ax_d, "d")

    # e — conditional pooled uncertainty and structure-preserving exact permutation result, copied from S051.
    y = np.arange(len(STRATEGY_ORDER))
    for yi, row in enumerate(pooled.itertuples(index=False)):
        color = STRATEGY_COLOR[row.strategy]
        ax_e.hlines(yi, row.AUC_CI95_low, row.AUC_CI95_high, color=color, linewidth=1.6, zorder=1)
        ax_e.scatter(row.ROC_AUC, yi, s=24, color=color, edgecolor="white", linewidth=0.35, zorder=2)
        ax_e.text(min(row.AUC_CI95_high + 0.025, 1.0), yi, f"AUC {row.ROC_AUC:.2f}\nP = {row.external_permutation_P:.3g}", va="center", fontsize=3.65)
    ax_e.set_xlim(0.35, 1.02)
    ax_e.set_ylim(-0.55, 2.55)
    ax_e.set_yticks(y)
    ax_e.set_yticklabels([STRATEGY_LABEL[s] for s in STRATEGY_ORDER], fontsize=4.0)
    ax_e.invert_yaxis()
    ax_e.set_xlabel("Conditional pooled AUC (95% interval)", fontsize=4.35, labelpad=1)
    ax_e.tick_params(axis="x", labelsize=3.9, length=2, pad=1)
    ax_e.tick_params(axis="y", length=0, pad=1)
    ax_e.set_title("Conditional pooled uncertainty", fontsize=5.8, fontweight="bold", loc="left", pad=2)
    ax_e.text(0.0, -0.23, "n = 20; structure-preserving exact permutation, 28,000 permutations.", transform=ax_e.transAxes,
              fontsize=3.55, color="#5F5F5F", va="top")
    panel_label(ax_e, "e", x=-0.20)

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "Figure5_partial_transportability_and_heterogeneity"
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
        "status": "PASS", "figure": "Figure 5", "backend": "Python/matplotlib",
        "input_manifest_validation": "PASS", "input_sha256": verified,
        "assertions": {
            "heldout_predictions": 60, "heldout_samples_per_strategy": 20,
            "pooled_auc": {row.strategy: row.ROC_AUC for row in data["pooled"].itertuples(index=False)},
            "gene_rank_dataset_auc": data["per_dataset"].loc[data["per_dataset"]["strategy"].eq("gene_rank"), "ROC_AUC"].tolist(),
            "feature_display_count": 20,
            "feature_display_strategy_counts": data["top_features"]["strategy"].value_counts().to_dict(),
            "exact_permutations": 28000,
        },
        "output_dimensions_mm": [183, 135],
        "interpretation_boundary": "Pooled and per-dataset estimates are conditional on the assembled small datasets. Dataset-level heterogeneity prevents diagnostic, biomarker, external-clinical-validation, or broad transportability claims.",
    }
    (OUT / "Figure5_build_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
