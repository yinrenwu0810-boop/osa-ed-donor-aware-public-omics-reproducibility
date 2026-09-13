from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.8,
    "legend.frameon": False,
})

OUTPUT_ROOT = Path(__file__).resolve().parents[1]
ROOT = OUTPUT_ROOT.parent
OUT = OUTPUT_ROOT / "figures" / "Figure7_20260912"
SOURCE_DIR = OUT / "source_data"
MANIFEST = OUTPUT_ROOT / "figure_contracts" / "input_manifest.sha256.tsv"
SOURCES = {
    "S070": "revision_v2/07_virtual_perturbation/07_enrichment/v1/hypoxia_primary_summary.tsv",
    "S071": "revision_v2/07_virtual_perturbation/10_exploratory_hypoxia_sensitivity/attempt_20260910_03/06_analysis_multiverse.tsv",
    "S072": "revision_v2/07_virtual_perturbation/10_exploratory_hypoxia_sensitivity/attempt_20260910_02/03_target_donor_seed_hypoxia.tsv",
    "S073": "revision_v2/07_virtual_perturbation/10_exploratory_hypoxia_sensitivity/attempt_20260910_02/04_leave_one_seed_out.tsv",
    "S074": "revision_v2/07_virtual_perturbation/10_exploratory_hypoxia_sensitivity/attempt_20260910_02/04_leave_one_donor_out.tsv",
    "S075": "revision_v2/07_virtual_perturbation/10_exploratory_hypoxia_sensitivity/attempt_20260910_02/05_target_vs_matched_controls.tsv",
    "S076": "revision_v2/07_virtual_perturbation/10_exploratory_hypoxia_sensitivity/attempt_20260910_02/05_leading_edge_stability.tsv",
    "S077": "revision_v2/07_virtual_perturbation/10_exploratory_hypoxia_sensitivity/attempt_20260910_03/06_hypoxia_mapping_qc.tsv",
}
TARGETS = ["TYMS", "EFNB2", "LRRC17"]
COLORS = {"TYMS": "#2A9D8F", "EFNB2": "#E69F00", "LRRC17": "#8A6BB8"}
SPEC_COLORS = {"S0": "#4D4D4D", "S1": "#729ECE", "S2": "#C2A5CF"}


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def verify_inputs() -> dict[str, str]:
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str).set_index("source_id")
    hashes = {}
    for source_id, relative_path in SOURCES.items():
        path = ROOT / relative_path
        actual = digest(path)
        if manifest.loc[source_id, "relative_path"] != relative_path or actual.lower() != manifest.loc[source_id, "sha256"].lower():
            raise ValueError(f"input drift: {source_id}")
        hashes[source_id] = actual
    return hashes


def panel_label(ax, text: str, x: float = -0.15) -> None:
    ax.text(x, 1.06, text, transform=ax.transAxes, fontsize=8.5, fontweight="bold")


def jaccard(left: str, right: str) -> float:
    a, b = set(left.split(";")), set(right.split(";"))
    return len(a & b) / len(a | b)


def main() -> None:
    input_hashes = verify_inputs()
    formal = pd.read_csv(ROOT / SOURCES["S070"], sep="\t")
    multi = pd.read_csv(ROOT / SOURCES["S071"], sep="\t")
    seed = pd.read_csv(ROOT / SOURCES["S072"], sep="\t")
    leave_seed = pd.read_csv(ROOT / SOURCES["S073"], sep="\t")
    leave_donor = pd.read_csv(ROOT / SOURCES["S074"], sep="\t")
    controls = pd.read_csv(ROOT / SOURCES["S075"], sep="\t")
    leading = pd.read_csv(ROOT / SOURCES["S076"], sep="\t")
    mapping = pd.read_csv(ROOT / SOURCES["S077"], sep="\t")
    if (multi.shape[0], seed.shape[0], leave_seed.shape[0], leave_donor.shape[0], controls.shape[0]) != (12, 45, 45, 9, 3):
        raise ValueError("unexpected frozen source dimensions")
    s0 = multi.loc[multi["specification"].eq("S0")].set_index("target_gene").reindex(TARGETS)
    f0 = formal.set_index("target_gene").reindex(TARGETS)
    if not np.allclose(s0["NES"], f0["NES"]) or not np.allclose(s0["FDR_within_Hallmark"], f0["padj"]):
        raise ValueError("S0 formal reconciliation failed")
    if not multi.loc[multi["specification"].eq("S3"), ["p_value", "FDR_within_Hallmark"]].isna().all().all():
        raise ValueError("S3 inferential boundary failed")
    if not (multi["NES"] > 0).all() or not (multi.loc[multi["specification"].isin(["S0", "S1", "S2"]), "FDR_within_Hallmark"] >= 0.05).all():
        raise ValueError("direction or FDR boundary failed")
    source_mapping = mapping.loc[mapping["row_type"].eq("SOURCE_MAPPING")].iloc[0]
    if int(source_mapping["universe_member_count"]) != 191 or int(source_mapping["source_member_count_unique"]) != 200:
        raise ValueError("mapping audit failed")

    spec_pairs = []
    for target, group in multi.loc[multi["specification"].isin(["S0", "S1", "S2"])].groupby("target_gene"):
        rows = group.set_index("specification")
        for left, right in itertools.combinations(["S0", "S1", "S2"], 2):
            spec_pairs.append({"comparison_type": "between_specification_rank_level", "target_gene": target,
                               "comparison_id_a": left, "comparison_id_b": right,
                               "leading_edge_jaccard": jaccard(rows.loc[left, "leading_edge"], rows.loc[right, "leading_edge"]),
                               "leading_edge_size_a": rows.loc[left, "leading_edge_size"],
                               "leading_edge_size_b": rows.loc[right, "leading_edge_size"],
                               "status": "DISPLAY_DERIVED_FROM_FROZEN_S071"})
    spec_pairs = pd.DataFrame(spec_pairs)
    leading_display = pd.concat([leading, spec_pairs], ignore_index=True)
    stability = pd.concat([
        seed.assign(analysis="Seed-level (n=45)"),
        leave_seed.assign(analysis="Leave-one-seed (n=45)"),
        leave_donor.assign(analysis="Leave-one-donor (n=9)"),
    ], ignore_index=True, sort=False)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    for frame, name in [
        (multi, "Figure7a_specification_NES_source_data.tsv"),
        (multi.loc[multi["specification"].isin(["S0", "S1", "S2"])], "Figure7b_inferential_FDR_source_data.tsv"),
        (stability, "Figure7c_direction_stability_source_data.tsv"),
        (controls, "Figure7d_matched_control_source_data.tsv"),
        (leading_display, "Figure7e_leading_edge_stability_source_data.tsv"),
        (mapping.loc[mapping["row_type"].eq("SOURCE_MAPPING")], "Figure7f_mapping_source_data.tsv"),
    ]:
        frame.to_csv(SOURCE_DIR / name, sep="\t", index=False)

    fig = plt.figure(figsize=(7.205, 6.102))
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.1])
    fig.subplots_adjust(left=.095, right=.975, top=.91, bottom=.09, wspace=.66, hspace=.78)
    fig.text(.095, .965, "Exploratory sensitivity analysis of HALLMARK_HYPOXIA", fontsize=9, fontweight="bold")
    fig.text(.095, .942, "Robust positive rank direction, but FDR-negative and nonspecific relative to matched genes", fontsize=5.3, color="#555555")
    a, b, c = [fig.add_subplot(grid[0, i]) for i in range(3)]
    d, e, f = [fig.add_subplot(grid[1, i]) for i in range(3)]

    # a: NES across frozen and exploratory specifications
    offsets = {"TYMS": -.09, "EFNB2": 0, "LRRC17": .09}
    for target in TARGETS:
        q = multi.loc[multi["target_gene"].eq(target)].set_index("specification").reindex(["S0", "S1", "S2", "S3"])
        x = np.arange(4) + offsets[target]
        a.plot(x, q["NES"], color=COLORS[target], lw=.9, alpha=.85)
        a.scatter(x[:3], q["NES"].iloc[:3], s=20, color=COLORS[target], zorder=3, label=target)
        a.scatter(x[3], q["NES"].iloc[3], s=26, facecolor="white", edgecolor=COLORS[target], marker="s", lw=1.1, zorder=3)
    a.axhline(0, color="#777777", lw=.65)
    a.set_xlim(-.35, 3.35); a.set_ylim(0, 1.05); a.set_xticks(range(4)); a.set_xticklabels(["S0\nformal", "S1\nexplor.", "S2\nexplor.", "S3\ndescr."], fontsize=4.6)
    a.set_ylabel("HALLMARK_HYPOXIA NES", fontsize=5.1); a.tick_params(axis="y", labelsize=4.7)
    a.set_title("Direction across specifications", loc="left", fontsize=6.7, fontweight="bold")
    a.legend(fontsize=4.3, loc="lower left", bbox_to_anchor=(0, .09), handlelength=.8)
    a.text(.02, .29, "S3: descriptive only; no P/FDR", transform=a.transAxes, fontsize=3.9, color="#666666")
    panel_label(a, "a")

    # b: only inferential FDR values
    q = multi.loc[multi["specification"].isin(["S0", "S1", "S2"])].copy()
    ybase = np.arange(len(TARGETS))
    for i, spec in enumerate(["S0", "S1", "S2"]):
        values = q.loc[q["specification"].eq(spec)].set_index("target_gene").reindex(TARGETS)["FDR_within_Hallmark"]
        y = ybase + (i - 1) * .18
        b.scatter(values, y, s=19, color=SPEC_COLORS[spec], label=spec, zorder=3)
    b.axvspan(0, .05, color="#F4CCCC", alpha=.65, zorder=0)
    b.axvline(.05, color="#B22222", lw=.7, ls="--")
    b.set_xlim(0, 1.03); b.set_yticks(ybase); b.set_yticklabels(TARGETS, fontsize=5.1); b.invert_yaxis()
    b.set_xlabel("FDR within Hallmark", fontsize=5.1); b.tick_params(axis="x", labelsize=4.7)
    b.set_title("Inferential FDR remains > 0.05", loc="left", fontsize=6.7, fontweight="bold")
    b.legend(fontsize=4.2, ncol=1, loc="upper left", columnspacing=.7, handletextpad=.25)
    b.text(.06, .04, "S3 has no inferential P/FDR", transform=b.transAxes, fontsize=3.8, color="#666666")
    panel_label(b, "b")

    # c: direction stability distributions
    analyses = ["Seed-level (n=45)", "Leave-one-seed (n=45)", "Leave-one-donor (n=9)"]
    rng = np.random.default_rng(20260912)
    for ai, analysis in enumerate(analyses):
        for ti, target in enumerate(TARGETS):
            values = stability.loc[(stability["analysis"].eq(analysis)) & (stability["target_gene"].eq(target)), "NES"].to_numpy()
            pos = ai + (ti - 1) * .22
            c.scatter(pos + rng.uniform(-.045, .045, len(values)), values, s=8, color=COLORS[target], alpha=.52, edgecolor="none")
            c.hlines(np.median(values), pos-.075, pos+.075, color=COLORS[target], lw=1.4)
    c.axhline(0, color="#777777", lw=.65)
    c.set_xlim(-.48, 2.48); c.set_ylim(0, 1.08); c.set_xticks(range(3)); c.set_xticklabels(["Seed", "LOO\nseed", "LOO\ndonor"], fontsize=4.6)
    c.set_ylabel("HALLMARK_HYPOXIA NES", fontsize=5.1); c.tick_params(axis="y", labelsize=4.7)
    c.set_title("Positive direction under resampling", loc="left", fontsize=6.7, fontweight="bold")
    c.text(.02, .04, "points: model resamples; bars: median", transform=c.transAxes, fontsize=3.8, color="#666666")
    panel_label(c, "c")

    # d: matched-gene calibration
    calibration = controls.set_index("target_gene").reindex(TARGETS)
    yd = np.arange(3)
    d.hlines(yd, calibration["signed_NES_empirical_percentile"], calibration["conservative_absolute_upper_tail"], color="#B0B0B0", lw=1.1)
    d.scatter(calibration["signed_NES_empirical_percentile"], yd, s=29, c=[COLORS[t] for t in TARGETS], label="signed percentile", zorder=3)
    d.scatter(calibration["conservative_absolute_upper_tail"], yd, s=29, facecolors="white", edgecolors="#333333", lw=.9, label="conservative abs. upper-tail", zorder=3)
    for y, (_, row) in enumerate(calibration.iterrows()):
        d.text(min(row["conservative_absolute_upper_tail"] + .03, .95), y, f"{row['conservative_absolute_upper_tail']:.2f}", va="center", fontsize=4.1)
    d.set_xlim(0, 1.02); d.set_yticks(yd); d.set_yticklabels(TARGETS, fontsize=5.1); d.invert_yaxis()
    d.set_xlabel("Empirical position vs 10 matched genes", fontsize=5.0); d.tick_params(axis="x", labelsize=4.7)
    d.set_title("No target-specific calibration signal", loc="left", fontsize=6.7, fontweight="bold")
    d.legend(fontsize=3.9, loc="upper left", handletextpad=.3)
    d.text(.02, .04, "positions, not P values", transform=d.transAxes, fontsize=3.8, color="#666666")
    panel_label(d, "d")

    # e: leading-edge stability, including deterministic S0-S2 display summary
    labels = ["Within\ndonor-seed", "Between\ndonor", "Between\ntarget", "S0–S2\nspecs"]
    groups = ["within_target_donor_seed", "between_donor_five_seed_consensus", "between_target_formal_consensus", "between_specification_rank_level"]
    rng = np.random.default_rng(20260913)
    data = [leading_display.loc[leading_display["comparison_type"].eq(group), "leading_edge_jaccard"].to_numpy() for group in groups]
    e.boxplot(data, positions=np.arange(4), widths=.48, patch_artist=True,
              boxprops={"facecolor": "#E5E5E5", "edgecolor": "#555555", "linewidth": .7},
              medianprops={"color": "#222222", "linewidth": 1}, whiskerprops={"linewidth": .7}, capprops={"linewidth": .7},
              flierprops={"marker": "", "markersize": 0})
    for i, values in enumerate(data):
        e.scatter(i + rng.uniform(-.12, .12, len(values)), values, s=7, color="#4C78A8", alpha=.55, edgecolor="none")
    e.set_ylim(.45, 1.02); e.set_xticks(range(4)); e.set_xticklabels(labels, fontsize=4.4); e.set_ylabel("Leading-edge Jaccard", fontsize=5.0); e.tick_params(axis="y", labelsize=4.7)
    e.set_title("Leading-edge overlap is stable", loc="left", fontsize=6.7, fontweight="bold")
    e.text(.02, .04, "S0–S2 only; S3 is a descriptive union", transform=e.transAxes, fontsize=3.6, color="#666666")
    panel_label(e, "e")

    # f: mapping adequacy and exact unavailable list
    mapped, unavailable = int(source_mapping["universe_member_count"]), int(source_mapping["source_member_count_unique"]) - int(source_mapping["universe_member_count"])
    f.barh([0], [mapped], color="#5B8DB8", height=.38, label="exactly mapped")
    f.barh([0], [unavailable], left=[mapped], color="#D99058", height=.38, label="unavailable")
    f.text(mapped / 2, 0, f"{mapped} mapped", ha="center", va="center", fontsize=5.1, color="white", fontweight="bold")
    f.text(mapped + unavailable / 2, 0, f"{unavailable}", ha="center", va="center", fontsize=5.1, color="white", fontweight="bold")
    f.set_xlim(0, 200); f.set_ylim(-.85, .42); f.set_yticks([]); f.set_xlabel("Original HALLMARK_HYPOXIA members", fontsize=5.0); f.tick_params(axis="x", labelsize=4.7)
    f.set_title("Exact-symbol mapping quality", loc="left", fontsize=6.7, fontweight="bold")
    f.legend(fontsize=3.9, loc="upper left", bbox_to_anchor=(0, .88), ncol=2)
    missing = source_mapping["missing_members"].replace(";", ", ")
    f.text(0, -.31, "191/200 (95.5%) exact mapping\nUnavailable: " + missing, transform=f.transData, fontsize=3.7, va="top", color="#555555", wrap=True)
    panel_label(f, "f")

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "Figure7_exploratory_hypoxia_sensitivity"
    for ext, kwargs in [(".svg", {}), (".pdf", {}), (".png", {"dpi": 300}), (".tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
        fig.savefig(base.with_suffix(ext), facecolor="white", **kwargs)
    plt.close(fig)
    audit = {
        "status": "PASS", "figure": "Figure 7", "backend": "Python/matplotlib", "input_manifest_validation": "PASS",
        "input_sha256": input_hashes,
        "assertions": {"all_12_NES_positive": True, "inferential_FDR_count": 9, "inferential_FDR_all_ge_0_05": True,
                       "seed_level_rows": 45, "leave_one_seed_rows": 45, "leave_one_donor_rows": 9,
                       "matched_control_rows": 3, "leading_edge_rows": int(leading_display.shape[0]),
                       "mapping_exact": f"{mapped}/200", "S3_no_P_or_FDR": True},
        "output_dimensions_mm": [183, 155],
        "interpretation_boundary": "Exploratory sensitivity analysis only: positive NES is rank direction, not expression upregulation; robust direction does not override formal FDR-negative and nonspecific matched-control results.",
    }
    (OUT / "Figure7_build_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "output": str(OUT)}))


if __name__ == "__main__":
    main()
