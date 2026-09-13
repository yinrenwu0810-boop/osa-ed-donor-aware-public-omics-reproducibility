#!/usr/bin/env python3
"""Create a bounded Figure 5 from sealed VP-G05--G07 outputs."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


VP = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
G05 = VP / "validation" / "VP_G05_stability_COMPLETE.json"
G06 = VP / "validation" / "VP_G06_enrichment_COMPLETE.json"
G07 = VP / "validation" / "VP_G07_external_scope_COMPLETE.json"


def label_panel(ax, letter: str) -> None:
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def box(ax, x, y, width, height, title, body, color):
    patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.02,rounding_size=0.025", facecolor=color, edgecolor="#4B5563", linewidth=0.8)
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height * 0.66, title, ha="center", va="center", fontsize=8.4, fontweight="bold")
    ax.text(x + width / 2, y + height * 0.32, body, ha="center", va="center", fontsize=7.1, linespacing=1.2)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if any((OUT / f"Figure5_virtual_perturbation_hypothesis_prioritization.{suffix}").exists() for suffix in ("png", "pdf", "tiff")):
        raise SystemExit("Refusing to overwrite an existing Figure 5 export")
    g05 = json.loads(G05.read_text(encoding="utf-8"))
    g06 = json.loads(G06.read_text(encoding="utf-8"))
    g07 = json.loads(G07.read_text(encoding="utf-8"))
    scores = {row["target_gene"]: row for row in g06["ED_signature_summary"]["ED_UP_FDR_ranked_GSEA"]}
    targets = ["EFNB2", "LRRC17", "TYMS"]
    nes = [scores[target]["NES"] for target in targets]
    fdr = [scores[target]["FDR_within_target"] for target in targets]
    target_background = ["nonDM-ED\n3 donors", "nonDM-ED\n3 donors", "normal\n3 donors"]

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Microsoft YaHei", "Arial"], "font.size": 8, "axes.unicode_minus": False})
    fig = plt.figure(figsize=(7.1, 6.4), layout="constrained")
    grid = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.32)

    ax_a = fig.add_subplot(grid[0, :])
    ax_a.set_axis_off(); ax_a.set_xlim(0, 1); ax_a.set_ylim(0, 1)
    label_panel(ax_a, "a")
    steps = [
        ("G04", "45 个预声明\n网络虚拟KO", "#DCE6F1"),
        ("G05", "3/3 供体、5 种子\n稳定性通过", "#D9EAD3"),
        ("G06", "完整有序网络排名\n与ED签名关联", "#FFF2CC"),
        ("G07", "无精确T1外部\n真实扰动记录", "#FCE4D6"),
    ]
    for idx, (name, detail, color) in enumerate(steps):
        x = 0.02 + idx * 0.245
        box(ax_a, x, 0.31, 0.205, 0.39, name, detail, color)
        if idx < len(steps) - 1:
            ax_a.annotate("", xy=(x + 0.238, 0.505), xytext=(x + 0.210, 0.505), arrowprops={"arrowstyle": "->", "lw": 1, "color": "#6B7280"})
    ax_a.text(0.5, 0.11, "所有步骤为计算/审计证据；G08 为后续可行性设计，尚无细胞实验数据。", ha="center", va="center", fontsize=8.2, color="#374151")

    ax_b = fig.add_subplot(grid[1, 0])
    label_panel(ax_b, "b")
    colors = ["#0072B2", "#009E73", "#E69F00"]
    bars = ax_b.bar(targets, nes, color=colors, width=0.62, edgecolor="#374151", linewidth=0.5)
    ax_b.set_ylim(0, 2.15); ax_b.set_ylabel("网络排名 GSEA NES")
    ax_b.set_title("与 ED-FDR 上调基因集的关联", fontsize=9, pad=7)
    ax_b.spines[["top", "right"]].set_visible(False)
    ax_b.set_axisbelow(True); ax_b.yaxis.grid(True, color="#E5E7EB", linewidth=0.6)
    for bar, value, q in zip(bars, nes, fdr):
        ax_b.text(bar.get_x() + bar.get_width() / 2, value + 0.06, f"FDR={q:.2g}", ha="center", va="bottom", fontsize=7)
    ax_b.text(0.5, -0.28, "NES 反映有序网络排名关联，\n不是 KO 后表达方向或效应量。", transform=ax_b.transAxes, ha="center", va="top", fontsize=7.1, color="#4B5563")

    ax_c = fig.add_subplot(grid[1, 1])
    label_panel(ax_c, "c")
    ax_c.set_axis_off(); ax_c.set_xlim(0, 1); ax_c.set_ylim(0, 1)
    ax_c.text(0.02, 0.98, "跨供体稳定性与外部边界", ha="left", va="top", fontsize=9, fontweight="bold")
    for row, target, background, color in zip([0.70, 0.48, 0.26], targets, target_background, colors):
        ax_c.scatter(0.08, row, s=62, marker="o", color=color, edgecolors="#374151", linewidths=0.5)
        ax_c.text(0.16, row + 0.045, target, ha="left", va="center", fontsize=8.6, fontweight="bold")
        ax_c.text(0.16, row - 0.045, f"{background}；5 seeds；PASS", ha="left", va="center", fontsize=7.1)
        ax_c.text(0.96, row, "3/3", ha="right", va="center", fontsize=8.2, fontweight="bold", color="#166534")
    ax_c.plot([0.03, 0.97], [0.14, 0.14], color="#D1D5DB", linewidth=0.8)
    ax_c.text(0.03, 0.08, f"精确T1外部参照：{g07['exact_T1_external_reference']}", ha="left", va="center", fontsize=7.0, color="#9A3412")

    for suffix, dpi in (("png", 600), ("tiff", 600), ("pdf", None)):
        target = OUT / f"Figure5_virtual_perturbation_hypothesis_prioritization.{suffix}"
        fig.savefig(target, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
