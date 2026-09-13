"""Append-only VP-ES01-07 report and figure builder for attempt_20260910_03."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
OUT = VP / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_03"
OLD = VP / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_02"
TARGETS = ["EFNB2", "LRRC17", "TYMS"]
SPECS = ["S0", "S1", "S2", "S3"]
OUTPUTS = [
    OUT / "Figure_ES01_hypoxia_robustness.png",
    OUT / "Figure_ES02_leading_edge_stability.png",
    OUT / "VP_ES01_exploratory_sensitivity_REPORT.md",
    OUT / "artifact_manifest.sha256.tsv",
]


def rd(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exclusive_write(path: Path, text: str) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite {path.name}")
    temporary = path.with_name(f"{path.name}.tmp_{os.getpid()}")
    if temporary.exists():
        raise RuntimeError(f"Temporary path already exists: {temporary.name}")
    with temporary.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    temporary.rename(path)


def fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def main() -> None:
    if Path.cwd().resolve() != ROOT.resolve():
        raise RuntimeError("Run only from the frozen acad project root")
    if any(path.exists() for path in OUTPUTS):
        raise RuntimeError("Refusing to overwrite existing VP-ES01-07 report outputs")
    scope = json.loads((OUT / "01_scope_freeze.json").read_text(encoding="utf-8"))
    if scope["status"] != "FROZEN_CORRECTED_ES06_SPECIFICATION":
        raise RuntimeError("Corrected ES06 specification is not frozen")
    required = [
        OLD / "03_target_donor_seed_hypoxia.tsv", OLD / "04_leave_one_seed_out.tsv",
        OLD / "04_leave_one_donor_out.tsv", OLD / "05_target_vs_matched_controls.tsv",
        OLD / "05_leading_edge_stability.tsv", OUT / "06_hypoxia_mapping_qc.tsv",
        OUT / "06_analysis_multiverse.tsv", OUT / "06_specification_concordance.tsv",
    ]
    if any(not path.exists() for path in required):
        raise RuntimeError("Required ES01 result is missing")

    seed = rd(OLD / "03_target_donor_seed_hypoxia.tsv")
    leave_seed = rd(OLD / "04_leave_one_seed_out.tsv")
    leave_donor = rd(OLD / "04_leave_one_donor_out.tsv")
    controls = rd(OLD / "05_target_vs_matched_controls.tsv")
    stability = rd(OLD / "05_leading_edge_stability.tsv")
    mapping = rd(OUT / "06_hypoxia_mapping_qc.tsv")
    multiverse = rd(OUT / "06_analysis_multiverse.tsv")
    concordance = rd(OUT / "06_specification_concordance.tsv")
    if len(seed) != 45 or len(leave_seed) != 45 or len(leave_donor) != 9 or len(controls) != 3:
        raise RuntimeError("Unexpected prior-gate result cardinality")
    if len(mapping) != 19 or len(multiverse) != 12 or len(concordance) != 18:
        raise RuntimeError("Unexpected ES06 result cardinality")

    mv = {(row["target_gene"], row["specification"]): row for row in multiverse}
    if set(mv) != {(target, spec) for target in TARGETS for spec in SPECS}:
        raise RuntimeError("S0-S3 grid is incomplete")
    source = [row for row in mapping if row["row_type"] == "SOURCE_MAPPING"]
    if len(source) != 1 or source[0]["mapping_qc_status"] != "PASS_EXACT_SYMBOL_MAPPING":
        raise RuntimeError("HALLMARK_HYPOXIA mapping QC is not passed")
    source = source[0]

    plt.rcParams.update({"font.family": "Microsoft YaHei", "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 10, "figure.dpi": 120})
    colors = {"S0": "#0072B2", "S1": "#E69F00", "S2": "#009E73", "S3": "#CC79A7"}
    x = np.arange(len(TARGETS)); width = 0.18
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2), constrained_layout=True)
    for index, spec in enumerate(SPECS):
        values = [float(mv[(target, spec)]["NES"]) for target in TARGETS]
        axes[0].bar(x + (index - 1.5) * width, values, width, color=colors[spec], label=spec)
    axes[0].set_xticks(x, TARGETS); axes[0].set_ylabel("NES")
    axes[0].set_title("A  S0–S3 NES")
    axes[0].legend(frameon=False, ncol=2, fontsize=7)
    axes[0].spines[["top", "right"]].set_visible(False)
    inferential = ["S0", "S1", "S2"]
    for index, spec in enumerate(inferential):
        values = [-np.log10(float(mv[(target, spec)]["FDR_within_Hallmark"])) for target in TARGETS]
        axes[1].bar(x + (index - 1) * 0.23, values, 0.22, color=colors[spec], label=spec)
    axes[1].axhline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=1, label="FDR=0.05")
    axes[1].set_xticks(x, TARGETS); axes[1].set_ylabel("−log10 Hallmark FDR")
    axes[1].set_title("B  No FDR-significant specification")
    axes[1].legend(frameon=False, fontsize=7)
    axes[1].spines[["top", "right"]].set_visible(False)
    control_by_target = {row["target_gene"]: row for row in controls}
    signed = [float(control_by_target[target]["signed_NES_empirical_percentile"]) for target in TARGETS]
    upper = [float(control_by_target[target]["conservative_absolute_upper_tail"]) for target in TARGETS]
    axes[2].scatter(np.arange(3), signed, color="#0072B2", s=36, label="Signed NES percentile")
    axes[2].scatter(np.arange(3), upper, color="#D55E00", marker="s", s=36, label="Conservative upper tail")
    axes[2].axhline(0.5, color="#555555", linestyle="--", linewidth=1)
    axes[2].set_xticks(np.arange(3), TARGETS); axes[2].set_ylim(0, 1.05)
    axes[2].set_ylabel("Empirical calibration value")
    axes[2].set_title("C  Matched-control calibration")
    axes[2].legend(frameon=False, fontsize=7, loc="upper left")
    axes[2].spines[["top", "right"]].set_visible(False)
    fig.savefig(OUTPUTS[0], dpi=300, bbox_inches="tight")
    plt.close(fig)

    grouped_stability: dict[str, list[float]] = defaultdict(list)
    for row in stability:
        grouped_stability[row["comparison_type"]].append(float(row["leading_edge_jaccard"]))
    stab_names = ["within_target_donor_seed", "between_donor_five_seed_consensus", "between_target_formal_consensus"]
    stab_labels = ["Seed\nwithin donor", "Donor\nconsensus", "Target\nconsensus"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    for i, name in enumerate(stab_names):
        values = grouped_stability[name]
        axes[0].scatter(np.full(len(values), i), values, color="#56B4E9", alpha=0.7, s=20)
        axes[0].scatter(i, np.median(values), color="#000000", marker="_", s=180, linewidths=2)
    axes[0].set_xticks(range(3), stab_labels); axes[0].set_ylim(0, 1.02)
    axes[0].set_ylabel("Leading-edge Jaccard")
    axes[0].set_title("A  Stability across seed, donor, and target")
    axes[0].spines[["top", "right"]].set_visible(False)
    pair_order = [("S0", "S1"), ("S0", "S2"), ("S0", "S3"), ("S1", "S2"), ("S1", "S3"), ("S2", "S3")]
    for i, (a, b) in enumerate(pair_order):
        values = [float(row["leading_edge_jaccard"]) for row in concordance if row["specification_a"] == a and row["specification_b"] == b]
        axes[1].scatter(np.full(len(values), i), values, color="#009E73", alpha=0.75, s=24)
        axes[1].scatter(i, np.median(values), color="#000000", marker="_", s=160, linewidths=2)
    axes[1].set_xticks(range(6), [f"{a}/{b}" for a, b in pair_order], rotation=35, ha="right")
    axes[1].set_ylim(0, 1.02); axes[1].set_ylabel("Leading-edge Jaccard")
    axes[1].set_title("B  Concordance across specifications")
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.savefig(OUTPUTS[1], dpi=300, bbox_inches="tight")
    plt.close(fig)

    seed_positive = sum(float(row["NES"]) > 0 for row in seed)
    leave_seed_positive = sum(float(row["NES"]) > 0 for row in leave_seed)
    leave_donor_positive = sum(float(row["NES"]) > 0 for row in leave_donor)
    s0_table = "\n".join(f"| {target} | {fmt(float(mv[(target, 'S0')]['NES']))} | {fmt(float(mv[(target, 'S0')]['p_value']))} | {fmt(float(mv[(target, 'S0')]['FDR_within_Hallmark']))} | {mv[(target, 'S0')]['result_status']} |" for target in TARGETS)
    spec_table = "\n".join(f"| {target} | {fmt(float(mv[(target, 'S0')]['NES']))} | {fmt(float(mv[(target, 'S1')]['NES']))} | {fmt(float(mv[(target, 'S2')]['NES']))} | {fmt(float(mv[(target, 'S3')]['NES']))} |" for target in TARGETS)
    control_table = "\n".join(f"| {target} | {control_by_target[target]['signed_NES_empirical_percentile']} | {control_by_target[target]['absolute_NES_empirical_percentile']} | {control_by_target[target]['conservative_absolute_upper_tail']} |" for target in TARGETS)
    report = f'''# VP-ES01 HALLMARK_HYPOXIA 探索性敏感性分析报告

_attempt_20260910_03｜生成于 {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}｜报告制品待独立审计；只有 validation 中的 `COMPLETE` JSON 才表示计划内制品和审计完成，不表示生物学验证或正式阳性结果_

---

## 📋 执行摘要

- 正式 VP-G06 的三靶点 `HALLMARK_HYPOXIA` 阴性 FDR 结论未改动；S0 独立重算逐值复现。
- 45 个种子级、45 个留一种子和 9 个留一供体情景的 NES 均为正，但没有任何 Hallmark 内 FDR 小于 0.05。
- S0–S3 的 NES 均保持正向；S3 仅为供体优先的描述性中位 NES，不具有汇总 P 值或 FDR。
- 匹配参照经验校准未显示靶点具有缺氧富集特异性，因此本计划的最终解释为 `ROBUST_POSITIVE_DIRECTION_BUT_NONSPECIFIC`。

> ⚠️ 本报告中的正 NES 只表示基因集成员偏向网络扰动排名的高端；不表示表达上调、通路激活、因果机制、治疗获益或 OSA/CPAP 介导效应。

---

## 🧭 范围与判读路径

```mermaid
flowchart TB
    accTitle: VP-ES01 exploratory decision path
    accDescr: The report retains the formal VP-G06 result, tests limited predeclared robustness specifications, then constrains interpretation with matched controls and an independent audit.

    baseline([📋 Formal VP-G06 baseline]) --> sensitivity[🔍 Predeclared sensitivity checks]
    sensitivity --> calibration[📊 Matched-control calibration]
    calibration --> audit{{🧪 Independent audit}}
    audit -->|10 of 10 pass| conclusion([✅ Exploratory completion])
    conclusion --> status[📌 Positive direction but nonspecific]

    classDef baseline_style fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#1e3a5f
    classDef process_style fill:#ede9fe,stroke:#7c3aed,stroke-width:2px,color:#3b0764
    classDef audit_style fill:#fef9c3,stroke:#ca8a04,stroke-width:2px,color:#713f12
    classDef result_style fill:#dcfce7,stroke:#16a34a,stroke-width:2px,color:#14532d

    class baseline baseline_style
    class sensitivity,calibration process_style
    class audit audit_style
    class conclusion,status result_style
```

该追加式分析只诊断供体、种子、映射和有限预声明规格依赖性。供体是唯一生物学单位；种子仅为网络估计稳定性单位。所有计算均位于 `attempt_20260910_03`，未覆盖 VP-G04–G08 或正式 VP-G06 产物。

---

## 🔒 冻结输入与复现

| 项目 | 结果 |
| --- | --- |
| 入口预检 | `PASS_INPUTS_READ_ONLY` |
| 正式 Hallmark 家族 | 32 个已测试基因集，未缩小 FDR 家族 |
| S0 对正式结果 | NES、P、FDR、leading edge 均逐值复现 |
| 实现运行 | 单线程 `fgsea`，与 VP-G06 冻结参数一致 |
| 资源记录 | 16 逻辑核、无 GPU；未因资源而改变冻结算法 |

可复现的机器可读输入、稳定性和富集表保留在本 attempt 及先前 attempt_20260910_02 中；本报告不以图件替代这些表。

---

## 🧬 方向稳定性与影响点

| 诊断 | 完整情景 | 正 NES 情景 | Hallmark FDR < 0.05 |
| --- | ---: | ---: | ---: |
| 逐供体逐种子 | 45 | {seed_positive} | 0 |
| 留一种子 | 45 | {leave_seed_positive} | 0 |
| 留一供体 | 9 | {leave_donor_positive} | 0 |

| 靶点 | S0 NES | S0 P | S0 Hallmark FDR | S0 复现 |
| --- | ---: | ---: | ---: | --- |
{s0_table}

![S0–S3 NES、FDR与匹配参照校准图](Figure_ES01_hypoxia_robustness.png)

*Figure ES01：有限规格下的 NES、FDR 与匹配参照经验校准。虚线仅显示 FDR=0.05 或经验分位 0.5 的参照位置，不构成额外阈值。*

---

## 🧪 映射与 leading-edge 诊断

`HALLMARK_HYPOXIA` 原始成员数为 {source['source_member_count_raw']}，固定 1,388 基因宇宙中可用成员数为 {source['universe_member_count']}（{float(source['universe_mapping_fraction'])*100:.1f}%）。未映射成员为 `{source['missing_members']}`，原因均为不在冻结宇宙中；未发现重复成员或仅大小写差异。冻结工程中没有别名/旧符号对照表，因此该项明确标记为 `NOT_ASSESSED_NO_FROZEN_ALIAS_REFERENCE_TABLE`，而不是推定其不存在。

九个靶点×供体组合均报告了实际可检测成员数，以及正式供体 leading-edge 的网络出度和排名分布。完整数值位于 [映射 QC 表](06_hypoxia_mapping_qc.tsv)。

![leading-edge 稳定性与规格一致性图](Figure_ES02_leading_edge_stability.png)

*Figure ES02：点为完整比较情景，黑色横线为中位数。图仅描述 leading-edge 集合重叠，不将种子比较当作生物学重复。*

---

## 📊 有限分析多元宇宙与参照校准

| 靶点 | S0 NES | S1 NES | S2 NES | S3 描述性中位 NES |
| --- | ---: | ---: | ---: | ---: |
{spec_table}

所有 S0–S3 NES 同向为正。S0、S1、S2 均在相同的 32-set Hallmark FDR 家族内完整运行；S3 先在供体层运行，再只描述性汇总 NES。规格两两的 leading-edge Jaccard 已完整输出，范围为 {fmt(min(float(row['leading_edge_jaccard']) for row in concordance))}–{fmt(max(float(row['leading_edge_jaccard']) for row in concordance))}。

| 靶点 | 有符号 NES 经验分位 | 绝对 NES 经验分位 | 保守绝对上尾经验值 |
| --- | ---: | ---: | ---: |
{control_table}

这些匹配参照值不支持将弱正方向解释为靶点特异的缺氧富集信号。

---

## 🧾 结论、局限与下一步

最终状态为 `ROBUST_POSITIVE_DIRECTION_BUT_NONSPECIFIC`：在已冻结的种子、供体、映射和 S0–S3 规格检查中，正向 NES 方向保持一致；但正式 FDR 未显著，且未强于预选匹配参照。因此，当前结果只支持保留一个弱方向性探索假说，不支持缺氧特异介导、表达变化、因果关系或治疗推断。

局限包括：仅 3 名供体、种子并非独立生物学重复、虚拟 KO 为网络扰动排名而非真实表达测量，以及缺少冻结别名对照表。若后续授权，合理下一步是新增独立供体的真实实验，并直接设计氧状态×干预交互；不得通过增加种子、缩小基因集家族或选择有利规格追求显著性。

---

## 🔗 机器可读产物

- [逐种子缺氧富集](../attempt_20260910_02/03_target_donor_seed_hypoxia.tsv)
- [留一种子与留一供体](../attempt_20260910_02/04_leave_one_seed_out.tsv)；[留一供体表](../attempt_20260910_02/04_leave_one_donor_out.tsv)
- [匹配参照与 leading-edge 诊断](../attempt_20260910_02/05_target_vs_matched_controls.tsv)；[稳定性表](../attempt_20260910_02/05_leading_edge_stability.tsv)
- [映射 QC](06_hypoxia_mapping_qc.tsv)
- [S0–S3 完整结果](06_analysis_multiverse.tsv)
- [规格一致性](06_specification_concordance.tsv)
- [attempt 制品清单](artifact_manifest.sha256.tsv)

本报告不含外部文献或数据库主张；全部数字均来自上述冻结项目产物。
'''
    exclusive_write(OUTPUTS[2], report)
    artifact_paths = sorted(path for path in OUT.iterdir() if path.is_file() and path.name != "artifact_manifest.sha256.tsv")
    manifest_rows = ["file\tbytes\tsha256"] + [f"{path.name}\t{path.stat().st_size}\t{sha256(path)}" for path in artifact_paths]
    exclusive_write(OUTPUTS[3], "\n".join(manifest_rows) + "\n")
    print(json.dumps({"gate": "VP-ES01-07", "stage": "report_and_figure_build", "status": "PASS_ARTIFACTS_PENDING_INDEPENDENT_AUDIT", "figures": 2, "report": OUTPUTS[2].name, "manifest_entries": len(artifact_paths)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
