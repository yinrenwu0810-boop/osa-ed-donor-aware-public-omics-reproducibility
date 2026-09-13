#!/usr/bin/env python3
"""Apply the frozen strict Hallmark bridge rule for endothelial IH, OSA, and CPAP."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"
PROTOCOL = REVISION / "00_protocol"
RESULTS = REVISION / "02_results"
VALIDATION = REVISION / "04_validation"
SOURCE_RELATIVE = "analysis/results/integration/Hallmark_cross_dataset_evidence.tsv"


def write_numbered_attempt(prefix: str, status: str, payload: str) -> Path:
    number = 1
    while True:
        if any(VALIDATION.glob(f"{prefix}_{number:02d}_*.json")):
            number += 1
            continue
        target = VALIDATION / f"{prefix}_{number:02d}_{status}.json"
        try:
            with target.open("x", encoding="utf-8", newline="\n") as handle:
                handle.write(payload)
            return target
        except FileExistsError:
            number += 1


def atomic_replace_json(path: Path, payload: str) -> None:
    counter = 0
    while True:
        temporary = path.with_name(f".{path.name}.{os.getpid()}.{counter}.tmp")
        try:
            handle = temporary.open("x", encoding="utf-8", newline="\n")
            break
        except FileExistsError:
            counter += 1
    try:
        with handle:
            handle.write(payload)
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source = ROOT / SOURCE_RELATIVE
    if not source.is_file():
        raise FileNotFoundError(source)

    input_manifest = pd.DataFrame(
        [
            {
                "relative_path": SOURCE_RELATIVE,
                "bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
        ]
    )
    input_manifest.to_csv(
        PROTOCOL / "gate_02_input_manifest.sha256.tsv", sep="\t", index=False
    )

    table = pd.read_csv(source, sep="\t")
    required = [
        "pathway",
        "endo243_NES",
        "endo243_FDR",
        "endo205_NES",
        "endo205_FDR",
        "osaAM_NES",
        "osaAM_FDR",
        "cpapAM_NES",
        "cpapAM_FDR",
    ]
    missing = [column for column in required if column not in table]
    if missing:
        raise RuntimeError(f"Missing required columns: {missing}")

    table["endo243_FDR_pass"] = table["endo243_FDR"].lt(0.05)
    table["endo205_FDR_pass"] = table["endo205_FDR"].lt(0.05)
    table["osaAM_FDR_pass"] = table["osaAM_FDR"].lt(0.05)
    table["cpapAM_FDR_pass"] = table["cpapAM_FDR"].lt(0.05)
    estimate_columns = [
        "endo243_NES",
        "endo243_FDR",
        "endo205_NES",
        "endo205_FDR",
        "osaAM_NES",
        "osaAM_FDR",
        "cpapAM_NES",
        "cpapAM_FDR",
    ]
    table["all_four_estimated"] = table[estimate_columns].notna().all(axis=1)
    table["all_four_NES_nonzero"] = table[
        ["endo243_NES", "endo205_NES", "osaAM_NES", "cpapAM_NES"]
    ].ne(0).all(axis=1)
    table["all_four_FDR_pass"] = table[
        ["endo243_FDR_pass", "endo205_FDR_pass", "osaAM_FDR_pass", "cpapAM_FDR_pass"]
    ].all(axis=1)
    table["endothelial_direction_concordant_v2"] = (
        np.sign(table["endo243_NES"]) == np.sign(table["endo205_NES"])
    )
    table["OSA_matches_endothelial_direction"] = (
        np.sign(table["osaAM_NES"]) == np.sign(table["endothelial_NES_mean"])
    )
    table["CPAP_reverses_OSA_direction"] = (
        np.sign(table["cpapAM_NES"]) == -np.sign(table["osaAM_NES"])
    )
    table["complete_directional_pattern"] = table[
        [
            "endothelial_direction_concordant_v2",
            "OSA_matches_endothelial_direction",
            "CPAP_reverses_OSA_direction",
        ]
    ].all(axis=1) & table["all_four_estimated"] & table["all_four_NES_nonzero"]
    table["strict_clinical_bridge_pass"] = (
        table["all_four_estimated"]
        & table["all_four_FDR_pass"]
        & table["complete_directional_pattern"]
    )
    table["revised_bridge_class"] = np.select(
        [
            table["strict_clinical_bridge_pass"],
            table["complete_directional_pattern"],
        ],
        ["STRICT_BRIDGE_PASS", "DIRECTIONAL_PATTERN_FDR_INCOMPLETE"],
        default="OTHER",
    )

    first = [
        "pathway",
        "revised_bridge_class",
        "strict_clinical_bridge_pass",
        "complete_directional_pattern",
        "all_four_FDR_pass",
        "all_four_estimated",
        "all_four_NES_nonzero",
        "endothelial_direction_concordant_v2",
        "OSA_matches_endothelial_direction",
        "CPAP_reverses_OSA_direction",
    ]
    remaining = [column for column in table.columns if column not in first]
    table = table[first + remaining].sort_values(
        ["strict_clinical_bridge_pass", "complete_directional_pattern", "endothelial_max_FDR", "pathway"],
        ascending=[False, False, True, True],
    )
    table.to_csv(RESULTS / "Hallmark_clinical_bridge_v2.tsv", sep="\t", index=False)

    strict_pathways = table.loc[table["strict_clinical_bridge_pass"], "pathway"].tolist()
    summary = pd.DataFrame(
        [
            ("hallmark_pathways_total", len(table)),
            ("all_four_estimated", int(table["all_four_estimated"].sum())),
            ("complete_directional_pattern", int(table["complete_directional_pattern"].sum())),
            ("all_four_FDR_pass", int(table["all_four_FDR_pass"].sum())),
            ("strict_clinical_bridge_pass", int(table["strict_clinical_bridge_pass"].sum())),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(RESULTS / "Hallmark_clinical_bridge_summary_v2.tsv", sep="\t", index=False)

    checks = {
        "source_rows_50": len(table) == 50,
        "all_four_estimated_49": int(table["all_four_estimated"].sum()) == 49,
        "strict_pass_count_1": len(strict_pathways) == 1,
        "strict_pass_is_HALLMARK_HYPOXIA": strict_pathways == ["HALLMARK_HYPOXIA"],
        "strict_pathway_all_four_FDR_pass": bool(
            table.loc[table["strict_clinical_bridge_pass"], "all_four_FDR_pass"].all()
        ),
        "strict_pathway_complete_directional_pattern": bool(
            table.loc[
                table["strict_clinical_bridge_pass"], "complete_directional_pattern"
            ].all()
        ),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "strict_pathways": strict_pathways,
        "rule": (
            "endo243, endo205, OSA-AM, and CPAP-AM each FDR < 0.05; "
            "the two endothelial NES directions agree; OSA-AM matches the mean endothelial direction; "
            "CPAP-AM reverses the OSA-AM direction"
        ),
        "evidence_boundary": (
            "This is a result-known audit sensitivity rule applied to an existing observational/longitudinal bridge. "
            "It is supportive and does not establish clinical mediation or treatment efficacy."
        ),
    }
    payload = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    write_numbered_attempt(
        "gate_02_validation_attempt", report["status"], payload
    )
    atomic_replace_json(
        VALIDATION / "gate_02_validation_latest.json", payload
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
