#!/usr/bin/env python3
"""Freeze VP-G01 inputs and the result-blind virtual-KO run design.

This script reads only the seven declared authoritative inputs.  It never
opens the H5AD, creates a matrix, installs a package, or runs a model.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_INPUTS = (
    {
        "relative_path": "revision_v2/01_work/GSE206528_replay/attempt_07/GSE206528_annotated_v2.h5ad",
        "bytes": 1228957045,
        "sha256": "4390e11cda2732bd171ee259cf6f61cdb60373117248edc9864e1ace1dadb293",
    },
    {
        "relative_path": "revision_v2/01_work/GSE206528_doublet_sensitivity/attempt_01/doublet_calls_by_cell.tsv",
        "bytes": 7941574,
        "sha256": "5145e59877a32df9d49fef3d7e1dcd9492dbbaf4182a1cf9389570f45b3c04b8",
    },
    {
        "relative_path": "revision_v2/01_work/GSE206528_doublet_sensitivity/attempt_02/edgeR_Hallmark/TYMS_EFNB2_LRRC17_sensitivity.tsv",
        "bytes": 5030,
        "sha256": "34b149ae5343412e88bb9e6f4c383f11f610450821c8aa423a183be9e291d976",
    },
    {
        "relative_path": "revision_v2/02_results/IH_ED_candidate_classification_v2.tsv",
        "bytes": 107737,
        "sha256": "a43aea8f977bd2cd084793b628d2ed5a9c3d2ff329a8c442c7140865ba0648fe",
    },
    {
        "relative_path": "analysis/data/processed/gene_sets/Hallmark.gmt",
        "bytes": 45517,
        "sha256": "569848c1ac4530b108cbec5b8d1e41a583066b09fb0fdbc61188493271b4c3b7",
    },
    {
        "relative_path": "analysis/data/processed/gene_sets/Reactome.gmt",
        "bytes": 746995,
        "sha256": "8d086721afbb33e9ff69eb5a465be841ad30b0ed0074c5eb682ace0e316333cc",
    },
    {
        "relative_path": "analysis/data/processed/gene_sets/GO_BP.gmt",
        "bytes": 4266646,
        "sha256": "a800139aa6798de4b6eac2b681a921331416336c0bdf598896a467c8f723ec2b",
    },
)

SEEDS = (2026082701, 2026082702, 2026082703, 2026082704, 2026082705)
MODEL_PARAMETERS = {
    "package": "scTenifoldKnk",
    "package_version_required": "1.0.3",
    "qc": False,
    "nc_lambda": 0,
    "nc_nNet": 10,
    "nc_nCells": 500,
    "nc_nComp": 3,
    "nc_scaleScores": True,
    "nc_symmetric": False,
    "nc_q": 0.90,
    "td_K": 3,
    "td_maxIter": 1000,
    "td_maxError": 1e-5,
    "td_nDecimal": 3,
    "ma_nDim": 2,
    "nCores": 4,
    "small_donor_rule": "Only if a donor has fewer than 500 fibroblasts: nc_nCells=min(500,floor(0.8*cell_count)).",
}
TARGETS = (
    {
        "gene": "TYMS",
        "condition": "normal",
        "donors": ("Normal_1", "Normal_2", "Normal_3"),
        "disease_direction": "decreased_in_nonDM_ED_vs_normal",
        "virtual_ko_context": "normal_fibroblast",
        "interpretation_limit": "Network effect of simulated TYMS loss; not disease causality.",
        "expected_new_logFC": -2.37146446967605,
        "expected_new_FDR": 0.0393705997310064,
    },
    {
        "gene": "EFNB2",
        "condition": "organic_ED_nonDM",
        "donors": ("non-DM_1", "non-DM_2", "non-DM_3"),
        "disease_direction": "increased_in_nonDM_ED_vs_normal",
        "virtual_ko_context": "nonDM_ED_fibroblast",
        "interpretation_limit": "Potential ED-network dependency; not predicted reversal of ED.",
        "expected_new_logFC": 1.81594592717166,
        "expected_new_FDR": 0.0393705997310064,
    },
    {
        "gene": "LRRC17",
        "condition": "organic_ED_nonDM",
        "donors": ("non-DM_1", "non-DM_2", "non-DM_3"),
        "disease_direction": "increased_in_nonDM_ED_vs_normal",
        "virtual_ko_context": "nonDM_ED_fibroblast",
        "interpretation_limit": "Potential ED-network dependency; not predicted reversal of ED.",
        "expected_new_logFC": 2.15123356693459,
        "expected_new_FDR": 0.0449741882936955,
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="\n", dir=path.parent, delete=False) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    temporary.replace(path)


def write_json_atomic(path: Path, payload: object) -> None:
    write_text_atomic(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def tsv_text(fieldnames: list[str], rows: list[dict[str, object]]) -> str:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def program_version(command: list[str]) -> str:
    try:
        completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace", check=False)
        return (completed.stdout or completed.stderr).strip()
    except OSError as error:
        return f"UNAVAILABLE: {error}"


def source_target_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = {
        row["gene_symbol"]: row
        for row in rows
        if row["file"] == "differential_fibroblast_nonDM_ED_vs_normal.tsv" and row["gene_symbol"] in {item["gene"] for item in TARGETS}
    }
    if set(selected) != {item["gene"] for item in TARGETS}:
        raise ValueError("The frozen sensitivity input lacks exactly one fibroblast nonDM-ED row for every target.")
    return selected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--rscript", required=True, type=Path)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    vp_root = project_root / "revision_v2" / "07_virtual_perturbation"
    protocol_dir = vp_root / "01_protocol"
    validation_dir = vp_root / "validation"
    timestamp = utc_now()

    manifest_rows: list[dict[str, object]] = []
    input_ok = True
    for specification in EXPECTED_INPUTS:
        path = project_root / specification["relative_path"]
        exists = path.is_file()
        observed_bytes = path.stat().st_size if exists else None
        observed_sha256 = sha256(path) if exists else None
        status = "PASS" if exists and observed_bytes == specification["bytes"] and observed_sha256 == specification["sha256"] else "HOLD_INPUT_DRIFT"
        input_ok = input_ok and status == "PASS"
        manifest_rows.append(
            {
                "relative_path": specification["relative_path"],
                "bytes_expected": specification["bytes"],
                "bytes_observed": observed_bytes if observed_bytes is not None else "MISSING",
                "sha256_expected": specification["sha256"],
                "sha256_observed": observed_sha256 if observed_sha256 is not None else "MISSING",
                "status": status,
            }
        )
    write_text_atomic(
        protocol_dir / "input_manifest.sha256.tsv",
        tsv_text(list(manifest_rows[0]), manifest_rows),
    )

    sensitivity_path = project_root / EXPECTED_INPUTS[2]["relative_path"]
    source_rows = source_target_rows(sensitivity_path)
    target_rows: list[dict[str, object]] = []
    target_ok = True
    for target in TARGETS:
        source = source_rows[target["gene"]]
        observed_logfc = float(source["new_logFC"])
        observed_fdr = float(source["new_FDR"])
        matches = abs(observed_logfc - target["expected_new_logFC"]) < 1e-12 and abs(observed_fdr - target["expected_new_FDR"]) < 1e-12
        target_ok = target_ok and matches
        target_rows.append(
            {
                "gene": target["gene"],
                "source_contrast": "fibroblast nonDM-ED vs normal",
                "new_logFC": format(observed_logfc, ".14g"),
                "new_FDR": format(observed_fdr, ".14g"),
                "disease_direction": target["disease_direction"],
                "primary_virtual_KO_context": target["virtual_ko_context"],
                "primary_condition": target["condition"],
                "interpretation_limit": target["interpretation_limit"],
                "source_value_check": "PASS" if matches else "HOLD_INPUT_DRIFT",
            }
        )
    write_text_atomic(protocol_dir / "target_direction_freeze.tsv", tsv_text(list(target_rows[0]), target_rows))

    run_rows: list[dict[str, object]] = []
    order = 0
    for target in TARGETS:
        for donor in target["donors"]:
            for seed in SEEDS:
                order += 1
                run_rows.append(
                    {
                        "run_order": order,
                        "run_id": f"VP_MAIN_{order:02d}",
                        "gene": target["gene"],
                        "condition": target["condition"],
                        "donor": donor,
                        "seed": seed,
                        "matrix_directory": f"04_prepared/fibroblast/{target['condition']}/{donor}",
                        "output_attempt_directory": f"05_runs/main/{target['gene']}/{target['condition']}/{donor}/{seed}/attempt_01",
                        "status": "FROZEN_NOT_RUN",
                    }
                )
    run_design_ok = len(run_rows) == 45 and len({row["run_id"] for row in run_rows}) == 45 and len({(row["gene"], row["condition"], row["donor"], row["seed"]) for row in run_rows}) == 45
    write_text_atomic(protocol_dir / "run_matrix.tsv", tsv_text(list(run_rows[0]), run_rows))

    freeze_status = "PASS" if input_ok and target_ok and run_design_ok else "HOLD_INPUT_DRIFT"
    freeze = {
        "gate": "VP-G01",
        "status": freeze_status,
        "frozen_at": timestamp,
        "protocol_source": {
            "path": "revision_v2/07_virtual_perturbation/00_manual/VIRTUAL_PERTURBATION_WORKBOOK_zh.md",
            "version": "1.0",
            "manual_freeze_date": "2026-08-27",
        },
        "project_root": str(project_root),
        "input_manifest": "01_protocol/input_manifest.sha256.tsv",
        "input_count": len(EXPECTED_INPUTS),
        "target_direction_freeze": "01_protocol/target_direction_freeze.tsv",
        "run_matrix": "01_protocol/run_matrix.tsv",
        "main_run_count": len(run_rows),
        "random_seeds": list(SEEDS),
        "model_parameters": MODEL_PARAMETERS,
        "primary_endpoints": [
            "Full-gene differential-regulation ranking",
            "Donor-level consensus perturbation genes",
            "HALLMARK_HYPOXIA enrichment in the perturbation ranking",
            "Empirical perturbation strength versus expression- and degree-matched negative genes",
        ],
        "secondary_endpoints": [
            "Full Hallmark, Reactome, and GO Biological Process ranked enrichment",
            "ED fibroblast FDR-gene and predeclared top-ranked ED signature overlap",
            "EFNB2/LRRC17 normal-background conditional sensitivity",
            "Post-main-run descriptive extensions only",
        ],
        "result_blind_boundary": {
            "h5ad_contents_opened": False,
            "prepared_matrices_generated": False,
            "virtual_knockout_run": False,
            "biological_result_interpreted": False,
        },
        "environment": {
            "python_version": sys.version.split()[0],
            "python_executable": sys.executable,
            "rscript_path": str(args.rscript),
            "rscript_version": program_version([str(args.rscript), "--version"]),
            "platform": platform.platform(),
        },
    }
    write_json_atomic(protocol_dir / "analysis_freeze.json", freeze)

    validation = {
        "gate": "VP-G01",
        "status": "VP_G01_freeze_PASS" if freeze_status == "PASS" else "HOLD_INPUT_DRIFT",
        "validated_at": timestamp,
        "checks": {
            "exactly_seven_authoritative_inputs": len(manifest_rows) == 7,
            "all_authoritative_input_hashes_match": input_ok,
            "all_target_source_values_match": target_ok,
            "five_random_seeds_frozen": len(SEEDS) == 5,
            "forty_five_unique_main_runs_frozen": run_design_ok,
            "model_parameter_nCores_is_four": MODEL_PARAMETERS["nCores"] == 4,
            "result_blind_boundary_maintained": True,
        },
        "artifacts": [
            "01_protocol/input_manifest.sha256.tsv",
            "01_protocol/analysis_freeze.json",
            "01_protocol/target_direction_freeze.tsv",
            "01_protocol/run_matrix.tsv",
        ],
    }
    write_json_atomic(validation_dir / "VP_G01_freeze_PASS.json", validation)
    print(json.dumps({"gate": "VP-G01", "status": validation["status"], "main_run_count": len(run_rows)}, ensure_ascii=False))
    return 0 if freeze_status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
