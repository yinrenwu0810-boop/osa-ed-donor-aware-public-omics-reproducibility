"""Independent read-only audit for VP-ES01 attempt_20260910_03.

This script deliberately does not import the report/figure builder.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image

ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
OUT = VP / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_03"
OLD = VP / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_02"
VALIDATION = VP / "validation"
AUDIT = VALIDATION / "VP_ES01_exploratory_sensitivity_AUDIT01.json"
COMPLETE = VALIDATION / "VP_ES01_exploratory_sensitivity_COMPLETE.json"
TARGETS = {"EFNB2", "LRRC17", "TYMS"}
SPECS = {"S0", "S1", "S2", "S3"}


def rd(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rd_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(name: str, passed: bool, evidence: str) -> dict[str, object]:
    return {"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence}


def write_new(path: Path, payload: dict[str, object]) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite {path.name}")
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    if Path.cwd().resolve() != ROOT.resolve():
        raise RuntimeError("Run only from the frozen acad project root")
    if AUDIT.exists() or COMPLETE.exists():
        raise RuntimeError("Refusing to overwrite existing ES01 audit or complete state")
    mandatory = [OUT / "VP_ES01_exploratory_sensitivity_REPORT.md", OUT / "artifact_manifest.sha256.tsv", OUT / "Figure_ES01_hypoxia_robustness.png", OUT / "Figure_ES02_leading_edge_stability.png"]
    if any(not path.exists() for path in mandatory):
        raise RuntimeError("Report, manifest, or figures are missing")

    checks: list[dict[str, object]] = []
    manifest_input = rd(OUT / "input_manifest.sha256.tsv")
    hash_mismatches = []
    for row in manifest_input:
        path = VP / row["path"]
        if not path.exists() or str(path.stat().st_size) != row["bytes"] or sha256(path) != row["sha256"]:
            hash_mismatches.append(row["path"])
    checks.append(check("01_frozen_input_hashes", not hash_mismatches, f"verified {len(manifest_input)} manifest entries; mismatches={hash_mismatches}"))

    registry = rd_csv(VP / "05_runs" / "run_registry.tsv")
    control_tasks = rd(VP / "01_protocol" / "control_execution_v2" / "control_task_matrix.tsv")
    main_ok = len(registry) == 45 and all(row["status"] == "PASS_TECHNICAL" for row in registry) and len({(r["gene"], r["condition"], r["donor"], r["seed"]) for r in registry}) == 45
    control_ok = len(control_tasks) == 90 and len({(r["control_gene"], r["condition"], r["donor"], r["seed"]) for r in control_tasks}) == 90
    checks.append(check("02_main_and_control_runs_unchanged", main_ok and control_ok, f"main={len(registry)} unique PASS_TECHNICAL; control tasks={len(control_tasks)} unique"))

    structure = {(r["gene"], r["condition"], r["donor"]) for r in registry}
    donor_counts = {target: len({donor for gene, condition, donor in structure if gene == target}) for target in TARGETS}
    seed_counts = {(gene, condition, donor): sum(1 for r in registry if (r["gene"], r["condition"], r["donor"]) == (gene, condition, donor)) for gene, condition, donor in structure}
    checks.append(check("03_target_donor_seed_structure", set(r["gene"] for r in registry) == TARGETS and all(n == 3 for n in donor_counts.values()) and all(n == 5 for n in seed_counts.values()), f"donors={donor_counts}; seed counts={sorted(set(seed_counts.values()))}"))

    mv = rd(OUT / "06_analysis_multiverse.tsv")
    formal = rd(VP / "07_enrichment" / "v1" / "hypoxia_primary_summary.tsv")
    formal_by_target = {r["target_gene"]: r for r in formal}
    s0 = [r for r in mv if r["specification"] == "S0"]
    s0_ok = len(s0) == 3
    for row in s0:
        ref = formal_by_target.get(row["target_gene"], {})
        s0_ok = s0_ok and row["result_status"] == "PASS_NUMERIC_RECONCILIATION" and row["S0_leading_edge_exact_match"] == "TRUE" and all(float(row[key]) <= 1e-12 for key in ["rank_input_max_abs_difference", "S0_ES_abs_difference", "S0_NES_abs_difference", "S0_p_value_abs_difference", "S0_FDR_abs_difference"]) and all(abs(float(row[key]) - float(ref[source])) <= 1e-12 for key, source in [("NES", "NES"), ("p_value", "pval"), ("FDR_within_Hallmark", "padj")])
    checks.append(check("04_S0_numeric_reproduction", s0_ok, "S0 compared independently to formal hypoxia_primary_summary.tsv with zero-tolerance fields"))

    leave_seed = rd(OLD / "04_leave_one_seed_out.tsv")
    leave_donor = rd(OLD / "04_leave_one_donor_out.tsv")
    checks.append(check("05_complete_leave_one_analyses", len(leave_seed) == 45 and len(leave_donor) == 9 and len({(r["target_gene"], r["donor"], r["excluded_seed"]) for r in leave_seed}) == 45 and len({(r["target_gene"], r["excluded_donor"]) for r in leave_donor}) == 9, f"leave-one-seed={len(leave_seed)}; leave-one-donor={len(leave_donor)}"))

    grid_ok = {(r["target_gene"], r["specification"]) for r in mv} == {(target, spec) for target in TARGETS for spec in SPECS}
    non_sig_present = sum(r["FDR_within_Hallmark"] != "NA" and float(r["FDR_within_Hallmark"]) >= 0.05 for r in mv if r["specification"] in {"S0", "S1", "S2"}) == 9
    checks.append(check("06_all_predeclared_specifications_reported", grid_ok and non_sig_present, f"S0-S3 grid={len(mv)} rows; all nine rank-level specifications retain non-significant FDR"))

    family_ok = all(r["Hallmark_tested_set_count"] == "32" for r in mv if r["specification"] in {"S0", "S1", "S2"})
    checks.append(check("07_formal_Hallmark_FDR_family_preserved", family_ok, "S0-S2 each record the unchanged 32-set tested Hallmark family"))

    controls = rd(OLD / "05_target_vs_matched_controls.tsv")
    selection = rd(VP / "06_consensus" / "negative_control_selection_v2" / "negative_control_selection.tsv")
    selected_counts = {target: sum(r["target_gene"] == target for r in selection) for target in TARGETS}
    control_ok = len(controls) == 3 and all(r["matched_controls"] == "10" for r in controls) and all(n == 30 for n in selected_counts.values())
    checks.append(check("08_matched_controls_not_reselected", control_ok, f"reported controls=10 per target; frozen donor-level selection rows={selected_counts}"))

    report = (OUT / "VP_ES01_exploratory_sensitivity_REPORT.md").read_text(encoding="utf-8")
    biology_ok = "供体是唯一生物学单位" in report and "种子仅为网络估计稳定性单位" in report and "n=45" not in report
    checks.append(check("09_biological_unit_boundary", biology_ok, "report explicitly separates donor biological units from seed stability units"))

    prohibited = ["证实三个基因调控缺氧通路", "证明OSA通过三个基因导致ED", "证明CPAP通过三个基因逆转ED", "虚拟KO显示通路被上调", "虚拟KO显示通路被下调", "增加种子相当于增加样本量", "验证完成"]
    bad_language = [phrase for phrase in prohibited if phrase in report]
    source = [r for r in rd(OUT / "06_hypoxia_mapping_qc.tsv") if r["row_type"] == "SOURCE_MAPPING"]
    fig_ok = True
    for figure in mandatory[2:]:
        with Image.open(figure) as image:
            fig_ok = fig_ok and image.width >= 1000 and image.height >= 700
    report_structure_ok = len(re.findall(r"^# ", report, flags=re.M)) == 1 and "accTitle:" in report and "accDescr:" in report
    checks.append(check("10_report_boundary_and_artifact_quality", not bad_language and len(source) == 1 and source[0]["mapping_qc_status"] == "PASS_EXACT_SYMBOL_MAPPING" and fig_ok and report_structure_ok, f"prohibited_phrases={bad_language}; figure_QA={fig_ok}; report_structure={report_structure_ok}"))

    manifest = rd(OUT / "artifact_manifest.sha256.tsv")
    manifest_ok = True
    for row in manifest:
        path = OUT / row["file"]
        manifest_ok = manifest_ok and path.exists() and str(path.stat().st_size) == row["bytes"] and sha256(path) == row["sha256"]
    if not manifest_ok:
        checks[0]["status"] = "FAIL"
        checks[0]["evidence"] += "; attempt artifact manifest verification failed"
    passed = all(item["status"] == "PASS" for item in checks)
    audit_payload = {"gate": "VP-ES01-07", "attempt": "attempt_20260910_03", "status": "PASS_10_OF_10" if passed else "FAIL_AUDIT", "created_at_utc": datetime.now(timezone.utc).isoformat(), "checks": checks, "report_sha256": sha256(OUT / "VP_ES01_exploratory_sensitivity_REPORT.md"), "manifest_entries": len(manifest), "builder_imported": False}
    write_new(AUDIT, audit_payload)
    if not passed:
        raise RuntimeError("Independent audit failed; COMPLETE was not written")
    complete_payload = {"gate": "VP-ES01", "attempt": "attempt_20260910_03", "status": "COMPLETE", "completion_meaning": "Exploratory analysis artifacts and independent audit are complete; this is not biological validation or a replacement of formal VP-G06.", "final_interpretation": "ROBUST_POSITIVE_DIRECTION_BUT_NONSPECIFIC", "formal_result_protection": "VP-G06 HALLMARK_HYPOXIA formal results remain unchanged and not FDR-significant.", "audit_file": AUDIT.name, "audit_sha256": sha256(AUDIT), "created_at_utc": datetime.now(timezone.utc).isoformat()}
    write_new(COMPLETE, complete_payload)
    print(json.dumps({"gate": "VP-ES01-07", "status": "COMPLETE", "audit": "PASS_10_OF_10", "final_interpretation": complete_payload["final_interpretation"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
