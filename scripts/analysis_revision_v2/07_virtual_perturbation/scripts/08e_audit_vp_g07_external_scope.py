#!/usr/bin/env python3
"""Independently audit VP-G07 evidence artifacts and, only on PASS, seal completion."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
OUT = VP / "08_external_reference"
VAL = VP / "validation"
AUDIT = VAL / "VP_G07_external_scope_AUDIT01.json"
COMPLETE = VAL / "VP_G07_external_scope_COMPLETE.json"
REQUIRED = [
    OUT / "public_perturbation_target_coverage.tsv", OUT / "search_log.tsv", OUT / "source_evidence.tsv",
    OUT / "gears_applicability.tsv", OUT / "unipert_g2cp_applicability.json", OUT / "VP_G07_external_scope_report.md",
    OUT / "artifact_manifest.sha256.tsv",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if AUDIT.exists() or COMPLETE.exists():
        raise SystemExit("Refusing to overwrite existing VP-G07 audit or completion record")
    checks: list[dict[str, object]] = []
    checks.append({"name":"required_artifacts_exist", "pass":all(p.is_file() for p in REQUIRED), "detail":[str(p) for p in REQUIRED if not p.is_file()]})
    if not all(p.is_file() for p in REQUIRED):
        payload = {"gate":"VP-G07", "status":"HOLD", "checks":checks, "generated_at_utc":datetime.now(timezone.utc).isoformat()}
        AUDIT.write_text(json.dumps(payload, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        raise SystemExit("HOLD: missing required artifacts")
    coverage = rows(REQUIRED[0])
    covered = {row["target_gene"] for row in coverage}
    checks.append({"name":"all_targets_explicit", "pass":set(["EFNB2","LRRC17","TYMS"]).issubset(covered), "detail":sorted(covered)})
    checks.append({"name":"lrrc17_bounded_negative_retained", "pass":any(r["target_gene"]=="LRRC17" and r["evidence_status"]=="NO_VERIFIED_HIT" for r in coverage), "detail":"required explicit negative row"})
    checks.append({"name":"required_coverage_fields_nonempty", "pass":all(all(value.strip() for value in row.values()) for row in coverage), "detail":f"{len(coverage)} rows"})
    checks.append({"name":"no_unjustified_T1_claim", "pass":not any(r["cell_context_tier"]=="T1_EXACT" or r["matches_corpus_cavernosum_fibroblast"]=="YES" for r in coverage), "detail":"no exact corpus-cavernosum fibroblast record"})
    gears = rows(REQUIRED[3])
    checks.append({"name":"gears_not_cross_cell_validation", "pass":all(r["is_direct_external_reference"]=="NO" and r["gears_adjudication"] not in {"DIRECT_EXTERNAL_REFERENCE","RELATED_CONTEXT_REFERENCE"} for r in gears), "detail":[r["gears_adjudication"] for r in gears]})
    g2cp = json.loads(REQUIRED[4].read_text(encoding="utf-8"))
    checks.append({"name":"g2cp_all_conditions_adjudicated", "pass":g2cp.get("g2cp_status")=="NOT_APPLICABLE_YET" and len(g2cp.get("g2cp_start_conditions", []))==5 and all(x.get("status")=="UNMET" for x in g2cp["g2cp_start_conditions"]), "detail":g2cp.get("g2cp_status")})
    search = rows(REQUIRED[1])
    checks.append({"name":"sealed_ncbi_query_provenance_retained", "pass":sum(1 for r in search if r["search_id"].startswith("NCBI_"))==18, "detail":f"{sum(1 for r in search if r['search_id'].startswith('NCBI_'))} NCBI rows"})
    evidence = rows(REQUIRED[2])
    checks.append({"name":"source_urls_and_exclusions_retained", "pass":all(r["source_url"].strip() for r in evidence) and any(r["inclusion_decision"]=="SEARCH_HIT_NOT_COVERAGE" for r in evidence), "detail":f"{len(evidence)} source rows"})
    manifest = rows(REQUIRED[6])
    manifest_ok = True
    for entry in manifest:
        p = VP / entry["relative_path"]
        manifest_ok = manifest_ok and p.is_file() and digest(p)==entry["sha256"] and str(p.stat().st_size)==entry["bytes"]
    checks.append({"name":"artifact_manifest_reconciles", "pass":manifest_ok, "detail":f"{len(manifest)} entries"})
    passed = all(bool(c["pass"]) for c in checks)
    timestamp = datetime.now(timezone.utc).isoformat()
    payload = {"gate":"VP-G07", "audit_id":"AUDIT01", "status":"PASS" if passed else "HOLD", "generated_at_utc":timestamp, "independence":"This script does not import the builder and re-reads artifact contents, manifest hashes, frozen labels, and prohibited claims.", "checks":checks, "errors":[] if passed else [c["name"] for c in checks if not c["pass"]]}
    AUDIT.write_text(json.dumps(payload, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    if not passed:
        raise SystemExit("HOLD: independent audit failed; completion record not written")
    completion = {"gate":"VP-G07", "status":"PASS", "completed_at_utc":timestamp, "technical_completion":"All frozen expected G07 evidence artifacts were built and AUDIT01 independently passed.", "exact_T1_external_reference":"NO_VERIFIED_T1_EXACT_RECORD", "gears_status":"No direct external validation; no cross-cell GEARS prediction or combination claim was made.", "unipert_g2cp_status":"NOT_APPLICABLE_YET", "scientific_interpretation":"T2-T5 records document public perturbation coverage only. They do not validate corpus-cavernosum fibroblast biology, ED reversal, therapeutic benefit, or causality.", "next_gate":"VP-G08 remains NOT_STARTED and requires separate explicit user authorization.", "audit_sha256":digest(AUDIT), "artifact_manifest_sha256":digest(OUT / "artifact_manifest.sha256.tsv"), "implementation_sha256":{"scripts/08d_build_vp_g07_external_scope.py":digest(VP / "scripts" / "08d_build_vp_g07_external_scope.py"), "scripts/08e_audit_vp_g07_external_scope.py":digest(VP / "scripts" / "08e_audit_vp_g07_external_scope.py")}}
    COMPLETE.write_text(json.dumps(completion, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status":"PASS", "audit":str(AUDIT), "complete":str(COMPLETE)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
