#!/usr/bin/env python3
"""VP-ES01-02: read-only seed-rank reconstruction and VP-G05 reconciliation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP_ROOT = PROJECT_ROOT / "revision_v2" / "07_virtual_perturbation"
ATTEMPT = VP_ROOT / "10_exploratory_hypoxia_sensitivity" / "attempt_20260910_01"
PREFLIGHT = ATTEMPT / "00_preflight.json"
SCOPE = ATTEMPT / "01_scope_freeze.json"
REGISTRY = VP_ROOT / "05_runs" / "run_registry.tsv"
GENE_UNIVERSE = VP_ROOT / "04_prepared" / "gene_universe.tsv"
FROZEN_RANKINGS = VP_ROOT / "06_consensus" / "all_virtual_KO_rankings.tsv.gz"
OUTPUTS = [
    ATTEMPT / "02_seed_rank_qc.tsv",
    ATTEMPT / "02_seed_pair_stability.tsv",
    ATTEMPT / "02_seed_aggregate_reconciliation.tsv",
]
FAILURE = ATTEMPT / "02_seed_rank_reconstruction_FAILURE_RETAINED.json"

EXPECTED_DIFF_FIELDS = ["gene", "distance", "Z", "FC", "p.value", "p.adj"]
EXPECTED_RANKING_FIELDS = [
    "target_gene", "condition", "donor", "gene", "consensus_rank",
    "median_standardized_rank", "in_donor_top5pct",
]
EXPECTED_SEEDS = {"2026082701", "2026082702", "2026082703", "2026082704", "2026082705"}
TOP_K = 200
NUMERIC_TOLERANCE = 1e-15


class GateFailure(RuntimeError):
    """A retained VP-ES01-02 failure that blocks downstream biology."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_delimited(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    if path.exists():
        raise GateFailure(f"refusing to overwrite existing output: {path.name}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_frozen_gate_state() -> dict[str, str]:
    if Path.cwd().resolve() != PROJECT_ROOT.resolve():
        raise GateFailure(f"working directory is not the frozen project root: {Path.cwd()}")
    if any(path.exists() for path in OUTPUTS) or FAILURE.exists():
        raise GateFailure("VP-ES01-02 output or retained failure already exists; a new attempt is required")
    if not PREFLIGHT.is_file() or not SCOPE.is_file():
        raise GateFailure("required VP-ES01-00/01 records are missing")
    preflight = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    scope = json.loads(SCOPE.read_text(encoding="utf-8"))
    if preflight.get("status") != "PASS_INPUTS_READ_ONLY" or preflight.get("errors"):
        raise GateFailure("VP-ES01-00 is not a clean PASS_INPUTS_READ_ONLY record")
    if scope.get("status") != "FROZEN_EXPLORATORY_SPECIFICATION":
        raise GateFailure("VP-ES01-01 scope is not frozen")
    if scope.get("preflight_record", {}).get("sha256") != sha256(PREFLIGHT):
        raise GateFailure("VP-ES01-00 preflight hash differs from the hash embedded in the scope freeze")
    if scope.get("primary_ranking", {}).get("metric") != "distance":
        raise GateFailure("scope does not retain distance as the primary ranking metric")
    if scope.get("primary_ranking", {}).get("tie_break") != "gene symbol ascending":
        raise GateFailure("scope does not retain the frozen gene-symbol tie break")
    manifest = read_delimited(ATTEMPT / "input_manifest.sha256.tsv", "\t")
    manifest_hashes = {row["path"]: row["sha256"] for row in manifest if row.get("exists") == "True"}
    needed = {
        "05_runs/run_registry.tsv": REGISTRY,
        "04_prepared/gene_universe.tsv": GENE_UNIVERSE,
        "06_consensus/all_virtual_KO_rankings.tsv.gz": FROZEN_RANKINGS,
    }
    for relative, path in needed.items():
        if not path.is_file() or manifest_hashes.get(relative) != sha256(path):
            raise GateFailure(f"HOLD_INPUT_DRIFT: {relative} differs from the VP-ES01-00 input manifest")
    return {"preflight_sha256": sha256(PREFLIGHT), "scope_sha256": sha256(SCOPE)}


def load_gene_universe() -> tuple[list[str], set[str]]:
    rows = read_delimited(GENE_UNIVERSE, "\t")
    if not rows or set(rows[0]) != {"gene_order", "gene", "source_var_index_0based", "selection_reason"}:
        raise GateFailure("fixed gene universe has an unexpected schema")
    ordered = [row["gene"] for row in rows]
    if len(ordered) != 1388 or len(set(ordered)) != 1388:
        raise GateFailure(f"fixed gene universe has {len(ordered)} genes / {len(set(ordered))} unique genes, expected 1,388")
    if [int(row["gene_order"]) for row in rows] != list(range(1, 1389)):
        raise GateFailure("fixed gene universe order is not the frozen 1..1,388 sequence")
    return ordered, set(ordered)


def load_registry() -> list[dict[str, str]]:
    rows = read_delimited(REGISTRY, ",")
    required = {"run_order", "run_id", "gene", "condition", "donor", "seed", "status", "output_attempt"}
    if len(rows) != 45 or not rows or required - set(rows[0]):
        raise GateFailure("main run registry does not have the frozen 45-row schema")
    if any(row["status"] != "PASS_TECHNICAL" for row in rows):
        raise GateFailure("main run registry contains a non-PASS_TECHNICAL run")
    run_ids = [row["run_id"] for row in rows]
    identities = [(row["gene"], row["condition"], row["donor"], row["seed"]) for row in rows]
    if len(set(run_ids)) != 45 or len(set(identities)) != 45:
        raise GateFailure("main run registry contains duplicate run identity")
    groups: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for row in rows:
        groups[(row["gene"], row["condition"], row["donor"])].add(row["seed"])
    if len(groups) != 9 or any(seeds != EXPECTED_SEEDS for seeds in groups.values()):
        raise GateFailure("main run registry is not the frozen 3-target x 3-donor x 5-seed design")
    return rows


def reconstruct_seed_rank(row: dict[str, str], universe: set[str]) -> tuple[dict[str, int], dict[str, object]]:
    source = Path(row["output_attempt"]) / "differential_regulation.tsv"
    if not source.is_file() or source.stat().st_size == 0:
        raise GateFailure(f"missing or empty differential_regulation.tsv: {row['run_id']}")
    records = read_delimited(source, "\t")
    fields = list(records[0]) if records else []
    if fields != EXPECTED_DIFF_FIELDS:
        raise GateFailure(f"unexpected differential-regulation schema in {row['run_id']}: {fields}")
    genes = [record["gene"] for record in records]
    if len(records) != 1388 or len(set(genes)) != len(genes) or set(genes) != universe:
        raise GateFailure(f"gene universe mismatch in {row['run_id']}")
    try:
        distances = {record["gene"]: float(record["distance"]) for record in records}
    except ValueError as exc:
        raise GateFailure(f"non-numeric distance in {row['run_id']}: {exc}") from exc
    if not all(math.isfinite(value) for value in distances.values()):
        raise GateFailure(f"non-finite distance in {row['run_id']}")
    ordered = sorted(distances, key=lambda gene: (-distances[gene], gene))
    ranks = {gene: index for index, gene in enumerate(ordered, start=1)}
    qc = {
        "target_gene": row["gene"], "condition": row["condition"], "donor": row["donor"], "seed": row["seed"],
        "run_id": row["run_id"], "output_attempt": row["output_attempt"], "gene_universe_size": len(genes),
        "fixed_gene_universe_size": len(universe), "target_present": row["gene"] in ranks,
        "finite_distance": True, "distance_ranking": "descending", "tie_break": "gene_symbol_ascending",
        "status": "PASS_SEED_RANK_QC",
    }
    return ranks, qc


def spearman_no_ties(left: dict[str, int], right: dict[str, int]) -> float:
    genes = sorted(set(left) & set(right))
    n = len(genes)
    if n < 2:
        return math.nan
    sum_sq = sum((left[gene] - right[gene]) ** 2 for gene in genes)
    return 1.0 - (6.0 * sum_sq) / (n * (n * n - 1))


def jaccard_top_k(left: dict[str, int], right: dict[str, int]) -> float:
    a = {gene for gene, rank in left.items() if rank <= TOP_K}
    b = {gene for gene, rank in right.items() if rank <= TOP_K}
    return len(a & b) / len(a | b) if a | b else math.nan


def load_frozen_rankings() -> dict[tuple[str, str, str, str], tuple[int, float, bool]]:
    with gzip.open(FROZEN_RANKINGS, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if fields != EXPECTED_RANKING_FIELDS:
            raise GateFailure(f"all_virtual_KO_rankings schema differs from frozen VP-G05 schema: {fields}")
        source: dict[tuple[str, str, str, str], tuple[int, float, bool]] = {}
        for row in reader:
            key = (row["target_gene"], row["condition"], row["donor"], row["gene"])
            if key in source:
                raise GateFailure("duplicate row in frozen all_virtual_KO_rankings")
            source[key] = (int(row["consensus_rank"]), float(row["median_standardized_rank"]), row["in_donor_top5pct"] == "True")
    if len(source) != 12492:
        raise GateFailure(f"frozen all_virtual_KO_rankings has {len(source)} rows, expected 12,492")
    return source


def write_failure(error: Exception, state: dict[str, str] | None) -> None:
    if FAILURE.exists():
        return
    record = {
        "gate": "VP-ES01-02", "stage": "per_seed_ranking_reconstruction_and_reconciliation",
        "status": "FAIL_RECONCILIATION" if not str(error).startswith("HOLD_INPUT_DRIFT") else "HOLD_INPUT_DRIFT",
        "created_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "error": str(error), "frozen_gate_record_hashes": state or {},
        "downstream_action": "Do not enter VP-ES01-03 or interpret biology. Retain this failure and investigate in a new append-only attempt.",
    }
    FAILURE.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    state: dict[str, str] | None = None
    try:
        state = validate_frozen_gate_state()
        _, universe = load_gene_universe()
        registry = load_registry()
        frozen = load_frozen_rankings()
        rank_maps: dict[str, dict[str, int]] = {}
        qc_rows: list[dict[str, object]] = []
        grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        for row in registry:
            ranks, qc = reconstruct_seed_rank(row, universe)
            rank_maps[row["run_id"]] = ranks
            qc_rows.append(qc)
            grouped[(row["gene"], row["condition"], row["donor"])].append(row)

        pair_rows: list[dict[str, object]] = []
        aggregate_rows: list[dict[str, object]] = []
        mismatch_count = 0
        for group, runs in sorted(grouped.items()):
            target_gene, condition, donor = group
            runs = sorted(runs, key=lambda item: int(item["seed"]))
            if len(runs) != 5:
                raise GateFailure(f"expected five seeds for {group}, found {len(runs)}")
            maps = [rank_maps[row["run_id"]] for row in runs]
            for left in range(5):
                for right in range(left + 1, 5):
                    pair_rows.append({
                        "target_gene": target_gene, "condition": condition, "donor": donor,
                        "seed_a": runs[left]["seed"], "seed_b": runs[right]["seed"], "n_common_genes": len(universe),
                        "spearman_rho": spearman_no_ties(maps[left], maps[right]), "top_k": TOP_K,
                        "top200_jaccard": jaccard_top_k(maps[left], maps[right]), "status": "PASS_RECOMPUTED",
                    })
            n = len(universe)
            scores = {gene: float(statistics.median([(n - rank_map[gene] + 1) / n for rank_map in maps])) for gene in universe}
            ranked = sorted(scores, key=lambda gene: (-scores[gene], gene))
            reconstructed_ranks = {gene: index for index, gene in enumerate(ranked, start=1)}
            top_n = math.ceil(n * 0.05)
            for gene in sorted(universe):
                key = (target_gene, condition, donor, gene)
                if key not in frozen:
                    raise GateFailure(f"frozen ranking lacks reconstructed key: {key}")
                frozen_rank, frozen_score, frozen_top = frozen[key]
                rank_match = reconstructed_ranks[gene] == frozen_rank
                score_difference = abs(scores[gene] - frozen_score)
                score_match = score_difference <= NUMERIC_TOLERANCE
                top_match = (reconstructed_ranks[gene] <= top_n) == frozen_top
                comparison = "PASS_EXACT_RECONCILIATION" if rank_match and score_match and top_match else "FAIL_RECONCILIATION"
                mismatch_count += comparison != "PASS_EXACT_RECONCILIATION"
                aggregate_rows.append({
                    "target_gene": target_gene, "condition": condition, "donor": donor, "gene": gene,
                    "reconstructed_rank": reconstructed_ranks[gene], "frozen_rank": frozen_rank,
                    "rank_match": rank_match, "reconstructed_median_standardized_rank": scores[gene],
                    "frozen_median_standardized_rank": frozen_score, "absolute_score_difference": score_difference,
                    "reconstructed_in_donor_top5pct": reconstructed_ranks[gene] <= top_n,
                    "frozen_in_donor_top5pct": frozen_top, "top5pct_match": top_match,
                    "comparison_status": comparison,
                })
        if len(qc_rows) != 45 or len(pair_rows) != 90 or len(aggregate_rows) != 12492:
            raise GateFailure("VP-ES01-02 output cardinality mismatch")
        if mismatch_count:
            for row in pair_rows:
                row["status"] = "NOT_INTERPRETED_RECONCILIATION_FAILED"
        write_tsv(OUTPUTS[0], qc_rows, list(qc_rows[0]))
        write_tsv(OUTPUTS[1], pair_rows, list(pair_rows[0]))
        write_tsv(OUTPUTS[2], aggregate_rows, list(aggregate_rows[0]))
        if mismatch_count:
            raise GateFailure(f"{mismatch_count} of 12,492 rows failed VP-G05 reconciliation")
        print(json.dumps({"gate": "VP-ES01-02", "status": "PASS_RECONCILIATION", "seed_qc_rows": 45, "seed_pair_rows": 90, "aggregate_rows": 12492, "mismatches": 0}, ensure_ascii=False))
        return 0
    except Exception as exc:
        write_failure(exc, state)
        print(json.dumps({"gate": "VP-ES01-02", "status": "FAILED_OR_HELD", "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
