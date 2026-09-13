"""VP-G05 step 1: frozen seed stability and donor consensus aggregation.

This script consumes only retained VP-G04 differential-regulation artifacts.
It performs no pathway analysis, no negative-control selection, and no
biological interpretation.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import shutil
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
REGISTRY = VP / "05_runs" / "run_registry.tsv"
G04_COMPLETE = VP / "validation" / "VP_G04_main_runs_COMPLETE.json"
OUT = VP / "06_consensus"
VALIDATION = VP / "validation" / "VP_G05_stability_PRECONTROL.json"
TOP_K = 200
TOP_FRACTION = 0.05
MIN_SPEARMAN = 0.70
MIN_JACCARD = 0.40
SEEDS_PER_DONOR = 5
MIN_DONOR_FRACTION = 2 / 3


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def stable_ranks(attempt: Path) -> dict[str, int]:
    rows = read_tsv(attempt / "differential_regulation.tsv")
    if not rows or set(rows[0]) != {"gene", "distance", "Z", "FC", "p.value", "p.adj"}:
        raise RuntimeError(f"Unexpected differential-regulation schema: {attempt}")
    ordered = sorted(((row["gene"], float(row["distance"])) for row in rows), key=lambda pair: (-pair[1], pair[0]))
    if len({gene for gene, _ in ordered}) != len(ordered):
        raise RuntimeError(f"Duplicate genes: {attempt}")
    return {gene: index for index, (gene, _) in enumerate(ordered, start=1)}


def spearman_no_ties(left: dict[str, int], right: dict[str, int]) -> float:
    genes = sorted(set(left) & set(right))
    n = len(genes)
    if n < 2:
        return math.nan
    sum_sq = sum((left[gene] - right[gene]) ** 2 for gene in genes)
    return 1.0 - (6.0 * sum_sq) / (n * (n * n - 1))


def jaccard_top_k(left: dict[str, int], right: dict[str, int], k: int) -> float:
    a = {gene for gene, rank in left.items() if rank <= k}
    b = {gene for gene, rank in right.items() if rank <= k}
    return len(a & b) / len(a | b) if a | b else math.nan


def median(values: list[float]) -> float:
    return float(statistics.median(values)) if values else math.nan


def main() -> None:
    if not G04_COMPLETE.is_file() or json.loads(G04_COMPLETE.read_text(encoding="utf-8")).get("status") != "PASS":
        raise RuntimeError("VP-G05 requires a PASS VP-G04 complete audit.")
    if OUT.exists() or VALIDATION.exists():
        raise RuntimeError("VP-G05 output path already exists; refusing overwrite.")
    registry = read_csv(REGISTRY)
    if len(registry) != 45 or {row["status"] for row in registry} != {"PASS_TECHNICAL"}:
        raise RuntimeError("VP-G05 requires exactly 45 retained PASS_TECHNICAL main runs.")

    temporary = OUT.with_name(OUT.name + ".tmp")
    temporary.mkdir(parents=True, exist_ok=False)
    try:
        ranking_freeze = {
            "gate": "VP-G05",
            "status": "FROZEN_BEFORE_AGGREGATE_OUTPUTS",
            "ranking_metric": "distance",
            "ranking_direction": "descending; larger network perturbation distance ranks higher",
            "tie_break": "gene symbol ascending",
            "seed_aggregation": "within donor, median of per-seed standardized ranks (n-rank+1)/n",
            "top_k": TOP_K,
            "top_fraction": TOP_FRACTION,
            "seed_stability_threshold_spearman": MIN_SPEARMAN,
            "seed_stability_threshold_jaccard": MIN_JACCARD,
            "donor_consensus_rule": "gene is consensus when it enters top 5% in at least 2/3 donors",
            "prohibited": ["pathway inspection", "significance-only filtering", "result-driven seed selection"],
        }
        (temporary / "VP_G05_ranking_rule_freeze.json").write_text(json.dumps(ranking_freeze, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

        grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        for row in registry:
            grouped[(row["gene"], row["condition"], row["donor"])].append(row)
        rank_maps: dict[int, dict[str, int]] = {}
        for row in registry:
            rank_maps[int(row["run_order"])] = stable_ranks(Path(row["output_attempt"]))

        seed_spearman: list[dict[str, object]] = []
        seed_jaccard: list[dict[str, object]] = []
        donor_stability: list[dict[str, object]] = []
        donor_scores: dict[tuple[str, str, str], dict[str, float]] = {}
        donor_ranks: dict[tuple[str, str, str], dict[str, int]] = {}

        for key, runs in sorted(grouped.items()):
            gene, condition, donor = key
            runs = sorted(runs, key=lambda item: int(item["seed"]))
            if len(runs) != SEEDS_PER_DONOR:
                raise RuntimeError(f"Expected {SEEDS_PER_DONOR} seeds for {key}, found {len(runs)}")
            maps = [rank_maps[int(row["run_order"])] for row in runs]
            universe = set(maps[0])
            if any(set(rank_map) != universe for rank_map in maps[1:]):
                raise RuntimeError(f"Gene universe mismatch among seeds for {key}")
            pair_rhos, pair_jaccards = [], []
            for left in range(len(runs)):
                for right in range(left + 1, len(runs)):
                    rho = spearman_no_ties(maps[left], maps[right])
                    jac = jaccard_top_k(maps[left], maps[right], TOP_K)
                    pair_rhos.append(rho)
                    pair_jaccards.append(jac)
                    common = len(set(maps[left]) & set(maps[right]))
                    seed_spearman.append({"target_gene": gene, "condition": condition, "donor": donor, "seed_a": runs[left]["seed"], "seed_b": runs[right]["seed"], "n_common_genes": common, "spearman_rho": rho})
                    seed_jaccard.append({"target_gene": gene, "condition": condition, "donor": donor, "seed_a": runs[left]["seed"], "seed_b": runs[right]["seed"], "top_k": TOP_K, "jaccard": jac})
            med_rho, med_jaccard = median(pair_rhos), median(pair_jaccards)
            stable = med_rho >= MIN_SPEARMAN and med_jaccard >= MIN_JACCARD
            donor_stability.append({"target_gene": gene, "condition": condition, "donor": donor, "n_seeds": len(runs), "pair_count": len(pair_rhos), "median_seed_spearman": med_rho, "median_top200_jaccard": med_jaccard, "seed_stability_status": "PASS" if stable else "UNSTABLE_VIRTUAL_PERTURBATION"})
            n = len(universe)
            scores = {symbol: median([(n - rank_map[symbol] + 1) / n for rank_map in maps]) for symbol in universe}
            ranked = sorted(scores, key=lambda symbol: (-scores[symbol], symbol))
            donor_scores[key] = scores
            donor_ranks[key] = {symbol: index for index, symbol in enumerate(ranked, start=1)}

        donor_pairwise: list[dict[str, object]] = []
        consensus_rows: list[dict[str, object]] = []
        full_rankings: list[dict[str, object]] = []
        target_status: list[dict[str, object]] = []
        target_groups: dict[tuple[str, str], list[tuple[str, str, str]]] = defaultdict(list)
        for key in donor_ranks:
            target_groups[key[:2]].append(key)
        for target_key, donor_keys in sorted(target_groups.items()):
            gene, condition = target_key
            donor_keys = sorted(donor_keys, key=lambda item: item[2])
            pairwise_rhos = []
            for left in range(len(donor_keys)):
                for right in range(left + 1, len(donor_keys)):
                    rho = spearman_no_ties(donor_ranks[donor_keys[left]], donor_ranks[donor_keys[right]])
                    pairwise_rhos.append(rho)
                    donor_pairwise.append({"target_gene": gene, "condition": condition, "donor_a": donor_keys[left][2], "donor_b": donor_keys[right][2], "n_common_genes": len(donor_ranks[donor_keys[left]]), "spearman_rho": rho})
            stable_donors = [row for row in donor_stability if row["target_gene"] == gene and row["condition"] == condition and row["seed_stability_status"] == "PASS"]
            target_stable = len(stable_donors) >= math.ceil(len(donor_keys) * MIN_DONOR_FRACTION)
            target_status.append({"target_gene": gene, "condition": condition, "donors": len(donor_keys), "stable_donors": len(stable_donors), "median_donor_spearman": median(pairwise_rhos), "status": "PASS_STABILITY" if target_stable else "UNSTABLE_VIRTUAL_PERTURBATION"})
            genes = sorted(donor_ranks[donor_keys[0]])
            top_n = math.ceil(len(genes) * TOP_FRACTION)
            for symbol in genes:
                donor_top_count = sum(donor_ranks[key][symbol] <= top_n for key in donor_keys)
                median_score = median([donor_scores[key][symbol] for key in donor_keys])
                consensus_rows.append({"target_gene": gene, "condition": condition, "gene": symbol, "donors_top5pct": donor_top_count, "donors_total": len(donor_keys), "is_2of3_top5pct_consensus": donor_top_count >= math.ceil(len(donor_keys) * MIN_DONOR_FRACTION), "median_standardized_rank": median_score})
                for donor_key in donor_keys:
                    full_rankings.append({"target_gene": gene, "condition": condition, "donor": donor_key[2], "gene": symbol, "consensus_rank": donor_ranks[donor_key][symbol], "median_standardized_rank": donor_scores[donor_key][symbol], "in_donor_top5pct": donor_ranks[donor_key][symbol] <= top_n})

        write_tsv(temporary / "seed_pairwise_spearman.tsv", seed_spearman, ["target_gene", "condition", "donor", "seed_a", "seed_b", "n_common_genes", "spearman_rho"])
        write_tsv(temporary / "seed_pairwise_top200_jaccard.tsv", seed_jaccard, ["target_gene", "condition", "donor", "seed_a", "seed_b", "top_k", "jaccard"])
        write_tsv(temporary / "donor_stability.tsv", donor_stability, ["target_gene", "condition", "donor", "n_seeds", "pair_count", "median_seed_spearman", "median_top200_jaccard", "seed_stability_status"])
        write_tsv(temporary / "donor_pairwise_spearman.tsv", donor_pairwise, ["target_gene", "condition", "donor_a", "donor_b", "n_common_genes", "spearman_rho"])
        write_tsv(temporary / "candidate_donor_consensus.tsv", consensus_rows, ["target_gene", "condition", "gene", "donors_top5pct", "donors_total", "is_2of3_top5pct_consensus", "median_standardized_rank"])
        with gzip.open(temporary / "all_virtual_KO_rankings.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["target_gene", "condition", "donor", "gene", "consensus_rank", "median_standardized_rank", "in_donor_top5pct"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(full_rankings)
        report = {
            "gate": "VP-G05",
            "stage": "seed_and_donor_stability_before_controls",
            "status": "PASS_STABILITY_PRECONTROL",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "target_status": target_status,
            "output_contract": "Full rankings retained; pathway analysis and negative-control KO not performed in this stage.",
            "boundary": "Network perturbation ranking only; no pathway direction, expression reversal, or causal inference.",
        }
        (temporary / "VP_G05_stability_PRECONTROL.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        temporary.rename(OUT)
        VALIDATION.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"status": report["status"], "targets": len(target_status), "output": str(OUT)}, ensure_ascii=False))
    except Exception:
        failed = OUT.with_name(OUT.name + ".attempt_01_FAIL")
        if failed.exists():
            raise
        temporary.rename(failed)
        raise


if __name__ == "__main__":
    main()
