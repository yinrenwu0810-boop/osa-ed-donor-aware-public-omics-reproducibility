#!/usr/bin/env python3
"""Run the frozen bounded NCBI GEO and PubMed searches for VP-G07."""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


VP = Path(__file__).resolve().parents[1]
OUT_DIR = VP / "08_external_reference"
OUT = OUT_DIR / "ncbi_search_results.json"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TARGETS = ["EFNB2", "LRRC17", "TYMS"]


def get_json(endpoint: str, parameters: dict[str, str]) -> dict:
    url = f"{BASE}/{endpoint}?{urllib.parse.urlencode(parameters)}"
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "VP-G07-audit/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.load(response)
    time.sleep(0.4)
    return payload


def search(db: str, term: str) -> tuple[dict, list[dict]]:
    payload = get_json("esearch.fcgi", {"db": db, "term": term, "retmode": "json", "retmax": "100", "retstart": "0"})
    result = payload["esearchresult"]
    ids = result.get("idlist", [])
    records: list[dict] = []
    if ids:
        summary = get_json("esummary.fcgi", {"db": db, "id": ",".join(ids), "retmode": "json"})["result"]
        for uid in summary.get("uids", []):
            item = summary[uid]
            if db == "gds":
                records.append({
                    "uid": uid,
                    "accession": item.get("accession", ""),
                    "title": item.get("title", ""),
                    "summary": item.get("summary", ""),
                    "taxon": item.get("taxon", ""),
                    "entrytype": item.get("entrytype", ""),
                    "gdstype": item.get("gdstype", ""),
                    "n_samples": item.get("n_samples", ""),
                    "pubmedids": item.get("pubmedids", []),
                    "pdat": item.get("pdat", ""),
                    "ftp_link": item.get("ftplink", ""),
                })
            else:
                article_ids = {entry.get("idtype", ""): entry.get("value", "") for entry in item.get("articleids", [])}
                records.append({
                    "uid": uid,
                    "title": item.get("title", ""),
                    "authors": [author.get("name", "") for author in item.get("authors", [])],
                    "fulljournalname": item.get("fulljournalname", ""),
                    "pubdate": item.get("pubdate", ""),
                    "doi": article_ids.get("doi", ""),
                    "pmcid": article_ids.get("pmc", ""),
                })
    meta = {
        "db": db,
        "endpoint": f"{BASE}/esearch.fcgi and esummary.fcgi",
        "term": term,
        "server_count": int(result["count"]),
        "retmax": 100,
        "retrieved_count": len(ids),
        "complete_within_bound": int(result["count"]) == len(ids),
        "query_translation": result.get("querytranslation", ""),
    }
    return meta, records


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite existing search result: {OUT}")
    OUT_DIR.mkdir(parents=True, exist_ok=False)
    searches = []
    for target in TARGETS:
        concepts = {
            "single_cell_perturbation": f'"{target}"[All Fields] AND ("perturb-seq"[All Fields] OR "CROP-seq"[All Fields] OR ("single cell"[All Fields] AND "CRISPR"[All Fields]))',
            "related_cell_genetic": f'"{target}"[All Fields] AND (fibroblast[All Fields] OR endothelial[All Fields] OR "smooth muscle"[All Fields] OR pericyte[All Fields] OR "corpus cavernosum"[All Fields] OR penile[All Fields]) AND (CRISPR[All Fields] OR knockdown[All Fields] OR knockout[All Fields] OR overexpression[All Fields])',
            "transcriptomic_genetic": f'"{target}"[All Fields] AND (CRISPR[All Fields] OR knockdown[All Fields] OR knockout[All Fields] OR overexpression[All Fields]) AND (transcriptome[All Fields] OR transcriptomic[All Fields] OR "RNA-seq"[All Fields])',
        }
        for family, base_term in concepts.items():
            for db in ["gds", "pubmed"]:
                term = f"({base_term}) AND gse[ETYP]" if db == "gds" else base_term
                meta, records = search(db, term)
                searches.append({"target_gene": target, "query_family": family, "search": meta, "records": records})
    payload = {
        "gate": "VP-G07",
        "accessed_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Bounded targeted NCBI search; genome-wide studies may not mention individual targets in indexed metadata and are cross-checked through catalogues and guide/condition annotations.",
        "rate_limit": "Unauthenticated NCBI E-utilities, serialized at no more than 2.5 requests per second.",
        "searches": searches,
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"path": str(OUT), "searches": len(searches), "server_total": sum(entry["search"]["server_count"] for entry in searches), "retrieved_total": sum(entry["search"]["retrieved_count"] for entry in searches)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
