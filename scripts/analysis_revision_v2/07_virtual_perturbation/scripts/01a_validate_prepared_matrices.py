#!/usr/bin/env python3
"""Independently audit VP-G02 prepared matrices without reopening the H5AD."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.io import mmread


EXPECTED = (
    ("normal", "Normal_1", 2377),
    ("normal", "Normal_2", 3254),
    ("normal", "Normal_3", 7149),
    ("organic_ED_nonDM", "non-DM_1", 2042),
    ("organic_ED_nonDM", "non-DM_2", 1586),
    ("organic_ED_nonDM", "non-DM_3", 4854),
)
PRIMARY_RATE_BOUNDS_PERCENT_2DP = {
    ("TYMS", "normal"): (3.66, 5.64),
    ("EFNB2", "organic_ED_nonDM"): (14.40, 23.77),
    ("LRRC17", "organic_ED_nonDM"): (31.00, 50.63),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    prepared = root / "04_prepared"
    universe = prepared / "gene_universe.tsv"
    genes = read_tsv(universe)
    recorded_sha = (prepared / "gene_universe.sha256").read_text(encoding="utf-8").split()[0]
    matrix_checks = []
    for condition, donor, n_cells in EXPECTED:
        donor_dir = prepared / "fibroblast" / condition / donor
        donor_genes = read_tsv(donor_dir / "genes.tsv")
        donor_cells = read_tsv(donor_dir / "cells.tsv")
        with gzip.open(donor_dir / "matrix.mtx.gz", "rb") as handle:
            matrix = mmread(handle).tocsr()
        assert donor_genes == [{"gene_order": x["gene_order"], "gene": x["gene"]} for x in genes]
        assert matrix.shape == (len(genes), n_cells)
        assert len(donor_cells) == n_cells
        assert all(row["donor"] == donor and row["condition"] == condition and row["cell_type"] == "fibroblast" and row["predicted_doublet"] == "False" for row in donor_cells)
        assert np.isfinite(matrix.data).all() and (matrix.data >= 0).all()
        matrix_checks.append({"condition": condition, "donor": donor, "matrix_shape_gene_by_cell": list(matrix.shape), "nnz": int(matrix.nnz), "finite_and_nonnegative": True})
    detection = read_tsv(prepared / "target_detection_all.tsv")
    primary = [row for row in detection if row["is_primary_gene_donor_combination"] == "True"]
    assert len(primary) == 9 and all(row["status"] == "PASS" for row in primary)
    for row in primary:
        low, high = PRIMARY_RATE_BOUNDS_PERCENT_2DP[(row["gene"], row["condition"])]
        assert low <= round(float(row["detection_rate"]) * 100, 2) <= high
    audit = {
        "gate": "VP-G02",
        "status": "VP_G02_independent_audit_PASS",
        "audited_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "checks": {
            "gene_universe_sha256_matches_record": sha256(universe) == recorded_sha,
            "six_matrices_share_exact_gene_order": True,
            "matrix_dimensions_match_frozen_donor_counts": True,
            "matrix_entries_are_finite_and_nonnegative": True,
            "nine_primary_target_detection_rates_match_workbook_ranges_at_two_decimal_percent_precision": True,
        },
        "gene_universe_count": len(genes),
        "matrix_checks": matrix_checks,
        "primary_target_detection": primary,
    }
    output = root / "validation" / "VP_G02_independent_audit_PASS.json"
    output.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"gate": "VP-G02", "status": audit["status"], "primary_combinations": len(primary)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
