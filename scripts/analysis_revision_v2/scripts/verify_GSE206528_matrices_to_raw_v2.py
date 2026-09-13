#!/usr/bin/env python3
"""Chunked, value-level audit from the eight deposited matrices to raw h5ad."""

from __future__ import annotations

import hashlib
import json
import os
import tarfile
import gc
import csv
import gzip
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
MATRIX_DIR = ROOT / "analysis/data/processed/GSE206528_matrices"
ARCHIVE = ROOT / "analysis/data/raw/GSE206528/GSE206528_RAW.tar"
RAW = ROOT / "analysis/data/processed/GSE206528_raw.h5ad"
VALIDATION = ROOT / "revision_v2/04_validation"


def write_numbered_attempt(
    validation_dir: Path, prefix: str, status: str, payload: str
) -> Path:
    """Create the first unused numbered attempt without overwriting history."""
    attempt_number = 1
    while True:
        number = f"{attempt_number:02d}"
        if any(validation_dir.glob(f"{prefix}_{number}_*.json")):
            attempt_number += 1
            continue
        attempt = validation_dir / f"{prefix}_{number}_{status}.json"
        try:
            with attempt.open("x", encoding="utf-8", newline="\n") as handle:
                handle.write(payload)
            return attempt
        except FileExistsError:
            attempt_number += 1


def atomic_replace_json(path: Path, payload: str) -> None:
    """Replace a latest pointer atomically via a PID-tagged sibling file."""
    counter = 0
    while True:
        suffix = f".{os.getpid()}.{counter}.tmp"
        temporary = path.with_name(f".{path.name}{suffix}")
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


def sha256_stream(handle, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    while block := handle.read(chunk_size):
        digest.update(block)
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return sha256_stream(handle)


def verify_archive_members(files: list[Path]) -> list[dict[str, object]]:
    expected = {path.name: path for path in files}
    rows: list[dict[str, object]] = []
    with tarfile.open(ARCHIVE, "r") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        if {Path(member.name).name for member in members} != set(expected):
            raise RuntimeError("Archive member names do not equal the eight extracted matrix names")
        for member in sorted(members, key=lambda item: item.name):
            name = Path(member.name).name
            handle = archive.extractfile(member)
            if handle is None:
                raise RuntimeError(f"Cannot read archive member {member.name}")
            archive_hash = sha256_stream(handle)
            extracted_hash = sha256_file(expected[name])
            rows.append(
                {
                    "file": name,
                    "archive_bytes": member.size,
                    "extracted_bytes": expected[name].stat().st_size,
                    "archive_sha256": archive_hash,
                    "extracted_sha256": extracted_hash,
                    "byte_and_hash_match": member.size == expected[name].stat().st_size
                    and archive_hash == extracted_hash,
                }
            )
    return rows


def main() -> None:
    files = sorted(MATRIX_DIR.glob("*.csv.gz"))
    if len(files) != 8:
        raise RuntimeError(f"Expected eight matrices, observed {len(files)}")
    archive_rows = verify_archive_members(files)
    raw = ad.read_h5ad(RAW, backed="r")
    raw_shape = tuple(raw.shape)
    if raw.shape != (64993, 33694):
        raise RuntimeError(f"Unexpected raw h5ad shape: {raw.shape}")

    per_sample: list[dict[str, object]] = []
    global_matrix_mismatches = 0
    global_layer_mismatches = 0
    global_max_difference = 0.0
    global_source_total = 0
    global_raw_total = 0
    global_layer_total = 0
    all_genes_exact = True
    all_cells_exact = True
    raw_integer_nonnegative = True

    for path in files:
        accession = path.name.split("_", 1)[0]
        cell_mask = raw.obs["sample_accession"].astype(str).eq(accession).to_numpy()
        cell_indices = np.flatnonzero(cell_mask)
        if not len(cell_indices) or not np.array_equal(
            cell_indices, np.arange(cell_indices[0], cell_indices[-1] + 1)
        ):
            raise RuntimeError(f"Cells for {accession} are absent or non-contiguous")
        gene_offset = 0
        sample_matrix_mismatches = 0
        sample_layer_mismatches = 0
        sample_max_difference = 0.0
        sample_source_total = 0
        row_slice = slice(cell_indices[0], cell_indices[-1] + 1)
        observed_full = raw.X[row_slice, :]
        layer_full = raw.layers["counts"][row_slice, :]
        if not sparse.issparse(observed_full):
            observed_full = sparse.csr_matrix(observed_full)
        if not sparse.issparse(layer_full):
            layer_full = sparse.csr_matrix(layer_full)
        observed_full = observed_full.T.tocsr()
        layer_full = layer_full.T.tocsr()
        sample_raw_total = int(observed_full.sum())
        sample_layer_total = int(layer_full.sum())
        raw_integer_nonnegative &= bool(
            np.all(observed_full.data >= 0)
            and np.all(observed_full.data == np.floor(observed_full.data))
            and np.all(layer_full.data >= 0)
            and np.all(layer_full.data == np.floor(layer_full.data))
        )
        sample_genes_exact = True
        with gzip.open(path, "rb") as handle:
            header = next(csv.reader([handle.readline().decode("utf-8").rstrip("\r\n")]))
            source_cells = [str(value) for value in header[1:]]
            raw_cells = raw.obs_names[cell_indices].astype(str).tolist()
            cells_exact = source_cells == raw_cells
            all_cells_exact &= cells_exact
            for line in handle:
                comma = line.find(b",")
                if comma < 1:
                    raise RuntimeError(f"Malformed row {gene_offset + 2} in {path.name}")
                source_gene = line[:comma].decode("utf-8").strip('"')
                expected_gene = str(raw.var_names[gene_offset])
                genes_exact = source_gene == expected_gene
                sample_genes_exact &= genes_exact
                all_genes_exact &= genes_exact
                values = np.fromstring(line[comma + 1 :], dtype=np.int64, sep=",")
                if len(values) != len(source_cells):
                    raise RuntimeError(
                        f"{path.name} row {gene_offset + 2} has {len(values)} values; "
                        f"expected {len(source_cells)}"
                    )
                observed = observed_full.getrow(gene_offset)
                source_nonzero = np.flatnonzero(values)
                row_exact = np.array_equal(source_nonzero, observed.indices) and np.array_equal(
                    values[source_nonzero], observed.data
                )
                if not row_exact:
                    observed_dense = np.zeros(len(values), dtype=np.float64)
                    observed_dense[observed.indices] = observed.data
                    difference = values.astype(np.float64) - observed_dense
                    sample_matrix_mismatches += int(np.count_nonzero(difference))
                    sample_max_difference = max(
                        sample_max_difference, float(np.max(np.abs(difference)))
                    )
                sample_source_total += int(values.sum())
                gene_offset += 1
        if (
            np.array_equal(observed_full.indptr, layer_full.indptr)
            and np.array_equal(observed_full.indices, layer_full.indices)
            and np.array_equal(observed_full.data, layer_full.data)
        ):
            sample_layer_mismatches = 0
        else:
            for start in range(0, raw.n_vars, 1024):
                stop = min(start + 1024, raw.n_vars)
                layer_difference = (
                    observed_full[start:stop, :].astype(np.float64)
                    - layer_full[start:stop, :].astype(np.float64)
                )
                sample_layer_mismatches += int(layer_difference.count_nonzero())
        if gene_offset != raw.n_vars:
            raise RuntimeError(f"{path.name} contained {gene_offset} genes, expected {raw.n_vars}")
        global_matrix_mismatches += sample_matrix_mismatches
        global_layer_mismatches += sample_layer_mismatches
        global_max_difference = max(global_max_difference, sample_max_difference)
        global_source_total += sample_source_total
        global_raw_total += sample_raw_total
        global_layer_total += sample_layer_total
        per_sample.append(
            {
                "sample_accession": accession,
                "source_file": path.name,
                "cells": len(source_cells),
                "genes": gene_offset,
                "cell_order_exact": cells_exact,
                "gene_order_exact": sample_genes_exact,
                "source_vs_raw_mismatches": sample_matrix_mismatches,
                "raw_X_vs_counts_layer_mismatches": sample_layer_mismatches,
                "max_absolute_difference": sample_max_difference,
                "source_total_umi": sample_source_total,
                "raw_total_umi": sample_raw_total,
                "counts_layer_total_umi": sample_layer_total,
            }
        )
        print(
            json.dumps(
                {
                    "sample_accession": accession,
                    "cells": len(source_cells),
                    "source_vs_raw_mismatches": sample_matrix_mismatches,
                    "raw_X_vs_counts_layer_mismatches": sample_layer_mismatches,
                }
            ),
            flush=True,
        )
        del observed_full, layer_full
        gc.collect()
    raw.file.close()

    checks = {
        "archive_has_eight_matching_members": len(archive_rows) == 8
        and all(row["byte_and_hash_match"] for row in archive_rows),
        "raw_shape_exact": raw_shape == (64993, 33694),
        "eight_accessions_unique": len({row["sample_accession"] for row in per_sample}) == 8,
        "cell_order_exact_all_samples": bool(all_cells_exact),
        "gene_order_exact_all_samples": bool(all_genes_exact),
        "matrix_counts_exact": global_matrix_mismatches == 0
        and global_max_difference == 0,
        "raw_X_and_counts_layer_exact": global_layer_mismatches == 0,
        "counts_nonnegative_integers": bool(raw_integer_nonnegative),
        "total_umi_conserved": global_source_total == global_raw_total == global_layer_total,
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "gate": "03_matrix_to_raw_semantic_lineage",
        "evaluated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks": checks,
        "archive_members": archive_rows,
        "per_sample": per_sample,
        "totals": {
            "cells": sum(row["cells"] for row in per_sample),
            "genes_per_matrix": raw_shape[1],
            "source_vs_raw_mismatches": global_matrix_mismatches,
            "raw_X_vs_counts_layer_mismatches": global_layer_mismatches,
            "max_absolute_difference": global_max_difference,
            "source_total_umi": global_source_total,
            "raw_total_umi": global_raw_total,
            "counts_layer_total_umi": global_layer_total,
        },
        "evidence_boundary": (
            "This proves byte identity of the deposited archive members and value-level "
            "identity of the eight author-processed matrices with raw h5ad. The deposited "
            "matrices are not FASTQ-level raw data, and this check does not replay Scanpy."
        ),
    }
    VALIDATION.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    write_numbered_attempt(
        VALIDATION,
        "gate_03_matrix_to_raw_validation_attempt",
        report["status"],
        payload,
    )
    atomic_replace_json(
        VALIDATION / "gate_03_matrix_to_raw_validation_latest.json",
        payload,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["status"] == "PASS" else 2)


if __name__ == "__main__":
    main()
