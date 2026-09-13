#!/usr/bin/env python3
"""Verify that the v1 inputs frozen for revision v2 have not changed."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REVISION = ROOT / "revision_v2"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest_path = REVISION / "00_protocol/input_manifest.sha256.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    rows = []
    for item in manifest.itertuples(index=False):
        path = ROOT / item.relative_path
        exists = path.is_file()
        current_bytes = path.stat().st_size if exists else None
        current_sha256 = sha256(path) if exists else None
        rows.append(
            {
                "relative_path": item.relative_path,
                "exists": exists,
                "bytes_match": exists and current_bytes == item.bytes,
                "sha256_match": exists and current_sha256 == item.sha256,
            }
        )
    status = pd.DataFrame(rows)
    report = {
        "status": "PASS" if status[["exists", "bytes_match", "sha256_match"]].all().all() else "FAIL",
        "manifest_entries": len(status),
        "missing": int((~status["exists"]).sum()),
        "byte_mismatch": int((~status["bytes_match"]).sum()),
        "hash_mismatch": int((~status["sha256_match"]).sum()),
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "PASS":
        print(status.loc[~status["sha256_match"]].to_string(index=False))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
