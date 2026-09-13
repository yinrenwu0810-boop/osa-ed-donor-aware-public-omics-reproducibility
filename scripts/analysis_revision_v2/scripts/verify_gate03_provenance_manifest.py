#!/usr/bin/env python3
"""Read-only verifier for the Gate 03 provenance SHA256 manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = (
    WORKSPACE / "revision_v2" / "00_protocol" / "gate_03_provenance_manifest.sha256.tsv"
)
DEFAULT_VALIDATION_DIR = WORKSPACE / "revision_v2" / "04_validation"
REQUIRED_COLUMNS = {"path_kind", "path", "role", "bytes", "sha256"}
EXPECTED_MANIFEST_SHA256 = "c832788e13681be8c2d0d73e3d46cbbaa8e7b1eb4f579a1faaff70b90671043b"
EXPECTED_MANIFEST_ENTRIES = 31


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def resolve_manifest_path(path_kind: str, recorded_path: str) -> Path:
    if path_kind == "workspace":
        candidate = (WORKSPACE / Path(recorded_path)).resolve()
        try:
            candidate.relative_to(WORKSPACE.resolve())
        except ValueError as exc:
            raise ValueError(f"workspace path escapes workspace: {recorded_path}") from exc
        return candidate
    if path_kind == "absolute":
        candidate = Path(recorded_path)
        if not candidate.is_absolute():
            raise ValueError(f"absolute path is not absolute: {recorded_path}")
        return candidate
    raise ValueError(f"unsupported path_kind: {path_kind}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--write-validation", action="store_true")
    parser.add_argument("--validation-dir", type=Path, default=DEFAULT_VALIDATION_DIR)
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    manifest_sha256 = sha256(manifest_path)
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed_columns = reader.fieldnames
        entries = list(reader)

    required_columns_exact = (
        observed_columns is not None
        and len(observed_columns) == len(REQUIRED_COLUMNS)
        and set(observed_columns) == REQUIRED_COLUMNS
    )

    seen: set[tuple[str, str]] = set()
    details: list[dict[str, object]] = []
    if required_columns_exact:
        for entry in entries:
            identity = (entry["path_kind"], entry["path"])
            duplicate = identity in seen
            seen.add(identity)
            path_valid = True
            path_error = None
            try:
                target = resolve_manifest_path(*identity)
            except ValueError as exc:
                path_valid = False
                path_error = str(exc)
                target = None
            try:
                expected_bytes = int(entry["bytes"])
            except (TypeError, ValueError):
                expected_bytes = None
            expected_hash = entry["sha256"].lower()
            exists = bool(path_valid and target is not None and target.is_file())
            current_bytes = target.stat().st_size if exists and target is not None else None
            current_hash = sha256(target) if exists and target is not None else None
            details.append(
                {
                    "path_kind": entry["path_kind"],
                    "path": entry["path"],
                    "path_valid": path_valid,
                    "path_error": path_error,
                    "exists": exists,
                    "duplicate": duplicate,
                    "bytes_match": exists
                    and expected_bytes is not None
                    and current_bytes == expected_bytes,
                    "sha256_match": exists and current_hash == expected_hash,
                }
            )

    failures = [
        row
        for row in details
        if not (
            row["path_valid"]
            and row["exists"]
            and not row["duplicate"]
            and row["bytes_match"]
            and row["sha256_match"]
        )
    ]
    checks = {
        "manifest_identity_match": manifest_sha256 == EXPECTED_MANIFEST_SHA256,
        "manifest_entry_count_match": len(entries) == EXPECTED_MANIFEST_ENTRIES,
        "required_columns_exact": required_columns_exact,
        "no_path_escape": required_columns_exact
        and len(details) == len(entries)
        and all(row["path_valid"] for row in details),
        "all_targets_exist": required_columns_exact
        and len(details) == len(entries)
        and all(row["exists"] for row in details),
        "all_bytes_match": required_columns_exact
        and len(details) == len(entries)
        and all(row["bytes_match"] for row in details),
        "all_sha256_match": required_columns_exact
        and len(details) == len(entries)
        and all(row["sha256_match"] for row in details),
        "no_duplicates": required_columns_exact
        and len(details) == len(entries)
        and all(not row["duplicate"] for row in details),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "expected_manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "manifest_identity_match": checks["manifest_identity_match"],
        "manifest_entries": len(entries),
        "expected_manifest_entries": EXPECTED_MANIFEST_ENTRIES,
        "manifest_entry_count_match": checks["manifest_entry_count_match"],
        "observed_columns": observed_columns,
        "checks": checks,
        "workspace_entries": sum(row["path_kind"] == "workspace" for row in details),
        "absolute_entries": sum(row["path_kind"] == "absolute" for row in details),
        "path_errors": sum(not row["path_valid"] for row in details),
        "missing": sum(not row["exists"] for row in details),
        "duplicates": sum(row["duplicate"] for row in details),
        "byte_mismatch": sum(not row["bytes_match"] for row in details),
        "hash_mismatch": sum(not row["sha256_match"] for row in details),
        "failures": failures,
        "evidence_boundary": (
            "Confirms current bytes against this manifest only; it does not prove that the "
            "external Scanpy scripts are byte-identical to those used historically, recreate "
            "the historical software environment, or constitute an end-to-end replay."
        ),
    }
    if args.write_validation:
        args.validation_dir.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        write_numbered_attempt(
            args.validation_dir,
            "gate_03_provenance_manifest_validation_attempt",
            report["status"],
            payload,
        )
        atomic_replace_json(
            args.validation_dir / "gate_03_provenance_manifest_validation_latest.json",
            payload,
        )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
