#!/usr/bin/env python3
"""Fail-closed preflight for the portable GSE206528 Gate 03 package."""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "revision_v2/00_protocol/gate_03_provenance_manifest.sha256.tsv"
SNAPSHOT = ROOT / "tools/scanpy_snapshot"
VALIDATION = ROOT / "revision_v2/04_validation"
PACKAGE_MANIFEST = ROOT / "TRANSFER_PACKAGE_MANIFEST.sha256.tsv"
EXPECTED_MANIFEST_SHA256 = "c832788e13681be8c2d0d73e3d46cbbaa8e7b1eb4f579a1faaff70b90671043b"
EXPECTED_PACKAGES = {
    "scanpy": "1.12.1",
    "anndata": "0.13.1",
    "harmonypy": "0.0.10",
    "igraph": "1.0.0",
    "leidenalg": "0.12.0",
    "numpy": "2.4.6",
    "pandas": "3.0.3",
    "scipy": "1.18.0",
    "scikit-learn": "1.9.0",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def available_memory_gb() -> float | None:
    if sys.platform != "win32":
        return None

    class MemoryStatusEx(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    status = MemoryStatusEx()
    status.dwLength = ctypes.sizeof(status)
    if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        return status.ullAvailPhys / 1024**3
    return None


def package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def resolve_entry(path_kind: str, recorded_path: str) -> Path:
    if path_kind == "workspace":
        candidate = (ROOT / recorded_path).resolve()
        try:
            candidate.relative_to(ROOT)
        except ValueError as exc:
            raise RuntimeError(f"Workspace path escapes package: {recorded_path}") from exc
        return candidate
    if path_kind == "absolute":
        return (SNAPSHOT / Path(recorded_path).name).resolve()
    raise RuntimeError(f"Unsupported path kind: {path_kind}")


def verify_package_manifest() -> bool:
    if not PACKAGE_MANIFEST.is_file():
        return False
    with PACKAGE_MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    required_columns = {"relative_path", "bytes", "sha256"}
    if not entries or set(entries[0]) != required_columns:
        return False
    for entry in entries:
        target = (ROOT / entry["relative_path"]).resolve()
        try:
            target.relative_to(ROOT)
        except ValueError:
            return False
        if (
            not target.is_file()
            or target.stat().st_size != int(entry["bytes"])
            or sha256(target) != entry["sha256"]
        ):
            return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()

    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    artifact_failures = []
    for entry in entries:
        target = resolve_entry(entry["path_kind"], entry["path"])
        current_hash = sha256(target) if target.is_file() else None
        current_size = target.stat().st_size if target.is_file() else None
        if current_hash != entry["sha256"] or current_size != int(entry["bytes"]):
            artifact_failures.append(entry["path"])

    memory = available_memory_gb()
    disk = shutil.disk_usage(ROOT).free / 1024**3
    packages = {name: package_version(name) for name in EXPECTED_PACKAGES}
    checks = {
        "package_root_has_expected_name": ROOT.name == "acad",
        "manifest_identity_exact": sha256(MANIFEST) == EXPECTED_MANIFEST_SHA256,
        "manifest_entry_count_exact": len(entries) == 31,
        "transfer_package_manifest_exact": verify_package_manifest(),
        "all_31_frozen_inputs_exact": not artifact_failures,
        "scanpy_snapshot_present": SNAPSHOT.is_dir() and len(list(SNAPSHOT.glob("*.py"))) == 4,
        "python_is_3_12": sys.version_info[:2] == (3, 12),
        "packages_exact": packages == EXPECTED_PACKAGES,
        "available_memory_at_least_10_gb": memory is not None and memory >= 10.0,
        "available_disk_at_least_40_gb": disk >= 40.0,
    }
    report = {
        "status": "READY" if all(checks.values()) else "HOLD",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "available_memory_gb": memory,
        "available_disk_gb": disk,
        "packages": packages,
        "artifact_failures": artifact_failures,
        "evidence_boundary": (
            "READY only authorizes the technical replay attempt. It does not make "
            "Gate 03 PASS, does not complete Gate 03B doublet sensitivity, and "
            "does not establish independent biological validation."
        ),
    }
    if args.write_report:
        VALIDATION.mkdir(parents=True, exist_ok=True)
        path = VALIDATION / "target_machine_preflight_latest.json"
        path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    raise SystemExit(0 if report["status"] == "READY" else 2)


if __name__ == "__main__":
    main()
