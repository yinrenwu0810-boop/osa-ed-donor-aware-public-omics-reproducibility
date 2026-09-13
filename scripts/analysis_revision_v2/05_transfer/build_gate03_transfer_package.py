#!/usr/bin/env python3
"""Create a self-contained, hash-verifiable Gate 03 transfer package."""

from __future__ import annotations

import csv
import hashlib
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "revision_v2/00_protocol/gate_03_provenance_manifest.sha256.tsv"
PACKAGE_PARENT = ROOT.parent / "transfer_packages"
PACKAGE_NAME = "acad_gate03_handoff_2026-08-13"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def ignored(_: str, names: list[str]) -> set[str]:
    return {name for name in names if name == "__pycache__" or name.endswith(".pyc")}


def write_package_manifest(package_root: Path) -> Path:
    output = package_root / "TRANSFER_PACKAGE_MANIFEST.sha256.tsv"
    rows = []
    for path in sorted(package_root.rglob("*")):
        if not path.is_file() or path == output:
            continue
        rows.append(
            (path.relative_to(package_root).as_posix(), path.stat().st_size, sha256(path))
        )
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write("relative_path\tbytes\tsha256\n")
        for relative, size, digest in rows:
            handle.write(f"{relative}\t{size}\t{digest}\n")
    return output


def write_archive_sidecar(archive: Path) -> Path:
    sidecar = archive.with_suffix(archive.suffix + ".sha256.txt")
    text = (
        "Archive created at: "
        + datetime.now(timezone.utc).isoformat()
        + "\nArchive: "
        + str(archive)
        + "\nArchive SHA-256: "
        + sha256(archive)
        + "\nArchive bytes: "
        + str(archive.stat().st_size)
        + "\n"
    )
    sidecar.write_text(text, encoding="utf-8")
    return sidecar


def main() -> None:
    package_parent = PACKAGE_PARENT
    package_root = package_parent / PACKAGE_NAME / "acad"
    archive = package_parent / f"{PACKAGE_NAME}.zip"
    if package_root.exists() or archive.exists():
        raise RuntimeError(
            f"Refusing to overwrite existing transfer artifact: {package_root} or {archive}"
        )
    package_parent.mkdir(parents=True, exist_ok=True)
    package_root.mkdir(parents=True)
    try:
        shutil.copytree(ROOT / "revision_v2", package_root / "revision_v2", ignore=ignored)
        with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
            entries = list(csv.DictReader(handle, delimiter="\t"))
        for entry in entries:
            if entry["path_kind"] == "workspace":
                source = ROOT / entry["path"]
                target = package_root / entry["path"]
            elif entry["path_kind"] == "absolute":
                source = Path(entry["path"])
                target = package_root / "tools/scanpy_snapshot" / source.name
            else:
                raise RuntimeError(f"Unsupported path kind: {entry['path_kind']}")
            if not source.is_file():
                raise RuntimeError(f"Frozen source is absent: {source}")
            copy_file(source, target)
        manifest_path = write_package_manifest(package_root)
        with zipfile.ZipFile(
            archive, "w", compression=zipfile.ZIP_STORED, allowZip64=True
        ) as handle:
            for path in sorted(package_root.rglob("*")):
                if path.is_file():
                    handle.write(path, path.relative_to(package_root.parent))
        sidecar = write_archive_sidecar(archive)
        print(
            {
                "package_root": str(package_root),
                "archive": str(archive),
                "archive_bytes": archive.stat().st_size,
                "archive_sha256": sha256(archive),
                "archive_sha256_sidecar": str(sidecar),
                "package_manifest": str(manifest_path),
            }
        )
    except BaseException:
        raise


if __name__ == "__main__":
    main()
