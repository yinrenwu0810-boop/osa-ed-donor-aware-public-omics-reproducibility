from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure8_20260912"
MANIFEST = FIGURE_DIR / "Figure8_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure8_statistical_resolution_and_external_gap.svg",
    "Figure8_statistical_resolution_and_external_gap.pdf",
    "Figure8_statistical_resolution_and_external_gap.png",
    "Figure8_statistical_resolution_and_external_gap.tiff",
    "Figure8_build_audit.json",
    "Figure8_session_record.json",
    "Figure8_QA.md",
    "Figure8_legend.md",
    "source_data/Figure8a_effect_distribution_source_data.tsv",
    "source_data/Figure8a_top15_absolute_delta_source_data.tsv",
    "source_data/Figure8b_exact_allocation_schematic_source_data.tsv",
    "source_data/Figure8c_n3v3_MDE_source_data.tsv",
    "source_data/Figure8c_idealized_sample_size_source_data.tsv",
    "source_data/Figure8d_external_coverage_records_source_data.tsv",
    "source_data/Figure8d_coverage_matrix_source_data.tsv",
    "source_data/Figure8e_gap_map_source_data.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure8_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure8.py",
    "scripts/qa_figure8.py",
]


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


rows = []
for relative_path in RELATIVE_PATHS:
    path = FIGURE_DIR / relative_path
    if not path.is_file():
        raise FileNotFoundError(path)
    rows.append({"scope": "output", "relative_path": relative_path, "size_bytes": path.stat().st_size, "sha256": digest(path)})
for relative_path in UPSTREAM_PATHS:
    path = OUTPUT_ROOT / relative_path
    if not path.is_file():
        raise FileNotFoundError(path)
    rows.append({"scope": "upstream", "relative_path": relative_path, "size_bytes": path.stat().st_size, "sha256": digest(path)})

with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["scope", "relative_path", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

print(f"artifact_count={len(rows)}")
print(f"manifest={MANIFEST}")
