from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure6_20260912"
MANIFEST = FIGURE_DIR / "Figure6_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure6_virtual_perturbation_stability_and_ED_association.svg",
    "Figure6_virtual_perturbation_stability_and_ED_association.pdf",
    "Figure6_virtual_perturbation_stability_and_ED_association.png",
    "Figure6_virtual_perturbation_stability_and_ED_association.tiff",
    "Figure6_build_audit.json",
    "Figure6_QA.md",
    "Figure6_legend.md",
    "source_data/Figure6a_design_source_data.tsv",
    "source_data/Figure6b_donor_stability_source_data.tsv",
    "source_data/Figure6c_ED_UP_GSEA_source_data.tsv",
    "source_data/Figure6d_top5_overlap_source_data.tsv",
    "source_data/Figure6e_matched_calibration_source_data.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure6_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure6.py",
    "scripts/qa_figure6.py",
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
