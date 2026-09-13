from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure7_20260912"
MANIFEST = FIGURE_DIR / "Figure7_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure7_exploratory_hypoxia_sensitivity.svg",
    "Figure7_exploratory_hypoxia_sensitivity.pdf",
    "Figure7_exploratory_hypoxia_sensitivity.png",
    "Figure7_exploratory_hypoxia_sensitivity.tiff",
    "Figure7_build_audit.json",
    "Figure7_QA.md",
    "Figure7_legend.md",
    "source_data/Figure7a_specification_NES_source_data.tsv",
    "source_data/Figure7b_inferential_FDR_source_data.tsv",
    "source_data/Figure7c_direction_stability_source_data.tsv",
    "source_data/Figure7d_matched_control_source_data.tsv",
    "source_data/Figure7e_leading_edge_stability_source_data.tsv",
    "source_data/Figure7f_mapping_source_data.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure7_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure7.py",
    "scripts/qa_figure7.py",
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
