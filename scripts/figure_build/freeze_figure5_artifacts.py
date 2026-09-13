from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure5_20260912"
MANIFEST = FIGURE_DIR / "Figure5_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure5_partial_transportability_and_heterogeneity.svg",
    "Figure5_partial_transportability_and_heterogeneity.pdf",
    "Figure5_partial_transportability_and_heterogeneity.png",
    "Figure5_partial_transportability_and_heterogeneity.tiff",
    "Figure5_build_audit.json",
    "Figure5_QA.md",
    "Figure5_legend.md",
    "source_data/Figure5a_heldout_probability_source_data.tsv",
    "source_data/Figure5b_pooled_ROC_source_data.tsv",
    "source_data/Figure5c_per_dataset_AUC_source_data.tsv",
    "source_data/Figure5d_feature_stability_source_data.tsv",
    "source_data/Figure5e_pooled_metric_source_data.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure5_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure5.py",
    "scripts/qa_figure5.py",
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


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
