from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure1_20260911"
MANIFEST = FIGURE_DIR / "Figure1_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure1_parallel_evidence_and_strict_bridge.svg",
    "Figure1_parallel_evidence_and_strict_bridge.pdf",
    "Figure1_parallel_evidence_and_strict_bridge.png",
    "Figure1_parallel_evidence_and_strict_bridge.tiff",
    "Figure1_build_audit.json",
    "Figure1_QA.md",
    "Figure1_legend.md",
    "source_data/Figure1a_evidence_layers.tsv",
    "source_data/Figure1b_endothelial_filters.tsv",
    "source_data/Figure1c_clinical_bridge_filters.tsv",
    "source_data/Figure1d_strict_context_heatmap.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure1_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure1.py",
    "scripts/qa_figure1.py",
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

