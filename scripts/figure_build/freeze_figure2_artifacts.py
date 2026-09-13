from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure2_20260911"
MANIFEST = FIGURE_DIR / "Figure2_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure2_endothelial_convergence_and_heterogeneity.svg",
    "Figure2_endothelial_convergence_and_heterogeneity.pdf",
    "Figure2_endothelial_convergence_and_heterogeneity.png",
    "Figure2_endothelial_convergence_and_heterogeneity.tiff",
    "Figure2_build_audit.json",
    "Figure2_QA.md",
    "Figure2_legend.md",
    "source_data/Figure2a_volcano_source_data.tsv",
    "source_data/Figure2b_reported_correlation.tsv",
    "source_data/Figure2b_shared_effect_source_data.tsv",
    "source_data/Figure2c_filter_count_source_data.tsv",
    "source_data/Figure2d_strict_lodo_heatmap_source_data.tsv",
    "source_data/Figure2e_pathway_dotplot_source_data.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure2_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure2.py",
    "scripts/qa_figure2.py",
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
