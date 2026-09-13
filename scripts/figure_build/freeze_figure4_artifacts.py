from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure4_20260912"
MANIFEST = FIGURE_DIR / "Figure4_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure4_candidate_hierarchy_and_doublet_sensitivity.svg",
    "Figure4_candidate_hierarchy_and_doublet_sensitivity.pdf",
    "Figure4_candidate_hierarchy_and_doublet_sensitivity.png",
    "Figure4_candidate_hierarchy_and_doublet_sensitivity.tiff",
    "Figure4_build_audit.json",
    "Figure4_QA.md",
    "Figure4_legend.md",
    "source_data/Figure4a_candidate_hierarchy_source_data.tsv",
    "source_data/Figure4b_candidate_evidence_matrix_source_data.tsv",
    "source_data/Figure4c_doublet_gene_sensitivity_source_data.tsv",
    "source_data/Figure4d_doublet_pathway_sensitivity_source_data.tsv",
    "source_data/Figure4e_pathway_state_count_source_data.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure4_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure4.py",
    "scripts/qa_figure4.py",
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
