from __future__ import annotations

import csv
import hashlib
from pathlib import Path


OUTPUT_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = OUTPUT_ROOT / "figures" / "Figure3_20260912"
MANIFEST = FIGURE_DIR / "Figure3_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "Figure3_ED_corpus_cavernosum_localization.svg",
    "Figure3_ED_corpus_cavernosum_localization.pdf",
    "Figure3_ED_corpus_cavernosum_localization.png",
    "Figure3_ED_corpus_cavernosum_localization.tiff",
    "Figure3_build_audit.json",
    "Figure3_QA.md",
    "Figure3_legend.md",
    "source_data/Figure3a_UMAP_descriptive_source_data.tsv.gz",
    "source_data/Figure3b_donor_celltype_composition_source_data.tsv",
    "source_data/Figure3c_donor_pseudobulk_pathway_source_data.tsv",
    "source_data/Figure3d_gene_celltype_effect_source_data.tsv",
    "source_data/Figure3e_fibroblast_donor_dot_source_data.tsv",
]
UPSTREAM_PATHS = [
    "figure_contracts/Figure3_contract.md",
    "figure_contracts/input_manifest.sha256.tsv",
    "scripts/build_figure3.py",
    "scripts/qa_figure3.py",
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
