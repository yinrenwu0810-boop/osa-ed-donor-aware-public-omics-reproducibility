from __future__ import annotations

import csv
import hashlib
from pathlib import Path


BASE = Path(__file__).resolve().parents[1]
OUTPUT = BASE / "W2_artifact_manifest.sha256.tsv"
RELATIVE_PATHS = [
    "paper_spine_config.json",
    "confirmed_motivation.md",
    "confirmed_contribution.md",
    "source_inventory.md",
    "evidence_bank.md",
    "figure_asset_map.md",
    "claim_register.md",
    "section_blueprints.md",
    "writing_rationale_matrix.md",
    "figure_contracts/FIGURE_CONTRACTS_MASTER.md",
    "figure_contracts/Figure1_contract.md",
    "figure_contracts/Figure2_contract.md",
    "figure_contracts/Figure3_contract.md",
    "figure_contracts/Figure4_contract.md",
    "figure_contracts/Figure5_contract.md",
    "figure_contracts/Figure6_contract.md",
    "figure_contracts/Figure7_contract.md",
    "figure_contracts/Figure8_contract.md",
    "figure_contracts/figure_source_data_inventory.tsv",
    "figure_contracts/input_manifest.sha256.tsv",
    "figure_contracts/SOURCE_FREEZE_STATUS.md",
    "W2_FIGURE_CONTRACTS_COMPLETE.md",
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


rows = []
for rel in RELATIVE_PATHS:
    path = BASE / rel
    if not path.is_file():
        raise FileNotFoundError(path)
    rows.append({"relative_path": rel, "size_bytes": path.stat().st_size, "sha256": digest(path)})

with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["relative_path", "size_bytes", "sha256"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"artifact_count={len(rows)}")
print(f"manifest={OUTPUT}")

