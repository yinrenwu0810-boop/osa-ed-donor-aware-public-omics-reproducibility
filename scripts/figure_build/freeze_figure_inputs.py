from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = OUTPUT_DIR.parent
INVENTORY = OUTPUT_DIR / "figure_contracts" / "figure_source_data_inventory.tsv"
MANIFEST = OUTPUT_DIR / "figure_contracts" / "input_manifest.sha256.tsv"
STATUS = OUTPUT_DIR / "figure_contracts" / "SOURCE_FREEZE_STATUS.md"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


rows: list[dict[str, str]] = []
with INVENTORY.open("r", encoding="utf-8-sig", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        path = PROJECT_ROOT / row["relative_path"]
        exists = path.is_file()
        rows.append(
            {
                **row,
                "exists": str(exists).upper(),
                "size_bytes": str(path.stat().st_size) if exists else "",
                "sha256": sha256_file(path) if exists else "",
            }
        )

fieldnames = list(rows[0].keys())
with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

missing = [row for row in rows if row["exists"] != "TRUE"]
status = "PASS" if not missing else "HOLD_MISSING_INPUT"
timestamp = datetime.now(timezone.utc).astimezone().isoformat()
lines = [
    "# Figure source freeze status",
    "",
    f"- Generated: `{timestamp}`",
    f"- Inventory rows: **{len(rows)}**",
    f"- Missing inputs: **{len(missing)}**",
    f"- Status: **{status}**",
    "",
    "This gate verifies only file presence and SHA-256 identity. It does not upgrade scientific evidence or authorize new analysis.",
]
if missing:
    lines.extend(["", "## Missing", ""])
    lines.extend(f"- `{row['source_id']}`: `{row['relative_path']}`" for row in missing)
STATUS.write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"status={status}")
print(f"inventory_rows={len(rows)}")
print(f"missing_inputs={len(missing)}")
print(f"manifest={MANIFEST}")

