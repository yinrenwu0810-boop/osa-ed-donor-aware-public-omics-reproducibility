"""Create an additive, hash-bound erratum for the VP-G04 label defect.

Historical run artifacts are intentionally never modified.  This document
binds each retained historical run_parameters.json hash to its corrected G04
metadata and records the pre/post source-script hashes.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad")
VP = ROOT / "revision_v2" / "07_virtual_perturbation"
PRE_CORRECTION_RUNNER_SHA256 = "16caccfbf1c5ac2314d092639d4e78bbafe0bc6ef50ef26ea664286aa90a9632"
RAW_GATE = "VP-G03"
RAW_PURPOSE = "Technical pilot only; biological interpretation prohibited."
CORRECTED_GATE = "VP-G04"
CORRECTED_PURPOSE = "Main frozen virtual-KO; biological interpretation prohibited."


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def main() -> None:
    registry_path = VP / "05_runs" / "run_registry.tsv"
    runner_path = VP / "scripts" / "02_run_sctenifoldknk.R"
    rows = read_csv(registry_path)
    if len(rows) != 45 or {row["status"] for row in rows} != {"PASS_TECHNICAL"}:
        raise RuntimeError("Erratum requires exactly 45 retained PASS_TECHNICAL G04 rows.")
    entries = []
    for row in rows:
        parameters_path = Path(row["output_attempt"]) / "run_parameters.json"
        parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
        if parameters.get("gate") != RAW_GATE or parameters.get("purpose") != RAW_PURPOSE:
            raise RuntimeError(f"Unexpected raw metadata for {row['run_id']}; refusing generic correction.")
        entries.append(
            {
                "run_order": int(row["run_order"]),
                "run_id": row["run_id"],
                "historical_run_parameters_sha256": sha256(parameters_path),
                "raw_metadata": {"gate": RAW_GATE, "purpose": RAW_PURPOSE},
                "corrected_metadata": {"gate": CORRECTED_GATE, "purpose": CORRECTED_PURPOSE},
            }
        )
    report = {
        "document": "VP-G04 additive metadata erratum v1",
        "status": "PASS_ADDITIVE_METADATA_ERRATUM",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "45 frozen VP-G04 main virtual-KO tasks only",
        "policy": "Historical result artifacts and their manifests remain unmodified. This erratum is an additive hash-bound provenance layer.",
        "defect": "The historical runner hard-coded VP-G03 / Technical pilot metadata while executing frozen VP-G04 main tasks.",
        "runner_source_sha256_before_correction": PRE_CORRECTION_RUNNER_SHA256,
        "runner_source_sha256_after_correction": sha256(runner_path),
        "registry_sha256": sha256(registry_path),
        "entries": entries,
    }
    destination = VP / "validation" / "VP_G04_metadata_erratum_v1.json"
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "entries": len(entries), "output": str(destination)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
