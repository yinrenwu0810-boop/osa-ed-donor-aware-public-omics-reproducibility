#!/usr/bin/env python3
"""Audit the versioned Chinese-manuscript integration of VP-G04--G07 results."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from docx import Document


ACADEMIC_ROOT = Path(__file__).resolve().parents[3]
REVISION = ACADEMIC_ROOT / "revision_v2"
VP = REVISION / "07_virtual_perturbation"
HERE = Path(__file__).resolve().parent
V2 = REVISION / "06_delivery" / "IH_OSA_ED_中文修订工作稿_v2.docx"
V3 = REVISION / "06_delivery" / "IH_OSA_ED_中文修订工作稿_v3_含虚拟扰动结果.docx"
SOURCE = HERE / "manuscript_zh.v3_virtual_perturbation.src.md"
FIGURE = HERE / "figures" / "Figure5_virtual_perturbation_hypothesis_prioritization.png"
AUDIT = HERE / "VP_manuscript_zh_v3_AUDIT02.json"
G04 = VP / "validation" / "VP_G04_main_runs_COMPLETE.json"
G05 = VP / "validation" / "VP_G05_stability_COMPLETE.json"
G06 = VP / "validation" / "VP_G06_enrichment_COMPLETE.json"
G07 = VP / "validation" / "VP_G07_external_scope_COMPLETE.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if AUDIT.exists():
        raise SystemExit("Refusing to overwrite manuscript-integration audit")
    source = SOURCE.read_text(encoding="utf-8")
    doc_text = "\n".join(p.text for p in Document(V3).paragraphs)
    g04 = json.loads(G04.read_text(encoding="utf-8"))
    g05 = json.loads(G05.read_text(encoding="utf-8"))
    g06 = json.loads(G06.read_text(encoding="utf-8"))
    g07 = json.loads(G07.read_text(encoding="utf-8"))
    checks = [
        {"name": "sealed_gates_pass", "pass": g04.get("status") == "PASS" and g05.get("status") == "PASS" and g06.get("status") == "PASS" and g07.get("status") == "PASS"},
        {"name": "new_v3_does_not_overwrite_v2", "pass": V2.is_file() and V3.is_file() and sha256(V2) != sha256(V3)},
        {"name": "figure_and_results_section_present", "pass": FIGURE.is_file() and "## 2.8. 冻结的虚拟扰动分析" in source and "图5｜虚拟扰动" in source and "虚拟扰动的跨供体稳定性" in doc_text},
        {"name": "methods_added_and_software_renumbered", "pass": "## 4.10. 虚拟扰动、稳定性评估和外部覆盖审计" in source and "## 4.11. 软件与可重复性" in source},
        {"name": "no_overclaim_language", "pass": all(token in source for token in ["不构成外部生物学验证", "不表示KO后基因或通路", "G08不属于本研究的已获得数据"])},
        {"name": "limits_and_nonresults_retained", "pass": "没有识别到人海绵体成纤维细胞中三个靶点的真实遗传扰动T1外部参照" in source and "未获得Hallmark Hypoxia的FDR显著网络富集" in source and "ED下调FDR基因集仅含3个基因" in source},
    ]
    passed = all(item["pass"] for item in checks)
    payload = {
        "artifact": "IH_OSA_ED_中文修订工作稿_v3_含虚拟扰动结果.docx",
        "status": "PASS_VERSIONED_MANUSCRIPT_INTEGRATION" if passed else "FAIL_MANUSCRIPT_INTEGRATION_AUDIT",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "v2_docx_sha256": sha256(V2),
        "v3_docx_sha256": sha256(V3),
        "v3_source_sha256": sha256(SOURCE),
        "figure_sha256": sha256(FIGURE),
        "boundary": "The v3 manuscript adds only completed VP-G04--G07 computational results. VP-G08 is described as a prospective feasibility design with no wet-lab data, not as experimental validation."
    }
    AUDIT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "audit": str(AUDIT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
