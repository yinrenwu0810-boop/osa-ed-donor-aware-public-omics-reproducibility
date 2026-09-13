#!/usr/bin/env python3
"""Audit v4: the versioned VP manuscript with Figure 5 physically embedded."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from docx import Document


HERE = Path(__file__).resolve().parent
REVISION = HERE.parents[1]
VP = REVISION / "07_virtual_perturbation"
V2 = REVISION / "06_delivery" / "IH_OSA_ED_中文修订工作稿_v2.docx"
V3 = REVISION / "06_delivery" / "IH_OSA_ED_中文修订工作稿_v3_含虚拟扰动结果.docx"
V4 = REVISION / "06_delivery" / "IH_OSA_ED_中文修订工作稿_v4_含虚拟扰动结果与图5.docx"
SOURCE = HERE / "manuscript_zh.v3_virtual_perturbation.src.md"
FIGURE = HERE / "figures" / "Figure5_virtual_perturbation_hypothesis_prioritization.png"
AUDIT = HERE / "VP_manuscript_zh_v4_AUDIT03.json"
GATES = [
    VP / "validation" / "VP_G04_main_runs_COMPLETE.json",
    VP / "validation" / "VP_G05_stability_COMPLETE.json",
    VP / "validation" / "VP_G06_enrichment_COMPLETE.json",
    VP / "validation" / "VP_G07_external_scope_COMPLETE.json",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if AUDIT.exists():
        raise SystemExit(f"Refusing to overwrite existing audit: {AUDIT}")
    if not all(path.is_file() for path in (V2, V3, V4, SOURCE, FIGURE, *GATES)):
        raise SystemExit("A required sealed input is missing")

    source = SOURCE.read_text(encoding="utf-8")
    doc = Document(V4)
    paragraphs = doc.paragraphs
    doc_text = "\n".join(p.text for p in paragraphs)
    caption_index = next((i for i, p in enumerate(paragraphs) if p.text.startswith("图5｜虚拟扰动的跨供体稳定性")), None)
    gate_statuses = [json.loads(path.read_text(encoding="utf-8")).get("status") for path in GATES]
    checks = [
        {"name": "sealed_g04_to_g07_pass", "pass": gate_statuses == ["PASS", "PASS", "PASS", "PASS"]},
        {"name": "v2_and_v3_preserved_as_prior_versions", "pass": len({sha256(V2), sha256(V3), sha256(V4)}) == 3},
        {"name": "vp_results_methods_and_boundaries_retained", "pass": all(token in source for token in ["## 2.8. 冻结的虚拟扰动分析", "## 4.10. 虚拟扰动、稳定性评估和外部覆盖审计", "不构成外部生物学验证", "G08不属于本研究的已获得数据"])},
        {"name": "figure5_is_physically_embedded", "pass": doc.inline_shapes.__len__() >= 1 and caption_index is not None},
        {"name": "figure5_precedes_its_caption", "pass": caption_index is not None and any(paragraph._p.xpath('.//a:blip') for paragraph in paragraphs[max(0, caption_index - 1):caption_index])},
        {"name": "reader_facing_limitations_retained", "pass": all(token in doc_text for token in ["未获得Hallmark Hypoxia的FDR显著网络富集", "ED下调FDR基因集仅含3个基因", "没有识别到人海绵体成纤维细胞中三个靶点的真实遗传扰动T1外部参照", "尚无细胞实验数据"])},
    ]
    passed = all(item["pass"] for item in checks)
    payload = {
        "artifact": V4.name,
        "status": "PASS_VERSIONED_MANUSCRIPT_INTEGRATION_WITH_EMBEDDED_FIGURE" if passed else "FAIL_MANUSCRIPT_INTEGRATION_AUDIT",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "gate_statuses": gate_statuses,
        "v2_docx_sha256": sha256(V2),
        "v3_docx_sha256": sha256(V3),
        "v4_docx_sha256": sha256(V4),
        "source_sha256": sha256(SOURCE),
        "figure_sha256": sha256(FIGURE),
        "boundary": "Only completed VP-G04--G07 computational evidence is added. G08 remains a prospective feasibility design without wet-lab results.",
    }
    AUDIT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "audit": str(AUDIT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
