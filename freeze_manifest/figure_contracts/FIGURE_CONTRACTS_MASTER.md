# Figure 1–8 Contract Index

_Frozen after confirmation of Option A | Target: Scientific Reports | Backend: Python only_

## Common export contract

- Eight main figures and no main-text tables; the eight-display-item limit is therefore preserved.
- Full-width figures use 183 mm; height is figure-specific and must remain legible on an A4/Letter manuscript page.
- Required deliverables per figure: editable SVG, vector PDF, 600-dpi TIFF, review PNG, one Python build script, panel-level source data, dependency/session record, and QA report.
- Statistical unit must be stated in every legend. Donors or study datasets are biological units where applicable; cells, random seeds, leave-one-out recomputations, pathways, and matched genes are not biological replicates.
- A plot may mechanically filter, reshape, label, or summarize a sealed table only under the rule frozen in its contract. Any new statistical test, changed FDR family, changed candidate set, or changed matching rule triggers `HOLD_NEW_ANALYSIS_REQUIRED`.

## Readiness

| Figure | Contract status | New inferential analysis? | Main reason |
|---|---|---:|---|
| 1 | READY_AFTER_MECHANICAL_DERIVATION | No | Existing bridge tables fully support the claim. |
| 2 | READY_AFTER_RULED_DERIVATION | No | Representative pathway display needs a deterministic selection rule; no statistics are recomputed. |
| 3 | READY_AFTER_RULED_DERIVATION | No | Donor composition and candidate/pathway panels exist; optional donor dots require descriptive normalization only. |
| 4 | READY_FROM_SEALED_TABLES | No | Candidate class and doublet sensitivity tables are complete. |
| 5 | READY_FROM_SEALED_TABLES | No | Predictions, per-dataset metrics, intervals, permutation results, and feature stability are present. |
| 6 | READY_FROM_SEALED_TABLES | No | Stability, ED enrichment, overlap, eligibility, and controls are present. |
| 7 | READY_WITH_PANEL_B_REVISION | No | Replace the unavailable all-32-set S1/S2 panel with the available Hypoxia-only S0–S2 FDR panel. |
| 8 | READY_AFTER_MECHANICAL_DERIVATION | No | Exact-P resolution, power, and external coverage sources are present. |

## Gate result

All 54 registered inputs were present when the freeze script ran. The SHA-256 manifest is the authority for later drift checks. No main figure is currently blocked by a missing source, provided the Figure 7 panel revision and the non-inferential derivation rules below are retained.

