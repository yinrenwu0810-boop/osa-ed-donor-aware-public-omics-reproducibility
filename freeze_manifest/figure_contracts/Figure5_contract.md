# Figure 5 Contract — Partial transportability and marked heterogeneity

**Core conclusion:** IH-related signals show partial leave-one-dataset-out transportability, but dataset-level AUC and selected features are heterogeneous and do not constitute diagnostic validation.

- Figure archetype: quantitative grid.
- Target/output: Scientific Reports main Figure 5; Python; 183 × 135 mm; SVG/PDF/TIFF/PNG.
- Hero evidence: pooled prediction/ROC shown beside the 0.11–1.00 dataset AUC range.
- Status: `READY_FROM_SEALED_TABLES`; no new test.

## Panel map

| Panel | Unique question | Frozen display | Source IDs |
|---|---|---|---|
| a | What probabilities were assigned to held-out samples? | Jittered held-out probabilities grouped by dataset and true IH label. | S050 |
| b | How do the three frozen strategies perform in pooled predictions? | Pooled ROC curves and AUCs for gene-rank, Hallmark-rank, and adaptive strategies. | S050–S051 |
| c | Is performance consistent across held-out datasets? | Per-dataset AUC point/range display, prominently showing gene-rank AUC 0.11, 1.00, and 0.94. | S052 |
| d | Are selected features stable? | Outer-fold selection frequency and median absolute coefficient for a capped, deterministic feature set. | S053–S054 |
| e | How much uncertainty/calibration is visible? | Pooled AUC with conditional bootstrap interval and structure-preserving exact permutation P value; n=20 and 28,000 permutations shown. | S051 |

## Display rule

Panel d displays all features selected in at least two outer folds; if more than 20, retain the 20 highest selection counts then median absolute coefficient, with alphabetical ties. No model refitting is allowed.

## Integrity and reviewer risk

The terms “transportability” and “heterogeneity” are required. Avoid “diagnostic model”, “biomarker panel”, and “external clinical validation”. Confidence intervals and permutation P values are conditional on the assembled small datasets and do not remove dataset heterogeneity.

