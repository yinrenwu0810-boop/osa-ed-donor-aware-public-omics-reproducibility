# Figure 6 Contract — Virtual-perturbation stability and ED signature association

**Core conclusion:** Predeclared virtual perturbations of TYMS, EFNB2, and LRRC17 yield stable network rankings and ED_UP signature associations in the tested model, but these are predictive network effects rather than biological knockout validation or expression reversal.

- Figure archetype: schematic-led composite.
- Target/output: Scientific Reports main Figure 6; Python; 183 × 145 mm; SVG/PDF/TIFF/PNG.
- Hero evidence: ED_UP ranked enrichment and Top-5-percent overlap, supported by stability and matched controls.
- Status: `READY_FROM_SEALED_TABLES`; no new test.

## Panel map

| Panel | Unique question | Frozen display | Source IDs |
|---|---|---|---|
| a | What was the repeated model design? | Three targets × three donor backgrounds × five seeds = 45 main runs, with 90 matched controls; seeds shown as technical/model repeats. | S003, S060–S066 |
| b | Are rankings stable across seeds and donors? | Donor-level median seed Spearman and Top-200 Jaccard, with donor-consensus summary. | S060–S062 |
| c | Do consensus rankings associate with ED_UP? | NES/FDR point display for the 18-gene ED_UP_FDR signature. | S063 |
| d | Do top network-ranked genes overlap ED signatures? | Top-5-percent overlap sizes and within-target FDR for ED_UP/ED_DOWN. | S064 |
| e | Are target effects unusual relative to matched genes? | Matched-control empirical calibration for frozen distance metrics. | S066 |

## Statistics and source-data rule

Use the reported GSEA and hypergeometric results without re-running. ED_DOWN ranked GSEA is labelled `NOT_TESTED_SIZE_3` using the frozen eligibility rule; it is not plotted as a null enrichment result.

## Integrity and reviewer risk

Never label the output as gene upregulation/downregulation, ED reversal, causal mechanism, target validation, or wet-laboratory knockout. Donors are distinct cellular backgrounds; seeds are not biological replicates. Matched controls calibrate network-score specificity only.

