# Figure 7 Contract — Exploratory Hypoxia robustness and nonspecificity

**Core conclusion:** Across frozen and exploratory specifications, all three targets retain a positive Hypoxia direction, but none achieves Hallmark FDR<0.05 and none is specific relative to matched controls.

- Figure archetype: quantitative grid.
- Target/output: Scientific Reports main Figure 7; Python; 183 × 155 mm; SVG/PDF/TIFF/PNG.
- Hero evidence: S0–S3 NES beside matched-control empirical calibration.
- Status: `READY_WITH_PANEL_B_REVISION`; explicitly exploratory; no new analysis.

## Panel map

| Panel | Unique question | Frozen display | Source IDs |
|---|---|---|---|
| a | Is the target-level direction stable to specification? | EFNB2, LRRC17, and TYMS S0–S3 NES; S3 marked descriptive with no P/FDR. | S070–S071 |
| b | Does any inferential specification pass FDR? | Hypoxia FDR for S0, S1, and S2 only, with 0.05 line; S3 shown as not inferential. | S071 |
| c | Is direction stable to seeds and donors? | Distributions of 45 seed-level, 45 leave-one-seed, and 9 leave-one-donor NES values. | S072–S074 |
| d | Are target NES values unusual versus matched genes? | Signed percentile and conservative absolute upper-tail empirical value for each target. | S075 |
| e | Is the leading edge stable? | Jaccard summaries across specifications and within target/donor/seed comparisons. | S071, S076 |
| f | Was the gene set adequately mapped? | 191/200 exact mapping with the nine unavailable symbols listed compactly. | S077 |

## Frozen revision from merged V3

The V3 plan proposed a full 32-Hallmark S0–S2 FDR panel. The registered append-only outputs contain target-level S0–S3 Hypoxia results but do not contain a complete 32-set S1/S2 table. Panel b is therefore narrowed to the available Hypoxia FDR values. Recreating all 32 sets would be a new analysis and is not authorized in this stage.

## Integrity and reviewer risk

Figure title and legend must contain “exploratory sensitivity analysis”. Positive NES is a rank-direction result, not expression upregulation. The fixed interpretation is `robust positive direction but nonspecific`; stability cannot replace the formal FDR-negative result.

