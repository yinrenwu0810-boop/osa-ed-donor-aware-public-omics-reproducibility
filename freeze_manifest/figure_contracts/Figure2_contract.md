# Figure 2 Contract — Endothelial convergence within cross-study heterogeneity

**Core conclusion:** Endothelial IH responses converge in a strictly prioritized subset and at pathway level, while genome-wide effect agreement remains limited and individual datasets retain null FDR findings.

- Figure archetype: quantitative grid with a dominant cross-study effect panel.
- Target/output: Scientific Reports main Figure 2; Python; 183 × 155 mm; SVG/PDF/TIFF/PNG.
- Hero evidence: shared-gene effect density plus the 258 strict LODO-FDR subset.
- Status: `READY_AFTER_RULED_DERIVATION`; no new inferential analysis.

## Panel map

| Panel | Unique question | Frozen display | Source IDs |
|---|---|---|---|
| a | Do all three studies yield the same single-gene signal? | Three aligned volcano plots with identical axes where feasible and explicit FDR-significant counts, including zero-result studies. | S020–S022 |
| b | How similar are shared effects? | Pairwise hexbin/density plots for shared genes, with the already reported Spearman coefficient annotated. | S023–S024 |
| c | What survives progressively stricter rules? | Count/upset-style summary for same direction, stability-screen, full-meta FDR, and all-three LODO-FDR sets. | S010–S011 |
| d | Are strict genes driven by one study? | Effect heatmap for a deterministic subset of strict LODO-FDR genes across all three studies; no forest plot because sealed standard errors are not registered. | S011 |
| e | Does pathway-level evidence add interpretable structure? | Dot plot for frozen Hallmark and Reactome subsets, encoding NES/effect direction and FDR without merging correction families. | S025–S027 |

## Mechanical derivation rules

Panel d takes the 12 lowest full-meta FDR genes among `strict_lodo_FDR_candidate=True`, with six positive and six negative when available; ties resolve alphabetically. Panel e includes HALLMARK_HYPOXIA plus the five Hallmark pathways with lowest `endothelial_max_FDR`, and up to six Reactome pathways with FDR<0.05 selected by absolute NES then alphabetical tie-break. Hallmark and Reactome are faceted and never pooled for FDR.

## Integrity and reviewer risk

The figure must not imply universal gene-level agreement. Dataset is the replication level for cross-study stability; genes are tested features, not replicates. Zero FDR hits are reported as results rather than hidden.

