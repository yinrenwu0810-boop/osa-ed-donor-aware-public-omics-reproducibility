# Figure 1 Contract — Parallel evidence layers and strict clinical bridge

**Core conclusion:** IH exposure, OSA/CPAP clinical context, and ED target tissue are parallel evidence layers, and HALLMARK_HYPOXIA is the only pathway satisfying the complete prespecified strict bridge.

- Figure archetype: schematic-led composite.
- Target/output: Scientific Reports main Figure 1; Python; 183 × 120 mm; SVG/PDF/TIFF/PNG.
- Hero evidence: the strict bridge screen and the four-context Hypoxia pattern.
- Status: `READY_AFTER_MECHANICAL_DERIVATION`; no new statistical test.

## Panel map

| Panel | Unique question | Frozen display | Source IDs |
|---|---|---|---|
| a | What evidence layers were kept distinct? | Left-to-right parallel study map with datasets, groups, and independent units; no causal arrows. | S001–S003 |
| b | How strongly did endothelial evidence narrow? | Nested/count flow: 5,139 common → 1,282 full-meta FDR; alongside 1,766 same-direction → 1,035 stability screen → 258 strict LODO-FDR. | S010–S011 |
| c | How did the clinical bridge narrow? | 50 total → 49 estimable → 12 complete directional → 3 all-four FDR → 1 strict bridge. | S012–S013 |
| d | Which pathway survives and why? | NES heatmap for the 12 complete-directional pathways across the four prespecified strict contexts; FDR<0.05 marked separately. | S013–S014 |

## Statistics and source-data rule

No recomputation. Counts and NES/FDR values are copied exactly. Panel d includes all 12 pathways with `complete_directional_pattern=True`, ordered with the sole strict bridge first and the remainder by pathway name; this prevents favorable post-hoc ranking.

## Integrity and reviewer risk

The figure must say that OSA/CPAP support is a clinical-context bridge, not proof that Hypoxia mediates OSA-associated ED. “Normal” is replaced by “author-designated normal-erection tumour-margin reference” wherever relevant. Empty cells mean not estimated, never zero.

