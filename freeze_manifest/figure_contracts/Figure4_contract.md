# Figure 4 Contract — Candidate hierarchy and doublet sensitivity

**Core conclusion:** TYMS, EFNB2, and LRRC17 are the only ED-FDR-supported candidates under the frozen hierarchy, and their evidence classes persist after doublet exclusion within the same donor atlas.

- Figure archetype: quantitative grid.
- Target/output: Scientific Reports main Figure 4; Python; 183 × 140 mm; SVG/PDF/TIFF/PNG.
- Hero evidence: candidate evidence matrix plus before/after doublet results.
- Status: `READY_FROM_SEALED_TABLES`; no new test.

## Panel map

| Panel | Unique question | Frozen display | Source IDs |
|---|---|---|---|
| a | How many candidates have FDR versus nominal-P ED support? | Two-level count display: 3 ED-FDR-supported and 563 nominal-P exploratory. | S040–S041 |
| b | Which evidence dimensions support each candidate? | Matrix for endothelial direction/meta FDR, ED cell-type effect/FDR, and doublet-retained class. Historical composite ranking and centrality are excluded. | S034, S041, S044 |
| c | Are three-gene effects altered after doublet exclusion? | Paired before/after log2FC with FDR encoded separately, faceted by cell type/contrast. | S042 |
| d | Are pathway estimates globally stable after exclusion? | Old versus new Hallmark NES scatter with identity line; color indicates significance-state transition. | S043 |
| e | How many pathway states are retained, gained, or lost? | Deterministic counts from old/new FDR<0.05 states. | S043 |

## Statistics and source-data rule

All P values and FDRs are copied from the sealed before/after analyses. Panel e is a contingency count only. The full 49-state-change list remains supplementary/source data.

## Integrity and reviewer risk

Doublet exclusion is a sensitivity analysis in the same dataset, not independent replication. Multiple cell-type readouts from the same donors are nested observations and must not be counted as separate confirmations.

