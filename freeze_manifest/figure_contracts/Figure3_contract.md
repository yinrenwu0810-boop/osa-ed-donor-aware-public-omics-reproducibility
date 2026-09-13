# Figure 3 Contract — ED corpus-cavernosum localization

**Core conclusion:** IH-linked pathways and the prioritized genes localize across multiple ED corpus-cavernosum compartments, but inference remains donor-level and group is inseparable from procurement context in this atlas.

- Figure archetype: asymmetric mixed-modality figure.
- Target/output: Scientific Reports main Figure 3; Python; 183 × 155 mm; SVG/PDF/TIFF/PNG.
- Hero evidence: donor-pseudobulk pathway and three-gene compartment effects.
- Status: `READY_AFTER_RULED_DERIVATION`; no new test.

## Panel map

| Panel | Unique question | Frozen display | Source IDs |
|---|---|---|---|
| a | What cell states are represented? | Descriptive UMAP of 64,993 cells, clearly labelled non-inferential. | S030 |
| b | Are cell proportions dominated by particular donors? | Donor × cell-type stacked composition; donor labels grouped by condition without cell-level P values. | S031 |
| c | Which compartments carry IH-linked programmes in ED? | Donor-pseudobulk Hallmark NES/FDR heatmap for prespecified structural and immune compartments. | S033 |
| d | Where do TYMS, EFNB2, and LRRC17 show ED effects? | Three-gene × cell-type log2FC/FDR matrix with missing estimates shown as blank. | S034 |
| e | What do the donor-level observations look like? | Optional descriptive donor dots for the three genes only, generated from pseudobulk counts under the rule below. | S035–S036 |

## Descriptive derivation rule for panel e

Panel e is allowed only if the existing pseudobulk sample table identifies all plotted donor columns unambiguously. Display log2 counts-per-million with a fixed 0.5 prior count, facet by cell type, and do not add a new hypothesis test or significance symbol. If mapping is ambiguous, omit panel e and record `HOLD_PANEL_E_SOURCE_MAPPING`; the figure can still proceed with panels a–d.

## Integrity and reviewer risk

The legend must state 3 tumour-margin reference donors, 3 non-diabetic ED implantation donors, and 2 diabetic ED donors where shown. Cells are not independent replicates. UMAP proximity and cell abundance are descriptive. The normal-reference label must preserve the author-designated tumour-margin context.

