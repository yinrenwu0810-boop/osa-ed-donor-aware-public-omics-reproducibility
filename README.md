# OSA-ED donor-aware public-omics reproducibility archive

This repository archives the versioned analysis and figure-generation scripts, panel-level source-data tables, and frozen input manifest used for the OSA-ED public-data manuscript package dated 2026-09-13.

## Included

- `scripts/figure_build/`: Figure 1-8 build, QA, and release-freeze scripts.
- `scripts/analysis_revision_v2/`: versioned analysis scripts preserved from `revision_v2`, including virtual-perturbation workflow scripts.
- `figures/Figure*/source_data/`: panel-level source data used to render the eight main figures.
- `freeze_manifest/figure_contracts/`: Figure 1-8 contracts, the 54-item SHA-256 frozen-input manifest, and the freeze-status record.
- `RELEASE_MANIFEST.sha256.tsv`: SHA-256 inventory for every file in this release.

## Scope and scientific boundaries

This is a reproducibility archive for a public-data, donor-aware, cross-dataset hypothesis-prioritization analysis. It does not establish OSA-to-ED causality, diagnostic performance, therapeutic targets, experimental validation, or a biological mechanism. Biological independent units are donors or datasets, not cells, random seeds, paths, or matched controls.

The source data tables are figure-level derived data. The repository deliberately excludes original sequencing files, individual-level clinical records, raw single-cell objects, author-controlled submission metadata, manuscript drafts, and rendered manuscript figures. Obtain primary public data from the originating repositories cited by their GEO accessions and comply with their terms.

Some archived build scripts verify paths and SHA-256 values against the original sealed project layout; this archive preserves them verbatim for provenance. The panel-level source tables enable inspection of plotted values but do not substitute for omitted upstream raw inputs.

## Integrity check

On PowerShell, run:

```powershell
Get-Content RELEASE_MANIFEST.sha256.tsv | Select-Object -Skip 1 | ForEach-Object {
  $parts = $_ -split "`t", 2
  $actual = (Get-FileHash -LiteralPath $parts[1] -Algorithm SHA256).Hash.ToLower()
  if ($actual -ne $parts[0]) { throw "Hash mismatch: $($parts[1])" }
}
```

## Licence

This repository uses a dual licence. Source code in `scripts/` is available under the MIT License. Derived figure-level source data and documentation/frozen manifests are available under CC BY 4.0. Read [LICENSE.md](LICENSE.md) for the full scope, attribution conditions, and exclusions for upstream public data.
