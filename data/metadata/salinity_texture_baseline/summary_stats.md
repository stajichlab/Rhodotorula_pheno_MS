# Summary Statistics: salinity_texture_baseline

<!-- Generated: 2026-09-11 -->
<!-- Script: manual (pandas, interactive) -->

## Overview

| Property | Value |
|----------|-------|
| Rows | 997 (colony-level; 312 distinct `Strain ID`, 300 distinct `Strain` codes) |
| Columns | 87 |
| File size | 497 KB (gzip, from 1.4 MB CSV) |
| Date range | N/A (single timepoint: Hours=90) |
| Format | CSV.gz |

## Column summaries

| Column | Type | Non-null | Unique | Min | Max | Mean | Top values |
|--------|------|----------|--------|-----|-----|------|------------|
| Strain ID | numeric | 997 | 312 | 1 | 325 | — | — |
| Species | categorical | 964 | 21 | — | — | — | R. mucilaginosa (712), R. paludigena (45), R. dairenensis (35) |
| Salinity (%w/v) | numeric | 997 | 1 | 0 | 0 | 0 | fixed by construction |
| Hours | numeric | 997 | 1 | 90 | 90 | 90 | fixed by construction |
| Shape_Area | numeric | 997 | — | 10 | 160,297 | 39,320 | px (10-px minimum is likely a segmentation artifact, not a real colony) |
| Shape_Solidity | numeric | 997 | — | see notes | see notes | — | macro-scale outline irregularity |
| Texture_AngularSecondMoment-avg-scale05 | numeric | 991 | — | 0.039 | 0.239 | 0.123 | std 0.027 |
| Texture_Contrast-avg-scale05 | numeric | 991 | — | 0.503 | 9.715 | 0.910 | std 0.535 (right-skewed: a few very high-contrast colonies) |
| Texture_Correlation-avg-scale05 | numeric | 991 | — | -0.306 | 0.860 | 0.657 | std 0.094 |
| Texture_Entropy-avg-scale05 | numeric | 991 | — | 2.671 | 5.121 | 3.701 | std 0.268 |
| (9 remaining `-avg-scale05` Haralick metrics) | numeric | 991 | — | — | — | — | see `schema.yaml` for definitions; full describe() not reproduced here for brevity |
| (52 directional `-degXXX-scale05` columns) | numeric | 991 | — | — | — | — | same 13 metrics x 4 angles; use for anisotropy checks, not as the primary summary |

## Missing data summary

| Column | Missing count | Missing % | Pattern / notes |
|--------|---------------|-----------|-----------------|
| Environment | 112 | 11.2% | Carried over from source; likely unscored/unknown for some strains |
| Origin | 34 | 3.4% | Carried over from source |
| Species | 33 | 3.3% | Same 33(ish) rows overlap with missing Strain |
| Strain | 32 | 3.2% | Rows with no strain assignment ("unidentified spots"-type, consistent with `control_phenotype_90_110h`) |
| All 65 `Texture_*` columns | 6 | 0.6% | Same 6 rows across the whole texture-column family — a colony-segmentation edge case in the source pipeline |

## Quality flags

- `Shape_Area` minimum of 10 px is almost certainly a segmentation
  artifact (a speck, not a colony) — any strain-level aggregation should
  either filter tiny objects or use median/robust statistics across
  replicates, not the raw mean.
- `Texture_Contrast-avg-scale05` and several other metrics are
  right-skewed (max >10x the median) — a handful of colonies with
  unusually high local contrast. Worth flagging as either real rough
  outliers or imaging artifacts (debris, condensation) before treating as
  biology.
- 300 distinct `Strain` codes but 312 distinct numeric `Strain ID` values
  — 12 IDs have no resolvable strain code in this subset (consistent with
  the 32 rows missing `Strain`).
- Multiple replicate rows exist per strain (up to several) — no
  strain-level collapse has been done in this raw file.

## Notes

Only the `-avg-scale05` (rotation-invariant) columns are recommended as
the starting point for any strain-level morphology summary; the 4
directional angle columns are retained for anisotropy diagnostics (e.g.
distinguishing radial ridging from isotropic roughness) but are unlikely
to be the primary comparison axis.
