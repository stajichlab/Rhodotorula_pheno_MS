# Summary Statistics: control_phenotype_90_110h

<!-- Generated: 2026-08-15 -->
<!-- Script: manual (python3/pandas one-off, see data/raw/control_phenotype_90_110h/CONTROL_PHENOTYPE_90_110H.md for pipeline this table came from) -->

## Overview

| Property | Value |
|----------|-------|
| Rows | 314 (303 distinct `strain_code` after excluding 10 NULL-strain_code rows and collapsing 1 duplicate code) |
| Columns | 26 |
| File size | 112,639 bytes (~110 KB) |
| Date range | N/A (not temporal; `timepoint_h` is hours-since-plate-start, not a calendar date) |
| Format | CSV |

## Column summaries

| Column | Type | Non-null | Unique | Min | Max | Mean | Top values |
|--------|------|----------|--------|-----|-----|------|------------|
| strain_id | numeric | 314 | 314 | 1 | ~320 | — | — |
| strain_code | categorical | 304 | 303 | — | — | — | `TFCN_48D-10` (2 rows) |
| genus | categorical | 302 | 3 | — | — | — | Rhodotorula (300), Cystobasidium (1), Pseudomicrostroma (1) |
| species | categorical | 302 | 18 | — | — | — | R. mucilaginosa (216), R. paludigena (17), R. toruloides (10), R. diobovata (10) |
| n_colonies | numeric | 314 | — | 1 | 8 | 3.69 (sd 0.80) | median 4 |
| n_replicate_wells | numeric | 314 | — | 1 | 8 | 3.65 (sd 0.81) | median 4 |
| timepoint_h | numeric | 314 | 2 | 105.0 | 108.0 | 107.20 | 108h (main cadence, majority), 105h (alternate cadence) |
| area_median | numeric | 314 | — | 801.0 | 69,820.5 | 32,837.7 (sd 8,673.3) | px |
| l_median | numeric | 314 | — | 70.18 | 79.82 | 74.87 (sd 1.42) | CIELAB L* |
| a_median | numeric | 314 | — | 0.60 | 14.21 | 9.81 (sd 2.27) | CIELAB a* |
| b_median | numeric | 314 | — | -1.37 | 9.07 | 3.27 (sd 1.78) | CIELAB b* |

(`*_mean`, `*_var`, `*_sd` columns for area/l/a/b omitted here for brevity —
same ranges as their `*_median` counterparts; see `schema.yaml` for full
column list.)

## Missing data summary

| Column | Missing count | Missing % | Pattern / notes |
|--------|---------------|-----------|-----------------|
| strain_code | 10 | 3.2% | "Unidentified spots" per upstream doc — currently dropped by this repo's groupby ingestion, see Known Issues in `provenance.md` |
| genus | 12 | 3.8% | Same 10 strain_code-NULL rows + 2 more with a strain_code but no species/genus call |
| species | 12 | 3.8% | Same as genus |
| origin | 12 | 3.8% | Same as genus/species |
| environment | 34 | 10.8% | Broader gap than genus/species — some identified strains still lack an environment label |
| area_var / area_sd | 14 each | 4.5% | Rows with n_colonies < 2 (variance undefined for a single replicate) |
| l_var / l_sd / a_var / a_sd / b_var / b_sd | 14 each | 4.5% | Same n_colonies < 2 rows as area_var/sd |

## Quality flags

- 28/314 rows have `n_colonies < 3` — upstream doc flags these as having
  unreliable variance/SD estimates; downstream analyses relying on
  within-strain spread (not just the median/mean) should filter or flag
  these explicitly.
- 10/314 rows have no `strain_code` at all and are effectively invisible
  to every downstream table in this repo (silently dropped by
  `pandas.groupby`) — real colony measurements exist for these spots but
  aren't currently usable without a corrected/matched strain_code. Worth
  a manual look at the source imaging records if these strains matter for
  final analyses (a few appear to correspond to strains that also lack a
  species call in the legacy YPD2 table, e.g. "17-334P-4", "25-325P-1" —
  possibly the same underlying strain-ID gaps recurring across both
  phenotyping pipelines).
- Colony area and a\* show a systematic (not just noisy) offset vs. the
  legacy YPD2 table even at this latest timepoint — see cross-reference in
  `provenance.md`'s Known Issues and
  `analysis/examine_phenotype_calling/RECOMMENDATION.md`.

## Notes

- This table is the **canonical phenotype source** for
  `analysis/scripts/build_strain_phenotype_table.py --source control_90_110`
  as of 2026-08-15 (see `.living/decisions.md`).
- Values were cross-validated against two other timepoint windows
  (`control_70_80`, `control_80_90`) and the legacy YPD2 table in
  `analysis/examine_phenotype_calling/` before being chosen as canonical —
  see that directory's `RECOMMENDATION.md` for the full comparison
  (correlation plots, agreement statistics).
