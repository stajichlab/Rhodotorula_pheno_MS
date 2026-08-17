# control_phenotype_90_110h

Canonical per-strain colony color + colony size phenotype table for this
project, chosen 2026-08-15 over the legacy YPD2 table and two other
timepoint windows — see
`analysis/examine_phenotype_calling/RECOMMENDATION.md` in this repo for the
comparison and reasoning (best agreement with YPD2 on every trait: L\*, a\*,
b\*, C\*, colony area; most complete species coverage: 303/304 strains with
a species call vs. 301/318 in YPD2).

## Where this came from

Copied verbatim (immutable, unmodified) from a sibling project on
2026-08-15:

```
/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/analysis/control_late_timepoint_phenotype/results/phenotype_control_timepoint_90_110.csv
```

SHA256: `6379c1f199d2a3f4da219210dd505c4ecd3ce65af909b67e37921d961fe69f6e`
(matches source file at copy time — reverify if the upstream file changes).

That sibling project's own documentation
(`control_late_timepoint_phenotype/CONTROL_LATE_PHENOTYPE.md`, read
2026-08-15) describes the pipeline as follows — reproduced here since this
repo depends on the file but does not own the pipeline that produced it:

> **Purpose.** A per-strain phenotype table for the control media only
> (Copper concentration = 0 mM, i.e. growth on plain YPD), taken at a late
> timepoint (hours_since_plate_start in [90, 110]) that has data for
> essentially all strains. The phenotype emphasis is colony color (CIELAB
> L\*, a\*, b\*) plus colony size, aggregated across every replicate colony
> of a strain.
>
> **Method.**
> 1. Control filter — only plates with Copper concentration = 0 mM.
> 2. Late timepoint — the imaging rig runs on two interleaved ~3h cadences,
>    so a fixed single hour would drop ~85 strains imaged only on the
>    alternate cadence. Instead, hours are rounded to the nearest integer
>    (the imaging "pass"), and each strain gets its latest pass within the
>    window [90, 110]h (passes 105h / 108h); all colonies from that pass
>    are used as the strain's replicates.
> 3. Per-colony trait values (one row of the source view = one colony at
>    one timepoint): colony size = `Shape_Area` (px); L\*/a\*/b\* =
>    `ColorLab_{L*,a*,b*}Median` (CIELAB, pixel-median of the colony's Lab
>    histogram, robust to intra-colony outliers).
> 4. Per-strain aggregation across all replicate colonies at the strain's
>    chosen pass: median, mean, sample variance, and sample SD for each of
>    area/L\*/a\*/b\*, plus `n_colonies` and `n_replicate_wells` (distinct
>    run:well spots).
>
> **Coverage (2026-08-15 re-run).** This window covers 314/320 strains
> (303 distinct `strain_code` after collapsing plate/replicate rows in
> this repo's ingestion — see `data/metadata/control_phenotype_90_110h/summary_stats.md`);
> 286 strains with ≥3 colonies at their latest pass, 28 flagged
> (`n_colonies < 3`, unreliable variance/SD); 1 strain has no control
> (Cu=0) image in any window and 5 have no late control image at all.
> Range sanity: L\* 70-80, a\* ~0.6-14, b\* -1.4-9.1, colony area
> ~0.8k-70k px.
>
> **Caveats (from upstream doc).** Variance/SD computed on small samples
> (median 4 colonies) — treat as within-strain replicate spread, not
> population variance. Strains on the alternate cadence measured ~3h
> earlier than the main group (105h vs 108h) — recorded per strain in
> `timepoint_h`. `strain_id = NULL` rows (unidentified spots) excluded.
> Sample variance (`var_samp`/`stddev_samp`) used throughout.

## Open questions for the PI (not yet answered here)

The above is what I (the agent) could reconstruct from reading the
upstream project's own documentation and duckdb schema secondhand — it
covers the *computational* pipeline but not the underlying *wet-lab*
experimental design. Please fill in (edit this section directly, or tell
me the answers and I'll write them up) in
`data/metadata/control_phenotype_90_110h/provenance.md`:

- Imaging rig / instrument make & model, and where it's physically located.
- Plate/well layout: how many colonies per strain per plate, randomization
  scheme, whether replicate colonies are biological (independent
  inoculations) or technical (same inoculation, multiple wells).
- Media recipe details beyond "YPD, Cu=0 mM" (e.g. exact YPD formulation,
  agar %, incubation temperature — YPD2's `Incubation Temp (°C)` column
  suggested 30°C for the legacy table; confirm whether this control
  pipeline used the same).
- Whether/how strains were randomized or blocked across plates (relevant
  for the "Library Plate" batch-effect consideration already flagged in
  `analysis/INTEGRATED_ANALYSIS_STRATEGY.md`'s Key Considerations).
- Contact person(s) for the `Rhodotorula_phenotypes` pipeline if different
  from the PI, for future questions.

## Known ingestion-time issues (found while building this repo's phenotype
table from this file — see `analysis/scripts/build_strain_phenotype_table.py`)

- 1 of 304 rows has no `species` value (`TFCN_17-331Y-3`).
- 1 duplicate `strain_code` in the raw file (`TFCN_48D-10`, 2 rows) —
  collapsed to the mean across rows in this repo's derived tables, per
  `build_strain_phenotype_table.py`'s documented convention.
- Colony area and a\* both show a small but persistent systematic offset
  vs. the legacy YPD2 table even at this latest timepoint window — see
  `analysis/examine_phenotype_calling/RECOMMENDATION.md`'s caveats section.
  Likely a pipeline/methodology difference between the two imaging
  systems, not a data error; don't mix absolute area/a\* values between
  YPD2- and this-source-derived tables in the same model.

## How this repo uses it

`analysis/scripts/build_strain_phenotype_table.py --source control_90_110`
now points at `data/raw/control_phenotype_90_110h/phenotype_control_timepoint_90_110.csv`
(previously read directly from the sibling project's mutable path — see
`.living/decisions.md` 2026-08-15 entry for why this was ingested).
