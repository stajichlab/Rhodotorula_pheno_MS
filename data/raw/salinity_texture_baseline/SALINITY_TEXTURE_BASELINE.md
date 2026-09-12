# salinity_texture_baseline

Colony-image GLCM (gray-level co-occurrence matrix / Haralick) texture
features for the 0% w/v salinity control condition, Hours=90, YPDN media,
subsetted from the sibling project's salinity-stress screen:

`/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Phenotyping/Salinity/analysis/master_measurements_combined.csv`

(216 MB, 129,627 rows across all salinity levels 0/1/3/5% w/v and
timepoints 0-100h; NOT copied here in full — see below).

## Why this dataset

Brought in 2026-09-11 to test the PI's hypothesis that AHL-autoinducer-
like MS features correlate with colony morphology (smooth vs. rough).
This project (`Rhodotorula_pheno_MS`) had no colony-texture/morphology
phenotype of any kind — only CIELAB color and shape-area
(`control_phenotype_90_110h`, `EXFAB_UCR-005`). The Salinity screen's
per-colony image analysis pipeline computes standard Haralick texture
descriptors (contrast, correlation, entropy, angular second moment,
inverse difference moment, etc., at 4 GLCM angles + an averaged
`-avg-scale05` summary), which are a legitimate quantitative proxy for
colony surface texture/roughness. This is **not** a PI-validated
categorical smooth/rough call — it is the best available quantitative
stand-in, adopted for lack of one.

## What was extracted vs. left out

Only the **0% w/v salinity, Hours=90, Media=YPDN** subset was extracted
(997 colony-level rows, 312 distinct numeric `Strain ID`, 300 distinct
`Strain` codes after dropping rows with missing strain metadata) — the
1%/3%/5% salinity conditions and other timepoints in the source file are
out of scope for this hypothesis and were not copied in, to avoid
duplicating an unrelated (and much larger) experiment inside this repo.
Columns retained: strain/species identity columns, `Salinity (%w/v)`,
`Hours`, `Media`, replicate/plate columns, `Shape_*` (colony size/shape,
kept for QC and as a decoy/confound check — colony area is already a
known confound elsewhere in this project, see
`.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md`),
and all `Texture_*` GLCM columns. Full per-pixel image data, other
salinity levels, and other timepoints remain only in the source project.

## Why Hours=90

Chosen to match this project's existing `control_phenotype_90_110h`
color-phenotype convention (90-110h post-inoculation timepoint window).
Hours=90 gave the best strain coverage among the available discrete
timepoints near that window (90/95/96/100h all present; 90h alone covers
312/318 strains vs. ~150-160 at the other three). This is a single
timepoint, not a 90-110h window average like the color table -- the
Salinity screen's imaging cadence doesn't have per-strain data across
that whole window at every strain.

Stored as `salinity0_hours90_texture.csv.gz` (gzip -9, ~497 KB, matching
this repo's existing `.csv.gz` convention for `EXFAB_UCR-005`) rather than
plain CSV; `pandas.read_csv` reads it transparently.

## Strain-code crosswalk

The source file's `Strain ID` is a numeric ID from the same ID space as
`data/metadata/EXFAB_UCR-005/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz`'s
`Strain ID` column, which maps to this project's `strain_code`
(`Strain` column in that crosswalk file, matching
`control_phenotype_90_110h`'s `strain_code`). Verified: 274/303 strains
in the canonical color-phenotype panel have a Hours=90/Salinity=0 texture
record here. The `Strain` column in this dataset is carried over directly
from the source file's own crosswalk and was NOT independently
re-verified strain-by-strain beyond the ID-overlap check above.

## Known caveats (see also `data/metadata/salinity_texture_baseline/provenance.md`)

1. **Proxy, not ground truth**: this is GLCM texture computed from
   photographs under the Salinity project's imaging rig, not a
   PI-scored or validated smooth/rough phenotype.
2. **Cross-project imaging batch effect, unchecked**: this data was
   collected under a different imaging setup/pipeline
   (`Rhodotorula_Phenotyping/Salinity`) than this project's own color
   phenotyping (`control_phenotype_90_110h`, from the sibling
   `Rhodotorula_phenotypes` pipeline). Whether the two imaging rigs
   introduce systematic differences that could confound
   texture-vs-anything comparisons has not been tested.
3. **Repurposed control condition, not a dedicated morphology assay**:
   0% salinity is the *unstressed baseline* of a salinity-stress
   screen, not an experiment designed to characterize morphology per se.
   Whatever biological state the strains were in at the time of that
   screen (which may differ in age of stock culture, plate batch, etc.
   from the MS/color phenotyping campaign) is what's captured here.
4. **6 rows have missing texture values** (colony segmentation/edge
   cases in the source pipeline, not addressed here) and 32-112 rows
   have missing Strain/Species/Origin/Environment metadata (likely the
   same "unidentified spots" issue documented for
   `control_phenotype_90_110h`).
5. Multiple replicate colony rows per strain (up to several per Strain
   ID) -- any strain-level use requires an explicit aggregation
   decision (mean/median across replicates), not yet made in this raw
   table.
