# Which phenotype table to use going forward

## Verdict: `control_90_110` (105/108h window)

Use `phenotype_control_timepoint_90_110.csv` from the sibling
`Rhodotorula_phenotypes` project as the canonical strain-level color/size
phenotype source going forward, via
`analysis/scripts/build_strain_phenotype_table.py --source control_90_110`
(now the script's default). This is a provisional call pending formal
`mycelium:ingest` of that file into `data/metadata/` -- see the open item
below.

## Why

1. **Best agreement with the original YPD2 table, on every trait, of the
   three timepoint windows.** From `timepoint_comparison_summary.txt`
   (YPD2 vs. each window, n=298 matched strains):

   | trait | vs 70-80h | vs 80-90h | vs 90-110h |
   |---|---|---|---|
   | L\* | r=0.911 | r=0.928 | **r=0.944** |
   | a\* | r=0.955 | r=0.963 | **r=0.969** |
   | b\* | r=0.973 | r=0.984 | **r=0.990** |
   | C\* | r=0.947 | r=0.958 | **r=0.967** |
   | area | r=0.909 | r=0.927 | **r=0.952** |

   90-110h wins on every single trait, not just on average -- a consistent
   signal, not noise.

2. **Real, monotonic timepoint drift, not just measurement noise.**
   `timepoint_progression_scatter.png` compares the three control windows
   against each other: adjacent windows (70-80 vs 80-90, 80-90 vs 90-110)
   agree more tightly than distant ones (70-80 vs 90-110) for every trait.
   If the differences between windows were just per-image noise, all three
   pairwise comparisons would look equally scattered regardless of how far
   apart the windows are -- they don't. Colonies are still visibly
   developing between 70h and 110h, and 90-110h is the closest of the three
   to both plateau and to the original (also late-timepoint) YPD2
   measurement.

3. **Best species/strain-code coverage.** 303 of 304 strains have a species
   call (vs. 301/316 in YPD2 -- 15 missing), and several strain codes that
   were typo'd/unlabeled in YPD2 were corrected (see
   `analysis/INTEGRATED_ANALYSIS_STRATEGY.md`'s Phase 1 status section for
   the full strain-membership diff).

## Caveats to carry forward, not to ignore

- **Colony area has a persistent systematic offset vs. YPD2 at every
  window, including 90-110h** (`area (control)` sits below `area (YPD2)`
  for most strains in `ypd2_vs_control_windows_scatter.png` -- not just
  scatter, a real shift). This did not shrink to zero even at the latest
  window, which points to a **methodology difference between the two
  imaging/measurement pipelines**, not purely a growth-stage effect. Treat
  absolute colony-area values as pipeline-specific; area *comparisons
  within* one source (e.g. across control_90_110 strains) are fine, but
  don't mix area values from YPD2 and control_* sources in the same model.
- **a\* has a smaller version of the same systematic offset** (control a\*
  runs slightly below YPD2 a\* even at 90-110h, most visible in the
  low-to-mid a\* range) -- worth remembering since a\* is a direct input to
  the orange_score composite. It doesn't change the species ranking
  finding (reproduced independently on both YPD2 and control-derived
  sources in Phase 1), but absolute orange_score values are not directly
  comparable across sources.
- **`DBVPG_3771`** (a genome-sequenced *R. dairenensis* strain) is absent
  from all three control_* windows (present only in YPD2) -- still a real
  strain-coverage loss for the focal species, independent of which window
  is chosen.
- **Not yet formally ingested.** All three control_* files are still read
  directly from `/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/...`
  in `build_strain_phenotype_table.py`'s `SOURCES` dict. If `control_90_110`
  is confirmed as the standing choice, run it through `mycelium:ingest` into
  `data/metadata/` with provenance (what "control"/timepoint window mean,
  how it supersedes YPD2) so the dependency isn't on a sibling project's
  mutable file layout.

## Supporting files in this directory

- `compare_phenotype_sources.py` -- first pass, YPD2 vs. the (now removed)
  single `control_late` file.
- `compare_timepoint_sources.py` -- the full comparison behind this
  recommendation: all 3 control windows vs. each other, and vs. YPD2, for
  color (L\*/a\*/b\*/C\*) and colony size (area).
- `timepoint_progression_scatter.png/.pdf`, `ypd2_vs_control_windows_scatter.png/.pdf`,
  `timepoint_comparison_summary.txt` -- the plots and correlation tables
  this recommendation is based on.
