# Data Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### control_phenotype_90_110h
```yaml
name: control_phenotype_90_110h
type: other                 # per-strain colony color (CIELAB) + colony size phenotype table
source: Sibling project Rhodotorula_phenotypes (Jason Stajich lab), derived from db/rhodotorula_phenotypes.duckdb
date_acquired: 2026-08-15
format: CSV (1 file)
rows: 314 (303 distinct strain_code)
columns: 26
size: 112,639 bytes (~110 KB)
raw_path: data/raw/control_phenotype_90_110h/
processed_path: (none yet -- consumed directly by analysis/scripts/build_strain_phenotype_table.py --source control_90_110)
metadata_path: data/metadata/control_phenotype_90_110h/
status: validated
known_issues:
  - 10/314 rows have no strain_code ("unidentified spots") -- currently silently dropped by this repo's groupby-based ingestion
  - 2/304 strain_code-bearing rows have no species call
  - 1 duplicate strain_code (TFCN_48D-10, 2 rows) -- collapsed to mean in derived tables
  - 28/314 strains have n_colonies < 3 -- unreliable variance/SD
  - Colony area and a* show a small systematic offset vs. the legacy YPD2 table even at this latest timepoint (pipeline/methodology difference, not error)
  - Wet-lab experimental design details (imaging rig, plate randomization, biological vs. technical replicates) not yet documented -- pending PI input, see provenance.md
access_restrictions: institutional-only
tags: [phenotype, color, CIELAB, colony-size, rhodotorula, canonical-source]
```

The canonical strain-level color/colony-size phenotype table for this
project (chosen 2026-08-15 after comparing three timepoint windows against
each other and against the legacy YPD2 table — see
`analysis/examine_phenotype_calling/RECOMMENDATION.md`). Supersedes direct
cross-project reads of the sibling `Rhodotorula_phenotypes` repo for new
analysis work; `analysis/scripts/build_strain_phenotype_table.py --source control_90_110`
now reads from `data/raw/control_phenotype_90_110h/` in this repo instead.
Wet-lab provenance (imaging/plate design) is still pending PI review — see
open questions in `data/raw/control_phenotype_90_110h/CONTROL_PHENOTYPE_90_110H.md`.

### salinity_texture_baseline
```yaml
name: salinity_texture_baseline
type: imaging                # colony-image GLCM/Haralick texture features, morphology proxy
source: Sibling project Rhodotorula_Phenotyping/Salinity (Jason Stajich lab), derived from master_measurements_combined.csv
date_acquired: 2026-09-11
format: CSV.gz (1 file)
rows: 997 (312 distinct Strain ID, 300 distinct strain_code)
columns: 87
size: 497 KB (gzip)
raw_path: data/raw/salinity_texture_baseline/
processed_path: (none yet)
metadata_path: data/metadata/salinity_texture_baseline/
status: raw
known_issues:
  - This is a quantitative TEXTURE PROXY for colony morphology, not a PI-validated smooth/rough categorical phenotype -- no such phenotype exists anywhere reachable from this project as of 2026-09-11
  - Collected under a different imaging pipeline/rig (Rhodotorula_Phenotyping/Salinity) than this project's own color phenotyping -- cross-pipeline batch effects not checked
  - 0% salinity is the unstressed baseline of a salinity-STRESS screen, not a dedicated morphology assay
  - Hours=90 is a single timepoint (not a 90-110h window average like control_phenotype_90_110h), chosen only for nominal timepoint consistency
  - 6/997 rows missing all Texture_* values; 32-112/997 rows missing Strain/Species/Origin/Environment
  - Strain ID -> strain_code crosswalk (via EXFAB_UCR-005) checked only at ID-set-overlap level (274/303 canonical panel strains matched), not spot-checked
  - Multiple replicate colony rows per strain -- no strain-level aggregation applied yet
access_restrictions: institutional-only
tags: [phenotype, morphology, texture, GLCM, Haralick, colony-imaging, rhodotorula, proxy]
```

Brought in specifically to test the PI's 2026-09-11 hypothesis that
AHL-autoinducer-like MS features correlate with colony morphology
(smooth vs. rough) — see `analysis/ahl_autoinducer_search/`. GLCM
(Haralick) texture metrics (contrast, correlation, entropy, angular
second moment, etc.) computed per colony image are the best available
quantitative stand-in for colony surface morphology in this project's
data ecosystem; treat any downstream use as testing a proxy, not a
validated trait, until/unless the PI confirms a real smooth/rough score
exists elsewhere. See `SALINITY_TEXTURE_BASELINE.md` for the full
extraction rationale and caveats.

### EXFAB_UCR-005 (legacy — pre-mycelium, not fully retrofitted)
```yaml
name: EXFAB_UCR-005
type: other                 # legacy phenotype (YPD2) + copper-AUC + MS2 sample crosswalk metadata
source: UCR EXFAB facility / lab phenotyping pipeline (pre-dates this repo's mycelium adoption)
date_acquired: unknown (files dated 2026-07-02 / 2026-08-11 in filenames)
format: CSV.gz / TSV.gz (multiple files)
rows: 318 (YPD2_phenotypic, legacy)
columns: 45
size: unknown
raw_path: (not placed under data/raw/ -- files live directly in data/metadata/EXFAB_UCR-005/, pre-dating the raw/metadata split convention)
processed_path: (none)
metadata_path: data/metadata/EXFAB_UCR-005/
status: raw
known_issues:
  - 15-16 strains missing Species call (see analysis/examine_phenotype_calling/ for comparison against control_phenotype_90_110h, which resolves most of this)
  - Does not follow the data/raw/ + data/metadata/ split convention -- files sit directly under data/metadata/EXFAB_UCR-005/ with no data/raw/EXFAB_UCR-005/ counterpart
access_restrictions: institutional-only
tags: [phenotype, color, CIELAB, copper-resistance, legacy, pre-mycelium]
```

Pre-dates this repo's mycelium adoption (2026-08-15) and does not follow
the `data/raw/` + `data/metadata/` split convention (files were placed
directly under `data/metadata/EXFAB_UCR-005/` with no immutable raw copy
elsewhere). Retained as-is rather than retroactively restructured, since
it's actively read by multiple existing scripts by that exact path
(`scripts/pcoa_color_phenotype.py`, `analysis/copper/scripts/common.py`,
etc.) — moving it would require updating every reader. Superseded by
`control_phenotype_90_110h` as the canonical phenotype source for new work
as of 2026-08-15, but kept for cross-validation
(`analysis/examine_phenotype_calling/`) and by scripts not yet migrated.
