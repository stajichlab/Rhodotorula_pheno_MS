# Provenance: salinity_texture_baseline

## Source

**Type**: derived (filtered subset of a sibling project's dataset)

**Origin**:
- Derived from: `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Phenotyping/Salinity/analysis/master_measurements_combined.csv`
  (Jason Stajich lab, `Rhodotorula_Phenotyping/Salinity` project — a
  salinity-stress colony-imaging screen; 216 MB, 129,627 rows, all
  salinity levels 0/1/3/5% w/v and timepoints 0-100h, YPDN media).
  Transformation applied: filtered to rows where
  `Salinity (%w/v) == 0 AND Hours == 90 AND Media == "YPDN"`, and
  restricted to the identity/shape/texture columns listed in
  `schema.yaml` (dropped per-image intensity/bounding-box/grid columns
  not relevant to morphology).

**Citation / accession**: N/A (internal lab project, not published)

## Acquisition details

**Date acquired**: 2026-09-11

**Obtained by**: Jason Stajich (jason.stajich@ucr.edu), via Claude Code
session in `Rhodotorula_Metabolites/Rhodotorula_pheno_MS`

**Method**: `pandas.read_csv` subset + column filter, run interactively;
no dedicated ingestion script committed (single one-off subset — see
decision log entry 2026-09-11 in `.living/decisions.md`).

**Checksum**: not recorded for either the source file or this subset.

## Access restrictions

**Restriction level**: institutional-only

**Details**: Same restriction level as this project's other strain
phenotype data (`control_phenotype_90_110h`, `EXFAB_UCR-005`) — internal
lab strain collection data.

## Known issues

- This is a **proxy** for colony morphology (smooth/rough), not a
  PI-validated categorical morphology call — none exists anywhere
  reachable from this project as of 2026-09-11 (see
  `.living/learnings.md` same date).
- Collected under a **different imaging pipeline/rig** than this
  project's own color phenotyping (`control_phenotype_90_110h`, sourced
  from the sibling `Rhodotorula_phenotypes` pipeline). Cross-pipeline
  batch effects have not been checked.
- 0% salinity is the **unstressed baseline of a salinity-stress
  screen**, not a dedicated morphology assay — the strains' physiological
  state at the time of this screen may differ from the state during the
  MS/color phenotyping campaign (different stock-culture age, plate
  batch, imaging date).
- `Hours=90` is a **single timepoint**, chosen only for nominal
  consistency with `control_phenotype_90_110h`'s 90-110h window; it is
  not itself a window average.
- 6/997 rows have missing texture values across the full `Texture_*`
  column family (same rows); 32-112/997 rows have missing
  Strain/Species/Origin/Environment identity fields.
- `Strain ID` -> `strain_code` crosswalk relies on
  `data/metadata/EXFAB_UCR-005/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz`
  and was checked only at the level of ID-set overlap (274/303 canonical
  panel strains matched), not spot-checked strain-by-strain.
- Multiple replicate colony rows exist per strain — no strain-level
  aggregation has been applied in this raw table.

## Contact

**Primary contact**: Jason Stajich, UCR (jason.stajich@ucr.edu)

**Backup contact**: None.

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-09-11 | Initial ingestion: 0% salinity / Hours=90 / YPDN subset for the AHL-autoinducer/morphology hypothesis |
