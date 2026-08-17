# Provenance: control_phenotype_90_110h

## Source

**Type**: derived

**Origin**:
- Derived from: `db/rhodotorula_phenotypes.duckdb` (`v_phenotype` view joined
  with `condition_plate_factor` and `strain` tables), in the sibling
  project `Rhodotorula_phenotypes`
  (`/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/`). Built
  by that project's `analysis/control_late_timepoint_phenotype/scripts/build_phenotype_table.py`,
  reproducible there via `bash analysis/control_late_timepoint_phenotype/run.sh`.
  This repo ingests the resulting CSV output, not the upstream duckdb
  itself.

**Citation / accession**: N/A (internal lab database, not a public accession)

<!-- PI: is there an upstream sample/imaging accession or LIMS ID this should reference? -->

## Acquisition details

**Date acquired**: 2026-08-15

**Obtained by**: Claude (agent), on request from Jason Stajich (PI), copied
verbatim from the sibling project's `results/` output directory.

**Method**: `cp` from
`/bigdata/stajichlab/jstajich/projects/Rhodotorula_phenotypes/analysis/control_late_timepoint_phenotype/results/phenotype_control_timepoint_90_110.csv`
to `data/raw/control_phenotype_90_110h/` in this repo, no transformation.
File mtime at time of copy: 2026-08-15 16:32.

**Checksum**: SHA256 `f1070e6e0b01b7e1f2eaf869d078dae981a502ecd2f168e3d1d73a173e945c47`

> Note: The `species` column for 4 strains was corrected on 2026-08-16 to
> match the whole-genome fastANI reassignment recorded in
> `analysis/integrated_analysis/phase_siderophore/ANI_check/SPECIES_REASSIGNMENT.md`
> (strains `TFCN_1A-1-3`, `TFCN_1B-1-2` R. pacifica→R. mucilaginosa,
> `TFCN_1A-1-2` R. paludigena→R. mucilaginosa, `TFCN_152C-6`
> R. toruloides→R. taiwanensis). All other values unchanged; no rows
> added/removed.

## Access restrictions

**Restriction level**: institutional-only

**Details**: Lab-internal data (Stajich lab, UCR). No formal DUA on file
that the agent is aware of. <!-- PI: confirm/correct restriction level. -->

## Known issues

- **10/314 raw rows have no `strain_code`** ("unidentified spots" per
  upstream doc) and are currently silently dropped by this repo's
  `groupby(strain_code)` ingestion step (pandas drops NaN group keys) —
  their color/area measurements exist in the raw file but do not
  currently reach any strain-level table in this repo. See
  `data/raw/control_phenotype_90_110h/CONTROL_PHENOTYPE_90_110H.md` for
  the affected `strain_id`/`strain_name` values.
- 2/304 strain-code-bearing rows have no `species` call (`TFCN_17-331Y-3`,
  and one of the two `TFCN_48D-10` rows — the duplicate strain_code, see
  next point).
- 1 duplicate `strain_code` (`TFCN_48D-10`, 2 rows, different
  plates/replicates) — collapsed to the mean in this repo's derived
  tables (`build_strain_phenotype_table.py`), not in the raw file itself.
- 28/314 strains have `n_colonies < 3` at their chosen imaging pass —
  upstream doc flags these as having unreliable variance/SD estimates
  (small-sample within-strain replicate spread, not population variance).
- Colony area and a\* both show a small, persistent systematic offset vs.
  the legacy YPD2 phenotype table (`data/metadata/EXFAB_UCR-005/`), even
  at this latest timepoint window — see
  `analysis/examine_phenotype_calling/RECOMMENDATION.md`. Likely a
  pipeline/methodology difference (different imaging systems), not a data
  error in either table — don't mix absolute area/a\* values between the
  two sources in the same statistical model.
- **Wet-lab experimental design details not yet documented here** —
  imaging rig identity, plate layout/randomization, biological vs.
  technical replicate structure, exact media recipe. See "Open questions
  for the PI" in `CONTROL_PHENOTYPE_90_110H.md`. <!-- PI: please fill in
  or dictate answers; this provenance record is incomplete without them. -->

## Contact

**Primary contact**: Jason Stajich (PI), jasonst@ucr.edu

**Backup contact**: <!-- PI: who else knows the Rhodotorula_phenotypes imaging/plate pipeline, if not you? -->

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-08-15 | Initial ingestion of the 90-110h timepoint window as this repo's canonical phenotype source, superseding direct cross-project reads of the (now-removed) single `control_late` file and the legacy YPD2 table for new analysis work. |
| 1.1 | 2026-08-16 | Corrected `species` for 4 strains per whole-genome fastANI reassignment (see note above the checksum). New SHA256 saved. |
