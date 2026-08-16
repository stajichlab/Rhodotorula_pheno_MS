# Within-species color↔metabolome sweep: all species with ≥5 strains (2026-08-16)

Script: `analysis/scripts/phase2_within_species_association.py --min-strains 5`
(same design as `WITHIN_SPECIES_MUCILAGINOSA.md` and
`ROBUSTNESS_AND_MULTIVARIATE.md`'s 2-species robustness check, extended to
every remaining species with at least 5 strains having both color and MS
data). a\* only (primary predictor); C\* not run for these small species
(see Recommendation below).

## Species covered

| Species | n strains (phenotype+MS) | Cell fraction n | Supernatant n | Result (a\*) |
|---|---|---|---|---|
| *R. mucilaginosa* | 206 | — | — | Null — see `WITHIN_SPECIES_MUCILAGINOSA.md` (not rerun here) |
| *R. paludigena* | 10 | — | — | Null — see `ROBUSTNESS_AND_MULTIVARIATE.md` (not rerun here) |
| *R. toruloides* | 10 | — | — | Null — see `ROBUSTNESS_AND_MULTIVARIATE.md` (not rerun here) |
| *R. dairenensis* | 8 | 7 | 8 | **Null**, 0/10,164-10,416 hits (2 features excluded as constant: 4,909 cell / 2,768 supernatant) |
| *R. diobovata* | 8 (1 strain has no genome-tree tip) | 8 | 8 | **Null**, 0 hits (3,283 / 2,735 constant-excluded) |
| *R. taiwanensis* | 6 | 6 | 6 | **Null**, 0 hits (5,619 / 3,081 constant-excluded) |
| *R. sp. clade I* | 5 | 5 | 5 | **Null**, 0 hits (5,289 / 3,197 constant-excluded) |
| *R. sphaerocarpa* | 5 (1 strain has no genome-tree tip) | 5 | 5 | **Null**, 0 hits (4,759 / 3,112 constant-excluded) |

All 8 species tested to date (2 large + 6 small) come back null for
a\*↔metabolome at BH-FDR<0.05.

## Result files
Each species produced 2 CSVs in this directory (`within_species_<Species>_association_a.csv`,
`within_species_<Species>_association_area_decoy.csv`), one row per
deduplicated feature group × fraction:

- `within_species_Rhodotorula_dairenensis_association_{a,area_decoy}.csv`
- `within_species_Rhodotorula_diobovata_association_{a,area_decoy}.csv`
- `within_species_Rhodotorula_taiwanensis_association_{a,area_decoy}.csv`
- `within_species_Rhodotorula_sp_clade_I_association_{a,area_decoy}.csv`
- `within_species_Rhodotorula_sphaerocarpa_association_{a,area_decoy}.csv`

## Important caveat: the negative control is uninformative at this n

For *R. mucilaginosa* (n=206), the area decoy showed a real, strong signal
(1,524/10,164 cell-fraction hits) — proof the pipeline can detect a true
effect at that power, which is what makes the a\*/C\* null trustworthy
there. **For all 5 new small species (n=5-8), the area decoy also came
back at 0 hits** — but this is expected purely from low power (constant-
feature exclusion alone removes 48-51% of cell-fraction features at these
sample sizes), not evidence the pipeline is well-calibrated at this n.
**The hard gate in the script only checks `n_perm >= 100`, not that the
decoy actually found something** — so it passed formally, but a 0-hit
decoy at n=5-8 does not carry the same evidential weight as it does at
n=206 or n=10 (where `ROBUSTNESS_AND_MULTIVARIATE.md`'s decoys, while
also null, at least had more features survive the constant-abundance
filter). **These 6 small-species results should be read as "no signal
detected, in a test that may not have been able to detect one," not as
independently confirmed nulls.** This mirrors this project's established
pattern (`.living/findings/phenotype-metabolome-association-statistical-power.md`)
of power-ceiling nulls recurring at every scale tried so far.

## Interpretation

Combined with the 2 already-tested larger species, this sweep now covers
every species in the panel with ≥5 strains and MS data (8 of 17-18
species; the remaining 9-10 species have 1-3 strains each and are not
statistically testable at all). No species shows a color↔metabolome
signal. This strengthens the whole-panel Phase 2 null somewhat (a
species-specific signal confined to one of the untested tiny-n species
can't be ruled out, but nothing in the well-sampled tier shows it), while
the small-species caveat above means this is additional *consistent*
evidence, not additional *independently powered* evidence.

## Recommendation (not yet done)
- C\* was not run for the 5 new small species — low expected value given
  a\*'s null and the already-thin power, but cheap to add if useful.
- If any of these small species is of specific biological interest (e.g.
  *R. sphaerocarpa*/*R. taiwanensis* — Idea 5's top convergent-color
  candidates, see `../phase5_genome_linkage/idea5_regime_shift/RESULTS.md`),
  consider whether it's worth waiting for more strains before trusting a
  null there specifically.

## Reproduce
```
for SP in "Rhodotorula dairenensis" "Rhodotorula diobovata" "Rhodotorula taiwanensis" \
          "Rhodotorula sp. clade I" "Rhodotorula sphaerocarpa"; do
  python3 analysis/scripts/phase2_within_species_association.py --species "$SP" --predictor area --min-strains 5
  python3 analysis/scripts/phase2_within_species_association.py --species "$SP" --predictor a     --min-strains 5
done
```
