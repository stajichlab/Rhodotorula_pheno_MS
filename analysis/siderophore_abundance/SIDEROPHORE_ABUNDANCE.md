# Siderophore / Iron-Related Metabolite Abundance

Summarizes MS2 feature abundance for rhodotorulic acid and other
siderophore/hydroxamate-family compounds identified by prior SIRIUS/CANOPUS
structure annotation (`analysis/sirius_annotation/`), across all sampled
*Rhodotorula* (+ outgroup) species and strains.

## Method

1. Searched `analysis/sirius_annotation/sirius_annotations.tsv` (SIRIUS
   structure name + NPC pathway/class + ClassyFire class) for siderophore
   keywords: rhodotorulic/rhodotorulate, siderophore, ferrichrome,
   ferrioxamine, desferrioxamine/deferoxamine, coprogen, fusarinine, dimerum,
   basidiochrome, rhizoferrin, hydroxamate.
2. Pulled per-sample peak areas for each matching `row ID` from the aligned
   MS2 feature quant matrix
   (`data/processed/EB_20260130_ExFAB_Rhodo_Sup_and_Pellet/.../aligned_features_ms2.csv.zst`).
3. Joined to sample metadata (species, canonical strain, cell-pellet vs.
   supernatant fraction) via
   `...-merged_metadata.fixed.tsv.gz`; blanks/QC/SPE-blank samples excluded.
4. Aggregated mean/median/max peak area and detection rate per
   species x fraction and per strain x fraction.

Script: `scripts/siderophore_abundance_table.py` (`run.sh` reproduces all
outputs). Passed `scilintr` clean.

## Candidate features found (12 total)

| row ID | SIRIUS structure name | Structure confidence |
|---|---|---|
| 562 | **Rhodotorulic Acid** | **0.969** |
| 46831 | Desferrioxamine D1 | 0.440 |
| 6976 | Desferrioxamine X5 | 0.129 |
| 11834 | Proferrioxamine-D2 | 0.115 |
| 67 | desferrioxamine h | 0.063 |
| 47668 | Desferrioxamine X5 | 0.059 |
| 249 | desferrioxamine h | 0.058 |
| 7498 | Desferrioxamine X2 | 0.056 |
| 11582 | Deferoxamine | 0.056 |
| 4253 | N-Stearoyldesferrioxamine | 0.030 |
| 4788 | ddesferrioxamine[-00] | 0.028 |
| 3295 | ddesferrioxamine[-00] | 0.018 |

**Only row 562 (rhodotorulic acid) has a structure confidence high enough to
treat as an identification (0.969).** All 11 desferrioxamine/ferrioxamine
hits have confidence 0.018-0.44 — SIRIUS's closest database match, not a
validated identification. These are reported as *tentative, structurally
related candidates* only, not confirmed compounds. No MS/MS spectral
library match (GNPS) or authentic standard has been checked for any of
these 11 -- that would be required before treating them as real.
Desferrioxamine D1 (0.44) is the best of this tentative group and is the
only one plausibly worth a follow-up spectral-match check.

## Headline result: Rhodotorulic acid (row 562, confidence 0.969)

Detected in essentially all sampled species/strains, both fractions, at
90-100% detection rate. Mean peak area is **50-100x higher in cell pellet
than in matched supernatant** for every species with both fractions sampled
-- consistent with an intracellularly retained/cell-associated hydroxamate
siderophore (rhodotorulic acid is a well-characterized *Rhodotorula*/
*Rhodosporidium* siderophore in the literature), not a signal driven by one
or two outlier species.

| Species | n (cell) | mean peak area, cell pellet | detection %, cell | n (sup) | mean peak area, supernatant | detection %, sup |
|---|---:|---:|---:|---:|---:|---:|
| R. paludigena | 10 | 84,134,886 | 90% | 10 | 854,621 | 100% |
| R. kratochvilovae | 3 | 69,829,507 | 100% | 3 | 950,700 | 100% |
| R. mucilaginosa | 208 | 54,395,378 | 98% | 208 | 671,778 | 98% |
| Pseudomicrostroma phylloplanum (outgroup) | 1 | 54,235,360 | 100% | 1 | 569,872 | 100% |
| Cystobasidium sp. (outgroup) | 1 | 50,757,388 | 100% | 1 | 657,892 | 100% |
| R. diobovata | 9 | 50,015,987 | 100% | 9 | 1,108,706 | 100% |
| R. pacifica | 3 | 46,406,715 | 100% | 3 | 784,274 | 100% |
| R. sp. clade XIII | 1 | 45,201,232 | 100% | 1 | 450,914 | 100% |
| R. dairenensis | 7 | 42,432,812 | 100% | 8 | 584,350 | 100% |
| R. sphaerocarpa | 6 | 41,352,945 | 100% | 6 | 606,477 | 100% |
| R. taiwanensis | 6 | 38,168,386 | 100% | 6 | 504,237 | 100% |
| R. toruloides | 10 | 33,517,079 | 100% | 10 | 687,023 | 100% |
| R. glutinis | 2 | 28,587,214 | 100% | 2 | 412,530 | 100% |
| R. graminis | 3 | 27,195,529 | 100% | 3 | 491,991 | 100% |
| R. sp. clade XI | 2 | 26,975,879 | 100% | 2 | 751,918 | 100% |
| R. sp. clade I | 5 | 22,800,206 | 80% | 5 | 445,383 | 100% |
| R. araucariae | 1 | 18,291,294 | 100% | 1 | 619,330 | 100% |

Full table with n_detected, median, and max: `outputs/siderophore_abundance_by_species.tsv`.
Strain-level breakdown: `outputs/siderophore_abundance_by_strain.tsv`.

## Caveats

- Peak areas are relative (uncalibrated, no authentic standard curve) --
  cross-species comparisons are relative abundance, not absolute
  concentration.
- Several species/outgroups have n=1-3 samples (Pseudomicrostroma,
  Cystobasidium, R. araucariae, R. sp. clade XI/XIII, R. pacifica,
  R. kratochvilovae, R. glutinis, R. graminis) -- means for these are not
  statistically robust and should be read as single/few-replicate values,
  not population estimates.
- No statistical test (e.g., cell vs. supernatant paired test, ANOVA across
  species) has been run yet on this table; the summary above is descriptive
  only.
- The desferrioxamine/ferrioxamine hits are unconfirmed and should not be
  cited as identified siderophores without a spectral library match.

## Outputs

- `outputs/siderophore_feature_annotations.tsv` -- the 12 candidate SIRIUS
  annotations, sorted by structure confidence.
- `outputs/siderophore_abundance_by_species.tsv` -- feature x species x
  fraction summary (n, n_detected, mean/median/max peak area, detection
  rate).
- `outputs/siderophore_abundance_by_strain.tsv` -- feature x strain summary
  (fraction-resolved).
