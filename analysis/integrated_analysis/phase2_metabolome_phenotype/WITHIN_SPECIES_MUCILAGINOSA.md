# Within-species color <-> metabolome association: *R. mucilaginosa* (2026-08-15)

Script: `analysis/scripts/phase2_within_species_association.py`. Follows
directly from PHASE2_SUMMARY.md's "Recommended next steps" #4: the
whole-panel test's effective sample size is bounded by ~17-18 independent
species-level lineages; a single species with many strains sidesteps that
almost entirely. *R. mucilaginosa* is the obvious candidate here: 216
phenotyped strains, 206 with both color and MS data, 201 with genome data
for real strain-level phylogenetic blocking (not the coarser species-tree
blocks the whole-panel script uses).

## Why this is a meaningfully different test, not just "Phase 2 on a subset"

- **Much higher effective n.** 206 strains, blocked into 27 clades from
  the actual strain-level genome tree (vs. 6 species-level blocks for the
  whole panel) — the phylogenetic-independence tax is far smaller.
- **Real color variation to explain.** *a\** within *R. mucilaginosa*
  alone: mean 10.26, sd 1.78, range **0.60–14.21** — essentially the same
  range as the entire 303-strain panel (0.60–14.21). This species alone
  spans almost the whole orange-color spectrum in this dataset. A null
  result here is not explained by "not enough within-species color
  variation to correlate with anything."
- **CORRECTION (2026-08-15, later same day)**: this section originally
  reported the colony-area decoy as "0 hits" within *R. mucilaginosa*,
  based on a `--n-perm 20` smoke-test run that was mistakenly never
  rerun at the full `--n-perm 200` before the hard gate was treated as
  satisfied (the gate checks file existence/freshness, not permutation
  count — a real gap). **A properly-powered rerun (`--n-perm 200`,
  `--seed 0`, deterministic) shows the area confound DOES reproduce
  within *R. mucilaginosa***: 2,025/9,437 cell-fraction compound groups
  significant at BH-FDR<0.05 (0/10,255 supernatant) — not meaningfully
  different in magnitude from the whole-panel scan's 1,524/10,164. The
  "area confound is purely between-species" interpretation below is
  **wrong** and superseded — see the corrected table and interpretation
  further down. The a\*/C\* real-predictor results (both run at the full
  `--n-perm 200` from the start) are unaffected by this error and remain
  valid.

## Result: color null holds; area decoy corrected

| Predictor | Fraction | n strains | n tested | BH-FDR<0.05 hits | Null-expected |
|---|---|---|---|---|---|
| a\* (primary) | cell | 206 | 9,437 | **0** | ~0 |
| a\* (primary) | supernatant | 206 | 10,255 | **0** | ~0 |
| C\* (secondary) | cell | 206 | 9,437 | **0** | ~0 |
| C\* (secondary) | supernatant | 206 | 10,255 | **0** | ~0 |
| area (decoy), corrected | cell | 206 | 9,437 | **2,025** | ~0 |
| area (decoy), corrected | supernatant | 206 | 10,255 | **0** | ~0 |

The corrected decoy result actually **strengthens** confidence in the
color null rather than weakening it: color (a\*/C\*) still shows zero
association even though this exact pipeline, on this exact strain set, is
demonstrably capable of detecting a strong effect (colony area) when one
is present. The area confound is present both between- and within-species
— its cause remains an open question (see
`.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md`,
to be updated) rather than a purely between-species artifact as
originally (incorrectly) concluded.

Top nominal (uncorrected) hits for a\* have modest effect sizes
(|Spearman rho| ≈ 0.18-0.21, empirical FDR ≈ 0.48-0.48 — nowhere near
significant) — a genuine null, not a near-miss obscured by multiple
testing.

## Interpretation

This substantially strengthens the Phase 2 whole-panel null: with much
better power (206 strains, real strain-level phylogenetic blocking, and
color variation spanning almost the full observed range), **no
metabolome feature tracks color within the largest, best-sampled species
in this panel either.** Combined with the whole-panel result, the
evidence increasingly suggests that if colony a\*/C\* is linked to
specific MS2-detectable compounds at all, the relationship is either (a)
too weak to detect with current sample sizes even at within-species
power, (b) present only in specific other species/clades not tested here,
(c) driven by compounds not well captured by this untargeted LC-MS2 method
(e.g. compounds that don't ionize well, are below detection, or aren't
chromatographically resolved), or (d) not actually present as a simple
linear/monotonic univariate relationship (a more complex, multivariate, or
threshold-based relationship could still exist and evade a per-feature
correlation test).

## Recommended next steps

1. **Try one or two other species with enough strains** for a similar
   test as a robustness check (e.g. *R. paludigena* n=17, *R. toruloides*
   n=10, *R. diobovata* n=10) — smaller n, weaker power, but worth a quick
   look given how easy this script makes it (`--species` flag).
2. **Consider a multivariate approach**: instead of per-feature univariate
   correlation, test whether color axes are jointly predictable from a
   multivariate combination of features (e.g. sparse PLS/regularized
   regression) — could detect a real relationship that's diffuse across
   many weakly-correlated features rather than concentrated in one strong
   hit, which a per-feature test with FDR correction is specifically bad
   at detecting.
3. **Revisit whether untargeted LC-MS2 is even capturing the relevant
   compounds** — carotenoids are chemically difficult (poor ionization,
   often below LOD in standard reversed-phase LC-MS), and this whole
   analysis assumes the metabolome data contains a detectable pigment
   signal. Targeted analysis (the HPLC/LC-UV-Vis work, contingent on the
   SIRIUS re-run per the 2026-08-15 grilling session decision) would
   resolve this ambiguity directly rather than continuing to test an
   assumption the untargeted data may not be positioned to confirm.

## Reproduce

```
python3 analysis/scripts/phase2_within_species_association.py --species "Rhodotorula mucilaginosa" --predictor area   # negative control, must run first
python3 analysis/scripts/phase2_within_species_association.py --species "Rhodotorula mucilaginosa" --predictor a      # primary
python3 analysis/scripts/phase2_within_species_association.py --species "Rhodotorula mucilaginosa" --predictor C      # secondary
```
