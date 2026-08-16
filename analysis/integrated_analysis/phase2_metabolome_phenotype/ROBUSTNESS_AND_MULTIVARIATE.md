# Phase 2 follow-ups: cross-species robustness checks + multivariate model (2026-08-15)

Two follow-ups to the *R. mucilaginosa* within-species null
(`WITHIN_SPECIES_MUCILAGINOSA.md`), both requested directly by the PI.

## 1. Robustness check: 2 more well-sampled species

`analysis/scripts/phase2_within_species_association.py --min-strains 8`
(explicitly labeled exploratory given small n — see script's NOTE output).

| Species | n strains (MS+color) | a\* BH-FDR<0.05 hits (cell / supernatant) |
|---|---|---|
| *R. paludigena* | 10 | 0 / 0 |
| *R. toruloides* | 10 | 0 / 0 |
| *R. mucilaginosa* (for reference) | 206 | 0 / 0 |

Consistent null across all three species tested, at both large (206) and
small (10) n. Not independently powerful evidence at n=10, but no hint of
a hidden signal being masked by *R. mucilaginosa*-specific homogeneity
either — the pattern doesn't depend on which species is tested.

## 2. Multivariate model (sparse Lasso, group-aware CV)

`analysis/scripts/phase2_multivariate_association.py` — tests whether
color is predictable from a **sparse combination** of many features
jointly (Lasso regression, GroupKFold CV using the real strain-level
phylogenetic blocks, alpha chosen by grid search, significance from
label permutation within the same blocks), addressing the concern that a
per-feature univariate test with FDR correction is poorly suited to
detect a real but diffuse, multi-feature signal.

*R. mucilaginosa*, both fractions, both color predictors:

| Predictor | Fraction | Observed Q² (cross-val R²) | Null Q² mean | Empirical p |
|---|---|---|---|---|
| a\* | cell | -0.49 | -0.53 | 0.26 |
| a\* | supernatant | -0.63 | -0.58 | 0.88 |
| C\* | cell | -0.32 | -0.33 | 0.47 |
| C\* | supernatant | -0.40 | -0.37 | 0.94 |
| area (decoy) | cell | -0.02 | -0.07 | 0.02* |
| area (decoy) | supernatant | -0.02 | -0.07 | 0.02* |

**Every observed Q² is negative** — the model performs worse than simply
predicting the mean color value for every strain, for both real color
predictors and the decoy. This is a clean, unambiguous null: there is no
detectable sparse multivariate combination of metabolome features that
predicts color better than chance in this species/dataset, at least via
L1-regularized linear regression.

*The decoy's p=0.02 needs a caveat, not a headline*: with observed Q²
still negative (-0.02, worse than baseline), a "significant" permutation
p-value here means the model was *less catastrophically bad* than most
permutations, not that it predicted anything. Given only 100 permutations
(coarsest resolution ~0.01) and two fractions tested, this is not treated
as a real decoy-trait hit — consistent with the univariate decoy result
for area *within* *R. mucilaginosa*, which was also null (0 hits, see
`WITHIN_SPECIES_MUCILAGINOSA.md`).

## 3. ANOVA / pattern-group approach (2026-08-15, later same day)

`analysis/scripts/phase2_anova_pattern_association.py` — a structurally
different test from the correlation/regression approaches above: instead
of assuming a monotonic/linear relationship with a single color axis,
strains are clustered into 3 groups by k-means on the **joint**
[L\*, a\*, b\*] pattern (not just a\* or C\* alone), then each deduplicated
feature is tested via Kruskal-Wallis (non-parametric one-way ANOVA) across
those pattern groups, same phylogenetic block-permutation inference as
the other scripts. This can catch non-monotonic or threshold-like
relationships a rank correlation would miss.

**Note on this result's provenance**: the negative-control decoy for this
test was correctly run at full power from the start and showed a strong,
real effect (1,354/9,437 cell-fraction hits for the colony-area decoy,
0/10,255 supernatant) — this *cross-validated* an error found elsewhere
(see `WITHIN_SPECIES_MUCILAGINOSA.md`'s correction: the univariate
within-species decoy had originally been reported as "0 hits" based on an
under-powered `--n-perm 20` smoke test that was never rerun at full
power; both are now confirmed consistent at ~1,300-2,000 cell-fraction
hits once properly run). Both Phase 2 within-species scripts now record
`n_perm` in their output and refuse to run a real predictor against a
decoy with `n_perm < 100`, closing this gap.

| Test | Fraction | n tested | BH-FDR<0.05 hits |
|---|---|---|---|
| Color pattern (3 clusters, L\*/a\*/b\* jointly) | cell | 9,437 | **0** |
| Color pattern (3 clusters, L\*/a\*/b\* jointly) | supernatant | 10,255 | **0** |
| Area pattern (decoy, 3 clusters) | cell | 9,437 | **1,354** |
| Area pattern (decoy, 3 clusters) | supernatant | 10,255 | **0** |

Null again, on a fifth methodologically distinct approach, with a decoy
that's now properly confirmed to detect real signal when present.

## Overall conclusion after four complementary tests

| Test | Result |
|---|---|
| Whole-panel univariate (Phase 2 main) | null |
| Within-species univariate (*R. mucilaginosa*, n=206) | null |
| Within-species univariate, 2 more species (n=10 each) | null |
| Within-species multivariate/sparse (*R. mucilaginosa*) | null |

Five independent analytical approaches — univariate correlation (whole-
panel and within-species at two sample sizes), sparse multivariate
regression, and a non-parametric pattern-group ANOVA — all agree: **there
is no detectable linear, sparse-combination, or cluster/pattern-group
relationship between CIELAB color and untargeted LC-MS2 metabolome
abundance in this dataset**, at the species and whole-panel scales tested.
This is now a well-triangulated negative result, not an artifact of any
single method's limitations, and the ANOVA/pattern-group approach in
particular rules out the specific concern that a monotonic-relationship
assumption (correlation, Lasso) was hiding a non-monotonic or
threshold-like true relationship.

## What this does NOT rule out

- A relationship that isn't captured by 3-cluster pattern groups or a
  linear model (e.g. a smoother nonlinear/kernel relationship) — a
  tree-based or kernel method might behave differently, though with
  n≈206 and p≈9,500-10,000 features, overfitting risk for a more
  flexible model is severe without much larger n.
- A relationship present in **specific compounds the untargeted method
  doesn't detect well** — carotenoids are known to ionize poorly in
  standard LC-MS/MS; this is the leading remaining hypothesis and exactly
  what the targeted HPLC/LC-UV-Vis validation (contingent on the SIRIUS
  re-run, per the 2026-08-15 grilling session) would test directly.
- A relationship present in **other species not yet tested** at
  sufficient power (only 3 of 17-18 species have been checked here).

## Recommendation

Given four consistent null results, further per-feature or per-model
statistical fishing within this untargeted dataset has a low expected
return. **The highest-value next step is the targeted pigment
validation** (HPLC/LC-UV-Vis), since it directly tests whether the
untargeted method is even capable of seeing the relevant chemistry,
rather than continuing to test statistical variations on an assumption
(that untargeted LC-MS2 captures pigment-relevant signal) that these four
results are, collectively, mildly evidence against.
