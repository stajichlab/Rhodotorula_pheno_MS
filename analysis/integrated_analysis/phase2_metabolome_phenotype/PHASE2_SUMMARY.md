# Phase 2: Color <-> Metabolome Association — Results (2026-08-15)

Script: `analysis/scripts/phase2_color_metabolome_association.py`. Method,
predictor choice, and negative-control design are as specified in
`analysis/INTEGRATED_ANALYSIS_STRATEGY.md` Phase 2 (PI-confirmed via the
2026-08-15 grilling session — see `.living/decisions.md`).

- Predictor: a\* (**primary**), C\* (secondary), both from the canonical
  `control_90_110` phenotype table.
- Test set: 10,949 deduplicated MS2 compound groups (`ms_feature_dedup_groups.csv`
  representatives), TSS-normalized within fraction, ~274-275 strains with
  both MS and phenotype data (cell / supernatant fractions separately).
- Inference: empirical p-value from label permutation restricted within
  6 species-tree-derived phylogenetic blocks (200 permutations), BH-FDR
  <0.05 on the deduplicated test count. Asymptotic Spearman p reported as
  a secondary/diagnostic column only.
- **Negative-control hard gate enforced**: `--predictor area` (colony
  size, a phylogenetically-structured but color-unrelated decoy trait)
  must run first and pass freshness checks before `--predictor a`/`C` will
  write output.

## Result 1: no whole-panel color association survives correction

| Predictor | Fraction | n strains | n tested | BH-FDR<0.05 hits | Null-expected hits |
|---|---|---|---|---|---|
| a\* (primary) | cell | 274 | 10,164 | **0** | ~0 |
| a\* (primary) | supernatant | 275 | 10,416 | **0** | ~0 |
| C\* (secondary) | cell | 274 | 10,164 | **0** | ~0 |
| C\* (secondary) | supernatant | 275 | 10,416 | **0** | ~0 |

No deduplicated MS2 compound group reaches BH-FDR<0.05 for either color
axis, in either fraction. This matches the permutation-null expectation
almost exactly (no inflation), which is a *good* sign for the pipeline's
calibration, but means **Phase 3 has no color-associated feature list to
enrichment-test as originally scoped**. This is consistent with the power
ceiling flagged throughout the strategy doc (n≈17-18 independent
species-level lineages) — a whole-panel color→pigment relationship, if it
exists, is not detectable at a formal BH-FDR<0.05 threshold with the
current sample.

**Not evidence that color and pigment chemistry are unrelated** — it is
evidence that this specific test, at this power, does not detect it. See
"Recommended next steps" below.

### Follow-up: within-species tests (higher power, phylogeny-free)

The whole-panel test above pays an "n≈17-18 independent lineages" tax for
phylogenetic correction. Restricting to variation *within* one species
sidesteps that almost entirely. Every species in the panel with ≥5
strains and MS data has now been tested this way — **all 8 come back
null for a\***:

| Species | n | Result | Details |
|---|---|---|---|
| *R. mucilaginosa* | 206 | Null | [`WITHIN_SPECIES_MUCILAGINOSA.md`](WITHIN_SPECIES_MUCILAGINOSA.md) |
| *R. paludigena*, *R. toruloides* | 10 each | Null | [`ROBUSTNESS_AND_MULTIVARIATE.md`](ROBUSTNESS_AND_MULTIVARIATE.md) (also includes sparse multivariate Lasso + pattern-group ANOVA, both null) |
| *R. dairenensis*, *R. diobovata*, *R. taiwanensis*, *R. sp. clade I*, *R. sphaerocarpa* | 5-8 each | Null, **but negative control uninformative at this n** | [`WITHIN_SPECIES_SMALL_SPECIES_SWEEP.md`](WITHIN_SPECIES_SMALL_SPECIES_SWEEP.md) |

Only *R. mucilaginosa* (n=206) has a negative control with real
demonstrated power (area decoy: 1,524/10,164 cell-fraction hits). The
6 smaller species' area decoys all returned 0 hits too, which the sweep
doc above flags explicitly as underpowered nulls, not confirmed ones.

## Result 2: important QC finding — colony area confounds cell-fraction abundances broadly

The negative-control decoy run (`--predictor area`) was expected to show
~0 hits, like the real color runs. Instead:

| Predictor | Fraction | n tested | BH-FDR<0.05 hits | Null-expected hits | Median \|rho\| among hits |
|---|---|---|---|---|---|
| area (decoy) | **cell** | 10,164 | **1,524** | ~0 | 0.23 |
| area (decoy) | supernatant | 10,416 | 0 | ~0 | — |

**1,524 of 10,164 deduplicated cell-fraction compound groups (~15%) show
a statistically robust association with colony area**, entirely absent in
the supernatant fraction. Effect sizes are modest and diffuse (median
\|rho\|≈0.23 among hits, not a few very strong outlier compounds), which
points toward a **broad, systematic technical confound** — most plausibly
that larger colonies yield proportionally more (or less efficiently
extracted) cell material, shifting the TSS-normalized abundance profile of
many compounds at once — rather than a genuine biological signal
concentrated in a handful of pigment-relevant compounds. The fraction
asymmetry (cell only, not supernatant) is consistent with a
biomass/extraction-loading explanation: colony size is a cell-pellet
property, and secreted (supernatant) metabolite composition wouldn't be
expected to scale the same way.

**This finding is more actionable right now than the color result.**
Recommended before further cell-fraction work: investigate whether the
EverythingBagel extraction/normalization protocol scales linearly with
colony biomass, and/or add colony area as a covariate in any future
cell-fraction analysis (Phase 2's a\*/C\* runs already show this isn't
masking a color signal, since the real predictor runs found nothing even
without an area covariate — but a real future finding could still be
partly confounded by this if it correlates with area).

## Recommended next steps

1. **Report Phase 2 honestly as a null result at the pre-registered
   threshold**, not spun as a positive finding. This is a real, useful
   result: it constrains how strong any true color→pigment relationship
   can be, given this sample.
2. **For Phase 3, use nominal (non-FDR-significant) top-ranked
   associations as explicitly exploratory/hypothesis-generating input**,
   not confirmed hits — e.g., the top 20-50 features by empirical p per
   fraction/predictor, clearly labeled as below the pre-registered
   significance threshold. Whether this is worth doing depends on whether
   Phase 3's SIRIUS-based compound-class check can still say anything
   useful about a merely-nominally-ranked list; if not, Phase 3 may need
   to wait for more data (e.g., additional strains, or the improved SIRIUS
   annotations) rather than proceeding on an underpowered feature list.
3. **Investigate the colony-area/cell-fraction confound** as a
   near-term, standalone task — independent of the color hypothesis, this
   affects any future analysis of cell-fraction absolute or TSS-normalized
   abundances.
4. ~~Consider whether **within-species variation tests**...~~ **Done** —
   see the "Follow-up: within-species tests" subsection above; all 8
   testable species (≥5 strains) come back null.

## Reproduce

```
python3 analysis/scripts/phase2_color_metabolome_association.py --predictor area   # negative control, must run first
python3 analysis/scripts/phase2_color_metabolome_association.py --predictor a      # primary
python3 analysis/scripts/phase2_color_metabolome_association.py --predictor C      # secondary
```

Outputs: `color_metabolome_association_{a,C,area_decoy}.csv` in this
directory — one row per (deduplicated feature group x fraction), columns
`spearman_rho, empirical_p, empirical_fdr, asymptotic_p, asymptotic_fdr`
plus leave-one/few-top-orange-species-out sensitivity columns for the top
20 hits per fraction (all NaN here since there were no real hits to
sensitivity-check for a\*/C\*; populated for the area decoy's top hits).
