# Sterol/ergostane cluster vs. copper-resistance AUC (naive quick-check, 2026-08-16)

Follow-up to `EXPANDED_SEARCH_RESULTS.md`. SIRIUS's NPC-class breakdown of
the 73 "Terpenoids"-pathway features (asked by PI: "terpene composition of
the mass spec molecules?") surfaced 7 features independently classed as
Ergostane/Cholane/Cholestane steroids -- the single largest NPC subclass in
that bucket, and structurally consistent with the row-846 ergosterol ISF
candidate that already had a naive copper-AUC lead. Ran the same naive
Spearman quick-check (`analysis/scripts/idea1_auc_quickcheck.py`, promoted
from the earlier one-off inline script per that doc's recommendation) on
the 3 most chemically self-consistent SIRIUS ergostane calls plus row 846
as an anchor.

**Still explicitly a screen, not a validated result**: no phylogenetic
block permutation, no negative control, no BH-FDR -- same caveat as
`EXPANDED_SEARCH_RESULTS.md`.

## Candidates tested

| Row ID | SIRIUS call | Confidence | Dedup group |
|---|---|---|---|
| 846 | (ergosterol ISF, mass-search candidate, not itself SIRIUS-annotated) | -- | asid_31011 (anchor, previously tested) |
| 9852 | Peroxyergosterol | 0.694 | asid_9852 |
| 6682 | Ergost-3,5,7,9(11),22-pentaen | 0.834 | asid_6682 |
| 35014 | 7-Hydroxyergosterol | 0.337 | asid_3176 |

## Results (naive Spearman, TSS-normalized, canonical_strain-averaged, n~267-268)

| Row ID | Fraction | n | rho | p (naive) |
|---|---|---|---|---|
| 846 | cell | 267 | 0.218 | 0.0003 |
| 846 | supernatant | 268 | -0.076 | 0.21 |
| 9852 | cell | 267 | 0.169 | 0.0055 |
| 9852 | supernatant | 268 | -0.007 | 0.91 |
| 6682 | cell | 267 | **0.278** | **4.0e-6** |
| 6682 | supernatant | 268 | **-0.240** | **7.0e-5** |
| 35014 | cell | 267 | 0.226 | 0.0002 |
| 35014 | supernatant | 268 | 0.258 | 1.8e-5 |

(846's rho differs slightly from the earlier one-off number, 0.218 vs.
0.226 -- within noise of a minor strain-averaging implementation detail
between the old inline script and the new promoted one; not a discrepancy
worth chasing.)

## Interpretation

All four independently SIRIUS/mass-called ergostane-class features show
the **same direction** of naive correlation with copper-resistance AUC in
the cell fraction (rho 0.17-0.28, all p<0.006). This is materially
stronger circumstantial evidence than the single row-846 result alone: if
this were pure noise, four separate features drawn from the same
biosynthetic family would not be expected to agree in sign this
consistently. Two additional details are worth carrying forward:

- **Row 6682 flips sign between fractions** (cell +0.28, supernatant
  -0.24, both significant) -- a plausible intracellular-retention-vs-
  secretion pattern for a membrane sterol, not just fraction noise.
- **Row 35014 is positive in BOTH fractions** -- differs from the other
  three (cell-only signal), so not a perfectly uniform pattern across the
  cluster; worth keeping distinct rather than treating the cluster as one
  monolithic signal.

## Status / next step (done 2026-08-17)

The full-rigor test is now done: added `--predictor auc` (mean_auc_rate from
`sample_metadata.csv.gz`) to `phase2_within_species_association.py` and ran it
on *R. mucilaginosa* (207 strains with copper AUC), gated by a freshly rerun
`area` colony-size decoy (same hard-gate convention as the color tests).

**Result: the naive whole-panel sterol-AUC signal does NOT survive**
**within-species phylogenetic block permutation.** All four cluster features
collapse to |rho| <= 0.12 within *R. mucilaginosa* (all empirical_p > 0.10),
and the overall within-species AUC run finds 0 BH-FDR<0.05 hits in either
fraction (null mean 0.0 hits/permutation).

This is the third straight failure of a naive whole-panel correlation to
survive phylogenetic/within-species rigor (F-005 in
`.living/findings/phylogenetic-confounding-of-trait-molecular-associations.md`),
and the second for the copper-AUC phenotype specifically (the naive amino-acid
hits failed the *R. mucilaginosa*-only check first). Consistent with the
established finding, the whole-panel cell-fraction rho 0.17-0.28 was driven by
between-species phylogenetic structure (copper-resistant species carry more
cell sterol), not by within-species co-variation.

Outputs:
- `analysis/integrated_analysis/phase2_metabolome_phenotype/within_species_Rhodotorula_mucilaginosa_association_auc.csv`
- `.../within_species_Rhodotorula_mucilaginosa_association_area_decoy.csv` (refreshed decoy)
- Reproduce: `python3 analysis/scripts/phase2_within_species_association.py --species "Rhodotorula mucilaginosa" --predictor auc`

A species-level (whole-panel, species-tree-blocked) variant --
`phase2_color_metabolome_association.py` machinery with mean_auc_rate -- is
still available if we ever want the species-lineage-level view, but it is the
weaker design (effective sample size bounded by ~17 species lineages) and the
within-species result already tells the story; not run this pass.

## Reproduce

```
python3 analysis/scripts/idea1_auc_quickcheck.py --row-id 846 9852 6682 35014
```
Output: `sterol_cluster_auc_quickcheck.csv` in this directory.
