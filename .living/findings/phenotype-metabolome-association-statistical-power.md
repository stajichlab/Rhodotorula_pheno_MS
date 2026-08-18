---
topic: phenotype-metabolome-association-statistical-power
description: Whether whole-panel correlation between a quantitative organismal phenotype and untargeted metabolome features can reach formal significance given a strain panel with modest phylogenetic independence (small number of species-level lineages despite many strains).
created: 2026-08-15
last_updated: 2026-08-15
status: active
---

# Phenotype-metabolome association statistical power

## F-001: No whole-panel colony-color <-> metabolome association survives phylogenetically-aware correction at n~275 strains / ~17-18 species
**Status:** preliminary
**Claim:** Across ~275 fungal strains (17-18 species) with paired colony
color (CIELAB a*, C*) and untargeted LC-MS2 metabolome data, a
whole-panel Spearman correlation test (10,949 deduplicated compound
groups per fraction, phylogenetically-block-restricted permutation
p-values, BH-FDR<0.05) found **zero** significant color-metabolome
associations for either color axis, in either fraction (cell,
supernatant) — matching the ~0 hits expected under the permutation null
almost exactly (no inflation, but also no signal detected).
**Implications:** A whole-panel design across many strains does not
automatically confer high statistical power for phylogenetically-aware
tests if the number of *independent* lineages (species-level tree tips)
stays small — here, ~275 strains reduce to an effective sample closer to
17-18 for the block-permutation null once phylogenetic non-independence
is accounted for. This is a concrete demonstration that "more strains"
is not the same as "more power" for this class of test; the true limiting
factor is the number of independent evolutionary origins/lineages
sampled, not the total strain count. A null result at this threshold does
not rule out a real color-pigment relationship — it means this specific
test, at this power, cannot detect one weaker than what ~17-18
independent points can resolve at alpha=0.05 with FDR correction across
~10,000 tests.
**Tags:** statistical-power, phylogenetic-correction, metabolomics, phenotype-association, permutation-test, null-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | Phase 2 color-metabolome association (a\*, C\* runs) | EB_20260130_ExFAB_Rhodo_Sup_and_Pellet aligned_features_ms2 (10,949 deduplicated compound groups), control_90_110 phenotype table | Rhodotorula_pheno_MS | 0/10,164 (cell) and 0/10,416 (supernatant) compound groups significant at BH-FDR<0.05 for a\*; same for C\*; permutation null also ~0 hits (well-calibrated, not underpowered-looking-like-a-bug) | supports |

### Open Questions
- Would relaxing to a nominal/exploratory threshold (not FDR-corrected)
  surface a top-ranked candidate list worth pursuing with targeted
  validation, even without formal significance?
- Does this null result hold up if compound groups are aggregated further
  (e.g. to NPC-pathway-level summed abundance) rather than tested
  individually, trading specificity for a smaller, less severely
  FDR-penalized test set?
- Would a multivariate approach (sparse PLS / regularized regression
  jointly across features) detect a diffuse relationship a per-feature
  univariate test with FDR correction is poorly suited to find?

## F-006: Extreme-group (quartile) color contrast also null — 6th independent method to agree
**Status:** preliminary
**Claim:** `analysis/scripts/extreme_group_color_association.py` split
strains into top/bottom quartile groups on a\*/C\*/area and tested for
compound-abundance differences via a block-permuted rank-sum statistic —
a design built to catch threshold/nonlinear effects a continuous
correlation could miss, at the acknowledged cost of lumping species
together (PI-flagged caveat). Null for both a\* and C\*, both fractions
(0 BH-FDR<0.05 hits); the area decoy shows 1,723 cell-fraction hits,
confirming real power at this n.
**Implications:** A 6th independent statistical method (after the 5 in
F-002 through F-004) finds no color-metabolome signal, now including one
specifically designed to be more sensitive to threshold effects than
Spearman correlation. Strengthens confidence this is a genuine null, not
a methodological blind spot.
**Tags:** extreme-group, rank-sum, null-result, quartile-split

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Extreme-group test (extreme_group_color_association.py) | Same feature matrix as F-002/F-003 | Rhodotorula_pheno_MS | Null for a\*/C\*, both fractions; decoy well-calibrated | refines |

## F-005: All 8 species with ≥5 strains now tested within-species for color↔metabolome association — all null
**Status:** preliminary
**Claim:** Extended the within-species test (F-003's design) from 3 to
all 8 species in the panel with ≥5 strains and MS data:
*R. dairenensis* (n=8), *R. diobovata* (n=8), *R. taiwanensis* (n=6),
*R. sp. clade I* (n=5), *R. sphaerocarpa* (n=5), in addition to the
already-tested *R. mucilaginosa* (n=206), *R. paludigena* (n=10),
*R. toruloides* (n=10). All 8 come back null for a\*↔metabolome at
BH-FDR<0.05.
**Implications:** No species-specific color↔metabolome signal in the
well-sampled tier of the panel. **Caveat**: the negative-control decoy
(colony area) also returned 0 hits for all 5 new small species (n=5-8),
unlike *R. mucilaginosa*'s decoy which showed real signal (1,524/10,164
hits) -- at n=5-8, ~48-51% of cell-fraction features are excluded as
constant before testing even begins, so a 0-hit decoy at this n does not
carry the same evidential weight. These 5 results should be read as
"no signal detected in a test that may not have been able to detect one,"
not independently confirmed nulls, per this project's established
power-ceiling pattern (see F-002/F-003 in this file).
**Tags:** within-species, null-result, small-species, negative-control-caveat, power-ceiling

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Small-species sweep (phase2_within_species_association.py --min-strains 5) | Same feature matrix as F-002/F-003 | Rhodotorula_pheno_MS | 5 new species, all null; decoys uninformative at this n | refines |

## F-002: Within-species test (R. mucilaginosa, n=206, much higher power) also finds no color-metabolome association
**Status:** supported (2 consistent entries: whole-panel test F-001 above, and this species-restricted follow-up)
**Claim:** Restricting to *R. mucilaginosa* alone (216 phenotyped strains,
206 with paired color+MS data, 201 with genome data for real strain-level
phylogenetic block permutation across 27 blocks) — a design that sidesteps
almost all of the whole-panel test's "effective n bounded by ~17-18
species" limitation — still found **zero** BH-FDR<0.05 color-metabolome
associations for a\* or C\*, in either fraction. This is not a
low-power-due-to-low-variance artifact: a\* within this one species alone
spans nearly the full range seen across the entire 303-strain panel
(0.60-14.21 vs. panel-wide 0.60-14.21). Top nominal (uncorrected) hits had
modest effect sizes (|rho|~0.18-0.21) far from significance — a clean
null, not a near-miss.
**Implications:** The null result is not primarily explained by the
whole-panel test's species-level power ceiling — a much better-powered,
essentially phylogeny-light design on the single largest, most
color-diverse species in the panel reached the same conclusion. This
shifts the more likely explanations toward: the untargeted LC-MS2 method
may not be capturing the relevant pigment compounds well (carotenoids are
known to ionize poorly), the true relationship may be multivariate/diffuse
rather than concentrated in single features a per-feature test can detect,
or there may genuinely be little to no direct color-metabolome linkage at
the resolution this data can measure. **CORRECTION (2026-08-15, later
same day)**: the decoy-trait claim below was wrong — see corrected
Evidence Ledger row. The color null itself is unaffected and, if
anything, better supported: this exact pipeline demonstrably detects a
strong effect (colony area) on this exact strain set when one is
present, and still finds nothing for color.
**Tags:** statistical-power, within-species, null-result, mucilaginosa, phylogenetic-correction

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | Phase 2 within-species association (R. mucilaginosa, a\*/C\* runs) | Same feature/phenotype data as F-001, restricted to R. mucilaginosa strains | Rhodotorula_pheno_MS | 0/9,437 (cell) and 0/10,255 (supernatant) compound groups significant at BH-FDR<0.05 for a\*; same for C\* | supports |
| 2026-08-15 (correction, same day) | Rerun of the decoy area trait at proper `--n-perm 200` (original claim used an underpowered `--n-perm 20` smoke test never rerun before the hard gate was treated as satisfied) | Same as above | Rhodotorula_pheno_MS | Area decoy actually shows 2,025/9,437 (cell) and 0/10,255 (supernatant) hits — the confound DOES reproduce within *R. mucilaginosa*, contradicting the original "between-species only" claim, which is retracted | refines |

### Open Questions
- Untargeted LC-MS2 may simply not resolve carotenoids well (poor
  ionization) — would a targeted pigment method (HPLC/LC-UV-Vis, per the
  2026-08-15 grilling session's contingent-on-SIRIUS decision) find a
  relationship this untargeted approach cannot see?
- Does the same within-species null hold for other species with enough
  strains? (Partially answered by F-004 below: yes for 2 more species,
  small n.)

## F-004: Robustness check (2 more species) and a sparse multivariate model both also find no color-metabolome signal
**Status:** robust (4 consistent entries across different methods: whole-panel univariate F-001, within-species univariate F-002/F-003, 2 more species, and a group-aware Lasso multivariate model — same dataset/project but methodologically independent approaches)
**Claim:** Extending the *R. mucilaginosa* within-species null (F-002/F-003)
in two directions, both null: (a) the same univariate within-species test
in *R. paludigena* (n=10) and *R. toruloides* (n=10) found 0 hits each,
consistent with mucilaginosa despite much smaller n; (b) a sparse
multivariate model (Lasso regression, GroupKFold cross-validation using
real phylogenetic blocks, permutation significance) testing whether color
is predictable from a *combination* of many weakly-correlated features
(addressing the concern that per-feature FDR-corrected tests can miss a
diffuse signal) found cross-validated R² **negative** for both a\* and C\*
in both fractions within *R. mucilaginosa* (observed Q² -0.32 to -0.63,
i.e. worse than predicting the mean; empirical p 0.26-0.94) — no signal
detected by this method either.
**Implications:** Four methodologically distinct analyses (whole-panel
univariate, within-species univariate at two sample sizes, within-species
multivariate/sparse) now agree: no detectable linear or sparse-combination
relationship between CIELAB color and untargeted LC-MS2 metabolome
abundance in this dataset. This substantially de-risks the conclusion that
the null is a real property of this dataset/question rather than an
artifact of any single method (FDR-correction stringency, per-feature test
insensitivity to diffuse signal, or insufficient power). The most likely
remaining explanations are that untargeted LC-MS2 doesn't capture the
relevant pigment chemistry well (carotenoids ionize poorly), or that fewer
than 3-4 of the 17-18 species tested happen to carry the informative
signal. Recommendation: further statistical exploration of this untargeted
dataset has low expected marginal return; targeted pigment validation
(HPLC/LC-UV-Vis) is the higher-value next step.
**Tags:** statistical-power, multivariate, lasso, null-result, robustness-check, mucilaginosa, paludigena, toruloides

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | Within-species robustness check (R. paludigena, R. toruloides, n=10 each) | Same feature/phenotype data as F-002, restricted per species | Rhodotorula_pheno_MS | 0 BH-FDR<0.05 hits for a\* in either species, either fraction | supports |
| 2026-08-15 | Multivariate (Lasso + GroupKFold + permutation) association, R. mucilaginosa | Same as F-002/F-003 | Rhodotorula_pheno_MS | Observed cross-val R² negative for a\*/C\* in both fractions; none beat permutation null | supports |

### Open Questions
- Would a nonlinear/kernel method behave differently, or is overfitting
  risk at this n/p ratio prohibitive for anything more flexible than
  Lasso?
- Does the null hold in species not yet tested (only 3/17-18 checked)?
- Would the targeted HPLC/LC-UV-Vis validation (contingent on SIRIUS
  re-run) find a relationship this untargeted approach structurally
  cannot see?

## F-007: SIRIUS class aggregation and GSEA-style enrichment also find no color class association
**Status:** supported
**Claim:** `analysis/scripts/class_level_association.py` aggregates the
10,949 dedup features into SIRIUS chemical classes (NPC pathway, NPC
class, ClassyFire class — only the ~29%-annotated subset) and tests a\*/C\*
via (strategy 2) class-mean z-score Spearman and (strategy 3) GSEA-style
KS enrichment ES, both against block-restricted (species-tree-clade)
500-permutation nulls with BH-FDR. Null for both color axes: 0-1 of 502
classes survive FDR<0.05 (a\*: 1 NPC-pathway hit in supernatant; C\*: 0);
top a\* signals are single-feature classes (n_members=1) that don't
survive. This is the 7th/8th independent method to agree null (after the
5 in F-002/F-004, F-006, F-005's sweep).
**Why aggregation doesn't rescue signal:** SIRIUS classes in this data are
mostly tiny (typical n_members=1), so class aggregation neither pools weak
within-class signal nor materially reduces multiple testing.
**Decoy/calibration angle:** the area decoy saturates (129/502 classes,
43% of cell/NPC classes) — traced to the F-001 biomass/global-abundance
artifact, `spearman(mean z-abundance, area) = -0.30` (p=3e-7) in the cell
fraction — not a machine bug: a globally shuffled true-null phenotype
through the identical pipeline gives the expected ~3-5% FPR, and
per-species blocks barely reduce decoy hits (43%→41%), so area's signal is
within-species scalar abundance scaling, not species-tree structure. The
default 6-clade restriction is weak in the cell fraction (clade_1 holds
230/274 strains).
**Tag:** class-aggregation, gsea, null-result, biomass-artifact
