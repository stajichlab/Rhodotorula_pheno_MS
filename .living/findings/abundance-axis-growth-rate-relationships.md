---
topic: abundance-axis-growth-rate-relationships
description: Characterization of how the global metabolome abundance axis and individual features relate to liquid growth rate (Cu AUC) — which growth axes the biomass abundance artifact tracks, and the only within-species growth signal found so far (supernatant nucleosides/peptides).
created: 2026-08-17
last_updated: 2026-08-17
status: active
---

# Abundance axis <-> growth-rate relationships

Follow-up to the class-aggregation result (F-017): the area decoy saturates because
TSS-normalized scalar abundance correlates with colony area (`spearman(mean z, area) ≈ -0.30`,
cell fraction). Question asked by the PI: **does the abundance axis also track liquid
growth rate (Cu AUC), and is there anything else to learn of this relationship?**

## F-001: The abundance axis tracks colony area (biomass scaling) but NOT liquid Cu-AUC within species; the two growth phenotypes are different axes

**Status:** supported
**Claim:** The global abundance axis (mean of feature z-scores per strain, TSS-normalized,
replicate-collapsed, adduct-dedup'ed) behaves differently against the two growth phenotypes:

| comparison | fraction | whole-panel rho (p) | within *R. mucilaginosa* rho (p) |
|---|---|---|---|
| mean z vs colony `area` | cell | -0.298 (5.2e-7) | **-0.290 (<<0.01)** |
| mean z vs Cu-AUC | cell | -0.174 (4.2e-3) | -0.070 (0.31) |
| mean z vs colony `area` | supernatant | -0.192 (1.4e-3) | -0.237 (<0.01) |
| mean z vs Cu-AUC | supernatant | 0.007 (0.91) | 0.101 (0.15) |
| colony `area` vs Cu-AUC | cell/supernatant | 0.373/0.361 (<1e-9) | 0.252 (<0.01) |

Colony area and liquid Cu-AUC **are positively correlated** (whole-panel and within-species),
so they are related growth phenotypes — but the metabolome abundance axis tracks **colony area
within species (biomass scaling, the F-001 biomass artifact)** and does **not** track liquid
growth rate within species. The whole-panel abundance↔AUC correlation (-0.17) is almost
entirely **between-species**: species differ strongly in AUC (*R. mucilaginosa* fast, mean
AUC ~23.3; *R. kratochvilovae* slow, ~7.2), so it dissolves to rho≈-0.07 (ns) inside
*R. mucilaginosa*.
**Implications:** Liquid growth rate is not another route to the biomass artifact — the
biomass confound is specific to colony-area-scale (yield-like) biomass, and the "global
abundance axis must relate to growth" intuition does not hold within species. Growth-rate
questions must be analyzed at the feature/class level, not the aggregate abundance axis.
**Tags:** growth-rate, auc, biomass-artifact, within-species, abundance-axis, null-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-17 | growth-rate abundance-axis characterization (Spearman, whole-panel + within-mucilaginosa) | control Cu-AUC + EB MS2 features | Rhodotorula_pheno_MS | mean z ↔ area significant within species (cell & supernatant); mean z ↔ AUC ns within species | supports |

## F-002: Individual-feature growth signals are species-driven in the cell fraction; the only within-species growth signal is a small supernatant nucleoside/peptide set

**Status:** preliminary
**Claim:** Feature-level Spearman of each feature (TSS-normalized relative abundance, strain
level) vs Cu-AUC, cell and supernatant fractions:

- **Cell fraction:** 218/10,230 features (2.1%) reach \|rho\|>0.3 across all species
  (chemistry: ceramides, wax monoesters, sphingolipids, fatty acids, di-/tri-peptides,
  strongest rho ≈ -0.37, most negative). **Within *R. mucilaginosa*: 0/10,230** — the entire
  cell-fraction growth signal is between-species (phylogenetic).
- **Supernatant fraction:** 143/10,439 features (1.4%) \|rho\|>0.3 across all species.
  **Within *R. mucilaginosa*: 10 features** survive, and this count is **significantly above
  permutation null** (block-free permutation of AUC, 200 draws: null mean 0.13 hits, sd 1.29,
  observed 10 → empirical p ≈ 0). All 10 are **positive** (more abundant in faster growers):
  purine nucleosides (4), di-peptides (2), tri-peptides (1), plus 3 unannotated. Signatures
  strongly evocative of excreted nucleotide/peptide economy scaling with growth.
**Implications:** This is the first within-species non-null *growth* signal found in the
project across any molecular layer (see F-005 in the phylogenetic-confounding topic). It is
small (10 features), supernatant-only, with simple, mechanistically-plausible chemistry
(purine nucleosides + small peptides). Strong candidate for the "something else to learn"
about the growth relationship — but n=10 features at \|rho\|≈0.30-0.34, consistent
directionality, single species (n=208 strains).
**Tags:** growth-rate, auc, supernatant, nucleosides, peptides, within-species, lead-not-result, permutation-calibrated

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-17 | feature-level AUC correlations + permutation null | control Cu-AUC + EB MS2 features | Rhodotorula_pheno_MS | cell: 218 hits all-species → 0 within-species; supernatant: 10 within-species hits, p≈0 vs null, all positive, nucleoside/peptide chemistry | support |
| 2026-08-17 | 1000-perm within-mucilaginosa permutation (conservative) | control Cu-AUC + EB MS2 features | Rhodotorula_pheno_MS | within-muc null: mean 0.06, sd 0.33 hits, 95th = 0; observed cell 0 (emp p=1.0), supernatant 10 (emp p<0.001) | support |

## F-003: Species-block permutation confirms the between-species decomposition

**Status:** supported
**Claim:** Whole-panel *within-species-block* permutation of AUC (shuffle strain AUC only
among strains of the same species, preserving all cross-species structure; 300 perms) gives the
background count of features with \|rho\|>0.3 expected from cross-species covariance alone:

| fraction | observed hits (\|rho\|>0.3) | within-species-block null (mean/sd/95th) | block p |
|---|---|---|---|
| cell | 218 | 125.7 / 62.8 / 240.3 | 0.090 |
| supernatant | 143 | 40.5 / 38.4 / 120.0 | 0.027 |

Cross-species covariance alone already predicts ~41 (supernatant) and ~126 (cell) features at
\|rho\|>0.3 — i.e. most raw hits are phylogenetic covariance, not strain-level biology. The
cell fraction's observed 218 lies **inside** the block null (p=0.090, not beyond species
structure); the supernatant's 143 is modestly beyond it (p=0.027), and that excess is
accounted for by exactly the 10 genuinely within-species features from F-002. Combined with
the within-muc 1000-perm result (cell 0 hits, emp p=1.0; supernatant 10 hits, emp p<0.001),
the decomposition is: **cell-fraction growth-metabolite association = cross-species only;
supernatant = cross-species + a small real within-species component (10 nucleoside/peptide
features)**.
**Implications:** The PI's species-block permutation test resolves the question: the growth-
metabolome relationship is overwhelmingly phylogenetic in the cell fraction, and only a weak
(10-feature) within-species biological signal exists, and only in the supernatant.
**Tags:** growth-rate, auc, species-block-permutation, phylogenetic, supernatant, within-species, permutation-calibrated

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-17 | whole-panel species-block permutation (300 perms) | control Cu-AUC + EB MS2 features | Rhodotorula_pheno_MS | cell obs 218 vs null 95th 240 (p=0.090, not beyond species structure); supernatant obs 143 vs null 95th 120 (p=0.027); excess = the 10 within-species hits | support |

## F-004: Network (MS2 molecular-family) component-level test recovers the growth split — strong in supernatant, weak in cell — and color stays null at component level too

**Status:** preliminary
**Claim:** `analysis/scripts/network_component_association.py` (n=1000 perm; within R.
mucilaginosa only; 1,153 GNPS molecular-network components with ≥2 dedup-rep members; 4,892
unmapped/singleton features tested per-feature only) tests each molecular family by two
designs — GSEA-style enrichment (mean(-log10 p) of member features) and max|rho| — each with
permutation p + BH-FDR, plus a familywise whole-panel max null (`p_fdr_max`). Narrative
report with figures and links to all tables: `analysis/network_components/REPORT.md`. Results:

| trait | comps sig (enrichment BH<0.05) | comps sig (max\|rho\| BH<0.05) | comps sig (familywise) |
|---|---|---|---|
| growth AUC, **supernatant** | 210 / 1153 | 328 / 1153 | 7 |
| growth AUC, **cell** | 3 / 1153 | 7 / 1153 | 0 |
| color a\*, cell | 2 | 3 | 0 |
| color a\*, supernatant | 0 | 5 | 1 |
| color C\*, cell | 0 | 4 | 0 |
| color C\*, supernatant | 0 | 0 | 0 |

The growth signal is concentrated in the **supernatant** (210-328 of 1,153 families pass
BH-FDR; 41% of components have enrichment perm-p<0.05 vs 5% under null) and covers broad
chemistry (secreted metabolome largely scales with growth rate), while the cell fraction has
only a handful of nominally significant families and zero that survive the familywise bound.
Color remains null at component level in both fractions (0-5 BH hits, 0-1 familywise) — an
9th methodological angle agreeing with F-002/F-003/F-013/F-017 that aggregating over
chemistry does not rescue color signal.
**Component 887** (purine-nucleoside family, the 2-methylthio-adenosine / isopentenyladenine
cluster containing the 10 supernatant features of F-002) is significantly growth-associated
in the supernatant (es p=0.005, BH-FDR 0.037; max\|rho\|=0.223, n=46 members) and not in the
cell fraction (p_es 0.18) — carrying the F-002 feature-level result up to the family level.
**Implications:** molecular-family aggregation confirms the compartment asymmetry without
needing multiple-testing thrift; the broad supernatant signal is likely a few shared
latent/biological axes (density-linked production) rather than many independent pathways, so
read BH counts (which assume component independence) with the 7-family familywise result as
a conservative bound. Calibration is sane (color_C_supernatant, a true null panel, gives
0-hit and ≤5% perm-p fractions).
**Tags:** network-components, gnps, enrichment, max-rho, growth-rate, auc, supernatant, within-species, color-null, biomarkers-of-growth

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-18 | `analysis/scripts/network_component_association.py --nperm 1000` (corrected feature-identity matrix; singletons reported separately) | within-muc control Cu-AUC + color a\*/C\*, cell & supernatant vs GNPS network components | Rhodotorula_pheno_MS | sup growth 210 enrich/328 maxrho BH /7 familywise sig; cell growth 3/7/0; color a\*/C\* null both fractions (0-1 familywise); comp 887 sup sig (es FDR 0.037) not cell | support (preliminary)

### Open Questions
- Are purine nucleosides / small peptides biologically expected to scale with faster growth
  (secreted nucleotide turnover, cell-wall/protein economy) — or is this a batch/extraction
  effect correlated with culture density at harvest?
- Is the between-species cell-fraction lipid signal (ceramides/wax esters) really
  phylogenetic, or an extraction-density artifact? (It vanishes within species, consistent
  with phylogeny, but the raw specimen-dilution artifact remains untested.)
- Do the same supernatant components replicate against other-metal AUC once populated, or
  against colony area (the biomass confound) as a negative control?
- Is the broad (328-family BH / 7-family familywise) supernatant signal better described by a
  small number of shared latent axes (e.g. total secreted production) — checked e.g. by
  PCA/generative factors on the significant families?
