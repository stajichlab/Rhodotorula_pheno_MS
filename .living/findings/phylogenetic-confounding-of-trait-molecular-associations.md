---
topic: phylogenetic-confounding-of-trait-molecular-associations
description: Cross-cutting pattern across independent phenotype/molecular-layer pairs in this project — whole-panel/naive trait-molecular associations recur, but do not survive phylogenetic correction or within-species restriction.
created: 2026-08-15
last_updated: 2026-08-15
status: active
---

# Phylogenetic confounding of trait <-> molecular-feature associations

## F-001: Whole-panel trait-molecular "hits" have not survived phylogeny-aware or within-species testing in any pairing tried so far
**Status:** preliminary
**Claim:** Two independent phenotype/molecular-layer pairings in this
project show the same structural pattern: (1) **colony color (a\*/C\*)
vs. MS2 metabolome** — no formal whole-panel naive test was run (analysis
went straight to phylogeny-aware testing), but the phylogeny-aware
whole-panel test, a within-species test in the largest/most color-diverse
species (*R. mucilaginosa*, n=206), two smaller within-species robustness
checks (*R. paludigena*, *R. toruloides*, n=10 each), and a sparse
multivariate model all found **zero** detectable association, despite
color varying across nearly its full observed range within a single
species alone (see
`phenotype-metabolome-association-statistical-power.md`, F-001 through
F-004). (2) **Copper-resistance growth rate (AUC) vs. proteome-wide amino
acid composition** (`analysis/copper/` — a separate, earlier analysis) —
a *naive* whole-panel Pearson correlation found 5 amino acids (S, L, Q, W,
T) significant (BH-FDR q<0.05, `analysis/copper/outputs/naive_vs_pgls_comparison.csv`),
and PGLS (phylogeny-corrected, whole-panel) still called several of these
"significant" — but a **within-*R. mucilaginosa*-only sensitivity check
found none of them held up** (`sensitivity_mucilaginosa_only.csv`: p =
0.20-0.68 for all 6 tested amino acids, none close to significant).
**Implications:** Both cases show the same shape: apparent trait-molecular
associations are easy to find at the whole-panel scale (species differ
in many correlated ways simply by virtue of shared ancestry — geography,
ecology, everything), and even whole-panel PGLS correction does not
always catch this (case 2 above: PGLS still called several amino acids
significant at the whole-panel level despite the signal not holding up
within one species). **The within-species restriction is currently the
single most diagnostic test available in this project for distinguishing
a real trait-molecular link from a phylogenetic/ecological confound** —
more diagnostic than whole-panel PGLS alone. Any future whole-panel
"hit" in this project (Phase 3 SIRIUS enrichment, Phase 5 genome linkage,
or a revisited copper-AUC analysis) should be treated as provisional
until checked against a within-species (or few-species) restriction,
given the 100% "does not replicate within-species" track record so far
across the two pairings tested.
**Tags:** phylogenetic-confounding, within-species-validation, methodology, cross-cutting, color, copper-resistance, growth-rate

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | Phase 2 color-metabolome (whole-panel + within-species + multivariate) | control_90_110 phenotype + EB MS2 features | Rhodotorula_pheno_MS | No association at any scale tested (4 methods) | supports |
| (earlier, exact date not in this session's context) | copper/scripts 02_pgls_analysis.R + 04_sensitivity_species_balance.R | Cu_AUC phenotype + BFD aa_freq | Rhodotorula_pheno_MS | Naive + whole-panel-PGLS "hits" (S,L,Q,W,T + C,E) did not replicate in the R. mucilaginosa-only sensitivity check | supports |

### Open Questions
- Is this pattern specific to how strongly phylogenetically structured
  *Rhodotorula* species/strains are in this dataset (heavy sampling
  imbalance, e.g. 216/303 strains are *R. mucilaginosa*), or would it
  recur in a differently-structured panel?
- **PI plans to revisit the copper-resistance AUC data/analysis** (noted
  2026-08-15) — the naive/PGLS/sensitivity results above should be treated
  as a snapshot from the current data version, not a final word; update
  this entry when that revisit happens.
- Would a within-species design be worth adopting as the *default* primary
  test in this project going forward, given its track record, rather than
  whole-panel PGLS as primary with within-species as a secondary check (the
  current convention per the Species-Level Collapse procedure)?
