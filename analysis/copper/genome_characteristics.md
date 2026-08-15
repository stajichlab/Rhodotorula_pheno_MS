# Genome amino-acid composition vs. copper resistance: a first-principles re-test

## Goal

Previous work (`../../Rhodotorula_MS2_pheno_explore/`) reported relationships between
strain amino-acid composition and copper tolerance. Before extending or trusting
that result, this analysis re-derives it from scratch with the simplest defensible
model at each step, and explicitly tests whether the trend survives correction for
shared ancestry (phylogenetic non-independence) using the BFD `phyling` protein
tree. If a naive (species-as-independent-samples) correlation disappears under
phylogenetic correction, it was very likely a lineage/taxonomy artifact, not a
trait-level association with copper resistance.

## Data sources

| Data | Path | Role |
|---|---|---|
| Copper resistance | `data/metadata/EXFAB_UCR-005/Cu_AUC.20260811.fixed.csv.gz` | `mean_auc_rate` = growth-curve AUC in copper-stress medium, per strain (`SAMPLE_NAME`), 275 rows |
| Genome-wide AA composition | `BFD/results/genome_stats/aa_freq/*/<strain>.aa_freq.csv.gz` | Proteome-wide frequency of each of the 20 standard amino acids, one file per strain (150 strains with genomes/annotations in BFD) |
| Phylogeny | `BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/protein-Rhodotorula-taxa_276.fungi_odb10.fasttree.support.treefile` | FastTree ML tree from PHYling BUSCO (fungi_odb10) protein markers, 276 tips, used for phylogenetic correction |

## Strain identifier reconciliation

The three sources name strains differently (tree tip `Rhodotorula_mucilaginosa_DBVPG_3776.proteins`,
aa_freq file `DBVPG3776.aa_freq.csv.gz`, metadata `DBVPG_3776`). `scripts/common.py`
normalizes all three (strip `_`/`-`/whitespace, uppercase) and matches metadata/aa_freq
keys as substrings of the normalized tree-tip label. This is a pragmatic, verifiable
join — diagnostics are written to `outputs/join_diagnostics.txt` on every run so
silent mismatches are visible rather than assumed away.

**Coverage:** of 275 Cu_AUC strains, **134 (49%)** have both a phyling tree placement
and an aa_freq file and form the analysis set (`outputs/copper_aa_master_table.csv`).
13 strains have no tree placement (not sequenced/assembled or excluded from PHYling);
128 have a tree placement but no aa_freq file (BFD annotation lags the phylogeny —
see `BFD/TODO.md`). **The 134-strain set is heavily species-skewed: 99/134 (74%) are
_R. mucilaginosa_**, with the next largest group (_R. paludigena_) at only 9 strains.
This imbalance is carried through to the results below and is the single biggest
caveat on the conclusions.

## Assumptions

1. **Proteome-wide AA frequency is a meaningful genome-level trait.** We use the
   whole-proteome average composition per strain (not per-gene, not restricted to
   any gene family/pathway) as the simplest possible genomic phenotype. This is
   deliberately blunt — it cannot distinguish "more Cys because of a few
   metallothionein-like proteins" from "globally Cys-richer proteome" — and is meant
   as a null-model-style screen, not a mechanistic claim.
2. **`mean_auc_rate` is treated as the copper-resistance phenotype**, used as-is with
   no additional normalization (assumed already comparable across strains/plates).
3. **Each strain is an independent statistical unit** in the naive model, and a tip
   in a phylogenetic tree with among-strain covariance in the PGLS model. Multiple
   strains of the same species are kept separate (not averaged to one point per
   species) so within-species variation is available to the model — but see the
   sensitivity check below for why this matters.
4. **Phylogenetic correction uses Pagel's λ (ape `corPagel`, ML-estimated)**, not a
   fixed Brownian-motion (λ=1) or star-phylogeny (λ=0) model. λ is fit per-AA rather
   than assumed, so the degree of phylogenetic signal in the residuals is itself an
   output, not an input assumption.
5. **Ten of the tree's edges among the 134 retained tips have branch length 0**
   (near-identical/unresolved strain pairs in the BUSCO-marker tree). These are
   nudged to length 1e-6 before fitting so the phylogenetic covariance matrix is
   non-singular — a standard workaround, but it means near-zero-branch strain pairs
   contribute close to zero independent phylogenetic information to λ estimation.
6. **Primary vs. exploratory hypotheses are pre-declared and reported separately.**
   Cys, His, Asp, Glu, Met are flagged `primary_hypothesis=TRUE` a priori (thiol/
   imidazole/carboxylate/thioether side chains with plausible direct roles in metal
   binding or redox buffering); all 20 AAs are also scanned exploratorily with
   Benjamini-Hochberg FDR correction (q < 0.05) across the 20 tests.

## Statistical strategy

1. **`00_build_master_table.py`** — join Cu_AUC to aa_freq via the tree, write
   `outputs/copper_aa_master_table.csv` + `outputs/join_diagnostics.txt`.
2. **`01_naive_correlation.py`** — per-AA Pearson & Spearman correlation of AA
   frequency vs. `mean_auc_rate`, strains treated as independent, BH-FDR across 20
   AAs. This reproduces the kind of test the prior analysis likely started from.
3. **`02_pgls_analysis.R`** — prune the tree to the 134 matched tips, refit each
   AA's relationship as `auc ~ AA_freq` by GLS with a `corPagel` correlation
   structure (λ estimated by ML), BH-FDR across 20 AAs.
4. **`03_compare_naive_vs_pgls.py`** — merge both result tables and classify each
   AA as `survives_phylo_correction`, `naive_only_lineage_artifact`, `pgls_only`
   (masked by phylogenetic confounding in the naive test), or `not_significant`.
5. **`04_sensitivity_species_balance.R`** — given the 74% _R. mucilaginosa_
   imbalance, refit the top hits *within* the _R. mucilaginosa_ subset alone (its
   own pruned subtree) to check whether the signal holds within one species or is
   solely an inter-species effect riding on the imbalance.

Run all steps: `bash scripts/run.sh`.

## Results (n = 134 strains, 2026-08-15 run)

**Naive test:** 12/20 AAs are significant at BH q < 0.05 (strong signal), led by Ser,
Leu, Gln, Trp, Thr (all q < 1e-5). This alone is not informative — with 74% of
strains from one species, a naive correlation is exactly what you'd expect from
species-level clustering in both AA usage and copper tolerance, independent of any
real trait relationship.

**PGLS (full 134-strain tree):** 7/20 AAs remain BH-significant. Comparing to the
naive test (`outputs/naive_vs_pgls_comparison.csv`):

| Verdict | AAs |
|---|---|
| `survives_phylo_correction` (naive-sig AND pgls-sig) | S, L, Q, W, T |
| `pgls_only` (masked in naive test, revealed by PGLS) | **C** (primary hypothesis), E (primary hypothesis) |
| `naive_only_lineage_artifact` (naive-sig, not pgls-sig) | A, F, H, N, G, V, I, P, K, M, Y |
| `not_significant` in either | D, R |

Two things stand out:
- **Cys is the strongest primary-hypothesis hit** (PGLS q = 0.0028) and only
  emerges *after* phylogenetic correction — its naive association was suppressed
  (q = 0.059), i.e. phylogenetic structure was masking, not manufacturing, this one.
  That is the pattern expected if Cys usage has a real, lineage-independent
  relationship with copper handling (thiol-based metal chelation is directly
  plausible for Cu).
- **11 of 20 naive hits (55%) evaporate** under phylogenetic correction — most of
  the naive "signal" was indeed a taxonomy/lineage artifact, exactly the failure
  mode this re-test was designed to catch.

Estimated λ values cluster around 0.25-0.5 for all AAs (moderate, not full,
phylogenetic signal in the AA~AUC residuals) — see `outputs/pgls_correlation_results.csv`.

**Sensitivity check (within-species):** Refitting PGLS for the five non-primary
top hits (S, L, Q, W, T) *and* Cys within the 99-strain _R. mucilaginosa_-only
subtree, **none remain significant** (all p > 0.2), and several flip sign (L, Q, W, C
all change direction relative to the full dataset; see
`outputs/sensitivity_mucilaginosa_only.csv`). This is the key caveat: the
"survives_phylo_correction" result in the full 134-strain tree is not corroborated
by within-species variation in the dominant species — it is consistent with an
across-species (inter-clade) effect that Pagel's λ, estimated on a tree where most
mass sits in one clade, does not fully absorb.

## Bottom line

The naive AA-vs-copper-AUC correlations reported previously are **not simply
reproduced as-is**: a majority collapse once phylogeny is modeled, confirming they
were largely lineage artifacts. A handful of associations (led by Cys, plausible
on mechanistic grounds) survive PGLS on the full tree, but **do not replicate within
the one species that supplies 74% of the data**, so the current data cannot
distinguish "real cross-lineage trait covariation" from "an artifact of uneven
species sampling that Pagel's λ under-corrects for." **Treat this as
hypothesis-generating (Cys/thiol content as a candidate copper-tolerance
correlate worth targeted follow-up), not as a confirmed trend.**

## Limitations / next steps

- Proteome-wide average composition is a blunt instrument; a natural follow-up is
  restricting to specific functional categories already in BFD (`function/pfam_hmmscan`,
  `function/cazy`, metal-transport/metallothionein-like Pfam domains) rather than
  whole-proteome averages.
- Species imbalance (74% _R. mucilaginosa_) limits what PGLS can actually
  distinguish; balancing the sampled tree (e.g. one strain per species, or explicit
  species-level random effects / phylogenetic mixed model) would be a stronger next
  step before treating any AA as a real correlate.
- 128 strains with Cu_AUC and a tree placement still lack an aa_freq file — closing
  that BFD annotation gap (per `BFD/TODO.md`) would roughly double the usable n and
  directly improve species balance if it disproportionately adds non-mucilaginosa
  strains.
- This analysis used AUC as reported; it does not re-derive or QC the growth-curve
  fitting itself.
- See `../YPD/color_shape_growth/YPD_COLOR_SHAPE_GROWTH.md` for the parallel,
  non-stress-condition check (growth/color under plain YPD) using the same
  join + naive/PGLS approach, run as an independent baseline.
