# Integrated Analysis — Executive Summary

> Navigation + verdict document for `analysis/integrated_analysis/`. Full
> plan and rationale: `analysis/INTEGRATED_ANALYSIS_STRATEGY.md`. Full
> evidence ledgers: `.living/findings/*.md`
> (`.living/findings/FINDINGS_REGISTRY.md` is the master index). This file
> answers "what's been tried and what's the verdict" — read the linked
> `.md`/`.csv` files for methods and full numbers. No plots exist in these
> phase directories (all output is tabular/CSV) except where noted.

| Phase | Question | Status | Verdict |
|---|---|---|---|
| [1 — Phenotype](#phase-1--phenotype-characterization-done) | Characterize color phenotype, phylogenetic structure | Done | Reframed the whole project |
| [2 — Metabolome](#phase-2--color--metabolome-association-done-null) | Does color correlate with any of 10,949 MS2 compound groups? | Done | **Null**, 5 independent methods agree |
| [3 — Idea 1 targeted re-mining](#phase-3--idea-1-targeted-carotenoidsterol-mass-re-mining-done--live-lead) | Can targeted mass search find pigment features the untargeted pipeline missed? | Done | Found candidates; **no color link**, but a real copper-AUC lead (unvalidated) |
| [5 — Idea 5 regime shift](#phase-5--idea-5-bayesian-regime-shift-detection-done) | Where on the tree did color shift convergently? | Done | Diffuse posterior; top candidate corroborates Phase 1 ranking |
| [5 — Idea 3 pigment genes](#phase-5--idea-3-pigment-gene-genome-screen-done-for-screening-stage) | Are pigment-pathway genes present/variable across the genome panel? | Done (screening stage) | Melanin = laccase route (convergent 3-method finding); first genome↔color test not yet run |

---

## Phase 1 — Phenotype characterization (done)
**Path**: `phase1_phenotype/`
**Scripts**: `analysis/scripts/build_strain_phenotype_table.py`, `build_species_level_tables.py`, `prune_species_tree.R`, `phylogenetic_signal.R`, `convergent_color_test.R`

- `strain_phenotype_table.csv` / `species_phenotype_table.csv` — CIELAB L\*/a\*/b\*/C\*/h°/orange_score at strain and species level (canonical `control_90_110` source).
- `phenotype_phylogenetic_signal.csv` — Blomberg's K / Pagel's λ per trait.
- `species_tree.nwk` — 17-tip species-level pruned tree (278-taxa PHYling backbone).
- `convergent_color_candidates.csv` — heuristic convergence screen (superseded by Idea 5's formal Bayesian version, see Phase 5 below).

**Key finding, project-reframing**: *R. dairenensis* — the species the project was originally built around — ranks only **4th of 17-18 species** on the orange_score composite, reproduced across multiple phenotype sources (`analysis/examine_phenotype_calling/RECOMMENDATION.md`). This motivated the 2026-08-15 pivot from an *R.-dairenensis*-specific investigation to a whole-panel color→compound→genome framing (see strategy doc's "Reframing" section).

---

## Phase 2 — Color ↔ metabolome association (done, **null**)
**Path**: `phase2_metabolome_phenotype/`
**Scripts**: `analysis/scripts/phase2_color_metabolome_association.py`, `phase2_within_species_association.py`, `phase2_multivariate_association.py`, `phase2_anova_pattern_association.py`

| Method | File | n | Result |
|---|---|---|---|
| Whole-panel univariate + phylo block-permutation | `PHASE2_SUMMARY.md`, `color_metabolome_association_{a,C}.csv` | ~275 strains, 17-18 species blocks | **0/10,164-10,416** hits at BH-FDR<0.05, either color axis, either fraction |
| Within-species univariate, *R. mucilaginosa* | `WITHIN_SPECIES_MUCILAGINOSA.md`, `within_species_Rhodotorula_mucilaginosa_association_a.csv` | n=206, real strain-tree blocking | **Null**, despite near-full panel-wide color range within this one species |
| Within-species univariate, 2 more species | `ROBUSTNESS_AND_MULTIVARIATE.md`, `within_species_Rhodotorula_{paludigena,toruloides}_association_a.csv` | n=10 each | Null |
| Within-species univariate, remaining 5 species with ≥5 strains | `WITHIN_SPECIES_SMALL_SPECIES_SWEEP.md`, `within_species_Rhodotorula_{dairenensis,diobovata,taiwanensis,sp_clade_I,sphaerocarpa}_association_a.csv` | n=5-8 each | Null, but negative control uninformative at this n (flagged explicitly) |
| Sparse multivariate (Lasso, group-CV) | `ROBUSTNESS_AND_MULTIVARIATE.md`, `multivariate_Rhodotorula_mucilaginosa_a.csv` | *R. mucilaginosa* | Null — cross-val R² negative throughout |
| Pattern-group ANOVA (k-means + Kruskal-Wallis) | `anova_Rhodotorula_mucilaginosa_color.csv` | *R. mucilaginosa* | Null |

**Every species in the panel with ≥5 strains and MS data (8 of 17-18) has now been tested within-species — all null.**

**Verdict**: 5 independent statistical methods agree — **no detectable color↔metabolome association** in this dataset at current resolution. Hard-gated negative controls (`*_area_decoy.csv` files, colony area as a phylogenetically-structured but color-unrelated decoy) confirm each pipeline can detect a real effect when present (area decoy shows 1,354-2,025 hits depending on method/fraction), so this is not a power/calibration artifact of the pipeline itself.
**Real side-finding**: colony area broadly confounds cell-fraction (not supernatant) metabolite abundances — `.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md`.
**Full evidence ledger**: `.living/findings/phenotype-metabolome-association-statistical-power.md` (F-002 through F-004); cross-cutting pattern (color AND copper-AUC naive hits both fail within-species restriction): `.living/findings/phylogenetic-confounding-of-trait-molecular-associations.md` (F-005).

This 5-method null result triggered the `analysis/ideas/2026-08-15-color-metabolome-genome-null-brainstorm/` ideation session (see `analysis/ideas/EXECUTIVE_SUMMARY.md`) that produced Ideas 1/3/5/6 below.

---

## Phase 3 — Idea 1: targeted carotenoid/sterol mass re-mining (done, **live lead**)
**Path**: `phase3_metabolome_phenotype_idea1/`
**Script**: `analysis/scripts/idea1_targeted_mass_remining.py`, `idea1_auc_quickcheck.py`

- `RESULTS.md` — original 9-compound carotenogenesis pathway search (11 mass matches). Standout: a torularhodin-mass candidate (row 21315) that SIRIUS actively misannotated as "Polyamines" — confirms a concrete misannotation risk.
- `MS2_FRAGMENTATION_CHECK.md` — manual spectral inspection demoted one candidate (phytofluene, row 1735 — precursor-stability artifact, no diagnostic fragment), left row 21315 unconfirmed-but-plausible.
- `EXPANDED_SEARCH_RESULTS.md` — expanded to 19 compounds/4 categories (apocarotenoids, oxygenated carotenoids, sterol-pathway markers): 31 matches. Standout: an **ergosterol candidate (row 846)**, independently ISF-confirmed by EverythingBagel's own pipeline.
- `STEROL_CLUSTER_AUC_CHECK.md`, `sterol_cluster_auc_quickcheck.csv` — extended to a 4-feature ergostane cluster (SIRIUS-called Peroxyergosterol, Ergost-3,5,7,9(11),22-pentaen, 7-Hydroxyergosterol + row-846 anchor). **All 4 show the same positive-direction naive correlation with copper-resistance AUC in the cell fraction** (rho 0.17-0.28, p<0.006).

**Verdict on color**: no candidate from either pass correlates with a\* — this re-mining did not rescue a hidden color signal, consistent with Phase 2.
**Verdict on copper-AUC**: a real, multi-feature-corroborated **lead, not a validated result** — no phylogenetic block-permutation or negative control has been run on it yet, and this project has a track record of naive whole-panel hits (including this exact AUC phenotype) failing that check once tested.
**Full evidence ledger**: `.living/findings/carotenoid-pathway-detectability-in-untargeted-lcms.md` (F-001 through F-004; registry F-007/F-008/F-009).
**Next step (not started)**: phylogenetically-corrected test of the ergosterol/sterol-cluster↔copper-AUC signal (predictor-swapped variant of `phase2_within_species_association.py`, `mean_auc_rate` instead of color).

---

## Phase 5 — Idea 5: Bayesian regime-shift detection (done)
**Path**: `phase5_genome_linkage/idea5_regime_shift/`
**Scripts**: `analysis/scripts/idea5_regime_shift_detection.R`, `idea5_contrast_pairs.R`

- `RESULTS.md`, `regime_shift_amean_summary.txt` — formal `bayou` reversible-jump OU MCMC on a\* (species tree), replacing the earlier coarse heuristic (`convergent_color_test.R`, Phase 1). k (number of shifts): posterior mean 1.8, 95% HPD 0-4. **No branch reaches a confident posterior probability** (top 5 branches all 0.10-0.21) — a power-ceiling result, same story as Phase 1/2, now in a formal Bayesian framework.
- `contrast_pairs.csv` — top candidate branches paired with nearest non-candidate sister clade by patristic distance. Top candidates: *R. sphaerocarpa*+*R. taiwanensis* (pp=0.177) and *R. glutinis* (pp=0.176) — **these independently match Phase 1's original heuristic #1/#2 ranking**, despite the low absolute posterior support.

**Verdict**: diffuse/inconclusive on its own, but the cross-method agreement with Phase 1 is itself informative — these two clades are the most defensible convergent-color-evolution candidates in the project so far, worth using as the contrast pairs if/when Idea 3's genome data supports a convergence test.
**Full evidence ledger**: `.living/findings/convergent-color-evolution-in-rhodotorula.md` (F-006).

---

## Phase 5 — Idea 3: pigment gene genome screen (done for screening stage)
**Path**: `phase5_genome_linkage/idea3_pigment_hmm_search/`
**Scripts**: `scripts/swissprot_pigment_crossref.py`, `scripts/pfam_pigment_screen.py`, `scripts/parse_hmm_hits.py`

Three independent screens of pigment-pathway gene presence/copy-number across all 278 BFD genomes (2,188,032 proteins): a PI-curated 28-profile custom HMM panel, a coarse Pfam-domain screen, and a SwissProt best-hit keyword cross-reference. Outputs: `outputs/pigment_hmm_strain_summary.csv`, `outputs/pfam_pigment_strain_summary.csv`, `outputs/swissprot_pigment_strain_summary.csv` (all 278-strain × gene-family matrices), plus per-protein hit tables.

**Verdict, convergent finding**: all 3 methods independently agree this genus melanizes via the **laccase route, not tyrosinase** (laccase present ~278/278 strains in all 3 screens; tyrosinase near-absent in all 3).
**Verdict, caution**: 10/28 custom HMM profiles found zero hits genome-wide (likely lineage-specific paralog discrimination — needs a PI sanity check on 4 gene names: `ayg1`/`crtB`/`hppd`/`scd`); 5 HMM-panel + 2 Pfam families show implausibly high copy numbers (superfamily-level domain matches, not confirmed single-gene orthologs — not yet safe to use as a predictor).
**Full evidence ledger**: `.living/findings/pigment-gene-genomic-screen-rhodotorula.md` (F-001 through F-003; registry F-010/F-011).
**Next step (not started, this is the project's first real genome↔color test)**: test the sparse-presence families (`crtR` in 52/278 strains, `crtQ` in 48/278, `hgd` in 28/278) as presence/absence predictors against species-level color.

---

## Phase 4 (classifier) and Phase 6 (deferred) — not started
Per `analysis/INTEGRATED_ANALYSIS_STRATEGY.md`: Phase 4 (whole-panel classifier linking color→compound class→gene) has no work yet — it depends on Phase 3 producing an enrichment-worthy feature list, which Phase 2's null result removed. Phase 6 is explicitly deferred in the strategy doc.
