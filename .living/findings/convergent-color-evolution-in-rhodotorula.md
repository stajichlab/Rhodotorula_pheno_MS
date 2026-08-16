---
topic: convergent-color-evolution-in-rhodotorula
description: Whether darker/redder colony color has evolved independently more than once across the Rhodotorula species panel, and where on the phylogeny.
created: 2026-08-15
last_updated: 2026-08-15
status: active
---

# Convergent color evolution in Rhodotorula

## F-001: Formal Bayesian regime-shift detection (bayou) finds diffuse, not decisive, support for color-shift locations, but recovers the same top candidate as the earlier coarse heuristic
**Status:** preliminary
**Claim:** A reversible-jump OU regime-shift model (bayou 2.3.2, 20,000
MCMC generations, single chain) fit to species-mean a\* on the 17-species
tree found posterior support for ~1-2 regime shifts (k mean 1.8, 95% HPD
0-4) but **no single branch reached a decisive posterior probability**
(max ~0.21, below the conventional ~0.3 cutoff for a confident call). The
top-ranked candidate clade, *R. sphaerocarpa* + *R. taiwanensis* (pp
0.177, clade mean a\*=11.75 vs. nearest sister *R. sp. clade XI* a\*=7.58),
independently matches the #1/#2-ranked species by a\* from this project's
much earlier (Phase 1) simple ranking — an external consistency check
between two independent methods. *R. glutinis* (pp 0.176) also matches a
candidate flagged by the original coarse convergent-evolution heuristic
(`convergent_color_test.R`).
**Implications:** The diffuse posterior support (no branch clears a
confident threshold) is itself a finding, not a failed analysis — it's
the same "not enough independent phylogenetic data points" story that
recurs throughout this project (Phase 1 phylogenetic-signal test, Phase 2
species-level block permutation), now demonstrated in a formal Bayesian
framework rather than an ad hoc heuristic. With only 17 species tips,
this tree likely cannot support confident, publication-grade localization
of discrete color-gain events — a hard ceiling on this class of question
for this dataset, independent of which method is used. The two most
defensible candidate origins for downstream work (Idea 3's genome-side
follow-up, Idea 5's Step 3) are *R. sphaerocarpa*+*R. taiwanensis* and
*R. glutinis* — both corroborated by an independent method; the
*Pseudomicrostroma phylloplanum* and *R. mucilaginosa* candidates in the
same top-5 list are lower-confidence (likely tree-placement/branch-length
artifacts, see full writeup) and should not be treated as equally
supported.
**Tags:** convergent-evolution, bayou, regime-shift, phylogenetic-comparative-methods, power-ceiling, taiwanensis, sphaerocarpa, glutinis

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | bayou regime-shift MCMC (idea5_regime_shift_detection.R) on species_tree.nwk, a\*_mean | species_phenotype_table.csv (17 species) | Rhodotorula_pheno_MS | No branch reached decisive pp; top candidate (R. sphaerocarpa+R. taiwanensis) matches Phase 1's independent top-2 ranking | supports |

### Open Questions
- Would a multi-chain run with Gelman-Rubin convergence diagnostics and
  more generations sharpen the posterior, or is 17 tips a hard ceiling
  regardless of chain length/count?
- Does the same top-2 candidate set (R. sphaerocarpa+R. taiwanensis;
  R. glutinis) hold up when rerun on orange_score_mean instead of a\*
  (not yet done — exploratory tier)?
- Once Idea 3's candidate carotenoid-pathway-gene calls exist: do these
  two candidate clades independently show the same genomic signature
  relative to their paired sisters (Idea 5 Step 3)?
