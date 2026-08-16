# Idea 5, Steps 1-2: Formal regime-shift detection + contrast pairs (2026-08-15)

Scripts: `analysis/scripts/idea5_regime_shift_detection.R` (Step 1),
`analysis/scripts/idea5_contrast_pairs.R` (Step 2). Replaces the coarse
heuristic in `convergent_color_test.R` (above-mean orange_score AND
phylogenetically-distant-from-*R.-dairenensis*) with a formal Bayesian
reversible-jump OU model (`bayou` 2.3.2) fit to a\* (primary predictor)
on the 17-tip species tree.

## Step 1: bayou MCMC result

- 20,000 generations, 30% burnin, single chain (not yet run with multiple
  independent chains for a formal convergence check — see Caveats).
- Effective sample sizes: 51-206 across parameters (borderline for
  `alpha`, adequate for the rest) — fine for an exploratory pass, not
  yet publication-grade.
- **k (number of regime shifts): posterior mean 1.8, 95% HPD interval
  0-4.** The tree supports roughly 1-2 shifts but with substantial
  uncertainty — expected given only 17 tips.
- **No single branch reached a decisive posterior probability**: the
  top 5 branches all sit in the 0.10-0.21 range (see
  `regime_shift_amean_branch_posterior.csv`), well below what would
  normally be called a confident shift call (`shiftSummaries()` found
  nothing clearing pp>0.3). This is itself informative — it is the same
  "not enough independent data points" power-ceiling story that has
  recurred throughout this project (Phase 1's phylogenetic signal test,
  Phase 2's species-level block permutation), now showing up in a
  formal Bayesian framework instead of an ad hoc heuristic.

## Step 2: candidate contrast pairs

Top 5 branches by posterior probability, mapped to descendant clades and
paired against their nearest phylogenetic non-candidate sister
(`contrast_pairs.csv`):

| pp | Candidate clade | Clade mean a\* | Nearest non-candidate sister | Sister a\* | Patristic dist. |
|---|---|---|---|---|---|
| 0.177 | *R. sphaerocarpa* + *R. taiwanensis* | **11.75** | *R. sp. clade XI* | 7.58 | 0.567 |
| 0.176 | *R. glutinis* | 5.81 | *R. graminis* | 6.72 | 0.232 |
| 0.147 | *Pseudomicrostroma phylloplanum* (outgroup) | 10.90 | *R. toruloides* | 8.05 | 1.140 |
| 0.142 | *R. mucilaginosa* | 10.26 | *R. pacifica* | 10.56 | 0.043 |
| 0.104 | *R. glutinis* + *R. sp. clade I* | 7.91 | *R. graminis* | 6.72 | 0.261 |

**Consistency check with prior work**: the top candidate
(*R. sphaerocarpa* + *R. taiwanensis*) is exactly the pair that ranked
#1/#2 by a\* in Phase 1's species ranking (`convergent_color_candidates.csv`)
— an independent, formal method recovering the same signal the earlier
coarse heuristic and the raw ranking both pointed to. That's a real
consistency win for the analysis pipeline, not a coincidence to explain
away.

**Caveats on individual rows**:
- The *Pseudomicrostroma phylloplanum* row is a single-strain outgroup
  genus with substantial phylogenetic distance to its paired sister
  (1.14, the largest in the table) — plausibly a tree-placement/outgroup
  artifact rather than a genuine within-*Rhodotorula* carotenoid-pathway
  convergence event. Treat with more skepticism than the *Rhodotorula*-
  only rows.
- The *R. mucilaginosa* row has a suspiciously tiny patristic distance to
  its "nearest non-candidate sister" (0.043) — likely reflects
  *R. mucilaginosa*'s own branch-length placement (206 strains collapsed
  to one representative tip) rather than a meaningful evolutionary
  distance; also *R. mucilaginosa*'s a\* (10.26) barely differs from its
  paired sister's (10.56), so this row likely doesn't represent a real
  color-gain event regardless of its posterior probability — the pp
  ranking is picking something else up (possibly tree-topology/branch-
  length structure near the root of the *Rhodotorula* clade) that
  deserves a closer look before being treated as a "candidate origin."

## What this means for the rest of Idea 5 / Idea 3

Given the diffuse posterior support, Step 3 ("do independent origins share
molecular/genomic correlates") should proceed but with explicit weighting
toward the two cleanest rows — *R. sphaerocarpa*+*R. taiwanensis* (pp
0.177, biologically consistent with prior ranking) and *R. glutinis* (pp
0.176, consistent with the original heuristic's flagged candidates) — and
treat the *Pseudomicrostroma* and *R. mucilaginosa* rows as lower-
confidence / needing manual review before inclusion in any contrast
analysis.

This is exactly the kind of candidate list Idea 3's genome-side work
should eventually test against once the BFD rebuild + candidate-gene
ortholog calls are ready: do *R. sphaerocarpa*/*R. taiwanensis* and
*R. glutinis* independently show the same *crtS* (or other candidate
pathway gene) variant/copy-number pattern relative to their paired
sisters?

## Caveats on the method itself

- Single MCMC chain, 20,000 generations — adequate to see the qualitative
  story (diffuse support, no decisive shift) but not a publication-grade
  run. A rigorous version would run ≥2 independent chains, check Gelman-
  Rubin diagnostics (`gelman.R()`, available in bayou), and likely
  increase generations until all parameters' ESS clear ~200-500.
- Only a\* was tested (primary predictor per the 2026-08-15 grilling
  session decision). orange_score_mean can be rerun the same way
  (`--trait-col orange_score_mean`) as an exploratory-tier follow-up —
  not yet done.
- 17 tips is small for bayou by the standard of most published
  applications (typically 50-300+ tips) — the diffuse result should be
  read as "this tree probably doesn't have enough independent lineages to
  formally localize shifts with confidence," which is itself a valid and
  useful finding, not a failed analysis.

## Reproduce

```
Rscript analysis/scripts/idea5_regime_shift_detection.R
Rscript analysis/scripts/idea5_contrast_pairs.R
```
