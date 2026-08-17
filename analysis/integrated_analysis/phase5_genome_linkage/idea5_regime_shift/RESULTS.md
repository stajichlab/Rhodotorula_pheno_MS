# Idea 5, Steps 1-2: Formal regime-shift detection + contrast pairs (2026-08-15)

Scripts: `analysis/scripts/idea5_regime_shift_detection.R` (Step 1),
`analysis/scripts/idea5_contrast_pairs.R` (Step 2). Replaces the coarse
heuristic in `convergent_color_test.R` (above-mean orange_score AND
phylogenetically-distant-from-*R.-dairenensis*) with a formal Bayesian
reversible-jump OU model (`bayou` 2.3.2) fit to a\* (primary predictor)
on the 16-tip species tree.

## Step 1: bayou MCMC result

- 40,000 generations, 30% burnin, single chain (not yet run with multiple
  independent chains for a formal convergence check — see Caveats). Re-run
  on the post-ANI 16-tip species tree (R. pacifica removed; the two
  reclassified strains merged into R. mucilaginosa).
- Effective sample sizes: 141-388 across parameters (adequate for an
  exploratory pass, `alpha` borderline) — not yet publication-grade.
- **k (number of regime shifts): posterior mean 1.39, 95% HPD interval
  0-3.** The tree supports roughly 1-2 shifts but with substantial
  uncertainty — expected given only 16 tips.
- **A single branch now carries most of the signal**: the terminal branch
  to *R. taiwanensis* has posterior shift probability 0.256 (see
  `regime_shift_amean_branch_posterior.csv`), with its sister *R.
  sphaerocarpa* branch next at 0.119. Still below the pp>0.3 bar
  `shiftSummaries()` uses for a confident call, but a much more coherent
  distribution than the pre-reassignment run — the high-a\* signal is
  concentrated on the *R. sphaerocarpa*/*R. taiwanensis* subclade rather
  than spread diffusely across near-root branches.

## Step 2: candidate contrast pairs

Top 5 branches by posterior probability, mapped to descendant clades and
paired against their nearest phylogenetic non-candidate sister
(`contrast_pairs.csv`):

| pp | Candidate clade | Clade mean a\* | Nearest non-candidate sister | Sister a\* | Patristic dist. |
|---|---|---|---|---|---|
| 0.256 | *R. taiwanensis* | **12.73** | *R. sp. clade XI* | 7.58 | 0.562 |
| 0.119 | *R. sphaerocarpa* | 10.76 | *R. sp. clade XI* | 7.58 | 0.571 |
| 0.080 | *R. diobovata* | 8.15 | *R. graminis* | 6.72 | 0.471 |
| 0.073 | *R. glutinis* | 5.81 | *R. graminis* | 6.72 | 0.232 |
| 0.072 | *R. glutinis* + *R. sp. clade I* | 7.91 | *R. graminis* | 6.72 | 0.261 |

**Consistency check with prior work**: the top candidate is now a single
species, *R. taiwanensis* (a\* 12.73 — the highest in the panel), with
its sister *R. sphaerocarpa* second — i.e. both members of the
*R. sphaerocarpa*/*R. taiwanensis* subclade that ranked #1/#2 by a\* in
Phase 1's species ranking (`convergent_color_candidates.csv`). A formal
Bayesian method recovering the same subclade the earlier coarse heuristic
and the raw ranking both pointed to is a real consistency win for the
analysis pipeline.

**Caveats on individual rows**:
- The *R. taiwanensis* and *R. sphaerocarpa* rows are a shared subclade
  (both descend from a common ancestor) — not independent origins, so as
  contrast units they count as a single candidate event (the high-a\*
  signal on the *R. sphaerocarpa*/*R. taiwanensis* subclade as a whole).
  The post-ANI removal of *R. pacifica* (which previously sat between
  *R. mucilaginosa* and the sphaerocarpa/taiwanensis clade as part of the
  diffuse signal) sharpened this into the clean two-species cluster.
- The *R. diobovata* row (pp 0.080) has a modest a\* contrast (8.15 vs
  6.72) — lower priority than the top pair but a reasonable secondary
  candidate to carry into Step 3.
- *R. glutinis* (a\* 5.81, low) appears as a candidate shift clade at
  pp 0.073 — this reflects the model's freedom to place a downward
  shift on a low-trait branch, not a color-gain event, so it is of lower
  biological interest for the "convergent color gain" question.

## What this means for the rest of Idea 5 / Idea 3

Posterior support for any single branch is still below a decisive
threshold, but the signal is now coherently concentrated on the
*R. sphaerocarpa*/*R. taiwanensis* high-a\* subclade. Step 3 ("do
independent origins share molecular/genomic correlates") should proceed
with explicit weighting toward the *R. taiwanensis*/*R. sphaerocarpa*
row pair (pp 0.256/0.119, the same clade ranked #1/#2 by a\* in Phase 1)
as the primary candidate, with *R. diobovata* as a secondary. The
post-reassignment tree no longer confounds the signal with a
*R. pacifica*-adjacent branch.

This is exactly the kind of candidate list Idea 3's genome-side work
should eventually test against once the BFD rebuild + candidate-gene
ortholog calls are ready: do *R. sphaerocarpa*/*R. taiwanensis* and
*R. glutinis* independently show the same *crtS* (or other candidate
pathway gene) variant/copy-number pattern relative to their paired
sisters?

## Caveats on the method itself

- Single MCMC chain, 40,000 generations — adequate to see the qualitative
  story (signal concentrated on the sphaerocarpa/taiwanensis subclade but
  still below a decisive pp threshold) but not a publication-grade run. A
  rigorous version would run ≥2 independent chains, check Gelman-Rubin
  diagnostics (`gelman.R()`, available in bayou), and likely increase
  generations until all parameters' ESS clear ~200-500.
- Only a\* was tested (primary predictor per the 2026-08-15 grilling
  session decision). orange_score_mean can be rerun the same way
  (`--trait-col orange_score_mean`) as an exploratory-tier follow-up —
  not yet done.
- 16 tips is small for bayou by the standard of most published
  applications (typically 50-300+ tips) — posterior support for any single
  branch stays below a decisive threshold, but the concentration of the
  signal on the sphaerocarpa/taiwanensis subclade is itself a useful and
  consistent finding, not a failed analysis.

## Reproduce

```
Rscript analysis/scripts/idea5_regime_shift_detection.R --ngen 40000
Rscript analysis/scripts/idea5_contrast_pairs.R
```
