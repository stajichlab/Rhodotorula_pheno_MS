# Idea Session: Restructuring the data to find color-metabolome/genome patterns after 5 null results

**Date**: 2026-08-15
**Trigger**: Five independent statistical methods (whole-panel + within-species correlation, sparse Lasso, pattern-group ANOVA) all found zero association between CIELAB color and MS2 metabolome abundance, with negative controls confirming the pipelines aren't simply underpowered. Genome-side association hasn't been attempted at all yet. Asked for a 7-persona brainstorm (chemistry, metabolomics/MS, fungal genomics, quantitative genetics, evolutionary biology, causal inference, ecology) on how to restructure the data or approach.

**Personas used**: 7 (3 catalog: Quantitative Geneticist, Evolutionary Biologist, Causal Inference Researcher, Ecologist; 3 custom: Natural Products/Analytical Chemist, Metabolomics/MS Specialist, Fungal Genomics Bioinformatician). 2 ideas each = 14 ideas total.

## Index

| # | Persona | Idea | Effort | Key risk (one line) |
|---|---|---|---|---|
| 1a | Chemist | [Targeted carotenoid extraction + APCI/APPI-MS](01_natural-products-chemist.md) | Medium | Sample viability / instrument access; small-n |
| 1b | Chemist | [Predicted-mass targeted re-mining of existing raw MS data](01_natural-products-chemist.md) | **Low** | None new — reuses existing files |
| 2a | Metabolomics/MS | [Molecular-network-aware compound family scores](02_metabolomics-ms-specialist.md) | Medium | Clustering resolution is a hidden researcher DOF |
| 2b | Metabolomics/MS | [Dark-matter MS2 spectral-similarity kernel regression](02_metabolomics-ms-specialist.md) | High | Research-grade method, needs custom implementation |
| 3a | Bioinformatician | [Candidate carotenoid pathway genotyping (BFD panel)](03_fungal-genomics-bioinformatician.md) | **Low** | Pfam domains too coarse; needs curation |
| 3b | Bioinformatician | [Genome-wide orthogroup PAV/CNV scan](03_fungal-genomics-bioinformatician.md) | High | Large multiple-testing burden, OrthoFinder compute cost |
| 4a | Quant. Geneticist | [Variance-component/heritability decomposition of color](04_quantitative-geneticist.md) | **Low** | Needs replicate-level (not just final) phenotype data |
| 4b | Quant. Geneticist | [Metabolome polygenic score vs. color](04_quantitative-geneticist.md) | Medium | Still needs a defensible aggregation weighting scheme |
| 5a | Evolutionary Biologist | [Convergence-restricted association test (independent color-gain events)](05_evolutionary-biologist.md) | High | N=independent origins (3-4), not N=strains — inherently low power |
| 5b | Evolutionary Biologist | [Environment/ecology as covariate](05_evolutionary-biologist.md) | Low–Medium | — |
| 6a | Causal Inference | [Determine area's causal role (confounder/mediator/collider) before trusting any null](06_causal-inference-researcher.md) | **Low** | If area is a collider, "just adjust for it" is actively wrong |
| 6b | Causal Inference | [Stop treating species as pure nuisance — genotype-to-color path-specific reanalysis](06_causal-inference-researcher.md) | Medium | BFD genotype-color test is a heavier lift than the rest |
| 7a | Ecologist | [Niche-stratified photoprotection — does color track environment?](07_ecologist.md) | **Low** | Origin/Environment metadata may be sparse/inconsistent |
| 7b | Ecologist | [Generalist vs. specialist niche breadth as predictor](07_ecologist.md) | High | Needs a defensible niche-breadth metric |

## Cross-cutting theme

Three personas independently converged on the same structural point from different
angles: **the project has been testing "does color correlate with individual
metabolite features" when it may need to test "can the data even see the
relevant molecules, and is color even causally upstream of what we're
measuring" first.** The chemist (1b), metabolomics specialist (2a), and
causal-inference researcher (6a) all flag a *prerequisite* check that's
cheap and hasn't been done, rather than a new expensive analysis.

## Panel's prioritized recommendations (top 7)

Ranked by (a) low effort / no new data required, and (b) resolving a
precondition that the existing 5 null results depend on:

1. **6a — Determine colony area's causal role before trusting any null.**
   Every one of the 5 existing tests used area purely as a decoy/negative
   control. Nobody has formally tested whether area is a *confounder*,
   *mediator*, or *collider* relative to color and metabolite abundance —
   if it's a collider (e.g., the image-segmentation pipeline's color
   estimate depends on colony size), the whole negative-control framing
   could be subtly miscalibrated. Cheapest, most foundational check on this
   list; should arguably happen before more tests are run at all.

2. **1b — Predicted-mass targeted re-mining of the existing raw MS data.**
   Zero new data, zero new wet lab. Build an exact-mass/MS2-fragment-ladder
   target list for the full carotenogenesis pathway (not just torulene/
   torularhodin) and check it against the raw feature list directly —
   tests whether the untargeted pipeline's own deduplication/SIRIUS
   annotation step silently orphaned or miscollapsed the real pigment
   signal before concluding the null reflects true biology.

3. **3a — Candidate carotenoid pathway genotyping (BFD panel).**
   Skips the metabolome's 4.8%-annotation bottleneck entirely by going
   straight from known pathway genes (CrtYB/CrtI-family Pfam domains,
   already in the fully-populated BFD tables) to color, using the
   already-built PGLS/species-level infrastructure. This is also the
   **first-ever genome↔color test in this project** (Phase 5 of the
   strategy doc, never reached) — high information value for low effort.

4. **4a — Variance-component/heritability decomposition of color.**
   A basic precondition check nobody has run: how much of the color
   variance is even heritable/strain-consistent (vs. plate/batch/
   measurement noise) before spending more effort hunting for what
   explains it? Cheap, and directly informs whether continuing to chase
   metabolome/genome correlates is well-motivated at all.

5. **2a — Molecular-network-aware compound family scores.**
   The EverythingBagel pipeline already computed a full GNPS-style
   molecular network (`compound_network.graphml`, `compound_pairs.tsv`)
   that this project has never used. Aggregating to chemical-family/
   community level (flagging the carotenoid/terpenoid cluster
   specifically) and rerunning the *already-validated* permutation +
   decoy stack on cluster scores is a moderate-effort way to use the 95%
   of features SIRIUS never annotated.

6. **6b / 3b combined direction — genotype-to-color test, hypothesis-free
   backup.** Once 3a's candidate-gene test runs, 6b's point (species/clade
   currently only treated as a nuisance to block out, never as a genuine
   causal pathway from genotype→color) argues for *also* running an
   unbiased genome-wide scan (3b) as a complement, not instead of, the
   candidate-gene approach — same logic as this project's existing
   primary/exploratory split convention.

7. **5b / 7a (near-duplicate, pick one) — Environment/ecology as covariate.**
   Both the evolutionary biologist and the ecologist independently flagged
   the same unused Origin/Environment metadata fields and the same
   photoprotectant/desiccation-resistance hypothesis for carotenoid
   function. Low effort, reuses existing PGLS scaffolding, and is a
   genuinely untested axis — color may track ecological niche rather than
   phylogeny or metabolome at all.

**Not recommended to start with** (high effort / high risk / low near-term
payoff, per the panel's own honesty about feasibility): 2b (research-grade
spectral kernel), 5a (convergence test, inherently underpowered at
N=3-4 independent origins), 7b (needs a new niche-breadth metric to be
invented), 3b as a *starting* point (better as a follow-up to 3a).

## Next steps

- Promote items 1, 3, 4 (all Low effort, no new data) to `todo/` items?
- Item 2 (area causal role) is arguably a blocker worth resolving before
  any of the others are fully trusted — PI's call on sequencing.
- None of these require the SIRIUS re-run to complete except deeper
  compound-ID work under idea 2a/2b — everything else can start now.
