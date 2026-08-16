# Ideas — Executive Summary

> Navigation document for `analysis/ideas/`. This directory holds
> persona-based ideation sessions and their resulting development plans.
> For where each idea's actual results landed, see
> `analysis/integrated_analysis/EXECUTIVE_SUMMARY.md`.

## 2026-08-15-color-metabolome-genome-null-brainstorm
**Path**: `2026-08-15-color-metabolome-genome-null-brainstorm/`
**Trigger**: Phase 2's 5-independent-method null result for color↔metabolome association (see integrated-analysis summary, Phase 2).
**Format**: 7 personas (natural-products chemist, metabolomics/MS specialist, fungal genomics bioinformatician, quantitative geneticist, evolutionary biologist, causal-inference researcher, ecologist) — `01`–`07_*.md`, synthesized in `00_index.md` (14 ideas total, ranked).

### Idea status

| # | Idea | Status | Where it landed |
|---|---|---|---|
| 1 | Targeted re-mining of raw MS data against a carotenogenesis-pathway mass list | **Done, live lead** | `analysis/integrated_analysis/phase3_metabolome_phenotype_idea1/` — no color signal, but an unvalidated ergosterol/sterol-cluster↔copper-AUC lead |
| 3 | Candidate carotenoid/melanin-pathway-gene genotyping vs. color (first genome↔color test) | **Screening done, association test not started** | `analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/` — pigment genes broadly present; laccase-route melanin confirmed 3 ways; sparse-presence families (`crtR`/`crtQ`/`hgd`) identified as the first testable candidates but not yet run against color |
| 5 | Bayesian regime-shift detection (formal replacement for the heuristic convergence screen) | **Done** | `analysis/integrated_analysis/phase5_genome_linkage/idea5_regime_shift/` — diffuse posterior (power ceiling at 17 tips), but top candidate corroborates Phase 1's original ranking |
| 6 | Genome-wide GWAS-style scan (ambitious tail of Idea 3) | **Not started** | Would require real variant/SNP calling (assembly-derived protein calls are sufficient for Idea 3's gene-level tier, not for this); infrastructure requirements scoped in `DEVELOPMENT_PLAN.md` but no work done |
| 2, 4, 7-14 | Determine colony-area's causal role; compound-class/pathway-level aggregation; molecular-network reuse; others (see `00_index.md` full ranked list) | **Not started** | — |

**Also in this directory**:
- `DEVELOPMENT_PLAN.md` — detailed plan for Ideas 1/3/5/6, including the resolution that SNP calling is NOT needed for Idea 3's gene-level tier (only for Idea 6's genome-wide tier).
- `IDEA1_CHEMIST_FRAMEWORK.md` — two-phase framework for Idea 1 (in-silico mass re-mining first, wet-lab escalation only if empty) — Phase 1 of this framework is what's been executed; wet-lab escalation not started.
