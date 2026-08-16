<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-08-16

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 0 sections | 2026-08-15 | — |
| decisions.md | 0 entries (large — read selectively) | 2026-08-16 | — |
| learnings.md | 4 entries | 2026-08-15 | pandas groupby silently drops rows with a NULL group key, NaN correlation values silently produce spuriously significant permutation p-values, `genome_strain_species_busco_map.csv` has no generating script — species-level analyses can't pick up new BFD taxa on their own, "Hard gate" negative-control checks need to verify statistical adequacy, not just file freshness |
| log/ | 3 sessions | 2026-08-15 | rhodotorula-pheno-ms (3) |
| findings/ | 9 findings across 5 topics | 2026-08-16 | carotenoid-pathway-detectability-in-untargeted-lcms, convergent-color-evolution-in-rhodotorula, phenotype-metabolome-association-statistical-power, phylogenetic-confounding-of-trait-molecular-associations, biomass-scaling-artifacts-in-extraction-based-metabolomics |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-08-16 (heuristic)

## Tag clusters

- **silent-failure** (3 entries) — L-1, L-2, L-4
- **permutation-test** (2 entries) — L-2, L-4

## Most recent (10)

- [2026-08-15] L-1: pandas groupby silently drops rows with a NULL group key
- [2026-08-15] L-2: NaN correlation values silently produce spuriously significant permutation p-values
- [2026-08-15] L-3: `genome_strain_species_busco_map.csv` has no generating script — species-level analyses can't pick up new BFD taxa on their own
- [2026-08-15] L-4: "Hard gate" negative-control checks need to verify statistical adequacy, not just file freshness

## By tag

- `silent-failure`: L-1, L-2, L-4
- `permutation-test`: L-2, L-4
- `busco-map`: L-3
- `data-loss`: L-1
- `groupby`: L-1
- `hard-gate`: L-4
- `missing-script`: L-3
- `nan-handling`: L-2
- `negative-control`: L-4
- `numpy`: L-2
- `pandas`: L-1
- `phenotype-ingestion`: L-1
- `phylogeny`: L-3
- `reproducibility`: L-3
- `smoke-test`: L-4
- `species-level-analysis`: L-3
- `statistical-power`: L-4
- `statistics`: L-2

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
