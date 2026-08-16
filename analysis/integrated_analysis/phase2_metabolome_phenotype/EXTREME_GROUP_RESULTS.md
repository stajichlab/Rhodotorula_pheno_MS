# Extreme-group (high vs. low orange/red) color↔metabolome test (2026-08-16)

PI follow-up question: instead of a continuous whole-panel correlation
(null, see `PHASE2_SUMMARY.md`), does partitioning strains into
high-orange/red vs. low-orange/red groups and comparing compound
abundances directly find anything? Explicitly acknowledged by the PI as
lumping species together and ignoring phylogeny.

**Script**: `analysis/scripts/extreme_group_color_association.py` — top
vs. bottom quartile (25%) on a\*/C\*/area, rank-biserial effect size
(Mann-Whitney-U equivalent, vectorized rank-sum), empirical p-value from
species-tree-**block-restricted** label permutation (500 perms) — so this
design is not fully blind to phylogeny either, it's a different,
higher-power question layered on the same block-permutation machinery as
the rest of Phase 2. Cell and supernatant always tested separately.

## Group composition
High/low a\* groups (73 high, 65 low strains) and C\* groups (71/63) are
dominated by *R. mucilaginosa* (34-60 strains per group) with the
remainder spread across most other species — see per-run species-count
breakdown printed to stderr / reproducible via the script. This is the
phylogenetic-lumping the PI flagged: a real signal here could still be
driven almost entirely by *R. mucilaginosa*'s own within-species
variation rather than a genus-wide pattern.

## Results

| Predictor | Fraction | n high / n low | BH-FDR<0.05 hits | Null-expected |
|---|---|---|---|---|
| area (decoy) | cell | 66/65 | **1,723** | ~0 |
| area (decoy) | supernatant | 67/65 | 0 | ~0 |
| a\* (primary) | cell | 73/65 | **0** | ~0 |
| a\* (primary) | supernatant | 73/65 | **0** | ~0 |
| C\* (secondary) | cell | 71/63 | **0** | ~0 |
| C\* (secondary) | supernatant | 72/63 | **0** | ~0 |

**Null for color, in both fractions, at quartile extremes.** The decoy
confirms the design has real power to detect an effect at this n (1,723
cell-fraction hits for area, consistent with the already-known
biomass/extraction confound — `.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md`),
so this isn't underpowered the way the small within-species sweep was.
**This is a 6th independent method (on top of the 5 in `PHASE2_SUMMARY.md`)
that finds no color↔metabolome signal**, now including a design
specifically built to catch threshold/nonlinear effects a continuous
correlation could miss.

## Output files
`extreme_group_association_{a,C,area_decoy}.csv` in this directory — one
row per deduplicated feature group × fraction, columns `rank_biserial`,
`empirical_p`, `empirical_fdr`.

## Reproduce
```
python3 analysis/scripts/extreme_group_color_association.py --predictor area   # negative control, run first
python3 analysis/scripts/extreme_group_color_association.py --predictor a      # primary
python3 analysis/scripts/extreme_group_color_association.py --predictor C      # secondary
```
