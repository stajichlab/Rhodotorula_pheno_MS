# Genome amino-acid composition vs. non-stress growth/color phenotypes

## Goal

Independent baseline for `../../copper/genome_characteristics.md`: repeat the same
proteome-wide-AA-composition-vs-phenotype test using YPD2 (non-stress, plain-medium)
phenotypes instead of copper-stress AUC. Because these phenotypes are not
copper-related, any AA association here is either (a) a generic life-history/
growth-rate correlate not specific to metal stress, or (b) evidence that the
join/PGLS pipeline itself produces false positives on unrelated traits — both
useful to know before trusting the copper result.

## Data & phenotypes

`data/metadata/EXFAB_UCR-005/YPD2_phenotypic.20260702.fixed.csv.gz`, 316 unique
strains under plain YPD medium (`Media == "YPDN"`), four phenotypes taken as-is:

- `Mean_Shape_Area` — mean colony area across imaging timepoints, used as a
  growth-size proxy (not a fitted growth rate; see limitations).
- `Mean_ColorLab_L*Mean`, `_a*Mean`, `_b*Mean` — CIE-Lab colony color (L*=
  lightness, a*=green-red, b*=blue-yellow), mean across replicates.

A strain with multiple metadata rows (different plates) is collapsed to the mean
across rows before joining (only 2/316 strains had >1 row in practice).

## Method

Identical pipeline to the copper analysis: join to `BFD/results/genome_stats/aa_freq`
via the phyling tree using the shared `analysis/copper/scripts/common.py` matcher,
then naive Pearson correlation (BH-FDR across 20 AAs, per phenotype) followed by
PGLS with ML-estimated Pagel's λ (BH-FDR per phenotype). See
`../../copper/genome_characteristics.md` for the full rationale (same assumptions
and zero-branch-length handling apply here). Run: `bash scripts/run.sh`.

**Coverage:** 145/316 YPD2 strains (46%) have both a tree placement and an aa_freq
file (`outputs/join_diagnostics.txt`). Species composition is similarly skewed:
101/145 (70%) _R. mucilaginosa_.

## Results (n = 145 strains, 2026-08-15 run)

| Phenotype | Naive BH-significant AAs | Survive PGLS (BH q<0.05) |
|---|---|---|
| `Mean_Shape_Area` (growth proxy) | 14/20 | **0** |
| `Mean_ColorLab_L*Mean` (lightness) | 3/20 (D, V, E) | **1** (D) |
| `Mean_ColorLab_a*Mean` (green-red) | 15/20 | **0** |
| `Mean_ColorLab_b*Mean` (blue-yellow) | 7/20 | **0** |

Full tables: `outputs/naive_correlation_results.csv`, `outputs/pgls_correlation_results.csv`,
`outputs/naive_vs_pgls_comparison.csv`.

**This is a stronger and cleaner negative result than the copper case.** For
growth size and two of three color axes, essentially all of the naive AA
correlation vanishes under phylogenetic correction — i.e. AA composition tracks
species identity, and species identity tracks colony size/color (unsurprising:
these are taxonomically diagnostic traits), but there is no evidence of a
trait-level link between bulk proteome AA composition and growth size or hue
once lineage is accounted for. Lightness (L*) retains one weak hit (Asp), not
pursued further here.

**Comparison to copper:** the copper analysis found 5-7 AAs (led by Ser and,
after correction, Cys) surviving PGLS out of 20 — a much larger fraction than any
YPD phenotype here. That contrast is mildly reassuring that the copper signal is
not simply "PGLS always leaves some AAs significant by chance/model artifact" —
but the copper result's own within-species sensitivity check
(`../../copper/genome_characteristics.md`) still failed to replicate, so neither
result should be over-interpreted from this alone.

## Limitations

- `Mean_Shape_Area` is a snapshot mean, not a fitted growth-rate/AUC parameter
  (unlike the copper phenotype) — it is a size, not a rate, proxy.
- Same species imbalance and zero-branch-length caveats as the copper analysis
  apply; no within-species sensitivity check was run here (lower priority than
  copper given the near-uniformly null result).
- No multiple-testing correction was applied *across* the four phenotypes
  (each phenotype's 20-AA family is BH-corrected independently); treat the 4x20
  results as a screen, not a confirmatory family-wise test.
