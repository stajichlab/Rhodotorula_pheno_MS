# Differential Features Analysis

Pairwise species comparisons of MS2 features, annotated with
GNPS library matches and SIRIUS/CANOPUS predictions.

- **23,745** significant features (FDR < 5%)
- **11** comparisons with significant features (out of 110 total)
- **790** SIRIUS annotations (663 with structure predictions)

## Cross-comparison rollup

| View | Features | Link |
|------|----------|------|
| All significant features (every comparison concatenated) | 23,745 | [open](all_significant_features_summary.html) |

## Cell pellet comparisons

| Comparison | Features | Identified | Dashboard |
|------------|----------|------------|-----------|
| Cell pellet: R. mucilaginosa vs R. toruloides | 2,439 | 2,167 | [dashboard](./cell_mucilaginosa_vs_toruloides/dashboard.html) |
| Cell pellet: R. diobovata vs R. mucilaginosa | 2,072 | 1,928 | [dashboard](./cell_diobovata_vs_mucilaginosa/dashboard.html) |
| Cell pellet: R. mucilaginosa vs R. paludigena | 1,321 | 1,214 | [dashboard](./cell_mucilaginosa_vs_paludigena/dashboard.html) |
| Cell pellet: R. mucilaginosa vs R. taiwanensis | 1,081 | 1,003 | [dashboard](./cell_mucilaginosa_vs_taiwanensis/dashboard.html) |
| Cell pellet: R. mucilaginosa vs R. sphaerocarpa | 637 | 598 | [dashboard](./cell_mucilaginosa_vs_sphaerocarpa/dashboard.html) |
| Cell pellet: R. dairenensis vs R. mucilaginosa | 276 | 262 | [dashboard](./cell_dairenensis_vs_mucilaginosa/dashboard.html) |

## Supernatant comparisons

| Comparison | Features | Identified | Dashboard |
|------------|----------|------------|-----------|
| Supernatant: R. diobovata vs R. mucilaginosa | 1,705 | 1,417 | [dashboard](./supernatant_diobovata_vs_mucilaginosa/dashboard.html) |
| Supernatant: R. mucilaginosa vs R. taiwanensis | 1,261 | 1,056 | [dashboard](./supernatant_mucilaginosa_vs_taiwanensis/dashboard.html) |
| Supernatant: R. mucilaginosa vs R. paludigena | 1,079 | 863 | [dashboard](./supernatant_mucilaginosa_vs_paludigena/dashboard.html) |
| Supernatant: R. mucilaginosa vs R. sphaerocarpa | 748 | 692 | [dashboard](./supernatant_mucilaginosa_vs_sphaerocarpa/dashboard.html) |

## Supernatant vs cell pellet (paired, within species)

| Comparison | Features | Identified | Dashboard |
|------------|----------|------------|-----------|
| R. mucilaginosa: supernatant vs cell pellet | 11,126 | 6,421 | [dashboard](./mucilaginosa_sup_vs_cell/dashboard.html) |

## Individual comparison tables

Each comparison directory also contains:
- `compound_summary.html` -- sortable/filterable table of significant features
- `volcano.pdf` / `volcano.png` -- volcano plot
- `top_features.pdf` / `top_features.png` -- top features plot
- `differential_features.csv.gz` -- full differential features table
