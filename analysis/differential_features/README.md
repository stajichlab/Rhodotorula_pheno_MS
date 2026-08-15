# Differential Features Analysis

Pairwise species comparisons of MS2 features, annotated with
GNPS library matches and SIRIUS/CANOPUS predictions.

- **11,953** significant features (FDR < 5%)
- **10** comparisons with significant features (out of 110 total)
- **790** SIRIUS annotations (663 with structure predictions)

## Cross-comparison rollup

| View | Features | Link |
|------|----------|------|
| All significant features (every comparison concatenated) | 11,953 | [open](all_significant_features_summary.html) |

## Cell pellet comparisons

| Comparison | Features | Identified | Dashboard |
|------------|----------|------------|-----------|
| Cell pellet: R. mucilaginosa vs R. toruloides | 2,600 | 953 | [dashboard](./cell_mucilaginosa_vs_toruloides/dashboard.html) |
| Cell pellet: R. diobovata vs R. mucilaginosa | 2,085 | 676 | [dashboard](./cell_diobovata_vs_mucilaginosa/dashboard.html) |
| Cell pellet: R. mucilaginosa vs R. paludigena | 1,270 | 358 | [dashboard](./cell_mucilaginosa_vs_paludigena/dashboard.html) |
| Cell pellet: R. mucilaginosa vs R. taiwanensis | 961 | 276 | [dashboard](./cell_mucilaginosa_vs_taiwanensis/dashboard.html) |
| Cell pellet: R. mucilaginosa vs R. sphaerocarpa | 642 | 201 | [dashboard](./cell_mucilaginosa_vs_sphaerocarpa/dashboard.html) |
| Cell pellet: R. dairenensis vs R. mucilaginosa | 268 | 66 | [dashboard](./cell_dairenensis_vs_mucilaginosa/dashboard.html) |

## Supernatant comparisons

| Comparison | Features | Identified | Dashboard |
|------------|----------|------------|-----------|
| Supernatant: R. diobovata vs R. mucilaginosa | 1,706 | 512 | [dashboard](./supernatant_diobovata_vs_mucilaginosa/dashboard.html) |
| Supernatant: R. mucilaginosa vs R. taiwanensis | 859 | 264 | [dashboard](./supernatant_mucilaginosa_vs_taiwanensis/dashboard.html) |
| Supernatant: R. mucilaginosa vs R. paludigena | 782 | 238 | [dashboard](./supernatant_mucilaginosa_vs_paludigena/dashboard.html) |
| Supernatant: R. mucilaginosa vs R. sphaerocarpa | 780 | 231 | [dashboard](./supernatant_mucilaginosa_vs_sphaerocarpa/dashboard.html) |

## Individual comparison tables

Each comparison directory also contains:
- `compound_summary.html` -- sortable/filterable table of significant features
- `volcano.pdf` / `volcano.png` -- volcano plot
- `top_features.pdf` / `top_features.png` -- top features plot
- `differential_features.csv.gz` -- full differential features table
