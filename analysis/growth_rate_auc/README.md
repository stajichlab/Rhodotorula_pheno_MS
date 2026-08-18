# Growth-rate (Cu-AUC) <> supernatant metabolome: within-species signal

**Script**: `analysis/scripts/build_growth_auc_feature_table.py`
**Inputs**: `analysis/linked_data/sample_metadata.csv.gz`, `feature_abundance_matrix.csv.gz`,
`ms_feature_dedup_groups.csv`, `analysis/sirius_annotation/sirius_annotations.tsv`, EB/GNPS
library search results
**Outputs** (this directory):
- `outputs/within_mucilaginosa_cell_auc_features.tsv.gz`, `within_mucilaginosa_supernatant_auc_features.tsv.gz` —
  every feature (adduct-dedup'ed, replicate-collapsed) with its within-*R. mucilaginosa*
  Spearman rho vs liquid growth rate (Cu-AUC `mean_auc_rate`).
- `outputs/compound_summary.tsv` — the **10 supernatant features** with |rho|>0.3, annotated
  with the exact same identity pipeline + column schema as the pairwise differential
  comparisons (`scripts/build_compound_summary.py`), plus the per-feature association
  columns (`rho_auc`, `emp_p_perm_within_muc`).
- `outputs/compound_summary.html` — sortable/filterable interactive table rendered by the
  same generator used for the differential comparisons
  (`scripts/generate_compound_table_html.py`).

**Why only these 10**: feature-level Spearman of each feature vs `mean_auc_rate`, restricted to
*R. mucilaginosa* (n=208-212, the only species with adequate sample size). Cell-fraction
features show **0** hits at |rho|>0.3; supernatant shows 10, and 1000-permutation testing
gives an empirical p of 0.000999 per hit (null mean 0.06 hits) — significantly more than
chance. All association columns are within-species; the species-level/whole-panel signal is
entirely between-species (phylogenetic) drift, see
`.living/findings/abundance-axis-growth-rate-relationships.md` (F-002/F-003).

**Chemistry of the 10** (all more abundant in faster growers except row 4299, which is
negative and unidentified): purine nucleosides (rows 48, 1483, 1762 — methylthioadenosine /
guanosine-like), di-/tri-peptides (rows 1979, 4951), an analog match to Cholesterol (22697),
an analog match to a chlorinated mass-bank dimer (28139), and 3 unidentified (4299 neg, 5327,
17148).

Regenerate:

```bash
python3 analysis/scripts/build_growth_auc_feature_table.py
```
