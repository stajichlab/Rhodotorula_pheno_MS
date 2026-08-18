---
session_id: 2026-08-17-001
project: rhodotorula-pheno-ms
branch: "main"
started: 2026-08-17T12:08:33-0700
ended: 2026-08-17T13:05:00-0700
duration_minutes: 57
files_changed: analysis/scripts/class_level_association.py,analysis/class_level_aggregation/outputs/class_association_{area,a,C}.csv,analysis/class_level_aggregation/outputs/class_enrichment_{area,a,C}.csv,analysis/class_level_aggregation/CLASS_LEVEL_AGGREGATION.md,analysis/ANALYSIS_MANIFEST.md,.living/findings/FINDINGS_REGISTRY.md,.living/findings/phenotype-metabolome-association-statistical-power.md,.living/log/2026-08-17-001-rhodotorula-pheno-ms.md
---

## Session Log

### 12:08 — Session started
- Branch: `main`
- Resuming from: 2026-08-16-001-rhodotorula-pheno-ms.md
- Context: prior session designed + wrote `analysis/scripts/class_level_association.py` (Strategies 2+3), compile-clean but unrun.

### ~12:15 — class_level_association.py: smoke test + bug fixes
- First run (decoy, n-perm 10) surfaced 2 bugs:
  1. GSEA `es_over_ordering` used `n_all = n_features` (unfiltered count) instead of the valid-feature ranking length → broadcast error. Fixed inside the function (uses `order_desc.size`).
  2. `pd.unique(list)` deprecation → `np.unique`.
- Deeper bug: **constant-feature filter false negatives** — 6 features passed `za.std(axis=0) > 0` with std≈7e-18 (float noise), yielding NaN rho. Replaced with a span/scale tolerance filter (`col_span > eps*scale*100`) and reordered so the NaN-column mask is applied to `members` too (that fix required the new `members_hv` intermediate). Runtime ~7-15s for all 6 fraction×ontology combos at n-perm=500.
- Fixed gate/target naming mismatch: decoy run writes `class_association_area.csv` but gate checked `class_association_area_decoy.csv` → gate updated to the actual filename.

### ~12:35 — Rerun: full decoy + calibration check
- Run `--predictor area` (n-perm 500): **129/502 classes at FDR<0.05** (e.g. 43% of cell/NPC classes). Suspicious — verified NOT a machinery bug:
  - Globally-shuffled true-null phenotype through the identical pipeline → expected ~3-5% FPR (0-6.7% across the 6 cells), median |rho|≈0.05.
  - Per-species blocks barely reduce decoy hits (43%→41%) → the signal is within-species, not species-tree-clade removable.
  - Root cause: global biomass/abundance artifact — `spearman(mean z-abundance, area) = -0.30` (p=3e-7), consistent with the existing F-001 biomass-scaling finding. This is why *every* class shows similar rho≈0.29.
- Conclusion recorded: the area decoy saturates via the known abundance artifact; the permutation machinery is well-calibrated.

### ~12:45 — Full runs: predictors a\* and C\*
- `--predictor a` (n-perm 500): 1/502 classes at FDR<0.05 (NPC pathway, supernatant); top nominal hits are all single-member classes (n_members=1) → feature-level effects, don't survive BH-FDR.
- `--predictor C` (n-perm 500): 0/502 at FDR<0.05.
- Consistent with Phase 2 feature-level result (top features rho≈0.33-0.37 but FDR≈0.24-0.43) — class-level aggregation does **not** rescue color signal (SIRIUS classes mostly n_members=1, so little pooling / little multiple-testing reduction).

### ~12:55 — Writeup
- Wrote `analysis/class_level_aggregation/CLASS_LEVEL_AGGREGATION.md` (design, per-predictor result table, calibration-verification section, interpretation).
- Updated `analysis/ANALYSIS_MANIFEST.md` with a `class_level_aggregation` section.
- Added finding F-017 to `.living/findings/FINDINGS_REGISTRY.md` and F-007 section in `.living/findings/phenotype-metabolome-association-statistical-power.md`.
- Files changed: class_level_association.py, CLASS_LEVEL_AGGREGATION.md, ANALYSIS_MANIFEST.md, FINDINGS_REGISTRY.md, phenotype-metabolome-association-statistical-power.md, 6 output CSVs.
