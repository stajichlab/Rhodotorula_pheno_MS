---
session_id: 2026-08-17-002
project: rhodotorula-pheno-ms
branch: "main"
started: 2026-08-17T13:05:00-0700
ended: 2026-08-17T14:10:00-0700
duration_minutes: 65
files_changed: .living/findings/abundance-axis-growth-rate-relationships.md,.living/findings/FINDINGS_REGISTRY.md,.living/log/2026-08-17-002-rhodotorula-pheno-ms.md,analysis/class_level_aggregation/CLASS_LEVEL_AGGREGATION.md
---

## Session Log

### 13:05 — Continued from 2026-08-17-001
- Question from PI: the abundance axis must relate to liquid growth rate (Cu AUC); is there something else to learn of this relationship?

### 13:10 — Debugged NaN in growth-rate Spearman
- First attempt returned `nan` — two causes: (a) constant features → std=0 → NaN z (fixed with span>0 filter before z-scoring); (b) string-vs-int row-ID mismatch in the annotation join (fixed with `.astype(int)`).
- Clean read via `sample_metadata.csv.gz` (has `mean_auc_rate` + `canonical_strain` + fraction) + TSS-normalized rep-only feature matrix collapsed by `canonical_strain`.

### 13:20 — Abundance axis vs growth: negative/weak
- Whole-panel: `spearman(mean z, Cu-AUC)` cell = -0.174 (p=4e-3), supernatant = 0.007 (ns).
- Within *R. mucilaginosa*: cell = -0.070 (p=0.31), supernatant = 0.101 (p=0.15) — the whole-panel signal is between-species (species differ hugely in AUC: mucilaginosa ~23.3 vs kratochvilovae ~7.2).
- Colony `area` vs Cu-AUC: rho=0.37 whole-panel, 0.25 within-muc — related phenotypes, but the abundance axis tracks area (within-species) and NOT AUC. Recorded as F-018.

### 13:30 — Feature-level: cell hits are species-driven, supernatant has a real within-species set
- Cell: 218/10,230 features \|rho\|>0.3 vs AUC all-species (ceramides, wax esters, sphingolipids, di/tripeptides; rho≈-0.3-0.37). Within-muc: **0/10,230**.
- Supernatant: 143/10,439 all-species; within-muc **10 features**, vs permutation null mean 0.13 (200 draws, empirical p≈0) — significantly above null, ALL positive (more abundant in faster growers): puient nucleoside x4, di/tripeptides x3, 3 unannotated. Recorded as F-019.

### 13:45 — Writeup
- Created `.living/findings/abundance-axis-growth-rate-relationships.md` (F-018, F-019 with evidence ledgers + open questions).
- Updated FINDINGS_REGISTRY (F-018, F-019).

### 14:05 — Species-block + high-power permutation resolution
- No within-strain replication exists in this dataset (1 sample/strain), so per-sample treatment ≡ strain-level for the muc-supernatant set.
- Whole-panel within-species-block perm (300 perms, shuffle AUC only among strains of the same species, preserving cross-species structure):
  - cell: obs 218 hits vs null 125.7±62.8 (95th 240), p=0.090 → NOT beyond species structure
  - supernatant: obs 143 vs null 40.5±38.4 (95th 120), p=0.027 → modestly beyond; excess = the 10 within-species hits
- Within-muc 1000-perm: cell 0 hits (null mean 0.06, emp p=1.0); supernatant 10 (null mean 0.06, emp p<0.001).
- Conclusion: cell-fraction growth-metabolome association is cross-species only; supernatant has a small real within-species component (10 nucleoside/peptide features). Added F-003 to topic doc, upgraded F-019 to supported in FINDINGS_REGISTRY.

### Caveat logged
- Early top-feature rho≈-0.5 came from a near-constant-feature artifact (tiny-norm columns); scipy-verified real top is rho≈-0.37. The 10 supernatant hits were computed via vectorized rank correlation and verified against scipy spearmanr on spot checks.
