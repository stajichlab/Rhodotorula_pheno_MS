# Class-Level Aggregation: Strategies 2 + 3

**Script**: `analysis/scripts/class_level_association.py` (one shared code path for all three strategies)
**Inputs**: `analysis/linked_data/feature_abundance_matrix.csv.gz`, `sample_metadata.csv.gz`, `ms_feature_dedup_groups.csv`, `analysis/integrated_analysis/phase1_phenotype/strain_phenotype_table.csv`, `species_tree.nwk`, `analysis/sirius_annotation/sirius_annotations.tsv`
**Outputs**: `outputs/class_association_{area,a,C}.csv` (strategy 2) and `class_enrichment_{area,a,C}.csv` (strategy 3)
**Date**: 2026-08-17

## What this is

A second-pass analysis that moves up one level from the Phase 2 feature-level scans
(`phase2_metabolome_phenotype/`) by aggregating the 16,332 raw features (10,949 after
adduct/isotopologue de-duplication) into **chemical classes** using SIRIUS annotations, so we
can ask whether *whole compound classes* associate with the color phenotypes a\*/C\* — which is
more interpretable and statistical more powerful than the feature level when within-class
signals are weak but pervasive.

Only the *annotated* subset is used (~29% of features): 3,217 features with an NPC pathway and
NPC class label, 3,131 with a ClassyFire class. The background universe for each ontology is
therefore "features annotated with that ontology".

## Design

- **Predictors**: a\* (primary), C\* (secondary), area (negative control / decoy; gated).
- **Fractions**: `cell` (whole-cell extract) and `supernatant`, analyzed separately.
- **Replicates**: replicate samples of the same strain are collapsed by tolerance; the input
  sample matrix is TSS-normalized per sample; features are z-scored per feature.
- **Strategy 1 (basis)**: feature de-duplication already handled upstream —
  `ms_feature_dedup_groups.csv`. Each dedup representative is a feature row here
  (`is_group_representative`).
- **Strategy 2 (class score)**: for each class, score each strain as the mean of its member
  features' z-scores, then Spearman-rank-correlate the class score against the ranked phenotype.
- **Strategy 3 (GSEA-style)**: rank all features by their signed correlation with the phenotype,
  compute a Kolmogorov-Smirnov running-sum enrichment score (ES) for each class, and test the ES
  against the same permutation null.
- **Null / significance**: block-restricted permutation of the phenotype within species-tree
  clades (6 clades by default, seed 0, 500 perms), giving an empirical fraction `(excc+1)/(nperm+1)`
  for each test, then BH-FDR across classes within each fraction/ontology.
- **Hard gate** (mirrors Phase 2): a\*/C\* runs refuse to start unless a fresh `area` decoy
  output exists and is newer than the input data. Verified working in this run.

## Results

| predictor | strategy 2 hits (FDR<0.05) | strategy 3 mean-rho hits | strategy 3 GSEA-ES hits |
|-----------|----------------------------|--------------------------|--------------------------|
| area (decoy) | 129 / 502 classes | many | many |
| a\* | 1 / 502 (1 NPC pathway, supernatant) | 0 | 0 |
| C\* | 0 / 502 | 0 | 0 |

### a\* top signals (not FDR-significant)

The strongest a\* associations are overwhelmingly **single-feature classes** (class size = 1),
i.e. the class label is just wrapped around one feature's rho; after BH-FDR they do not survive.
Best nominal p-values: Epothilones (NPC class, cell, mean-rho emp. p=0.004 unadjusted),
Pyrroloquinoline alkaloids (0.008), Stilbenes/Monomeric stilbenes (0.018). No coherent
pigment-chemistry story at class level.

### Area decoy saturates — expected, not a calibration failure

The area decoy returns 129/502 classes significant at FDR<0.05 and a median |class rho| of
0.09-0.21. Three checks establish this is real biology, not a bug:

1. **Permutation machinery is calibrated.** Substituting a globally-shuffled (truly null)
   phenotype through the identical pipeline yields the expected ~3-5% of classes at
   p<0.05 (range 0-6.7% across the 6 fraction×ontology cells, median |rho|≈0.05).
2. **The association is not removable by species structure.** Restricting the permutation to
   per-species blocks instead of 6 clades barely changes the hit rate (43%→41% in cell/NPC) —
   the area signal is *within*-species, not a species-tree artifact. (Note: cell-fraction
   clade_1 holds 230/274 strains, so the default 6-clade restriction is weak in the cell
   fraction regardless.)
3. **Global abundance confound.** `spearman(mean z-abundance, area) = -0.30` (p=3×10⁻⁷) across
   strains in the cell fraction. Growth rate (AUC in Copper media) tracks overall measured
   metabolome scalar magnitude, so most classes pick up the same signed rho. This is the
   documented biomass-scaling artifact for cell-extract data
   (`.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md`, Phase 2).

### Bottom line

Partially consistent with the project's repeated null: a\* and C\* show **no class-level
association that survives block-permuted BH-FDR**, matching the Phase 2 feature-level result
(even the top Phase 2 features — rho ≈ 0.33-0.37 — had FDR ≈ 0.24). The class-level aggregation
does **not** rescue color signal, in part because SIRIUS classes are mostly tiny (typical
n_members=1) so aggregation neither pools weak signals nor reduces the multiple-testing burden
effectively. The decoy again demonstrates the study *has* sensitivity (129 classes), the
sensitivity is aimed at the global biomass axis in cell extracts, and color is not riding that
axis in any interpretable fraction.

## Follow-up: is the abundance axis the same thing as growth rate?

Since the decoy saturates via scalar abundance (`mean z ≈ -0.30` vs colony `area`), we asked
whether the same axis also tracks liquid growth (Cu AUC). It does **not** within species:
`spearman(mean z, Cu-AUC)` = -0.17 whole-panel (cell) but -0.07 within *R. mucilaginosa*; the
whole-panel value is between-species drift (species differ ~3× in AUC). Colony `area` and
Cu-AUC do correlate (rho≈0.37 / 0.25 within-species) — they are related phenotypes but recruit
the metabolome on different axes. At the feature level, the cell-fraction growth signal
(218/10,230 features |rho|>0.3, lipid chemistry) is entirely between-species (0/10,230 within
*mucilaginosa*); the only within-species growth signal found so far is a small supernatant set
of 10 features (purine nucleosides + di/tripeptides, all positive, p≈0 vs permutation null).
Full characterization in
`.living/findings/abundance-axis-growth-rate-relationships.md` (F-018, F-019).

## Interpretation / recommendation

- No new class-level color–metabolome association to follow up directly.
- The dominant remaining signal in this data (cell fraction) is the global biomass/abundance
  artifact, which continues to argue for per-sample abundance normalization/calibration before
  any further correlational mining, consistent with earlier findings.
- For discovery going forward the higher-value paths remain the targeted ones already flagged
  (HPLC/LC-UV-Vis pigment validation; validated candidate gene ↔ color tests) rather than
  further aggregation over the untargeted SIRIUS labels.

## Re-runs

```bash
python3 analysis/scripts/class_level_association.py --predictor area --n-perm 500 --seed 0
python3 analysis/scripts/class_level_association.py --predictor a    --n-perm 500 --seed 0
python3 analysis/scripts/class_level_association.py --predictor C    --n-perm 500 --seed 0
```
