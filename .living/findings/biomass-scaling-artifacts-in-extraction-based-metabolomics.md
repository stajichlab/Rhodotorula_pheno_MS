---
topic: biomass-scaling-artifacts-in-extraction-based-metabolomics
description: Whether sample biomass/size (e.g. colony area, culture density) introduces broad, systematic bias into extraction-based untargeted metabolomics abundances even after standard normalization (TSS).
created: 2026-08-15
last_updated: 2026-08-15
status: active
---

# Biomass scaling artifacts in extraction-based metabolomics

## F-001: Colony size broadly confounds TSS-normalized cell-fraction (not supernatant) metabolite abundances
**Status:** preliminary
**Claim:** In an untargeted LC-MS2 metabolomics dataset of ~275 fungal
strains (cell pellet + supernatant fractions), colony area (a
biomass/growth proxy unrelated to color or pigment chemistry) showed a
statistically robust, broad association with TSS-normalized cell-fraction
metabolite abundances (1,524 of 10,164 deduplicated compound groups,
~15%, BH-FDR<0.05 under a phylogeny-aware permutation null that showed
~0 expected hits) but essentially no association in the paired
supernatant fraction from the same strains (0/10,416 hits). Effect sizes
among "hits" were modest and diffuse (median |Spearman rho| ≈ 0.23), not
concentrated in a few strongly-correlated compounds — consistent with a
systematic extraction/loading artifact scaling with cell biomass, rather
than a genuine biological signal.
**Implications:** Standard total-sum-scaling (TSS) normalization does not
fully remove a biomass-scaling confound in cell-pellet (intracellular)
extraction-based metabolomics, at least in this dataset/pipeline. Any
downstream analysis correlating a strain-level trait against cell-fraction
abundances should check whether that trait itself correlates with colony
size/biomass before trusting the result — a real trait-metabolite
association could be partly or wholly a biomass-scaling artifact if the
trait and colony size are themselves correlated (e.g., in the originating
project, growth rate/AUC phenotypes are also measured for these same
strains). Supernatant (secreted) metabolite fractions appear less
susceptible to this specific artifact, which is itself a useful,
reportable methodological observation for anyone measuring both fractions.
**Tags:** metabolomics, normalization, TSS, biomass, confound, extraction, cell-fraction, methodology

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | Phase 2 color-metabolome negative-control decoy run | EB_20260130_ExFAB_Rhodo_Sup_and_Pellet aligned_features_ms2 (10,949 deduplicated compound groups, ~275 strains) | Rhodotorula_pheno_MS | Colony area (decoy trait) associated with 1,524/10,164 cell-fraction compound groups at BH-FDR<0.05 (permutation null ~0 expected); 0/10,416 in supernatant | supports |

### Open Questions
- Is this a TSS-normalization insufficiency specifically, or would a
  different normalization (e.g. probabilistic quotient normalization,
  or explicit biomass-based scaling using the colony area itself)
  eliminate it?
- Does this artifact also appear in other extraction-based metabolomics
  datasets from the same EverythingBagel pipeline / same lab, or is it
  specific to this colony-plate-based sample-prep protocol?
- Which specific compound classes are most affected — is it uniform
  across all NPC pathways, or concentrated in particular chemistries
  (e.g. lipids, which would be consistent with a membrane/biomass-content
  scaling explanation)?
