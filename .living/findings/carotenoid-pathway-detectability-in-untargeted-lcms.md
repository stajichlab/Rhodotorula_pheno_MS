---
topic: carotenoid-pathway-detectability-in-untargeted-lcms
description: Whether the carotenogenesis pathway (phytoene through torularhodin) is detectable at all in this project's untargeted LC-MS2 data, and whether SIRIUS's compound annotations for candidate features are trustworthy.
created: 2026-08-15
last_updated: 2026-08-15
status: active
---

# Carotenoid pathway detectability in untargeted LC-MS2

## F-001: Exact-mass search finds plausible carotenoid-pathway feature matches, including one SIRIUS actively misannotated, but none correlate with color
**Status:** preliminary
**Claim:** A targeted exact-mass search (20 ppm, [M+H]+/[M+Na]+/[M+NH4]+/
[M+H-H2O]+ adducts) of the full carotenogenesis pathway (phytoene through
torularhodin) against the raw EB feature table (53,040 features, not the
16,332-feature has-MS2 subset) found 11 raw-feature matches across 6 of 9
pathway compounds (no matches for neurosporene/lycopene/γ-carotene/
β-carotene). Two stand out: a phytofluene [M+NH4]+ candidate (row 1735,
-1.0 ppm, 220 total scans -- the highest-intensity match by far,
unannotated by SIRIUS) and a torularhodin [M+H]+ candidate (row 21315,
+2.5 ppm, 6 total scans, **which SIRIUS annotated as "Polyamines" / a
chemically implausible cyclopentane-dicarboxamide structure** -- a
concrete, not hypothetical, instance of SIRIUS misannotating a
carotenoid-mass feature). Both candidates' dedup groups were already
inside Phase 2's tested set and show no meaningful a\* correlation
(Spearman rho 0.03-0.12, FDR 0.48-0.94 in both fractions).
**Implications:** This is genuinely double-sided evidence. On one hand,
it strengthens the case that this project's ESI-based untargeted method
*can* detect signal in the carotenoid mass range at all (ruling out the
strongest form of "ESI can't ionize this compound class") and confirms a
concrete SIRIUS misannotation risk for at least one candidate feature
(not just a theoretical concern). On the other hand, neither standout
candidate's abundance correlates with color, so even if these ARE real
carotenoid-pathway intermediates, they don't explain the color-metabolome
null result -- this re-mining pass did not rescue a hidden signal.
Structural confirmation (MS2 fragmentation check for polyene-
characteristic neutral losses) has not yet been done for either
candidate; without it, "plausible mass match" should not be upgraded to
"confirmed carotenoid feature."
**Tags:** carotenoids, untargeted-metabolomics, SIRIUS-misannotation, exact-mass-search, phytofluene, torularhodin, null-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | Targeted exact-mass re-mining (idea1_targeted_mass_remining.py) | EB aligned_features_ms2.csv (53,040 raw features), cross-referenced against Phase 2's color association results | Rhodotorula_pheno_MS | 11 mass matches found (2 high-confidence), neither correlates with a\* | refines |

### Open Questions
- Why zero matches for neurosporene/lycopene/γ/β-carotene specifically --
  genuinely absent, or a limitation of the 20 ppm / 4-adduct search
  scope?
- Does the RT pattern of the 11 matches (5.0-7.6 min, not a clean
  monotonic desaturation-order series) support or undermine treating them
  as a real biosynthetic pathway signal?
- Does a genuine reference-spectrum comparison (e.g. GNPS library search)
  for row 21315 confirm or refute the torularhodin-family identification
  its hand-inspected fragmentation pattern is consistent with?

## F-003: Expanded search (apocarotenoids, additional carotenoids, sterol pathway) finds a candidate ergosterol feature correlated with copper-resistance AUC (naive, not yet validated)
**Status:** preliminary
**Claim:** Expanding the exact-mass search to 19 compounds across 4
categories (original carotenoid pathway, apocarotenoids, additional
oxygenated carotenoids, sterol-pathway precursor-competition markers)
found 31 raw-feature matches (20 ppm). The most chemically self-
consistent new match is an ergosterol [M+H-H2O]+ candidate (row 846,
-4.3 ppm) that EB's own pipeline independently flagged as an in-source
fragment -- exactly the relationship expected for a sterol's dominant
dehydration ion, two independent lines of evidence agreeing. This
candidate (dedup group asid_31011) shows no color correlation (a\*
Spearman rho 0.02/-0.06) but a naive, uncorrected correlation with
copper-resistance growth rate (mean_auc_rate) in the cell fraction:
Spearman rho=0.226, p=0.0002, n=267. An astaxanthin [M+H]+ candidate
(row 9384, +2.4 ppm, 37 scans) was also found but shows no correlation
with either phenotype.
**Implications:** The ergosterol-AUC correlation is mechanistically
plausible (ergosterol is the principal fungal membrane sterol; membrane
sterol content is an established determinant of metal-ion tolerance) and
worth formal testing, but per this project's established pattern
(phylogenetic-confounding-of-trait-molecular-associations F-005: 2-for-2
track record of naive whole-panel hits -- including this *same* AUC
phenotype's amino-acid correlations -- failing phylogenetic block-
permutation or within-species restriction), **this is explicitly a lead,
not a result**, until run through the same rigor (phylogenetic block
permutation, negative control) already standard elsewhere in this
project. Not yet done.
**Tags:** ergosterol, copper-resistance, sterol-pathway, precursor-competition, naive-correlation, needs-validation

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Expanded exact-mass search + AUC quick-check | Same feature table as F-001/F-002, cross-referenced against sample_metadata's mean_auc_rate | Rhodotorula_pheno_MS | Ergosterol candidate: naive rho=0.226, p=0.0002 (cell fraction) with copper-AUC; no color correlation for any new candidate | refines |

### Open Questions
- Does the ergosterol-AUC correlation survive phylogenetic block
  permutation (species-level primary, within-*R. mucilaginosa* secondary)
  and a hard-gated negative control, or does it join the copper-AUC
  amino-acid "hits" that didn't survive the same check?
- Is astaxanthin biologically plausible for *Rhodotorula* specifically, or
  should that candidate be deprioritized on genus-level chemistry grounds
  regardless of the mass match?

## F-004: A cluster of 4 independently-called ergostane-class features all show the same-direction naive copper-AUC correlation, strengthening (but not validating) the ergosterol lead
**Status:** preliminary
**Claim:** SIRIUS's NPC-class breakdown of the 73 "Terpenoids"-pathway
features found an Ergostane/Cholane/Cholestane steroid cluster (7
features) -- the single largest NPC subclass in that bucket, structurally
consistent with F-003's row-846 ergosterol ISF candidate. Naive Spearman
quick-check (`analysis/scripts/idea1_auc_quickcheck.py`) of the 3 most
chemically self-consistent SIRIUS ergostane calls (Peroxyergosterol,
Ergost-3,5,7,9(11),22-pentaen, 7-Hydroxyergosterol) plus the row-846
anchor against copper-resistance AUC found **all four show the same
positive-direction correlation in the cell fraction** (rho 0.17-0.28, all
p<0.006 naive). Row 6682 additionally flips to a significant negative
correlation in supernatant (rho=-0.24); row 35014 is positive in both
fractions.
**Implications:** Four independently-called features from the same
biosynthetic family agreeing in sign is stronger circumstantial support
than any single feature alone -- raises the prior that F-003's ergosterol
signal reflects real membrane-sterol/copper-tolerance biology rather than
one mass coincidence. **Still explicitly not validated**: no phylogenetic
block permutation or negative control has been run on any of these four
features yet, and this project has a 2-for-2 (soon 3-for-3?) track record
of naive whole-panel hits failing that check. See
`analysis/integrated_analysis/phase3_metabolome_phenotype_idea1/STEROL_CLUSTER_AUC_CHECK.md`.
**Tags:** ergosterol, sterol-pathway, copper-resistance, naive-correlation, needs-validation, multi-feature-corroboration

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Ergostane-cluster AUC quick-check (idea1_auc_quickcheck.py) | Same feature table as F-001/F-003, cross-referenced against sample_metadata's mean_auc_rate | Rhodotorula_pheno_MS | 4/4 ergostane-class features positive in cell fraction (rho 0.17-0.28, p<0.006 naive) | refines |

### Open Questions
- Does the row-846/cluster ergosterol-AUC signal survive phylogenetic
  block permutation + hard-gated negative control (the full test still
  not yet run)?
- Why does row 6682 flip sign between fractions while row 35014 doesn't --
  real differential retention/secretion biology, or feature-specific
  noise?

## F-002: MS2 fragmentation check demotes one standout candidate, leaves the other genuinely ambiguous
**Status:** preliminary
**Claim:** Manual inspection of the raw MS2 spectra (pulled from
`aligned_features_filled.mgf` by FEATURE_ID) for F-001's two standout
candidates gives opposite verdicts from what intensity/mass-accuracy
alone suggested. **Row 1735 (phytofluene candidate) is demoted**: its
spectrum is a largely undissociated precursor (base peak = precursor ion
itself) with only scattered low-intensity fragments and, critically, no
peak near the expected ~543.5 m/z NH3 neutral-loss diagnostic for a real
[M+NH4]+ adduct -- its earlier "highest intensity in the list" status
reflected precursor stability, not confirmatory fragmentation. **Row
21315 (torularhodin candidate) remains a live, unconfirmed candidate**:
its spectrum shows two clean, paired water-loss fragments (299.30->281.29
and 327.34->309.31, both Δ18.01) consistent with a carboxylic-acid/
hydroxyl-bearing structure (torularhodin has a terminal -COOH), plus two
substantial (not noise-level) mid-mass fragments (282.28, 310.31)
plausible for mid-backbone cleavage of a C40 polyene. Not confirmed
against any reference spectrum.
**Implications:** Overall spectral intensity/total_scans is not a
reliable proxy for identification confidence in either direction here --
the highest-intensity candidate turned out to be the weaker
identification, and the lower-intensity candidate the more chemically
coherent one. This is a useful general caution for any future targeted
re-mining pass in this project. Net effect on Idea 1: still no rescued
color signal (row 21315 doesn't correlate with a\* either, per F-001), but
row 21315 is now the sole remaining candidate worth a genuine reference-
spectrum comparison (e.g. GNPS library search) before any wet-lab
escalation decision.
**Tags:** carotenoids, ms2-fragmentation, torularhodin, phytofluene, spectral-quality-caution

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-15 | Manual MS2 inspection (rows 1735, 21315) from aligned_features_filled.mgf | Same as F-001 | Rhodotorula_pheno_MS | Row 1735 demoted (no NH3 loss, noise-dominated); row 21315 chemically coherent but unconfirmed | refines |
