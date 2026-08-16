# Idea 1, Phase 1: Predicted-mass targeted re-mining (2026-08-15)

Script: `analysis/scripts/idea1_targeted_mass_remining.py`. Searches the
fuller EB feature table (53,040 raw features, not just the 16,332-feature
has-MS2 working subset) directly by exact mass for the full
carotenogenesis pathway (phytoene → phytofluene → ζ-carotene →
neurosporene → lycopene → γ/β-carotene → torulene → torularhodin),
independent of SIRIUS's own compound calling — per
`IDEA1_CHEMIST_FRAMEWORK.md` Phase 1.

## Result: 11 raw-feature matches within 20 ppm, across 9/9 compounds' adduct sets except neurosporene/lycopene/gamma/beta-carotene (0 hits for those 4)

| Compound | Adduct | row ID | ppm error | RT (min) | total_scans | SIRIUS call |
|---|---|---|---|---|---|---|
| phytofluene | [M+NH4]+ | 1735 | **-1.0** | 5.59 | **220** | none (unannotated) |
| torulene | [M+Na]+ | 16228 | -1.2 | 5.57 | 9 | none |
| torularhodin | [M+H]+ | 21315 | **+2.5** | 5.70 | 6 | "Polyamines" — implausible structure (a cyclopentane dicarboxamide) |
| torularhodin | [M+H]+ | 26444 | -11.5 | 6.36 | 11 | "Open-chain polyketides" — implausible steroid-like structure |
| torularhodin | [M+H-H2O]+ | 3297 | +14.7 | 5.70 | 65 | none |
| ζ-carotene | [M+Na]+ | 4076, 10220 | +6.3 to +6.7 | 6.97, 7.17 | 19, 20 | none |
| ζ-carotene | [M+H-H2O]+ | 30947 | +8.2 | 7.56 | 14 | none |
| phytoene | [M+H]+ | 35469 | -16.6 | 7.28 | 12 | none |
| phytoene | [M+NH4]+ | 44251 | -17.5 | 5.60 | 13 | none |
| phytofluene | [M+Na]+ | 15654 | -6.5 | 5.05 | 8 | none |

**No matches for neurosporene, lycopene, γ-carotene, or β-carotene** at
any adduct within 20 ppm.

## Two standout candidates worth follow-up

1. **Row 1735 (phytofluene, [M+NH4]+, -1.0 ppm)** — the tightest mass
   match in the list AND by far the highest intensity (220 total scans,
   vs. single digits to low tens for most others). This is the most
   credible individual candidate here. Unannotated by SIRIUS (part of the
   95%+ "dark matter" fraction) — exactly the kind of feature Idea 1's
   framework predicted the untargeted pipeline could have missed.
2. **Row 21315 (torularhodin, [M+H]+, +2.5 ppm)** — a tight match for
   the single most diagnostic compound in this pathway (torularhodin is
   the *crtS* product that differentiates orange from darker-orange/red
   coloration). **SIRIUS annotated this exact feature as "Polyamines" /
   a chemically implausible cyclopentane-dicarboxamide structure** — a
   concrete instance of the SIRIUS-misannotation risk the chemist persona
   flagged, not a hypothetical one. Low intensity (6 scans) tempers
   confidence, but the mass accuracy and the SIRIUS mismatch both make
   this worth a manual MS2-fragmentation check.

## Cross-check against Phase 2's existing color test: no rescued signal

Both standout candidates' dedup groups were already inside Phase 2's
tested set (`is_group_representative=True`) and show **no meaningful a\*
correlation** (`color_metabolome_association_a.csv`):

| dedup_group_id | fraction | Spearman rho | empirical p | empirical FDR |
|---|---|---|---|---|
| asid_462 (phytofluene candidate) | cell | 0.12 | 0.08 | 0.48 |
| asid_462 | supernatant | 0.03 | 0.68 | 0.94 |
| asid_21315 (torularhodin candidate) | cell | 0.08 | 0.27 | 0.65 |
| asid_21315 | supernatant | 0.12 | 0.06 | 0.80 |

**This re-mining pass did not recover a hidden color signal** — even if
these features are real carotenoid-pathway intermediates, they don't
correlate with a\* in this dataset either. Consistent with the rest of
this session's results (5 statistical methods, all null).

## Honest caveats on the whole approach

- **Weak evidence at 20 ppm without independent verification.** Most
  matches are low-intensity (single-digit to low-tens total_scans) with
  ppm errors spread across the full ±20 ppm window — some of these are
  plausibly coincidental mass matches in a large search space (36
  targets × ~53,000 features), not confirmed carotenoid detections.
  Exact mass alone is not structural confirmation, especially for
  isomeric/isobaric carotenoid intermediates.
- **No MS2 fragmentation check performed yet** — per the framework's
  Phase 1 Step 4, a plausible mass match should have its MS2 spectrum
  manually inspected for polyene-characteristic neutral losses before
  being trusted. Not done in this pass (would need to pull individual
  scan spectra from the raw mzML/MGF files, not just the aligned feature
  table used here).
- **Zero matches for neurosporene/lycopene/gamma/beta-carotene** could
  mean those specific intermediates genuinely aren't present/detectable,
  or could reflect the coarse 20 ppm / adduct-list limits of this quick
  pass — a systematic RT-based co-elution check (expecting the whole
  desaturation series to elute in a tight RT window under one
  chromatographic method) was not performed; the RTs observed here (5.0-7.6
  min) don't show an obviously clean monotonic desaturation-order pattern,
  which itself is a mild caution flag against over-interpreting these as a
  real biosynthetic series.

## Recommendation

Given the color-correlation cross-check came back null for both
standout candidates, and given the wet-lab escalation (IDEA1_CHEMIST_FRAMEWORK.md
Phase 2) is the natural next step regardless of this pass's outcome per
the decision gate already documented there — **this Phase 1 pass is
inconclusive but not encouraging**: it found plausible mass matches
(strengthening the case that ESI *can* pick up some signal in this mass
range) but none of them explain color variation. Before investing in
Phase 2's wet-lab APCI/APPI escalation, a cheap next step would be pulling
the raw MS2 spectra for row 1735 and row 21315 specifically and checking
fragmentation plausibility (polyene neutral-loss ladder) — if those don't
look like carotenoids on inspection, that's a fast way to rule out this
lane without further instrument time.

## Reproduce

```
python3 analysis/scripts/idea1_targeted_mass_remining.py --ppm 20
```
