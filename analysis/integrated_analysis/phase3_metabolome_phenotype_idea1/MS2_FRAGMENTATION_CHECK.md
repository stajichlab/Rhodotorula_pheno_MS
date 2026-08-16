# Idea 1, Phase 1 follow-up: MS2 fragmentation check for the two standout candidates (2026-08-15)

Per `RESULTS.md`'s recommendation, pulled the raw MS2 spectra for the two
standout mass matches from `aligned_features_filled.mgf`
(`FEATURE_ID` = row ID) and inspected fragmentation plausibility by hand
— not yet a formal reference-spectrum comparison (no torularhodin/
phytofluene reference spectra on hand in this repo), but enough to assess
whether these look like real, informative MS2 spectra of anything, versus
noise-dominated precursor survival.

## Row 1735 (phytofluene [M+NH4]+ candidate) — WEAKENS on inspection

```
PEPMASS=560.5187  RT=336.71s (5.61 min)  Source: SUP_106.mzML
Fragments (m/z, intensity):
  53.02 (2954)  62.28 (33352)  68.04 (2945)  73.46 (2607)  88.11 (13671)
  105.14 (3347)  117.22 (3170)  346.21 (2894)  371.54 (3210)
  560.29 (44545)  560.52 (671657)  <- base peak, essentially the precursor itself
```

- The base peak is the (largely unfragmented) precursor. Every other
  fragment is low-intensity (all ≤44,545 vs. base peak 671,657 — under
  7% relative intensity, several under 1%).
- **No peak near ~543.5 m/z** — the expected [M+NH4]+ → [M+H]+ ammonia
  neutral-loss (Δ17.03), which is the standard diagnostic for a real
  ammonium-adduct CID spectrum. Its absence is a real negative signal,
  not just "no evidence either way."
- The scattered low-mass fragments (53-117) don't form a recognizable
  polyene neutral-loss ladder (no clean, intensity-consistent series of
  losses).
- **Conclusion: this spectrum does not support the phytofluene
  identification.** The earlier "220 total scans" that made this look
  like the most credible candidate reflects precursor-ion intensity/
  stability, not confirmatory fragmentation — a useful reminder that
  scan count / intensity is not the same signal as structural evidence.

## Row 21315 (torularhodin [M+H]+ candidate) — more interesting, genuinely ambiguous

```
PEPMASS=565.4054  RT=342.76s (5.71 min)  Source: C_64.mzML
Fragments (m/z, intensity):
  52.53 (3299)  58.49 (2967)  62.82 (26475)  65.58 (3484)
  211.09 (4196)  214.54 (3263)
  281.29 (28171)  282.28 (118042)  299.30 (26728)
  309.33 (27951)  310.31 (132619)  310.36 (4756)  327.34 (24755)
  565.23-565.58 (precursor region, base peak 565.41 = 559530)
```

- **Two clean water-loss pairs**: 299.30 → 281.29 (Δ18.01, H2O) and
  327.34 → 309.31 (Δ18.01, H2O). Water loss from a carboxylic acid or
  hydroxyl group under CID is a textbook fragmentation, and torularhodin
  specifically has a terminal -COOH (it's the carboxylic-acid form of
  torulene) — this is at least *consistent* with, not contradicting, the
  proposed identity.
- 282.28 and 310.31 are both substantial, well-formed peaks (118,042 and
  132,619 — the two largest fragment peaks in the spectrum, ~20-24% of
  base peak intensity, well above noise level) roughly half the
  precursor mass (565.41 - 282.28 ≈ 283; 565.41 - 310.31 ≈ 255) —
  plausible for a mid-backbone cleavage of a C40 polyene chain, though
  this is a coarse plausibility check, not a confirmed fragment
  assignment (would need exact-formula fragment matching against a
  polyene fragmentation model to confirm).
- **Conclusion: this spectrum is not noise, and its fragmentation pattern
  (paired water losses, substantial mid-mass fragments) is chemically
  coherent with an oxygenated, roughly-C40-backbone structure.** This is
  more credible than row 1735's spectrum, despite its much lower total
  intensity/scan count — a useful reminder in the other direction, that
  overall intensity alone doesn't rank spectral quality either.
- **Still not a confirmed identification.** No reference spectrum
  comparison was done (none available in this repo); the water-loss
  pattern is consistent with many oxygenated C40 isoprenoids, not
  uniquely diagnostic of torularhodin specifically over other candidate
  structures at a similar mass (e.g. other carotenoid-related
  acids/alcohols, or an unrelated oxygenated lipid of similar mass that
  happens to fall in this m/z window).

## Net effect on Idea 1's standing recommendation

- **Row 1735 is no longer a credible lead** — demote it. The intensity
  that made it look promising was a property of an undissociated
  precursor, not evidence for the proposed identity.
- **Row 21315 remains a live, if unconfirmed, candidate.** Its
  fragmentation is chemically coherent with an oxygenated carotenoid-like
  structure, it's the single feature SIRIUS most visibly misannotated
  (flagging it as "Polyamines"), and the mass match itself was tight
  (+2.5 ppm). Still, **it does not correlate with color** (a\* Spearman
  rho 0.08/0.12, FDR 0.65/0.80 — see `RESULTS.md`), so even a confirmed
  identification here would not, on its own, explain the color-metabolome
  null.
- **Recommendation unchanged from RESULTS.md**: before any wet-lab
  escalation (`IDEA1_CHEMIST_FRAMEWORK.md` Phase 2), a genuine reference-
  spectrum comparison for row 21315 (against a public spectral library
  entry for torularhodin/torulene-family compounds, if one exists, e.g.
  via GNPS library search) would be the next cheap, no-new-data step —
  this hand inspection is suggestive, not confirmatory.

## Reproduce

Spectra extracted by `FEATURE_ID` from
`data/processed/EB_20260130_ExFAB_Rhodo_Sup_and_Pellet/b773ffa18c2b41e5a3484526293a54f9/b773ffa18c2b41e5a3484526293a54f9/nf_output/feature_finding/aligned_features_filled.mgf`
(row IDs 1735 and 21315). No script was written for this one-off
extraction (a simple linear MGF `BEGIN IONS`/`END IONS` block scan by
`FEATURE_ID`); worth turning into a small reusable
`analysis/scripts/extract_mgf_spectrum.py --feature-id N` utility if this
kind of spot-check becomes routine.
