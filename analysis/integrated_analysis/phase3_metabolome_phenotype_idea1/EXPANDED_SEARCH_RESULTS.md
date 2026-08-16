# Idea 1, Phase 1 (expanded): Broader target list + phenotype cross-check (2026-08-16)

Follow-up to `RESULTS.md`/`MS2_FRAGMENTATION_CHECK.md`, per PI request to
expand the mass search beyond the core carotenogenesis pathway and check
candidates against phenotype more broadly (not just color).

## Expanded target list

`analysis/scripts/idea1_targeted_mass_remining.py` now searches 19
compounds (76 targets with 4 adducts) across 4 categories — see the
script's module docstring for full formula/rationale citations:

| Category | Compounds | Why added |
|---|---|---|
| carotenoid_pathway (original) | phytoene → torularhodin (9) | Core biosynthetic pathway |
| apocarotenoid | β-ionone, β-cyclocitral, dihydroactinidiolide | Smaller, more polar cleavage fragments — more ESI-friendly than intact C40 backbones; SIRIUS already called one existing feature "Apocarotenoids" |
| additional_carotenoid | astaxanthin, canthaxanthin, echinenone, zeaxanthin | Oxygenated variants from related carotenogenic yeasts — **not strain-confirmed for *Rhodotorula***, broader net |
| sterol_pathway | ergosterol, squalene, lanosterol | Shares the same upstream isoprenoid/MVA precursor pool as carotenoids — tests a precursor-competition hypothesis, not a pigment search |

31 raw-feature matches found (20 ppm), up from 11 in the original pass.

## Two new standout candidates

1. **Ergosterol [M+H-H2O]+, row 846** (-4.3 ppm, 71 total scans, **flagged
   by EB's own pipeline as an in-source fragment**). This is the most
   chemically self-consistent match found in either search pass: sterol
   [M+H-H2O]+ (loss of the 3-hydroxyl) is the textbook dominant ion for
   this compound class, and the untargeted pipeline's own independent
   ISF-detection algorithm flagged exactly that relationship — two
   separate pieces of evidence agreeing, not one mass coincidence.
2. **Astaxanthin [M+H]+, row 9384** (+2.4 ppm, 37 scans). Tight mass
   match; astaxanthin is not classically reported in *Rhodotorula*
   specifically (more an *Xanthophyllomyces*/*Phaffia* pigment), so this
   should be read as "plausible mass, uncertain biology" pending a
   literature check, not an expected hit.

Also of note: the earlier torularhodin candidate (row 21315) is
**isobaric with canthaxanthin** (both C40H52O2) — mass alone cannot
distinguish them; `MS2_FRAGMENTATION_CHECK.md`'s water-loss pattern
(consistent with a hydroxyl/carboxyl loss) is arguably more consistent
with torularhodin's -COOH than canthaxanthin's diketone structure, but
this is not a confident structural call either way.

## Phenotype cross-check: color (a\*) — null for all three

| Candidate | Fraction | a\* Spearman rho | empirical FDR (from Phase 2's already-run pipeline) |
|---|---|---|---|
| Ergosterol (asid_31011) | cell | 0.02 | 0.99 |
| Ergosterol (asid_31011) | supernatant | -0.06 | 0.85 |
| Astaxanthin (asid_9384) | cell | -0.06 | 0.73 |
| Astaxanthin (asid_9384) | supernatant | -0.05 | 0.89 |
| Torularhodin/canthaxanthin (asid_21315) | cell | 0.08 | 0.65 |
| Torularhodin/canthaxanthin (asid_21315) | supernatant | 0.12 | 0.80 |

No color signal for any of the new candidates either — same story as
before.

## Phenotype cross-check: copper-resistance AUC — one real lead, NOT yet validated

Quick TSS-normalized naive Spearman correlation (mean_auc_rate from
`analysis/linked_data/sample_metadata.csv.gz`, **not** the full
phylogenetically-corrected pipeline — see caveat below):

| Candidate | Fraction | n | AUC Spearman rho | p (naive, uncorrected) |
|---|---|---|---|---|
| **Ergosterol (asid_31011)** | **cell** | **267** | **0.226** | **0.0002** |
| Ergosterol (asid_31011) | supernatant | 268 | -0.072 | 0.24 |
| Astaxanthin (asid_9384) | cell | 267 | constant/undetected | — |
| Astaxanthin (asid_9384) | supernatant | 268 | 0.006 | 0.92 |
| Torularhodin/canthaxanthin (asid_21315) | cell | 267 | 0.064 | 0.30 |
| Torularhodin/canthaxanthin (asid_21315) | supernatant | 268 | -0.018 | 0.77 |

**Ergosterol's cell-fraction abundance correlates with copper-resistance
growth rate (rho=0.23, p=0.0002 naive, n=267)** — a real, nominally
significant, and mechanistically plausible signal (ergosterol is the
principal fungal plasma-membrane sterol; membrane sterol content is a
well-established determinant of metal-ion tolerance in fungi generally).

**This is explicitly NOT a validated finding yet.** Per this project's
now-established pattern (`.living/findings/phylogenetic-confounding-of-trait-molecular-associations.md`,
F-005: 2-for-2 track record of naive whole-panel correlations failing to
survive phylogenetic block-permutation or within-species restriction —
including the *exact same copper-AUC phenotype*, whose naive amino-acid
"hits" did not hold up in the *R. mucilaginosa*-only sensitivity check),
this naive rho/p should be treated as a lead requiring the same rigor
before being trusted, not a result. It was computed as a quick screen
(single Spearman correlation, no phylogenetic blocking, no negative
control, no BH-FDR context) specifically to decide whether investing in
the full pipeline is worthwhile — and on that basis, **yes, this looks
worth the full test**.

## Recommendation / next step (not yet done)

Run ergosterol (asid_31011, cell fraction) through the same rigor as
everything else in this project before trusting it: phylogenetic
block-permutation empirical p-value (species-level primary,
*R. mucilaginosa*-only as secondary), with the same hard-gated
negative-control convention (e.g. colony area as decoy) already
established in `phase2_color_metabolome_association.py`/
`phase2_within_species_association.py`. This would need a small new
script variant (predictor = mean_auc_rate from sample_metadata rather
than the color phenotype table) — not yet written.

## Reproduce

```
python3 analysis/scripts/idea1_targeted_mass_remining.py --ppm 20
```
The AUC quick-check above was a one-off inline script, not yet saved as
a reusable file — worth promoting to
`analysis/scripts/idea1_auc_quickcheck.py` if this lane gets pursued
further.
