# N-acyl amino acid search

## Why this search exists

Expansion of the AHL autoinducer search (2026-09-16, PI request) to a
chemically related molecule class: N-acyl amino acids (a fatty acyl chain
condensed onto an amino acid's free amine). These are made by the same
acyl-CoA/acyl-ACP donor chemistry as AHLs, and some LuxI-family
acyltransferases are documented to make N-acyl amino acids instead of (or
alongside) AHLs. Unlike AHLs, N-acyl amino acids are NOT an
exclusively-bacterial signal class — they are a broadly distributed
natural product family (documented in bacteria as "lipoamino acids"/
biosurfactants, and in other contexts including host-microbiome lipid
signaling). A quorum-sensing role is documented for specific bacterial
systems; it is not established as a general property of the class, and
not established in fungi at all.

**This search answers "is this chemical class present," not "does it
function as a quorum signal here."** Those are different questions —
see the caveats section below.

## Method

`analysis/scripts/nacyl_amino_acid_mass_remining.py` — same design as the
AHL search: exact-mass re-mining of the raw EB feature table
(`aligned_features_ms2.csv`), independent of SIRIUS's own compound calls,
at 20 ppm tolerance.

Target list, built systematically rather than cherry-picked: all 20
standard proteinogenic amino acids, plus homoserine (included because
this project already established *R. mucilaginosa* has lactonase activity
— see `AHL_AUTOINDUCER_SEARCH.md`'s literature-check section — so the
open-chain hydrolysis product of any bacterial AHL contaminant, N-acyl-
homoserine, is a directly motivated additional target). Each amino acid
crossed with saturated acyl chain length Cn, n=4-18 (same range as the AHL
search), crossed with 3 adducts ([M+H]+, [M+Na]+, [M+NH4]+). 21 amino
acids x 15 chain lengths x 3 adducts = 945 targets.

Formula derivation (amide condensation, losing H2O):
`N-acyl-AA = C(n+a) H(2n+h-2) N(1) O(o+1)` for a saturated Cn fatty acid
and amino acid C(a)H(h)N(1)O(o). Verified against a known reference
compound: N-lauroylglycine (n=12, glycine) predicts C14H27NO3, matching
its published formula.

**Not searched**: negative-mode ions. N-acyl amino acids carry a free
carboxylic acid group and would be expected to ionize well as [M-H]- in
ESI negative mode. Checked this project's raw feature table directly
(the `adduct` column) — it contains **only positive-mode adducts**.
Negative-mode LC-MS data does not exist for this sample set. This is a
real, unclosable gap with the current data, not an oversight. Also not
searched: unsaturated/hydroxylated acyl chains, non-standard amino acids
beyond homoserine.

## Result

918 raw feature matches across 541 distinct feature rows (a row can match
more than one amino-acid/chain-length/adduct combination when those
combinations happen to be isobaric). This raw count is not evidence of
anything by itself — like the AHL search, no decoy/permutation null has
been run, and this formula space (CxHyNOz, a mix of common amino-acid-
and lipid-like elemental ratios) is if anything less distinctive than
the AHL search's, so an even higher false-positive rate should be
expected here without a null control.

## Structural cross-check: stronger than the AHL search, still not proof of function

228 of the 541 matched rows have an independent SIRIUS structural call.
Filtering to rows where **SIRIUS's own assigned formula exactly matches
the target formula** (the same discriminator that ruled out the AHL
search's one candidate) narrows this to 25 rows. Most of these 25 are
chemically incompatible on inspection (dipeptides, tripeptides — wrong
structure class despite matching formula, the same coincidental-overlap
pattern seen throughout the AHL search). But a genuinely different result
appears in the long-chain, arginine/lysine-conjugate subset:

| row_id | target compound | SIRIUS class | SIRIUS structure name | scans |
|---|---|---|---|---|
| 4109 | N-C16:0-arginine | N-acyl amines | **Palmitoyl arginine** | 26 |
| 51126 | N-C14:0-arginine | N-acyl amines | **N2-(1-Oxotetradecyl)-L-arginine** (= N-myristoyl-arginine) | 9 |
| 51152 | N-C16:0-lysine | N-acyl amines | **N6-Palmitoyl lysine** | 6 |
| 24998 / 40740 | N-C16:0-histidine | N-acyl amines | (class match, no specific structure name) | 7 / 6 |

For rows 4109, 51126, and 51152, SIRIUS's own independent structure
assignment — not just the elemental formula — names the exact compound
the target list was built to find. This is a materially stronger
correspondence than anything found in the AHL search, where the one
formula-adjacent candidate turned out to be a mismatched, unrelated
compound.

**A second candidate cluster was flagged and then set aside as likely
non-biological**: rows 3254/13159 (N-C13:0-lysine) and 6227/41488
(N-C15:0-lysine/N-C16:0-proline) also have exact SIRIUS-formula matches
and an "N-acyl amines" class call, but the specific structures SIRIUS
names are betaine-type quaternary-ammonium lipids (e.g.
"3-(Tetradecanoylamino)propyl(carboxymethyl)dimethylammonium"), one of
them explicitly a **deuterium-labeled** compound name
("...tricosadeuteriododecanoylamino..."). A deuterated compound would
not share the same nominal mass as its non-deuterated formula at this
precision, so this specific match most likely reflects a spectral-
library reference entry (e.g. an isotope-labeled internal standard used
in this or another LC-MS run) rather than the true identity of the
detected ion. Flagged, not treated as supporting evidence.

## Caveats — what this does and does not show

1. **This is not decoy-null-controlled.** Same open item as the AHL
   search. Until a permutation/decoy control is run, the false-positive
   rate of a 945-target, 20 ppm search in this feature table is unknown.
2. **N-acyl amino acids, unlike AHLs, are chemically unremarkable for a
   fungus to make.** Fatty-acid amide conjugation is ordinary lipid
   chemistry (related to ceramide/sphingolipid biosynthesis and other
   acyltransferase activity fungi are known to have). Detecting one is
   therefore much less surprising, on priors, than detecting an AHL — but
   it is also much weaker evidence of a *signaling* function. The
   detection of Palmitoyl-arginine-like compounds says "this chemical
   class exists in the extract," not "this strain uses it to signal."
3. **Fungal-vs-bacterial origin is not resolved.** N-acyl-arginine and
   N-acyl-lysine are documented bacterial lipoamino acids (membrane/
   stress-response lipids in several bacterial genera). The same trace-
   contamination caveat that applies to the AHL search applies here,
   though less severely, since fungi/other eukaryotes are not excluded
   from making this compound class the way they effectively are for
   AHLs (no known LuxI-type synthase in any fungal genome, see
   `AHL_AUTOINDUCER_SEARCH.md`).
4. **No strain-level detection-vs-morphology analysis has been run yet**
   for these candidates. The AHL search's Phase 2 (detection vs.
   colony-texture proxy) has not been repeated here. If this line is
   pursued further, that is the natural next step for rows 4109, 51126,
   and 51152 specifically.
5. **MS2 fragment confirmation has not been attempted** on any row here
   either.

## Bottom line

This is a more promising exploratory lead than the AHL search produced —
real, chemically coherent SIRIUS structure calls for N-acyl-arginine and
N-acyl-lysine at long acyl chain lengths (C14, C16), not just formula
coincidences. It is not evidence of quorum sensing, and not yet evidence
of anything strain- or morphology-specific. The next concrete steps, if
pursued, are: (1) the same decoy/permutation null owed to the AHL search,
(2) strain-level detection of rows 4109/51126/51152 against the texture
proxy (the Phase 2 pattern already built for AHLs), and (3) MS2 spectral
inspection of these three rows specifically, since they are the first
candidates in this whole investigation to clear a basic structural
plausibility bar.

## Files
- `analysis/scripts/nacyl_amino_acid_mass_remining.py` — target-list builder + search
- `nacyl_amino_acid_target_list.csv` — 945 target m/z values (21 amino acids x 15 chain lengths x 3 adducts)
- `nacyl_amino_acid_mass_matches.csv` — 918 raw, unfiltered feature matches at 20 ppm (541 distinct rows)
