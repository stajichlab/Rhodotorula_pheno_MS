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

## What is known about these specific compounds (literature, 2026-09-16)

Background research (not specific to this dataset) on Palmitoyl arginine,
N-myristoyl-arginine, and N6-Palmitoyl lysine as natural products:

- **No paper documents these exact compounds as characterized natural
  products from a named organism.** No fungal biosynthetic route is
  documented for any of them.
- **The closest well-characterized relatives are ornithine lipids (OL)
  and lysine lipids (LL)** — a genuine, established bacterial "aminolipid"
  class (OlsB→OlsA biosynthetic pathway: N-acylation of ornithine or
  lysine, then a second acylation). One review states these lipids are
  **"phosphorus-free and found exclusively in bacteria."** Their
  documented biological role is a **phosphate-starvation response**:
  bacteria substitute these phosphorus-free lipids for phospholipids
  when phosphate is scarce, and separately they have been linked to
  cationic-antimicrobial-peptide/colistin resistance. This is a
  **membrane-lipid role, not a signaling role.**
- **A separate, unrelated documented role exists for N-acyl glycines**
  specifically (e.g. "commendamide," N-3-hydroxypalmitoyl-glycine): human
  gut bacteria produce these as **host-GPCR-mimicking metabolites**
  (agonists of human GPR132/GPR119) — a form of host-microbiome chemical
  communication, but not a bacterium-to-bacterium density-dependent
  quorum signal (Cohen, Fischbach et al., *Nature* 2017, "Commensal
  bacteria make GPCR ligands that mimic human signalling molecules").
- **No literature source directly supports a quorum-sensing role for
  N-acyl-arginine or N-acyl-lysine specifically.** The two documented
  roles found (phosphate-stress membrane lipid; host-GPCR ligand) are
  both real and published, but neither is "quorum sensing" in the sense
  this investigation set out to test.
- **Bottom line on origin**: as with AHLs, the literature-documented
  biology of this compound class is bacterial, with no fungal precedent.
  Compound identity alone cannot distinguish genuine fungal biosynthesis
  from bacterial contamination — but nothing found makes fungal
  production more likely than the contamination explanation already
  reached for AHLs.

## Compartment (cell vs. supernatant) and species-variation follow-up (2026-09-16)

`analysis/scripts/nacyl_amino_acid_compartment_species_analysis.py` pulled
raw per-sample peak areas for the 3 strong candidates (rows 4109, 51126,
51152) plus the 2 supplementary histidine rows, mapped to strain_code,
species, and fraction (cell/supernatant) via the Cu_AUC crosswalk (ID
lookup only).

**All 5 candidates are overwhelmingly cell-associated, not
supernatant-associated.** Paired Wilcoxon signed-rank tests (log10 peak
area + 1, paired by strain, n=265 strains with both fractions) are
significant for every candidate (p<0.0001), all in the same direction:

| row | detected in cell (of 265) | detected in supernatant (of 266) |
|---|---|---|
| 4109 (Palmitoyl arginine) | 252 | 10 |
| 51126 (N-myristoyl-arginine) | 164 | 2 |
| 51152 (N6-Palmitoyl lysine) | 155 | 1 |
| 24998 (histidine, supplementary) | 139 | 7 |
| 40740 (histidine, supplementary) | 94 | 3 |

**This is directly relevant to the quorum-sensing question, and it argues
against it.** A quorum signal has to leave the cell to be sensed by
neighbors — it should be detectable extracellularly (supernatant), not
almost exclusively intracellular/cell-pellet-associated. This compartment
pattern is, however, exactly what the aminolipid (ornithine-lipid/
lysine-lipid) literature above predicts: a phosphorus-free **membrane**
lipid stays with the cell fraction. The MS compartment data and the
literature-documented biology of the closest-known relative compound
class point the same direction, independently of each other.

**Production does vary across species, in the cell fraction.**
Kruskal-Wallis tests across species (n>=3 strains/species, descriptive —
NOT phylogenetically corrected) are significant for all 3 strong
candidates (p=0.001–0.009) and both supplementary rows (p<0.0001 and
p=0.049). *R. taiwanensis* and *R. paludigena* consistently rank highest;
*R. mucilaginosa* (this panel's dominant species by strain count) sits
mid-to-low; some species show near-total absence for some compounds —
e.g. row 51152 (N6-Palmitoyl lysine) is essentially undetected in *R.
toruloides* and *R. sp. clade I* (median log-peak = 0) while present in
most other species tested. Supernatant-fraction species comparisons are
not meaningful for any candidate — almost every strain's supernatant
value is 0, so there is essentially nothing to compare there.

Full per-strain/species/fraction data:
`nacyl_amino_acid_compartment_species.csv`; full test output:
`nacyl_amino_acid_compartment_species_diagnostics.txt`.

## Caveats — what this does and does not show

1. **This is not decoy-null-controlled.** Same open item as the AHL
   search. Until a permutation/decoy control is run, the false-positive
   rate of a 945-target, 20 ppm search in this feature table is unknown.
2. **N-acyl amino acids, unlike AHLs, are chemically unremarkable for a
   fungus to make** in the sense that fatty-acid amide conjugation is
   ordinary lipid chemistry fungi are known to do (e.g. ceramide/
   sphingolipid biosynthesis) — but no fungal precedent exists for THIS
   specific compound family (see literature section above). Detecting a
   member of this class is less surprising on priors than detecting an
   AHL, but it is not evidence of a signaling function, and the
   compartment result above (below) argues specifically against a
   signaling role.
3. **Fungal-vs-bacterial origin is not resolved by compound identity
   alone.** The literature-documented biology of the closest-known
   relative class (ornithine/lysine aminolipids) is described as
   exclusive to bacteria. Combined with the cell-restricted compartment
   pattern, this is more consistent with either (a) a bacterial
   membrane-lipid contaminant, or (b) a structurally analogous but
   previously undocumented fungal membrane lipid, than with a fungal
   signaling molecule.
4. **Species-level variation is descriptive, not phylogenetically
   corrected.** The Kruskal-Wallis tests above do not account for shared
   ancestry between species; a real phylogenetic signal or a lineage
   effect could look like "species variation" without being trait-
   specific. This project's established block-permutation framework
   (used elsewhere, e.g. `phase2_metabolome_phenotype`) has not been
   applied here.
5. **MS2 fragment confirmation has not been attempted** on any row here.

## Bottom line

This remains a more promising exploratory lead than the AHL search
produced — real, chemically coherent SIRIUS structure calls for
N-acyl-arginine and N-acyl-lysine at long acyl chain lengths (C14, C16),
not just formula coincidences, and real, statistically clear variation
across strains and species in the cell fraction. **But the new
compartment result points away from a quorum-sensing interpretation
specifically**: these compounds are almost entirely cell-associated, not
released into the supernatant, which is what a membrane lipid looks like,
not what a diffusible signal looks like — and that is exactly the role
documented in the literature for the closest known relative compound
class (bacterial phosphate-stress aminolipids). The next concrete steps,
if pursued, are: (1) the same decoy/permutation null owed to the AHL
search, (2) MS2 spectral inspection of rows 4109/51126/51152 to confirm
structure, and (3) if a signaling role is still of interest, checking
whether ANY signal is detectable in the supernatant fraction specifically
(the 1-10 strains where it was) rather than the cell fraction, since that
is where a real quorum signal would have to be.

## Files
- `analysis/scripts/nacyl_amino_acid_mass_remining.py` — target-list builder + search
- `nacyl_amino_acid_target_list.csv` — 945 target m/z values (21 amino acids x 15 chain lengths x 3 adducts)
- `nacyl_amino_acid_mass_matches.csv` — 918 raw, unfiltered feature matches at 20 ppm (541 distinct rows)
- `analysis/scripts/nacyl_amino_acid_compartment_species_analysis.py` — cell-vs-supernatant + species-variation follow-up
- `nacyl_amino_acid_compartment_species.csv` / `_diagnostics.txt` — full per-strain/species/fraction data and test output
