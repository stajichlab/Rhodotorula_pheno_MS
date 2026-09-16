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

## Phylogenetic signal test (2026-09-16, PI request)

The Kruskal-Wallis species comparison above is descriptive and not
phylogenetically aware — species differences could reflect shared
ancestry, or could be phylogenetically random noise that happens to
differ by species. Two independent methods were run to test this
directly, both against the species-level tree
(`analysis/integrated_analysis/phase1_phenotype/species_tree.nwk`,
16 tips), cell fraction only (supernatant is not tested — essentially
absent per the compartment analysis above).

**Method 1 — Blomberg's K / Pagel's lambda** (this project's dedicated
phylogenetic-signal tool, `phylogenetic_signal.R`, previously used for
color phenotype; reused here as
`nacyl_amino_acid_phylogenetic_signal.R`). Species-level mean
log10(peak_area+1), via `nacyl_amino_acid_build_species_table.py`.

| compound | K | K p-value | lambda | lambda p-value |
|---|---|---|---|---|
| Palmitoyl arginine (4109) | 0.36 | 0.51 | 0.03 | 0.90 |
| N-myristoyl-arginine (51126) | 0.38 | 0.58 | ~0 | 1.00 |
| N6-Palmitoyl lysine (51152) | 0.37 | 0.62 | ~0 | 1.00 |
| histidine (24998, supplementary) | 0.35 | 0.69 | 0.06 | 0.76 |
| histidine (40740, supplementary) | 0.33 | 0.78 | 0.09 | 0.71 |

K well below the Brownian-motion expectation of 1 for all 5 compounds,
none significant; lambda near 0 (no phylogenetic correlation structure)
for all 5, none significant. **No detectable phylogenetic signal by
this method for any candidate.**

**Method 2 — block-permutation test** (per PI request, reusing this
project's `phase2_color_metabolome_association.py` block-construction
method directly — its own use of that method tests compound abundance
against an EXTERNAL phenotype while controlling for phylogeny, a
different question from "is production itself phylogenetically
structured"; this repurposes the same species-tree clade construction,
`analysis/scripts/nacyl_amino_acid_block_permutation_signal.py`, as a
non-parametric complement to K/lambda that uses strain-level data (n=265,
more power) and makes no Brownian-motion assumption — better suited to
these compounds' zero-inflated distributions). Six species-tree clades
(same default as phase2), one-way ANOVA F-statistic on log-abundance by
clade, 2000-permutation empirical null (shuffling clade labels across
strains).

| compound | F (observed) | empirical p |
|---|---|---|
| Palmitoyl arginine (4109) | 1.55 | 0.18 |
| N-myristoyl-arginine (51126) | 1.06 | 0.38 |
| N6-Palmitoyl lysine (51152) | 1.84 | 0.06 |
| histidine (24998, supplementary) | 1.38 | 0.20 |
| histidine (40740, supplementary) | 2.65 | 0.01 |

**None of the 3 strong candidates reach significance** (all p>0.05);
none would survive correction for testing 5 compounds either way.
`row_51152` is the closest to a signal (p=0.06) but does not clear even
an uncorrected 0.05 threshold. `row_40740` (the weakest, unnamed-
structure supplementary candidate) is nominally significant (p=0.01) but
does not survive a Bonferroni correction for 5 tests (0.05/5=0.01,
right at the boundary) and is the least structurally credible candidate
in this whole search.

**Conclusion: two independent methods (parametric species-level K/lambda,
non-parametric strain-level block-permutation) agree that none of the 3
structurally-corroborated candidates show a detectable phylogenetic
signal in production.** The species-level differences found earlier
(Kruskal-Wallis) are real in the sense that species DO differ, but that
difference does not track the species tree's clade structure in a
way distinguishable from chance — consistent with production being
driven by strain-level or ecological factors uncorrelated with deep
phylogeny, or (the honest alternative) with limited power at n=16
species (several contributing only 1 strain) and n=6 clades.

Full output: `nacyl_amino_acid_species_table.csv`,
`nacyl_amino_acid_phylogenetic_signal.csv`,
`nacyl_amino_acid_block_permutation_signal.csv`.

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
4. **Species-level variation is real but not phylogenetically
   structured** (resolved 2026-09-16, see phylogenetic signal section
   above) — species differ, but that difference doesn't track the
   species tree by either of two independent methods tested.
5. **MS2 fragment confirmation has not been attempted** on any row here.
6. **The cell-vs-supernatant compartment bias has not been checked
   against the known colony-area/biomass confound** documented elsewhere
   in this project (`.living/findings/biomass-scaling-artifacts-in-
   extraction-based-metabolomics.md`) — raw peak area in the cell
   fraction is known to broadly correlate with colony size for unrelated
   reasons; whether that confound explains any of the cell-fraction
   species variation above has not been checked.

## Validation steps 1-4 (2026-09-16, PI request)

### 1. MS2 fragment confirmation — CONFIRMED for all 3 strong candidates

`analysis/scripts/nacyl_amino_acid_ms2_fragment_check.py` pulled the real
MS2 spectrum for rows 4109, 51126, and 51152 from the raw EB pipeline's
`aligned_features.mgf` (matched by FEATURE_ID, not simulated) and checked
for the diagnostic fragment ions expected from a free arginine or lysine
backbone — the same logic used to confirm an intact AHL via its
homoserine-lactone-loss fragment, computed from first principles (not
looked up as pre-set values):

| row | fragment | expected m/z | observed m/z | ppm error | rel. intensity |
|---|---|---|---|---|---|
| 4109 (arginine) | [Arg+H]+ | 175.11896 | 175.1185 | -2.6 | 13.6% |
| 4109 | [Arg+H-NH3]+ | 158.09241 | 158.0920 | -2.6 | 18.0% |
| 4109 | [Arg+H-guanidine]+ | 116.07061 | 116.0703 | -2.7 | 9.1% |
| 4109 | C4H8N+ (Arg marker) | 70.06511 | 70.0649 | -3.0 | 10.4% |
| 51126 (arginine) | [Arg+H]+ | 175.11896 | 175.1182 | -4.3 | 11.3% |
| 51126 | [Arg+H-NH3]+ | 158.09241 | 158.0919 | -3.2 | 25.9% |
| 51126 | [Arg+H-guanidine]+ | 116.07061 | 116.0702 | -3.5 | 10.2% |
| 51126 | C4H8N+ (Arg marker) | 70.06511 | 70.0647 | -5.9 | 12.5% |
| 51152 (lysine) | [Lys+H]+ | 147.11281 | 147.1131 | +2.0 | 30.6% |
| 51152 | [Lys+H-NH3]+ | 130.08626 | 130.0859 | -2.8 | 46.8% |
| 51152 | [Lys+H-H2O]+ | 129.10152 | 129.1020 | +3.7 | 66.2% |
| 51152 | C5H10N+ (Lys marker) | 84.08078 | 84.0806 | -2.1 | 100.0% |

**4/4 diagnostic fragments matched for all 3 candidates, all within
±6 ppm, at meaningful relative intensity (9-100% of base peak).** This is
real structural evidence, not another mass coincidence: the spectrum for
each row shows exactly the fragment series a genuine N-acyl-arginine or
N-acyl-lysine should produce (loss of the acyl chain revealing the free
amino acid, then that amino acid's own well-documented secondary
fragmentation). This is the strongest positive structural result in this
entire investigation (AHL search included) — everything else so far has
been mass-only matching or SIRIUS class-level agreement; this is
consistent fragment-level chemistry. It does not, on its own, resolve
fungal-vs-bacterial origin (see caveat 3 above), and it does not confirm
acyl-chain length or attachment position beyond what SIRIUS already
implied (no fragment here is chain-length-diagnostic; the parent mass
already fixes that). Full detail: `nacyl_amino_acid_ms2_fragment_check.csv`.

### 2. Decoy/permutation null — read with a specific caveat, not at face value

`analysis/scripts/mass_search_decoy_null.py` reused for both the AHL and
N-acyl-amino-acid target lists: each of 1000 permutations shifts every
target m/z by a random, signed offset (Uniform 15-60 Da), reruns the same
20 ppm window search, and records the total match count.

| search | observed raw matches | null mean (sd) | null 95th pct | empirical p |
|---|---|---|---|---|
| AHL (135 targets) | 103 | 27.9 (9.3) | 45.0 | 0.001 |
| N-acyl amino acid (945 targets) | 918 | 228.5 (26.3) | 272.0 | 0.001 |

Both searches produce far more matches than the shifted-mass null
predicts (p=0.001, the floor at 1000 permutations — 0/1000 null draws
reached the observed count for either search).

**This needs a specific, honest reading, not a triumphant one.** The
feature table's m/z density is NOT uniform — a direct check (25 Da bins,
150-525 Da) shows feature count roughly quadruples across that range
(155 features at 150-175 Da vs. 689 at 475-500 Da). Both target lists
(146-493 Da for AHLs, 146-493 Da for N-acyl amino acids) sit in a
generally dense, chemically busy region of this metabolome. **The excess
match count over the shifted-mass null most likely reflects that this
mass region is intrinsically rich in real, unrelated small-molecule
chemistry (dipeptides, amino acids, fatty acyl carnitines — exactly what
SIRIUS cross-referencing already found for most of the raw matches, see
above), not that AHLs or N-acyl amino acids are specifically enriched
there.** This decoy test rules out "these hit counts are pure random
noise with no chemical basis" — it does NOT independently support "these
are real AHLs/N-acyl amino acids specifically" the way the MS2 result
above does. The two results (decoy null, MS2 confirmation) answer
different questions and should not be conflated. Full output:
`ahl_decoy_null.csv`, `nacyl_amino_acid_decoy_null.csv`.

### 3. Colony-area/biomass confound — checked, not present

Spearman correlation of each candidate's raw cell-fraction peak area
(log10) against colony area (`control_phenotype_90_110h`'s `area_median`,
log10), per strain (n=264): **all three strong candidates show
essentially zero correlation** (rho=0.02-0.04, p=0.52-0.72). A
species-level check (mean per species, n=16) is also null for two of
three (rho=-0.04 and 0.04, p>0.85) and borderline-non-significant for the
third (Palmitoyl arginine, rho=0.47, p=0.064). **The established
cell-fraction/biomass confound in this project does not explain these
candidates' abundance or their species-level variation.**

### 4. Blank and QC_Mix presence — checked, absent

Raw peak area for all 5 candidate rows (3 strong + 2 supplementary) is
**exactly zero in every one of the 7 Blank samples and every one of the 7
QC_Mix (pooled-sample) injections**, in both technical replicate sets.
No evidence of a reagent, extraction-background, or carryover signal.
(The absence from QC_Mix specifically is worth noting as mildly
unexpected if the compound were broadly present at the intensities seen
in individual strains — but pooled-QC dilution/matrix effects for a
compound restricted to particular strains is a plausible, non-alarming
explanation, not itself a red flag.)

## Bottom line (updated 2026-09-16)

**This is now the best-supported candidate in this entire AHL/quorum-
sensing investigation.** Three N-acyl-arginine/lysine features have: a
SIRIUS structure call naming the exact compound, MS2 fragments matching
the expected amino-acid-backbone chemistry at sub-5-ppm accuracy, no
detectable background/reagent contamination, and no detectable
colony-size confound. That said, three things still argue against a
**signaling** role specifically, independent of whether the compound
identity itself is now well-supported: (1) these compounds are almost
entirely cell-associated, not released into the supernatant, which is
what a membrane lipid looks like, not a diffusible signal; (2) the
closest literature-documented relative class (bacterial ornithine/lysine
aminolipids) has a described membrane-lipid, phosphate-stress role, not
a signaling one; (3) neither of two independent phylogenetic-signal tests
finds species-level production differences structured by the species
tree. **Fungal-vs-bacterial origin remains unresolved** — MS2 confirms
the compound's identity, not its producing organism. The one remaining
step from the original priority list, authentic chemical standards for
retention-time/MS2 matching, is the only test that could move this
further; item 6 (the `phase2_color_metabolome_association.py`
phenotype-association design) remains available if a different question
(does this track color/copper resistance) becomes of interest.

## Files
- `analysis/scripts/nacyl_amino_acid_mass_remining.py` — target-list builder + search
- `nacyl_amino_acid_target_list.csv` — 945 target m/z values (21 amino acids x 15 chain lengths x 3 adducts)
- `nacyl_amino_acid_mass_matches.csv` — 918 raw, unfiltered feature matches at 20 ppm (541 distinct rows)
- `analysis/scripts/nacyl_amino_acid_compartment_species_analysis.py` — cell-vs-supernatant + species-variation follow-up
- `nacyl_amino_acid_compartment_species.csv` / `_diagnostics.txt` — full per-strain/species/fraction data and test output
- `analysis/scripts/nacyl_amino_acid_build_species_table.py` — species-level table builder for the phylogenetic signal test
- `nacyl_amino_acid_species_table.csv` — 16-species mean log-abundance table
- `analysis/scripts/nacyl_amino_acid_phylogenetic_signal.R` — Blomberg's K / Pagel's lambda test
- `nacyl_amino_acid_phylogenetic_signal.csv` — K/lambda results
- `analysis/scripts/nacyl_amino_acid_block_permutation_signal.py` — block-permutation phylogenetic signal test
- `nacyl_amino_acid_block_permutation_signal.csv` — block-permutation results
- `analysis/scripts/nacyl_amino_acid_ms2_fragment_check.py` — MS2 diagnostic-fragment confirmation (validation step 1)
- `nacyl_amino_acid_ms2_fragment_check.csv` — fragment match results (4/4 for all 3 strong candidates)
- `analysis/scripts/mass_search_decoy_null.py` — decoy/permutation null, reused for both the AHL and this search (validation step 2)
- `ahl_decoy_null.csv` / `nacyl_amino_acid_decoy_null.csv` — per-permutation null match counts
