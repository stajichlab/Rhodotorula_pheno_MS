# AHL and Quorum-Sensing Molecule Search: Comprehensive Report

Date: 2026-09-16
Project: Rhodotorula_pheno_MS
Analysis folder: `analysis/ahl_autoinducer_search/`

## Scope statement

This report covers the MS-data search for **N-acyl-homoserine lactones
(AHLs)** and, as of 2026-09-16, an expansion to **N-acyl amino acids** in
the *Rhodotorula* MS2 metabolomics dataset. It also covers a cross-project
GWAS follow-up and the literature check that shaped the interpretation.

**The MS mass search covers only these two chemical families.** It did
not search for AI-2 (a furanosyl borate diester, structurally unrelated to
both), diffusible signal factor (DSF), or the aromatic-alcohol quorum-
sensing molecules native to fungi (farnesol, tyrosol, phenylethanol,
tryptophol). Those molecule classes are mentioned below only as literature
context, not as things this project searched for in the MS data. If a
search for those classes is wanted, it has not been done yet.

## Bottom line

**AHLs**: no positive evidence was found that *Rhodotorula* strains
produce them. Four independent checks were run. All four are consistent
with "no AHL production," none with "AHL production confirmed."

| Check | Result |
|---|---|
| Exact-mass search of the MS2 feature table | 103 raw matches, but the detection pattern does not discriminate between strains |
| Structural cross-check (SIRIUS) on the best candidate | Argues against AHL identity |
| Cross-project GWAS for a genetic locus | No locus survives scrutiny |
| Literature basis for the hypothesis | Rests on a functional bioreporter assay, not a confirmed chemical structure; the same species is separately documented to degrade AHLs, not make them |

This does not prove *Rhodotorula* cannot make AHLs. It means nothing found
so far supports that it does. The one remaining un-run check is MS2
fragment confirmation (below).

**N-acyl amino acids**: a real, variable compound class, but a follow-up
argues against a signaling role. Three long-chain (C14-C16) N-acyl-
arginine/lysine candidates have an independent SIRIUS structure call that
names the exact target compound (e.g. "Palmitoyl arginine"), not merely a
matching mass — stronger correspondence than anything found for AHLs.
Literature on the closest known relative class (bacterial ornithine/
lysine "aminolipids") describes a membrane-lipid, phosphate-stress role,
not signaling. Consistent with that, all candidates are almost entirely
cell-associated, not present in the supernatant (p<0.0001) — the opposite
compartment pattern from what a diffusible quorum signal would need.
Production does vary significantly by species in the cell fraction. See
section 8.

## 1. Where the hypothesis came from, and why it is weak on its own terms

The PI traced the likely source of the "*Rhodotorula* makes AHLs"
hypothesis to:

> Wilson et al. 2025, "Characterization of virulence-related phenotypes of
> *Candida parapsilosis* and *Rhodotorula mucilaginosa* isolated from the
> International Space Station (ISS)," *Life Sciences in Space Research*
> 45:16–24. doi:10.1016/j.lssr.2025.01.002. PMID 40280638.

The paper's own abstract states: "autoinducer (AI) production was
detected by activation of a reporter fluorescent gene present in
biosensor bacterial strains." This is a **functional bioassay**, not a
chemical measurement. The abstract reports no mass spectrometry and no
NMR. The paper's stated result is "increased ... long-chain autoinducer
production" in ISS-isolated *R. mucilaginosa* versus ATCC controls,
alongside increased capsule production, biofilm formation, antifungal
resistance, and nematode virulence — an overall "enhanced virulence"
narrative, not a targeted AHL-versus-morphology test. The paper's own
morphology-adjacent readouts are capsule production and filamentation,
not a smooth/rough colony score.

Bacterial LuxR-type AHL bioreporters are documented in the
quorum-sensing literature to respond to non-AHL lipophilic compounds,
including some fatty acids. *Rhodotorula* produces abundant lipids and
fatty-acid derivatives. A positive bioreporter signal from a
*Rhodotorula* extract is therefore also consistent with (a) cross-
reactive fungal lipid chemistry, or (b) residual bacterial AHL in a
non-axenic culture. The bioreporter result alone cannot distinguish these
from genuine fungal AHL synthesis. The paper itself was not accessed
(paywalled); this reading is based on the indexed abstract only.

## 2. Is fungal AHL biosynthesis known to exist at all?

A background literature and database check, independent of this
project's own MS data, found:

- **No LuxI-type AHL synthase ortholog exists in any fungal genome per
  KEGG.** KEGG ortholog K18096 (the CoA-utilizing LuxI homolog) is
  restricted to alpha-Proteobacteria. No Pfam/InterPro entry for a
  homoserine-lactone synthase domain family in fungi was found.
- **No natural fungal AHL production has been cloned or characterized.**
  The one closely related paper (*Commun. Biol.* 2025) engineers a
  bacterial LuxI gene into *S. cerevisiae* to make it produce AHLs — an
  explicitly heterologous, engineered system, not a natural yeast
  pathway.
- **Directly relevant counter-finding: *Rhodotorula mucilaginosa* is
  published as an AHL-degrading yeast**, with lactonase activity
  confirmed against C6-HSL, 3-oxo-C6-HSL, and 3-hydroxy-C6-HSL (Tan et
  al. 2014, *Sensors* 14:6463, PMC4029656). AHL-inactivating
  (lactonase/acylase) activity is reported as widespread across yeasts
  generally, and the basidiomycetous yeast *Trichosporon loubieri* is
  shown to consume AHLs as a carbon/nitrogen source (PMC3859043). This is
  the opposite of the hypothesis: the genus is documented to destroy
  AHLs, not make them.
- Fungi have their own, chemically unrelated quorum-sensing chemistry
  (farnesol, tyrosol, phenylethanol, tryptophol — aromatic
  alcohols/sesquiterpenes that regulate yeast-to-hypha switching;
  Wongsuk et al. 2016, *J. Basic Microbiol.*). These are not homoserine
  lactones and were not searched for in this project's MS data.

**Implication**: if an AHL-mass signal appears in the MS data, trace
bacterial contamination (possibly since degraded by the strain's own
lactonase) is a more likely explanation than novel fungal biosynthesis.
A genuinely novel fungal AHL pathway would not yet be in KEGG or Pfam, so
this does not rule out the hypothesis. But there is no known fungal AHL
synthase gene to anchor a genomic search against.

## 3. MS mass search (Phase 1)

**Method**: `analysis/scripts/ahl_targeted_mass_remining.py` searched the
raw EB feature table (`aligned_features_ms2.csv`, 53,040+ features),
independent of SIRIUS's own compound calls. The target list was built
from the general AHL formula, not hand-picked, across acyl chain length
n = C4 to C18:

- unsubstituted Cn-HSL: C(n+4)H(2n+5)NO3
- 3-oxo-Cn-HSL: C(n+4)H(2n+3)NO4
- 3-hydroxy-Cn-HSL: C(n+4)H(2n+5)NO4

45 compounds × 3 adducts ([M+H]+, [M+Na]+, [M+NH4]+) = 135 targets,
searched at 20 ppm mass tolerance. n≥12 was tagged `long_chain`, matching
the PI's stated interest in long-chain AHLs; n<12 was kept as an in-search
comparison range.

**Result**: 103 raw feature matches (41 long-chain, 62 short/medium-chain)
within 20 ppm. Files: `ahl_target_list.csv` (targets), `ahl_mass_matches.csv`
(matches).

**What this result does not mean**: an exact-mass match is not a
structure confirmation. The AHL formula family (CxHyNO3 or CxHyNO4, MW
roughly 170–340) is a common, non-distinctive elemental-composition
space. No decoy or permutation control has been run to measure how many
matches a random mass list of the same size would produce in this same
53,040-feature table at the same tolerance. Until that control is run,
the false-positive rate of this search is unknown, and the 103 matches
should be read as unfiltered candidates, not as detections.

**Structural cross-check**: 39 of the 103 matched rows have an
independent SIRIUS structural call. The single closest label to AHL
chemistry is row 28552, matched as `C5-HSL` (short/medium-chain, not
long-chain), which SIRIUS calls "N-acyl amines." But SIRIUS's own
assigned formula for that row, C11H13NO3, does not match the target AHL
formula, C9H15NO3, that the mass search aimed at. This means the mass
match at that row is very likely a coincidental overlap with an unrelated
compound (4-(cyclopropanecarbonylamino)-3-methylbutanoic acid per
SIRIUS), not corroborating evidence for an AHL. No other matched row has
a SIRIUS class resembling an acyl-lactone or acyl-amide signaling
molecule. The rest are called dipeptides, amino acids, fatty acyl
carnitines, chalcones, and terpenoids — chemical classes structurally
incompatible with AHLs.

## 4. Strain-level detection versus a colony-texture proxy (Phase 2)

This step asked a narrower question: among the 41 long-chain candidate
features, do any strains show a clear detection signal, and does that
signal track with a colony-surface-texture proxy (built separately, see
`AHL_AUTOINDUCER_SEARCH.md` for the texture proxy's own provenance and
caveats — it stands in for "smooth vs. rough," itself thought to be a
proxy for capsule production, and is not a validated phenotype).

**Binary detection is uninformative.** "Detected" was defined as a
feature's peak area exceeding its own blank-sample floor. Under this
definition, 266 of 266 strains with MS data (100%) show at least one of
the 41 long-chain candidates as "detected." A blank-floor-only threshold
has no power to discriminate in a table this large: with 41 narrow mass
windows checked against a 16,000+-feature background, almost any sample
will clear a blank floor somewhere by chance. This result should not be
read as evidence of widespread AHL production.

**Continuous intensity check is null.** Using the raw (not abundance-
normalized) peak area of each strain's single most intense long-chain
candidate feature: Spearman rho(intensity, texture-smoothness) = −0.097,
p = 0.12, n = 260. Not significant. The sign runs opposite the
hypothesis (higher candidate intensity trends with a rougher-reading
texture score, though this is within noise). The 20 strains with the
highest candidate intensity have a lower median smoothness score (−0.39)
than the full tested panel (+0.16) — again descriptive, not a formal
test, and not adjusted for the known confound between raw peak area and
colony size documented elsewhere in this project.

Full output: `ahl_strain_detection_vs_morphology.csv` and its
diagnostics file.

## 5. Cross-project GWAS follow-up (genotype, not MS data)

This is not an MS-data result, but it bears on the same hypothesis, so it
is included for completeness.

A dedicated GWAS was run in the sibling `Rhodotorula_phenotypes` project
(which has an existing, validated GEMMA-based GWAS pipeline for *R.
mucilaginosa*) to ask whether any genetic locus explains colony-texture
variation. Thirteen raw Haralick/GLCM texture measurements, already
present in that project's own imaging data, were tested as separate GWAS
traits across 213 and 182-strain panels.

Several traits initially showed large counts of statistically significant
SNPs (Contrast: 3,885; SumVariance: 3,056; SumAverage: 1,518). Checking
each hit's carrier strains against the existing kinship matrix showed
every one of these traces to a small group of closely related strains
(pairwise relatedness in the top 2–4% of the whole panel), not an
independent genetic locus. One more common-variant hit (`SumAverage`)
turned out to sit at the same location as an already-known locus for
colony lightness — expected, since `SumAverage` is literally average
image brightness, not a distinct texture signal.

**Result: no credible genetic locus for colony texture was found.** Full
detail: `~/projects/Rhodotorula_phenotypes/analysis/gwas/GWAS.md` section
24, decisions D-32 and D-33 in that project's `.living/decisions.md`.

## 6. What has not been done

Stated explicitly, per this project's standard for not overclaiming:

1. **No decoy/permutation null has been run** on the 103 AHL matches or
   the 918 N-acyl amino acid matches. The false-positive rate of a 20 ppm
   search over either formula family, in this specific feature table, is
   not known.
2. **No MS2 fragment confirmation has been attempted** on any candidate
   in either search. The diagnostic AHL fragmentation pattern
   (homoserine-lactone ring loss) has not been checked against any of the
   103 AHL rows' MS2 spectra; nor has any N-acyl amino acid candidate's
   MS2 spectrum been inspected.
3. **No search was run for AI-2, DSF, or the native fungal quorum-sensing
   molecules** (farnesol, tyrosol, phenylethanol, tryptophol). These are
   mentioned in this report only as literature background.
4. **No PI-validated smooth/rough or capsule-production phenotype has
   been located** for this strain panel. The texture proxy used in
   Phase 2 is a stand-in, not a direct measurement.
5. **The cross-pipeline batch effect** between the texture-proxy imaging
   rig and this project's own color-phenotyping rig has not been
   checked.
6. **No strain-level detection analysis has been run for the N-acyl
   amino acid candidates** (the Phase-2-style step already done for
   AHLs). This has not been done for rows 4109, 51126, or 51152.

## 7. Recommended next step, if this is pursued further

For AHLs, the single highest-value next step is the decoy/permutation
null (item 1 above). Given four independent lines of evidence already
point away from the AHL hypothesis, this is offered as due diligence, not
as an investigation expected to reverse the current conclusion.

For N-acyl amino acids, the priority order is different, because this
search actually produced structurally corroborated candidates: (a) pull
MS2 spectra for rows 4109, 51126, and 51152 and check for the expected
fatty-acyl and amino-acid-backbone fragment ions; (b) run the same
strain-level detection-vs-texture analysis already built for AHLs on
these three rows; (c) the decoy/permutation null still applies and would
help calibrate how surprised to be by finding these matches at all.

## 8. N-acyl amino acid expansion (2026-09-16, PI request)

N-acyl amino acids share the same acyl-CoA/acyl-ACP-donor chemistry as
AHLs. Some LuxI-family acyltransferases are documented to make N-acyl
amino acids instead of, or alongside, AHLs. Full detail, method, and all
caveats: `NACYL_AMINO_ACID_SEARCH.md`. Summary:

**Method**: same design as the AHL search (exact-mass re-mining,
independent of SIRIUS, 20 ppm). Target list built systematically: all 20
standard proteinogenic amino acids plus homoserine (motivated by this
project's own finding that *R. mucilaginosa* has lactonase activity — see
section 2 — so the open-chain hydrolysis product of a bacterial AHL
contaminant is a directly motivated target), each crossed with acyl chain
length C4-C18 and 3 positive-mode adducts. 945 targets total. Negative-
mode ions (which the free carboxylic-acid group would favor) were not
searched — checked directly, and this project's raw feature table
contains only positive-mode data.

**Result**: 918 raw matches across 541 distinct feature rows. Same
"unfiltered, no decoy null" caveat as the AHL search applies, and if
anything a higher expected false-positive rate, since this formula space
is even less distinctive than the AHL family's.

**Structural cross-check — the one meaningfully different result in this
whole investigation**: filtering to rows where SIRIUS's own independently
assigned formula exactly matches the target formula narrows 541 rows to
25. Most of the 25 are chemically incompatible on inspection (wrong
structure class despite matching formula — the same coincidental-overlap
pattern seen throughout the AHL search). But three rows are genuinely
different: SIRIUS's own structure name, not just its formula, names the
exact target compound —

- row 4109, target N-C16:0-arginine, SIRIUS: **"Palmitoyl arginine"**
- row 51126, target N-C14:0-arginine, SIRIUS: **"N2-(1-Oxotetradecyl)-L-arginine"** (= N-myristoyl-arginine)
- row 51152, target N-C16:0-lysine, SIRIUS: **"N6-Palmitoyl lysine"**

This is a materially stronger correspondence than the AHL search
produced (where the one formula-adjacent candidate turned out to be an
unrelated, mismatched compound). A second small cluster of exact-formula
matches was inspected and set aside as likely non-biological — one of the
proposed SIRIUS structures is named as a deuterium-labeled compound,
consistent with a spectral-library reference/internal-standard entry
rather than something present in the sample.

**What this does and does not show**: N-acyl amino acids are ordinary
lipid chemistry for a fungus (related to ceramide/sphingolipid
acyltransferase activity), so finding one is far less surprising than
finding an AHL — but that also makes it much weaker evidence of a
*signaling* function specifically. This result says "this chemical class
plausibly exists in the extract," not "this strain uses it to signal."

**Follow-up literature check and compartment/species analysis
(2026-09-16)**: no paper documents Palmitoyl arginine, N-myristoyl-
arginine, or N6-Palmitoyl lysine as characterized natural products from
any organism, fungal or bacterial. The closest documented relatives are
ornithine lipids and lysine lipids — a real, established bacterial
"aminolipid" class, described in the literature as **exclusively
bacterial** and produced as a **membrane-lipid substitute under phosphate
starvation**, not as a signal. No quorum-sensing role is documented for
this compound class in any organism.

Pulling per-sample data for these three rows (plus two weaker histidine
candidates) and mapping to strain, species, and fraction (cell vs.
supernatant) found: **all 5 candidates are overwhelmingly cell-
associated**, detected in 155-252 of 265 cell samples but only 1-10 of
266 supernatant samples (paired Wilcoxon signed-rank p<0.0001 for every
candidate, same direction each time). A diffusible quorum signal has to
leave the cell to work; a membrane lipid does not. This compartment
pattern matches the membrane-lipid role documented for the closest known
relative class, not a signaling role. Production does vary significantly
across species in the cell fraction (Kruskal-Wallis p=0.001-0.009,
descriptive, not phylogenetically corrected) — *R. taiwanensis* and
*R. paludigena* consistently highest, some species near-absent for some
compounds. Full detail: `NACYL_AMINO_ACID_SEARCH.md`.

## Files in this analysis

| File | Contents |
|---|---|
| `AHL_AUTOINDUCER_SEARCH.md` | Working log of the AHL investigation, in the order it happened |
| `ahl_target_list.csv` | 135 AHL mass targets (45 compounds × 3 adducts) |
| `ahl_mass_matches.csv` | 103 raw, unfiltered AHL mass matches at 20 ppm |
| `NACYL_AMINO_ACID_SEARCH.md` | Full method, results, and caveats for the N-acyl amino acid expansion |
| `nacyl_amino_acid_target_list.csv` | 945 N-acyl amino acid mass targets |
| `nacyl_amino_acid_mass_matches.csv` | 918 raw, unfiltered matches at 20 ppm (541 distinct rows) |
| `strain_texture_table.csv` | 298-strain colony-texture proxy table |
| `ahl_strain_detection_vs_morphology.csv` | Phase 2 detection-vs-texture output |
| `figures/fig1-4_*` | Texture variance by species and by metric |
| `../scripts/ahl_targeted_mass_remining.py` | Phase 1 mass search |
| `../scripts/build_strain_texture_table.py` | Texture-proxy build |
| `../scripts/ahl_strain_detection_vs_morphology.py` | Phase 2 detection-vs-texture analysis |
| `../scripts/ahl_morphology_variance_figures.py` | Figure generation |
