# AHL Autoinducer Search

## Hypothesis under test
PI hypothesis (2026-09-11): long-chain N-acyl-L-homoserine lactones (AHLs
— canonical bacterial quorum-sensing autoinducers, e.g. 3-oxo-C12-HSL) vary
across *Rhodotorula* strains/species and their production is linked to
smooth vs. rough colony morphology.

## Status: Phase 1 (mass search) + Phase 2 (strain-level detection vs. morphology proxy) done. Current evidence is null / non-discriminating. See "Phase 2 results" below. The likely source publication for the hypothesis (Wilson et al. 2025, ISS *R. mucilaginosa*) and a broader fungal-AHL-biosynthesis literature check (both below) independently weaken the premise further.

## Primary source of the hypothesis, and why it should be discounted (2026-09-11)

PI identified the likely origin of the "Rhodotorula produces AHL" premise:
**Wilson et al. 2025, "Characterization of virulence-related phenotypes of
*Candida parapsilosis* and *Rhodotorula mucilaginosa* isolated from the
International Space Station (ISS)", *Life Sciences in Space Research*
45:16-24, doi:10.1016/j.lssr.2025.01.002** (PMID 40280638). Full abstract
retrieved (PubMed); paper itself not accessed (paywalled).

**What that paper actually did, per its own abstract**: "autoinducer (AI)
production was detected by activation of a reporter fluorescent gene
present in biosensor bacterial strains" -- i.e. a **bacterial bioreporter
functional assay**, not a chemical/structural identification (no MS, no
NMR reported in the abstract). It reports "increased ... long-chain
autoinducer production ... in the ISS isolates" alongside increased
capsule production, biofilm formation, antifungal resistance, and
nematode virulence in ISS *R. mucilaginosa* vs. ATCC controls -- an
overall "enhanced virulence phenotype" narrative, not a specific
AHL-vs-morphology test. No smooth/rough colony-texture score is reported
in the abstract either; its morphology-adjacent readouts are capsule
production and filamentation.

**Why this weakens the premise, independent of this project's own null MS
result**: a bacterial LuxR-type AHL bioreporter is a functional assay, and
such reporters are well documented in the quorum-sensing literature to
cross-react with non-AHL lipophilic ligands (certain fatty acids and
fatty-acid-like compounds can activate LuxR-family receptors
non-specifically). *Rhodotorula* is a prolific lipid/fatty-acid-derivative
producer, so a positive bioreporter signal from a fungal extract does not
by itself demonstrate a genuine AHL was present -- it is at least as
consistent with cross-reactive fungal lipid chemistry, or residual
bacterial AHL carried over in a non-axenic culture, as with true fungal
AHL biosynthesis. Combined with the KEGG/Pfam absence of any fungal LuxI
ortholog and the specific published finding that *R. mucilaginosa* is an
AHL-**degrading** (lactonase-positive) organism (both below), the
bioreporter signal in the ISS paper is better explained by one of those
two alternatives than by novel fungal AHL biosynthesis. This does not
mean the ISS finding is wrong on its own terms (a real bioreporter signal
was presumably observed) -- it means "long-chain autoinducer production"
in that paper's specific, narrow (functional-bioassay) sense should not be
read as "this fungus makes AHLs" in the structural-chemistry sense this
project's MS mass-search was built to test.

## Background literature check (2026-09-11): is fungal AHL biosynthesis even known to exist?

Before pursuing this further, a background research pass checked whether
any fungus has ever been shown to biosynthesize AHLs at all (independent of
this project's own MS data). Summary (full citations at the end of this
section):

- **No LuxI-type AHL synthase ortholog exists in any fungal genome per
  KEGG.** The relevant KEGG ortholog (K18096, the CoA-utilizing LuxI
  homolog) is restricted to alpha-Proteobacteria; no fungal/eukaryotic
  genome is annotated with it, and no Pfam/InterPro entry for a
  "homoserine-lactone synthase" domain family in fungi was found.
- **No natural fungal AHL production has been cloned or characterized.**
  The one directly relevant paper (*Commun. Biol.* 2025, engineering an
  AHL-based quorum-sensing circuit in *S. cerevisiae*) is explicitly an
  **engineered, heterologous** system (bacterial LuxI expressed in yeast)
  — the paper itself frames this as something yeast does not do naturally.
- **Directly relevant counter-finding: *Rhodotorula mucilaginosa* — this
  panel's dominant species (213/298 strains in the texture table) — is
  specifically published as an AHL-DEGRADING (quorum-quenching) yeast**,
  with lactonase activity confirmed against C6-HSL, 3-oxo-C6-HSL, and
  3-hydroxy-C6-HSL (Tan et al. 2014, *Sensors* 14:6463, PMC4029656).
  AHL-inactivating (lactonase/acylase) activity is reported as widespread
  across yeasts generally, and the basidiomycetous yeast *Trichosporon
  loubieri* was shown to consume AHLs as a carbon/nitrogen source (PMC3859043).
  This is the **opposite** direction from the hypothesis: this genus is
  documented to destroy AHLs, not make them.
- Fungi do have their own, chemically unrelated quorum-sensing chemistry
  (farnesol, tyrosol, phenylethanol, tryptophol — aromatic
  alcohols/sesquiterpenes controlling yeast-hyphal morphogenesis;
  Wongsuk et al. 2016, *J. Basic Microbiol.*) — not homoserine lactones.

**Implication**: this reframes any AHL-mass signal in the existing MS data
as more parsimoniously explained by trace bacterial contamination
(possibly since degraded by the strain's own lactonase activity, which
would also help explain why no clean, discriminating detection signal was
found in Phase 2) than by novel fungal biosynthesis. It does not rule out
the hypothesis — a genuinely novel fungal AHL pathway would by definition
not yet be in KEGG/Pfam — but there is currently no positive-control gene
to anchor a genomic search against, and the literature gives a more
likely alternative explanation for the null Phase 2 result than "wrong
detection method."

**Key sources**: KEGG KO K18096; Tan et al., *Rhodotorula mucilaginosa*, a
Quorum Quenching Yeast Exhibiting Lactonase Activity, *Sensors* 2014,
14:6463 (PMC4029656); Ren et al., *Trichosporon loubieri* AHL degradation
(PMC3859043); "Inactivation of bacterial quorum sensing signals N-acyl
homoserine lactones is widespread in yeasts" (2018); Engineering an
AHL-based quorum-sensing circuit in *S. cerevisiae*, *Commun. Biol.* 2025
(PMC12717128, engineered/heterologous); Wongsuk et al. 2016, *J. Basic
Microbiol.* (fungal QS chemistry: farnesol/tyrosol/phenylethanol/tryptophol).

## Morphology-variance figures (2026-09-11)

`analysis/scripts/ahl_morphology_variance_figures.py` -> `figures/`:

- `fig1_smoothness_by_species.png/pdf` — smoothness_z by species (boxplot +
  jittered strip for n>=3 species; diamonds for the 8 singleton/n=2
  species, where a box is not meaningful). One pathological outlier
  (`TFCN_43A-4`, smoothness_z=-43, driven by an extreme raw Contrast
  value) is excluded from the plotted axis only (named in the title; the
  underlying data is untouched).
- `fig2_smoothness_distribution.png/pdf` — strain-level smoothness_z
  histogram (n=297 plotted, same outlier excluded and named), roughly
  unimodal, median 0.13.
- `fig3_raw_texture_components_by_species.png/pdf` — the 4 raw Haralick
  components feeding smoothness_z (Contrast, Entropy, AngularSecondMoment,
  InverseDifferenceMoment), by species, so the composite's apparent
  variance isn't taken on faith.
- `fig4_all_texture_metrics_distribution.png/pdf` — all 13 raw
  `-avg-scale05` Haralick metrics' distributions across all 298 strains
  pooled (independent of species) — most are unimodal and reasonably
  symmetric; Contrast is heavily right-skewed (the same outlier strain
  drives its long tail).

Species-level spread in `smoothness_z` (n>=3 species) is modest and
overlapping (interquartile ranges mostly span roughly -2 to +2, medians
clustered near 0) — there is no visually obvious species with a
dramatically different texture profile from the rest, *R. mucilaginosa*'s
enormous sample size (213/298) aside.

## What has been done
`analysis/scripts/ahl_targeted_mass_remining.py` — exact-mass re-mining of
the raw EB feature table (`aligned_features_ms2.csv`, 53,040+ rows, same
table used by `idea1_targeted_mass_remining.py`), independent of SIRIUS's
own class calls. Target list generated programmatically from the general
AHL formula, swept over acyl chain length n = C4-C18:

- unsubstituted Cn-HSL: C(n+4)H(2n+5)NO3
- 3-oxo-Cn-HSL: C(n+4)H(2n+3)NO4
- 3-hydroxy-Cn-HSL: C(n+4)H(2n+5)NO4

45 compounds x 3 adducts ([M+H]+, [M+Na]+, [M+NH4]+) = 135 targets,
searched at 20 ppm tolerance. `n>=12` tagged `long_chain` per the stated
hypothesis; `n<12` retained in the same search as a comparison range.

**Result**: 103 raw feature matches (41 `long_chain`, 62
`short_medium_chain`) across the 20 ppm window. See
`outputs? -> ahl_mass_matches.csv` (this dir) and `ahl_target_list.csv`.

## Important caveat — not yet a finding
AHLs are canonical **bacterial** signals (LuxI/LuxR family), not a
previously-reported fungal metabolite class. A hit here has at least three
non-exclusive explanations, and this search cannot distinguish them:

1. genuine fungal AHL production or AHL-mimic biosynthesis (would be
   novel and publishable if confirmed),
2. trace bacterial contamination of the axenic culture/extract,
3. coincidental isobaric overlap with an unrelated fungal metabolite —
   the AHL formula family (CxHyNO3/NO4, MW ~170-340) is a common,
   non-distinctive elemental-composition space, so 20 ppm mass-only
   matching has a real false-positive rate here that has **not yet been
   quantified** (no permutation/decoy-mass control has been run against
   this target list, unlike this project's other targeted-mass work,
   e.g. `phase3_metabolome_phenotype_idea1` and
   `phase_siderophore/siderophore_mass_remining.py`).

**Before treating any of the 103 hits as a lead**, per this repo's
robust-analysis convention, the following are required, not optional:
- a decoy/permutation null (random mass targets in the same m/z range,
  same ppm tolerance, same adduct set) to estimate the expected hit count
  under no biological signal,
- MS2 spectral inspection of any surviving candidate against the
  diagnostic AHL fragmentation pattern (homoserine lactone ring loss,
  m/z 102.055 immonium-type fragment), since exact mass alone does not
  confirm structure,
- cross-check against SIRIUS's own class call for the same row ID (already
  joined into `ahl_mass_matches.csv` where available) — a structurally
  incompatible SIRIUS call (e.g. a lipid/glycerophospholipid class) is
  evidence against the AHL interpretation for that row, not for it.

None of this has been done yet. The 103 rows are unfiltered raw matches.

## Known gap: no colony-morphology (smooth/rough) phenotype found
Searched this repo's `data/` and `analysis/` trees, and the sibling
phenotyping projects reachable from this session
(`Rhodotorula_Phenotyping`, `ExRhodotorula_Phenotypes`, `Phenotypes_Study1`)
for `smooth|rough|morpholog|texture|mucoid|colony.?(shape|form|surface|
margin|elevation)`. No categorical colony-morphology phenotype (smooth vs.
rough, or any texture/surface trait) was found anywhere accessible in this
session. The phenotype data currently available in this project
(`control_phenotype_90_110h`, `EXFAB_UCR-005`) covers only CIELAB color
(L*/a*/b*, chroma, hue), colony area/shape-area, and copper-media growth
AUC — no surface-texture trait. **This is stated as a gap, not assumed to
be resolved or absent** — the PI may have this scored elsewhere (image
archive not yet OCR'd/scored, a spreadsheet not yet ingested, or a call
that needs to be made from existing colony images). This blocks any
morphology-side of the hypothesis until sourced.

## Phase 2 results (2026-09-11): strain-level detection vs. colony-texture morphology proxy

Per PI direction, this phase looked for direct positive evidence (which
strains show the candidate features, and whether they skew toward
smoother texture) rather than running the project's full
phylogenetically-aware block-permutation framework. That heavier
framework is deferred, not abandoned -- everything below is descriptive/
exploratory and NOT adjusted for species/clade structure.

**Morphology proxy**: `strain_texture_table.csv` (298 strains, from
`salinity_texture_baseline`, see below) reduced to a single
`smoothness_z` composite = z(AngularSecondMoment) + z(InverseDifference-
Moment) − z(Contrast) − z(Entropy) (same z-score-sum convention as this
project's existing `orange_score`).

**AHL-candidate detection**: `analysis/scripts/ahl_strain_detection_vs_morphology.py`
pulled per-sample peak areas for the 41 long-chain (n≥12) candidate row
IDs from the raw EB feature table, mapped MS sample columns to
`strain_code` via the `Cu_AUC` crosswalk (ID lookup only, not a phenotype
source), and defined "detected" as peak area above that feature's own
blank-sample floor.

**Result 1 — binary detection is uninformative**: 266/266 strains
(100%) with MS data show at least one of the 41 long-chain candidates
above its blank floor. This is not evidence of widespread AHL production
-- it means a blank-floor-only threshold has no discriminating power in
a 16k+-feature untargeted table: essentially any sample has *some* peak
above blank noise somewhere across 41 narrow mass windows. This binary
check is retained in the output for transparency but should not be read
as a finding.

**Result 2 — continuous intensity check is null**: using the raw (non-
abundance-normalized) peak area of each strain's single most intense
long-chain candidate feature, `Spearman rho(max_peak_area_long_chain,
smoothness_z) = -0.097, p = 0.12, n = 260`. Not significant, and the
sign runs opposite the hypothesis (higher candidate intensity trending
with *lower*, i.e. rougher-reading, smoothness_z, though this is well
within noise). The top-20-by-intensity strains have a lower median
smoothness_z (-0.39) than the full tested panel (+0.16) -- again
descriptive, not a formal test, and not adjusting for the known
raw-peak-area/colony-size confound documented elsewhere in this project.

**Result 3 — the single "best" SIRIUS cross-reference argues against, not for, AHL identity**: of the 103 raw mass hits, 39 have an independent SIRIUS structural call; the closest thing to an AHL-adjacent class label is row 28552 (`C5-HSL`, short/medium-chain, not even long-chain) called "N-acyl amines" -- but SIRIUS's own assigned formula for that row (`C11H13NO3`) does not match the target AHL formula (`C9H15NO3`) the mass search was aiming at, meaning the match is very likely a coincidental isobaric overlap with an unrelated compound (4-(cyclopropanecarbonylamino)-3-methylbutanoic acid), not corroborating evidence. No other candidate row has any SIRIUS class resembling an acyl-lactone/amide signaling molecule; the rest are dipeptides, amino acids, fatty acyl carnitines, chalcones, terpenoids, etc. -- all structurally incompatible with AHLs.

**Bottom line as of 2026-09-11**: no positive evidence for the hypothesis has been found. The mass-search hits do not show a discriminating detection pattern, the one available structural cross-check argues against AHL identity for its best candidate, and the continuous intensity-vs-smoothness relationship is null. Independently, the literature source that likely motivated the hypothesis (Wilson et al. 2025, ISS *R. mucilaginosa*, see above) rests on a bacterial bioreporter functional assay, not a chemically-confirmed AHL structure, and bioreporters of this type are known to cross-react with fungal lipid chemistry -- so the premise itself is weaker than "a paper detected AHLs in this species" would suggest. This does not rule out the hypothesis outright -- MS2 spectral confirmation (the diagnostic AHL fragment) has still not been attempted on any candidate -- but combined with the biosynthesis-pathway absence and the degradation-not-production literature, the balance of evidence currently weighs against it.

## Next steps (open)
1. Run the decoy/permutation null on the 103 raw hits, to formally
   quantify the false-positive rate of a 20 ppm search over this AHL
   formula family (script not yet written) -- would explain, and put a
   number on, the 100% binary-detection saturation found in Phase 2.
2. If a decoy-adjusted candidate still stands out, pull its MS2 scan and
   check against the AHL diagnostic fragment (homoserine-lactone ring
   loss) before trusting it structurally -- not yet attempted for any
   row, since none has cleared a basic plausibility bar yet.
3. Cross-pipeline batch-effect check between the Salinity-screen imaging
   rig and this project's own color-phenotyping rig (flagged as unchecked
   at ingestion) -- relevant if the texture proxy is used further.
4. If a future step needs a rigorous (not descriptive) test of any
   surviving candidate against morphology, use this project's established
   phylogenetically-aware block-permutation framework (same pattern as
   `phase2_metabolome_phenotype`) -- deferred in this pass per PI
   direction, not because it's unnecessary in general.
5. **Still open**: no PI-validated smooth/rough categorical phenotype
   has been located; `strain_texture_table.csv` remains a proxy.

## Files
- `analysis/scripts/ahl_targeted_mass_remining.py` — Phase 1 target-list builder + search
- `ahl_target_list.csv` — 135 target m/z values (45 compounds x 3 adducts)
- `ahl_mass_matches.csv` — 103 raw, unfiltered feature matches at 20 ppm
- `analysis/scripts/build_strain_texture_table.py` — strain-level texture (morphology proxy) aggregation
- `strain_texture_table.csv` / `strain_texture_table_diagnostics.txt` — 298 strains, 13 avg-scale05 Haralick metrics
- `analysis/scripts/ahl_strain_detection_vs_morphology.py` — Phase 2: strain-level AHL detection vs. smoothness_z
- `ahl_strain_detection_vs_morphology.csv` / `_diagnostics.txt` — Phase 2 full output + log (null result, see above)
