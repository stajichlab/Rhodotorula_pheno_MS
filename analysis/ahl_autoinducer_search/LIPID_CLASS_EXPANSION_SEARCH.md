# Structurally-related lipid class expansion (2026-09-16)

Follow-up to `NACYL_AMINO_ACID_SEARCH.md`'s confirmed N-acyl-arginine/
lysine candidates: PI asked whether other structurally similar lipid
classes were worth considering. Three were identified and searched, with
different motivations and different results.

## 1. Diacylated ornithine/lysine lipids — NEGATIVE

**Motivation**: the confirmed mono-acyl candidates match the
INTERMEDIATE of the documented bacterial ornithine-lipid/lysine-lipid
pathway (OlsB→OlsA). The mature, biologically active form of that
pathway is diacylated — a second fatty acid esterified onto the first
acyl chain's 3-hydroxyl.

**Method**: `analysis/scripts/diacyl_aminolipid_mass_remining.py`.
Formula built from first principles (headgroup + 3-hydroxy fatty acid 1
− H2O + fatty acid 2 − H2O), ornithine and lysine headgroups, both acyl
chains swept systematically over C10–C18 (294 targets).

**Result**: 53 raw matches, 8 distinct rows, **zero** with a SIRIUS
annotation of any kind (let alone a matching one), and none share a
molecular-network component with the confirmed mono-acyl candidates.
Several ppm errors are also mediocre (up to ±19 ppm, near the search
tolerance boundary). **No support for the mature diacylated form in
this dataset.** This is itself informative: the mono-acyl form has real,
multi-method support; the predicted next biosynthetic step does not,
which is a real constraint on how far along the pathway (if bacterial)
these strains' contaminant/lipid signal actually goes, if that
interpretation is correct at all.

## 2. Fungal ceramides (N-acyl sphingoid bases) — POSITIVE, and structurally distinct from the aminolipid story

**Motivation**: unlike the aminolipids (documented as bacterial), N-acyl-
sphingoid-base lipids are textbook, well-documented FUNGAL membrane
chemistry — this checks whether the broader "N-acyl amines" SIRIUS-class
signal in this dataset includes real fungal sphingolipid chemistry,
not only (presumptively bacterial) amino acid conjugates.

**Method**: `analysis/scripts/fungal_ceramide_mass_remining.py`.
Phytosphingosine and dihydrosphingosine long-chain bases, plain and
2-hydroxy fatty acids, chain lengths C14–C26 (fungal ceramides commonly
carry very-long-chain fatty acids) — 84 targets.

**Result**: 23 raw matches, 21 distinct rows; **7 rows have a SIRIUS
formula exactly matching the target** (a script bug was caught and fixed
here — SIRIUS omits the "1" subscript in single-atom counts, e.g. "NO3"
vs. this project's "N1O3"; the corrected, element-count-normalized
comparison found 7 real matches where the first pass had reported 0).
**6 of the 7 are independently classed "Ceramides" by SIRIUS**, three
with real structure names:

- **"Armillaramide"** (row 31296) — a genuine, previously published
  natural-product name, isolated from *Armillaria* (honey fungus). A
  named fungal ceramide, not a generic/systematic label.
- **"2-Octadecanoylaminohexadecane-1,3-diol"** (rows 6274, 14537, 37800,
  three features at different retention times/adducts) — systematic
  name for a C16-sphingoid-base/C18-fatty-acid ceramide (same gross
  formula as this project's "C18 base + C16 acyl" target, an expected
  base-vs-acyl ambiguity in ceramide mass search that doesn't affect the
  formula-level confirmation).
- **"N-(1,3-dihydroxyoctadecan-2-yl)tetradecanamide"** (rows 41699,
  42657) — systematic name for a C18-sphingoid-base/C14-fatty-acid
  ceramide.

**MS2/network check**: rows 41699 and 42657 (the C14 ceramide pair, two
different adducts) sit in the SAME molecular-network component (1301) —
internally consistent. The C16 ceramide's three instances (6274, 14537,
37800) do NOT all share a network component, a minor inconsistency (could
mean genuinely different retention-time isomers, or just imperfect
network-edge thresholding) — not fatal to the finding, but not perfectly
clean either.

**Reading**: this is a real, independently-corroborated finding for a
genuinely fungal-plausible lipid class — unlike the arginine/lysine
aminolipids (bacterial-documented, fungal origin questionable), ceramide/
sphingolipid biosynthesis is normal, expected chemistry for a
basidiomycete yeast. This does not bear on the AHL/quorum-sensing
question at all (ceramides are not documented signaling molecules in
that sense) — it answers a narrower question (is some of the broader
"N-acyl amines" signal genuine fungal lipid biology) with yes.

## 3. Broader N-acyl amide family — mostly noise, but one credible, fungal-plausible hit

**Motivation**: least specifically motivated of the three — a systematic
sweep of the wider N-acyl-amide natural-product space documented in gut-
microbiome literature (Cohen et al. 2017: N-acyl ethanolamines, taurines,
GABA conjugates, and biogenic-amine conjugates), rather than a targeted
follow-up on an existing lead.

**Method**: `analysis/scripts/nacyl_amide_broad_mass_remining.py`.
Ethanolamine, taurine, GABA, dopamine, serotonin, tryptamine headgroups,
acyl chains C8–C22 (270 targets).

**Result**: 191 raw matches, 180 distinct rows, 29 with an exact-matching
SIRIUS formula. **Most of these 29 are the same isobaric-coincidence
pattern seen throughout this investigation** — e.g. "N-C15:0-GABA"
matched by SIRIUS's "Palmitoyl sarcosine" (a different, unrelated
lipoamino acid that happens to share the same gross formula), several
GABA-target rows matched by unrelated ceramide names, several dopamine-
target rows matched by unrelated heterocycles. **Two rows are a real,
credible, headgroup-consistent exception**:

- **Row 19360, target N-C16:0-ethanolamine, SIRIUS: "Palmitoylethanolamide"
  (class: "N-acyl ethanolamines (endocannabinoids)")** — exact name and
  class match. **MS2-checked**: the real spectrum's BASE PEAK is m/z
  62.0599 (calc. protonated ethanolamine [C2H8NO]+ = 62.0600, <2 ppm) —
  the free-headgroup fragment, the same diagnostic logic used to confirm
  the arginine/lysine candidates — plus a clear [M+H−H2O]+ loss-of-water
  fragment at m/z 282.2788 (calc. 282.2791), the other classic NAE
  fragmentation. Both expected fragments present, strong intensity.
  **This is a genuine, MS2-confirmed compound identity, independent of
  the arginine/lysine family.**
- **Rows 15924 and 32948, target N-C18:0-ethanolamine, SIRIUS:
  "N-(2-hydroxyethyl)octadecanamide"** (= stearoylethanolamide, the C18
  homolog of the row-19360 compound) at two separate features — not yet
  MS2-checked, but the same exact-name match pattern.

**Reading**: N-acylethanolamines (NAEs) are a well-documented,
broadly-conserved lipid class across plants, fungi, and animals (not
primarily bacterial the way the aminolipids are) — this is the first
candidate in this whole investigation where genuine fungal biosynthesis
is chemically the MORE parsimonious explanation, not merely
un-excluded. NAEs do have documented signaling-adjacent roles in some
organisms (membrane-stress response, developmental signaling in plants;
a role in fungal growth/development has been proposed in some
literature) — **this specific claim has not been freshly verified this
session** and should be checked before being treated as established, the
same way the AHL/aminolipid literature was verified earlier in this
investigation.

## Overall bottom line

Of three additional lipid classes considered: one (diacyl aminolipids) is
negative, one (fungal ceramides) is positive and adds a real, separate,
fungal-plausible finding unrelated to quorum sensing, and one (the broad
N-acyl amide sweep) is mostly noise but surfaces a single credible,
MS2-confirmed, fungal-plausible compound (Palmitoylethanolamide) that is
chemically distinct from the bacterial-aminolipid story. None of these
three additional classes are documented quorum-sensing molecules — they
answer "what real lipid chemistry exists in this data" (a question worth
asking given how much of the earlier work concerned distinguishing real
compounds from mass coincidences), not "is there a quorum signal here."

## Files
- `analysis/scripts/diacyl_aminolipid_mass_remining.py` / `diacyl_aminolipid_target_list.csv` / `diacyl_aminolipid_mass_matches.csv`
- `analysis/scripts/fungal_ceramide_mass_remining.py` / `fungal_ceramide_target_list.csv` / `fungal_ceramide_mass_matches.csv`
- `analysis/scripts/nacyl_amide_broad_mass_remining.py` / `nacyl_amide_broad_target_list.csv` / `nacyl_amide_broad_mass_matches.csv`
