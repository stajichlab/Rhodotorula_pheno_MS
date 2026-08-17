# Strategy: Predicting Pigment Compounds and Genomic Correlates from Colony Color Across *Rhodotorula*

## Motivation (revised 2026-08-15 — see "Reframing" below for what changed and why)

**Core question**: can colony color (CIELAB L\*/a\*/b\*), measured across the
full panel of ~300 phenotyped strains spanning 17-18 *Rhodotorula* (+ 2
outgroup) species, be used to build a classifier that predicts or supports
which metabolome compounds — specifically pigment-related compounds
(carotenoids for orange/yellow; other red pigments possibly related to
torulene/torularhodin) — are enriched in a strain's chemical profile? And
can the compounds/pathways implicated this way be linked to genome-level
differences (gene/protein composition, pathway completeness, sequence
variants) across strains?

**Working hypothesis**: carotenoids (a well-established *Rhodotorula*
pigment class — torulene, torularhodin, β-carotene) underlie orange/yellow
coloration; other, chemically distinct red pigments may also contribute and
are not yet identified in this dataset.

**Why this is a whole-panel prediction/classification problem, not a
single-species comparison**: every phenotyped strain has a color, and (for
the MS-covered subset) a metabolome profile. The question is not "what
makes species X different" but "does color predict pigment chemistry
across the whole diversity of strains we have" — the natural unit of
replication is the ~300-strain panel (18 species), not one clade. This
matters for statistical power: a whole-panel regression/classification
uses far more of the data's information than a single-clade contrast does.
Phylogenetic non-independence among strains within a species (and among
closely related species) still needs correcting for — the machinery built
for that (species-level collapse, PGLS, block permutation) carries forward
unchanged from before — but the "explain one extreme clade" power ceiling
that constrained the original framing no longer applies in the same way.

---

## Reframing (2026-08-15): what changed from the original plan

The project originally started from the observation that *R. dairenensis*
looked like a color outlier and asked "what genomic features make
*R. dairenensis* uniquely darker orange." Phase 1 work (see
`analysis/integrated_analysis/phase1_phenotype/`, and the git history of
this document) surfaced two problems with that framing:

1. **Statistical**: *R. dairenensis* is a single monophyletic clade
   (n=8-9 strains), so any genomic/metabolomic feature distinguishing it
   from the rest of the tree is confounded with phylogeny — a classic
   "n=1 clade" comparative-biology problem that PGLS/PIC/block permutation
   can formally correct for but cannot resolve in principle (they can't
   manufacture a second independent origin of the trait).
2. **Empirical**: on the actual color data (reproduced across the original
   YPD2 phenotype table and an independently-rebuilt, more complete
   `control_90_110` table — see "Phenotype data provenance" below),
   *R. dairenensis* is only the **4th-highest of 17-18 species** on a
   darker-orange composite score, well behind *R. taiwanensis*,
   *R. sphaerocarpa*, and *R. glutinis*. The species the motivation section
   was built around is not obviously the phenotypic outlier.

Rather than force the analysis to explain one species that may not even be
the extreme case, **the question has been generalized**: use color across
the *whole* panel to predict pigment compound chemistry and its genomic
correlates. This is both more statistically defensible (whole-panel
regression instead of single-clade contrast) and more directly useful (a
classifier/predictive tool, not a one-species just-so story). *R.
dairenensis* remains one strain group of interest within this larger
analysis, but is no longer the sole justification for it.

All of the machinery built under the old framing — phenotype table
construction, species-level collapse, phylogenetic tree pruning,
phylogenetic-signal testing, block-permutation testing, the BFD genome
functional-annotation database — carries forward unchanged; only the
*framing* and the specific downstream phases change. The full prior
discussion (the "n=1 clade" analysis, the *R. dairenensis*-specific
enrichment questions, the original 5-phase plan) is preserved in
**Appendix: Superseded *R. dairenensis*-Specific Framing** at the bottom of
this document for provenance — it documents real, reusable methodological
review even though the framing it was written for has changed.

### Second review of the reframed plan (Fable, 2026-08-15)

The reframing was itself sent back to Fable for review before implementation
started. Verdict: **a real improvement, but it shrinks the n=1-clade
problem rather than eliminating it** — with only 17-18 species-level
lineages, a trait with strong phylogenetic signal (a\*/C\* already qualify)
will still be dominated by a handful of species at the color extremes
(*R. taiwanensis*, *R. sphaerocarpa*, *R. glutinis*), so a "whole-panel"
result could still substantially be a 2-3-species story. The review also
flagged that Phase 4's planned LOSO-CV classifier does not by itself rule
out the same circularity risk in the original Phase 4 ML stack — a
classifier can hit good LOSO accuracy by learning phylogeny *through*
color rather than any color-specific relationship — and that Phase 3's
compound-class enrichment test is only as strong as SIRIUS's 4.8% feature
annotation rate allows. The specific fixes from this review (leave-one/
few-species-out sensitivity checks, a concrete two-part negative-control
design, a phylogenetically-structured-decoy-trait baseline for the
classifier, an explicit annotation-coverage caveat, and additional
compound-ID methods beyond SIRIUS) are folded directly into Phases 2-4
and "Key Considerations" below, not kept as a separate list, so they are
read as part of the plan rather than an addendum easy to skip.

---

## Data Assets Inventory

| Data Layer | Source | Key Dimensions |
|---|---|---|
| **Phenotype (canonical)** | `control_90_110` window, read via `analysis/scripts/build_strain_phenotype_table.py --source control_90_110` (from the sibling `Rhodotorula_phenotypes` project — not yet ingested into `data/metadata/`, see below) | 303 strains × CIELAB L\*, a\*, b\* (+ derived C\*, h°, orange_score), colony area |
| **Phenotype (legacy)** | `data/metadata/EXFAB_UCR-005/YPD2_phenotypic.20260702.fixed.csv.gz` | 318 strains; used for cross-validation against the canonical source, see `analysis/examine_phenotype_calling/` |
| **Metabolome** | `data/processed/EB_20260130_ExFAB_Rhodo_Sup_and_Pellet/...aligned_features_ms2.csv.zst` | 16,332 MS2 features × 594 samples (cell pellet + supernatant) |
| **Compound annotation** | `analysis/sirius_annotation/sirius_annotations.tsv` | 790 of 16,332 features SIRIUS-annotated with structure guess + NPC pathway/class + ClassyFire class. **73 features under the Terpenoids NPC pathway, including 3-4 explicitly "Carotenoids (C40, β-β)" / "Apocarotenoids(ε-)" NPC-class calls** — confirms carotenoid-class chemistry is present and detectable in this dataset, though specific structure-name calls (not literally torulene/torularhodin/β-carotene) should be treated as class-level evidence, not confirmed structure IDs, per SIRIUS's known unreliability on this compound class |
| **Genome – Functional** | `BFD/db/BFD.duckdb` | 276 strains, **fully populated** across pfam (3.27M rows), merops (402K), cazy (58K), signalp (103K), targetp (240K), tmhmm (401K), wolfpsort (2.17M), predgpi (2.17M), idp (3.0M) — checked 2026-08-15, all 9 tables have `count(distinct species_prefix) = 276`; the earlier "only 1-2 strains populated" blocker is resolved |
| **Genome – Gene/Protein** | Same BFD database | gene_proteins (2.1M), gene_info, gene_transcripts, codon_freq, aa_freq |
| **Genome – Assembly Stats** | Same BFD database | asm_stats, busco_genome (`complete_pct` etc., all 276 strains), telomere_summary |
| **Phylogeny** | `BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile` (PHYling protein tree, 278 strain-level tips) | Pruned to one tip per species via `analysis/scripts/prune_species_tree.R` → `analysis/integrated_analysis/phase1_phenotype/species_tree.nwk` (16 species with genome data) |
| **Existing Analyses** | `analysis/` directory | Color PCoA (`color_phenotype_ordination/`), MS feature PCoA (`ms_feature_ordination/`), differential features by species pair (`differential_features/`), SIRIUS annotations |

### Strain Overlap
- **Phenotype (canonical, control_90_110)**: 303 strains, 1 missing species call
- **Phenotype (legacy YPD2)**: 318 strains, 15 missing species call
- **MS data**: 594 samples from ~318 strains (cell + supernatant pairs)
- **BFD genome**: 276 strains, fully functionally annotated, 17-18 species represented
- Cross-referencing genome coverage against strain lists is naming-convention-sensitive (strain codes were corrected between YPD2 and `control_90_110`) — use `analysis/integrated_analysis/phase1_phenotype/genome_strain_species_busco_map.csv` (built from BFD directly) as the authoritative genome-side strain/species list rather than the older, now-stale `BFD/strain_coverage_summary.tsv`.

### Phenotype data provenance
Three timepoint-window rebuilds of the control-media phenotype
(`control_70_80`, `control_80_90`, `control_90_110`) superseded a single
now-removed `control_late` file from the sibling `Rhodotorula_phenotypes`
project. `control_90_110` (latest imaging window) was chosen as canonical
after comparing all three windows against each other and against YPD2 for
color and colony size — see `analysis/examine_phenotype_calling/RECOMMENDATION.md`
for the full reasoning; in short, it has the best agreement with YPD2 on
every trait and the most complete species coverage. **Not yet formally
ingested into `data/metadata/`** — still read cross-project; ingest via
`mycelium:ingest` if it remains the standing choice.

---

## Strategy: Whole-Panel Color → Compound → Genome Analysis

### Phase 1: Color Phenotype Characterization — DONE

**Goal**: build a clean, per-strain and per-species color phenotype table,
usable as either a predictor or response variable, and characterize its
distribution and phylogenetic structure across the whole panel (not framed
around any one species).

**Implemented** (`analysis/scripts/build_strain_phenotype_table.py` →
`build_species_level_tables.py` → `prune_species_tree.R` →
`phylogenetic_signal.R`, outputs in `analysis/integrated_analysis/phase1_phenotype/`):
- Strain-level and species-level tables of L\*, a\*, b\*, C\* (chroma),
  h° (hue angle), and an orange_score composite (z-sum of -L\*, a\*, b\*).
- Species-level phylogenetic signal (Blomberg's K, Pagel's λ): a\* and C\*
  show significant species-level signal (K p≈0.02-0.04); the orange_score
  composite does not clearly (K p≈0.09-0.15 depending on source) at n=17
  species — a useful caution that a single composite score may wash out
  signal present in its individual components (a\* in particular).
- `convergent_color_test.R`: ranks species by orange_score and flags
  phylogenetically-distant, above-average candidates. Retained as a
  diagnostic (useful for understanding phylogenetic structure of the trait)
  but no longer a hard *gate* the way it was under the old framing, since
  the new analysis doesn't depend on any one species being the unique
  outlier.

**Carried forward, still needed**:
- **Species-level collapse is still the right default unit** for any
  phylogenetically-corrected test (PGLS/PIC/Blomberg's K/Pagel's λ) — see
  "Species-Level Collapse: Procedure" below, unchanged from before.
- Whether to use a\*, C\*, or the orange_score composite as *the* color
  predictor for Phase 2 onward is an open call — given the composite's weak
  phylogenetic signal, consider running Phase 2 on **a\* and C\*
  individually as well as** (not only) orange_score, and letting the data
  say which is more predictive of pigment chemistry rather than committing
  to the composite by default.

---

### Phase 2: Color ↔ Metabolome Association (whole panel)

**Goal**: for every MS2 feature (not just the 790 SIRIUS-annotated ones),
test whether its abundance across strains correlates with color — using
the *entire* strain panel as the sample, not a species contrast.

**Steps**:
0. **Feature de-duplication (adducts/isotopologues/in-source fragments) —
   DONE, implemented 2026-08-15 (`analysis/scripts/dedupe_ms_features.py`).**
   16,332 MS2 features are not 16,332 independent chemical entities: the
   same underlying compound routinely produces multiple aligned features
   (different adducts — [M+H]+, [M+Na]+, [2M+H]+ — and isotopologues —
   M+1, M+2 — plus in-source fragments, ions that are decomposition
   products of a co-eluting parent, not independent compounds at all).
   BH-FDR assumes roughly independent tests; without this step, a compound
   that happens to be color-associated would otherwise contribute several
   "independent" hits instead of one, both inflating the apparent number
   of significant features and double-counting evidence in Phase 3's
   enrichment test.

   **EverythingBagel (EB) already computes this grouping** during feature
   finding — it does not need to be re-derived from scratch via RT/mass-
   difference heuristics. The fuller per-run EB output
   (`nf_output/feature_finding/feature_finding_results/aligned_features_ms2.csv`,
   53,040 total detected features before the has-MS2 filter down to this
   project's working 16,332) carries `isotope_source_id`,
   `adduct_source_id`, `is_default_adduct`, `is_isf`, `isf_parent_id`
   columns that this project's simplified working matrix
   (`analysis/linked_data/feature_abundance_matrix.csv.gz`) had dropped.
   `dedupe_ms_features.py` recovers these from the fuller EB output and
   applies them to the working 16,332-feature set: ISF features are folded
   into their (possibly transitive) non-ISF parent's group, then features
   are grouped by `adduct_source_id` (EB's own isotopologue/same-adduct-
   type grouping), with the `is_default_adduct==True` row (else the
   highest-`total_scans` row) as each group's representative.

   **Result: 16,332 raw features → 10,949 deduplicated groups (33%
   reduction)**, 1,774 features flagged as in-source fragments and folded
   into their parent's group. Output:
   `analysis/linked_data/ms_feature_dedup_groups.csv` (row ID →
   dedup_group_id, group_size, is_group_representative, is_isf_member) +
   `ms_feature_dedup_summary.txt`. **10,949, not 16,332, is the real Phase
   2 test count for BH-FDR purposes** — merge `ms_feature_dedup_groups.csv`
   onto the color-correlation results and correct on the deduplicated
   group count, either testing only `is_group_representative==True` rows
   or aggregating (e.g. mean/max abundance) within each group first.

   **Known limitation, not yet applied**: cross-adduct-type compound
   merging (the same compound detected as, separately, [M+H]+ and [M+Na]+
   — the same underlying molecule but a *different* `adduct_source_id`,
   since that column only groups isotopologues within one adduct type).
   EB's `aligned_compounds.tsv` has this information (`n_adducts` per
   compound, ~4.9% of all EB-detected compounds span >1 adduct type) but
   only as a free-text `members` field keyed by m/z, not row IDs — cross-
   referencing it would need a separate, noisier m/z-matching join, not
   attempted here. This means the 10,949-group count is conservative
   (slightly too high, never too low) relative to the true number of
   independent chemical entities — acceptable for FDR purposes (erring
   toward more tests, not fewer, is the safe direction) but worth revisiting
   if a specific compound's cross-adduct identity matters later (e.g. a
   Phase 3 carotenoid candidate that might have both [M+H]+ and [M+Na]+
   forms split across two groups here).
   Report how many of the 16,332 raw features collapse into how many
   deduplicated compound groups — this number, not 16,332, is the real
   Phase 2 test count for FDR purposes.
1. **Continuous color–feature correlation**: for each **deduplicated**
   feature/compound group, Spearman correlation of TSS-normalized
   abundance vs. color, across all strains with both color and MS data
   (~270-300 depending on overlap). BH-FDR correction on the deduplicated
   test count. Stratify by fraction (cell pellet vs. supernatant) —
   pigments are more likely cell-pellet-enriched (intracellular), so this
   stratification itself is informative, not just a nuisance split.
   **Predictor choice (PI-confirmed, 2026-08-15 grilling session): a\* is
   the sole pre-registered PRIMARY predictor** (most direct,
   mechanistically-motivated axis for carotenoid-driven color — the
   green-red CIELAB axis). **C\* is secondary** (run and reported
   alongside, not part of the primary/exploratory split).
   **orange_score is exploratory-only**, not even secondary — Phase 1
   found it has weaker species-level phylogenetic signal than a\*/C\*
   individually (K p≈0.09-0.15 vs. p≈0.02-0.04), and pre-registering three
   co-equal "primary" predictors would just move the multiple-comparisons
   problem up a level rather than resolve it. Report exploratory-tier
   results as exploratory, not folded into the headline result.
2. **Phylogenetic correction**: species-level random effect (mixed model:
   `feature_abundance ~ color + (1|species)`) as the primary test, plus
   `block_permutation.py` (already implemented) permuting within
   phylogenetic blocks, both at species-level and as strain-level
   sensitivity checks per the "Species-Level Collapse" procedure. With
   only 17-18 species-level groups, the `(1|species)` variance component
   in a mixed model is poorly estimated by asymptotic ML/REML — **treat
   the block-permutation empirical p-value as primary**, the mixed-model
   p-value as a secondary/diagnostic check, not the reverse. Also add a
   **strain-level random effect** for the pellet/supernatant pairing
   within a strain (flagged as a gap in the original review, still
   applies here) to avoid pseudoreplicating paired samples.
3. **Leave-one-species-out / leave-clade-out sensitivity check** (Fable
   review, 2026-08-15): for the top color-associated features, recompute
   the correlation/mixed-model result after dropping each of the top-3
   orange species (*R. taiwanensis*, *R. sphaerocarpa*, *R. glutinis*, per
   Phase 1's ranking) one at a time, and again after dropping all three
   together. A "whole-panel" result that disappears when 1-3 species are
   removed is substantially a small-clade story wearing a panel-wide
   label — report this explicitly for every headline hit carried into
   Phase 3, not just run once and left implicit.
4. **Negative control — concrete design, and a HARD PIPELINE GATE**
   (Fable review; enforcement mechanism PI-confirmed, 2026-08-15 grilling
   session: this is not a documented-but-manual step — the Phase 2 script
   must refuse to write a "final" results file unless both null runs below
   have already completed and their output is newer than the input data).
   Two complementary nulls, both run through the *exact* Phase 2 pipeline
   (dedup → correlation → phylogenetic correction) before trusting real
   hits:
   (a) **Label permutation within phylogenetic blocks**: shuffle each
   strain's color value among strains in the same block-permutation clade
   (reuses `block_permutation.py`'s block structure), rerun the full
   feature-correlation scan, repeat ≥200x, and report the empirical
   distribution of "number of BH-FDR<0.05 hits" under the null — this is
   the empirical false-positive-rate calibration.
   (b) **Matched phylogenetically-structured decoy trait**: run the same
   pipeline using colony area (already in the phenotype table, not
   causally tied to pigment chemistry) in place of color, and compare its
   hit count/identity to the real color scan's. If area "predicts" a
   similar number/identity of MS2 features as color does, the color scan
   isn't demonstrating anything color-specific — it's just detecting
   general phylogenetic structure in the metabolome, which any
   phylogenetically-patterned strain trait would also pick up.

**Status: IMPLEMENTED AND RUN, 2026-08-15** —
`analysis/scripts/phase2_color_metabolome_association.py`, full results in
`analysis/integrated_analysis/phase2_metabolome_phenotype/PHASE2_SUMMARY.md`.
**Result: null at the pre-registered threshold** — 0/10,164 (cell) and
0/10,416 (supernatant) deduplicated compound groups reach BH-FDR<0.05 for
either a\* or C\*, matching the permutation-null expectation almost
exactly (well-calibrated, not a bug — see `.living/findings/phenotype-metabolome-association-statistical-power.md`).
Phase 3 therefore has no formally-significant color-associated feature
list to enrichment-test as originally scoped; see PHASE2_SUMMARY.md's
"Recommended next steps" for options (explicit nominal/exploratory ranking,
within-species tests, or waiting for more data). **Important side finding**:
the negative-control decoy run (colony area) uncovered a broad,
previously-unknown confound — colony size associates with 1,524/10,164
cell-fraction (not supernatant) compound groups at BH-FDR<0.05 (permutation
null ~0 expected) — see `.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md`.
This is exactly the kind of result the hard-gated negative control (PI-
confirmed, grilling session) was designed to catch before trusting a real
scan; it was.

**Output**: a ranked list of deduplicated MS2 compound groups by
color-association strength, annotated with which top-species removals (if
any) it survives — this is the direct input to Phase 3.

---

### Phase 3: Compound-Class Enrichment (hypothesis test)

**Goal**: directly test the carotenoid hypothesis — are terpenoid/carotenoid-
class compounds enriched among the color-associated features from Phase 2,
more than expected by chance?

**Caveat that governs how to read this whole phase** (Fable review,
2026-08-15): only 790/16,332 (4.8%) of MS2 features are SIRIUS-annotated at
all, and that 790-feature background is not a random sample of the
metabolome — SIRIUS annotates what matches spectral libraries well, which
skews toward common/well-characterized compound classes. An enrichment
test against this background answers "among features SIRIUS could
confidently call, are terpenoids overrepresented in color-associated ones,"
which is narrower and more fragile than "is carotenoid chemistry driving
color" — do not let Phase 3's result be read as the latter. With only 3-4
carotenoid-class calls total in the whole annotated set, **the enrichment
test as scoped has essentially no power to detect anything short of near-
total concentration of those specific hits in the color-associated list**
— compute and report the minimum detectable enrichment (e.g. via a
hypergeometric power calculation for this exact background size and hit
count) *before* running the test, so a null result isn't misread as
"carotenoids ruled out" when it may just mean "underpowered."

**Steps**:
1. Join Phase 2's deduplicated, color-associated feature list to
   `sirius_annotations.tsv`'s NPC pathway/class calls.
2. Enrichment test (Fisher's exact or hypergeometric): is the Terpenoids
   NPC pathway (and specifically the Carotenoids/Apocarotenoids NPC
   classes) over-represented among color-associated features vs. the
   background rate in all SIRIUS-annotated features (73/790 ≈ 9.2%
   Terpenoids background rate)? Report alongside the minimum-detectable-
   enrichment power calculation above.
3. Repeat for other candidate pigment-related classes flagged in the
   original motivation (Amino acids and Peptides pathway — cyclic
   peptides/diketopiperazines — is the single largest NPC pathway in this
   dataset at 359/790 features, so it's worth checking whether *it*, not
   terpenoids, is what actually tracks color; this would be a genuine,
   reportable surprise rather than a hypothesis failure).
4. **Report both directions honestly**: if terpenoids/carotenoids are NOT
   enriched among color-associated features, that's a real, useful result
   (rules out the straightforward hypothesis, motivates targeted carotenoid
   MS/MS or LC-UV/Vis work instead of relying on untargeted SIRIUS calls)
   — but see the power caveat above before drawing that conclusion.
5. **Compound-ID methods beyond SIRIUS structure calls** (Fable review;
   the PI is adding data for this now — SIRIUS is re-running with
   improvements as of 2026-08-15). Recommended, in rough priority order:
   - **MS/MS spectral networking** (GNPS-style cosine-similarity molecular
     networking) across all 16,332 (or deduplicated) features, independent
     of SIRIUS structure ID. Carotenoids/apocarotenoids fragment in
     recognizable, related patterns even without a confident library
     match — an unannotated feature cluster that co-clusters with the 3-4
     SIRIUS carotenoid hits is much stronger evidence than the SIRIUS
     calls alone, and would surface color-relevant compound families
     SIRIUS missed entirely (it annotated only 4.8% of features).
   - **Retention-time/mass matching against literature values** for known
     *Rhodotorula* carotenoids (torulene, torularhodin, β-carotene,
     γ-carotene all have published RT/MS characteristics under comparable
     LC-MS conditions) — a targeted look-up against the full feature list,
     not reliant on untargeted SIRIUS calls at all. Cheap to run once
     reference values are compiled; should be done regardless of what
     SIRIUS's re-run produces.
   - **Isotope-pattern / characteristic-loss screening**: carotenoids have
     distinctive polyene chromophores with recognizable neutral-loss and
     isotope-pattern signatures separate from general SIRIUS scoring —
     worth a dedicated, targeted screen of the full feature list rather
     than relying only on whatever SIRIUS's confidence threshold passed.

**Key question**: does the color-metabolome association data support
carotenoids (vs. some other compound class) as the chemistry underlying
orange/yellow/red coloration in this panel?

---

### Phase 4: Classifier — Predict Pigment Compound Class from Color

**Goal**: build and validate a model that predicts pigment-related compound
presence/abundance from color phenotype alone — the actual "classifier"
requested.

**Steps**:
1. **Target definition (PI-confirmed, 2026-08-15 grilling session):
   CONTINUOUS**, not binary — the summed/aggregated abundance of whatever
   compound group Phase 3 flags as color-associated (not restricted in
   advance to the 3-4 SIRIUS carotenoid-class hits, since Phase 3 may find
   a different class is actually color-associated). Avoids an arbitrary
   presence/absence threshold given only 3-4 confidently-annotated
   carotenoid features to calibrate one against, and preserves more signal
   at the already-limited ~17-18 species-level effective sample size.
   Thresholded binary accuracy may still be reported as a secondary,
   interpretability-only metric, not the primary target.
2. **Model**: start simple — logistic/linear regression on L\*/a\*/b\*/C\*
   as predictors of the target, before reaching for Random Forest/XGBoost;
   given the panel is ~300 strains (much better than the old single-clade
   n), a tree ensemble is defensible here but should be benchmarked against
   the simple model, not assumed better.
3. **Validation**: **leave-one-species-out cross-validation** (LOSO-CV) is
   still the right choice — it tests whether the color→compound
   relationship generalizes across species boundaries, not just within one
   species' strains. Report LOSO-CV performance alongside a naive
   random-fold CV to show how much of the apparent performance is
   phylogeny-driven.
   **Important caveat (Fable review, 2026-08-15): LOSO-CV alone does not
   rule out the same circularity flagged in the original plan's Phase 4 ML
   stack, it just relocates it.** If color and the compound-class target
   are both phylogenetically structured (they are — Phase 1 found
   significant species-level signal in a\*/C\*), a classifier can hit good
   LOSO accuracy purely by learning "high-a\*/C\* species are in the
   terpenoid-high clade" — i.e. recovering phylogeny *through* color,
   without any causal color→chemistry link. LOSO-CV protects against
   memorizing individual strains, not against this. **Required
   diagnostic**: benchmark LOSO-CV performance against (a) the same
   phylogenetically-structured-decoy-trait null used in Phase 2 (does
   colony area predict the same compound-class target about as well as
   color does?), and (b) a within-tree label-permutation null (shuffle
   color among strains within the same phylogenetic block, refit, repeat
   ≥100x, compare real LOSO accuracy to the null distribution). Only
   report the classifier as demonstrating a color-specific relationship if
   it clears both — beating naive random-fold CV is not sufficient on its
   own.
4. **Feature importance**: which color axis (L\*, a\*, or b\*/C\*) is most
   predictive — directly informs whether "darker" (L\*), "redder" (a\*), or
   "more chroma" (C\*) is the operative signal, refining the orange_score
   composite definition itself if useful.

**Key question**: can color phenotype alone predict pigment compound class
presence/abundance well enough to be useful as a low-cost proxy/classifier?

---

### Phase 5: Genome Linkage

**Goal**: for whichever compound(s)/classes Phase 3-4 establish as
color-predictive, test genomic correlates across the whole panel.

**Steps**:
1. **Candidate pathway genes first**: known carotenoid biosynthesis Pfam
   domains (phytoene synthase / CrtB, phytoene desaturase / CrtI, lycopene
   cyclase / CrtY-type, and fungal-specific carotenogenesis genes) —
   search presence/copy-number/sequence variation across all 276 genomes.
   This is the pre-registered, primary hypothesis test (small number of
   genes, decided in advance, not a fishing expedition).
2. **Broader genome–compound association**: Spearman/PGLS of Pfam/CAZy/
   MEROPS feature counts (proteome-size-normalized, per
   "Species-Level Collapse" below) against the color-predictive compound's
   abundance, species-level primary + strain-level sensitivity check, BH-FDR.
   This part **is** exploratory (large number of tests) and should be
   reported as such, separate from the pre-registered pathway-gene test.
3. **Negative-control calibration**, same rationale as Phase 2/old-Phase-3.

**Key question**: are known or novel genomic features associated with the
compounds identified as color-predictive?

---

### Phase 6: Mechanistic Follow-Up (deferred, gated on Phase 5)

Same content as the old Phase 5 (pathway enrichment, antiSMASH/BGC search,
OrthoFinder/CAFE5 gene-family evolution, targeted metabolite/genetic
validation) — retained as the final stage, gated on Phase 5 producing a
short candidate gene/compound list, not scoped in detail until then. See
Appendix for the original text (still applicable, just later in sequence
and no longer *R. dairenensis*-specific — apply to whichever
species/strains carry the color-predictive compounds).

**Exception, moved earlier (Fable review, 2026-08-15)**: a lightweight
**targeted HPLC/LC-UV-Vis pigment quantification** on a representative
strain subset spanning the color range would ideally run in parallel with
Phase 2-3, not deferred here. So much of Phases 3-5's interpretation rests
on compound identity being right, and SIRIUS is known to be unreliable
specifically for carotenoids (class-level calls only, 3-4 hits total) — a
modest wet-lab check would be cheap insurance against building Phases 4-5
on a shaky compound call.

**Status (PI-confirmed, 2026-08-15 grilling session)**: this is **not on
a fixed timeline** — it's contingent on what the currently-in-progress
improved SIRIUS re-run shows first, not scheduled independently. Phases
2-4 proceed now on SIRIUS + spectral-networking evidence without waiting
for it. Consequence: **any Phase 3-5 compound-identity claim must be
explicitly flagged as provisional** (e.g. "provisionally identified as
carotenoid-class, pending confirmation") in every write-up until/unless
this validation happens — not stated as confirmed.

---

## Species-Level Collapse: Procedure

*(Unchanged from the original plan — still the correct approach for any
phylogenetically-corrected test in Phases 2, 3, or 5 above.)*

### Rationale
Phylogenetic comparative methods (PGLS, PIC, Blomberg's K, Pagel's λ) model
trait evolution along the branches of a species tree, where each **tip
represents one evolutionary lineage**. When multiple strains of the same
species are entered as separate tips (or as replicate rows against a
single species-tree tip), the within-species covariance among those
strains is not part of the Brownian-motion/OU model being fit, and the
effective degrees of freedom used by the test are silently inflated.
Collapsing to one row per species (one per tip in the PHYling species tree)
aligns the unit of analysis with the unit the statistical model assumes.

### Step 1 — Build species-level phenotype table
Implemented: `analysis/scripts/build_species_level_tables.py`. Group
strain-level rows by species; for each trait compute central tendency
(mean/median) + within-species SD/n. Flag species with n_strains = 1 (no
within-species variance estimate).

### Step 2 — Build species-level genome/metabolome feature matrices
For count features (Pfam/CAZy/MEROPS/codon/AA-frequency/MS2 abundance):
normalize per strain first (e.g. by total protein count for genome counts,
TSS for MS abundance), *then* aggregate to species mean/median. For
binary/presence features: species-level **prevalence**, not a hard
collapse. For assembly-quality covariates (BUSCO completeness, N50): species
mean, retained as PGLS covariates, not just predictors.

### Step 3 — Prune/collapse the phylogeny to species-level tips
Implemented: `analysis/scripts/prune_species_tree.R`. Picks one
representative strain per species (best BUSCO completeness via
`--quality-col`) and drops the rest; `ape::drop.tip`'s default
edge-collapsing correctly attaches the representative at the species' stem
branch length, not its own arbitrary depth. Output:
`analysis/integrated_analysis/phase1_phenotype/species_tree.nwk`.

### Step 4 — Run comparative tests at the species level (primary)
Merge species-level phenotype/genome/metabolome tables on `species`,
matched to tree tips. Include assembly-quality/proteome-size covariates in
the model formula. Run the negative-control calibration pass through the
same pipeline before trusting real-trait results.

### Step 5 — Strain-level results as secondary robustness check
Re-run naive (non-phylogenetic) tests at the strain level as a
sensitivity/robustness comparison, not the primary result. Species-level +
strain-level agreement is the strongest evidence; strain-level-only
significance should be flagged as likely pseudoreplication-driven.

### Step 6 — Within-species variation tests remain strain-level
Any question that is explicitly about within-species variation (e.g. do
the 8-9 *R. dairenensis* strains' internal color variation track their
internal compound variation?) stays at full strain resolution — this is
the one part of the analysis that is genuinely phylogeny-free at the tip
level, and remains useful as a complementary check wherever a species has
enough strains (n≥5ish) to support it, not just for *R. dairenensis*.

### Software note
`ape::drop.tip`, `phytools`, `caper::comparative.data` (R) for pruning/PGLS;
`build_species_level_tables.py` / `prune_species_tree.R` already implemented
and generic (accept any strain-level input table via CLI flags), reusable
across Phases 2, 3, and 5 without modification.

---

## Key Considerations (carried forward)

### Species-level power ceiling (Fable review, 2026-08-15)
State this plainly wherever the plan is summarized externally: **the
whole-panel reframing raises statistical power relative to the old
single-clade contrast, but the hard ceiling for any phylogenetically
corrected (species-level) test in this project is n≈17-18 independent
species lineages** — not the ~300 strains or 16,332 MS2 features headline
numbers suggest. Any trait with strong phylogenetic signal (a\*/C\* already
qualify) will in practice be shaped by a handful of species at the
extremes (Phase 1: *R. taiwanensis*, *R. sphaerocarpa*, *R. glutinis*), so
"whole-panel, n=17" is a real improvement over "n=1 clade" but is not
high-powered in an absolute sense — don't let downstream write-ups imply
otherwise. Every headline Phase 2-5 result should report a rough minimum-
detectable-effect-size or the leave-one/few-species-out sensitivity result
(see Phase 2 step 3) alongside its p-value/FDR, not the p-value alone.

### Multiple Testing
Phase 2 (16,332 raw MS2 features → **10,949 deduplicated groups**, see
Phase 2 step 0) and Phase 5's exploratory arm (~5,000+ Pfam domains, etc.)
are both large test batteries — BH-FDR < 0.05 throughout, on the
**10,949-group** count for Phase 2, with effect sizes reported alongside
FDR (an FDR-significant hit at this sample size
can still be a small, uninteresting effect).

### Pre-registration discipline
Given the number of features/tests across Phases 2-5, keep the same
primary/exploratory split discipline flagged in the original review:
pathway/class hypotheses decided in advance (carotenoid pathway genes,
Terpenoids NPC class) are primary; everything else is exploratory and
reported as such.

### Compound-ID annotation coverage (Fable review, 2026-08-15)
Only 4.8% (790/16,332) of MS2 features currently have any SIRIUS
annotation, and the carotenoid-specific calls within that set are
class-level, not structure-confirmed (see Phase 3). Any enrichment or
association result that depends on SIRIUS's NPC pathway/class calls
inherits this coverage/detection bias — treat Phase 3's result as
conditional on "among features SIRIUS could annotate," not as a claim
about the whole metabolome, until the spectral-networking/RT-matching
methods in Phase 3 step 5 (or the PI's improved SIRIUS re-run, in progress
as of 2026-08-15) materially raise annotation coverage.

### Within-species restriction as the most diagnostic check (added 2026-08-15)
Phase 2's results (color vs. metabolome — null at whole-panel,
within-species x3, and multivariate scales; see PHASE2_SUMMARY.md,
WITHIN_SPECIES_MUCILAGINOSA.md, ROBUSTNESS_AND_MULTIVARIATE.md) turn out
to have a direct precedent in this project: the earlier, separate
`analysis/copper/` analysis (copper-resistance growth rate/AUC vs.
proteome-wide amino acid composition) found 5 amino acids (S, L, Q, W, T)
significant by naive whole-panel correlation, and **PGLS — whole-panel
phylogenetic correction — still called several of these significant**
(`analysis/copper/outputs/naive_vs_pgls_comparison.csv`). But a
within-*R. mucilaginosa*-only sensitivity check
(`sensitivity_mucilaginosa_only.csv`) found **none of them held up**
(p = 0.20-0.68 for all 6 tested). See
`.living/findings/phylogenetic-confounding-of-trait-molecular-associations.md`
for the full cross-cutting writeup.

**Two independent phenotype/molecular-layer pairings in this project now
show the same pattern**: whole-panel "hits" — even PGLS-corrected ones —
have not survived within-species restriction. This makes within-species
testing the single most diagnostic check available in this project so
far, more so than whole-panel PGLS alone: **treat any future whole-panel
"hit" (Phase 3 enrichment, Phase 5 genome linkage, or a revisited
copper-AUC analysis) as provisional until checked against a within-species
restriction.**

**Caveat (PI, 2026-08-15): the copper-AUC result above is not a closed
finding** — the PI plans to revisit that phenotype's underlying data
further. The pattern it currently illustrates (naive/PGLS "hit" not
surviving within-species testing) is worth keeping as a methodological
precedent regardless, but the specific amino acids/numbers should be
treated as provisional pending that revisit, not as a settled result to
build on.

### Biological Context
*Rhodotorula* species produce carotenoid pigments (torulene, torularhodin,
β-carotene) broadly associated with orange-red coloration. Candidate
explanations for color variation across the panel, all testable within
this plan once Phase 2-3 identify which compounds actually track color:
1. Differential carotenoid production (pathway gene expression/copy number)
2. Different carotenoid composition (torulene vs. β-carotene ratio)
3. Lipid droplet storage capacity (more/larger droplets concentrating pigment)
4. Degradation pathway differences
5. Regulatory differences (transcription factors controlling biosynthesis)
6. A chemically distinct, non-carotenoid red pigment class not yet
   identified in this dataset (open possibility per the working hypothesis
   — Phase 3's enrichment test is exactly what would surface this if true)

---
---

## Appendix: Superseded *R. dairenensis*-Specific Framing (preserved for provenance)

Everything below this line is the original plan and its subsequent
statistical review, written when the project was framed around explaining
why *R. dairenensis* specifically is darker orange. It is superseded by the
whole-panel framing above, but preserved because (a) it contains real,
reusable methodological analysis (the "n=1 clade" problem, phylogenetic
signal results, phenotype-source comparison work) that informed the design
above, and (b) *R. dairenensis* remains a strain group of interest within
the new analysis even though it's no longer the sole justification for it.

### Original Motivation
*R. dairenensis* strains (n=10 phenotyped, n=9 with full genome + MS data)
exhibit a **darker orange color** that separates them from other
*Rhodotorula* species in both the CIELAB color phenotype ordination (PCoA
of L\*, a\*, b\*) and the Bray–Curtis PCoA of MS2 metabolome features. The
goal is to identify **genomic features** (gene family expansions/
contractions, functional domain gains/losses, pathway-level differences)
that correlate with—and potentially drive—this phenotypic difference, by
building an integrated framework that connects genome → metabolome →
phenotype.

**Note (2026-08-15)**: this framing's premise was checked in Phase 1 and
did not hold up as stated — *R. dairenensis* was found to be only the 4th
of 17-18 species on a darker-orange composite score, not the outlier. The
original PCoA cited above may show something different than the composite
score (that specific comparison was never completed) — worth a quick check
if anyone wants to understand exactly why the original framing arose, but
it does not block the whole-panel plan above, which doesn't depend on the
answer.

### Original 5-Phase Plan (Phases 2-5 as originally scoped)

#### Phase 2: Metabolome–Phenotype Association (Extend Existing) [original]
**Goal**: Identify which MS2 metabolomic features track with the darker
orange color phenotype **in R. dairenensis specifically**.
- Differential features by species pair (`scripts/differential_features_by_species.py`)
  already identifies MS features differentially abundant between
  *R. dairenensis* and other species; SIRIUS annotations provide compound
  identity for some features.
- Continuous phenotype–metabolome correlation (Spearman + BH-FDR,
  stratified by fraction); MWAS-color mixed model
  (`feature_abundance ~ orange_score + species_random_effect`) with block
  permutation; compound class enrichment via SIRIUS NPC pathway/class
  (terpenoids/carotenoid pathway, fatty acids/lipid droplets, peptides/
  cyclic peptides were the classes flagged as worth watching).

*(This phase's methodology — correlation, mixed model, block permutation,
compound class enrichment — is essentially identical to the new Phase 2-3
above; only the framing, "in R. dairenensis" vs. "across the whole panel,"
differs. The new phases are a direct generalization of this one.)*

#### Phase 3: Genome–Phenotype Association [original]
**Goal**: Identify genomic features that correlate with the darker orange
color phenotype **in R. dairenensis vs. the rest of the tree**.
- Build strain × functional feature matrices from BFD (Pfam ~5,000 IDs,
  MEROPS ~300, CAZy ~100, SignalP, TargetP, TMHMM, WolfPSort, PredGPI, IDP,
  codon usage, AA frequency, genome stats) for all 276 strains.
- Per-feature association tests (Fisher's/Mann-Whitney for binary, Spearman
  for count features), BH-FDR < 0.05.
- Phylogenetic correction: PGLS, PIC, block permutation; naive vs.
  phylogenetically-corrected comparison.

*(Superseded by the new Phase 5, which runs the same feature-matrix/
association-test machinery against whichever compounds Phase 2-4 establish
as color-predictive, across the whole panel rather than one clade.)*

#### Phase 4: Integrated Genome–Metabolome–Phenotype Model [original]
**Goal**: multi-omics model (sPLS/MOFA, Random Forest/XGBoost + SHAP,
LOSO-CV, bipartite network analysis) predicting *R. dairenensis*'s color
phenotype from combined genome + metabolome features.

*(Superseded by the new Phase 4's simpler, panel-wide classifier — the
original's sPLS/MOFA/network-analysis stack was flagged in the review below
as high risk of overfitting/circularity at single-clade sample sizes; the
whole-panel version starting from a simple model per new-Phase-4 step 2 is
the more defensible version of the same idea.)*

#### Phase 5: Mechanistic Hypothesis Generation & Validation [original]
Pathway-level analysis (KEGG/GO enrichment; carotenoid biosynthesis, lipid
metabolism, ROS/stress response, secreted peptide pathways flagged as
priorities), BGC identification (antiSMASH) comparing *R. dairenensis* to
close relatives, OrthoFinder/CAFE5 gene-family evolution analysis, and a
targeted validation strategy (metabolite structure confirmation, genetic
knockout/complementation, comparative validation against convergent
darker-orange lineages if found).

*(Content carried forward essentially unchanged as the new Phase 6, applied
to whichever species/strains turn out to carry the color-predictive
compounds rather than assumed to be *R. dairenensis*.)*

### Original Implementation Plan scaffolding
```
Software Dependencies:
  Python: pandas, numpy, scipy, scikit-learn, statsmodels, matplotlib, seaborn
  R (phylogenetic models): ape, phytools, phylolm, caper
  Multi-omics: mixOmics (R) or pyMOFA (Python) for sPLS/MOFA
  Network: networkx, pyvis
  SHAP: shap library for feature importance

Output structure (original, R. dairenensis-framed; superseded by whatever
directory layout the new Phases 2-6 actually produce, but the phase-by-phase
script/output naming convention below is still a reasonable template):
analysis/integrated_analysis/
├── phase1_phenotype/        # done, still valid, see above
├── phase2_metabolome_phenotype/
│   ├── metabolome_color_correlation.csv
│   ├── mwas_color_results.csv
│   └── compound_class_enrichment.csv
├── phase3_genome_phenotype/
│   ├── genome_feature_matrices/
│   ├── genome_color_association.csv
│   ├── pgls_results.csv
│   └── pathway_enrichment.csv
├── phase4_integrated_model/
│   ├── spls_results/
│   ├── ml_prediction_results.csv
│   ├── shap_feature_importance.csv
│   └── network_analysis/
└── phase5_mechanistic/
    ├── hypothesis_table.csv
    ├── gene_cluster_comparison.csv
    └── validation_strategy.md
```

### Statistical Review & Revisions (external review by Fable, 2026-08-15)

An independent statistical/scientific review of the original *R.
dairenensis*-specific plan raised the following, most of which motivated
the reframing above and/or still applies to the whole-panel version:

**The "n=1 clade" problem** (the core finding that motivated reframing):
*R. dairenensis* is a single monophyletic clade. PGLS, PIC, and block
permutation correct for expected covariance under a phylogenetic model of
trait evolution, but they do not manufacture independent replicates of a
single clade's color/genome co-evolution event — no statistical correction
can distinguish "this feature drives color" from "this feature and color
coincidentally changed on the same branch" without a second, independent
origin of the trait. This is the reason the project no longer treats
*R. dairenensis* alone as sufficient grounds for a genome-association
study; the whole-panel reframing sidesteps it by using the full diversity
of independent species/lineages as the sample instead of one clade vs. the
rest.

**Statistical power**: naive (non-phylogenetic) tests on a clade contrast
show many nominal hits because the clade is a strong, cohesive signal —
largely illusory for causal purposes. Phylogenetically corrected, the
effective sample size collapses toward the number of independent
trait-origin events (as few as 1 for a single clade). This is the
quantitative version of the n=1 problem above, and part of why the
whole-panel version (mixed models / PGLS across 17-18 independent species
lineages, not 1 clade vs. rest) is preferred.

**Additional gaps identified** (all still apply to the new Phases 2-5):
- Technical/batch covariates (assembly quality, proteome size, sequencing
  batch, isolation source/geography) should be covariates in genome
  association models, not just predictors.
- No pre-registered primary vs. exploratory hypothesis split across the
  large number of feature-matrix types and tests — still needed (see "Key
  Considerations" above).
- No negative-control calibration run (phylogenetically neutral or
  permuted trait through the full pipeline to check empirical false
  positive rate) — still needed, now specified in the new Phase 2/5.
- MS pellet/supernatant pairing within a strain needs a strain-level random
  effect in addition to species-level, to avoid pseudoreplication — still
  needed, now specified in the new Phase 2.
- BGC/OrthoFinder work (old Phase 5, new Phase 6) is premature to scope in
  detail before a filtered candidate list exists from the association
  phases — still true, now gated on the new Phase 5 instead of old Phase 3.
- Effect sizes should be reported alongside FDR, not FDR alone — still true.

**Species-level collapse vs. strain-level analysis**: the recommendation to
use species-level collapse as the primary unit for phylogenetically
corrected tests, with strain-level as a secondary robustness check, is
unchanged and carried forward as-is into the new plan's "Species-Level
Collapse: Procedure" section above.

### Phase 1 implementation history (preserved log, 2026-08-15)

Phase 1 was implemented and run end-to-end under the original framing:
`analysis/scripts/build_strain_phenotype_table.py` → `build_species_level_tables.py`
→ `prune_species_tree.R` → `phylogenetic_signal.R` → `convergent_color_test.R`,
outputs in `analysis/integrated_analysis/phase1_phenotype/`.

**Result that motivated the reframing**: on the original YPD2 phenotype
data (316 of 318 strains usable), *R. dairenensis* was only the
**4th-highest of 17 species** on the orange_score composite, behind
*R. taiwanensis* (~7x higher), *R. sphaerocarpa*, and *R. glutinis*. This
held regardless of genus-filtering (ruling out an outgroup-composition
artifact). `phenotype_phylogenetic_signal.csv` showed a\*/C\* with
significant species-level phylogenetic signal (K p=0.04, p=0.035) but the
orange_score composite itself did not clearly (K p=0.087) at n=17 species.
`convergent_color_candidates.csv` flagged *R. glutinis*, *R. sp. clade I*,
and *R. diobovata* as phylogenetically-distant, above-average-orange_score
candidates; *R. taiwanensis* and *R. sphaerocarpa* scored higher still but
were too closely related to *R. dairenensis* on the tree to count as
independent evidence under that filter.

**Cross-checked against an independently rebuilt, more complete phenotype
table** (`control_late`, later split into the three `control_*` timepoint
windows described in "Phenotype data provenance" above): reproduced the
same ranking (*R. dairenensis* still 4th of 17-18 species) — not an
artifact of YPD2's missing-species rows or outgroup composition.

This result is what prompted the 2026-08-15 reframing at the top of this
document. It does not need to be re-litigated to proceed with the
whole-panel plan above, but is preserved here since it's the empirical
trigger for the pivot.
