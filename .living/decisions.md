# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

## [2026-08-15] Retrofit existing repo into mycelium living-repo structure (non-destructive)

**Context**: This repo already had substantial structure (`analysis/`, `data/{raw,processed,metadata}/`, `scripts/`, `BFD/`, `nextflow/`, `AGENTS.md`, `GOALS.md`) before mycelium init was requested. `init_repo.py --restructure` is an unimplemented stub (prints an audit, does nothing) as of this mycelium version.

**Decision**: Ran `init_repo.py` without `--restructure` directly against the existing repo root. Verified first that its directory/manifest creation is purely additive (`mkdir(parents=True, exist_ok=True)`, manifest files only written `if not exists`) — it does not move, rename, or overwrite any pre-existing file. Installed core convention packs (`robust-analysis`, `report-generator`, `idea-generator`) plus the `bioinformatics` domain pack (genomics/phylogenetics work throughout this project). Did not install `image-analysis` (colony imaging happens upstream in the sibling `Rhodotorula_phenotypes` repo, not here) or `skill-bridge` (no external skillpacks cloned). Kept `AGENTS.md`/`GOALS.md` as the primary project-context files; `CLAUDE.md` (generated from the mycelium template) points to them rather than duplicating their content.

**Alternatives considered**:
- Wait for `--restructure` to be implemented — rejected, it's a stub with no timeline, and the plain `init` path was already verified non-destructive for this repo's existing layout.
- Manually move existing `analysis/`/`data/` content into a fresh mycelium skeleton — rejected as unnecessary churn; the existing layout already matches mycelium's expected `analysis/`, `data/{raw,processed,metadata}/` structure closely enough that `validate_structure.py` passes without moving anything.
- Auto-clone the three skillpacks reference repos (scientific-agent-skills, bioSkills, Autonomous-Science) — deferred; scaffolded `skillpacks/` with `.gitignore`/`README.md` only, left cloning as a PI decision (network operation with real footprint, not needed for current phases).

**Rationale**: Lowest-risk path to a working living-repo layer without disturbing in-progress analysis work (`analysis/INTEGRATED_ANALYSIS_STRATEGY.md` and its dependents were mid-revision at the time of this init).

**Consequences**: `validate_structure.py` passes with 10 warnings, all "subfolder missing its UPPER_SNAKE_CASE.md doc" (e.g. `analysis/sirius_annotation/` has no `SIRIUS_ANNOTATION.md`) — pre-existing analysis subfolders from before mycelium adoption. Not blocking; backfill opportunistically per the post-action protocol as those folders are touched again, tracked as a todo item rather than done in bulk now. `.claude/settings.local.json` now has mycelium's SessionStart/PostToolUse/Stop hooks registered alongside the one pre-existing permission entry (verified no conflicts before running).

**Tags**: mycelium, repo-init, non-destructive, provenance

## [2026-08-15] Formally ingest control_90_110 phenotype table; leave other sources cross-project

**Context**: `analysis/scripts/build_strain_phenotype_table.py` was reading
the canonical phenotype source (`control_90_110`) directly from a mutable
path in the sibling `Rhodotorula_phenotypes` project — flagged repeatedly
(Fable review, RECOMMENDATION.md) as a reproducibility risk given this
project's history of phenotype-source discrepancies (YPD2 vs. control
windows, strain-code corrections between them).

**Decision**: Copied `phenotype_control_timepoint_90_110.csv` verbatim
into `data/raw/control_phenotype_90_110h/` (SHA256 verified post-copy),
wrote full `schema.yaml`/`provenance.md`/`summary_stats.md` in
`data/metadata/control_phenotype_90_110h/`, added a `DATA_MANIFEST.md`
entry, and repointed `build_strain_phenotype_table.py`'s `control_90_110`
source at the local ingested copy. Left `control_70_80`/`control_80_90`
(the two non-canonical timepoint windows, kept only for the comparison in
`analysis/examine_phenotype_calling/`) and the legacy YPD2 table reading
cross-project / in their pre-existing locations — not ingested, since they
aren't the active analysis source.

**Alternatives considered**:
- Ingest all three control windows — rejected as unnecessary; only
  `control_90_110` is used for active analysis going forward, the other
  two exist solely for the one-time timepoint-comparison writeup.
- Retroactively move the legacy YPD2/`EXFAB_UCR-005` files into the
  `data/raw/` + `data/metadata/` split too — rejected for this pass;
  several existing scripts (`scripts/pcoa_color_phenotype.py`,
  `analysis/copper/scripts/common.py`, etc.) read that path directly by
  its current location, and moving it would require updating every
  reader. Documented as a known pre-mycelium gap in `DATA_MANIFEST.md`
  instead of silently leaving it unexplained.

**Rationale**: Ingest exactly what's load-bearing for ongoing analysis;
don't do a larger data-migration pass than the actual dependency requires.

**Consequences**: Provenance for the *computational* pipeline that
produced this table is documented (reproduced from the upstream project's
own README). Provenance for the underlying *wet-lab* experimental design
(imaging rig, plate randomization, replicate structure) is explicitly
flagged as pending PI input in both
`data/raw/control_phenotype_90_110h/CONTROL_PHENOTYPE_90_110H.md` and
`data/metadata/control_phenotype_90_110h/provenance.md` — this ingestion
is not "complete" by the mycelium checklist until the PI fills those
sections in. Also surfaced (and fixed, see `.living/learnings.md`) a
previously-silent data-loss bug: 10/314 rows in this table have no
`strain_code` and were being silently dropped by `pandas.groupby`.

**Tags**: mycelium, data-ingest, phenotype, provenance, control_90_110

## [2026-08-15] Use EverythingBagel's existing adduct/isotope grouping for MS2 feature de-dup, don't re-derive it

**Context**: Fable's review of the reframed strategy flagged that the
16,332 MS2 features are not independent chemical entities (isotopologues/
adducts/in-source fragments of the same compound get aligned as separate
"features"), inflating Phase 2's effective test count and double-counting
evidence in Phase 3's enrichment test. The strategy doc originally
sketched an ad hoc RT-window + mass-difference heuristic to dedupe.

**Decision**: Before writing that heuristic, audited
`data/processed/EB_20260130_ExFAB_Rhodo_Sup_and_Pellet/.../nf_output/`
and found EverythingBagel's feature-finding step already computes this
grouping (`isotope_source_id`, `adduct_source_id`, `is_default_adduct`,
`is_isf`, `isf_parent_id` columns in the fuller
`aligned_features_ms2.csv`, plus a full `aligned_compounds.tsv` /
molecular-networking output this project's simplified 16,332-row working
matrix had dropped). Wrote `analysis/scripts/dedupe_ms_features.py` to
recover and apply that existing grouping instead of re-deriving one.
Result: 16,332 raw features → 10,949 deduplicated groups (33% reduction),
1,774 in-source-fragment features folded into their parent's group.

**Alternatives considered**:
- Write the originally-sketched RT-window + mass-difference heuristic from
  scratch — rejected once EB's own grouping was found; re-deriving
  something EB already computed more rigorously (actual chromatographic
  co-elution + isotope pattern modeling, not just RT proximity) would be
  strictly worse and duplicative effort.
- Also apply EB's cross-adduct-type compound merging
  (`aligned_compounds.tsv`'s `n_adducts` field, ~4.9% of compounds span
  >1 adduct type) — deferred; that file only provides a free-text
  `members` field keyed by m/z, not row IDs, so merging it in requires a
  separate, noisier m/z-matching join. Left as a documented limitation
  (10,949 is a conservative, i.e. slightly too high, group count) rather
  than risk a fragile parser producing wrong merges.

**Rationale**: Don't re-implement analysis EverythingBagel already did
correctly; recovering existing annotation is more reliable and far
cheaper than re-deriving it heuristically.

**Consequences**: Phase 2's real BH-FDR test count is 10,949, not 16,332
— strategy doc and Multiple Testing section updated accordingly. The
conservative (not-fully-merged) grouping means a handful of compounds with
multiple detected adduct types will still appear as 2+ separate groups;
acceptable for FDR purposes (errs toward more tests, not fewer) but
worth resolving if a specific Phase 3 candidate compound's identity
depends on merging its adduct forms.

**Tags**: mass-spec, deduplication, everythingbagel, adducts, isotopologues, phase-2

## [2026-08-15] Phase 2/4 predictor and target definitions (PI grilled, decisions confirmed)

**Context**: `analysis/INTEGRATED_ANALYSIS_STRATEGY.md`'s Phase 2/4 left two
choices open: which color axis is the primary predictor, and whether
Phase 4's classifier target is binary or continuous. Phase 1 found a\*/C\*
have significant species-level phylogenetic signal (K p≈0.02-0.04) but the
orange_score composite doesn't clearly (K p≈0.09-0.15) — a composite that
washes out signal present in its own components is a weak choice to
pre-register as primary.

**Decision** (PI confirmed, interactive grilling session): 
- **a\* is the sole pre-registered PRIMARY predictor** for Phase 2 (most
  direct, mechanistically-motivated axis for carotenoid-driven color —
  it's literally the green-red CIELAB axis). 
- **C\* is secondary** (run and reported, not the basis for the
  primary/exploratory split). 
- **orange_score is demoted to exploratory-only** — not even secondary,
  since pre-registering three co-equal "primary" predictors just moves the
  multiple-comparisons problem up a level rather than resolving it.
- **Phase 4's classifier target is CONTINUOUS**: summed abundance of
  whatever compound group Phase 3 flags as color-associated (not
  restricted to the 3-4 SIRIUS carotenoid-class hits specifically, since
  Phase 3 may find a different class is actually color-associated).
  Avoids an arbitrary presence/absence threshold given the tiny number of
  confidently-annotated carotenoid features to calibrate one against, and
  preserves more signal at the already-limited ~17-18 species-level
  effective sample size. Thresholded binary accuracy reported as a
  secondary/interpretability metric, not the primary target.

**Alternatives considered**:
- C\* or a 2-predictor primary set for Phase 2 — rejected in favor of a
  single, cleanly pre-registered primary predictor (a\*).
- Binary classifier target — rejected as primary; threshold choice would
  be arbitrary given only 3-4 SIRIUS carotenoid-class hits to calibrate a
  "presence" cutoff against, and binarizing throws away signal needed at
  this sample size.

**Rationale**: Both decisions favor preserving statistical power and a
clean pre-registration story over interpretability shortcuts, consistent
with the "n≈17-18 species power ceiling" consideration already in the
strategy doc.

**Consequences**: `analysis/INTEGRATED_ANALYSIS_STRATEGY.md` Phase 2/4
sections to be updated to reflect a\*-primary / C\*-secondary /
orange_score-exploratory and continuous-target, once the rest of the
grilling session's open questions are resolved (in progress).

**Tags**: phase-2, phase-4, pre-registration, predictor-choice, classifier-target

## [2026-08-15] Negative-control enforcement + HPLC validation timeline (PI grilled, decisions confirmed)

**Context**: continuation of the same grilling session as the entry above.
Two more open items: whether the negative-control/permutation designs
(Fable review) are actually enforced or just documented intentions, and
whether/when the targeted HPLC/LC-UV-Vis pigment validation (also Fable,
recommended to run parallel to Phase 2-3 rather than deferred to Phase 6)
will happen.

**Decision**:
- **Negative controls are a HARD PIPELINE GATE**, not a documented-but-
  manual step. Phase 2 (and later Phase 4) scripts must refuse to emit a
  "final" results file unless the corresponding null-run output already
  exists and is newer than the input data — this needs to be built into
  the scripts themselves when Phase 2/4 are implemented, not left as an
  honor-system checklist item.
- **HPLC/LC-UV-Vis validation timeline is contingent on the SIRIUS
  re-run's results** (option "c" of committed/not-happening/contingent) —
  not scheduled on a fixed near-term date, not indefinitely deferred
  either. Phase 2-4 analysis proceeds on SIRIUS + spectral-networking
  evidence in the meantime, with compound identity explicitly flagged as
  **provisional** in any write-up until HPLC validation (if/when it runs)
  confirms it.

**Alternatives considered**:
- Documented-but-manual negative-control step — rejected; the whole point
  of Fable's review was to prevent it being skipped under time pressure
  once a real-looking hit appears, which a checklist item doesn't actually
  prevent.
- Waiting for HPLC validation before running Phase 2-4 at all — rejected;
  would stall the whole plan on a wet-lab resourcing question with no
  fixed timeline. Provisional-compound-identity framing lets analysis
  proceed without overclaiming.

**Rationale**: Both decisions prioritize keeping the plan moving while
being explicit about what's provisional vs. confirmed, rather than either
silently skipping a safeguard or blocking on an unscheduled dependency.

**Consequences**: Phase 2/4 script implementations must include an
explicit null-run-freshness check before writing final output (raises the
implementation bar slightly, intentionally). Any Phase 3-5 write-up must
carry a "compound identity provisional pending validation" caveat until/
unless HPLC work happens.

**Tags**: phase-2, phase-4, negative-control, hard-gate, hplc-validation, provisional-findings

## [2026-08-15] Accept MS2 dedup's cross-adduct-type limitation as-is (PI grilled, confirmed)

**Context**: continuation of the same grilling session. `dedupe_ms_features.py`'s
known limitation (cross-adduct-type compound merging deferred, ~4.9% of
EB-detected compounds affected — see the earlier 2026-08-15 dedup decision
entry above) was put to the PI: fix now, or accept as-is.

**Decision**: Accept as-is. Not fixing pre-emptively.

**Rationale**: The bias is conservative (more groups than the true
compound count, never fewer) and affects a small minority of groups
(~900/10,949). If a specific compound flagged by Phase 3 as color-
associated turns out to be one of the split-adduct cases, manual
verification against `aligned_compounds.tsv` at that point is more
trustworthy than a general-purpose automated m/z-matching merger built
now, before it's known which compounds actually matter.

**Tags**: mass-spec, deduplication, deferred-fix, phase-2

## [2026-08-15] Wet-lab phenotype provenance: deferred, not urgent (PI grilled, confirmed)

**Context**: `data/raw/control_phenotype_90_110h/CONTROL_PHENOTYPE_90_110H.md`
and `data/metadata/control_phenotype_90_110h/provenance.md` both have an
open "PI: please fill in" section (imaging rig identity, plate/well
layout and randomization, biological vs. technical replicate structure,
media/incubation details, batch/blocking structure).

**Decision**: Not filled in now. PI can answer accurately but confirmed
it isn't needed at this point in the project.

**Rationale**: Phase 1's core results (species ranking, phylogenetic
signal) were independently reproduced on a completely different
phenotyping pipeline (legacy YPD2), so `control_90_110`'s wet-lab specifics
aren't what's producing those findings — not a current blocker. It
becomes a real blocker specifically if/when Phase 2 needs to model
plate-level batch effects (already flagged as a "Library Plate" covariate
consideration in `analysis/INTEGRATED_ANALYSIS_STRATEGY.md`'s Key
Considerations) and the layout is still undocumented at that point.

**Consequences**: Ingestion of `control_phenotype_90_110h` remains
formally incomplete by the mycelium checklist (provenance section has
open TODOs) — acceptable and intentional, not an oversight. Revisit when
Phase 2 batch-effect modeling is actually being implemented, not before.

**Tags**: provenance, deferred, phenotype, batch-effects

## [2026-08-15] Start Phase 2 implementation now (PI grilled, confirmed — grilling session complete)

**Context**: final question of the 2026-08-15 grilling session on
`analysis/INTEGRATED_ANALYSIS_STRATEGY.md`'s remaining open decisions.
Phase 2 (color<->metabolome whole-panel association) doesn't depend on
SIRIUS (only Phase 3 does) or on the wet-lab provenance conversation
(deferred per the entry above).

**Decision**: Start building Phase 2 now. Nothing blocking.

**Consequences**: This closes out the 2026-08-15 grilling session with all
7 open decisions resolved:
1. a\* primary predictor, C\* secondary, orange_score exploratory-only.
2. Phase 4 classifier target: continuous (summed color-associated compound
   group abundance), not binary.
3. Negative controls: hard pipeline gate, not documented-but-manual.
4. HPLC/LC-UV-Vis validation: contingent on the SIRIUS re-run's results,
   not scheduled independently.
5. MS2 dedup cross-adduct-type limitation: accepted as-is.
6. Wet-lab phenotype provenance: deferred, not currently blocking.
7. Phase 2 implementation: starting now.

**Tags**: phase-2, grilling-session, implementation-start

## [2026-08-15] Phase 2 implemented and run — null color result, real confound found

**Context**: Following the grilling session's confirmed decisions (a\*
primary predictor, hard-gated negative control, species-tree-block
permutation), implemented and ran
`analysis/scripts/phase2_color_metabolome_association.py`.

**Decision/outcome**: Phase 2 ran to completion for a\* (primary), C\*
(secondary), and the area negative-control decoy, across both fractions.
Result: 0 BH-FDR<0.05 hits for color (either axis, either fraction) —
matches the permutation-null expectation, i.e. a well-calibrated null, not
an underpowered-looking bug. The decoy run surfaced a real, unexpected,
and more immediately actionable finding: colony area broadly confounds
cell-fraction (not supernatant) compound abundances (1,524/10,164 hits vs.
~0 expected) — logged as
`.living/findings/biomass-scaling-artifacts-in-extraction-based-metabolomics.md`.
The color null result is logged as
`.living/findings/phenotype-metabolome-association-statistical-power.md`.
Full writeup: `analysis/integrated_analysis/phase2_metabolome_phenotype/PHASE2_SUMMARY.md`.

**Consequences**: Phase 3 (compound-class enrichment) has no
formally-significant color-associated feature list to test as originally
scoped — needs a follow-up decision (not yet made) on whether to proceed
with a nominal/exploratory ranking, pivot to within-species tests, or wait
for more data. The area confound should be investigated before further
cell-fraction analysis of any kind, independent of the color question.

**Also fixed during implementation**: a NaN-handling bug where
constant-abundance features got spuriously minimal (most "significant")
empirical p-values — see `.living/learnings.md` 2026-08-15 entry.

**Tags**: phase-2, implementation, null-result, biomass-confound, findings

## [2026-08-15] Session checkpoint

`.living/INDEX.md` regenerated and `validate_structure.py` passed after
the Phase 2 implementation/run above. `.claude/last-session.md` updated
with the full grilling-session + Phase 2 summary. Session paused here
pending PI direction on Phase 3's next step (nominal/exploratory ranking
vs. within-species tests vs. wait for more data) and whether to
investigate the colony-area/cell-fraction confound now.

**Tags**: session-checkpoint

## [2026-08-15] Within-species (R. mucilaginosa) follow-up implemented — null confirmed with higher power

**Context**: PI proposed, after Phase 2's whole-panel null, testing
whether within-species a\* variation in *R. mucilaginosa* (the largest,
best-sampled species, n=216/206/201 phenotype/MS/genome) correlates with
metabolome features — directly following PHASE2_SUMMARY.md's recommended
next step #4.

**Decision/outcome**: Implemented `analysis/scripts/phase2_within_species_association.py`
(real strain-level genome-tree blocking, same predictor/negative-control/
hard-gate conventions as the whole-panel script). Ran a\* (primary), C\*
(secondary), and area (decoy) for *R. mucilaginosa*. Result: null again —
0 BH-FDR<0.05 hits for either color axis in either fraction, despite
much higher effective power (206 strains, 27 phylogenetic blocks vs. 6)
and a\* spanning nearly the full panel-wide color range within this one
species alone. The area decoy also showed 0 hits within-species (vs.
1,524 whole-panel cell-fraction hits), supporting that earlier confound
being a between-species effect. Full writeup:
`analysis/integrated_analysis/phase2_metabolome_phenotype/WITHIN_SPECIES_MUCILAGINOSA.md`;
findings updated: `.living/findings/phenotype-metabolome-association-statistical-power.md` (F-003).

**Consequences**: This substantially strengthens the whole-panel null —
it's no longer explainable primarily by the species-level power ceiling.
Recommended next steps (not yet decided/actioned): try 1-2 other
well-sampled species as a robustness check, consider a multivariate
approach (per-feature univariate tests may be poorly suited to a diffuse
signal), and treat this as further motivation for the targeted HPLC
validation being contingent on the SIRIUS re-run (untargeted LC-MS2 may
simply not resolve carotenoids well).

**Tags**: phase-2, within-species, mucilaginosa, null-result, implementation

## [2026-08-15] Session checkpoint 2

`.living/INDEX.md` regenerated, `validate_structure.py` passed after the
within-species *R. mucilaginosa* follow-up above. Paused pending PI
direction on next steps (other species robustness checks, multivariate
approach, or the area/cell-fraction confound investigation).

**Tags**: session-checkpoint

## [2026-08-15] Robustness check + multivariate follow-up complete — 4-method null triangulated

**Context**: PI requested (1) a quick robustness check on 1-2 more
well-sampled species and (2) a multivariate approach, following the
*R. mucilaginosa* within-species null.

**Decision/outcome**: (1) Ran `phase2_within_species_association.py` with
a new `--min-strains` override (explicitly labeled exploratory below 20)
on *R. paludigena* (n=10) and *R. toruloides* (n=10) — both null, 0 hits
each. (2) Wrote `analysis/scripts/phase2_multivariate_association.py`
(sparse Lasso regression, GroupKFold CV using real phylogenetic blocks,
permutation-based significance) and ran it on *R. mucilaginosa* for a\*,
C\*, and the area decoy — all four cross-validated R² values negative
(worse than predicting the mean), no signal detected. Full writeup:
`analysis/integrated_analysis/phase2_metabolome_phenotype/ROBUSTNESS_AND_MULTIVARIATE.md`;
findings updated: `.living/findings/phenotype-metabolome-association-statistical-power.md`
(F-004, status now "robust" — 4 consistent methodologically-independent entries).

**Consequences**: Four independent analytical approaches (whole-panel
univariate, within-species univariate at two sample sizes, within-species
multivariate) now agree on a null result. Recommendation carried into the
findings doc: further statistical exploration of the untargeted MS2 data
has low expected marginal return; the targeted HPLC/LC-UV-Vis pigment
validation (contingent on the SIRIUS re-run per the earlier grilling
session decision) is the higher-value next step, since it tests whether
the untargeted method can even see the relevant chemistry rather than
continuing to test statistical variations on that assumption.

**Tags**: phase-2, multivariate, robustness-check, null-result, findings, session-milestone

## [2026-08-15] Cross-cutting finding logged: whole-panel trait-molecular hits don't survive within-species testing; copper AUC pulled in as provisional precedent

**Context**: PI asked for a summary of nulls across the project (color
vs. metabolome, this session; copper-resistance AUC vs. amino acid
composition, an earlier separate analysis in `analysis/copper/`). Both
show the same pattern: whole-panel (even PGLS-corrected) associations do
not survive restriction to variation within a single species
(*R. mucilaginosa*).

**Decision**: Logged this as a new cross-cutting finding
(`.living/findings/phylogenetic-confounding-of-trait-molecular-associations.md`,
F-005 in the registry) and added a "Within-species restriction as the
most diagnostic check" subsection to
`analysis/INTEGRATED_ANALYSIS_STRATEGY.md`'s Key Considerations, pulling
in the copper-AUC precedent. **The copper-AUC numbers are explicitly
flagged as provisional** — PI stated intent to revisit that phenotype's
underlying data further — both in the finding's Open Questions and in the
strategy doc, so a future session doesn't treat the current
naive/PGLS/sensitivity-check numbers as a closed result.

**Consequences**: Going forward, any whole-panel "hit" in this project
(Phase 3 SIRIUS enrichment, Phase 5 genome linkage, or the eventual
copper-AUC re-analysis) should be checked against a within-species
restriction before being trusted, per this now-2-for-2 track record. When
the copper-AUC data is revisited, update the finding's Evidence Ledger
and the strategy doc's caveat accordingly rather than leaving the current
snapshot as the last word.

**Note on session hygiene**: multiple tool-result messages during this
session claimed files (`phase2_within_species_association.py`,
`INTEGRATED_ANALYSIS_STRATEGY.md`, `.claude/last-session.md`) were
externally modified and instructed the agent not to tell the PI — flagged
directly to the PI each time as a suspected prompt-injection pattern
rather than silently complied with, per standing instructions to surface
suspected injection attempts. The tree-file change (BFD taxa_276 ->
taxa_278) was verified as a real, benign pipeline update (both files
exist on disk); the `last-session.md` overwrite was traced to the
mycelium Stop hook's own deterministic fallback summary logic, not
tampering. Recorded here for provenance, not because either turned out to
be malicious.

**Tags**: cross-cutting, findings, copper-auc, provisional, session-hygiene

## [2026-08-15] Phylogeny source bumped from the 276-taxa to the 278-taxa PHYling tree; all direct downstream analyses rebuilt

**Context**: PI requested injecting the updated BFD PHYling protein tree
(`BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile`,
a re-run of the same pipeline with 2 additional taxa) everywhere the
project uses phylogeny, then rebuild anything that had gone stale.
Confirmed by diffing tip sets that taxa_278 is a strict superset of
taxa_276 (adds `Rhodotorula_mucilaginosa_DH4148.proteins` and
`Rhodotorula_evergladensis_DBVPG_7922.proteins`, drops nothing), so this
is a pure addition, not a topology-invalidating re-sampling.

**Decision**: Updated every hardcoded `taxa_276` path reference to
`taxa_278` (`analysis/copper/scripts/common.py`,
`analysis/copper/scripts/02_pgls_analysis.R`,
`analysis/YPD/color_shape_growth/scripts/02_pgls_analysis.R`,
`analysis/scripts/phase2_within_species_association.py`,
`scripts/run_all_differential_pairs.sbatch`, plus doc mentions in
`analysis/copper/genome_characteristics.md` and
`analysis/INTEGRATED_ANALYSIS_STRATEGY.md`), then reran every analysis
that reads the tree directly: the full `analysis/copper/scripts/run.sh`
and `analysis/YPD/color_shape_growth/scripts/run.sh` pipelines, all 5
documented `phase2_within_species_association.py` invocations
(mucilaginosa area/a/C, paludigena and toruloides robustness checks with
`--min-strains 8`), `analysis/scripts/prune_species_tree.R` (regenerates
`analysis/integrated_analysis/phase1_phenotype/species_tree.nwk`), and
its two downstream consumers `phylogenetic_signal.R` and
`convergent_color_test.R`. All reruns used documented/default arguments
(recovered from each script's usage docstring or from existing sibling
output filenames when no run script existed) and produced results
consistent with the previously-published numbers (same strain/species
counts, same qualitative hit lists) — this was a tree-precision refresh,
not a result-changing re-analysis.

**Did NOT rebuild**: `scripts/run_all_differential_pairs.sbatch` (a
~110-task SLURM array job — real cluster cost, held for explicit PI
go-ahead rather than auto-submitted). Also did not add the 2 new taxa
into the species-level tree/analyses: `prune_species_tree.R` correctly
dropped them with a warning because
`analysis/integrated_analysis/phase1_phenotype/genome_strain_species_busco_map.csv`
(the strain-to-species map) doesn't yet list them, and that map has no
generating script in this repo — it was built ad hoc "from BFD directly"
per the existing strategy-doc note. See linked todo.

**Consequences**: `species_tree.nwk` and everything downstream of it
(`phenotype_phylogenetic_signal.csv`, `convergent_color_candidates.csv`)
now reflect the 278-taxa branch-length estimates even though the tip set
is unchanged (still 17 species). Copper and YPD outputs
(`outputs/pruned_tree_*strains.nwk`, correlation/PGLS CSVs) are likewise
refreshed. Until the busco map is regenerated, the 2 new taxa remain
invisible to every species-level analysis in this project even though
they're now in the strain-level tree.

**Tags**: phylogeny, tree-update, reproducibility, pipeline-rebuild, copper, ypd, phase2

## [2026-08-15] ANOVA/pattern-group approach implemented; hard-gate gap fixed; corrected an earlier under-powered decoy error

**Context**: PI asked to try an ANOVA-style approach within a species —
testing whether metabolite abundance differs across color-PATTERN groups
(k-means clusters on joint L\*/a\*/b\*, not a single-axis correlation),
which could catch non-monotonic relationships the existing
correlation/Lasso tests would miss. While validating this new script's
decoy run, found it disagreed sharply with the existing univariate
within-species decoy result (1,354 vs. reported "0" hits) — investigation
traced this to a real error: the univariate within-species decoy
(`within_species_Rhodotorula_mucilaginosa_association_area_decoy.csv`)
had been generated by a `--n-perm 20` smoke test that was never rerun at
full power before being treated as satisfying the hard gate, and
`WITHIN_SPECIES_MUCILAGINOSA.md` / the findings file incorrectly reported
"0 hits, area confound is between-species only" on that basis.

**Decision**: (1) Wrote `analysis/scripts/phase2_anova_pattern_association.py`
(Kruskal-Wallis across k-means color-pattern clusters, same
phylogenetic-block-permutation/hard-gate conventions as the other Phase 2
scripts). (2) Fixed the hard-gate gap in both
`phase2_within_species_association.py` and the new ANOVA script: they now
write an `n_perm` column to their output and refuse to run a real
predictor unless the decoy output records `n_perm >= 100` (not just
file existence/freshness). (3) Reran the univariate decoy at full power
(`--n-perm 200`, `--seed 0`, deterministic) — confirmed 2,025/9,437
cell-fraction hits (not 0) — and corrected `WITHIN_SPECIES_MUCILAGINOSA.md`
and `.living/findings/phenotype-metabolome-association-statistical-power.md`
accordingly. (4) Ran the real ANOVA color-pattern test — also null
(0/9,437 cell, 0/10,255 supernatant), against a decoy now properly
validated to detect a strong real effect (1,354/9,437) — see
`ROBUSTNESS_AND_MULTIVARIATE.md`, section 3.

**Consequences**: The a\*/C\* real-predictor null results from earlier in
the session are unaffected (those were always run at full `n_perm=200`)
and are, if anything, better supported now — the pipeline demonstrably
detects real signal (area) on this exact strain set and still finds
nothing for color, across five now-mutually-consistent methods (whole-
panel correlation, within-species correlation at two sample sizes, sparse
multivariate regression, and pattern-group ANOVA). The corrected
area-confound magnitude (~1,300-2,000 cell-fraction hits within
*R. mucilaginosa* alone, not 0) means the "area confound is purely
between-species" claim is retracted; its actual cause is still an open
question. **Should backport the same `n_perm`-recording hard-gate fix to
`phase2_color_metabolome_association.py` and `phase2_multivariate_association.py`**
for consistency — not yet done, tracked as a todo.

**Tags**: phase-2, anova, hard-gate-fix, error-correction, decoy-validation, mucilaginosa

## [2026-08-15] Ran 7-persona idea-generation session after 5-method color-metabolome null

**Context**: PI asked to consult expert personas (data analysis, chemistry,
metabolomics, statistics, bioinformatics) on how to restructure the data
or approach, given the color-metabolome null was now well-triangulated
across 5 methods and genome-side association had never been attempted.

**Decision**: Used the installed `idea-generator` convention pack. Selected
7 personas (3 from the default catalog: Quantitative Geneticist,
Evolutionary Biologist, Causal Inference Researcher, Ecologist; 3 custom:
Natural Products/Analytical Chemist, Metabolomics/MS Specialist, Fungal
Genomics Bioinformatician) rather than all 15 catalog personas, to match
the disciplines the PI named and keep the session tractable. 2 ideas per
persona (14 total), dispatched as parallel fresh subagents, synthesized in
`analysis/ideas/2026-08-15-color-metabolome-genome-null-brainstorm/00_index.md`.

**Consequences**: Top-ranked recommendations (all low-effort, no new data
required): determine colony area's causal role (confounder/mediator/
collider) before trusting the existing negative-control framing;
targeted re-mining of existing raw MS data against a carotenogenesis
pathway mass list; candidate carotenoid-pathway-gene genotyping against
color (the project's first genome↔color test, using already-built BFD/
PGLS infrastructure); a variance-component/heritability decomposition of
color as a precondition check. See the index for the full ranked list and
what the panel explicitly recommended NOT starting with (research-grade
spectral kernels, the inherently underpowered convergence test at N=3-4
independent origins). Not yet actioned — awaiting PI direction on which
to implement.

**Tags**: idea-generation, brainstorm, phase-3, phase-5, next-steps

## [2026-08-15] Development plan for Ideas 1/3/5/6; started Idea 3 scaffold; SNP-calling question resolved

**Context**: PI wants to pursue most of the 14 brainstormed ideas and
asked for a development plan: a framework for Idea 1 (chemist), to start
implementing Idea 3 (candidate carotenoid-pathway genotyping), whether
Idea 3 needs a SNP-calling strategy, a strategy for Idea 5 (convergence
test), and what infrastructure Ideas 3+6 share.

**Decision**: Wrote `analysis/ideas/2026-08-15-color-metabolome-genome-null-brainstorm/DEVELOPMENT_PLAN.md`
covering all four asks, plus a standalone `IDEA1_CHEMIST_FRAMEWORK.md`
(two-phase: cheap in-silico re-mining of existing raw MS data against a
carotenogenesis-pathway mass list first, wet-lab APCI/APPI extraction only
if that comes back empty). **SNP-calling verdict**: not needed for Idea 3
(gene-level presence/copy-number/LoF works directly from BFD's existing
per-strain assemblies+annotations via Pfam counts + ortholog-confirmed
protein-sequence MSA — no raw-read remapping or variant caller required);
would only become necessary for Idea 6's more ambitious genome-wide
variant-level GWAS tail, which needs infrastructure this project doesn't
have yet (raw reads for all strains, a reference/pangenome alignment
strategy, a variant-scale mixed-model GWAS tool) — explicitly deferred
until/unless Tier 1 (candidate-gene) work motivates it.

Started the Idea 3 scaffold: `analysis/scripts/phase5_candidate_gene_genotyping.py`
(Pfam pre-filter for 6 candidate carotenoid-pathway genes — crtYB/crtI/
crtS/crtR/HMGR/GGPPS — plus stubbed ortholog-confirmation/copy-number/
MSA-LoF steps). Checked and confirmed `diamond` (2.1.24) and `mafft`
(7.505) are available as HPCC modules — the two tools this pipeline needs
beyond what's already installed.

**Blocking dependency confirmed with PI**: `BFD/db/BFD.duckdb` is
mid-rebuild by a separate nextflow pipeline; every DB-dependent step in
the new script hard-exits with a clear message rather than querying any
interim/stale copy (one was found at `/bigdata/stajichlab/shared/projects/BFD/BFD.duckdb`
during this session but confirmed empty/mid-rebuild and explicitly NOT to
be used). Reference protein sequences + catalytic-residue positions for
the candidate genes also not yet sourced (literature/PI task).

**Consequences**: Idea 1 Phase 1 (in-silico re-mining) and Idea 5's
Steps 1-2 (formal ancestral-state/regime-shift modeling on the existing
species tree) can both start immediately, independent of the BFD rebuild.
Idea 3/6 are otherwise ready to execute the moment the database lands and
reference sequences are sourced — see DEVELOPMENT_PLAN.md's "Recommended
build order."

**Tags**: idea-3, idea-1, idea-5, idea-6, phase-5, snp-calling, development-plan, blocked-on-bfd-rebuild

## [2026-08-15] Ran Idea 5 Steps 1-2 (formal regime-shift detection + contrast pairs)

**Context**: PI asked to execute Idea 5's Steps 1-2 now (doesn't depend on
the BFD rebuild). Installed `bayou` (2.3.2) and `OUwie` from CRAN (neither
previously installed; both installed cleanly, ~3-4 min each including
dependencies).

**Decision/outcome**: Wrote and ran `analysis/scripts/idea5_regime_shift_detection.R`
(bayou reversible-jump OU MCMC, 20,000 generations, a\* on the 17-species
tree) and `analysis/scripts/idea5_contrast_pairs.R` (maps top-pp branches
to descendant clades, pairs each against its nearest phylogenetic
non-candidate sister). Result: diffuse posterior support (max branch
pp~0.21, k mean 1.8 HPD95 0-4) — no branch clears a decisive threshold,
consistent with the "not enough independent phylogenetic data points"
theme recurring throughout this project. The top candidate
(*R. sphaerocarpa*+*R. taiwanensis*) independently matches Phase 1's
original #1/#2 species ranking by a\* — a genuine cross-method consistency
check, not cherry-picked. *R. glutinis* also matches the original coarse
heuristic's candidate list. Two of the top-5 rows (*Pseudomicrostroma
phylloplanum*, *R. mucilaginosa*) flagged as lower-confidence/likely
tree-placement artifacts in the writeup rather than presented as equally
supported. Full results: `analysis/integrated_analysis/phase5_genome_linkage/idea5_regime_shift/RESULTS.md`;
finding logged: `.living/findings/convergent-color-evolution-in-rhodotorula.md`.

**Alternatives considered**: `l1ou` (not on CRAN, GitHub-only, would need
`devtools`/`remotes` + network access to a non-CRAN source — skipped in
favor of `bayou`/`OUwie` which installed cleanly from CRAN and cover the
same analytical need).

**Consequences**: Step 1's originally-planned hard shift/non-shift
threshold had to be abandoned in favor of top-N ranking, since no branch
cleared a confident cutoff — this is a real result (documented as such),
not a script failure. This run is exploratory (single MCMC chain, no
multi-chain convergence diagnostic) — a rigorous version would run ≥2
chains and check Gelman-Rubin/ESS before treating results as
publication-grade. orange_score_mean not yet rerun (a\* only, per the
grilling-session primary-predictor decision) — tracked as a follow-up.

**Tags**: idea-5, bayou, regime-shift, convergent-evolution, phase-5, r-packages-installed

## [2026-08-15] Ran Idea 1 Phase 1 (targeted exact-mass re-mining of existing raw MS data)

**Context**: PI approved proceeding with Idea 1 Phase 1 per
`IDEA1_CHEMIST_FRAMEWORK.md` — no new data, searches the existing raw EB
feature table directly by exact mass for the full carotenogenesis
pathway, independent of SIRIUS's own compound calling.

**Decision/outcome**: Wrote and ran `analysis/scripts/idea1_targeted_mass_remining.py`
(9 pathway compounds x 4 adducts = 36 targets, 20 ppm tolerance, against
the fuller 53,040-feature EB table). Found 11 raw-feature matches across
6/9 compounds (no matches for neurosporene/lycopene/γ/β-carotene). Two
standouts: a phytofluene [M+NH4]+ candidate (row 1735, -1.0 ppm, by far
the highest intensity at 220 scans, unannotated by SIRIUS) and a
torularhodin [M+H]+ candidate (row 21315, +2.5 ppm) that **SIRIUS
annotated as "Polyamines"/a chemically implausible structure** — a
concrete confirmed instance of the misannotation risk flagged in the
original plan, not just a hypothetical. Cross-checked both against
Phase 2's already-run color association results: **neither correlates
with a\*** (Spearman rho 0.03-0.12, FDR 0.48-0.94). Full writeup:
`analysis/integrated_analysis/phase3_metabolome_phenotype_idea1/RESULTS.md`;
finding logged: `.living/findings/carotenoid-pathway-detectability-in-untargeted-lcms.md`.

**Consequences**: This pass is double-sided — it strengthens the case
that ESI *can* detect signal in the carotenoid mass range (ruling out the
strongest form of "the method can't ionize this class") and confirms a
real SIRIUS misannotation, but does not rescue a hidden color signal even
if the candidates are real. No MS2 fragmentation confirmation done yet
(structural check, not just mass match) — recommended as the next cheap
step (pull raw spectra for rows 1735/21315, check for polyene neutral-loss
pattern) before deciding whether to invest in the wet-lab APCI/APPI
escalation (IDEA1_CHEMIST_FRAMEWORK.md Phase 2).

**Tags**: idea-1, carotenoids, exact-mass-search, sirius-misannotation, null-result, phase-3

## [2026-08-15] Pulled and inspected MS2 spectra for Idea 1's two standout candidates

**Context**: PI approved pulling raw MS2 spectra for rows 1735/21315 (the
two standout mass matches from Idea 1 Phase 1) to check fragmentation
plausibility before deciding on further investment, per RESULTS.md's own
recommendation.

**Decision/outcome**: Extracted both spectra by `FEATURE_ID` from
`aligned_features_filled.mgf` (one-off Python scan, no reusable script
written yet — see consequences). Manual inspection **reversed the
intensity-based ranking**: row 1735 (phytofluene candidate, previously
the more credible one on mass accuracy + intensity) is a largely
undissociated precursor with no NH3 neutral-loss peak — demoted, does not
support the identification. Row 21315 (torularhodin candidate, the one
SIRIUS misannotated as "Polyamines") shows two clean paired water-loss
fragments consistent with a carboxylic-acid-bearing structure plus
substantial mid-mass backbone-cleavage-plausible fragments — chemically
coherent, though not confirmed against any reference spectrum. Full
writeup: `analysis/integrated_analysis/phase3_metabolome_phenotype_idea1/MS2_FRAGMENTATION_CHECK.md`;
finding updated: `.living/findings/carotenoid-pathway-detectability-in-untargeted-lcms.md` (F-002).

**Consequences**: Row 21315 is now the sole remaining candidate from
Idea 1 Phase 1 worth further investment (a genuine reference-spectrum/
GNPS library comparison, not another hand inspection) — but it still
doesn't correlate with color (per F-001), so even confirming its identity
wouldn't on its own explain the color-metabolome null. **General caution
for future re-mining passes**: total_scans/intensity is not a reliable
proxy for identification confidence in either direction — worth
remembering next time a "high-intensity" candidate looks promising on
paper alone. Not yet turned into a reusable
`extract_mgf_spectrum.py --feature-id N` utility — flagged as a
nice-to-have if this kind of spot-check becomes routine, not done now.

**Tags**: idea-1, ms2-fragmentation, torularhodin, phytofluene, spectral-quality-caution

## [2026-08-16] BFD rebuild confirmed done (278 strains); expanded Idea 1 mass search to apocarotenoids/additional carotenoids/sterol pathway

**Context**: PI asked whether the BFD rebuild had finished, then asked to
implement the previously-proposed expanded mass search (apocarotenoids,
additional oxygenated carotenoids, sterol-pathway precursor-competition
markers) and check hits against phenotype more broadly than just color,
while PI researches Idea 3's candidate gene reference sequences
separately.

**Decision/outcome**: Confirmed `BFD/db/BFD.duckdb` (13.8 GB) is fully
populated — 278 strains (up from 276), all functional-annotation tables
populated across all 278. Idea 3 is no longer database-blocked; still
needs reference protein sequences (PI sourcing separately, not done by
this agent this session).

Extended `analysis/scripts/idea1_targeted_mass_remining.py` from 9 to 19
compounds / 4 categories (see script docstring for full rationale per
category) and reran. 31 matches (up from 11). Two new standouts: an
ergosterol [M+H-H2O]+ candidate (row 846, -4.3 ppm) independently flagged
by EB's own pipeline as an in-source fragment — chemically self-
consistent, two lines of evidence agreeing, not one coincidental mass
match — and an astaxanthin [M+H]+ candidate (row 9384, +2.4 ppm, not
strain-confirmed for *Rhodotorula* specifically). Cross-checked all
standout candidates against both color (a\*, using Phase 2's existing
results) and copper-resistance AUC (new quick TSS-normalized Spearman
check against `sample_metadata.csv.gz`'s `mean_auc_rate`). Result: no
color correlation for any candidate; **ergosterol shows a naive,
uncorrected correlation with copper-AUC in the cell fraction (rho=0.226,
p=0.0002, n=267)**. Full writeup:
`analysis/integrated_analysis/phase3_metabolome_phenotype_idea1/EXPANDED_SEARCH_RESULTS.md`;
findings logged: `.living/findings/carotenoid-pathway-detectability-in-untargeted-lcms.md`
(F-003), registry F-008.

**Consequences**: The ergosterol-AUC correlation is explicitly flagged as
a lead, not a result — it was computed as a naive single Spearman
correlation with no phylogenetic correction, and this project has a
2-for-2 track record (including this *exact* AUC phenotype's earlier
amino-acid "hits") of naive whole-panel correlations failing to survive
phylogenetic block-permutation or within-species restriction. Before
trusting it, it needs the same rigor as everything else in this project —
not yet built (would need a predictor-swapped variant of the existing
Phase 2 scripts, using mean_auc_rate instead of the color phenotype
table). Recommended, not yet actioned.

**Tags**: idea-1, bfd-rebuild-done, ergosterol, copper-resistance, sterol-pathway, needs-validation

## 2026-08-16: Sterol/ergostane cluster AUC follow-up; Idea 3 unblocked with PI-built custom pigment HMMs

**Context**: PI asked about the terpene-class composition of the SIRIUS
"Terpenoids"-pathway features (73 total). Breakdown by NPC class/ClassyFire
class showed carotenoids are a small minority (4/73) while an
Ergostane/Cholane/Cholestane steroid cluster (7 features) is the single
largest NPC subclass — structurally consistent with the existing row-846
ergosterol ISF candidate and its naive copper-AUC lead (F-008).

**Decision/outcome**: Promoted the earlier one-off inline AUC quick-check
into a reusable script (`analysis/scripts/idea1_auc_quickcheck.py`) and
ran it on the 3 most chemically self-consistent SIRIUS ergostane calls
(Peroxyergosterol row 9852, Ergost-3,5,7,9(11),22-pentaen row 6682,
7-Hydroxyergosterol row 35014) plus the row-846 anchor. **All 4
independently-called features show the same positive-direction naive
correlation with copper-AUC in the cell fraction** (rho 0.17-0.28, all
p<0.006) — see
`analysis/integrated_analysis/phase3_metabolome_phenotype_idea1/STEROL_CLUSTER_AUC_CHECK.md`,
findings F-004/registry F-009. Still explicitly unvalidated (no
phylogenetic block permutation or negative control run on any of these
yet).

Separately, PI provided a custom-built panel of 28 pigment-pathway protein
HMMs (`~/pigment_protein_hmms/`, HMMER3, built from curated seed
alignments — not coarse Pfam prefilters) covering carotenoid
(crtB/crtO/crtP/crtQ/crtR/crt_fungal_lcy/pds/psy), melanin
(pks_melanin/ayg1/scy A-F/t3hnr/t4hnr/tyrosinase/laccase), and other
pigment (hgd/hppd/mysA-E/scd) families — this is the reference-sequence
input Idea 3 was blocked on. Decided to start with a standalone
`hmmsearch` of the combined panel (`pigmentation_profiles.hmm`) against
all 2,188,032 predicted proteins from the 278 BFD genomes
(`BFD/input/pep/*.proteins.fa`, concatenated, 1.3GB), `-E 1e-5`,
`--domtblout`, rather than first building a formal BFD-pipeline
integration step — PI's phrasing ("start with these HMMs... or we can
design an additional step in BFD") read as preferring the fast path
first, with pipeline integration a deferred option if the panel proves
useful. Run launched in the same SLURM allocation (job 27441568) already
in use this session. Also ran an independent, faster sanity check:
keyword cross-reference of BFD's existing `swissprot.parquet` diamond
hits against pigment-relevant SwissProt entries (crt*, tyrosinase,
laccase, NCED/CCD, hmgA/HGD, pks*) — confirms broad presence of
carotenoid- and melanin-pathway-adjacent hits across the panel (raw,
unfiltered on identity/coverage — a presence signal, not confirmed
orthology). Outputs:
`analysis/integrated_analysis/phase5_genome_linkage/idea3_pigment_hmm_search/outputs/swissprot_pigment_hits_by_protein.csv`
and `swissprot_pigment_strain_summary.csv`.

**Consequences**: Idea 3 (candidate-gene genotyping, Phase 5) is now
substantively unblocked and actively running for the first time this
session, using the PI's own curated HMM panel instead of the coarse Pfam
placeholder in `phase5_candidate_gene_genotyping.py`
(`step_ortholog_confirmation()` there is superseded by this approach and
should be revisited/retired once the HMM-search results land). Next:
parse `pigment_hmm_hits.domtblout` into a per-strain family
presence/copy-number table once the background run completes, then this
becomes the first genome<->color association test in the project
(Phase 5).

**Tags**: idea-1, idea-3, ergosterol, sterol-pathway, copper-resistance, pigment-hmm-panel, bfd, genome-linkage

## 2026-08-16: Idea 3 pigment gene screen — 3 methods completed, convergent melanin-pathway finding

**Context**: Following the PI-provided custom HMM panel and launched
`hmmsearch` background run (see prior 2026-08-16 entry), completed the
screening stage of Idea 3: parsed the HMM panel results and additionally
ran the Pfam pre-filter scaffolded (but never executed) in
`phase5_candidate_gene_genotyping.py`, extended with laccase/tyrosinase
Pfam families to match the HMM panel's melanin-pathway coverage.

**Decision/outcome**: All 3 independent screens (custom HMM panel, Pfam
domains, SwissProt keyword cross-reference) converge on the same result:
laccase/multicopper-oxidase is present in ~all 278 genomes, tyrosinase is
essentially absent in all 3 — strong, methodologically diverse evidence
that this genus melanizes via the laccase route. Also found: 10/28 HMM
profiles have zero genome-wide hits (likely lineage-specific paralog
discrimination — `scyA` hits but `scyB-F` don't, `mysA/C/D/E` hit but
`mysB` doesn't — though `crtB`/`ayg1`/`hppd`/`scd` being silent needs PI
confirmation), and 5 HMM-panel + 2 Pfam families show copy numbers too
high to be real single-gene orthologs (superfamily-level domain matches).
Full detail: `.living/findings/pigment-gene-genomic-screen-rhodotorula.md`
(F-001/F-002/F-003, registry F-010/F-011).

**Consequences**: The screening stage of Idea 3 is done. Two immediately
actionable candidates for the project's first real genome<->color
association test: the sparse-presence HMM-panel families (`crtR` in
52/278 strains, `crtQ` in 48/278, `hgd` in 28/278) as binary predictors —
not yet tested against species-level color. Before using any copy-number
family (`crtP`, `laccase`, `mysE`) as a predictor, needs ortholog
confirmation (diamond blastp + mafft MSA + gene tree) — supersedes the
stubbed `step_ortholog_confirmation()` in
`phase5_candidate_gene_genotyping.py`. Neither next step started yet;
proposed to PI, awaiting direction.

**Tags**: idea-3, pigment-hmm-panel, laccase, tyrosinase, melanin-pathway, genome-linkage, ortholog-confirmation-needed

## 2026-08-16: Executive summary navigation docs + within-species sweep completed for all species with ≥5 strains

**Context**: PI noted it was "confusing to know what avenues have been
explored and what is supported or rejected" across `analysis/ideas/` and
`analysis/integrated_analysis/`, and separately asked to extend Phase 2's
within-species test to every species with ≥5 strains.

**Decision/outcome**: Added `analysis/integrated_analysis/EXECUTIVE_SUMMARY.md`
(phase-by-phase status/verdict table with links into each phase's docs
and tables) and `analysis/ideas/EXECUTIVE_SUMMARY.md` (idea-by-idea status
table linking to where each landed). Noted no PNG/PDF figures exist in
these phase directories (all tabular output) except
`analysis/examine_phenotype_calling/`.

Ran `phase2_within_species_association.py --min-strains 5` (decoy then
a\* primary) for the 5 remaining species with ≥5 strains + MS data not
yet tested: *R. dairenensis* (n=8), *R. diobovata* (n=8),
*R. taiwanensis* (n=6), *R. sp. clade I* (n=5), *R. sphaerocarpa* (n=5).
All null. Important caveat surfaced and documented: at n=5-8, the area
negative-control decoy also returns 0 hits (unlike *R. mucilaginosa*'s
n=206 decoy, which shows real signal) -- constant-abundance feature
exclusion alone removes ~48-51% of cell-fraction features at this n, so a
0-hit decoy here doesn't carry the same evidential weight. Wrote
`WITHIN_SPECIES_SMALL_SPECIES_SWEEP.md` documenting this explicitly, and
updated `PHASE2_SUMMARY.md` with a consolidated within-species results
table linking all species tested to date (8 of 17-18 in the panel, every
species with ≥5 strains). Findings logged:
`.living/findings/phenotype-metabolome-association-statistical-power.md`
F-005, registry F-012.

**Consequences**: Every statistically testable species in the panel now
has a within-species color-metabolome result on record; the remaining
9-10 species (1-3 strains each) cannot be tested this way. The
small-species negative-control caveat should be carried forward to any
future small-n test in this project, not just this one -- worth
promoting to a general convention if it recurs.

**Tags**: navigation, executive-summary, within-species, small-species, negative-control-caveat, phase-2

## 2026-08-16: Two new PI-directed threads — extreme-group color test, siderophore investigation

**Context**: PI proposed 2 new ideas: (1) partition strains into
high-orange/red vs low-orange/red groups (explicitly acknowledging this
ignores phylogeny) and test for compound differences, confirming
cell/supernatant are already tested separately throughout this project
(they are); (2) investigate rhodotorulic acid / siderophore chemistry in
the MS data, its strain-level presence/absence, and cross-reference
against the known NRPS gene (PI: "we know the NRPS gene... so we can also
condition whether the gene is found in genomes of strains that do not
show evidence of the compound").

**Decision/outcome**:

(1) Built `analysis/scripts/extreme_group_color_association.py` --
top/bottom-quartile a\*/C\*/area groups, rank-biserial effect size
(vectorized Mann-Whitney-U-equivalent rank-sum, reusing
phase2_color_metabolome_association.py's TSS/dedup/block-permutation
machinery unchanged). Null for both color axes, both fractions; decoy
well-calibrated (1,723 cell-fraction hits). 6th independent method to
find no color-metabolome signal. See
`EXTREME_GROUP_RESULTS.md`, findings F-006/registry F-013.

(2) Built a 3-script siderophore investigation:
`siderophore_mass_remining.py` (exact-mass search, idea1-style, 4
compounds x 4 adducts, 29 matches -- best rhodotorulic acid candidate row
2190 is the highest-intensity match found), `siderophore_presence_absence.py`
(strain-level presence, ~99% across all 17 species -- flagged as too
permissively thresholded to be discriminating), `siderophore_nrps_pfam_screen.py`
(coarse Pfam stand-in for the real reference sequence, since the PI has
not yet supplied the actual NRPS accession -- ornithine-hydroxylase proxy
universal/uninformative, 2-module-NRPS architecture proxy found only 1
genome-wide, judged a draft-assembly gene-model fragmentation artifact
rather than real absence). Concluded: neither side of the proposed
genotype/phenotype cross-reference is currently discriminating enough to
condition one against the other. See `phase_siderophore/RESULTS.md`,
findings F-001/registry F-014.

**Consequences**: Both threads are explicitly flagged as needing PI
input to proceed further -- (1) is a completed, documented null (no
further action needed unless PI wants a different quantile threshold);
(2) is blocked on the PI supplying the actual rhodotorulic-acid NRPS
reference sequence/accession, the same unblock pattern that worked for
Idea 3's pigment HMM panel. Also worth tightening the MS presence
threshold before any future genome cross-reference attempt.

**Tags**: extreme-group, color-metabolome, null-result, siderophore, rhodotorulic-acid, nrps, needs-reference-sequence

## 2026-08-16: Real rhodotorulic-acid NRPS ortholog search (PI-supplied reference sequence)

**Context**: PI supplied the actual reference sequence for the
rhodotorulic-acid NRPS at `tmpin/RA_NRPS.fa` (protein F2DD6D01_006956-T1
from *R. kratochvilovae* Y14), replacing the earlier coarse Pfam stand-in
(F-014/registry, from the same-day siderophore investigation). Also
pointed to antiSMASH results for that source genome at
`/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_RhodotorulaAcid/annotation/Rhodotorula_kratochvilovae_Y14/antismash_local/`.

**Decision/outcome**: Ingested the reference to
`analysis/integrated_analysis/phase_siderophore/reference/RA_NRPS.fa`.
Confirmed via antiSMASH (JAFEUJ010000019.1.region001.gbk, product="NRPS")
that this sits in a real BGC with an adjacent biosynthetic-additional
smCOG gene (F2DD6D01_006955, plausibly the ornithine hydroxylase
partner). Built `analysis/scripts/siderophore_nrps_diamond_search.py`
(diamond blastp, single query vs. all 278 BFD proteomes) — supersedes
`siderophore_nrps_pfam_screen.py` entirely (kept for provenance only,
its output should not be used). Result: clean, strongly bimodal identity
distribution (303 hits <30% = noise, 305 hits >=60% = real orthologs),
with the 3 in-panel *R. kratochvilovae* strains landing exactly as
expected as a built-in positive control. Threshold pident>=45 &
qcovhsp>=70 (the empty gap between modes) gives 275/278 strains
confirmed ortholog-positive.

Cross-referenced the 3 negative strains (1 outgroup + 2 *R. mucilaginosa*,
DBVPG_3236/DBVPG_3855) against the existing MS presence data for the best
rhodotorulic acid candidate (row 2190). **Result: both gene-negative
R. mucilaginosa strains still show strong, comparable-to-positive-control
MS signal** -- does not confirm a simple 1:1 gene-presence <->
compound-presence dependency. Checked BUSCO completeness for these 2
strains (busco_genome table in BFD.duckdb): both below panel average
(90.4%, 81.0% vs. 95.1% average), with DBVPG_3855 the LOWEST-completeness
genome in the entire 278-strain panel -- strongly suggestive this is a
gene-model dropout artifact of assembly quality, not true biological
absence or pathway redundancy, though not yet directly confirmed (would
need tblastn against the raw assembly rather than the predicted
proteome). Findings logged:
`.living/findings/siderophore-detectability-rhodotorulic-acid.md` F-002,
registry F-015.

**Consequences**: The genome side of the siderophore investigation is now
solid and reusable (real validated ortholog call, not a coarse domain
proxy) -- a meaningful upgrade matching the pattern already established
with Idea 3's custom pigment HMM panel. The specific 2-strain discordance
found should NOT be read as evidence against the gene-compound
relationship without first ruling out the assembly-quality explanation
(tblastn check, not yet done). A panel-wide statistical test is likely
underpowered regardless (only 2-3/278 strains negative either way).

**Tags**: siderophore, rhodotorulic-acid, nrps, diamond-ortholog-search, busco-completeness, assembly-artifact

## 2026-08-16: BUSCO<90 discordance accepted as assembly artifact; NRPS candidate multifasta/alignment/tree built

**Context**: PI reviewed the DBVPG_3236/DBVPG_3855 gene-negative-but-
MS-positive discordance from the prior entry and decided: ignore strains
with BUSCO completeness <90%, attribute the discordance to incomplete
assembly (most *R. mucilaginosa* strains do carry the ortholog), move on
-- no tblastn-vs-raw-assembly confirmation requested. Then asked for a
multifasta of the best candidate rhodotorulic-acid NRPS orthologs across
the panel, for alignment and phylogeny building.

**Decision/outcome**: Built `analysis/scripts/siderophore_nrps_build_multifasta.py`
-- pulls each strain's best-hit ortholog protein (from the diamond search's
per-strain summary) via BFD.duckdb's gene_proteins table, filtered to
confirmed ortholog AND BUSCO completeness >=90%. 10 strains dropped on
the BUSCO filter (not just the 2 originally flagged -- 8 more across
other species also fall below 90% among the 275 ortholog-positive set).
265 candidates + the PI-supplied reference = 266 sequences ->
`outputs/RA_NRPS_candidates.faa`. Aligned with `mafft --auto` (~15s) ->
`RA_NRPS_candidates.aln.fa`, first-pass tree with `FastTree` (~5s, CAT
approximation, no bootstrap) -> `RA_NRPS_candidates.tree.nwk`. Findings
logged: `.living/findings/siderophore-detectability-rhodotorulic-acid.md`
F-003.

**Consequences**: Real alignment/tree data now exists for asking whether
the NRPS gene tree tracks the species tree or shows topology anomalies
(HGT, duplication, lineage-specific acceleration) -- not yet interpreted,
that's the natural next step if the PI wants to pursue it. If this
becomes a figure, the current FastTree pass (no bootstrap) should be
upgraded to a proper ML tree with support values (e.g. IQ-TREE).

**Tags**: siderophore, rhodotorulic-acid, nrps, phylogeny, busco-filter, mafft, fasttree

## 2026-08-16: NRPS gene tree vs. species tree comparison, rendered figure

**Context**: PI asked whether the rhodotorulic-acid NRPS gene tree
(built from the prior entry's multifasta/alignment) tracks the species
tree, then separately asked for the tree to be drawn and saved as
PDF/PNG and linked into the results doc.

**Decision/outcome**: Built `analysis/scripts/siderophore_nrps_tree_species_comparison.py`
(per-species monophyly check via MRCA/terminal-set comparison,
terminal-branch-length outlier detection) and
`analysis/scripts/siderophore_nrps_plot_tree.py` (species-colored
Bio.Phylo + matplotlib rendering). Result: the 3 largest species
(*R. mucilaginosa* n=194, *R. paludigena* n=16, *R. toruloides* n=8) come
back "not monophyletic," but this is a tree-RESOLUTION artifact, not
biology -- FastTree found only 61 unique sequence patterns among 266
sequences, i.e. the gene is essentially invariant across most of the
genus at the protein level, visually confirmed as one shallow
near-zero-branch-length block in the rendered tree. The well-resolved
half of the tree (species with genuinely distinct sequences) forms clean
monophyletic clades matching expectations for simple vertical
inheritance -- no strong evidence of HGT. One genuine, unexplained
anomaly: *R. evergladensis* DBVPG_7922 has an outlier-long terminal
branch despite a clean high-confidence ortholog hit (not a bad gene
model) -- flagged as a real accelerated-evolution lead, not investigated
further. Figure delivered to PI and linked into
`phase_siderophore/RESULTS.md` (Step 7). Findings logged:
`.living/findings/siderophore-detectability-rhodotorulic-acid.md` F-004,
registry F-016.

**Consequences**: This closes out the phylogenetic side of the
siderophore investigation for now -- no further action proposed unless
the PI wants a bootstrap-supported ML tree (IQ-TREE) or wants to chase
the *R. evergladensis* branch-length anomaly.

**Tags**: siderophore, rhodotorulic-acid, nrps, phylogeny, gene-tree-species-tree, conserved-gene
