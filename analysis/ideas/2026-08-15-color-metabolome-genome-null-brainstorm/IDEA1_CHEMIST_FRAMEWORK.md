# Framework: Idea 1 — Targeted carotenoid chemistry (chemist persona)

Combines both chemist-persona ideas (`01_natural-products-chemist.md`)
into one phased framework: 1b is a cheap, no-new-data prerequisite check;
1a is the wet-lab escalation if 1b doesn't resolve the question.

## Why this track exists

Five statistical methods found zero color-metabolome association in the
untargeted LC-MS2 data. The chemist's core point: **that's ambiguous
between "no real biology" and "the method can't see the molecule."**
Carotenoids (torulene, torularhodin, β-carotene) are long, fully-conjugated
polyene chains — poor ESI ionization, prone to isomerization/degradation,
often better suited to APCI/APPI ionization or dedicated pigment
chromatography than the generic reversed-phase-ESI setup an untargeted
LC-MS2 pipeline uses. This track resolves that ambiguity before any more
statistical firepower gets spent on the existing untargeted data.

## Phase 1 (now): Predicted-mass targeted re-mining of existing raw data — no new data

**Goal**: check whether the real pigment signal is present in the raw MS
files but was orphaned/miscollapsed by the untargeted pipeline's own
feature-finding, deduplication, or SIRIUS annotation steps — cheapest
possible check, reuses files already on disk.

### Steps
1. **Build an exact-mass + MS2-fragment-ladder target list** for the full
   carotenogenesis pathway, not just the end products:
   - Phytoene, phytofluene, ζ-carotene, neurosporene, lycopene (the
     desaturation series *crtI* walks through)
   - γ-carotene, β-carotene, torulene (the cyclization/branch products)
   - Torularhodin (*crtS* product) and any known intermediate oxidation
     states
   - Common in-source fragments/adducts expected for this compound class
     specifically ([M+H-H2O]+, loss of terminal ring fragments, etc. —
     literature-derived, not the generic adduct list EB already tried)
2. **Search the raw/aligned feature tables directly by exact mass + RT
   plausibility** (not by SIRIUS annotation) against this target list —
   i.e., don't trust the pipeline's own compound-calling for this specific
   search; do an independent, targeted mass lookup.
3. **Cross-reference hits against the dedup groups**
   (`analysis/linked_data/ms_feature_dedup_groups.csv`) to see which
   `dedup_group_id` each mass-list hit belongs to, and whether that group
   was one of the ~10,949 tested in Phase 2 — if a real carotenoid feature
   got folded into an ISF group or lost in the ~4.9% un-merged
   cross-adduct case (see the dedup script's documented limitation), that
   would explain a false null directly.
4. **If a plausible carotenoid-pathway mass hit is found**: pull its
   MS2 spectrum and manually/semi-manually assess fragmentation pattern
   plausibility (polyene neutral losses, etc.) before trusting the mass
   match alone — exact mass alone is not confirmation at this resolution.
5. **Report either way**: "found candidate features, here's where they
   sit in the existing analysis" (informs a Phase 2 rerun using the
   right features) OR "no plausible mass matches anywhere in the raw
   data" (real evidence the untargeted method structurally can't see
   this pathway, strengthening the case for Phase 2 escalation below).

### What's needed
- Reference exact masses / adduct rules for the pathway intermediates
  (literature lookup, not project data — can be compiled without any
  BFD/database dependency).
- Access to the raw/full feature table with m/z + RT
  (`aligned_features_ms2.csv` in the EB nf_output, already used this
  session for the dedup work — no new data acquisition).
- No wet lab, no new instrument time. **This can start immediately,
  independent of the BFD rebuild.**

## Phase 2 (contingent on Phase 1's outcome): Targeted wet-lab extraction + APCI/APPI-MS

**Goal**: if Phase 1 finds nothing even at the raw-mass level, that's real
evidence the untargeted method can't see this chemistry — escalate to a
purpose-built extraction/ionization method that can.

### Steps
1. **Select a stratified strain subset** spanning the observed color
   range (not all ~300 strains — this is a targeted, smaller-n
   confirmatory experiment) — use the existing `strain_phenotype_table.csv`
   a\*/C\* values to pick e.g. 10-20 strains spanning low to high a\*,
   including some of the "convergent candidate" species from Phase 1's
   `convergent_color_candidates.csv` (*R. glutinis*, etc.) if archived
   material/viable cultures exist for them.
2. **Confirm sample/culture viability** — check whether archived pellet
   material from the original EB extraction is still usable (degradation
   risk given carotenoids are photo/thermolabile) or whether fresh
   fermentation runs are needed for the selected subset.
3. **Standard carotenoid extraction protocol** (organic solvent extraction
   optimized for polyene stability — literature methods exist for
   *Rhodotorula*/*Xanthophyllomyces* specifically) with authentic
   standards (torulene, torularhodin, β-carotene if commercially available)
   run alongside for retention-time and response-factor calibration.
4. **C30 reversed-phase LC (the standard carotenoid-separation column
   chemistry, better isomer resolution than generic C18) + APCI or APPI
   ionization** (not the generic ESI the untargeted pipeline used) —
   either via instrument access already available to the lab or a core
   facility/collaborator.
5. **Quantify torulene/torularhodin/β-carotene (and ratios between them)
   per strain**, then correlate directly against CIELAB a\*/C\* for that
   same strain subset — a small, targeted, confirmatory test, not another
   FDR-corrected genome-wide scan.

### What's needed
- PI decision on strain subset and whether archived material is usable
  vs. needs fresh growth.
- Authentic carotenoid standards (procurement).
- APCI/APPI-capable LC-MS instrument access (may need a core facility or
  collaborator if not already available in-house).
- Realistic timeline: weeks, not days, given wet-lab + instrument
  scheduling — this is the "escalation" tier, not the near-term deliverable.

## Decision gate between phases

Phase 2 should only proceed if Phase 1 comes back genuinely empty (no
plausible carotenoid-pathway mass matches anywhere in the raw data, even
outside the already-tested dedup groups). If Phase 1 *does* find candidate
features, the right next step is re-running Phase 2's statistical tests
(`phase2_color_metabolome_association.py` etc.) specifically against those
recovered features/groups — cheap, and doesn't require any new wet lab.

## Status
Phase 1 not yet started (needs the exact-mass target list compiled first —
literature task, can begin immediately). Phase 2 not started, contingent
on Phase 1's outcome per the decision gate above.
