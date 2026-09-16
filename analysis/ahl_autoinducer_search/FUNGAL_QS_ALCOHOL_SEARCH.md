# Fungal quorum-sensing alcohol search (farnesol/tyrosol/tryptophol/phenylethanol)

## Why this search exists

Follow-up to the biofilm-signaling question (2026-09-16): unlike AHLs and
the N-acyl-lipid classes searched so far (all of which turned out to be
cell-associated membrane lipids, not diffusible signals), farnesol and
tyrosol are the mechanistically best-established, verified fungal
density-dependent signals for the exact phenomenon the PI is asking
about — coordinating yeast-to-hyphal/biofilm morphogenesis. Verified this
session (not assumed): farnesol represses the yeast-to-hypha transition
in *Candida albicans* (Hornby et al. 2001, *Appl. Environ. Microbiol.*,
PMID 11425711); tyrosol does the reverse and biofilm cells secrete
substantially more of it than planktonic cells. **Every citation found
for this mechanism is Ascomycota (Candida-clade) — no Basidiomycota or
*Rhodotorula*-specific literature was found.** This search tests the
chemistry directly rather than assuming it does or doesn't apply here.

## Method

`analysis/scripts/fungal_qs_alcohol_mass_remining.py` — same exact-mass
re-mining design as the earlier searches, but a small, curated list of 5
specific named compounds (not a combinatorial homolog sweep): farnesol
(C15H26O), farnesoic acid (C15H24O2, a related farnesol-pathway oxidation
product also reported as QS-active in some studies, included as a bonus
target not verified to the same standard as farnesol/tyrosol this
session), tyrosol (C8H10O2), tryptophol (C10H11NO), phenylethanol
(C8H10O). 3 positive-mode adducts each, 15 targets total, 20 ppm.

## Result: clean null

**Tyrosol and phenylethanol: zero matches at any adduct.** Not detected
in this dataset at all, at 20 ppm.

**Tryptophol (9 matches), farnesol (11), farnesoic acid (9): matches
exist, but none show real structural corroboration.** Only one row
across all 29 matches has a SIRIUS formula exactly matching its target
(row 37002, farnesol [M+H]+) — and SIRIUS's own independent structure
call for that row is **"pentadeca-8,13-dien-2-one"** (class: "Fatty
alcohols"), a linear dienone with the same gross formula (C15H26O) as
farnesol but a completely different structural class (no terpenoid ring/
branching) — the same coincidental-isobaric-overlap pattern seen
throughout this entire investigation, arguing against, not for, farnesol
identity at that row.

**A clear tell**: 15 of the 29 raw matches (all the farnesol and
farnesoic-acid [M+Na]+ hits) are **the exact same feature rows** that
already matched in both the original AHL search and the N-acyl amino
acid search (row IDs 456, 536, 578, 2442, 3595, 3721, 3937, 8755, 13744,
23461, 29560, 30278, 34196, 42242, 44520 — verified by direct set
intersection against `ahl_mass_matches.csv` and
`nacyl_amino_acid_mass_matches.csv`). This is the signature of a generic,
recurring "hot" mass region (around m/z 245 and 259, specifically the
[M+Na]+ adduct) that coincidentally overlaps almost any reasonably-sized
target mass placed there — not a real, repeated confirmation of anything.
This is the same phenomenon the decoy-null analysis flagged earlier: some
regions of this feature table are just densely populated with real,
unrelated small-molecule chemistry.

## Bottom line

**No evidence for farnesol, tyrosol, tryptophol, or phenylethanol in this
dataset.** Two of the five targets are entirely absent; the other
three's matches are either isobaric coincidences with an unrelated
compound class or reuse of an already-flagged generic noisy mass region.
This does not prove the mechanism is absent from *Rhodotorula* — the
literature check found no prior Basidiomycota data either way — but this
project's own MS data gives no support for it, on the same terms every
other candidate in this investigation has been held to.

**What this means for the biofilm question**: the single most
mechanistically-motivated candidate signal class has come back null.
Combined with everything else found (AHLs null; the confirmed N-acyl-
arginine/lysine, ceramide, and Palmitoylethanolamide hits all being
cell-associated structural lipids, not diffusible signals), **this
investigation has not identified a plausible coordinating signal for
Rhodotorula biofilm formation** in the existing MS data. The two most
productive remaining directions, per the earlier synthesis, are (a)
testing whether sphingolipid-pathway activity correlates with
biofilm-relevant strain differences as a required-machinery hypothesis
rather than a signal hypothesis, and (b) the interkingdom-signaling
angle — testing whether *Rhodotorula*'s behavior responds to bacterial
lipids (AHLs, aminolipids) present in its environment, which is a
co-culture/exposure question, not a mass-search question.

## Files
- `analysis/scripts/fungal_qs_alcohol_mass_remining.py` — target-list builder + search
- `fungal_qs_alcohol_target_list.csv` — 15 target m/z values (5 compounds x 3 adducts)
- `fungal_qs_alcohol_mass_matches.csv` — 29 raw, unfiltered matches at 20 ppm
