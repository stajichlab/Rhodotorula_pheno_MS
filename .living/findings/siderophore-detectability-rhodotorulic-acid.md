---
topic: siderophore-detectability-rhodotorulic-acid
description: Whether rhodotorulic acid and related iron-sequestration chemistry is detectable in the untargeted LC-MS2 data, its strain-level presence/absence, and whether the known NRPS biosynthetic gene can be genotyped from current genome data.
created: 2026-08-16
last_updated: 2026-08-16
status: active
---

# Siderophore (rhodotorulic acid) detectability and genome cross-reference

## F-001: Rhodotorulic acid is plausibly detected and near-universally present, but the genome-side NRPS screen is not yet discriminating enough to cross-reference against it
**Status:** preliminary
**Claim:** Exact-mass search (20 ppm, 4 compounds x 4 positive-mode
adducts, same design as the carotenoid-pathway search) found 29 raw
feature matches, including a high-confidence rhodotorulic acid candidate
(row 2190, [M+NH4]+, -2.96 ppm, 164 total scans -- the highest-intensity
match in the entire search). Strain-level presence (raw abundance >0)
is ~99% across all 17 species, especially in the supernatant fraction --
directionally consistent with rhodotorulic acid's known role as a
broadly-conserved, secreted, genus-wide siderophore, but the threshold
used (bare nonzero abundance) does not reliably distinguish real
low-level signal from background, so this should be read as "detectable
at some level in nearly every strain," not "confirmed present/absent."
A coarse genome-side Pfam screen (ornithine-hydroxylase PF00743 +
NRPS-domain co-occurrence, standing in for the real reference sequence
the PI has not yet supplied) found the hydroxylase proxy universal
(278/278 strains, uninformative) and the expected 2-adenylation-domain
NRPS architecture in only 1 genome of 278 -- almost certainly a draft-
assembly gene-model fragmentation artifact (large NRPS genes commonly
split across predicted gene models at scaffold gaps), not real absence
in 277/278 strains of a genus-conserved trait.
**Implications:** Neither side of the PI's proposed cross-reference
("condition whether the gene is found in strains that do not show
evidence of the compound") is currently usable: MS presence has no
discriminating variance at this threshold, and the genome screen cannot
reliably call NRPS module architecture from current gene models. The
path forward needs (a) the actual reference sequence/accession for the
characterized rhodotorulic-acid NRPS -- PI has this -- to enable a real
ortholog search (same upgrade Idea 3 took with its custom pigment HMM
panel), and (b) a tighter, intensity-based MS presence threshold.
**Tags:** siderophore, rhodotorulic-acid, nrps, iron-sequestration, needs-reference-sequence, gene-model-fragmentation

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Mass search + presence/absence + coarse Pfam NRPS screen (siderophore_*.py) | EB raw feature table (53,040 features) + BFD.duckdb pfam table (278 strains) | Rhodotorula_pheno_MS | 29 mass matches; ~99% strain presence (uninformative threshold); genome screen inconclusive (1/278 gene-model artifact) | refines |

### Open Questions
- Does MS2 fragmentation confirm row 2190 as rhodotorulic acid
  specifically (vs. an isobaric alternative), mirroring the carotenoid
  candidates' fragmentation check?
- Does an intensity/scan-count-thresholded presence call (rather than
  bare nonzero) produce a genuinely variable strain-level pattern worth
  testing against anything panel-wide?

## F-002: Real ortholog search (PI-supplied reference sequence) works cleanly, but the one direct genotype↔phenotype cross-reference it enables does NOT confirm the hypothesized link
**Status:** preliminary
**Claim:** PI supplied the actual reference sequence for the
rhodotorulic-acid NRPS (`RA_NRPS.fa`, protein F2DD6D01_006956-T1 from
*R. kratochvilovae* Y14, antiSMASH-confirmed NRPS cluster gene with an
adjacent biosynthetic-additional smCOG gene). `diamond blastp` of this
single query against all 278 BFD-panel proteomes (2,188,032 proteins)
gives a clean, strongly bimodal identity distribution (303 hits <30%
identity = noise/generic NRPS-domain cross-hits, vs. 305 hits >=60%
identity = real orthologs), with the 3 in-panel *R. kratochvilovae*
strains landing exactly as expected as a built-in positive control
(96.4-99.8% identity, ~100% coverage). This is a MUCH cleaner,
genuinely discriminating result than the coarse Pfam screen it
supersedes (F-001) -- 275/278 strains confirmed positive, only 3
negative (1 outgroup + 2 specific *R. mucilaginosa* strains, DBVPG_3236
and DBVPG_3855, out of 202 *R. mucilaginosa* strains total).
**Cross-referencing those 2 gene-negative strains against the MS
presence data (the actual test the PI proposed): the compound (row 2190,
the best rhodotorulic acid mass-search candidate) is still clearly
present in both, at abundances comparable to or higher than
ortholog-positive *R. mucilaginosa* controls.** This does not confirm a
simple 1:1 gene-presence <-> compound-presence dependency. Both
gene-negative strains have below-panel-average BUSCO completeness
(DBVPG_3236: 90.4%; DBVPG_3855: 81.0%, the LOWEST completeness in the
entire 278-strain panel, vs. panel average 95.1%), making an
assembly-quality gene-model dropout the most likely explanation for the
apparent "gene absence" rather than true biological gene loss or a
compensating paralog -- not yet distinguished.
**Implications:** The genome side of this investigation is now solid and
reusable (a real, validated ortholog call, not a coarse domain proxy).
The specific discordance found is NOT strong evidence against the
gene-compound link -- it's much more likely an artifact of 2 specific
lower-quality assemblies -- but this needs to be checked (tblastn against
raw assembly, not just predicted proteome) before either concluding
discordance or re-attempting a panel-wide statistical test.
**Tags:** siderophore, rhodotorulic-acid, nrps, diamond-ortholog-search, busco-completeness, assembly-artifact, gene-model-dropout

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Diamond ortholog search + MS cross-reference (siderophore_nrps_diamond_search.py) | PI-supplied RA_NRPS.fa vs. 278 BFD proteomes; MS presence table from F-001 | Rhodotorula_pheno_MS | 275/278 strains ortholog-positive; both gene-negative *R. mucilaginosa* strains still show strong MS signal; both have below-average BUSCO completeness | refines |

### Open Questions
- Worth a panel-wide statistical test (gene presence/absence vs. compound
  abundance) once MS2-confirmed -- not yet attempted, would need only
  3/278 strains negative either way so is likely underpowered regardless.

**PI decision (2026-08-16): closed.** Strains with BUSCO completeness
<90% are excluded going forward and the DBVPG_3236/DBVPG_3855 discordance
is attributed to incomplete assembly, not pursued as biology (most
*R. mucilaginosa* strains do carry the ortholog). No tblastn-vs-raw-
assembly confirmation was run -- this is an accepted working assumption,
not a proven mechanism, should it matter again later.

## F-003: Candidate-ortholog multifasta, alignment, and first-pass phylogeny built (265 strains + reference)
**Status:** preliminary
**Claim:** `analysis/scripts/siderophore_nrps_build_multifasta.py` pulled
the best-hit rhodotorulic-acid-NRPS ortholog protein sequence per strain
(from F-002's diamond search), filtered to confirmed ortholog AND BUSCO
completeness >=90% (10 strains dropped on the BUSCO filter, not just the
2 originally flagged -- 8 more across other species also fall below 90%
among the 275 ortholog-positive set). 265 candidates + the PI-supplied
reference (F2DD6D01_006956-T1) = 266 sequences, aligned with `mafft
--auto` and a first-pass tree built with `FastTree` (CAT approximation,
no bootstrap -- fine for an initial topology look, not publication-grade).
FastTree found only 61 unique sequence patterns among 266 sequences --
expected for a conserved single-copy gene (many near-identical
*R. mucilaginosa* copies), not a run error.
**Implications:** This is now real, ready-to-inspect data for asking
whether the NRPS gene tree tracks the species tree (simple vertical
inheritance) or shows any topology anomalies (horizontal transfer,
gene duplication/loss, unusually fast/slow evolution in particular
lineages) -- none of that interpretation has been done yet, this entry
just records that the inputs exist.
**Tags:** siderophore, rhodotorulic-acid, nrps, phylogeny, mafft, fasttree

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Multifasta + mafft + FastTree (siderophore_nrps_build_multifasta.py) | F-002's diamond ortholog set, BUSCO-filtered | Rhodotorula_pheno_MS | 266-sequence alignment and tree built, not yet interpreted | supports |

### Open Questions
- Is a proper ML tree with bootstrap support (e.g. IQ-TREE) worth
  building if this becomes a figure, vs. the current fast FastTree pass?

## F-004: Gene tree broadly consistent with vertical inheritance, but most of it is unresolvable — the gene is extremely conserved
**Status:** preliminary
**Claim:** `analysis/scripts/siderophore_nrps_tree_species_comparison.py`
checked per-species monophyly and terminal branch-length outliers on the
266-tip candidate gene tree (F-003). FastTree found only 61 unique
sequence patterns among 266 sequences -- the gene is conserved enough
that *R. mucilaginosa* (n=194), *R. paludigena* (n=16), and
*R. toruloides* (n=8) collapse into one shallow, near-zero-branch-length
polytomy block with a few sequences from other species interleaved,
which the automated monophyly check flags as "not monophyletic." This is
almost certainly a tree-resolution artifact (insufficient variable sites
at the protein level), not evidence of HGT or contamination -- visually
confirmed in the rendered tree
(`outputs/RA_NRPS_candidates.tree.{png,pdf}`, see also
`.living/findings` topic RESULTS.md for the embedded figure). The
well-resolved lower half of the tree (species with genuinely distinct
sequences: *R. dairenensis*, *R. diobovata*, *R. graminis*,
*R. kratochvilovae*, *R. sp.* clades I/XI, *R. sphaerocarpa*) forms
clean, correctly-grouped monophyletic clades -- broadly consistent with
simple vertical inheritance, no strong anomaly. Two smaller-scale
exceptions flagged but not investigated further: *R. taiwanensis*
clusters with one *R. toruloides* tip, and *R. pacifica* (n=2) sits
inside the *R. mucilaginosa* block. Two terminal-branch-length outliers
found: *Pseudomicrostroma phylloplanum* (expected -- most divergent
outgroup) and **R. evergladensis DBVPG_7922** (bl=0.234, NOT trivially
explained -- its diamond hit is a clean 69.7%-identity/99.8%-coverage
ortholog, not a bad gene model, so this looks like genuine accelerated
evolution in this one strain's copy).
**Implications:** This gene is under strong purifying selection /
essentially invariant across most of the genus at the protein level --
consistent with it being a core, conserved biosynthetic function (matches
the near-universal MS presence and near-universal genome presence found
earlier in this investigation). Not useful for resolving fine-scale
strain relationships, but the *R. evergladensis* branch-length anomaly is
a genuine, unexplained lead if that species becomes relevant elsewhere.
**Tags:** siderophore, rhodotorulic-acid, nrps, phylogeny, gene-tree-species-tree, evergladensis, conserved-gene

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-16 | Monophyly + branch-length check (siderophore_nrps_tree_species_comparison.py) | F-003's gene tree | Rhodotorula_pheno_MS | Well-resolved species form clean clades; big-n species collapse into an unresolvable polytomy (resolution artifact); R. evergladensis flagged as a genuine branch-length outlier | refines |
