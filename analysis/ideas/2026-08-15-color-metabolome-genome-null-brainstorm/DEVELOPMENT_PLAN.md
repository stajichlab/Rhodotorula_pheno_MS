# Development plan: Ideas 1, 3, 5, 6 (2026-08-15)

Follow-up to `00_index.md`. PI wants to pursue most of the 14 brainstormed
ideas and asked specifically for: (a) a framework for Idea 1 (chemist), (b)
to start implementing Idea 3 (candidate carotenoid-pathway genotyping), (c)
whether Idea 3 needs a SNP-calling strategy, (d) a strategy for Idea 5
(convergence test), and (e) what infrastructure Ideas 3+6 share.

**Blocking dependency**: the BFD genome database is currently being
rebuilt by a separate nextflow pipeline (PI, 2026-08-15) and will land at
`BFD/db/BFD.duckdb` when finished. Nothing below queries the database
directly — everything genome-side is scaffolded/ready-to-run, not
executed, until that lands. Do not query any interim/stale copy in the
meantime (confirmed with PI).

---

## Part A — SNP calling: do we need it?

**Short answer: not for Idea 3 as scoped. Only for the more ambitious tail
of Idea 6.** Two genuinely different tiers of genomic resolution are in
play, and they have very different infrastructure costs:

### Tier 1 — gene-level presence/copy-number/LoF (Idea 3, no SNP calling needed)
Every BFD strain already has its own independent genome assembly *and*
gene annotation (Pfam/MEROPS/CAZy hits, called proteins in `gene_proteins`).
This means:
- **Copy number / presence-absence** of a candidate gene family = a
  `COUNT(*)` over the relevant Pfam ID(s) per strain. Already-computed
  data, zero new pipeline needed.
- **Loss-of-function / catalytic-residue variants** = pull each strain's
  already-called protein sequence(s) for the candidate gene, multiple-
  sequence-align them (MAFFT) across all 276-278 strains, and inspect the
  alignment directly for premature stop codons, frameshift-pattern
  truncations, or substitutions at known catalytic residues (from
  published CrtYB/CrtI/CrtS structural or mutagenesis literature). This
  works on the **assembly-derived protein calls that already exist** —
  no raw-read remapping, no reference-genome alignment, no variant
  caller (GATK/freebayes/etc.) required.

The only new bioinformatic step Tier 1 needs is **ortholog identification**
(BLAST/reciprocal-best-hit against a curated reference protein set),
because Pfam IDs alone are too coarse (e.g. PF00494 "SQS_PSY" covers both
phytoene synthase/CrtYB-type enzymes AND unrelated squalene synthases;
PF01593 "Amino_oxidase" covers CrtI-type phytoene desaturases AND many
unrelated FAD-dependent oxidoreductases). That curation step is described
in Part B below.

### Tier 2 — fine-grained/genome-wide variant GWAS (Idea 6's ambitious tail, DOES need real variant calling)
If the goal becomes "test every variable site genome-wide against color,"
not just "does this candidate gene show LoF variants," gene-level
presence/absence is too coarse a unit — you'd need an actual SNP/variant
matrix. That requires infrastructure this project does NOT currently have:
- A common coordinate system across strains this diverged (multiple
  species, not just strains of one species) — either a high-quality
  reference genome + whole-genome alignment (e.g. `mummer`/`nucmer` or
  `minimap2` pairwise, or a graph-genome/pangenome approach given how
  phylogenetically broad this panel is), or a k-mer-based
  reference-free approach (e.g. kmer-GWAS / pyseer with unitigs).
- Raw sequencing reads for every strain (assemblies alone aren't enough
  for standard GATK/freebayes-style calling) — need to confirm these
  are archived/accessible for all 276-278 BFD strains, not just the
  assemblies.
- A phylogenetic mixed-model or linear-mixed-model GWAS framework
  (e.g. `GEMMA`, `pyseer` with a kinship matrix, or a custom PGLS-per-
  variant loop) to correct for population structure at variant scale —
  the species-level PGLS infrastructure already built in this project
  handles gene-level predictors fine but would need adaptation (or a
  purpose-built tool) for a variant matrix with potentially tens of
  thousands of sites.

**Recommendation**: build Tier 1 first (Idea 3, ready now). Only invest in
Tier 2 infrastructure if Tier 1 turns up a real candidate gene/region
worth fine-mapping, or if Tier 1 comes back null and a hypothesis-free
genome-wide scan (Idea 3b, "Genome-Wide Orthogroup PAV/CNV Scan") also
comes back null — at that point Tier 2 becomes the natural escalation, not
before. This mirrors the project's existing primary-hypothesis-first,
exploratory-second convention.

---

## Part B — Idea 3 implementation plan (candidate carotenoid pathway genotyping)

### Step 1: Curate the candidate gene/Pfam list (can do now, no DB needed)
Rhodotorula's carotenogenesis pathway is homologous to the well-studied
*Xanthophyllomyces dendrorhous* (*Phaffia rhodozyma*) pathway (both
Sporidiobolales). Candidate genes, in priority order:

| Gene | Function | Pfam ID(s) (coarse, needs curation) | Reference protein to seed ortholog search |
|---|---|---|---|
| *crtYB* | Bifunctional phytoene synthase / lycopene cyclase | PF00494 (SQS_PSY) | *X. dendrorhous* CrtYB (UniProt) |
| *crtI* | Phytoene desaturase (builds conjugated polyene chain -> color intensity) | PF01593 (Amino_oxidase) | *X. dendrorhous* CrtI |
| *crtS* | Cytochrome P450, torulene -> torularhodin (the step most likely to differentiate "orange" from "red/darker orange") | PF00067 (p450) | *X. dendrorhous* CrtS / CYP60related |
| *crtR* | CrtS's cognate cytochrome P450 reductase (torularhodin pathway) | PF00970 + PF00258 (FAD-binding + flavodoxin domains) | *X. dendrorhous* CrtR |
| HMGR | HMG-CoA reductase (upstream isoprenoid precursor supply, rate-limiting in many fungi) | PF00368 | generic fungal HMGR |
| GGPPS | Geranylgeranyl pyrophosphate synthase (precursor supply) | PF00348 | generic fungal GGPPS |

*crtYB*, *crtI*, and *crtS* are the highest-priority three — *crtS*
specifically is the mechanistically best candidate for explaining
*variation in hue/darkness* (not just presence/absence of color) since it
controls the torulene-to-torularhodin branch point, i.e. exactly the kind
of "different carotenoid composition" hypothesis in the strategy doc's
Biological Context section.

**Action needed from PI/further work**: pull the actual reference protein
sequences (UniProt/NCBI accessions) for *X. dendrorhous* CrtYB/CrtI/CrtS/
CrtR to seed the BLAST search — I can look these up once network/database
access for reference sequences is confirmed as in-scope, or PI can supply
them directly if already on hand from prior literature review.

### Step 2: Ortholog identification (needs BFD protein sequences, i.e. needs the DB)
1. Extract candidate proteins from `gene_proteins` (or the underlying
   `.proteins.fa` files in `BFD/input_all/pep/`) for all 276-278 strains,
   pre-filtered by the coarse Pfam IDs above (fast, cheap filter).
2. Reciprocal-best-hit BLAST against the *X. dendrorhous* reference set to
   confirm true orthology (not just shared Pfam domain) and assign each
   strain's candidate hit to a gene name (*crtYB*/*crtI*/*crtS*/*crtR*).
3. Flag strains with zero orthology hits (candidate gene absent/not
   assembled/too divergent) vs. multiple hits (possible duplication —
   copy-number signal itself, worth testing).

### Step 3: Presence/copy-number test against color (Tier 1a)
- Species-level (primary) and strain-level (secondary/robustness) count
  of each candidate gene per strain, using the already-built
  `build_species_level_tables.py` infrastructure.
- PGLS: `a* ~ crtS_copy_number + busco_completeness + log(protein_count)`,
  same convention as the copper-AUC PGLS script
  (`analysis/copper/scripts/02_pgls_analysis.R`) — reuse that script's
  structure directly, swap the predictor.
- Negative-control calibration: same hard-gated decoy convention as Phase
  2 (colony area).

### Step 4: Sequence-level LoF/catalytic-residue test (Tier 1b)
- MAFFT MSA of each candidate gene's protein sequences across all strains
  with a confirmed ortholog.
- Programmatic scan for: premature stop codons (sequence truncated well
  before the alignment's consensus length), frameshift-pattern garbling,
  and substitutions at literature-known catalytic/functional residues
  (needs the specific residue numbers from *X. dendrorhous* structural/
  mutagenesis papers — to be sourced alongside the reference sequences in
  Step 1).
- Binary or ordinal "predicted functional / predicted LoF" call per
  strain per gene, tested against color the same way as Step 3.

### Script scaffold
`analysis/scripts/phase5_candidate_gene_genotyping.py` has been drafted
(see that file) with Steps 1-2's logic stubbed and clearly marked
`# BLOCKED: needs BFD.duckdb` at the exact query points — ready to fill in
and run the moment the database rebuild finishes.

---

## Part C — Strategy for Idea 5 (convergence-restricted association test)

This is a **design document, not yet implementable** — the evolutionary
biologist persona flagged it as High effort / inherently low power
(N=independent origins, not N=strains), so this should be scoped
carefully before committing engineering time, per that persona's own
honesty about feasibility.

### Step 1: Formal ancestral state reconstruction (not yet done)
The existing `convergent_color_test.R` (Phase 1) is a coarse heuristic
(above-mean orange_score AND phylogenetically-distant-from-focal-species).
A real convergence test needs:
- Fit an Ornstein-Uhlenbeck or Brownian-motion-with-regime-shift model
  (e.g. R packages `l1ou`, `bayou`, or `OUwie`) to orange_score/a* on the
  **species-level tree** (`species_tree.nwk`, already built).
- This identifies the ML/Bayesian-supported set of branches where a
  "regime shift" (color-gain event) most likely occurred — formally
  replacing the current heuristic's "above mean + far from R. dairenensis"
  proxy.

### Step 2: Define independent contrast pairs
- For each inferred regime-shift branch, pair the descendant "gained
  color" clade against its phylogenetically nearest "did not gain color"
  sister/outgroup clade — this is the correct evolutionary unit of
  replication (phylogenetically independent contrasts), not raw strain
  count.
- Given the current species tree has only 17-18 tips, expect **at most
  3-4 usable contrast pairs** (matching the existing heuristic's
  candidates: *R. glutinis*, *R. sp. clade I*, *R. diobovata* vs.
  *R. dairenensis*/*R. taiwanensis*/*R. sphaerocarpa*). This should be
  stated explicitly in any write-up as a fundamental power ceiling, not a
  fixable engineering problem — N=3-4 independent origins is the ceiling
  regardless of strain count.

### Step 3: Test whether independent origins share molecular correlates
- For metabolome: within each contrast pair, rank deduplicated compound
  groups by effect size (not full FDR — too few independent tests to
  support formal multiple-testing correction at N=3-4); check whether the
  same compound groups rank highly across *multiple independent* contrast
  pairs (a convergence signal a pooled/whole-panel test structurally
  cannot detect, since pooling averages across independent origins rather
  than looking for repeated signal).
- For genome: same design against Idea 3's candidate gene calls (natural
  integration point between Ideas 3 and 5) — do independent color-gain
  clades independently show the same *crtS* (or other candidate)
  variant/copy-number pattern?

### Recommendation
Treat this explicitly as **hypothesis-generation, not a confirmatory
test** (per the persona's own framing) — a small, ranked list of
candidates that appear in ≥2 independent contrasts is the deliverable,
not a p-value. Sequence this AFTER Idea 3's Tier 1 results exist, since
Step 3's genome side needs Idea 3's gene calls as input, and the metabolome
side is cheap to add once Step 1/2's contrast pairs are defined (which
doesn't need the database, so Steps 1-2 can start now against the
already-built species tree).

---

## Part D — Shared infrastructure for Ideas 3 + 6

Idea 6b ("stop treating species/clade as pure nuisance — genotype-to-color
path-specific reanalysis") and Idea 3 converge on the same underlying
need: **a genome-to-color association test that doesn't just block out
phylogeny as a nuisance, but treats it as a real information source.**

### What's already built and reusable
- Species-level collapse infrastructure (`build_species_level_tables.py`,
  `prune_species_tree.R`) — works for any strain-level predictor table,
  including gene copy-number/LoF calls once Idea 3 produces them.
- PGLS scripts (`analysis/copper/scripts/02_pgls_analysis.R` as a direct
  template) — Pagel's-lambda GLS, already handles the phylogenetic
  covariance structure properly rather than just "blocking it out,"
  which directly addresses Idea 6b's concern.
- Block-permutation infrastructure (`scripts/block_permutation.py`,
  extended this session in the Phase 2 scripts) for empirical/negative-
  control-calibrated significance.
- The species-level tree itself (`species_tree.nwk`).

### What's missing (needed before Idea 3/6 can run end-to-end)
1. **The rebuilt BFD database** (blocking, in progress).
2. **Reference protein sequences + catalytic-residue annotations** for
   the candidate carotenoid genes (Part B, Step 1) — literature lookup,
   not a data-pipeline gap.
3. **Ortholog-calling step** (Part B, Step 2) — not yet written; **checked
   2026-08-15: no `blastp` module on HPCC, but `diamond` (modules
   2.0.13-2.1.24) is available — use `diamond blastp` as the standard
   faster drop-in.** `module load diamond/2.1.24` (latest).
4. **MAFFT** for the MSA-based LoF scan (Part B, Step 4) — **checked
   2026-08-15: available, `module load mafft/7.505` (latest).**
5. **A species-mean vs. strain-mean sensitivity framework specific to
   genomic predictors** — the existing Species-Level Collapse procedure
   covers this in principle (Step 2's "prevalence" convention for binary
   features), but hasn't been exercised end-to-end on a real genomic
   predictor yet; Idea 3 will be the first real test of that machinery.
6. **(Tier 2 only, not needed yet)** — see Part A: raw reads + reference/
   pangenome alignment + a variant-scale mixed-model GWAS tool, only if
   Tier 1 motivates escalating.

### Recommended build order
1. ~~Confirm diamond/mafft availability~~ — **done 2026-08-15**: both
   available as HPCC modules (`diamond/2.1.24`, `mafft/7.505`).
2. Source reference sequences + catalytic residues for *crtYB*/*crtI*/
   *crtS* (literature/PI, doesn't need the DB) — not yet done.
3. The moment `BFD/db/BFD.duckdb` lands: run
   `phase5_candidate_gene_genotyping.py` Steps 1-2 (ortholog calling),
   then Steps 3-4 (presence/LoF tests against color).
4. Feed Idea 3's output into Idea 5's Step 3 (convergence check on
   genomic correlates) and into Idea 6b's genotype-to-color path-specific
   framing.
