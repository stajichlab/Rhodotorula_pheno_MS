# Candidate Carotenoid Pathway Genotyping Across the 276-Strain BFD Panel

## Persona
**Fungal Genomics / Comparative Genomics Bioinformatician** — I think in terms of pathway completeness, gene copy number, and loss-of-function variants, not whole-genome fishing expeditions; if you already suspect *which* pathway makes the pigment, go get that pathway's genes directly instead of waiting on metabolomics to find it for you.

## Motivation
The whole project has treated "genome" as a downstream validation step (Phase 5) rather than a primary lens, even though Rhodotorula pigmentation is a textbook case of a short, well-characterized fungal/yeast carotenoid pathway (CrtYB-type bifunctional phytoene synthase/lycopene cyclase, CrtI phytoene desaturase, CrtS/CrtC-type downstream oxidases in *Xanthophyllomyces*/*Rhodotorula* torularhodin-producing species). A comparative genomicist's instinct when a phenotype is null against an unbiased metabolome screen is to ask whether the causal genotype is even variable at the sequence/copy-number level in the pathway believed to produce the phenotype, rather than assuming the pathway is invariant and searching everywhere else.

## Connection to Existing Data
- BFD already has full Pfam annotation (~5000 IDs, 3.27M hits) across 276-278 strains — the relevant Pfam domains (PF00494 squalene/phytoene synthase, PF01593 FAD-binding monooxygenase covering phytoene desaturase-like enzymes, PF04116 FA_hydroxylase-like for potential downstream oxidases) can be pulled directly with no new annotation run.
- The strain-level phylogenomic tree (fungi_odb10 BUSCO/PHYling/FastTree) already used for PGLS on the color/metabolome side (`analysis/YPD/color_shape_growth/`, `analysis/copper/`) is reusable verbatim for a genotype-vs-phenotype PGLS/phylogenetic-signal test here — no new tree needed.
- CIELAB a*/C* from the `control_90_110` canonical phenotype source is already cleaned and strain-matched to BFD strain IDs from prior work.
- Critically, this sidesteps the metabolome entirely — 5 independent null results (whole-panel, within-species x2, Lasso, ANOVA) all tested color against **detected/annotated metabolite abundance**, which is bottlenecked by the 4.8% SIRIUS annotation rate. Genome sequence doesn't have that annotation ceiling.

## Approach
1. Build a curated candidate gene set from characterized *Rhodotorula*/*Xanthophyllomyces dendrorhous*/*Sporidiobolus* carotenogenesis literature (CrtYB, CrtI, CrtS/CYP450 oxidase, CrtC-like) and pull matching Pfam/InterPro domain hits per strain from existing BFD tables — no new domain scan needed if Pfam IDs already exist; otherwise run a targeted HMMER/InterProScan pass restricted to these ~4-6 gene models only (cheap, hours not days).
2. Tabulate per-strain presence/absence and copy number (paralog count) for each candidate gene across all 276-278 strains, flagging species-level vs strain-level (within-*mucilaginosa*) variation separately, since the existing null results show *mucilaginosa* alone spans nearly the full color range.
3. Extract and align the actual CDS/protein sequences for each candidate gene per strain (from the same assemblies BFD annotated) to screen for premature stop codons, frameshifts, or catalytic-residue substitutions (e.g., known catalytic Cys/His in phytoene desaturase) that would silently knock out function despite gene presence — presence/absence alone is insufficient for a single-copy essential pathway.
4. Test copy number, PAV, and (if enough protein variation exists) a simple non-synonymous substitution burden score against CIELAB a*/C* using the same PGLS framework already built for color-vs-metabolome (`02_pgls_analysis.R` pattern), using the existing pruned tree to control phylogenetic non-independence.
5. If within-*mucilaginosa* PGLS is underpowered due to low sequence diversity in a single-copy gene, fall back to a phylogenetically-aware ANOVA of color by haplotype/genotype cluster (mirroring the pattern-group test already run on the metabolome side) so the same non-parametric machinery is reused.

## Expected Insights
Either (a) a genuine positive: copy number or specific catalytic-residue variants in CrtYB/CrtI/downstream oxidase genes track with a*/C*, which would be the first real color-linked genomic signal in the project and would also motivate re-examining the null carotenoid metabolome features (only 3-4 SIRIUS carotenoid-class calls) against genotype directly; or (b) a informative negative: the pathway is genomically invariant (single copy, intact, near-identical sequence) across the panel, which would strongly argue the color variation is transcriptional/regulatory/post-translational rather than encoded in pathway-gene sequence — redirecting future work toward promoter/expression-level analysis instead of more genome or metabolome fishing.

## Feasibility
- Effort: Low
- Data ready: Mostly (BFD Pfam tables and tree exist; needs a short curated candidate-gene list and possibly one narrow HMMER pass)
- Methods available: Standard tools (HMMER/Pfam lookup, existing PGLS scripts, BLAST for CDS extraction)
- Key risk: Rhodotorula's specific carotenogenesis gene models are not as well curated as *X. dendrorhous*'s; if Pfam-level domain hits are too coarse (multi-gene-family Pfam IDs like PF00494 hit many non-carotenoid isoprenoid synthases), the candidate set will need manual curation via BLAST/reciprocal-best-hit against a reference CrtYB/CrtI protein rather than pure Pfam ID filtering.

---

# Genome-Wide Orthogroup PAV/CNV Scan Against Color, Independent of Any Pathway Hypothesis

## Persona
**Fungal Genomics / Comparative Genomics Bioinformatician** — I know that gene-family expansion/contraction (CAFE5-style) and orthogroup presence/absence variation routinely explain trait variation in fungi even when the causal gene isn't in the textbook pathway, so I'd want an unbiased genome-wide genotype-phenotype scan as a complement to (not instead of) the candidate-gene approach.

## Motivation
Five independent metabolome-side null results (whole-panel, 2 within-species, Lasso, ANOVA) plus a real unexplained cell-fraction/colony-area confound suggest the causal biology may not be a simple detected/annotated small-molecule signal at all — it could be a structural genomic difference (gene loss, pathway paralog expansion, a regulatory gene family) that never shows up as a measurable LC-MS2 feature, or shows up as a feature in the 95.2% of features SIRIUS never annotated. A comparative genomicist's response to "the targeted small-molecule search failed five times" is to stop assuming the causal locus is even a metabolite-producing enzyme and instead ask which orthogroups/gene families vary in copy number across the phylogeny and correlate that variation, unbiased, against the phenotype.

## Connection to Existing Data
- BFD's Pfam (~5000 IDs), MEROPS (peptidases), and CAZy (carbohydrate-active enzyme) tables already give per-strain functional-domain content across all 276-278 strains — this is effectively a pre-computed orthogroup-like matrix that can be mined directly without running OrthoFinder from scratch, though true single-copy orthogroup calls would still benefit from an OrthoFinder pass on the existing proteomes.
- The strain-level PHYling/FastTree phylogeny is already built and already used for PGLS elsewhere in this repo, making a CAFE5-style gene-family evolution analysis (which needs an ultrametric or at least topologically-correct tree) close to plug-and-play.
- The dataset's extreme imbalance (206-216 of ~300-310 strains being *R. mucilaginosa*) is exactly the setting where a genome-wide domain-content scan is *more* robust than metabolome Lasso, because Pfam/CAZy/MEROPS hit counts are much lower-dimensional and less noisy per-strain than 10,949 compound groups, and don't depend on annotation completeness the way metabolite ID does.
- Directly complements the real side finding that colony area confounds cell-fraction metabolite abundance: a genome-wide domain-content scan run against colony area as a second phenotype (in parallel with color) could help explain that confound structurally (e.g., growth-rate-associated gene family differences) rather than leaving it as an unexplained artifact.

## Approach
1. Run OrthoFinder (or reuse Pfam-domain-count matrices as a fast proxy) across all 276-278 BFD proteomes to build a per-strain gene-family copy-number matrix; this is standard, off-the-shelf, and the proteomes/annotations already exist in BFD.
2. Filter to gene families/Pfam domains with non-trivial variance across strains (drop invariant core-eukaryotic families) and correct for genome assembly quality/BUSCO completeness per strain (already in BFD) to avoid spurious PAV calls driven by incomplete assemblies rather than true gene loss — a defensive-analysis step given `robust-analysis` conventions already active in this repo.
3. Run a phylogenetically-corrected association scan (PGLS or phylogenetic ANOVA per gene family, same tree/framework as `02_pgls_analysis.R`) between each variable gene family's copy number and CIELAB a*/C*, with explicit multiple-testing correction (FDR) given the number of families tested — analogous in spirit to the whole-panel + within-species block-permutation framework already used for metabolome features.
4. Separately run the same scan against colony area, to directly chase the unexplained cell-fraction/area confound described in the metabolome results — a two-birds move that reuses the identical pipeline.
5. For any gene family/orthogroup that survives FDR correction and phylogenetic control, pull its Pfam/InterPro annotation and cross-reference against the candidate carotenoid list from Idea 1 and against MEROPS/CAZy functional categories to interpret biologically, and flag it for a targeted BGC-adjacency check (are the significant genes physically clustered in the assembly, suggestive of a secondary-metabolite gene cluster antiSMASH hasn't been run to find yet).

## Expected Insights
A ranked, FDR-corrected list of gene families whose copy number/presence tracks color and/or colony area independent of any pigment-pathway hypothesis, which either (a) recovers the expected carotenoid genes from Idea 1 as independent validation, (b) surfaces an unexpected gene family (e.g., a CAZy or MEROPS family, or a transcription-factor Pfam domain) that reframes the causal hypothesis entirely, or (c) comes back null across the full annotated gene space, which — combined with 5 metabolome nulls — would be strong combined evidence that colony color in this panel is not encoded in gross gene content at all and is more likely driven by fine-grained regulatory/epigenetic variation, redirecting the project away from more genotype-content fishing.

## Feasibility
- Effort: High
- Data ready: Mostly (proteomes and Pfam/CAZy/MEROPS tables exist in BFD; needs OrthoFinder run and BUSCO-completeness-aware filtering, which is new compute but standard)
- Methods available: Standard tools (OrthoFinder, CAFE5 optional for formal gene-family evolution rates, existing PGLS scripts adapted to gene-family copy number instead of metabolite abundance)
- Key risk: Multiple-testing burden across thousands of gene families with an extremely imbalanced phylogeny (206+/278 strains in one species) risks either no power (harsh FDR correction) or confounding by assembly-quality artifacts masquerading as true PAV — must explicitly test whether significant hits are assembly-completeness-associated as a null check before trusting any result, per the repo's active robust-analysis convention.
