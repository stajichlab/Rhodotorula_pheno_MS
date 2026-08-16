# Molecular-Network-Aware Compound Family Scores vs. Color

## Persona
**Metabolomics / Mass Spectrometry Specialist** — I don't trust per-feature statistics on untargeted LC-MS2 data until I've looked at the spectral network, because MS2 fragmentation similarity groups chemically related "dark matter" that no single feature-level test can see.

## Motivation
Every null result so far (whole-panel, within-species, Lasso, k-means/Kruskal-Wallis) operated on 10,949 individually deduplicated compound groups treated as statistically independent, exchangeable units. That is exactly the assumption GNPS-style molecular networking exists to break: structurally related analogs (in-source fragments, adduct variants, biosynthetic-family members) produce correlated but distinct feature intensities, diluting per-feature signal and inflating the multiple-testing burden. The project already computed `compound_network.graphml`, `network.graphml`, `compound_pairs.tsv`, `filtered_pairs.tsv`, and `merged_pairs.tsv` via EverythingBagel — a full GNPS-equivalent network — and it has literally never been opened for analysis. That's the first thing a metabolomics specialist would reach for after five null univariate/multivariate tests, because it converts 10,949 noisy features into a much smaller number of coherent chemical-family scores, which is both a multiple-testing fix and a biologically motivated aggregation (only 4.8% of features have any SIRIUS class label, but network topology doesn't require an ID to group chemically related signal).

## Connection to Existing Data
- Directly reuses `compound_network.graphml` / `compound_pairs.tsv` / `filtered_pairs.tsv` from the EB pipeline output — no new wet-lab or MS data needed.
- Directly addresses the stated gap: "this has NOT yet been used in any analysis in this project."
- Reduces the dimensionality problem that likely hurt the Lasso (negative CV R^2) by replacing 10,949 correlated features with ~50-500 network-cluster eigenfeatures.
- Reuses the already-validated pipeline architecture (phylogenetic block permutation, colony-area decoy) from tests 1-5, just applied to a different feature space, so the negative-control gating logic doesn't need to be re-derived.
- Also has a natural hook into the carotenoid signal: the 3-4 SIRIUS-confirmed carotenoid-class features and 73 "Terpenoids" NPC-pathway features can be used to *label* whichever network cluster(s) they fall in, turning an anonymous network module into "the carotenoid neighborhood" even without full structure ID — carotenoids are the single most plausible color-pigment class in a *Rhodotorula* system and deserve a targeted look independent of the panel-wide scan.

## Approach
1. Parse `compound_network.graphml` (or `filtered_pairs.tsv`, whichever has the QC-recommended edge threshold) and run community detection (e.g., Louvain/Leiden on cosine-similarity edges) to define molecular families / spectral clusters.
2. For each cluster, compute a family-level abundance summary per strain-fraction sample (e.g., sum or first-PC eigenfeature of the deduplicated compound-group intensities that fall in that cluster), separately for cell pellet and supernatant fractions given the known colony-area/fraction confound.
3. Flag and separately report the cluster(s) containing the SIRIUS carotenoid/Terpenoids-pathway features as an a priori "pigment-relevant" hypothesis set, tested with tighter correction (small family, one-shot test) rather than buried in the panel-wide scan.
4. Re-run the existing validated statistical stack (Spearman/PGLS-style phylogenetic block permutation at species level, within-mucilaginosa repeat, colony-area decoy as hard gate) on cluster-level scores instead of individual compound groups, at both panel-wide and within-R.-mucilaginosa scope.
5. Cross-tabulate cluster membership against node degree / connectivity metrics (highly connected "hub" clusters vs. singleton features) to check whether the previously-tested individual features were mostly isolated singletons (poor test power) vs. members of large families now aggregated away.
6. If any cluster survives, pull `aligned_isf.tsv` to check whether the surviving edges are true structural analogs or artifactual in-source-fragment relationships before over-interpreting.

## Expected Insights
Either (a) a compound family — plausibly the carotenoid/terpenoid neighborhood — shows a phylogenetically-controlled association with a*/C* that no individual feature reached significance on alone, which would be the first positive color-metabolome link in the project and directly testable against the pigment-biosynthesis prior; or (b) the network-aggregated scan is also null, which meaningfully strengthens the existing null conclusion by ruling out "buried in feature-level noise" as an explanation and shifting the burden of proof toward "color is not compositionally MS-detectable under current media/extraction conditions" or "genetic, not metabolic-flux, in origin."

## Feasibility
- Effort: Medium
- Data ready: Mostly (network files exist; need to check `filtered_pairs.tsv` vs `merged_pairs.tsv` edge-weight semantics and confirm graphml node IDs map cleanly onto the 10,949 deduplicated compound groups used in prior tests)
- Methods available: Standard tools (networkx/igraph community detection, existing R PGLS/permutation scripts from `analysis/copper/scripts/02_pgls_analysis.R` and `analysis/YPD/color_shape_growth/scripts/02_pgls_analysis.R` as templates)
- Key risk: Community detection resolution parameter is a hidden researcher degree of freedom — must pre-register or sensitivity-sweep the clustering resolution rather than picking whichever gives a hit; also cluster-level aggregation can wash out a real signal if only one feature in a large family truly carries pigment information

---

# Dark-Matter MS2 Spectral-Similarity Kernel Regression Against Color

## Persona
**Metabolomics / Mass Spectrometry Specialist** — untargeted MS2 is 95%+ "dark matter" by design; the fix is never to wait for library IDs, it's to use the raw fragmentation spectra themselves as a chemical similarity space and regress the phenotype against the whole spectral landscape rather than against annotation-gated features.

## Motivation
Only 790/16,332 features (4.8%) carry any SIRIUS annotation, and SIRIUS is currently being re-run with more reference data — meaning every test so far implicitly treated the unannotated 95.2% as either invisible (if analysis subset to annotated features) or as undifferentiated background noise (if analysis lumped them in as anonymous "compound groups" with no chemical structure imposed between them). A metabolomics specialist's default reflex when annotation coverage is this poor is to stop waiting on identification and instead use MS2 spectral similarity itself (cosine/modified-cosine spectral matching, or embedding-based spectral fingerprints) as a continuous chemical-similarity kernel — this is exactly the reasoning behind GNPS spectral networking and newer tools like Spec2Vec/MS2DeepScore, and it lets every one of the 16,332 features contribute regardless of whether SIRIUS ever names it.

## Connection to Existing Data
- Targets the exact blind spot named in the brief: "Are we throwing away 95% of the informative signal by only looking at the 4.8% SIRIUS-annotated slice?"
- Uses the raw MS2 spectra that underlie the aligned feature table (upstream of the EB pipeline's 16,332-feature output) — data that must already exist as the pipeline's raw `.mzML`/consensus-spectra input, just never used as *spectra* in the phenotype analyses (only their intensities were).
- Builds on, rather than duplicates, the molecular-network idea above: that approach uses the pipeline's *already-computed* pairwise network edges (fast, cheap); this approach recomputes a full pairwise MS2 similarity kernel (or fits an ML embedding) directly from spectra, which is more expensive but not gated by whatever edge-filtering thresholds EB's network construction already applied, and can pick up sub-threshold structural similarity the graphml discarded.
- Directly usable with the same phylogenetic block-permutation and colony-area decoy-gate framework validated in tests 1-5, and can specifically interrogate the unexplained colony-area/cell-fraction confound by asking whether the metabolites driving that confound cluster spectrally (i.e., are they one biosynthetic family, e.g., membrane lipids or storage compounds scaling with cell mass, rather than 500 independent hits).
- Ambitious sibling to Idea 1: instead of clustering with EB's pre-set network thresholds, this fits a kernel/embedding model (e.g., kernel ridge regression or MMR/Mantel test: spectral-similarity distance matrix vs. color-difference matrix) that can also flexibly incorporate colony area as a covariate kernel.

## Approach
1. Extract per-feature representative MS2 spectra (peak lists) from the EB pipeline's raw spectral output (not just the aligned intensity table) for all 16,332 features, or the 10,949 deduplicated groups as the working unit.
2. Compute a pairwise spectral similarity matrix using modified cosine (accounting for precursor mass shifts, standard in GNPS-style workflows) and/or fit an unsupervised spectral embedding (e.g., MS2DeepScore or simple Spec2Vec-style peak-fragment embedding) to get each feature/compound-group a low-dimensional chemical-structure-informed coordinate, independent of SIRIUS annotation status.
3. Use this embedding two ways: (a) as unsupervised structure to re-cluster features at finer resolution than the EB network, cross-checked against Idea 1's clusters for agreement; (b) as a kernel for a Mantel/RQV-style test — sample-to-sample metabolome dissimilarity computed with spectral-informed weighting vs. sample-to-sample color (deltaE00 or a*/b* Euclidean) dissimilarity, with the phylogenetic block permutation as the null model exactly as in test 1.
4. Explicitly build a parallel kernel test using colony area as the response variable (positive control channel, given it's the trait known to correlate with cell-fraction metabolite abundance) to confirm the kernel approach recovers the known area signal before trusting a color null.
5. If a color signal emerges only in dark-matter (unannotated) regions of the embedding, flag those specific spectral neighborhoods for prioritized manual MS2 inspection / targeted re-analysis once the PI's SIRIUS re-run with expanded references completes.
6. Report embedding-neighborhood-level effect sizes, not single p-values, given multiple embedding/kernel choices are being compared (sensitivity sweep required per repo convention).

## Expected Insights
If dark-matter spectral structure correlates with color where annotated-feature-level tests found nothing, it identifies exactly which unidentified chemical family to prioritize for the incoming SIRIUS re-run and gives the PI a concrete "spend annotation effort here" target instead of annotating uniformly. If the kernel test also comes up null for color but strongly positive for the colony-area positive-control channel (as expected from the existing area-vs-fraction confound), it is strong, method-agnostic confirmation that the color null is real and not an artifact of annotation coverage or feature-level test insensitivity — closing off the "we just haven't looked hard enough at the chemistry" objection before this goes to writeup.

## Feasibility
- Effort: High
- Data ready: Needs preprocessing (raw per-feature MS2 spectra/peak lists must be pulled from EB's intermediate outputs, not just the final aligned intensity table; format and availability need to be confirmed with pipeline outputs)
- Methods available: Research-grade (modified-cosine matching is standard/matchms-implementable; MS2DeepScore/Spec2Vec embeddings are available as pretrained or trainable tools but applying them well to this instrument's spectra, and interpreting a Mantel-style kernel test with proper phylogenetic block permutation, is a non-trivial custom pipeline)
- Key risk: Without a positive control succeeding first (colony-area kernel test), a null color result here is uninterpretable — could mean "no signal" or "kernel/embedding poorly calibrated for this data"; also compute cost of full pairwise MS2 similarity at 10,949+ compound groups across 594 samples could be substantial and needs early scoping
