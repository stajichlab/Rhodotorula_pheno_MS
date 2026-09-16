#!/usr/bin/env Rscript
# N-acyl amino acid candidates -- phylogenetic signal test (2026-09-16, PI
# request). Same method and tool as analysis/scripts/phylogenetic_signal.R
# (Phase 1, color phenotype): Blomberg's K and Pagel's lambda via
# phytools::phylosig, on the species-level tree
# (analysis/integrated_analysis/phase1_phenotype/species_tree.nwk).
#
# This answers "is compound production phylogenetically structured" --
# a DIFFERENT question from phase2_color_metabolome_association.py's
# block-permutation design, which tests whether compound abundance
# correlates with an EXTERNAL phenotype (color, growth) while controlling
# for phylogeny as a nuisance. That test has not been run here; see
# NACYL_AMINO_ACID_SEARCH.md for why K/lambda was judged the more direct
# tool for the literal question asked, and for a note on what the phase2
# association design would additionally test if wanted.
#
# Response: mean log10(peak_area+1) per species, CELL FRACTION (see
# nacyl_amino_acid_build_species_table.py -- supernatant is not tested,
# these compounds are essentially absent there per the compartment
# analysis). Both phylosig(K) and phylosig(lambda) run their own built-in
# permutation/likelihood-ratio null test (K: nsim tip-label permutations;
# lambda: LR test against lambda=0) -- no separate decoy gate is added
# here, unlike phase2_color_metabolome_association.py's external decoy
# requirement, because K/lambda's significance test is intrinsic to the
# method itself, not a separate design choice.
#
# CAVEAT (stated up front, not just in a footnote): several of these
# species have n_strains=1 (Cystobasidium sp., Pseudomicrostroma
# phylloplanum, R. araucariae, R. sp. clade XIII) -- their species-level
# "mean" is a single strain's value, not an average. Four compounds also
# show hard zeros (non-detection) in several species -- K/lambda assume a
# continuous Brownian-motion-like trait; a trait with structural zero-
# inflation (presence/absence mixed with a continuous magnitude among
# detected strains) is not a clean fit to that assumption. Results here
# should be read as an approximation, flagged explicitly, not a
# textbook-clean phylogenetic signal estimate.

suppressMessages({
  library(ape)
  library(phytools)
  library(optparse)
})

option_list <- list(
  make_option("--tree", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_tree.nwk"),
  make_option("--species-table", type = "character",
              default = "analysis/ahl_autoinducer_search/nacyl_amino_acid_species_table.csv"),
  make_option("--out", type = "character",
              default = "analysis/ahl_autoinducer_search/nacyl_amino_acid_phylogenetic_signal.csv"),
  make_option("--nsim", type = "integer", default = 1000)
)
opt <- parse_args(OptionParser(option_list = option_list))

tree <- read.tree(opt$tree)
dat <- read.csv(opt[["species-table"]], stringsAsFactors = FALSE)
dat$species <- gsub(" ", "_", dat$species)
rownames(dat) <- dat$species

missing_tips <- setdiff(tree$tip.label, dat$species)
if (length(missing_tips) > 0) {
  stop("Species tree has tips absent from --species-table: ", paste(missing_tips, collapse = ", "))
}
dat <- dat[tree$tip.label, ]

zero_edges <- which(tree$edge.length == 0)
if (length(zero_edges) > 0) {
  cat(sprintf("Zero-length edges in species tree: %d (nudged to 1e-6)\n", length(zero_edges)))
  tree$edge.length[zero_edges] <- 1e-6
}

traits <- c("row_4109_palmitoyl_arginine", "row_51126_myristoyl_arginine",
            "row_51152_palmitoyl_lysine", "row_24998_histidine_MH", "row_40740_histidine_MNa")
traits <- intersect(traits, colnames(dat))

results <- list()
for (tr in traits) {
  x <- setNames(dat[[tr]], dat$species)
  x <- x[!is.na(x)]
  n_zero <- sum(x == 0)
  sub_tree <- if (length(x) < Ntip(tree)) drop.tip(tree, setdiff(tree$tip.label, names(x))) else tree

  k_res <- tryCatch(
    phylosig(sub_tree, x, method = "K", test = TRUE, nsim = opt$nsim),
    error = function(e) NULL
  )
  lambda_res <- tryCatch(
    phylosig(sub_tree, x, method = "lambda", test = TRUE),
    error = function(e) NULL
  )

  results[[tr]] <- data.frame(
    compound = tr,
    n_species = length(x),
    n_species_zero = n_zero,
    K = if (!is.null(k_res)) k_res$K else NA,
    K_p_value = if (!is.null(k_res)) k_res$P else NA,
    lambda = if (!is.null(lambda_res)) lambda_res$lambda else NA,
    lambda_logL = if (!is.null(lambda_res)) lambda_res$logL else NA,
    lambda_p_value = if (!is.null(lambda_res)) lambda_res$P else NA
  )
}

res <- do.call(rbind, results)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(res, opt$out, row.names = FALSE)

cat(sprintf("\nWrote %s\n", opt$out))
cat("Phylogenetic signal (species-level tree, n=", Ntip(tree), " species, CELL FRACTION):\n", sep = "")
print(res, row.names = FALSE)
cat(
  "\nCaveat: several species contribute only 1 strain, and several\n",
  "compounds have hard zeros (non-detection) in multiple species --\n",
  "K/lambda assume a Brownian-motion-like continuous trait; a mixed\n",
  "presence/absence + magnitude trait is an approximation to that, not a\n",
  "clean fit. Read results as approximate, not textbook-standard.\n",
  sep = ""
)
