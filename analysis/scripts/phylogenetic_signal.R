#!/usr/bin/env Rscript
# Phase 1 (analysis/INTEGRATED_ANALYSIS_STRATEGY.md), step 2: test for
# phylogenetic signal in the color phenotype across the species tree, using
# Blomberg's K and Pagel's lambda (phytools::phylosig).
#
# Run at the SPECIES level (species_tree.nwk from prune_species_tree.R +
# species_phenotype_table.csv from build_species_level_tables.py), per the
# "Species-Level Collapse" section of the strategy doc -- comparative-method
# signal statistics assume tips are independent evolutionary lineages, so
# feeding them multiple non-independent strains per species as separate
# tips would misstate both K and lambda.
#
# K and lambda both answer "how much does phenotypic similarity track
# phylogenetic relatedness" but are not interchangeable: K compares the
# trait's variance to a strict Brownian-motion expectation (K=1), while
# lambda scales the phylogenetic correlation structure itself and is
# estimated by ML, so lambda is more robust to a trait that evolved under
# something other than pure Brownian motion. Both are reported; a large
# K/lambda alongside a clade that is *also* the outlier for the trait
# (R. dairenensis here) does not by itself distinguish one deep trait-origin
# event from many independent ones -- see the "n=1 clade" discussion in the
# strategy doc. This script establishes the magnitude of phylogenetic signal;
# it is not, on its own, the convergent-evolution test (that's
# convergent_color_test.py).
#
# h_deg (hue angle) is circular data; phylosig assumes a linear trait, so
# its K/lambda for h_deg should be read with that caveat -- included for
# completeness, not as a primary result.

suppressMessages({
  library(ape)
  library(phytools)
  library(optparse)
})

option_list <- list(
  make_option("--tree", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_tree.nwk"),
  make_option("--species-table", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_phenotype_table.csv"),
  make_option("--out", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/phenotype_phylogenetic_signal.csv"),
  make_option("--nsim", type = "integer", default = 1000,
              help = "Number of randomizations for Blomberg's K permutation p-value.")
)
opt <- parse_args(OptionParser(option_list = option_list))

tree <- read.tree(opt$tree)
dat <- read.csv(opt[["species-table"]], stringsAsFactors = FALSE)
# ape::write.tree/read.tree round-trips unquoted spaces in tip labels as
# underscores (standard Newick behavior) -- species names came out of
# prune_species_tree.R with spaces (e.g. "Rhodotorula dairenensis") and
# come back from disk as "Rhodotorula_dairenensis". Normalize the table's
# species key the same way rather than trying to re-quote tip labels.
dat$species <- gsub(" ", "_", dat$species)
rownames(dat) <- dat$species

missing_tips <- setdiff(tree$tip.label, dat$species)
if (length(missing_tips) > 0) {
  stop("Species tree has tips absent from --species-table: ", paste(missing_tips, collapse = ", "))
}
extra_rows <- setdiff(dat$species, tree$tip.label)
if (length(extra_rows) > 0) {
  message(sprintf(
    "%d species in --species-table have no tree tip (no genome data) and are excluded from this test: %s",
    length(extra_rows), paste(extra_rows, collapse = ", ")
  ))
}
dat <- dat[tree$tip.label, ]

# zero-length branches make some tree operations behave oddly for signal
# stats too (e.g. K's expected-variance ratio); nudge as in
# copper/scripts/02_pgls_analysis.R
zero_edges <- which(tree$edge.length == 0)
if (length(zero_edges) > 0) {
  cat(sprintf("Zero-length edges in species tree: %d (nudged to 1e-6)\n", length(zero_edges)))
  tree$edge.length[zero_edges] <- 1e-6
}

# read.csv() (check.names=TRUE default) mangles "L*_mean" -> "L._mean" etc.,
# since "*" isn't a valid R name character -- match against the mangled form.
traits <- c("L._mean", "a._mean", "b._mean", "C._mean", "h_deg_mean", "orange_score_mean")
traits <- intersect(traits, colnames(dat))
trait_display <- c(L._mean = "L*", a._mean = "a*", b._mean = "b*", C._mean = "C*",
                    h_deg_mean = "h_deg", orange_score_mean = "orange_score")

results <- list()
for (tr in traits) {
  x <- setNames(dat[[tr]], dat$species)
  x <- x[!is.na(x)]
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
    trait = trait_display[[tr]],
    n_species = length(x),
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
cat("Phylogenetic signal (species-level tree, n=", Ntip(tree), " species):\n", sep = "")
print(res, row.names = FALSE)
cat(
  "\nInterpretation note: high K/lambda for orange_score alongside R. dairenensis\n",
  "being the trait outlier does NOT by itself establish independent convergent\n",
  "origins -- see convergent_color_test.py and the 'n=1 clade' section of\n",
  "analysis/INTEGRATED_ANALYSIS_STRATEGY.md.\n",
  sep = ""
)
