#!/usr/bin/env Rscript
# Idea 5, Step 1 (analysis/ideas/2026-08-15-color-metabolome-genome-null-brainstorm/
# DEVELOPMENT_PLAN.md, Part C): formal ancestral-state / regime-shift
# detection on the species-level color trait, replacing the coarse
# heuristic in convergent_color_test.R (above-mean-orange_score AND
# phylogenetically-distant-from-R.-dairenensis).
#
# Method: bayou (Uyeda & Harmon; reversible-jump MCMC over an
# Ornstein-Uhlenbeck process with an unknown number/location of
# adaptive-regime shifts) on the species-level tree
# (species_tree.nwk, 17 tips -- built by prune_species_tree.R). This
# directly answers "which branches most likely carry an independent
# color-gain event" with posterior support, rather than a hand-tuned
# above/below-mean threshold.
#
# Trait: a*_mean (PRIMARY predictor per the 2026-08-15 grilling session
# decision -- see .living/decisions.md). orange_score_mean can be rerun
# the same way as a documented follow-up (exploratory tier, per the same
# decision) by changing TRAIT_COL below.
#
# Tree size caveat: 17 tips is small for bayou (most published examples
# use 50-300+ tips). This is run as an exploratory/preliminary pass --
# generations kept modest (fast on this tree size) and results should be
# read as "candidate shift locations for Step 2's contrast-pair
# definition," not as publication-grade posterior support without a
# convergence-diagnostic pass (multiple independent chains, effective
# sample size check) which this script does not yet include.

suppressMessages({
  library(ape)
  library(bayou)
  library(optparse)
})

option_list <- list(
  make_option("--tree", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_tree.nwk"),
  make_option("--species-table", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_phenotype_table.csv"),
  make_option("--trait-col", type = "character", default = "a._mean",
              help = "Column in --species-table to use as the trait (R read.csv mangles 'a*_mean' -> 'a._mean'). Default: a._mean (primary predictor)."),
  make_option("--ngen", type = "integer", default = 20000,
              help = "MCMC generations. Small tree (17 tips) -> fast; increase for a publication-grade run."),
  make_option("--kmax", type = "integer", default = 5,
              help = "Max number of regime shifts allowed by the prior (17 tips / ~32 edges -> keep small)."),
  make_option("--seed", type = "integer", default = 0),
  make_option("--out-dir", type = "character",
              default = "analysis/integrated_analysis/phase5_genome_linkage/idea5_regime_shift")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

dir.create(opt[["out-dir"]], recursive = TRUE, showWarnings = FALSE)
run_name <- paste0("regime_shift_", gsub("[^A-Za-z0-9]", "", opt[["trait-col"]]))

tree <- read.tree(opt$tree)
dat_table <- read.csv(opt[["species-table"]], stringsAsFactors = FALSE)
# ape's Newick round-trip turns spaces in species names into underscores (see
# phylogenetic_signal.R / convergent_color_test.R for the same normalization)
dat_table$species <- gsub(" ", "_", dat_table$species)

if (!(opt[["trait-col"]] %in% colnames(dat_table))) {
  stop("--trait-col '", opt[["trait-col"]], "' not found. Available columns: ",
       paste(colnames(dat_table), collapse = ", "))
}
missing_sp <- setdiff(tree$tip.label, dat_table$species)
if (length(missing_sp) > 0) stop("Species tree has tips with no phenotype row: ", paste(missing_sp, collapse = ", "))

dat_table <- dat_table[match(tree$tip.label, dat_table$species), ]
dat <- setNames(dat_table[[opt[["trait-col"]]]], dat_table$species)
if (any(is.na(dat))) stop("NA trait values for: ", paste(names(dat)[is.na(dat)], collapse = ", "))

cat(sprintf("Trait: %s | n species: %d | range: %.3f to %.3f\n", opt[["trait-col"]], length(dat), min(dat), max(dat)))

# zero-length branches break bayou's likelihood the same way they break
# corPagel (see copper/scripts/02_pgls_analysis.R) -- nudge them
zero_edges <- which(tree$edge.length == 0)
if (length(zero_edges) > 0) {
  cat(sprintf("Zero-length edges: %d (nudged to 1e-6)\n", length(zero_edges)))
  tree$edge.length[zero_edges] <- 1e-6
}

prior <- make.prior(
  tree,
  dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", dk = "cdpois", dtheta = "dnorm"),
  param = list(
    dalpha = list(scale = 0.1),
    dsig2 = list(scale = 0.1),
    dk = list(lambda = 1, kmax = opt$kmax),
    dtheta = list(mean = mean(dat), sd = 1.5 * sd(dat))
  ),
  plot.prior = FALSE  # this is a non-interactive/batch run -- no X11 display available
)

setwd_orig <- getwd()
mcmc_dir <- file.path(normalizePath(opt[["out-dir"]]), "mcmc_chain")
dir.create(mcmc_dir, recursive = TRUE, showWarnings = FALSE)

mcmc <- bayou.makeMCMC(tree, dat, prior = prior, file.dir = mcmc_dir, outname = run_name, plot.freq = NULL, ticker.freq = max(opt$ngen %/% 10, 1000))
mcmc$run(opt$ngen)

chain <- mcmc$load()
chain <- set.burnin(chain, 0.3)

sink(file.path(opt[["out-dir"]], paste0(run_name, "_summary.txt")))
cat("bayou regime-shift detection summary\n")
cat("=====================================\n")
cat(sprintf("Trait: %s | tree: %s | n species: %d\n", opt[["trait-col"]], opt$tree, length(dat)))
cat(sprintf("MCMC generations: %d | burnin: 30%%\n\n", opt$ngen))
print(summary(chain))
sink()

shiftsumm <- tryCatch(
  shiftSummaries(chain, mcmc, pp.cutoff = 0.3),
  error = function(e) {
    cat("shiftSummaries failed (likely no shifts cleared the pp.cutoff in this exploratory run):", conditionMessage(e), "\n")
    NULL
  }
)

# per-branch posterior probability of carrying a shift -- the actual
# "which branches most likely carry an independent color-gain event"
# answer Step 2 needs, regardless of whether shiftSummaries succeeded
branch_pp <- Lposterior(chain, tree)
branch_pp$species_tip_label <- NA
for (i in seq_len(nrow(branch_pp))) {
  edge_idx <- as.integer(rownames(branch_pp)[i])
  if (!is.na(edge_idx) && edge_idx <= nrow(tree$edge)) {
    child <- tree$edge[edge_idx, 2]
    if (child <= length(tree$tip.label)) {
      branch_pp$species_tip_label[i] <- tree$tip.label[child]
    }
  }
}
branch_pp <- branch_pp[order(-branch_pp$pp), ]
write.csv(branch_pp, file.path(opt[["out-dir"]], paste0(run_name, "_branch_posterior.csv")), row.names = TRUE)

cat("\nTop branches by posterior shift probability (terminal-branch rows show the species tip; internal rows are ancestral branches, no single species label):\n")
print(head(branch_pp[, c("pp", "species_tip_label")], 10))

cat(sprintf("\nWrote %s\n", file.path(opt[["out-dir"]], paste0(run_name, "_summary.txt"))))
cat(sprintf("Wrote %s\n", file.path(opt[["out-dir"]], paste0(run_name, "_branch_posterior.csv"))))
