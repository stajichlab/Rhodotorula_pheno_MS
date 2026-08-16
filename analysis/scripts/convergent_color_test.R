#!/usr/bin/env Rscript
# Phase 1 gating step (see "Statistical Review & Revisions" ->
# "Core finding: the n=1 clade problem" in
# analysis/INTEGRATED_ANALYSIS_STRATEGY.md): does the darker-orange
# phenotype appear in more than one independently-evolved lineage?
#
# R. dairenensis being a single monophyletic clade means every genomic
# feature that distinguishes it from the rest of the tree is confounded
# with phylogeny -- PGLS/PIC/block permutation correct for expected
# covariance under a tree model, they do not manufacture a second,
# independent instance of "orange evolved here." If another species, ideally
# one NOT closely related to R. dairenensis, independently shows an
# elevated orange_score, that is the one thing that actually buys a second
# data point for phase 3-5 causal claims. This script does not run that
# downstream analysis; it identifies which species (if any) are worth
# treating as candidate convergent lineages before investing further.
#
# Method (exploratory flagging, not a formal test -- see caveat printed at
# the end): for each species (species-level tree tip, from
# prune_species_tree.R + build_species_level_tables.py output), report its
# orange_score_mean rank and its patristic distance from R. dairenensis on
# the species tree. A species is flagged "candidate_convergent" if BOTH:
#   (a) orange_score_mean is above the cross-species mean (i.e. on the
#       "darker orange" side at all), and
#   (b) its patristic distance from R. dairenensis is above the median
#       pairwise distance across all species pairs (i.e. it is NOT one of
#       dairenensis's closest relatives, where a similar score is more
#       likely to just be shared ancestry/incomplete lineage sorting rather
#       than independent convergence).
# This is a coarse first-pass filter meant to prioritize which species (if
# any) deserve a closer, dedicated convergence analysis -- not a
# publication-ready convergence test.

suppressMessages({
  library(ape)
  library(optparse)
})

option_list <- list(
  make_option("--tree", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_tree.nwk"),
  make_option("--species-table", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_phenotype_table.csv"),
  make_option("--focal-species", type = "character", default = "Rhodotorula dairenensis",
              help = "The species whose phenotype motivated this analysis."),
  make_option("--out", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/convergent_color_candidates.csv")
)
opt <- parse_args(OptionParser(option_list = option_list))

tree <- read.tree(opt$tree)
dat <- read.csv(opt[["species-table"]], stringsAsFactors = FALSE)
dat$species <- gsub(" ", "_", dat$species)  # match ape's Newick round-trip, see phylogenetic_signal.R
focal <- gsub(" ", "_", opt[["focal-species"]])

if (!(focal %in% tree$tip.label)) {
  stop("--focal-species '", opt[["focal-species"]], "' not found as a tree tip (as '", focal, "')")
}
dat <- dat[dat$species %in% tree$tip.label, ]
rownames(dat) <- dat$species

pdist <- cophenetic.phylo(tree)
dat$patristic_dist_from_focal <- pdist[dat$species, focal]

score_col <- "orange_score_mean"
if (!(score_col %in% colnames(dat))) {
  # read.csv mangles "orange_score_mean" fine (no special chars), this
  # branch only guards against a renamed input table
  stop("--species-table missing '", score_col, "' column")
}

mean_score <- mean(dat[[score_col]], na.rm = TRUE)
median_dist <- median(dat$patristic_dist_from_focal[dat$species != focal], na.rm = TRUE)

dat$above_mean_orange_score <- dat[[score_col]] > mean_score
dat$distant_from_focal <- dat$patristic_dist_from_focal > median_dist
dat$candidate_convergent <- dat$species != focal & dat$above_mean_orange_score & dat$distant_from_focal

dat <- dat[order(-dat[[score_col]]), ]
dat$orange_score_rank <- seq_len(nrow(dat))

out_cols <- c("species", "orange_score_rank", score_col, "patristic_dist_from_focal",
              "above_mean_orange_score", "distant_from_focal", "candidate_convergent")
out <- dat[, out_cols]
out$species <- gsub("_", " ", out$species)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(out, opt$out, row.names = FALSE)

n_candidates <- sum(out$candidate_convergent, na.rm = TRUE)
n_excluded <- sum(is.na(out$candidate_convergent))
cat(sprintf("\nWrote %s\n", opt$out))
cat(sprintf(
  "Focal species: %s (rank %d of %d by orange_score_mean = %.3f)\n",
  opt[["focal-species"]], out$orange_score_rank[out$species == opt[["focal-species"]]],
  nrow(out), out[[score_col]][out$species == opt[["focal-species"]]]
))
cat(sprintf("Cross-species mean orange_score_mean: %.3f | median patristic distance from focal: %.4f\n",
            mean_score, median_dist))
print(out, row.names = FALSE)
cat(sprintf(
  "\n%d candidate convergent species flagged (above-average orange_score AND phylogenetically distant from %s).\n",
  n_candidates, opt[["focal-species"]]
))
if (n_excluded > 0) {
  cat(sprintf("(%d species excluded from flagging: NA orange_score_mean, e.g. outside --genus-filter used upstream.)\n",
              n_excluded))
}
cat(
  "CAVEAT: this is a coarse exploratory filter (mean/median thresholds on\n",
  "n=", nrow(out), " species), not a formal convergence test. A flagged species\n",
  "is a candidate for closer inspection (ancestral state reconstruction, a\n",
  "dedicated OU/regime-shift model e.g. via l1ou or bayou, or literature check\n",
  "on that species' known pigmentation) before being treated as a genuine\n",
  "second independent origin for Phase 3-5 causal claims.\n",
  sep = ""
)
