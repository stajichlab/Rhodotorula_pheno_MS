#!/usr/bin/env Rscript
# Idea 5, Step 2 (DEVELOPMENT_PLAN.md Part C): define independent
# phylogenetic contrast pairs from Step 1's bayou branch-posterior output
# (idea5_regime_shift_detection.R). This is the correct evolutionary unit
# of replication for the convergence test -- N=independent origins, not
# N=strains -- feeding both Idea 5's own metabolome-side follow-up and,
# once ready, Idea 3's genome-side candidate-gene calls.
#
# Adjustment from the original plan: Step 1's posterior support is
# strongest on the terminal branches to R. taiwanensis / R. sphaerocarpa
# (pp ~0.26 / 0.12 on the post-ANI 16-tip tree — see
# idea5_regime_shift/regime_shift_amean_summary.txt), still below a clean
# shift/non-shift partition. Rather than force a hard threshold that
# isn't supported by the data, this script takes the
# top-N ranked branches by posterior probability as CANDIDATE shift
# clades (explicitly labeled as such, not confirmed shifts) and pairs
# each against its phylogenetically nearest non-candidate clade.

suppressMessages({
  library(ape)
  library(optparse)
})

option_list <- list(
  make_option("--tree", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_tree.nwk"),
  make_option("--species-table", type = "character",
              default = "analysis/integrated_analysis/phase1_phenotype/species_phenotype_table.csv"),
  make_option("--branch-posterior", type = "character",
              default = "analysis/integrated_analysis/phase5_genome_linkage/idea5_regime_shift/regime_shift_amean_branch_posterior.csv"),
  make_option("--top-n", type = "integer", default = 5,
              help = "Number of top-pp branches to treat as candidate shift clades."),
  make_option("--out", type = "character",
              default = "analysis/integrated_analysis/phase5_genome_linkage/idea5_regime_shift/contrast_pairs.csv")
)
opt <- parse_args(OptionParser(option_list = option_list))

tree <- read.tree(opt$tree)
dat <- read.csv(opt[["species-table"]], stringsAsFactors = FALSE)
dat$species <- gsub(" ", "_", dat$species)
rownames(dat) <- dat$species

bp <- read.csv(opt[["branch-posterior"]], row.names = 1)
# branch_posterior.csv row names are edge indices (see idea5_regime_shift_detection.R)
bp$edge_index <- as.integer(rownames(bp))
bp <- bp[order(-bp$pp), ]
top <- head(bp, opt[["top-n"]])

# lightweight descendant-tip lookup without requiring the phangorn package
node_descendant_tips <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node <= n_tip) return(tree$tip.label[node])
  kids <- tree$edge[tree$edge[, 1] == node, 2]
  unlist(lapply(kids, function(k) node_descendant_tips(tree, k)))
}

pdist <- cophenetic.phylo(tree)

results <- list()
candidate_tip_sets <- list()
for (i in seq_len(nrow(top))) {
  edge_idx <- top$edge_index[i]
  child_node <- tree$edge[edge_idx, 2]
  tips_in_clade <- node_descendant_tips(tree, child_node)
  candidate_tip_sets[[as.character(edge_idx)]] <- tips_in_clade
}

all_candidate_tips <- unique(unlist(candidate_tip_sets))
noncandidate_tips <- setdiff(tree$tip.label, all_candidate_tips)

for (i in seq_len(nrow(top))) {
  edge_idx <- as.character(top$edge_index[i])
  clade_tips <- candidate_tip_sets[[edge_idx]]
  # nearest non-candidate tip (or nearest tip not in this specific clade if
  # all tips are somehow candidates) by patristic distance, averaged if the
  # clade has >1 tip
  if (length(noncandidate_tips) > 0) {
    dists_to_noncand <- sapply(noncandidate_tips, function(nc) mean(pdist[clade_tips, nc]))
    nearest <- names(sort(dists_to_noncand))[1]
  } else {
    other_tips <- setdiff(tree$tip.label, clade_tips)
    dists_to_other <- sapply(other_tips, function(o) mean(pdist[clade_tips, o]))
    nearest <- names(sort(dists_to_other))[1]
  }
  results[[edge_idx]] <- data.frame(
    edge_index = top$edge_index[i],
    posterior_pp = top$pp[i],
    candidate_clade_tips = paste(clade_tips, collapse = ";"),
    candidate_clade_n_species = length(clade_tips),
    candidate_clade_mean_a = mean(dat[clade_tips, "a._mean"], na.rm = TRUE),
    nearest_noncandidate_sister = nearest,
    sister_a = dat[nearest, "a._mean"],
    patristic_distance = mean(pdist[clade_tips, nearest])
  )
}

out <- do.call(rbind, results)
out <- out[order(-out$posterior_pp), ]
out$species_underscore_note <- "species names use underscores (Newick round-trip); convert back to spaces for display"

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(out, opt$out, row.names = FALSE)

cat(sprintf("Top %d candidate shift clades (by posterior branch probability) paired with nearest non-candidate sister:\n\n", opt[["top-n"]]))
print(out[, c("posterior_pp", "candidate_clade_tips", "candidate_clade_mean_a", "nearest_noncandidate_sister", "sister_a", "patristic_distance")], row.names = FALSE)

n_distinct_clades <- length(unique(out$candidate_clade_tips))
cat(sprintf(
  "\n%d distinct candidate clades identified (some top-N branches may be nested/overlapping -- see candidate_clade_tips for overlap).\n",
  n_distinct_clades
))
cat(sprintf("Wrote %s\n", opt$out))
cat(
  "\nCAVEAT (per DEVELOPMENT_PLAN.md Part C): this is hypothesis-generation,\n",
  "not a confirmatory test. Posterior support for individual branches was\n",
  "strongest on the taiwanensis/sphaerocarpa subclade but stayed below ~0.3\n",
  "-- treat this contrast-pair list as candidates for Step 3 (shared\n",
  "molecular correlates), not as confirmed independent origins.\n",
  sep = ""
)
