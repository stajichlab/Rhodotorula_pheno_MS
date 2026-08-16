#!/usr/bin/env Rscript
# Phylogenetically-corrected re-test of the naive AA-vs-copper-AUC
# correlations from 01_naive_correlation.py, using the phyling fungi_odb10
# protein tree pruned to the 134 strains with both genome and Cu_AUC data.
#
# Model: mean_auc_rate ~ AA_frequency, fit by GLS with a Pagel's-lambda
# correlation structure (ape::corPagel) on the pruned tree. Lambda is
# estimated by ML, not fixed at 1, so the model interpolates between the
# naive (lambda=0, star phylogeny / OLS-equivalent) and full Brownian-motion
# (lambda=1) cases rather than assuming either extreme.
#
# Zero-length branches (20 in the source tree, from unresolved/near-identical
# strain pairs) make the phylogenetic covariance matrix singular under
# corPagel; each is nudged to a small positive length (1e-6) before fitting,
# a standard workaround for near-zero polytomies in comparative methods.

suppressMessages({
  library(ape)
  library(nlme)
})

repo <- "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS"
treefile <- file.path(repo, "BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile")
master_csv <- file.path(repo, "analysis/copper/outputs/copper_aa_master_table.csv")
out_csv <- file.path(repo, "analysis/copper/outputs/pgls_correlation_results.csv")
pruned_tree_out <- file.path(repo, "analysis/copper/outputs/pruned_tree_134strains.nwk")

dat <- read.csv(master_csv, stringsAsFactors = FALSE)
stopifnot(!any(duplicated(dat$tree_tip)))
rownames(dat) <- dat$tree_tip

tree <- read.tree(treefile)
missing_tips <- setdiff(dat$tree_tip, tree$tip.label)
if (length(missing_tips) > 0) {
  stop("Tree tips missing for: ", paste(missing_tips, collapse = ", "))
}
drop_tips <- setdiff(tree$tip.label, dat$tree_tip)
pruned <- drop.tip(tree, drop_tips)
stopifnot(Ntip(pruned) == nrow(dat))

# nudge zero-length branches so the corPagel covariance matrix is non-singular
zero_edges <- which(pruned$edge.length == 0)
cat(sprintf("Zero-length edges in pruned tree: %d (nudged to 1e-6)\n", length(zero_edges)))
pruned$edge.length[zero_edges] <- 1e-6
write.tree(pruned, pruned_tree_out)

dat <- dat[pruned$tip.label, ]

aa20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
primary_aa <- c("C", "H", "D", "E", "M")

results <- list()
for (aa in aa20) {
  d <- data.frame(auc = dat$mean_auc_rate, x = dat[[aa]], tip = dat$tree_tip)
  rownames(d) <- d$tip

  fit_pagel <- tryCatch(
    gls(auc ~ x, data = d, correlation = corPagel(0.5, phy = pruned, form = ~tip, fixed = FALSE),
        method = "ML"),
    error = function(e) NULL
  )
  if (is.null(fit_pagel)) {
    results[[aa]] <- data.frame(amino_acid = aa, n = nrow(d), slope = NA, se = NA, p_value = NA,
                                 lambda = NA, aic = NA, converged = FALSE, primary_hypothesis = aa %in% primary_aa)
    next
  }
  s <- summary(fit_pagel)
  coef_row <- s$tTable["x", ]
  lambda_est <- as.numeric(coef(fit_pagel$modelStruct$corStruct, unconstrained = FALSE))

  results[[aa]] <- data.frame(
    amino_acid = aa,
    n = nrow(d),
    slope = coef_row["Value"],
    se = coef_row["Std.Error"],
    p_value = coef_row["p-value"],
    lambda = lambda_est,
    aic = AIC(fit_pagel),
    converged = TRUE,
    primary_hypothesis = aa %in% primary_aa
  )
}

res <- do.call(rbind, results)
res$p_value_fdr_bh <- p.adjust(res$p_value, method = "BH")
res <- res[order(res$p_value), ]
write.csv(res, out_csv, row.names = FALSE)

cat(sprintf("\nWrote %s\n", out_csv))
cat("PGLS results (Pagel's lambda, ML), sorted by raw p-value:\n")
print(res[, c("amino_acid", "slope", "p_value", "p_value_fdr_bh", "lambda", "primary_hypothesis")], row.names = FALSE)
