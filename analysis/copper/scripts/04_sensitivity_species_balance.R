#!/usr/bin/env Rscript
# Sensitivity check: 99/134 strains (74%) are R. mucilaginosa. Pagel's-lambda
# PGLS already down-weights shared ancestry, but a single dominant species
# can still drive an apparent trend if that species alone spans most of the
# AUC and AA-frequency range. Refit PGLS restricted to R. mucilaginosa only
# (tests whether the signal holds *within* one species, i.e. is not solely
# an inter-species effect) to see whether the top hits from 02 are
# consistent in sign/magnitude in that subset.

suppressMessages({
  library(ape)
  library(nlme)
})

repo <- "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS"
master_csv <- file.path(repo, "analysis/copper/outputs/copper_aa_master_table.csv")
pruned_tree_file <- file.path(repo, "analysis/copper/outputs/pruned_tree_134strains.nwk")
out_csv <- file.path(repo, "analysis/copper/outputs/sensitivity_mucilaginosa_only.csv")

dat <- read.csv(master_csv, stringsAsFactors = FALSE)
rownames(dat) <- dat$tree_tip
tree <- read.tree(pruned_tree_file)

muc <- dat[dat$species == "Rhodotorula mucilaginosa", ]
sub_tree <- drop.tip(tree, setdiff(tree$tip.label, muc$tree_tip))
muc <- muc[sub_tree$tip.label, ]
cat(sprintf("R. mucilaginosa subset: n=%d strains\n", nrow(muc)))

# top hits from the full-dataset comparison (03): S, L, Q, W, T, C
target_aa <- c("S", "L", "Q", "W", "T", "C")

results <- list()
for (aa in target_aa) {
  d <- data.frame(auc = muc$mean_auc_rate, x = muc[[aa]], tip = muc$tree_tip)
  rownames(d) <- d$tip
  fit <- tryCatch(
    gls(auc ~ x, data = d, correlation = corPagel(0.5, phy = sub_tree, form = ~tip, fixed = FALSE), method = "ML"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    results[[aa]] <- data.frame(amino_acid = aa, n = nrow(d), slope = NA, p_value = NA, lambda = NA, converged = FALSE)
    next
  }
  s <- summary(fit)
  coef_row <- s$tTable["x", ]
  lambda_est <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  results[[aa]] <- data.frame(amino_acid = aa, n = nrow(d), slope = coef_row["Value"],
                               p_value = coef_row["p-value"], lambda = lambda_est, converged = TRUE)
}
res <- do.call(rbind, results)
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
print(res, row.names = FALSE)
