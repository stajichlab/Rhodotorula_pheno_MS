#!/usr/bin/env Rscript
# Phylogenetically-corrected re-test of the naive AA-vs-YPD-phenotype
# correlations from 01_naive_correlation.py. Mirrors
# analysis/copper/scripts/02_pgls_analysis.R; see that script/doc for the
# corPagel / zero-branch-nudge rationale, identical here.

suppressMessages({
  library(ape)
  library(nlme)
})

repo <- "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS"
treefile <- file.path(repo, "BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile")
master_csv <- file.path(repo, "analysis/YPD/color_shape_growth/outputs/ypd_aa_master_table.csv")
out_csv <- file.path(repo, "analysis/YPD/color_shape_growth/outputs/pgls_correlation_results.csv")
pruned_tree_out <- file.path(repo, "analysis/YPD/color_shape_growth/outputs/pruned_tree_146strains.nwk")

dat <- read.csv(master_csv, stringsAsFactors = FALSE)
stopifnot(!any(duplicated(dat$tree_tip)))
rownames(dat) <- dat$tree_tip

tree <- read.tree(treefile)
drop_tips <- setdiff(tree$tip.label, dat$tree_tip)
pruned <- drop.tip(tree, drop_tips)
stopifnot(Ntip(pruned) == nrow(dat))

zero_edges <- which(pruned$edge.length == 0)
cat(sprintf("Zero-length edges in pruned tree: %d (nudged to 1e-6)\n", length(zero_edges)))
pruned$edge.length[zero_edges] <- 1e-6
write.tree(pruned, pruned_tree_out)

dat <- dat[pruned$tip.label, ]

aa20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
phenotypes <- c("Mean_Shape_Area", "Mean_ColorLab_L.Mean", "Mean_ColorLab_a.Mean", "Mean_ColorLab_b.Mean")
# read.csv mangles the literal "L*Mean" column name to "L.Mean" etc.; recover actual names
pheno_cols <- grep("^Mean_(Shape_Area|ColorLab_)", names(dat), value = TRUE)
cat("Phenotype columns found:", paste(pheno_cols, collapse = ", "), "\n")

all_results <- list()
for (pheno in pheno_cols) {
  for (aa in aa20) {
    d <- data.frame(y = dat[[pheno]], x = dat[[aa]], tip = dat$tree_tip)
    rownames(d) <- d$tip
    fit <- tryCatch(
      gls(y ~ x, data = d, correlation = corPagel(0.5, phy = pruned, form = ~tip, fixed = FALSE), method = "ML"),
      error = function(e) NULL
    )
    key <- paste(pheno, aa)
    if (is.null(fit)) {
      all_results[[key]] <- data.frame(phenotype = pheno, amino_acid = aa, n = nrow(d),
                                        slope = NA, p_value = NA, lambda = NA, converged = FALSE)
      next
    }
    s <- summary(fit)
    coef_row <- s$tTable["x", ]
    lambda_est <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
    all_results[[key]] <- data.frame(phenotype = pheno, amino_acid = aa, n = nrow(d),
                                      slope = coef_row["Value"], p_value = coef_row["p-value"],
                                      lambda = lambda_est, converged = TRUE)
  }
}
res <- do.call(rbind, all_results)
res$p_value_fdr_bh <- NA
for (pheno in pheno_cols) {
  idx <- res$phenotype == pheno
  res$p_value_fdr_bh[idx] <- p.adjust(res$p_value[idx], method = "BH")
}
res <- res[order(res$phenotype, res$p_value), ]
write.csv(res, out_csv, row.names = FALSE)

cat(sprintf("\nWrote %s\n", out_csv))
for (pheno in pheno_cols) {
  hits <- res[res$phenotype == pheno & !is.na(res$p_value_fdr_bh) & res$p_value_fdr_bh < 0.05, "amino_acid"]
  cat(sprintf("  %s: %s\n", pheno, paste(hits, collapse = ", ")))
}
