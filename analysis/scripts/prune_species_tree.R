#!/usr/bin/env Rscript
# Prune a strain-level phylogeny (one tip per strain, e.g. the PHYling
# fungi_odb10 tree in BFD/results/phyling_pep/) down to one tip per species,
# for use with the species-level PGLS/PIC tests described in the
# "Species-Level Collapse: Procedure" section of
# analysis/INTEGRATED_ANALYSIS_STRATEGY.md.
#
# Method: for each species, pick one representative strain (by default the
# one with the highest value of --quality-col, e.g. BUSCO completeness or
# N50; falls back to the first tip alphabetically, with a warning, if no
# quality column is given) and drop every other strain of that species from
# the tree. This is the "preferred" option in the strategy doc rather than
# an explicit branch-length-averaging collapse, because ape::drop.tip's
# default collapse.singles=TRUE already does the right thing here: when all
# but one tip of a monophyletic clade are dropped, the internal nodes left
# with a single child are merged, so the remaining tip's edge length becomes
# the distance from the species' crown node's parent down through the
# representative's original tip -- i.e. the within-species branch length is
# absorbed away and the representative ends up positioned at the species'
# stem, not at its own arbitrary original depth. A separate "average the
# branch lengths across strains" implementation would fight what drop.tip
# already gives you for free.
#
# Requires the input strain<->species map to list every tree tip; any tip in
# the tree that's absent from the map is dropped from the pruned tree with a
# warning (e.g. genomes without phenotype/species data yet).

suppressMessages({
  library(ape)
  library(optparse)
})

option_list <- list(
  make_option("--tree", type = "character", help = "Input strain-level Newick tree."),
  make_option("--strain-species-map", type = "character",
              help = "CSV with columns tree_tip,species[,<quality-col>]."),
  make_option("--quality-col", type = "character", default = NULL,
              help = "Numeric column in the map used to pick the representative strain per species (higher = better, e.g. busco_completeness). If omitted, the first tip alphabetically is used and a warning is printed."),
  make_option("--out-tree", type = "character", help = "Output species-level Newick tree (tip labels = species names)."),
  make_option("--out-map", type = "character", help = "Output CSV: species, representative_tip, n_strains_in_species, within_species_max_patristic_dist, quality_value.")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (req in c("tree", "strain-species-map", "out-tree", "out-map")) {
  if (is.null(opt[[req]])) stop("Missing required argument: --", req)
}

tree <- read.tree(opt$tree)
map <- read.csv(opt[["strain-species-map"]], stringsAsFactors = FALSE)
stopifnot(all(c("tree_tip", "species") %in% colnames(map)))
stopifnot(!any(duplicated(map$tree_tip)))

quality_col <- opt[["quality-col"]]
if (!is.null(quality_col) && !(quality_col %in% colnames(map))) {
  stop("--quality-col '", quality_col, "' not found in ", opt[["strain-species-map"]])
}

tips_in_tree <- tree$tip.label
map_in_tree <- map[map$tree_tip %in% tips_in_tree, ]
tree_tips_not_in_map <- setdiff(tips_in_tree, map$tree_tip)
if (length(tree_tips_not_in_map) > 0) {
  warning(sprintf(
    "%d tree tip(s) absent from --strain-species-map; dropping from pruned tree: %s",
    length(tree_tips_not_in_map),
    paste(head(tree_tips_not_in_map, 10), collapse = ", ")
  ))
}
map_tips_not_in_tree <- setdiff(map$tree_tip, tips_in_tree)
if (length(map_tips_not_in_tree) > 0) {
  warning(sprintf(
    "%d strain(s) in --strain-species-map have no matching tree tip; ignored: %s",
    length(map_tips_not_in_tree),
    paste(head(map_tips_not_in_tree, 10), collapse = ", ")
  ))
}

# within-species patristic distance diagnostic, computed BEFORE pruning
# (needs the full strain-level tree; not meaningful after collapse)
full_dist <- cophenetic.phylo(tree)

pick_representative <- function(species_rows) {
  if (nrow(species_rows) == 1) {
    return(list(tip = species_rows$tree_tip[1], quality = NA_real_))
  }
  if (!is.null(quality_col)) {
    q <- suppressWarnings(as.numeric(species_rows[[quality_col]]))
    if (all(is.na(q))) {
      warning(sprintf(
        "species '%s': all %d strains have missing/non-numeric %s; falling back to alphabetical first tip",
        species_rows$species[1], nrow(species_rows), quality_col
      ))
      ord <- order(species_rows$tree_tip)
      return(list(tip = species_rows$tree_tip[ord[1]], quality = NA_real_))
    }
    best <- which.max(q)
    return(list(tip = species_rows$tree_tip[best], quality = q[best]))
  }
  ord <- order(species_rows$tree_tip)
  list(tip = species_rows$tree_tip[ord[1]], quality = NA_real_)
}

species_list <- split(map_in_tree, map_in_tree$species)
if (is.null(quality_col)) {
  warning("--quality-col not given: representative strain per species chosen alphabetically, ",
          "not by assembly quality. Pass --quality-col (e.g. busco_completeness) for a defensible choice.")
}

records <- list()
representative_tips <- character(0)
drop_tips <- character(0)

for (sp in names(species_list)) {
  rows <- species_list[[sp]]
  rep_info <- pick_representative(rows)
  rep_tip <- rep_info$tip
  other_tips <- setdiff(rows$tree_tip, rep_tip)

  within_max_dist <- if (length(rows$tree_tip) > 1) {
    max(full_dist[rows$tree_tip, rows$tree_tip])
  } else {
    NA_real_
  }

  representative_tips <- c(representative_tips, rep_tip)
  drop_tips <- c(drop_tips, other_tips)
  records[[sp]] <- data.frame(
    species = sp,
    representative_tip = rep_tip,
    n_strains_in_species = nrow(rows),
    within_species_max_patristic_dist = within_max_dist,
    quality_value = rep_info$quality
  )
}

drop_tips <- c(drop_tips, tree_tips_not_in_map)
pruned <- drop.tip(tree, drop_tips)
stopifnot(Ntip(pruned) == length(representative_tips))

tip_to_species <- setNames(names(species_list), representative_tips)
# guard against duplicate species names created by relabeling (shouldn't
# happen given split() keys are unique, but fail loudly rather than emit an
# ambiguous tree if it ever does)
new_labels <- tip_to_species[pruned$tip.label]
stopifnot(!any(is.na(new_labels)))
stopifnot(!any(duplicated(new_labels)))
pruned$tip.label <- unname(new_labels)

write.tree(pruned, opt[["out-tree"]])

out_map <- do.call(rbind, records)
out_map <- out_map[order(out_map$species), ]
write.csv(out_map, opt[["out-map"]], row.names = FALSE)

cat(sprintf(
  "Pruned %d-strain tree to %d-species tree.\n  Species tree: %s\n  Representative-strain map: %s\n",
  Ntip(tree), Ntip(pruned), opt[["out-tree"]], opt[["out-map"]]
))
n_single <- sum(out_map$n_strains_in_species == 1)
cat(sprintf("%d/%d species represented by a single strain (no within-species distance diagnostic).\n",
            n_single, nrow(out_map)))
