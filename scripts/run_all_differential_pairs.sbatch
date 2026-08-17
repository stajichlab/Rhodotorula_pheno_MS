#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH -c 2 --mem 4gb
#SBATCH --time=00:40:00
#SBATCH --array=1-110
#SBATCH -J diff_features
#SBATCH -o analysis/differential_features/logs/%A_%a.log

# Runs scripts/differential_features_by_species.py for every species pair
# (both fractions) listed in analysis/differential_features/pairs_to_run.tsv
# -- one SLURM array task per pair, so a slow/failed pair doesn't block the
# rest and everything fits comfortably in the short queue's 2h limit even
# though there are ~110 pairs total. Each task runs the base Mann-Whitney
# test plus TWO block-permutation checks (see block_permutation.py and
# differential_features_by_species.py's docstrings):
#   - plate:  --block-by "Library Plate", always informative, catches an MS
#             batch confound.
#   - clade:  --tree $TREE_PATH, catches a phylogenetic-relatedness
#             confound -- but is a *degenerate* (zero-power, not zero-
#             signal) test whenever the two species being compared don't
#             share any clade at the chosen resolution, which is the
#             common case for two well-separated species (see the
#             DEGENERATE-flagged diobovata-vs-mucilaginosa run). Each
#             task's log states plainly whether its clade test was
#             degenerate -- check that before trusting a clade result,
#             don't assume "0 significant" means "no real difference."
#
# TREE_PATH below is the 276-taxon FastTree-with-support-values tree
# (BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/...).
# Checked directly across all 110 pairs at N_CLADES=30: only 2 (both
# fractions of R. mucilaginosa vs R. pacifica) have any clade containing
# both species -- the other 108 will report DEGENERATE for the clade test,
# which is expected (see block_permutation.py) and not a bug: most species
# pairs here are cleanly monophyletic and non-overlapping in this tree, so
# clade-blocking has nothing to test for them. If a future tree changes
# this, update TREE_PATH in place (and rerun
# scripts/generate_species_pairs.py + resubmit) -- nothing else here needs
# to change.
#
# Before submitting:
#   1. Regenerate the pair list if the linked sample metadata has changed:
#        python3 scripts/generate_species_pairs.py
#   2. Make sure --array matches the current line count, e.g.:
#        sbatch --array=1-$(wc -l < analysis/differential_features/pairs_to_run.tsv) \
#            scripts/run_all_differential_pairs.sbatch
#      (a command-line --array overrides the #SBATCH default above, so this
#      is the safest way to submit if the pair count has ever changed since
#      this file was written).
#
# Run from the repo root: sbatch scripts/run_all_differential_pairs.sbatch

set -euo pipefail

REPO_DIR=/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS
cd "$REPO_DIR"

TREE_PATH=BFD/results/phyling_pep/protein/buildtree/fungi_odb10/fasttree/protein-Rhodotorula-taxa_278.fungi_odb10.fasttree.support.treefile
N_CLADES=30

module load miniconda3/py39_4.12.0

PAIRS_FILE=analysis/differential_features/pairs_to_run.tsv
if [[ ! -s "$PAIRS_FILE" ]]; then
    echo "missing $PAIRS_FILE -- run scripts/generate_species_pairs.py first" >&2
    exit 1
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PAIRS_FILE")
if [[ -z "$LINE" ]]; then
    echo "no line $SLURM_ARRAY_TASK_ID in $PAIRS_FILE (array range doesn't match the file -- see header comment)" >&2
    exit 1
fi
IFS=$'\t' read -r FRACTION SPECIES_A SPECIES_B <<< "$LINE"

echo "task $SLURM_ARRAY_TASK_ID: fraction=$FRACTION  A='$SPECIES_A'  B='$SPECIES_B'"

TREE_ARGS=()
if [[ -s "$TREE_PATH" ]]; then
    TREE_ARGS=(--tree "$TREE_PATH" --n-clades "$N_CLADES")
else
    echo "note: $TREE_PATH not found -- running without clade-blocking" >&2
fi

python3 scripts/differential_features_by_species.py \
    --species-a "$SPECIES_A" --species-b "$SPECIES_B" --fraction "$FRACTION" \
    --block-by "Library Plate" "${TREE_ARGS[@]}"
