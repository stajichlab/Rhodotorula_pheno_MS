#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH -c 2 --mem 4gb
#SBATCH --time=00:40:00
#SBATCH --array=1-10
#SBATCH -J paired_sup_vs_cell
#SBATCH -o analysis/differential_features/logs/paired_%A_%a.log

# Runs scripts/differential_features_by_species.py --sup-vs-cell for every
# species listed in analysis/differential_features/paired_species_to_run.tsv
# -- one SLURM array task per species. Each species' run tests the paired
# (same-strain) cell-vs-supernatant contrast within that species only:
#   - Wilcoxon signed-rank on each strain's (cell, supernatant) TSS pair,
#     so phylogenetic/genotype identity within the species is held constant
#     and we're asking "what does this species release (supernatant) vs
#     retain (cell)".
#   - --block-by "Library Plate": within-plate sign-flip permutation (see
#     block_permutation.py's paired_block_permutation_pvalues), the paired
#     analogue of the unpaired runs' plate-blocking confound check.
#
# Before submitting:
#   1. Regenerate the species list (drops species with < --min-n paired
#      strains) if the linked sample metadata has changed:
#        python3 scripts/generate_species_pairs.py --paired --min-n 3
#   2. Make sure --array matches the current line count, e.g.:
#        sbatch --array=1-$(wc -l < analysis/differential_features/paired_species_to_run.tsv) \
#            scripts/run_all_paired_sup_vs_cell.sbatch
#
# Run from the repo root: sbatch scripts/run_all_paired_sup_vs_cell.sbatch

set -euo pipefail

REPO_DIR=/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS
cd "$REPO_DIR"

module load miniconda3/py39_4.12.0

SPECIES_FILE=analysis/differential_features/paired_species_to_run.tsv
if [[ ! -s "$SPECIES_FILE" ]]; then
    echo "missing $SPECIES_FILE -- run scripts/generate_species_pairs.py --paired first" >&2
    exit 1
fi

SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_FILE")
if [[ -z "$SPECIES" ]]; then
    echo "no line $SLURM_ARRAY_TASK_ID in $SPECIES_FILE (array range doesn't match the file -- see header comment)" >&2
    exit 1
fi

echo "task $SLURM_ARRAY_TASK_ID: species='$SPECIES'"

python3 scripts/differential_features_by_species.py \
    --sup-vs-cell --species-a "$SPECIES" \
    --block-by "Library Plate"
