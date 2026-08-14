#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH -c 2 --mem 4gb
#SBATCH --time=00:20:00
#SBATCH --array=1-110
#SBATCH -J diff_features
#SBATCH -o analysis/differential_features/logs/%A_%a.log

# Runs scripts/differential_features_by_species.py for every species pair
# (both fractions) listed in analysis/differential_features/pairs_to_run.tsv
# -- one SLURM array task per pair, so a slow/failed pair doesn't block the
# rest and everything fits comfortably in the short queue's 2h limit even
# though there are ~110 pairs total.
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

python3 scripts/differential_features_by_species.py \
    --species-a "$SPECIES_A" --species-b "$SPECIES_B" --fraction "$FRACTION"
