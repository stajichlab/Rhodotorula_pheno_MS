#!/usr/bin/env bash
# Submit the SIRIUS chain as a SLURM array job on the short queue.
# One shard per array task, strictly serial (%1 concurrency).
#
# Prerequisites:
#   python3 scripts/build_compound_summary.py
#   python3 scripts/differential_features_by_species.py
#   sirius login active (see run_sirius_container_serial_full.sbatch)
#
# Usage:
#   scripts/run_sirius_array.sh [spectra_per_shard]
#       (default spectra_per_shard=50, ~20 min each on short queue)
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

SPECTRA_PER_SHARD="${1:-50}"
SIF=/bigdata/stajichlab/shared/singularity/sirius-6.3.12-linux-x64.sif
PIPELINE_SCRIPTS="$PWD/scripts/sirius_container_pipeline"

[[ -f "$SIF" ]] || { echo "Container image not found: $SIF" >&2; exit 1; }

# Step 1: Generate comprehensive target set (all significant, minus
# those already in sirius_annotations.tsv)
python3 scripts/select_all_sirius_targets.py

# Step 2: Export MGF for all targets
python3 scripts/export_sirius_targets_mgf.py

# Step 3: Shard the MGF
MGF=analysis/sirius_annotation/sirius_targets.mgf
SHARD_DIR=analysis/sirius_annotation/shards
N_SPECTRA=$(grep -c '^BEGIN IONS' "$MGF")
N_SHARDS=$(( (N_SPECTRA + SPECTRA_PER_SHARD - 1) / SPECTRA_PER_SHARD ))
echo "Sharding $N_SPECTRA spectra into $N_SHARDS shards (~${SPECTRA_PER_SHARD}/shard)"
rm -rf "$SHARD_DIR"
python3 "$PIPELINE_SCRIPTS/shard_mgf.py" \
  --input "$MGF" --out-dir "$SHARD_DIR" --n-shards "$N_SHARDS"

# Step 4: Submit array job (one shard per task, %1 = strict serial)
OUT_DIR=analysis/sirius_annotation/sirius_results
mkdir -p analysis/sirius_annotation/logs

JOBID=$(sbatch --chdir="$PWD" \
  --array=0-$((N_SHARDS-1))%1 \
  --export=ALL,SHARD_DIR="$PWD/$SHARD_DIR",OUT_DIR="$PWD/$OUT_DIR",SIF="$SIF" \
  --output="analysis/sirius_annotation/logs/%A_%a.log" \
  --error="analysis/sirius_annotation/logs/%A_%a.log" \
  --parsable "$PIPELINE_SCRIPTS/run_sirius_array.sbatch")

echo "Submitted array job $JOBID ($N_SHARDS shards, 1 at a time, short queue)"
echo "After it completes, merge and import with:"
echo "  python3 $PIPELINE_SCRIPTS/merge_sirius_shards.py --shard-root $PWD/$OUT_DIR --out-dir $PWD/$OUT_DIR/merged"
echo "  python3 scripts/import_sirius_annotations.py"
