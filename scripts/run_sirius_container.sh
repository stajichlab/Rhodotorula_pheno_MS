#!/usr/bin/env bash
# Submit the SIRIUS chain (formula -> fingerprint -> canopus -> structures ->
# write-summaries) over analysis/sirius_annotation/sirius_targets.mgf, sharded
# for checkpointing but run STRICTLY SERIALLY (one shard, one `sirius` process,
# at a time -- no SLURM array). See knowledgebase in the sibling project,
# ../Rhodotorula_MS2_pheno_explore/knowledgebase/sirius.md, for the full
# operational history behind this design: in SIRIUS 6.3.12, EVERY subcommand
# (even bare `formula`) needs an active `sirius login`, and concurrent
# processes touching the shared ~/.sirius-6.3 login/refresh-token file break
# login for everyone -- reproduced multiple times in that project via
# parallel SLURM arrays. Full serialization is the only design confirmed safe.
#
# Before running, generate the inputs this expects:
#   python3 scripts/build_compound_summary.py      # compound_summary.tsv per comparison
#   python3 scripts/select_sirius_targets.py ...    # analysis/sirius_annotation/sirius_targets.csv
#   python3 scripts/export_sirius_targets_mgf.py    # analysis/sirius_annotation/sirius_targets.mgf
#
# And confirm login is alive first (this job will fail on its very first
# shard's `formula` call otherwise, not just canopus/structures):
#   module load singularity
#   singularity exec --bind /bigdata,/scratch \
#     /bigdata/stajichlab/shared/singularity/sirius-6.3.12-linux-x64.sif sirius login --show
#
# Usage (from the repo root):
#   scripts/run_sirius_container.sh [spectra_per_shard] [after_jobid]
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

SPECTRA_PER_SHARD="${1:-20}"
AFTER_JOBID="${2:-}"
SIF=/bigdata/stajichlab/shared/singularity/sirius-6.3.12-linux-x64.sif
MGF=analysis/sirius_annotation/sirius_targets.mgf
SHARD_DIR=analysis/sirius_annotation/shards
OUT_DIR=analysis/sirius_annotation/sirius_project_container
PIPELINE_SCRIPTS="$PWD/scripts/sirius_container_pipeline"

[[ -f "$SIF" ]] || { echo "Container image not found: $SIF" >&2; exit 1; }
[[ -f "$MGF" ]] || { echo "Input MGF not found: $MGF -- run scripts/export_sirius_targets_mgf.py first" >&2; exit 1; }

N_SPECTRA=$(grep -c '^BEGIN IONS' "$MGF")
N_SHARDS=$(( (N_SPECTRA + SPECTRA_PER_SHARD - 1) / SPECTRA_PER_SHARD ))
echo "Sharding $N_SPECTRA spectra into $N_SHARDS shards (~${SPECTRA_PER_SHARD}/shard)"

if [[ ! -d "$SHARD_DIR" ]] || [[ $(find "$SHARD_DIR" -maxdepth 1 -name 'shard_*.mgf' 2>/dev/null | wc -l) -ne "$N_SHARDS" ]]; then
  rm -rf "$SHARD_DIR"
  python3 "$PIPELINE_SCRIPTS/shard_mgf.py" --input "$MGF" --out-dir "$SHARD_DIR" --n-shards "$N_SHARDS"
else
  echo "Reusing existing shards in $SHARD_DIR"
fi

mkdir -p analysis/sirius_annotation/logs
DEP_ARGS=()
[[ -n "$AFTER_JOBID" ]] && DEP_ARGS=(--dependency="afterany:$AFTER_JOBID")

JOBID=$(sbatch --chdir="$PWD" "${DEP_ARGS[@]}" \
  --export=ALL,SHARD_DIR="$PWD/$SHARD_DIR",OUT_DIR="$PWD/$OUT_DIR",SIF="$SIF" \
  --output="analysis/sirius_annotation/logs/%x-%j.log" --error="analysis/sirius_annotation/logs/%x-%j.log" \
  --parsable "$PIPELINE_SCRIPTS/run_sirius_container_serial_full.sbatch")
echo "Submitted job $JOBID ($N_SHARDS shards, serial)${AFTER_JOBID:+, after $AFTER_JOBID}"
echo "After it completes, merge with:"
echo "  python3 $PIPELINE_SCRIPTS/merge_sirius_shards.py --shard-root $PWD/$OUT_DIR --out-dir $PWD/$OUT_DIR/merged"
