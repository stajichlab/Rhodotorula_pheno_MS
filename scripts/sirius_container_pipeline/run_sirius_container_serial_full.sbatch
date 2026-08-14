#!/usr/bin/env bash
# Copied verbatim (2026-08-14) from the sibling project's proven pipeline:
# ../Rhodotorula_MS2_pheno_explore/algorithms/sirius_container_pipeline/ --
# see that project's knowledgebase/sirius.md for the full operational
# history (login-serialization requirement, JVM AOT-cache crash, output
# layout gotchas) this design is distilled from. Not submitted directly --
# scripts/run_sirius_container.sh (this repo's own wrapper) sets
# SHARD_DIR/OUT_DIR/SIF and submits this.
#SBATCH --job-name=sirius_container
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.log
#
# Full SIRIUS chain (formula -> fingerprint -> canopus -> structures ->
# write-summaries), one shard at a time, STRICTLY SERIAL -- not an array. Not
# submitted directly; this repo's launcher (scripts/run_sirius_container.sh)
# sets SHARD_DIR/OUT_DIR/SIF and submits this.
#
# Why fully serial, not "offline steps in parallel + network steps serial"
# (an earlier, wrong design of this pipeline): in SIRIUS 6.3.12 (confirmed via
# this exact container), EVERY subcommand -- including bare `formula` alone --
# requires an active `sirius login`, not just canopus/structures. The
# project's older learning that formula/fingerprint are offline was true for
# SIRIUS 5.8.1 but does not hold for 6.3.12. Since every step touches the
# shared ~/.sirius-<ver> login/refresh-token file, ANY concurrent SIRIUS
# process (even 2, even different analyses sharing this project's $HOME)
# races on that file and breaks login for everyone -- reproduced twice in this
# project (see .living/learnings.md, 2026-08-11 and 2026-08-13 entries). Full
# serialization -- one `sirius` process touching the shared workspace at any
# moment, project-wide -- is the only design confirmed safe. Each analysis's
# launcher chains its job after the others via --dependency so this holds
# across analyses too, not just within one.
#
# ZODIAC is deliberately excluded (not in the chain below) -- it is SIRIUS's
# only <DATASET TOOL> (needs all compounds together for its cross-spectrum
# Gibbs-sampling network), so sharding would change its result. Every step run
# here is a <COMPOUND TOOL> computed independently per spectrum -- sharding
# itself (splitting which spectra land in which project space) is still
# correct, only the *execution* had to stop being concurrent.
#
# SIRIUS's default "preserve already-computed tool results" behavior means
# shards with partial progress from an earlier (broken) attempt resume from
# whatever tool last completed rather than recomputing from scratch.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"

: "${SHARD_DIR:?set via --export}"
: "${OUT_DIR:?set via --export}"
: "${SIF:?set via --export}"

module load singularity

SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-12}"
mkdir -p "$OUT_DIR"
N_SHARDS=$(find "$SHARD_DIR" -maxdepth 1 -name 'shard_*.mgf' | wc -l)
echo "Running full SIRIUS chain serially over $N_SHARDS shards in $SHARD_DIR"

for SHARD in "$SHARD_DIR"/shard_*.mgf; do
  SHARD_ID=$(basename "$SHARD" .mgf)
  SHARD_OUT="$OUT_DIR/$SHARD_ID"
  echo "=== $SHARD_ID ==="
  JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G -XX:AOTMode=off" singularity exec --bind /bigdata,/scratch "$SIF" \
    sirius --cores "${SLURM_CPUS_PER_TASK:-4}" \
    --input "$SHARD" --output "$SHARD_OUT" \
    formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
    fingerprint \
    canopus \
    structures \
    write-summaries
done

echo "All $N_SHARDS shards completed"
