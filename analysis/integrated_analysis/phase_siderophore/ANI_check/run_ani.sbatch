#!/bin/bash
#SBATCH --job-name=tfcn_ani
#SBATCH --output=ani_%j.log
#SBATCH --error=ani_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=stajichlab

set -euo pipefail
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/fastani/1.34/bin:$PATH
MODIR=/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/analysis/integrated_analysis/phase_siderophore/ANI_check

cd "${SCRATCH:?}"

# copy the three query genomes to node-local scratch for fast I/O
for g in Rhodotorula_paludigena_TFCN_1A-1-2.scaffolds.fa \
         Rhodotorula_pacifica_TFCN_1A-1-3.scaffolds.fa \
         Rhodotorula_pacifica_TFCN_1B-1-2.scaffolds.fa; do
    cp "$MODIR/$g" .
done
ls Rhodotorula_*.scaffolds.fa > queries.list

fastANI --ql queries.list --rl "$MODIR/references.list" --matrix -t 16 \
        -o "$MODIR/ANI.out" > "$MODIR/ANI.log" 2>&1
echo "fastANI exit: $?"