#!/bin/bash
#SBATCH --job-name=ani_t152
#SBATCH --output=ani_%j.log
#SBATCH --error=ani_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --partition=stajichlab

set -euo pipefail
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/fastani/1.34/bin:$PATH
cd "${SCRATCH:?}"
cp /bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/analysis/integrated_analysis/phase_siderophore/ANI_check/t152/Rhodotorula_toruloides_TFCN_152C-6.scaffolds.fa .
echo "Rhodotorula_toruloides_TFCN_152C-6.scaffolds.fa" > q.list
fastANI --ql q.list --rl /bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/analysis/integrated_analysis/phase_siderophore/ANI_check/t152/refs.list --matrix -t 16 -o /bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/analysis/integrated_analysis/phase_siderophore/ANI_check/t152/ANI.out
