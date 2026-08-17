#!/usr/bin/bash -l
#SBATCH -p short --mem 1gb

module load singularity
nextflow run stajichlab/nf_phyling -profile singularity_slurm,ucr_hpcc --seq_type cds \
    --input input/cds --prefix Rhodotorula --markerset basidiomycota_odb12 --outdir results/phyling_cds
