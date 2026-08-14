#!/usr/bin/bash -l
#SBATCH -p short --mem 1gb

module load singularity
nextflow run stajichlab/nf_phyling -profile singularity_slurm,ucr_hpcc --seq_type protein \
    --input input/pep --prefix Rhodotorula --markerset fungi_odb10 --outdir results/phyling_pep
