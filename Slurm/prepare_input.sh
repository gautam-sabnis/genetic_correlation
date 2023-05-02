#!/bin/bash
#SBATCH -J prepare_input
#SBATCH -q batch -p compute
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 00:30:00
#SBATCH --mem=300GB

module load singularity
singularity exec -B/projects ~/genetic_correlation/Singularity/mousegwas.sif \
Rscript ~/genetic_correlation/Code/prepare_input.R 

