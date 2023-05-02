#!/bin/bash
#SBATCH -J Genetic_correlation
#SBATCH -q batch -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:30:00
#SBATCH --mem=128GB

#Example
anno=anno_jabs_2023-04-30.txt
geno=geno_jabs_2023-04-30.txt.gz
pheno=pheno_jabs_2023-04-30.txt

module load singularity
singularity exec -B/projects ~/genetic_correlation/Singularity/Gemma.sif \
gemma -g ~/genetic_correlation/example/${geno} \
-p ~/genetic_correlation/example/${pheno} \
-a ~/genetic_correlation/example/${anno} \
-gk -o jabs_relatednessmat


