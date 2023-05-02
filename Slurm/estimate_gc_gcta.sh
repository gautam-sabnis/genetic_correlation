#!/bin/bash
#SBATCH -J Genetic_correlation
#SBATCH -q batch -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH --mem=256GB
#SBATCH --array=1-66

BATCH_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}{p;q;}" < "xy_array.txt") #xy_array is a txt file of dimension # of array jobs x 2 containing the indices$
x=${BATCH_LINE%,*} #Extracts the first element of BATCH_LINE
y=${BATCH_LINE#*,} #Extracts the last element of BATCH_LINE

module load singularity
singularity exec -B/projects ~/genetic_correlation/Singularity/gcta.sif \
/gcta/gcta-1.94.1-linux-kernel-3-x86_64/./gcta-1.94.1 --bfile jabs_2023-04-30 \
--autosome --make-grm --out ~/genetic_correlation/gcta_results/jabs 

module load singularity
singularity exec -B/projects ~/genetic_correlation/Singularity/gcta.sif \
/gcta/gcta-1.94.1-linux-kernel-3-x86_64/./gcta-1.94.1 --reml-bivar ${x} ${y} \
--grm ~/genetic_correlation/gcta_results/jabs --pheno ~/genetic_correlation/example/pheno_jabs_2023-04-30_plink --reml-bivar-lrt-rg 0 --out ~/genetic_correlation/gcta_results/jabs_${x}_${y}

