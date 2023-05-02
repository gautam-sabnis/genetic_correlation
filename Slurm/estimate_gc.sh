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

#Example
anno=anno_jabs_2023-04-30.txt
geno=geno_jabs_2023-04-30.txt.gz
pheno=pheno_jabs_2023-04-30.txt

module load singularity
singularity exec -B/projects ~/genetic_correlation/Singularity/Gemma.sif \
gemma -g ~/genetic_correlation/example/${geno} \
-p ~/genetic_correlation/example/${pheno} \
-n ${x} ${y}  -a ~/genetic_correlation/example/${anno} -k \
~/genetic_correlation/example/jabs_relatednessmat.cXX.txt \
-lmm -o ~/genetic_correlation/gemma_results/jabs_${x}_${y}

