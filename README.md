### Estimate Genetic Correlations
Estimate genetic correlations using individual-level genotype data via *embarrassingly parallel* array jobs. 

#### Requirements
The pipeline is developed and tested in Linux and Mac OS environments. The following software and packages are required:

1. [R](https://www.r-project.org/)
2. [GEMMA](https://github.com/genetics-statistics/GEMMA)
3. [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview)
4. [mousegwas](https://github.com/TheJacksonLaboratory/mousegwas)
5. [genio](https://cran.r-project.org/web/packages/genio/index.html)
6. [python](https://www.python.org/)
7. [numpy](https://numpy.org/)
8. [pandas](https://pandas.pydata.org/)
9. [corrplot](https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html)
10. [Singularity](https://sylabs.io/docs/)

`mousegwas` comes preinstalled with several packages. This pipeline needs to be run using a Slurm cluster. 

### Tutorial

You can download the pipeline by 

    cd ~
    git clone https://github.com/gautam-sabnis/genetic_correlation
    cd genetic_correlation

To create the singularity images, run, for example

    singularity build --fakeroot gcta.sif gcta.def


You'll need the following **input** files to execute the pipeline:

- **CSV file**: A CSV file containing the traits, Strain and MouseID as columns. 

- **yaml**: Yaml file similar to the one in the example folder.
- **xy_array.txt**: An array of indices encoding all possible trait pairs.

To create `xy_array.txt`, run the following

        cd ~/genetic_correlation/example
        python3 Code/create_pairs.py --p example/pheno_jabs_2023-04-30.txt


#### 1. Prepare input files for GEMMA.
    Submit prepare_input.sh on the cluster.   

This runs prepare_input.R using the input files (*csv* and *yaml*). 

	Specify the args$input, args$yaml and name variables in prepare_input.R

This will create the following files (see the **example folder**),

- anno_name.txt 
- geno_name.txt  
- pheno_name.txt  
- genotypes_name.csv  
- phenotypes_name.csv  

**GEMMA** takes .txt files (**BIMBAM** files) as input. The .csv files are required later for converting **BIMBAM** files to **Plink** files. 

#### 2. Using GEMMA

    Submit estimate_gc_kinship.sh

This will estimate the kinship matrix. The estimated matrix is used, along with others, as input to the next script. 

    Submit estimate_gc.sh

This will estimate the binary genetic correlations for all trait pairs. 

#### 3. Using GCTA

**GCTA** takes Plink files (.bim, .bed, .fam) as input.

##### 3.1 Convert BIMBAM to Plink files.

    Run make_plink.R on the local machine.
Use the earlier created .txt and .csv files as input to make_plink.R. This will create the Plink files (`.bim, .bed, .fam`) in a temporary folder that you will need to rename and copy to your current working directory. See the example below - 

    R output
    Writing: /tmp/Rtmp5xb04r/jabs52c26564d2c5e.bed
    Writing: /tmp/Rtmp5xb04r/jabs52c26564d2c5e.bim
    Writing: /tmp/Rtmp5xb04r/jabs52c26564d2c5e.fam

    cd /tmp/Rtmp5xb04r/
    mv jabs52c26564d2c5e.bed jabs_2023-04-30.bed
    mv jabs52c26564d2c5e.bim jabs_2023-04-30.bim
    mv jabs52c26564d2c5e.fam jabs_2023-04-30.fam
    mv jabs_* ~/genetic_correlation/example/.

The Plink files, along with `pheno_jabs_2023-04-30_plink`, are used as inputs to the next script. 

	Submit estimate_gc_gcta.sh

Estimate bivariate genetic correlations using **GCTA** for all trait pairs.  

#### 4. Post-processing

Use `genetic_correlation_plot.R` file to
- extract the estimated parameters
- check for pairs for which the algorithm diverged.
- create a heatmap of the estimated bivariate genetic correlation matrix using **GEMMA**. Label this matrix as `cor_gemma`. 

Use `genetic_correlation_plot_gcta.R` to 
- extract the estimated parameters
- check for pairs for which the algorithm diverged.
- create a heatmap of the estimated genetic correlation matrix using **GCTA**. Label this matrix as `cor_gcta`. 
- Combine `cor_gemma` and `cor_gcta` in the lower and upper triangular parts of the matrix for comparison using `genetic_correlation_plot_gcta.R`. 


