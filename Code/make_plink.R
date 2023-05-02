#This file shouldn't be run from the terminal since the code creates plink files that get deleted after exiting the R session!

setwd("~/genetic_correlation/")
libs <- c("genio", "argparse")
sapply(libs, require, character.only = TRUE)

geno <- read.table("example/geno_jabs_2023-04-30.txt", header = FALSE, sep = "\t")
anno <- read.table("example/anno_jabs_2023-04-30.txt", header = FALSE, sep = "\t")
pheno <- read.table("example/pheno_jabs_2023-04-30.txt", header = FALSE, sep = "\t")
data <- read.csv("example/phenotypes_jabs_2023-04-30.csv", header = TRUE, stringsAsFactors = TRUE)
genotypes <- read.csv("example/strains_genotypes_all.csv", header = TRUE)
name <- "jabs_2023-04-30"

names(data)[names(data) == c("MouseID")] <- "IID"
names(data)[names(data) == c("Strain")] <- "FID"
data$FID <- as.numeric(as.factor(data$FID))
data$IID <- as.numeric(as.factor(data$IID))
data$IID <- seq(1:nrow(data))
colnames(data)[3:ncol(data)] <- seq(1,ncol(data)-2)
data[,-c(1,2)] <- scale(data[,-c(1,2)], center = TRUE, scale = TRUE)
write.table(data, paste0("pheno_", name, "_plink"), row.names=FALSE, col.names=FALSE)

m <- nrow(geno) #number of loci
n <- ncol(geno) - 3 #nrow(pheno)

bim <- make_bim(n = m)
bim$chr <- anno[,1]
bim$id <- paste0(anno[,3])
bim$ref <- genotypes[,"minor"]
bim$alt <- genotypes[,"major"]
bim$pos <- genotypes[,"bp38"]

fam <- make_fam(n = n)
fam$fam <- as.numeric(as.factor(data$FID))
fam$id <- as.numeric(as.factor(data$IID)) 
fam$id <- seq(1:nrow(fam))
fam$sex <- sample(1:2, n, replace = TRUE)
fam$pheno <- as.factor(rep(-9, nrow(fam))) 


X <- ((as.matrix(geno[, -seq(1,3)]))) 
rownames(X) <- bim$id
colnames(X) <- fam$id

#setwd(getwd())
file_plink <- tempfile(paste0("jabs"))

# Write genotypes, along with the BIM and FAM files we created.
# Omiting them would result in writing the original dummy version of these tables, before we edited them.
time_write_genio <- system.time(
    write_plink(file_plink, X, bim, fam)
)

print("Done!")

#The files are written here "/tmp"!
#Rename and copy them to ~/genetic_correlation/example. 
#cd /tmp/Rtmp5xb04r/
#mv jabs52c26564d2c5e.bed jabs_2023-04-30.bed
#mv jabs52c26564d2c5e.bim jabs_2023-04-30.bim
#mv jabs52c26564d2c5e.fam jabs_2023-04-30.fam
#mv jabs_* ~/genetic_correlation/example/.





