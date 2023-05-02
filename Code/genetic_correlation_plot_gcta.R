setwd("~/genetic_correlation/")
libs <- c("ggplot2", "corrplot","yaml")
sapply(libs, require, character.only = TRUE)

data <- read.csv("example/phenotypes_jabs_2023-04-30.csv")
features <- names(data)[-which(names(data) %in% c("Strain", "MouseID"))]

ind <- read.csv("example/xy_array.txt", header = FALSE)
names(ind) <- c("x", "y")

#Make sure all the correlations have been computed 
for (index in seq(nrow(ind))){
	tryCatch( 
		expr = { f <- read.delim(paste0("example/gcta_results/jabs_", ind[index,]$x, "_", ind[index,]$y, ".hsq"))
			}, 
		error = function(e) {
		message("File missing!")
		print(e)
	})
}

#Parse the files & extract the covariances
cor <- matrix(0, 12, 12)
cor_p <- matrix(0, 12, 12)
cor_vector_p <- as.numeric()

for (index in seq(nrow(ind))){
	if (!file.exists(paste0("example/gcta_results/jabs_", ind[index,]$x, "_", ind[index,]$y, ".hsq"))){
		cor_vector_p[index] <- NA
		next
	}
	cat("Index: x, y = ", ind[index,]$x, ",", ind[index,]$y, "\n")
	f <- read.delim(paste0("example/gcta_results/jabs_", ind[index,]$x, "_", ind[index,]$y, ".hsq"))
	cor_vector_p[index] <- as.numeric(gsub(" \\(one-tailed test\\)", "", f[16,2]))
}

cor_vector_p <- p.adjust(cor_vector_p, method = "fdr")

for (index in seq(nrow(ind))){
	if (!file.exists(paste0("example/gcta_results/jabs_", ind[index,]$x, "_", ind[index,]$y, ".hsq"))){
		cor_vector_p[index] <- NA
		next
	}
	cat("Index: x, y = ", ind[index,]$x, ",", ind[index,]$y, "\n")
	f <- read.delim(paste0("example/gcta_results/jabs_", ind[index,]$x, "_", ind[index,]$y, ".hsq"))
	cor[ind[index,]$x, ind[index,]$y] <- as.numeric(f[11, "Variance"])
	cor[ind[index,]$y, ind[index,]$x] <- as.numeric(f[11, "Variance"])

	cor_p[ind[index,]$x, ind[index,]$y] <- cor_vector_p[index]
	cor_p[ind[index,]$y, ind[index,]$x] <- cor_vector_p[index]

}

diag(cor) <- rep(1,12)

colnames(cor) <- features
rownames(cor) <- features

corrplot(cor, is.corr = FALSE, type = "lower", method = "pie", tl.col = "black", tl.cex = 0.5)
corrplot(cor, is.corr = FALSE, method = "square", tl.col = "black", tl.cex = 1.1, diag = TRUE, col.lim = c(-1,1), order = 'hclust', rect.col = "black", rect.lwd = 3, cl.cex = 1.25)
dev.print(pdf, "jabs_gc_clustered.pdf")


#cor_gemma <- cor #cor comes from genetic_correlation_plot.R
cor_gcta <- cor #cor comes from above

cor_comparison <- matrix(0, 12, 12)
cor_comparison_pval <- matrix(0, 12, 12)

cor_comparison[lower.tri(cor_comparison)] <- cor_gemma[lower.tri(cor_gemma)]
cor_comparison[upper.tri(cor_comparison)] <- cor_gcta[upper.tri(cor_gcta)]

diag(cor_comparison) <- 1
rownames(cor_comparison) <- features
colnames(cor_comparison) <- features

cor_comparison_pval <- cor_p
cor_comparison_pval[lower.tri(cor_comparison_pval)] <- 1
rownames(cor_comparison_pval) <- features
colnames(cor_comparison_pval) <- features

cor_comparison_pval[cor_comparison_pval == 0] <- 1
diag(cor_comparison_pval) <- 1


corrplot(cor_comparison, is.corr = FALSE, method = "square", tl.col = "black", tl.cex = 1.20, diag = TRUE, col.lim = c(-1,1), rect.col = "black", rect.lwd = 3, cl.cex = 1.20, p.mat = cor_comparison_pval, sig.level = 0.05, pch.cex = 1.9, insig = 'label_sig')

dev.print(pdf, "gemma_lt_gcta_ut_gcorr.pdf", width = 8, height = 8)