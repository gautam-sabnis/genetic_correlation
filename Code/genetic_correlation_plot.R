setwd("~/genetic_correlation")
libs <- c("ggplot2", "corrplot","yaml")
sapply(libs, require, character.only = TRUE)

data <- read.csv("example/phenotypes_jabs_2023-04-30.csv")
features <- names(data)[-which(names(data) %in% c("Strain", "MouseID"))]

ind <- read.csv("example/xy_array.txt", header = FALSE)
names(ind) <- c("x", "y")

#Make sure all the correlations have been computed 
for (index in seq(nrow(ind))){
	tryCatch( 
		expr = { f <- read.delim(paste0("example/gemma_results/jabs_", ind[index,]$x, "_", ind[index,]$y, ".log.txt"))
			}, 
		error = function(e) {
		message("File missing!")
		print(e)
	})
}

#Parse the files & extract the covariances
cor <- matrix(0, 12, 12)
for (index in seq(nrow(ind))){
	if (index %in% c()){
		next
	}
	cat("Index: x, y = ", ind[index,]$x, ",", ind[index,]$y, "\n")
	f <- read.delim(paste0("example/gemma_results/jabs_", ind[index,]$x, "_", ind[index,]$y, ".log.txt"))$X..
	line <- grep("## REMLE estimate for Vg in the null model", f)
	cor[ind[index,]$x, ind[index,]$y] <- as.numeric(f[line + 2])	
	cor[ind[index,]$y, ind[index,]$x] <- as.numeric(f[line + 2])
}


cor[cor > 1] <- 1
diag(cor) <- rep(1,12)

colnames(cor) <- features
rownames(cor) <- features

corrplot(cor, is.corr = FALSE, type = "lower", method = "pie", tl.col = "black", tl.cex = 0.5)
corrplot::corrplot(cor, is.corr = FALSE, method = "square", tl.col = "black", tl.cex = 1.1, diag = TRUE, col.lim = c(-1,1), order = 'hclust', rect.col = "black", rect.lwd = 3, cl.cex = 1.25)
dev.print(pdf, "jabs_gc_clustered.pdf")

cor_gemma <- cor
