### PCA of population based on protein expression

rm(list=ls())

## Set up data and output directories
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  data_dir <- args[1]
  print(paste("Data directory:", data_dir))
  output_dir <- args[2]
  print(paste("Output directory:", output_dir))
} else {
  stop(paste0("Invalid number of arguments provided (", length(args),
              ", expected 2). Please specify input and output directories."))
}

setwd(data_dir)
figures_dir <- file.path(output_dir, "Figures")

# Load packages
library(tidyverse)
library(FactoMineR)
library(colorspace)

# Load data
proteins = readRDS("Proteins_nd_imputed.rds")
print(ncol(proteins))

mypca=PCA(proteins, graph = FALSE)

saveRDS(mypca, "Proteins_PCA.rds")

ev=mypca$eig[,2]/100

# Scree plot
pdf(file.path(figures_dir, "PCA_scree_plot.pdf"), width=5, height=5)
par(mar=c(5,5,1,1))
plot(ev, pch=19, col="#0019a8", las=1, type="b", ylim=c(0,1),
     ylab="Proportion of explained variance", xlab="PCs")
points(cumsum(ev), pch=19, col="#dc251f", type="b")
legend("right", pch=19, col=c("#0019a8", "#dc251f"),
       legend=c("Proportion of e.v.", "Cumulative proportion of e.v."))
dev.off()