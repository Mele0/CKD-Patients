# OLINK: Impute missing values (no QC information)

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

library(impute)
proteins = readRDS("Proteins_denoised.rds")

imp = impute.knn(as.matrix(proteins),k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=362436069)
proteins_imp = imp$data

saveRDS(proteins_imp, "Proteins_na_imputed.rds")