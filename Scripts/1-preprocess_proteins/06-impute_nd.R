# OLINK: Impute values below LOD

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

library(imputeLCMD)

proteins = readRDS("Proteins_nd_filtered.rds")

imp = impute.QRILC(proteins, tune.sigma = 1)
proteins_imp = imp[[1]]

saveRDS(proteins_imp, "Proteins_nd_imputed.rds")