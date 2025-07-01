# OLINK: Denoise data using linear mixed models

rm(list=ls())

# Activate renv for lme4
renv::activate() # re-activate renv with lme4
renv::restore() # reproduce environment (install same package versions)

library(lme4)
library(dplyr)

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


proteins = readRDS("Proteins_na_filtered.rds")

dict_qc = readRDS(file.path(output_dir, "Dictionaries/QC_information_UKB.rds"))
dict_sample = readRDS(file.path(output_dir, "Dictionaries/UKB_sample_information.rds"))

# Create sample key with plate ID, well ID and processing start date

tech = dict_sample %>% select(PlateID, WellID)
rownames(tech) = dict_sample$eid

# Remove nuisance variance (plate, position and date)

print("Denoising random effects")
pb <- utils::txtProgressBar(style = 3)

denoised = vcov_res = NULL

for (k in 1:ncol(proteins)){
  print(colnames(proteins)[k])
  plate = dict_qc[make.names(dict_qc$Assay)==colnames(proteins)[k],]
  tech_spec = left_join(tech, plate %>% select(PlateID, Processing_StartDate), by = "PlateID")
  rownames(tech_spec) = rownames(tech)
  model1=lmer(proteins[,k]~(1|PlateID)+(1|WellID)+(1|Processing_StartDate), data=tech_spec,
              REML = FALSE,
              control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  vcov=as.data.frame(VarCorr(model1))$vcov
  names(vcov)=as.data.frame(VarCorr(model1))[,1]
  vcov=vcov[-length(vcov)]
  vcov_res = cbind(vcov_res, vcov)
  
  res = residuals(model1)
  names(res) = rownames(proteins)[!is.na(proteins[,k])]
  res = res[rownames(proteins)]
  denoised = cbind(denoised, res)
  
  utils::setTxtProgressBar(pb, k/ncol(proteins))
}
cat("\n")

colnames(denoised) = colnames(proteins)
rownames(denoised) = rownames(proteins)

saveRDS(denoised, "Proteins_denoised.rds")
saveRDS(vcov_res, "Random_effects_vcov.rds")

renv::deactivate()
