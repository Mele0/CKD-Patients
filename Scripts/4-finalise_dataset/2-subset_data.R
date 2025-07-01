# save subset of datas for included sample

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
final_dataset_dir <- args[2]
ext_rec_output_dir <- args[3]

source("../03-extraction_and_recoding/parameters/constants.R")

setwd(final_dataset_dir)

sample_eids_filename <- paste0("sample_eids_", INCIDENCE_YEAR_THRESHOLD, "years.rds")
file.copy(file.path(data_dir, sample_eids_filename), sample_eids_filename, overwrite = TRUE)
sample_eids <- readRDS(sample_eids_filename)
print("Copied sample_eids to final_dataset/")

proteomics <- readRDS(file.path(data_dir, "Proteins_nd_imputed.rds"))
scaled_proteomics_subset <- scale(proteomics[sample_eids,], center = TRUE, scale = TRUE) # mean 0, sd 1
saveRDS(scaled_proteomics_subset, paste0("proteins_scaled_", INCIDENCE_YEAR_THRESHOLD, "years.rds"))
saveRDS(attributes(scaled_proteomics_subset), paste0("scaling_params.rds"))
print(paste("Scaled and saved proteomics to final_dataset/ -> final dims:", paste(dim(scaled_proteomics_subset), collapse = ", ")))

covars <- readRDS(file.path(ext_rec_output_dir, "ukb_covars_final.rds"))
status <- readRDS(file.path(data_dir, paste0("ckd_status_", INCIDENCE_YEAR_THRESHOLD, "years.rds")))
df <- cbind(status[sample_eids, c("recr_date", "ckd_date", "death_date")], covars[sample_eids,])
saveRDS(df, paste0("covars_", INCIDENCE_YEAR_THRESHOLD, "years.rds"))
print(paste("Saved covars to final_dataset/ -> final dims:", paste(dim(covars[sample_eids,]), collapse = ", ")))

comorbs <- readRDS(file.path(ext_rec_output_dir, "ukb_comorbidities.rds"))
saveRDS(comorbs[sample_eids,], paste0("comorbidities_", INCIDENCE_YEAR_THRESHOLD, "years.rds"))
print(paste("Saved comorbidities to final_dataset/ -> final dims:", paste(dim(comorbs[sample_eids,]), collapse = ", ")))

