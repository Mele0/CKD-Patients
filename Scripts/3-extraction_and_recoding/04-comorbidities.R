# for comorbidities in outcome_definition/definitions,
# - date of diagnosis
# - source of data (HES/self-report)
# - case status: TRUE if prevalent / incident within INCIDENCE_YEAR_THRESHOLD

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_dir <- args[2]
def_dir <- args[3]
hes_outcomes_dir <- args[4]


source("scripts/comorb_diag_functions.R")
source("parameters/constants.R")

setwd(output_dir)

library(tidyverse)

extracted <- readRDS("ukb_extracted.rds")
comorb_list <- list.dirs(path = def_dir, full.names = FALSE, recursive = FALSE)

comorbs_df <- data.frame(row.names = rownames(extracted))
for (comorb in comorb_list) {
    if (comorb == "CKD") next # skip CKD

    # get self-reported field codingName and codes
    comorb_def_dir <- file.path(def_dir, comorb)
    filename <- list.files(path = comorb_def_dir, pattern = "self_reported_")
    if (length(filename) > 0) {
        filename <- filename[1]
        field_name <- gsub(".txt", "", filename)
        codes <- read_lines(file=file.path(comorb_def_dir, filename), skip=1)
        print(paste(comorb, "- getting", field_name, "with codes:", paste(codes, collapse=', ')))
        self_reported_year <- min_self_rep_year(extracted, field_name, codes)
    } else {
        print(paste(comorb, "- no self-reported data, skipping"))
        self_reported_year <- NA
    }

    # get date of diagnosis identified by ICD10 codes in hospital records
    print(paste(comorb, "- getting HES ICD10 outcomes"))
    hes_outcomes <- readRDS(file.path(hes_outcomes_dir, comorb, "output_final.rds"))
    hes_year <- year(hes_outcomes$date_diagnosis)

    # get date of recruitment
    recr_year <- year(hes_outcomes$date_recr)

    # get min and update df
    print(paste(comorb, "- combining data to get diagnosis year, data source, and case status"))
    diag <- diag_info(self_reported_year, hes_year, recr_year, comorb)
    comorbs_df <- cbind(comorbs_df, diag)
    
}

print("=======")
print(summary(comorbs_df))

# save full comorbs dataset
saveRDS(comorbs_df, "ukb_comorbidities.rds")