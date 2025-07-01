# filter participants by inclusion/exclusion criteria

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_dir <- args[2]
ext_rec_output_dir <- args[3]
ckd_outcomes_path <- args[4]

source("../03-extraction_and_recoding/parameters/constants.R")

setwd(output_dir)
covars <- readRDS(file.path(ext_rec_output_dir, "ukb_covars_final.rds"))
n18_outcomes <- readRDS(ckd_outcomes_path)

# check valid merge: matching eids per row
print(paste("VALID MERGE? - check matching EIDs per row:", all(rownames(covars) == n18_outcomes$eid)))

library(tidyverse)

ckd_status_df <- data.frame(
  recr_date = n18_outcomes$date_recr,
  death_date = n18_outcomes$date_death, 
  ckd_date = n18_outcomes$date_diagnosis, 
  esrd_date = covars$esrd_date, 
  esrd_source = if_else(is.na(covars$esrd_source), NA, paste0("UKBB_", covars$esrd_source))
) 
rownames(ckd_status_df) <- n18_outcomes$eid

### fix column types
ckd_status_df <- ckd_status_df %>% mutate(
  esrd_source = as.factor(esrd_source),
)

# count valid HES ICD10 diagnosis dates
print(paste("Extracted", sum(!is.na(ckd_status_df$ckd_date)), "HES ICD10 diagnoses")) # 24986 

### inclusion/exclusion criteria for defining CKD cases
ckd_status_df <- ckd_status_df %>% mutate(
  ckd_flg = !is.na(ckd_date),
  eligible_flg = case_when(
    is.na(ckd_date) ~ FALSE, # exclude healthy controls
    ckd_date <= recr_date ~ FALSE, # exclude prevalent cases at recruitment
    year(ckd_date) - year(recr_date)  > INCIDENCE_YEAR_THRESHOLD ~ FALSE, # exclude incident cases >10 years after recruitment
    esrd_date == "1900-01-01" ~ FALSE, # has ESRD, but date is unknown (field 42026, coding 272)
    !is.na(esrd_date) & esrd_date < recr_date ~ FALSE, # exclude prevalent ESRD cases at recruitment
    !is.na(esrd_date) & esrd_date < ckd_date ~ FALSE, # exclude ESRD cases diagnosed before CKD (error in data)
    .default = TRUE
  )
)
print(paste("Identified", sum(ckd_status_df$eligible_flg == TRUE), "CKD cases based on inclusion criteria")) # 14013

### restrict sample based on proteomics availability 
proteomics <- readRDS(file.path(data_dir, "Proteins_nd_imputed.rds"))
sample_eids <- intersect(rownames(ckd_status_df[ckd_status_df$eligible_flg == TRUE,]), rownames(proteomics))
print(paste("Final sample with available proteomics data:", length(sample_eids))) # 2022
saveRDS(sample_eids, file.path(data_dir, paste0("sample_eids_", INCIDENCE_YEAR_THRESHOLD, "years.rds")))

# save recr/death/ckd year and sample inclusion status
ckd_status_df$final_sample_flg <- rownames(ckd_status_df) %in% sample_eids
saveRDS(ckd_status_df, file.path(data_dir, paste0("ckd_status_", INCIDENCE_YEAR_THRESHOLD, "years.rds")))

### Table1 - selection bias checks 
print("Creating Table1s")
library(table1)
library(knitr)
ckd_status_df <- ckd_status_df %>% mutate(
  subset = case_when(
    ckd_flg == FALSE ~ "No CKD diagnosis",
    eligible_flg == FALSE ~ "Excluded based on CKD/ESRD diagnosis date",
    final_sample_flg == FALSE ~ "Eligible but no proteomics data",
    final_sample_flg == TRUE ~ "Included in final sample",
  ),
  recr_year = as.numeric(year(recr_date)),
  ckd_year = as.numeric(year(ckd_date)),
  esrd_year = as.numeric(year(esrd_date)),
  death_year = as.numeric(year(death_date)),
  years_diag_to_esrd = as.numeric((esrd_date - ckd_date) / 365.25),
  years_diag_to_death = as.numeric((death_date - ckd_date) / 365.25),
)
ckd_status_df$subset <- factor(ckd_status_df$subset, levels = c("No CKD diagnosis", "Excluded based on CKD/ESRD diagnosis date", "Eligible but no proteomics data", "Included in final sample"))
outcomes_t1 <- table1::table1(~ recr_year + ckd_year + esrd_year + esrd_year + years_diag_to_esrd + death_year + years_diag_to_death | subset, ckd_status_df) 
sink("table1_outcomes_dist.md")
kable(outcomes_t1)
sink(NULL)
print("Saved table1_outcomes_dist.md")

covars$subset <- ckd_status_df$subset
covars_t1 <- table1::table1(~ sex + yob + ethnicity + household_income + smoking_status + alcohol_freq + diabetes | subset, covars) # TODO use significant covars from lasso
sink("table1_covars_dist.md")
kable(covars_t1)
sink(NULL)
print("Saved table1_covars_dist.md")


### Comparing different incidence year threshold cutoffs
print("Testing different INCIDENCE_YEAR_THRESHOLD values")
### inclusion/exclusion criteria for defining CKD cases
test_thresholds <- sapply(0:20, function(INCIDENCE_YEAR_THRESHOLD) {
  ckd_status_df <- ckd_status_df %>% mutate(
    ckd_flg = !is.na(ckd_date),
    eligible_flg = case_when(
      is.na(ckd_date) ~ FALSE, # exclude healthy controls
      ckd_date <= recr_date ~ FALSE, # exclude prevalent cases at recruitment
      year(ckd_date) - year(recr_date)  > INCIDENCE_YEAR_THRESHOLD ~ FALSE, # exclude incident cases >10 years after recruitment
      esrd_date == "1900-01-01" ~ FALSE, # has ESRD, but date is unknown (field 42026, coding 272)
      !is.na(esrd_date) & esrd_date < recr_date ~ FALSE, # exclude prevalent ESRD cases at recruitment
      !is.na(esrd_date) & esrd_date < ckd_date ~ FALSE, # exclude ESRD cases diagnosed before CKD (error in data)
      .default = TRUE
    )
  )
  incl_count <- sum(ckd_status_df$eligible_flg == TRUE)
  sample_eids <- intersect(rownames(ckd_status_df[ckd_status_df$eligible_flg == TRUE,]), rownames(proteomics))
  sample_size <- length(sample_eids)
  return(c(INCIDENCE_YEAR_THRESHOLD, incl_count, sample_size))
})

print("Saving plots for INCIDENCE_YEAR_THRESHOLD tests")
test_thresholds <- data.frame(t(test_thresholds))
colnames(test_thresholds) <- c("threshold", "incident_ckd", "sample_size")
ggplot(test_thresholds, aes(x = threshold, y = sample_size, label = sample_size)) + geom_point() + geom_text(hjust = -0.3, vjust=0.1, angle = -40) + ggtitle("Incidence year threshold vs final sample size")
ggsave(file.path(output_dir, "threshold_v_sample_size.png"))
ggplot(pivot_longer(test_thresholds, c("incident_ckd", "sample_size")), aes(x = threshold, colour = as.factor(name), y = value, label = value)) + geom_point() + geom_line() + geom_text(hjust = -0.3, vjust=0.1, angle = -40) + labs(colour = "# Participants") + ggtitle("Effects on eligible CKD incident cases\nand final sample size")
ggsave(file.path(output_dir, "threshold_v_eligible_and_sample_lineplot.png"))
ggplot(pivot_longer(test_thresholds, c("incident_ckd", "sample_size")), aes(x = threshold, fill = as.factor(name), y = value, label = value)) + geom_bar(stat = "identity") + labs(colour = "# Participants") + ggtitle("Effects on eligible CKD incident cases\nand final sample size")
ggsave(file.path(output_dir, "threshold_v_eligible_and_sample_barplot.png"))
