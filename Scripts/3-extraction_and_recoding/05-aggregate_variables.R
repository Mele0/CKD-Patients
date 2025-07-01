args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_dir <- args[2]

source("scripts/aggr_functions.R")

setwd(output_dir)

library(tidyverse)

recoded <- readRDS("ukb_recoded.rds")
aggregated <- data.frame(row.names = rownames(recoded))

keep_col_mask <- rep(TRUE, ncol(recoded))
names(keep_col_mask) <- colnames(recoded)

### handle columns of 100% NA values
# error in recoding process - recover data from ukb_extracted
print("Checking for empty columns (fully NA) in ukb_recoded")
all_na_mask_recoded <- apply(recoded, 2, function(x) sum(!is.na(x)) == 0)
print("Checking ukb_extracted for potential recovery")
extracted <- readRDS("ukb_extracted.rds")
all_na_mask_extracted <- apply(extracted, 2, function(x) sum(!is.na(x)) == 0)

# data wiped from ukb_recoded but still available in ukb_extracted
print(paste("Recovering data from", sum(all_na_mask_recoded & !all_na_mask_extracted), "fields"))
for (col in colnames(recoded)[all_na_mask_recoded & !all_na_mask_extracted]) {
  # recover from ukb_extracted and manually recode variables
  print(paste("Recovering", col, "from ukb_extracted"))
  if (col == "esrd_date.0.0") {
    recoded$esrd_date <- extracted[,col]
    keep_col_mask <- drop_col(keep_col_mask, "esrd_date.0.0")
    next
  }
  recoded[,col] <- case_when(
    # field 884/904/1160/1269/1279/3829/3839, coding 100291
    col %in% c("moderate_activity.0.0", "vigorous_activity.0.0", "sleep_dur.0.0", "smoking_home.0.0", "smoking_outside.0.0", "num_stillbirths.0.0", "num_miscarriages.0.0") ~ case_when(
      extracted[,col] == -1 ~ NA, # do not know
      extracted[,col] == -3 ~ NA, # prefer not to answer
      .default = extracted[,col]
    ),
    # field 2734, coding 100584
    col == "num_live_births.0.0" ~ case_when(
      extracted[,col] == -3 ~ NA, # prefer not to answer
      .default = extracted[,col]
    ),
    # field 20008/20006, coding 13
    startsWith(col, "self_reported_illness_year.") | startsWith(col, "self_reported_cancer_year.") ~ case_when(
      extracted[,col] == -1 ~ NA, # date uncertain or unknown
      extracted[,col] == -3 ~ NA, # prefer not to answer
      .default = extracted[,col]
    ),
    # field 22671/22672/22674/22675/22677/22678/22680/22681, coding 13
    startsWith(col, "mean_IMT_") | startsWith(col, "max_IMT_") ~ case_when(
      extracted[,col] == 0 ~ NA, # measure invalid
      .default = extracted[,col]
    ),
    .default = extracted[,col]
  )
}

# QC drop columns with no valid data extracted (100% NAs)
for (col in colnames(recoded)[all_na_mask_extracted]) {
  print(paste("Removing", col, "- no valid data found"))
  keep_col_mask[col] <- FALSE
}

### self-reported conditions -- drop these cols, alr processed in comorbidity_status.R
# self_reported_cancer	Cancer code, self-reported		0,1,2,3	Cancer code, self-reported	0,1,2,3,4,5	0
# self_reported_illness	Non-cancer illness code, self-reported		0,1,2,3	Non-cancer illness code, self-reported	0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33	0
# self_reported_cancer_year	Interpolated Year when cancer first diagnosed	years	0,1,2,3	Interpolated Year when cancer first diagnosed (years)	0,1,2,3,4,5	0
# self_reported_illness_year	Interpolated Year when non-cancer illness first diagnosed	years	0,1,2,3	Interpolated Year when non-cancer illness first diagnosed (years)	0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33	0
keep_col_mask <- drop_col(keep_col_mask, "self_reported_")
print("Removing self_reported_* columns (already aggregated in comorbidities.R)")

### causes of death
# death_cause_primary	Underlying (primary) cause of death: ICD10		0,1	Underlying (primary) cause of death: ICD10	0	0
# death_cause_secondary	Contributory (secondary) causes of death: ICD10		0,1	Contributory (secondary) causes of death: ICD10	0,1,2,3,4,5,6,7,8,9,10,11,12,13,14	0
# ICD10 codes recoded for data coding 19
print("Aggregating causes of death")
death_pri <- t(apply( 
    recoded[,startsWith(colnames(recoded), "death_cause_primary.")], 
    1, as.array
  )) # row = array of ICD10 codes for each participant
death_sec <- t(apply( 
    recoded[,startsWith(colnames(recoded), "death_cause_secondary.")], 
    1, as.array
  )) # row = array of ICD10 codes for each participant
# flag if cause of death is N18
aggregated <- aggregated %>% mutate(
  death_ckd_pri = apply(death_pri, 1, function(pt) any(startsWith(pt, "N18"), na.rm = TRUE)),
  death_ckd_sec = apply(death_sec, 1, function(pt) any(startsWith(pt, "N18"), na.rm = TRUE)),
  death_ckd_any = (death_ckd_pri | death_ckd_sec)
)
keep_col_mask <- drop_col(keep_col_mask, "death_cause_")


### measured at baseline, but list of answers/measurements provided
# pulse_rate	Pulse rate, automated reading	bpm	0	Pulse rate, automated reading (bpm)	0,1	
print("Aggregating pulse_rate")
aggregated$pulse_rate <- get_mean(recoded, "pulse_rate.")
keep_col_mask <- drop_col(keep_col_mask, "pulse_rate.")

# bp_diastolic	Diastolic blood pressure, automated reading	mmhg	0	Diastolic blood pressure, automated reading (mmhg)	0,1	
print("Aggregating bp_diastolic")
aggregated$bp_diastolic <- get_mean(recoded, "bp_diastolic.")
keep_col_mask <- drop_col(keep_col_mask, "bp_diastolic.")

# bp_systolic	Systolic blood pressure, automated reading	mmhg	0	Systolic blood pressure, automated reading (mmhg)	0,1	
print("Aggregating bp_systolic")
aggregated$bp_systolic <- get_mean(recoded, "bp_systolic.")
keep_col_mask <- drop_col(keep_col_mask, "bp_systolic.")

# TODO ????
# education	Qualifications		0	Qualifications	0,1,2,3,4,5	0
print("Aggregating education")
education_rankings <- c(
  "College or University degree",
  "Other professional qualifications eg: nursing, teaching",
  "NVQ or HND or HNC or equivalent",
  "A levels/AS levels or equivalent",
  "O levels/GCSEs or equivalent",
  "CSEs or equivalent",
  "Other"
)
get_highest_ed <- function(ed_list) {
  # array of recoded 100305 values for 1 participant
  if (sum(!is.na(ed_list)) == 0)
    return(NA)
  if (sum(!is.na(ed_list)) == 1)
    return(ed_list[!is.na(ed_list)][1])
  ed_list_ranks <- sapply(ed_list, function(ed) {
    if (is.na(ed))
      return(NA)
    else
      return(match(ed, education_rankings))
  })
  ed_highest_rank <- min(ed_list_ranks, na.rm = TRUE)
  return(education_rankings[ed_highest_rank])
}
aggregated$education <- apply( 
    recoded[,startsWith(colnames(recoded), "education.")], 
    1, function(x) get_highest_ed(as.array(x))
  )
aggregated$education <- as.factor(aggregated$education)
keep_col_mask <- drop_col(keep_col_mask, "education.")

# employment	Current employment status		0	Current employment status	0,1,2,3,4,5,6
print("Aggregating employment")
aggregated$paid_employment <- apply(
  recoded[startsWith(colnames(recoded), "employment.")], 1,
  function(x) {
    if ("In paid employment or self-employed" %in% as.array(x))
      return(TRUE)
    else 
      return(FALSE)
  }
)
keep_col_mask <- drop_col(keep_col_mask, "employment.")

# disability	Attendance/disability/mobility allowance		0	Attendance/disability/mobility allowance	0,1,2
print("Aggregating disability")
aggregated$disability_allowance <- get_valid_flg(recoded, "disability.", c("Attendance allowance", "Disability living allowance", "Blue badge"))
keep_col_mask <- drop_col(keep_col_mask, "disability.")

# vascular_diagnosis	Vascular/heart problems diagnosed by doctor		0	Vascular/heart problems diagnosed by doctor	0,1,2,3
print("Aggregating vascular_diagnosis")
aggregated$vascular_diagnosis <- get_valid_flg(recoded, "vascular_diagnosis.", c("Heart attack", "Angina", "Stroke", "High blood pressure"))
keep_col_mask <- drop_col(keep_col_mask, "vascular_diagnosis.")

# med_pain	Medication for pain relief, constipation, heartburn		0	Medication for pain relief, constipation, heartburn	0,1,2,3,4,5
print("Aggregating med_pain")
aggregated$med_pain <- apply(
  recoded[,startsWith(colnames(recoded), "med_pain.")], 1, 
  function(x) {
    arr <- as.array(x)
    has_nsaid <- "Ibuprofen" %in% arr
    has_ppi <- "Omeprazole" %in% arr
    return(case_when(
      has_nsaid & has_ppi ~ "NSAIDs + PPI",
      has_nsaid ~ "NSAIDs",
      has_ppi ~ "PPI",
      sum(!is.na(arr)) != 0 ~ "Others",
      .default = NA
    ))
  }
)
aggregated$med_pain <- as.factor(aggregated$med_pain)
keep_col_mask <- drop_col(keep_col_mask, "med_pain.")

# vitamin_supp	Vitamin and mineral supplements		0	Vitamin and mineral supplements	0,1,2,3,4,5,6
print("Aggregating vitamin_supp")
aggregated$vitamin_supp <- get_valid_flg(recoded, "vitamin_supp.", c("Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D", "Vitamin E", "Folic acid or Folate (Vit B9)", "Multivitamins +/- minerals"))
keep_col_mask <- drop_col(keep_col_mask, "vitamin_supp.")

# physactivity_type	Types of physical activity in last 4 weeks		0	Types of physical activity in last 4 weeks	0,1,2,3,4
print("Aggregating physactivity_type")
get_highest_phys <- function(phys_list) {
  if (sum(!is.na(phys_list)) == 0)
    return(NA)
  phys_list_ranks <- sapply(phys_list, function(phys) {
    if (is.na(phys))
      return(NA)
    else {
      if (phys %in% c("Walking for pleasure (not as a means of transport)", "Light DIY (eg: pruning, watering the lawn)"))
        return(1)
      else if (phys == "Other exercises (eg: swimming, cycling, keep fit, bowling)")
        return(2)
      else if (phys %in% c("Strenuous sports", "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)"))
        return(3)
      else
        return(NA)
    }
  })
  phys_highest_rank <- max(phys_list_ranks, na.rm = TRUE)
  return(c("Light", "Moderate", "Heavy")[phys_highest_rank])
}
aggregated$physactivity_type <- apply( 
    recoded[,startsWith(colnames(recoded), "physactivity_type.")], 
    1, function(x) get_highest_phys(as.array(x))
  )
aggregated$physactivity_type <- as.factor(aggregated$physactivity_type)
keep_col_mask <- drop_col(keep_col_mask, "physactivity_type.")

# med_cvddm	Medication for cholesterol, blood pressure or diabetes		0	Medication for cholesterol, blood pressure or diabetes	0,1,2
print("Aggregating med_cvddm")
aggregated$med_cvddm <- get_valid_flg(recoded, "med_cvddm.", c("Cholesterol lowering medication", "Blood pressure medication", "Insulin"))
keep_col_mask <- drop_col(keep_col_mask, "med_cvddm.")

# mineral_supp	Mineral and other dietary supplements		0	Mineral and other dietary supplements	0,1,2,3,4,5
print("Aggregating mineral_supp")
aggregated$med_cvddm <- get_valid_flg(recoded, "mineral_supp.", c("Fish oil (including cod liver oil)", "Glucosamine", "Calcium", "Zinc", "Iron", "Selenium"))
keep_col_mask <- drop_col(keep_col_mask, "mineral_supp.")

# special_diet	Type of special diet followed		0	Type of special diet followed	0,1,2,3,4,5
print("Aggregating special_diet")
aggregated$special_diet <- apply(
  recoded[,startsWith(colnames(recoded), "special_diet.")], 1, 
  function(x) {
    arr <- as.array(x)
    return(case_when(
      "Vegetarian" %in% arr | "Vegan" %in% arr ~ "Vegetarian/Vegan",
      sum(!is.na(arr)) != 0 ~ "Others",
      .default = NA
    ))
  }
)
aggregated$special_diet <- as.factor(aggregated$special_diet)
keep_col_mask <- drop_col(keep_col_mask, "special_diet.")

###  measured at baseline, and only for subset of participants
## get mean of means and mean of max
# mean_IMT_120	Mean carotid IMT (intima-medial thickness) at 120 degrees 	micrometres	2,3	Mean carotid IMT (intima-medial thickness) at 120 degrees  (micrometres)	0
# mean_IMT_150	Mean carotid IMT (intima-medial thickness) at 150 degrees 	micrometres	2,3	Mean carotid IMT (intima-medial thickness) at 150 degrees  (micrometres)	0
# mean_IMT_210	Mean carotid IMT (intima-medial thickness) at 210 degrees 	micrometres	2,3	Mean carotid IMT (intima-medial thickness) at 210 degrees  (micrometres)	0
# mean_IMT_240	Mean carotid IMT (intima-medial thickness) at 240 degrees 	micrometres	2,3	Mean carotid IMT (intima-medial thickness) at 240 degrees  (micrometres)	0
print("Aggregating mean_IMT_[120|150|210|240]")
aggregated$mean_IMT_avg <- get_mean(recoded, "mean_IMT_")
keep_col_mask <- drop_col(keep_col_mask, "mean_IMT_")

# max_IMT_120	Maximum carotid IMT (intima-medial thickness) at 120 degrees 	micrometres	2,3	Maximum carotid IMT (intima-medial thickness) at 120 degrees  (micrometres)	0
# max_IMT_150	Maximum carotid IMT (intima-medial thickness) at 150 degrees 	micrometres	2,3	Maximum carotid IMT (intima-medial thickness) at 150 degrees  (micrometres)	0
# max_IMT_210	Maximum carotid IMT (intima-medial thickness) at 210 degrees 	micrometres	2,3	Maximum carotid IMT (intima-medial thickness) at 210 degrees  (micrometres)	0
# max_IMT_240	Maximum carotid IMT (intima-medial thickness) at 240 degrees 	micrometres	2,3	Maximum carotid IMT (intima-medial thickness) at 240 degrees  (micrometres)	0
print("Aggregating max_IMT_[120|150|210|240]")
aggregated$mean_IMT_avg <- get_mean(recoded, "max_IMT_")
keep_col_mask <- drop_col(keep_col_mask, "max_IMT_")

### clean up 
print("Cleaning up and saving final dataset")
# drop columns that have been aggregated
unchanged_cols <- recoded[,keep_col_mask]
# combine all columns
final <- cbind(aggregated, unchanged_cols)
# remove instance/array indices from remaining column names
colnames(final) <- gsub(".[0-9]+.[0-9]+$", "", colnames(final)) 

# save final covars dataset
saveRDS(final, "ukb_covars_final.rds")