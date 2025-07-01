rm(list=ls()) 

library(nnet) 
library(glmnet)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(sharp)
library(tableone)

############## Prepare comorbidities data

comorbids <- readRDS("comorbidities_10years.rds")

dim(comorbids) # 2022 66

# Select column names ending with 'status'
disease_status <- comorbids %>% select(ends_with("status"))

# Replace NA with FALSE in all columns
disease_status[] <- lapply(disease_status, function(col) {
  col <- ifelse(col == "TRUE", "TRUE", "FALSE")
  col[is.na(col)] <- "FALSE"  # Replace NAs explicitly
  return(col)
})
disease_status[] <- lapply(disease_status, function(col) {
  factor(col, levels = c("TRUE", "FALSE"))
})

# Remove columns with no cases
disease_status <- subset(disease_status, select = -c(maternal_DM_status, maternal_HTN_eclampsia_status))

# Summarize counts and proportions of "TRUE" and "FALSE" for each column
summary_df <- data.frame(
  Column = names(disease_status),
  TRUE_Count = sapply(disease_status, function(x) sum(x == "TRUE")),
  FALSE_Count = sapply(disease_status, function(x) sum(x == "FALSE")),
  TRUE_Proportion = sapply(disease_status, function(x) mean(x == "TRUE")),
  FALSE_Proportion = sapply(disease_status, function(x) mean(x == "FALSE"))
)

print(summary_df) # total 20 comorbidities conditions

dim(disease_status) # 2022 20
disease_status <- disease_status[order(rownames(disease_status)), ]

################### Prepare covars data

covars <- readRDS("imputed_var.rds")

dim(covars) # 2022 93

colnames(covars)[1] <- "ID"

# Set row names and remove the ID column
covars <- as.data.frame(covars)
rownames(covars) <- covars$ID
covars <- covars[, -1]

# Remove irrelevant columns
# 'yob' not needed as there is 'age_recruit'
# 'diabetes' not needed as this is diabetes diagnosed by doctor (field ID 2443)
# and we already have diabetes_status (self-reported and/or ICD10) in comorbids dataset
covars <- covars[, !(colnames(covars) %in% c(
  "death_ckd_pri", "death_ckd_sec", "death_ckd_any", 
  "yob", "diabetes", "recr_date", "ckd_date", 
  "death_date", "esrd_source", "esrd_date"
))]

dim(covars) # 2022 82
covars <- covars[order(rownames(covars)), ]

############################

# Combine covars and disease_status datasets
identical(rownames(covars), rownames(disease_status)) # TRUE
covars <- cbind(covars, disease_status)
dim(covars) # 2022 102

################ Calculate and add CKD stages as a variable

# Convert serum creatinine from micromol/L to mg/dl and update column name
covars <- covars %>% rename(creatinine_mg.dl = creatinine)
covars$creatinine_mg.dl <- covars$creatinine_mg.dl/88.42

# Calculcate eGFR using CKD-EPI 2021 equation, and add new column 'eGFR'
calculate_eGFR <- function(creatinine_mg.dl, age_recruit, sex) {
  
  # Set A and B values based on sex and creatinine level
  if (sex == "Female") {
    A <- 0.7
    B <- ifelse(creatinine_mg.dl <= 0.7, -0.241, -1.2)
    sex_multiplier <- 1.012
  } else if (sex == "Male") {
    A <- 0.9
    B <- ifelse(creatinine_mg.dl <= 0.9, -0.302, -1.2)
    sex_multiplier <- 1
  } 
  
  # Calculate eGFR
  eGFR <- 142 * ((creatinine_mg.dl / A) ^ B) * (0.9938 ^ age_recruit) * sex_multiplier
  
  return(eGFR)
}

covars <- covars %>%
  mutate(eGFR = mapply(calculate_eGFR, creatinine_mg.dl, age_recruit, sex))

# Classify CKD stages based on eGFR and add new column 'CKDstage_recruit'
classify_ckd_stage <- function(eGFR) {
  if (eGFR >= 90) {
    return("Stage 1")
  } else if (eGFR >= 60) {
    return("Stage 2")
  } else if (eGFR >= 45) {
    return("Stage 3a")
  } else if (eGFR >= 30) {
    return("Stage 3b")
  } else if (eGFR >= 15) {
    return("Stage 4")
  } else {
    return("Stage 5")
  }
}

covars <- covars %>%
  mutate(CKDstage_recruit = mapply(classify_ckd_stage, eGFR))

dim(covars) # 2022 104
str(covars)

# Change these categorical variables to factors
vars_to_factor <- c("education", "med_pain", "physactivity_type", "sex", "household_income", 
                    "insomnia", "snoring", "smoking_household", "oily_fish", "processed_meat", 
                    "salt_added", "alcohol_freq", "bodysize_age10", "height_age10", 
                    "multiple_birth", "maternal_smoking", "overall_health", "longstd_illness", 
                    "cancer_diagnosis", "menopause", "preg_loss", "ever_hrt", "smoking_status", 
                    "ethnicity", "CKDstage_recruit")

covars[vars_to_factor] <- lapply(covars[vars_to_factor], factor)
  
str(covars)


# Heatmap to show correlation amongst the variables
categorical_vars <- sapply(covars, is.factor) 
numerical_vars <- !categorical_vars  
covars_numeric <- covars[, numerical_vars]  # select only those columns that are numeric
pheatmap(cor(covars_numeric), cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = seq(-1,1, length.out = 100), border = NA)

################# Align with train_select (50%)
train_select <- readRDS("train_select.rds")
dim(train_select) # 1011 2

# Subset covars to match IDs in train_select
covars_train <- covars[rownames(covars) %in% train_select$ID, ]
dim(covars_train) # 1011 104

# Add cluster labels
covars_train <- cbind(cluster = train_select$cluster, covars_train)
covars_train$cluster <- as.factor(covars_train$cluster)

# Problem with NAs in logistic regression, so identify those columns
na_cols <- colnames(covars_train)[colSums(is.na(covars_train)) > 0]
print(na_cols)
na_counts <- colSums(is.na(covars_train))
na_counts[na_counts > 0]
# These are the 4 female-specific columns with 523 NAs (correspond with no of males)

# Change NA to "NotApplicable" for categorical columns
# i.e "menopause", "preg_loss", "ever_hrt"
for (col in na_cols) {
  if (is.factor(covars_train[[col]])) {
    covars_train[[col]] <- factor(covars_train[[col]], levels = c(levels(covars_train[[col]]), "NotApplicable"))
    covars_train[[col]][is.na(covars_train[[col]])] <- "NotApplicable"
  }
}
# check distribution of numeric column "num_live_births"
hist(covars_train$num_live_births)
summary(covars_train$num_live_births)
# Change to categorical with levels "0-2", ">=3" and NAs as "NotApplicable"
covars_train$num_live_births <- factor(
  ifelse(is.na(covars_train$num_live_births), "NotApplicable",
         ifelse(covars_train$num_live_births <= 2, "0-2", ">=3")),
  levels = c("0-2", ">=3", "NotApplicable")
)

# Given that stability-Lasso only select 1 variable for each cluster (cystatinC for cluster 0 and 1, platelet for cluster 2), we will divide variables into 4 domains and run Lasso for each domain

demographics_cols <- c("age_recruit", "sex", "ethnicity", "education", "paid_employment", 
                       "household_income", "IMD_England")

clinical_cols <- c("BMI", "waist_circ", "pulse_rate", "bp_diastolic", "bp_systolic", 
                   "smoking_status", "smoking_household", "smoking_home", "smoking_outside", 
                   "oily_fish", "processed_meat", "salt_added", "alcohol_freq", 
                   "physactivity_type", "moderate_activity", "vigorous_activity", 
                   "sleep_dur", "insomnia", "snoring", "med_pain", "vitamin_supp", 
                   "med_cvddm", "treatment_medcount", "bodysize_age10", "height_age10", 
                   "multiple_birth", "maternal_smoking", "overall_health", "menopause", 
                   "num_live_births", "preg_loss", "ever_hrt", "ts_ratio_adjusted", 
                   "eGFR", "CKDstage_recruit")

comorbidities_cols <- c("longstd_illness", "vascular_diagnosis", "cancer_diagnosis", 
                        "angina_status", "arrhythmia_status", "CVA_status", 
                        "cystic_renal_status", "diabetes_status", "ESRD_status", 
                        "gout_status", "heart_failure_status", "hypertension_status", 
                        "KU_congen_malform_status", "MI_status", "nephritic_syn_status", 
                        "nephrolithistatus", "nephrotic_syn_status", "obs_uropathy_status", 
                        "PVD_status", "renal_cancer_status", "SLE_status", "TIN_status")

covars_train_demographics <- covars_train[, c("cluster", demographics_cols)]
covars_train_clinical <- covars_train[, c("cluster", clinical_cols)]
covars_train_comorbidities <- covars_train[, c("cluster", comorbidities_cols)]

# Capture remaining columns for labs
all_selected_cols <- c("cluster", demographics_cols, clinical_cols, comorbidities_cols)
remaining_cols <- setdiff(names(covars_train), all_selected_cols)

covars_train_labs <- covars_train[, c("cluster", remaining_cols)]

# Check dataset dimensions 
dim(covars_train_demographics)
dim(covars_train_clinical)
dim(covars_train_comorbidities)
dim(covars_train_labs)

# Save datasets as RDS files
saveRDS(covars_train_demographics, file = "covars_train_demographics.rds")
saveRDS(covars_train_clinical, file = "covars_train_clinical.rds")
saveRDS(covars_train_comorbidities, file = "covars_train_comorbidities.rds")
saveRDS(covars_train_labs, file = "covars_train_labs.rds")
saveRDS(covars_train, file = "covars_train.rds")






# Generate Table 1 by clusters, using all 2022 participants
train_select <- readRDS("train_select.rds") # 1011 
train_refit <- readRDS("train_refit.rds") # 607
test <- readRDS("test.rds") # 404

covars$cluster <- coalesce(
  train_select$cluster[match(rownames(covars), train_select$ID)],
  train_refit$cluster[match(rownames(covars), train_refit$ID)],
  test$cluster[match(rownames(covars), test$ID)]
)

table(covars$cluster)

vars <- setdiff(names(covars), "cluster")

# Create the table by cluster
table1 <- CreateTableOne(vars = vars, strata = "cluster", data = covars)

# Print the summary table
print(table1, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE)

# Print the table to the PDF
pdf("table1_output.pdf", width = 8, height = 10)
print(table1, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE)

# Turn off the PDF device to save the file
dev.off()


