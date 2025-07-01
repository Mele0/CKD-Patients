# combines diagnosis information from HES registries and self-reported condition in UKBB verbal interview (data field 20002)
# requires access to cases defined by ICD10 from HES data in outcome_definition output file, output_final.rds
# requires access to self-reported illnesses and respective diagnosis years from extraction and recoding output file, ukb_extracted.rds
# specify conditions in codeList, using UKBB data coding 6

source("constants.R")
library(tidyverse)

# search specified field_name for self-reported illnesses (eg field 20002 extracted as self_reported_illness) 
# match with any code in codeList 
# return earliest reported valid year of diagnosis (eg field 20008 extracted as self_reported_illness_year)
# return NA if participant has never reported this illness ("healthy"), or if no valid dates reported
min_self_rep_year <- function(ukb_extracted, field_name, codeList, instanceList=0:3) {
  # get earliest self-reported year for each instance
  years <- sapply(instanceList, function(i) {
    # self-reported illnesses, data field 20002, data coding 6
    illness_codes <- t(apply( 
        ukb_extracted %>% select(starts_with(paste0(field_name, ".", i, "."))), 
        1, as.array
      )) # row = array of illness codes for each participant
    
    # year of first diagnosis, data field 20008, data coding 13
    illness_years <- ukb_extracted %>% 
      select(starts_with(paste0(field_name, "_year.", i, ".")))
    illness_years[illness_years == -1 | illness_years == -3] <- NA # remove invalid dates
    # print(paste("removed", sum(any(illness_years == -1, illness_years == -3)), "invalid years for instance", i))
    illness_years <- floor(illness_years) # remove decimal (see field 20008 notes: year 1970 is recorded as 1970.5)
    illness_years <- t(apply( 
        illness_years, 
        1, as.array
      )) # convert each row to array of corresponding diagnosis years for each participant
    
    return(
      sapply(1:nrow(illness_codes), function(pt) { # for each participant
          pt_years <- sapply(codeList, function(code) # match for any illness code in the codeList
            illness_years[pt, match(code, illness_codes[pt,])]) # get corresponding year of diagnosis
          return(
            if_else(
              sum(!is.na(pt_years)) == 0, 
              NA, # no valid year of diagnoses, return NA
              min(pt_years, na.rm = TRUE) # return earliest valid year
            )
          )
        }
      )
    )
  }
  )
  # get earliest reported year of diagnosis across all instances
  min_year <- apply(years, 1, function(x) if_else(
      sum(!is.na(x)) == 0, 
      NA, # NA if no year of diagnosis
      min(x, na.rm=TRUE)
    )) 
  
  return(min_year)
}

# compare diagnosis year from HES ICD10 and self-reported illness
# return dataframe with 
# - earliest year of diagnosis 
# - source of information (HES or self-reported)
# - case status
diag_info = function(self_rep, hes, recr, disease_name="diagnosis") {
  df = data.frame(self_reported_year = as.numeric(self_rep), hes_year = as.numeric(hes), recr_year = as.numeric(recr))
  df <- df %>% mutate(
    year = case_when(
        is.na(self_reported_year) & is.na(hes_year) ~ NA,
        # !is.na(self_reported_year) & !is.na(hes_year) ~ min(self_reported_year, hes_year),
        !is.na(self_reported_year) & !is.na(hes_year) & self_reported_year < hes_year ~ self_reported_year,
        !is.na(self_reported_year) & !is.na(hes_year) & self_reported_year >= hes_year ~ hes_year,
        is.na(self_reported_year) & !is.na(hes_year) ~ hes_year,
        !is.na(self_reported_year) & is.na(hes_year) ~ self_reported_year
      ),
    source = case_when(
        is.na(self_reported_year) & is.na(hes_year) ~ NA,
        !is.na(self_reported_year) & !is.na(hes_year) & self_reported_year < hes_year ~ "Self-reported",
        !is.na(self_reported_year) & !is.na(hes_year) & self_reported_year >= hes_year ~ "HES ICD10",
        is.na(self_reported_year) & !is.na(hes_year) ~ "HES ICD10", # diagnosis in Hospital Episode Statistics
        !is.na(self_reported_year) & is.na(hes_year) ~ "Self-reported" # Self-reported medical condition in verbal interview, field 20002
      ),
    status = (year < recr_year | year - recr_year <= INCIDENCE_YEAR_THRESHOLD), # case status = TRUE if prevalent or incident within 10 years
  )
  
  res <- data.frame(year=as.numeric(df$year), source=as.factor(df$source), status=as.factor(df$status))
  colnames(res) <- sapply(colnames(res), function(col) paste0(disease_name, "_", col))
  return(res)
}
