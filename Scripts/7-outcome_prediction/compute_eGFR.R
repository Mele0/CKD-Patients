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
