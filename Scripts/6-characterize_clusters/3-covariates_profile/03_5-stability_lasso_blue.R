rm(list=ls())

library(nnet) 
library(glmnet)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(sharp)
library(tableone)
library(caret)
library(pROC)
library(ROCR)
library(fastDummies)

covars_train_demographics <- readRDS("covars_train_demographics.rds")
# Apply one-hot encoding to factor variables with >2 levels
covars_train_demographics <- dummy_cols(
  covars_train_demographics, 
  select_columns = c("ethnicity", "education", "household_income"), 
  remove_selected_columns = TRUE,  # Remove original categorical columns
  remove_first_dummy = FALSE       # Keep all dummy variables for completeness
)
# Rename "sex" to "male" and convert values
covars_train_demographics$male <- ifelse(covars_train_demographics$sex == "Male", 1, 0)
covars_train_demographics$sex <- NULL






covars_train_clinical <- readRDS("covars_train_clinical.rds")
covars_train_clinical <- covars_train_clinical[, -36] # remove CKDstage
categorical_vars <- c(
  "smoking_status", "smoking_household", "oily_fish", "processed_meat", "salt_added", 
  "alcohol_freq", "physactivity_type", "insomnia", "med_pain", "bodysize_age10", 
  "height_age10", "overall_health", "menopause", "num_live_births", 
  "preg_loss", "ever_hrt"
) # factor with more than 2 levels, for one-hot encoding
covars_train_clinical <- dummy_cols(
  covars_train_clinical, 
  select_columns = categorical_vars, 
  remove_first_dummy = FALSE,  # Keep all categories
  remove_selected_columns = TRUE  # Remove original categorical columns
)
# Convert yes/no to 1/0 
covars_train_clinical$multiple_birth <- ifelse(covars_train_clinical$multiple_birth == "Yes", 1, 0)
covars_train_clinical$maternal_smoking <- ifelse(covars_train_clinical$maternal_smoking == "Yes", 1, 0)
covars_train_clinical$snoring <- ifelse(covars_train_clinical$snoring == "Yes", 1, 0)






covars_train_comorbidities <- readRDS("covars_train_comorbidities.rds")
covars_train_comorbidities[, -1] <- as.data.frame(lapply(covars_train_comorbidities[, -1], function(col) {
  if (is.factor(col)) {
    return(as.numeric(col) - 1)  # Convert factor levels to numeric (0,1)
  } else if (is.logical(col)) {
    return(as.numeric(col))  # Convert logical (TRUE/FALSE) to numeric (1/0)
  } else {
    return(col)  # Keep numeric columns unchanged
  }
}))



covars_train_labs <- readRDS("covars_train_labs.rds")



################# Demographics sub-category 
# Logistic regression
covars_train_demographics$cluster_bin <- ifelse(covars_train_demographics$cluster == 0, 1, 0)
logistic <- glm(cluster_bin ~ ., data=covars_train_demographics[, -1], family="binomial")
summary(logistic)

# Lasso
X <- as.matrix(covars_train_demographics[, !(colnames(covars_train_demographics) %in% c("cluster", "cluster_bin"))])
Y <- covars_train_demographics$cluster_bin

# Cross-validation
set.seed(1)
cv_model <- cv.glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model)

# Visualize the regularisation path
cv_model.path <- glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model.path)

# Use lambda.1se to avoid overfitting in test set
best_lambda <- cv_model$lambda.1se

# Fit the logistic model with Lasso, using the best lambda
lasso_model <- glmnet(X, Y, family="binomial", alpha=1, lambda=best_lambda)

# Selected variables by Lasso
beta_lasso <- coef(lasso_model, s = "best_lambda")[2:(ncol(X) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(selected_lasso)

# Stability selection
out <- VariableSelection(xdata = X, ydata = Y, verbose = FALSE, family="binomial")
CalibrationPlot(out)
selprop <- SelectionProportions(out)
selprop = sort(selprop, decreasing = TRUE)
print(selprop)

# Calibrated parameters
hat_params <- Argmax(out)
print(hat_params)

# Visualisation of selection proportions
par(mar = c(12, 5, 2, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "",
     ylab = "Selection Proportion", xaxt = "n", bty = "l",
     main = "Stability-Lasso selected demographics variables for Cluster 0 vs 1/2",
     col = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"),
     cex.lab = 1.5, ylim = c(0, 1.2))

abline(h = hat_params[2], lty = 2, col = "darkred")

for (i in 1:length(selprop)) {
  axis(side = 1, at = i, labels = names(selprop)[i],
       las = 2,
       col = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       col.axis = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       cex.axis = 0.7)
  if (selprop[i] >= hat_params[2]) {
    text(i, selprop[i] + 0.05, labels = round(selprop[i], 2),
         col = "red", cex = 0.6, pos = 3, offset = 0.5)
  }
}

# Selected variables by stability analysis
stability_demographics <- names(selprop[selprop >= hat_params[2]])

# Dataframe to store selprop
selection_results <- as.data.frame(t(selprop))
selection_results <- cbind(cluster = 0, selection_results)


################# Clinical sub-category 
# Remove eGFR
covars_train_clinical <- covars_train_clinical[, !colnames(covars_train_clinical) %in% "eGFR"]
# Logistic regression
covars_train_clinical$cluster_bin <- ifelse(covars_train_clinical$cluster == 0, 1, 0)
logistic <- glm(cluster_bin ~ ., data=covars_train_clinical[, -1], family="binomial")
summary(logistic)

# Lasso
X <- as.matrix(covars_train_clinical[, !(colnames(covars_train_clinical) %in% c("cluster", "cluster_bin"))])
Y <- covars_train_clinical$cluster_bin

# Cross-validation
set.seed(1)
cv_model <- cv.glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model)

# Visualize the regularisation path
cv_model.path <- glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model.path)

# Use lambda.1se to avoid overfitting in test set
best_lambda <- cv_model$lambda.1se

# Fit the logistic model with Lasso, using the best lambda
lasso_model <- glmnet(X, Y, family="binomial", alpha=1, lambda=best_lambda)

# Selected variables by Lasso
beta_lasso <- coef(lasso_model, s = "best_lambda")[2:(ncol(X) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(selected_lasso)

# Stability selection
out <- VariableSelection(xdata = X, ydata = Y, verbose = FALSE, family="binomial")
CalibrationPlot(out)
selprop <- SelectionProportions(out)
selprop = sort(selprop, decreasing = TRUE)
print(selprop)

# Calibrated parameters
hat_params <- Argmax(out)
print(hat_params)

# Visualisation of selection proportions
par(mar = c(12, 5, 2, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "",
     ylab = "Selection Proportion", xaxt = "n", bty = "l",
     main = "Stability-Lasso selected clinical variables for Cluster 0 vs 1/2",
     col = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"),
     cex.lab = 1, cex.main = 0.9, ylim = c(0, 1.2))

abline(h = hat_params[2], lty = 2, col = "darkred")

for (i in 1:length(selprop)) {
  axis(side = 1, at = i, labels = names(selprop)[i],
       las = 2,
       col = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       col.axis = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       cex.axis = 0.7)
  if (selprop[i] >= hat_params[2]) {
    text(i, selprop[i] + 0.05, labels = round(selprop[i], 2),
         col = "red", cex = 0.6, pos = 3, offset = 0.5)
  }
}

# Selected variables by stability analysis
stability_clinical <- names(selprop[selprop >= hat_params[2]])

# Dataframe to store selprop
selection_results_clinical <- as.data.frame(t(selprop))
selection_results <- cbind(selection_results, selection_results_clinical)
selection_results$eGFR <- 0.95 # add eGFR which was removed


################# Comorbidities sub-category 
# Logistic regression
covars_train_comorbidities$cluster_bin <- ifelse(covars_train_comorbidities$cluster == 0, 1, 0)
logistic <- glm(cluster_bin ~ ., data=covars_train_comorbidities[, -1], family="binomial")
summary(logistic)

# Lasso
X <- as.matrix(covars_train_comorbidities[, !(colnames(covars_train_comorbidities) %in% c("cluster", "cluster_bin"))])
Y <- covars_train_comorbidities$cluster_bin

# Cross-validation
set.seed(1)
cv_model <- cv.glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model)

# Visualize the regularisation path
cv_model.path <- glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model.path)

# Use lambda.1se to avoid overfitting in test set
best_lambda <- cv_model$lambda.1se

# Fit the logistic model with Lasso, using the best lambda
lasso_model <- glmnet(X, Y, family="binomial", alpha=1, lambda=best_lambda)

# Selected variables by Lasso
beta_lasso <- coef(lasso_model, s = "best_lambda")[2:(ncol(X) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(selected_lasso)

# Stability selection
out <- VariableSelection(xdata = X, ydata = Y, verbose = FALSE, family="binomial")
CalibrationPlot(out)
selprop <- SelectionProportions(out)
selprop = sort(selprop, decreasing = TRUE)
print(selprop)

# Calibrated parameters
hat_params <- Argmax(out)
print(hat_params)

# Visualisation of selection proportions
par(mar = c(12, 5, 2, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "",
     ylab = "Selection Proportion", xaxt = "n", bty = "l",
     main = "Stability-Lasso selected comorbidities variables for Cluster 0 vs 1/2",
     col = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"),
     cex.lab = 1, cex.main = 0.9, ylim = c(0, 1.2))

abline(h = hat_params[2], lty = 2, col = "darkred")

for (i in 1:length(selprop)) {
  axis(side = 1, at = i, labels = names(selprop)[i],
       las = 2,
       col = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       col.axis = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       cex.axis = 0.7)
  if (selprop[i] >= hat_params[2]) {
    text(i, selprop[i] + 0.05, labels = round(selprop[i], 2),
         col = "red", cex = 0.6, pos = 3, offset = 0.5)
  }
}

# Selected variables by stability analysis
stability_comorbidities <- names(selprop[selprop >= hat_params[2]])

# Dataframe to store selprop
selection_results_comorbidities <- as.data.frame(t(selprop))
selection_results <- cbind(selection_results, selection_results_comorbidities)


################# Labs sub-category 
# Remove cystatin_C
covars_train_labs <- covars_train_labs[, !colnames(covars_train_labs) %in% "cystatin_C"]
# Logistic regression
covars_train_labs$cluster_bin <- ifelse(covars_train_labs$cluster == 0, 1, 0)
logistic <- glm(cluster_bin ~ ., data=covars_train_labs[, -1], family="binomial")
summary(logistic)

# Lasso
X <- as.matrix(covars_train_labs[, !(colnames(covars_train_labs) %in% c("cluster", "cluster_bin"))])
Y <- covars_train_labs$cluster_bin

# Cross-validation
set.seed(1)
cv_model <- cv.glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model)

# Visualize the regularisation path
cv_model.path <- glmnet(X, Y, family="binomial", alpha=1)
plot(cv_model.path)

# Use lambda.1se to avoid overfitting in test set
best_lambda <- cv_model$lambda.1se

# Fit the logistic model with Lasso, using the best lambda
lasso_model <- glmnet(X, Y, family="binomial", alpha=1, lambda=best_lambda)

# Selected variables by Lasso
beta_lasso <- coef(lasso_model, s = "best_lambda")[2:(ncol(X) + 1), ]
selected_lasso <- names(beta_lasso)[which(beta_lasso != 0)]
print(selected_lasso)

# Stability selection
out <- VariableSelection(xdata = X, ydata = Y, verbose = FALSE, family="binomial")
CalibrationPlot(out)
selprop <- SelectionProportions(out)
selprop = sort(selprop, decreasing = TRUE)
print(selprop)

# Calibrated parameters
hat_params <- Argmax(out)
print(hat_params)

# Visualisation of selection proportions
par(mar = c(12, 5, 2, 1))
plot(selprop, type = "h", lwd = 3, las = 1, xlab = "",
     ylab = "Selection Proportion", xaxt = "n", bty = "l",
     main = "Stability-Lasso selected lab variables for Cluster 0 vs 1/2",
     col = ifelse(selprop >= hat_params[2], yes = "red", no = "grey"),
     cex.lab = 1, cex.main = 0.9, ylim = c(0, 1.2))

abline(h = hat_params[2], lty = 2, col = "darkred")

for (i in 1:length(selprop)) {
  axis(side = 1, at = i, labels = names(selprop)[i],
       las = 2,
       col = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       col.axis = ifelse(selprop[i] >= hat_params[2], yes = "red", no = "grey"),
       cex.axis = 0.7)
  if (selprop[i] >= hat_params[2]) {
    text(i, selprop[i] + 0.05, labels = round(selprop[i], 2),
         col = "red", cex = 0.6, pos = 3, offset = 0.5)
  }
}

# Selected variables by stability analysis
stability_lab <- names(selprop[selprop >= hat_params[2]])


# Dataframe to store selprop
selection_results_labs <- as.data.frame(t(selprop))
selection_results <- cbind(selection_results, selection_results_labs)
selection_results$cystatin_C <- 1 # add cystatinC which was removed

# Save as rds file
saveRDS(selection_results, file = "cluster0_selection_results.rds")

### Load refit and test datasets
covars_refit <- readRDS("covars_refit.rds")
covars_test <- readRDS("covars_test.rds")

covars_refit$cluster_bin <- ifelse(covars_refit$cluster == 0, 1, 0)
covars_refit <- covars_refit[, !colnames(covars_refit) %in% "CKDstage_recruit"] # remove CKDstage
covars_test$cluster_bin <- ifelse(covars_test$cluster == 0, 1, 0)
covars_test <- covars_test[, !colnames(covars_test) %in% "CKDstage_recruit"] # remove CKDstage

# Apply one-hot encoding to factor variables with >2 levels
categorical_vars <- c("ethnicity", "education", "household_income", 
                      "smoking_status", "smoking_household", "oily_fish", "processed_meat", "salt_added", 
                      "alcohol_freq", "physactivity_type", "insomnia", "med_pain", "bodysize_age10", 
                      "height_age10", "overall_health", "menopause", "num_live_births", 
                      "preg_loss", "ever_hrt"
)
covars_refit <- dummy_cols(
  covars_refit, 
  select_columns = categorical_vars, 
  remove_first_dummy = FALSE,  
  remove_selected_columns = TRUE  
)
covars_test <- dummy_cols(
  covars_test, 
  select_columns = categorical_vars, 
  remove_first_dummy = FALSE,  
  remove_selected_columns = TRUE  
)

covars_refit$male <- ifelse(covars_refit$sex == "Male", 1, 0)
covars_refit$sex <- NULL

covars_test$male <- ifelse(covars_test$sex == "Male", 1, 0)
covars_test$sex <- NULL


# Fit logistic model using covars_refit, using the stability-Lasso selected variables
selected_variables <- c(stability_demographics,
                        stability_clinical,
                        stability_comorbidities,
                        stability_lab)
selected_variables <- paste0("`", selected_variables, "`")
formula_str <- paste("cluster_bin ~", paste(selected_variables, collapse = " + "))
formula_obj <- as.formula(formula_str)
logistic_refit <- glm(formula_obj, data=covars_refit[, -1], family="binomial")
summary(logistic_refit)

summary_obj <- summary(logistic_refit)
results_df <- data.frame(
  Variable = rownames(summary_obj$coefficients),
  Coefficient = round(summary_obj$coefficients[, "Estimate"], 4),  
  P_Value = round(summary_obj$coefficients[, "Pr(>|z|)"], 4),  
  OR = round(exp(summary_obj$coefficients[, "Estimate"]), 2)  
)
significant_vars <- results_df[results_df$P_Value < 0.05, ]
print(significant_vars)

# Find optimal threshold (default is Youden's index)
lr.pred.probs = predict(logistic_refit, type = "response")
covars_refit <- cbind(covars_refit, probs = lr.pred.probs)
roc.obj <- roc(covars_refit$cluster_bin, covars_refit$probs)
optimal_threshold <- coords(roc.obj, "best", ret = "threshold")
optimal_threshold <- optimal_threshold[1, "threshold"]

# Plot the ROC curve
plot(roc.obj, col = "blue", main = "Cluster 0 ROC")
auc_value <- auc(roc.obj)
text(0.2, 0.05, labels = paste0("AUC = ", round(auc_value, 3)),
     col = "blue", cex = 0.7)


# Fit logistic regression on test set
lr.pred.probs.test <- predict(logistic_refit, type = "response", newdata = covars_test)
covars_test <- cbind(covars_test, pred.probs = lr.pred.probs.test)
covars_test$pred.class <- ifelse(covars_test$pred.probs >= optimal_threshold, 1, 0)
cf.mat <- table(covars_test$cluster_bin, covars_test$pred.class)
cf.mat

CA <- sum(diag(cf.mat)/length(covars_test$cluster_bin))
CA

TPR <- cf.mat[4]/sum(covars_test$cluster_bin == 1)
TPR

FPR <- cf.mat[3]/sum(covars_test$cluster_bin == 0)
FPR









