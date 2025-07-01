rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
final_dataset_dir <- args[1]
split_dir <- args[2]
output_dir <- args[3]

# currently using sharp v1.4.6
if (!require("sharp")) install.packages(c("sharp"), repos = "https://cran.rstudio.com/", dependencies=TRUE)

source("sharp_lasso_plots.R")
setwd(output_dir)

library(dplyr)
library(sharp)
library(caret)
library(pROC)

### Load train/test datasets
# train_select refers to the 50% of training data used to select importnat proteins; 
# train_refit is the next 30% of training data used to finalise the model; 
# test dataset is the 20% used to validate model
train_select <- readRDS(file.path(split_dir, "train_select.rds"))
train_refit <- readRDS(file.path(split_dir, "train_refit.rds"))
test <- readRDS(file.path(split_dir, "test.rds"))
proteins <- readRDS(file.path(final_dataset_dir, "proteins_scaled_10years.rds"))

print("Train_select subsample cluster proportions:")
print(table(train_select$cluster))
print("Train_refit subsample cluster proportions:")
print(table(train_refit$cluster))
print("Test subsample cluster proportions:")
print(table(test$cluster))

### Analyse protein associations for each cluster using LASSO and sharp stability selection 
for (clusterID in unique(test$cluster)) {
  # Extract cluster labels from split dataset 
  # convert to binary outcome for 1-vs-all analysis
  Y_train <- as.factor(ifelse(train_select$cluster == clusterID, 1, 0))
  Y_refit <- as.factor(ifelse(train_refit$cluster == clusterID, 1, 0))
  Y_test <- as.factor(ifelse(test$cluster == clusterID, 1, 0))
  
  rownames(proteins) <- as.character(rownames(proteins))  # Ensure row names are character
  train_select$ID <- as.character(train_select$ID)  # Ensure IDs are character
  train_refit$ID <- as.character(train_refit$ID)  # Ensure IDs are character
  test$ID <- as.character(test$ID)  # Ensure IDs are character
  
  # Remove ID and cluster columns from input features
  X_train <- as.matrix(proteins[train_select$ID, ])
  X_refit <- as.matrix(proteins[train_refit$ID, ])
  X_test <- as.matrix(proteins[test$ID, ])
  
  ## Lasso variable selection using sharp package on train_select
  t0 <- Sys.time()
  lasso_selection <- VariableSelection(
    xdata = X_train, ydata = Y_train, 
    verbose = FALSE,
    penalty.factor = rep(1, ncol(X_train)), # apply penalty to all columns
    family = "binomial", # binary logistic regression 
    n_cat = 3, # stable inclusion, unstable, stable exclusion
    pi_list = seq(0.5, 0.95, by = 0.01), # between 0.5 and 1 if n_cat = 3
    K = 100, seed = 1,
  )
  t1 <- Sys.time()
  print(paste("Cluster", clusterID, "- Time taken for sharp stability selection:", t1 - t0))
  
  # Visualise using calibration plot in sharp
  save_calibration_plot(lasso_selection, clusterID, output_dir)
  
  # Calibrated selection proportions
  selprop <- SelectionProportions(lasso_selection)
  
  # Calibrated parameters
  hat_params <- Argmax(lasso_selection)
  print("Sharp-calibrated parameters:")
  print(hat_params)
    
  print("Number of selected proteins:")
  print(table(selprop >= hat_params[2]))
  
  # Visualisation of selection proportions
  save_selprop_plot(selprop, hat_params, clusterID, output_dir)
  save_selprop_plot(selprop, hat_params, clusterID, output_dir, top = 40)

  ## Use train_refit to assess strength of association for the stably selected proteins in lasso selection
  selprop_stable <- selprop[selprop >= hat_params[2]]
  protein_stable <- names(selprop_stable)

  # save stably selected proteins to file
  save_stable_proteins(selprop, hat_params, clusterID, output_dir) 
  
  # logistic regression for cluster membership vs stable proteins
  X_refit_cluster <- as.matrix(X_refit[,protein_stable])
  Y_refit_cluster <- as.numeric(ifelse(train_refit$cluster == 0, 1, 0))
  protein_assoc_cluster <- glm(Y_refit_cluster ~.,family=binomial(link = 'logit'), data = as.data.frame(cbind(X_refit_cluster, Y_refit_cluster)))


  # plot regression coefficients
  save_betas_sorted_plot(protein_assoc_cluster, clusterID, output_dir)

  ### Evaluate predictive performance of model using test dataset
  predicted_test <- predict(protein_assoc_cluster, newdata = as.data.frame(X_test), type = "response")
  conf_mat <- confusionMatrix(data =as.factor(ifelse(predicted_test>0.5, 1, 0)), reference = Y_test)
  
  # save ROC curve
  roc_cluster <- roc(Y_test~predicted_test, auc = TRUE)
  save_roc_curve(
    roc_curve = roc_cluster, 
    auc = auc(roc_cluster), 
    clusterID, 
    output_dir
  ) 

  # save model
  save_stable_lasso_model(protein_assoc_cluster, conf_mat, clusterID, output_dir)
}

print("DONE")
