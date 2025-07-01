rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
cluster_labels_path <- args[1]
split_output_dir <- args[2]

setwd(split_output_dir)
library(tidyverse)
library(caret)

df <- read_csv(cluster_labels_path)

# Cluster membership
if ("ID" %in% colnames(df)) {
  df_new <- cbind(df$ID, df$Cluster)
} else if ("index" %in% colnames(df)) {
  df_new <- cbind(df$index, df$Cluster)
} else {
  df_new <- cbind(df$`Unnamed: 0`, df$Cluster)
}

df_new <- as.data.frame(df_new)
colnames(df_new) <- c("ID", "cluster")
dim(df_new) # 2016 by 2 (6 participants removed)

# Stratified split 50/30/20, preserving cluster proportions across subsets
set.seed(100)

train_select_idx <- createDataPartition(df_new$cluster, p = 0.5, list = FALSE)
train_select <- df_new[train_select_idx, ]

remaining_data <- df_new[-train_select_idx, ]
train_refit_idx <- createDataPartition(remaining_data$cluster, p = 0.6, list = FALSE)
train_refit <- remaining_data[train_refit_idx, ]
test <- remaining_data[-train_refit_idx, ]

# Check the split
cat("Original proportion:\n")
print(prop.table(table(df_new$cluster)))

cat("\nTrain Select proportion (50%):\n")
print(prop.table(table(train_select$cluster)))

cat("\nTrain Refit proportion (30%):\n")
print(prop.table(table(train_refit$cluster)))

cat("\nTest proportion (20%):\n")
print(prop.table(table(test$cluster)))

saveRDS(train_select, file = "train_select.rds")
saveRDS(train_refit, file = "train_refit.rds")
saveRDS(test, file = "test.rds")
