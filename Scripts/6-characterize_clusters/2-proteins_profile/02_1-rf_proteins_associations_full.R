rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
final_dataset_dir <- args[1]
split_dir <- args[2]
output_dir <- args[3]

# ---------------------------------------------------------
# Set working directory
setwd(output_dir)
# ---------------------------------------------------------

# Load required packages
library(randomForest)
library(pROC)
library(caret)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot) 
library(extrafont)
library(ggtext)
library(gtable)
library(patchwork)
library(tidyverse)


# Load scaled proteins, training and test datasets
train_select  <- readRDS(file.path(split_dir, "train_select.rds"))
train_refit  <- readRDS(file.path(split_dir, "train_refit.rds"))
test       <- readRDS(file.path(split_dir, "test.rds"))
proteins   <- readRDS(file.path(final_dataset_dir, "proteins_scaled_10years.rds"))

# Combine train_select and train_refit and ensure protein names are character
train <- rbind(train_select, train_refit)
train$ID <- as.character(train$ID)
test$ID  <- as.character(test$ID)
rownames(proteins) <- as.character(rownames(proteins))





############ Random Forest model tuning and training ############

# Default RF hyperparameter values
default_parameters <- list( mtry = floor(sqrt(ncol(proteins))), 
                            ntree = 500,
                            nodesize = 1,
                            maxnodes = 500)

# Loop across cluster IDs
for (clust_id in sort(unique(test$cluster))) {
  
  # Binary classification labels
  Y_train <- factor(ifelse(train$cluster == clust_id, "Cluster", "Non_Cluster"),
                    levels = c("Non_Cluster", "Cluster"))
  
  # Feature matrix
  X_train <- as.matrix(proteins[train$ID, ])
  
  # Train RF model with default hyperparameter values
  rf_default <- randomForest(x = X_train,
                             y = Y_train,
                             mtry     = default_parameters$mtry,
                             ntree    = default_parameters$ntree,
                             nodesize = default_parameters$nodesize,
                             maxnodes = default_parameters$maxnodes,
                             importance = TRUE)
  
  # Save default RF model
  saveRDS(rf_default,
          file = paste0("default_rf_cluster_", clust_id, ".rds"))

}


# RF Model Tuning
# Define hyperparameter grid search values
mtry_values <- c(floor(sqrt(ncol(proteins))),
                  floor(ncol(proteins)/3),
                  floor(ncol(proteins)/2))
ntree_values <- c(500, 1000)
nodesize_values <- c(1, 5)
maxnodes_values <- c(100, 500)


# Set up tuning for RF models independently for clusters based on optimal hyperparameter values
tuned_cluster <- function(clust_id) {
  # Binary classifications 
  Y_train <- factor(ifelse(train$cluster == clust_id, "C", "N"), levels = c("N","C"))
  Y_test <- factor(ifelse(test$cluster  == clust_id, "C", "N"), levels = c("N","C"))
  
  # Feature matrices
  X_train <- as.matrix(proteins[ train$ID, ])
  X_test <- as.matrix(proteins[ test$ID, ])
  
  optimal_auc   <- 0
  optimal_combo <- NULL
  
  for (m in mtry_values)
    for (nt in ntree_values)
      for (ns in nodesize_values)
        for (mx in maxnodes_values) {
          rf <- randomForest(x = X_train, y = Y_train,
                             mtry = m,
                             ntree = nt,
                             nodesize = ns,
                             maxnodes = mx,
                             importance = FALSE)
          p <- predict(rf, X_test, type = "prob")[,2]
          auc_values <- auc(roc(Y_test, p))
          if (auc_values >  optimal_auc) {
            optimal_auc   <- auc_values
            optimal_combo <- list(mtry=m, ntree=nt, nodesize=ns, maxnodes=mx)
          }
        }
  
  list(cluster = clust_id,
       cat(sprintf("Cluster %s | mtry=%d, ntree=%d, nodesize=%d, maxnodes=%d | AUC=%.4f\n",
                   clust_id, m, nt, ns, mx, auc_values)),
       optimal_params = optimal_combo,
       optimal_auc    =  optimal_auc)
}


# Tune each cluster and save
results <- lapply(sort(unique(train$cluster)), tuned_cluster)

# Save the optimal hyperparameter values as list
optimal_params_list <- setNames( lapply(results, `[[`, " optimal_params"),
  paste0("Cluster_", sapply(results, `[[`, "cluster"))
)
saveRDS(optimal_params_list, "optimal_hyperparameters.rds")

# Print summary of optimal hyperparameter values per cluster
for (r in results) {
  cat(sprintf("Cluster %s → AUC=%.3f → mtry=%d, ntree=%d, nodesize=%d, maxnodes=%d\n",
              r$cluster, r$optimal_auc,
              r$optimal_params$mtry, r$optimal_params$ntree,
              r$optimal_params$nodesize, r$optimal_params$maxnodes))
}


# Define optimal hyperparameter values per cluster
optimal_parameters <- list(
  "0" = list(mtry = sqrt(ncol(proteins)), ntree = 1000,  nodesize = 1, maxnodes = 500),
  "1" = list(mtry = sqrt(ncol(proteins)), ntree = 500, nodesize = 1, maxnodes = 100),
  "2" = list(mtry = sqrt(ncol(proteins)), ntree = 500,  nodesize = 1, maxnodes = 500))

# Train and save RF models per cluster
for (num in names(optimal_parameters)) {
  parameters <- optimal_parameters[[num]]
  
  # Binary classifications: "Cluster" vs. "Non_Cluster"
  Y_train <- factor(ifelse(train$cluster == as.integer(num), "Cluster", "Non_Cluster"),
                    levels = c("Non_Cluster", "Cluster"))
  
  # Feature matrix
  X_train <- as.matrix(proteins[train$ID, ])
  
  # Train RF using optimal hyperparameter values
  rf_model <- randomForest(x = X_train, y = Y_train,
                           mtry = parameters$mtry,
                           ntree = parameters$ntree,
                           nodesize = parameters$nodesize,
                           maxnodes = parameters$maxnodes,
                           importance = TRUE)
  
  # Save tuned RF model for each cluster
  saveRDS(rf_model, paste0("tuned_rf_cluster_", num, ".rds"))
}









############ Plotting barplots to compare top 20 predictive proteins per cluster based on MDA vs Gini ############

# Define cluster configurations
clusters <- c("0", "1", "2")
importance_types <- c("MDA", "Gini")
model_colors <- c("Default RF" = "cornsilk3", "Tuned RF" = "firebrick3")
cluster_labels <- c("0" = "Cluster Blue", "1" = "Cluster Orange", "2" = "Cluster Green")
panel_labels <- LETTERS[2:7]

# Extract top cluster-predictive proteins using MDA and Gini-Impurity
extract_top_importance <- function(rf_model, type = "MeanDecreaseAccuracy", top_n = 20) {
  column <- ifelse(type == "MDA", "MeanDecreaseAccuracy", "MeanDecreaseGini")
  imp <- importance(rf_model)
  imp_df <- data.frame(Protein = rownames(imp), Importance = imp[, column])
  imp_df <- imp_df[order(-imp_df$Importance), ][1:top_n, ]
  return(imp_df)
}

# Tabulating hyperparameter values table
param_table <- data.frame(
  Cluster = c("Cluster Blue", "Cluster Orange", "Cluster Green"),
  AUC = c(0.987, 0.98, 0.984),
  mtry = c(36, 36, 36),
  ntree = c(1000, 500, 500),
  nodesize = c(1, 1, 1),
  maxnodes = c(500, 100, 500))

# Format table
table_grob <- tableGrob(param_table, rows = NULL, theme = ttheme_default(base_size = 12))
title <- textGrob("A. Optimal Hyperparameter Values of Tuned RF Models",
                  gp = gpar(fontface = "bold", fontsize = 14, fontfamily = "Arial"), just = "center")

table_with_title <- gtable_add_rows(table_grob, heights = grobHeight(title) + unit(5, "mm"), pos = 0)
table_with_title <- gtable_add_grob(table_with_title, title, 1, 1, 1, ncol(table_with_title))

# Set up panels for plotting protein rankings side by side for MDA and Gini per cluster
plot_list <- list()
label_idx <- 1

for (cluster in clusters) {
  rf_default <- readRDS(paste0("default_rf_cluster_", cluster, ".rds"))
  rf_tuned   <- readRDS(paste0("tuned_rf_cluster_", cluster, ".rds"))
  cl_label <- cluster_labels[[cluster]]
  
  for (imp_type in c("MDA", "Gini")) {
    df_default <- extract_top_importance(rf_default, imp_type) %>% mutate(Model = "Default RF")
    df_tuned   <- extract_top_importance(rf_tuned, imp_type) %>% mutate(Model = "Tuned RF")
    combined <- bind_rows(df_default, df_tuned)
    
    p <- ggplot(combined, aes(x = reorder(Protein, Importance), y = Importance, fill = Model)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      scale_fill_manual(values = model_colors) +
      labs(
        title = paste0(cl_label, " - Top 20 Proteins (", imp_type, ")"),
        x = "Protein", y = "Importance") +
      theme_minimal(base_family = "Arial") +
      theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11))
    
    p <- ggdraw(p) +
      draw_label(panel_labels[label_idx], x = 0.02, y = 0.98, hjust = 0, vjust = 1,
                 fontface = "bold", size = 14, fontfamily = "Arial")
    
    plot_list[[paste0(cluster, "_", imp_type)]] <- p
    label_idx <- label_idx + 1
  }
}

# Save final combined panels
final_panel <- plot_grid(
  table_with_title,
  plot_grid(plot_list[["0_MDA"]], plot_list[["0_Gini"]], ncol = 2),
  plot_grid(plot_list[["1_MDA"]], plot_list[["1_Gini"]], ncol = 2),
  plot_grid(plot_list[["2_MDA"]], plot_list[["2_Gini"]], ncol = 2),
  ncol = 1,
  rel_heights = c(0.55, 1.1, 1.1, 1.1)
)

ggsave("rf_importance_combined_panel.png", plot = final_panel, width = 16, height = 26, dpi = 300)










############ Computing and plotting combined ROC Curves and confusion matrices ############

# Load trained RF models
rf_models <- list(
  "Cluster 0" = readRDS("rf_tuned_cluster_0.rds"),
  "Cluster 1" = readRDS("rf_tuned_cluster_1.rds"),
  "Cluster 2" = readRDS("rf_tuned_cluster_2.rds")
)

# Define thresholds computed using Youden's J
thresholds <- c("Cluster 0" = 0.385, "Cluster 1" = 0.4645, "Cluster 2" = 0.457)

# Manual color map for clusters
cluster_colors <- c(
  "Cluster 0" = "#1f77b4",  
  "Cluster 1" = "#ff7f0e",  
  "Cluster 2" = "#2ca02c" 
  )

# Change cluster names to colours
plot_labels <- c(
  "Cluster 0" = "Cluster Blue",
  "Cluster 1" = "Cluster Orange",
  "Cluster 2" = "Cluster Green"
)

# Match IDs
matched_ids <- intersect(test$ID, rownames(proteins))
test_filtered <- test[test$ID %in% matched_ids, ]
proteins_filtered <- proteins[matched_ids, ]
factor_levels <- c("Non-cluster", "Cluster")

# Plot ROC and confusion matrix panel
roc_plots <- list()
conf_plots <- list()

# Loop through all clusters
for (clusterID in 0:2) {
  cname <- paste0("Cluster ", clusterID)
  colour_name <- plot_labels[[cname]]
  color <- cluster_colors[[cname]]
  
  # Labels and input
  test_labels <- factor(ifelse(test_filtered$cluster == clusterID, "Cluster", "Non-cluster"), levels = factor_levels)
  X_test <- as.matrix(proteins_filtered)
  
  # Predict probability of belonging to cluster vs not
  probs <- predict(rf_models[[cname]], X_test, type = "prob")[, 2]
  pred <- factor(ifelse(probs > thresholds[[cname]], "Cluster", "Non-cluster"), levels = factor_levels)
  
  # Define performace metrics from ROC data
  roc_curve <- roc(test_labels, probs)
  auc_val <- round(auc(roc_curve), 3)
  roc_df <- data.frame(
    Sensitivity = rev(roc_curve$sensitivities),
    Specificity = rev(roc_curve$specificities),
    Cluster = cname
  )
  roc_df$AUC <- auc_val
  roc_plots[[cname]] <- roc_df
  
  # Plot confusion matrix with actual values and corresponding %
  cm <- confusionMatrix(pred, test_labels)
  cm_table <- as.data.frame(cm$table)
  cm_table$FreqRaw <- cm_table$Freq
  cm_table$Freq <- round(cm_table$Freq / sum(cm_table$Freq) * 100, 1)
  
  conf_plot <- ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(FreqRaw, " (", Freq, "%)")), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = color) +
    theme_minimal(base_size = 16) +
    ggtitle(paste0(colour_name, " Confusion Matrix")) +
    theme(
      plot.title = element_text(color = color, hjust = 0.5, size = 20, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold")
    )
  
  
  conf_plots[[cname]] <- conf_plot
}

# Combine ROC curves and AUCs into one plot
roc_all_df <- bind_rows(roc_plots)

# AUC values legend
auc_labels <- sapply(names(roc_plots), function(cname) {
  colour_name <- plot_labels[[cname]]
  auc_val <- unique(roc_plots[[cname]]$AUC)
  paste0(colour_name, " (AUC = ", auc_val, ")")
})
names(auc_labels) <- names(roc_plots)

# Final combined ROC Plot
ggplot(roc_all_df, aes(x = 1 - Specificity, y = Sensitivity, color = Cluster)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = cluster_colors, labels = auc_labels) +
  labs(
    title = "Evaluation of Cluster-specific RF Predictive Accuracy",
    x = "1 - Specificity", y = "Sensitivity", color = NULL
  ) +
  theme_minimal() +
  theme( plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.position = c(0.99, 0.01),
    legend.justification = c("right", "bottom"))

# Combine plots to panel and format
conf_plots[["Cluster 0"]] + conf_plots[["Cluster 1"]] + conf_plots[["Cluster 2"]] +
  plot_layout(nrow = 1) +
  plot_annotation(theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))

# Save the combined 1x3 confusion matrix panel as PDF
pdf("confusion_matrices_panel.pdf", width = 16, height = 6)
conf_plots[["Cluster 0"]] + 
  conf_plots[["Cluster 1"]] + 
  conf_plots[["Cluster 2"]] +
  plot_layout(nrow = 1) +
  plot_annotation(theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
dev.off()







############ Plot circular barplot for top 20 proteins per cluster and group by Olink functional panels ############

# Load importance rankings
top_rf_files <- list.files(output_dir, pattern = "top_rf_proteins_cluster_.*\\.rds", full.names = TRUE)
top_rf_list <- lapply(top_rf_files, readRDS)
names(top_rf_list) <- c("Cluster 0", "Cluster 1", "Cluster 2")
top_rf_df <- bind_rows(top_rf_list, .id = "Cluster")

# Load Olink functional panel annotations
panel_annotations <- readRDS(file.path(final_dataset_dir, "../Panel_annotation.rds"))
annotation_df <- data.frame(Protein = names(panel_annotations), FunctionalGroup = panel_annotations)

# Normalise and join annotations
top20 <- top_rf_df %>%
  group_by(Cluster) %>%
  slice_max(Importance, n = 20, with_ties = FALSE) %>%
  mutate(NormImportance = Importance / max(Importance)) %>%
  left_join(annotation_df, by = "Protein") %>%
  ungroup()

# Add empty spacing between bars of functional groups
empty_bar <- 2
clusters <- unique(top20$Cluster)
empty_df <- expand.grid(
  FunctionalGroup = unique(top20$FunctionalGroup),
  Cluster = clusters,
  Protein = NA,
  Importance = 0,
  NormImportance = 0
)
empty_df <- empty_df[rep(1:nrow(empty_df), each = empty_bar), ]
long_df <- bind_rows(top20, empty_df)

# Assign circular position and protein name per bar, including duplicates across clusters
long_df <- long_df %>%
  arrange(FunctionalGroup, Protein) %>%
  mutate(id = row_number())

# Label data with one label per bar
label_data <- long_df %>%
  filter(!is.na(Protein)) %>%
  mutate(angle = 90 - 360 * (id - 0.5) / nrow(long_df),
         hjust = ifelse(angle < -90, 1, 0),
         angle = ifelse(angle < -90, angle + 180, angle))


# Separator lines for functional group
base_data <- long_df %>%
  filter(!is.na(Protein)) %>%
  group_by(FunctionalGroup) %>%
  summarise(start = min(id), end = max(id), .groups = "drop") %>%
  mutate(title = (start + end) / 2)

# Define grid lines and cluster colours 
grid_lines <- tibble(y = c(0.25, 0.5, 0.75, 1.0))
cluster_colors <- c("Cluster 0" = "#377eb8", "Cluster 1" = "#ff7f00", "Cluster 2" = "#4daf4a")

# Plot circular barplot
p <- ggplot(long_df) +
  geom_bar(aes(x = as.factor(id), y = NormImportance, fill = Cluster),
           stat = "identity", width = 1, color = "black", linewidth = 0.1, na.rm = TRUE) +
  geom_hline(data = grid_lines, aes(yintercept = y),
             color = "grey70", linetype = "dotted", linewidth = 0.3) +
  geom_text(data = label_data,
            aes(x = id, y = NormImportance + 0.05, label = Protein, angle = angle, hjust = hjust),
            size = 3, fontface = "bold", inherit.aes = FALSE) +
  geom_segment(data = base_data,
               aes(x = start, xend = end, y = -0.05, yend = -0.05),
               color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(values = cluster_colors) +
  ylim(-0.4, 1.1) +
  coord_polar() +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  labs(title = "Circular Barplot of Top 20 RF Proteins per Cluster (MDA)",
       fill = "Cluster")

# Save ciruclar barplot as PDF
ggsave("circular_rf_proteins_mda.pdf", plot = p, width = 10, height = 10)


