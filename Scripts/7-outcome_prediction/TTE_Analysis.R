rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
imputed_covars_path <- args[1]
split_dir <- args[2]
output_dir <- args[3]

source("calculate_eGFR.R")
source("censored_dates.R")
setwd(output_dir)

library(dplyr)
library(ggplot2)
library(patchwork)
library(survival)
library(survminer)
library(caret)

# use binary cluster membership to cluster 1 (orange)
USE_BINARY_CLUSTER <- TRUE

### Load the imputed covariates data
covariates <- readRDS(imputed_covars_path)

# convert date types
covariates$esrd_date <- as.Date(covariates$esrd_date)
covariates$ckd_date <- as.Date(covariates$ckd_date)
covariates$death_date <- as.Date(covariates$death_date)

# convert serum creatinine from micromol/L to mg/dl and update column name
colnames(covariates)[colnames(covariates) == "creatinine"] <- "creatinine_mg.dl"
covariates$creatinine_mg.dl <- covariates$creatinine_mg.dl/88.42
# calculate eGFR
covariates <- covariates %>% mutate(eGFR = mapply(calculate_eGFR, creatinine_mg.dl, age_recruit, sex))

### Create dataframe for training
train_clusters <- rbind(
  readRDS(file.path(split_dir, "train_select.rds")),
  readRDS(file.path(split_dir, "train_refit.rds"))
)
train_clusters$ID <- as.character(train_clusters$ID)
rownames(train_clusters) <- train_clusters$ID

train_df <- cbind(train_clusters, covariates[train_clusters$ID,])

### Create dataframe for testing
test_clusters <- readRDS(file.path(split_dir, "test.rds"))
test_clusters$ID <- as.character(test_clusters$ID)
rownames(test_clusters) <- test_clusters$ID

test_df <- cbind(test_clusters, covariates[test_clusters$ID,])

get_surv_cols <- function(df, start_date, end_date, censored_date) {
  censored_date <- as.Date(censored_date, format = "%Y-%m-%d")
  
  df$event <- ifelse(!is.na(df[[end_date]]), 1, 0)
  df$stime <- as.Date(ifelse(!is.na(df[[end_date]]), df[[end_date]], censored_date)) - df[[start_date]]
  
  print(paste(">>", start_date, "to", end_date))
  print(paste("censored_date:", censored_date))
  print(paste("end_date max:", max(df[[end_date]], na.rm = TRUE)))
  print(paste("stime max:", max(df$stime)))
  print(paste("event #:", sum(df$event)))
  print(paste("final # at risk:", sum(df$stime == max(df$stime) & df$event == 0)))

  return(df)
}

cluster_km_curves <- function(data, start_date, end_date, censored_date, label, res_dir, cluster_bin=USE_BINARY_CLUSTER) {
  if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
  
  data <- get_surv_cols(data, start_date, end_date, censored_date)
  
  if (cluster_bin){
    surv_fit <- survfit(Surv(time = data$stime, event = data$event, type = "right") ~ cluster_bin, data = data) 
    colours <- c('#60594d', '#ff7f0e')
    labels <- c("Blue + Green", "Orange")
    filename <- paste0(end_date, "_", start_date, "_km_curves_binary.pdf")
  } else {
    surv_fit <- survfit(Surv(time = data$stime, event = data$event, type = "right") ~ cluster, data = data) 
    colours <- c('#1f77b4', '#ff7f0e', '#2ca02c')
    labels <- c("Blue", "Orange", "Green")
    filename <- paste0(end_date, "_", start_date, "_km_curves_cluster.pdf")
  }
  p <- ggsurvplot(surv_fit, data = data, pval = TRUE, conf.int = FALSE,
                  risk.table = TRUE,
                  ggtheme = theme_light(),
                  palette = colours,
                  # title = paste("Time to", label, "by Cluster"),
                  xlab = paste("Time (days) since CKD diagnosis"),
                  ylab = paste(label, "survival probability"),
                  legend = "none",
                  legend.title = " ",
                  legend.labs = labels,
                  title.position = "top",
                  title.hjust = 1)
  ggexport(filename = file.path(res_dir, filename), plot = p, device = "pdf")
  return(p)
}

save_cox_model <- function(model, res_dir, model_name, num_plots) {
  sink(file.path(res_dir, paste0(model_name, ".txt")))
  cat(paste("Cox regression", model_name))
  cat("\n\n========\n")
  # training concordance
  conc <- concordance(model)
  print(conc) 
  cat("\n\n========\n")
  # save model summary - HRs
  print(summary(model))
  cat("\n\n========\n")
  # HRs and 95% CIs
  cat("Hazard ratios and 95% CIs:\n")
  haz_ratios <- tryCatch(exp(coef(model)), error = function(e) {
    cat("Error getting HRs:", conditionMessage(e), "\n")
    NA
  })
  conf_ints <- tryCatch(exp(confint(model)), error = function(e) {
    cat("Error getting CIs:", conditionMessage(e), "\n")
    NA
  })
  tryCatch(print(cbind(haz_ratios, conf_ints)), error = function(e) {
    cat("Error combining HR and CIs:", conditionMessage(e), "\n")
    print(haz_ratios)
    print(conf_ints)
  })
  cat("\n\n========\n")
  cat("Proportional hazards assumptions test:\n")
  # check assumptions
  print(cox.zph(model))
  sink(NULL)

  # schoenfeld plots
  png(file.path(res_dir, paste0(model_name, "_scaled_schoenfeld_residuals.png")), width = 1000, height = 400)
  par(oma = c(1, 1, 5, 1), mfrow=c(1, num_plots), pty="s")
  plot(cox.zph(model, terms = TRUE))
  mtext(paste("Scaled Schoenfeld Residual plots for", model_name), side = 3, line = 1, outer = TRUE, cex = 2)
  dev.off()
}

evaluate_test_preds <- function(model, test_surv, test_data, cols) {
  # TODO bootstrap test and repeat 
  conc <- concordance(test_surv ~ predict(model, newdata = test_data[,cols], type = "risk"), data = test_data[,cols])
  return(conc)
}

compare_cox_reg <- function(train_data, test_data, start_date, end_date, censored_date, covars_list, res_dir, cluster_bin=USE_BINARY_CLUSTER) {
  if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
  
  train_data <- get_surv_cols(train_data, start_date, end_date, censored_date)
  test_data <- get_surv_cols(test_data, start_date, end_date, censored_date)
  
  # cluster_km_curves(train_data, paste0(end_date, "_"))

  # Model 0: demographic covariates
  train_surv <- Surv(time = train_data$stime, event = train_data$event, type = "right")
  model_0 <- coxph(train_surv ~ ., data = train_data[,covars_list])
  save_cox_model(model_0, res_dir, "model_0", length(covars_list))

  # Model 1: demographic covariates + clusters
  if (cluster_bin) {
    model_1 <- coxph(train_surv ~ ., data = train_data[,c(covars_list, "cluster_bin")])
  } else {
    model_1 <- coxph(train_surv ~ ., data = train_data[,c(covars_list, "cluster")])
  }
  save_cox_model(model_1, res_dir, "model_1", length(covars_list)+1)

  # Comparing Fits
  sink(file.path(res_dir, "test_comparison.txt"))
  cat("Likelihood ratio test:\n")
  print(anova(model_0, model_1, test = "LRT"))
  cat("\n\n========\n")

  # Evaluate prediction on test data
  test_surv <- Surv(time = test_data$stime, event = test_data$event, type = "right")

  concordance_model_0 <- evaluate_test_preds(model_0, test_surv, test_data, covars_list)
  if (cluster_bin) {
    test_data$cluster_bin <- as.factor(ifelse(test_data$cluster == 1, 1, 0))
    concordance_model_1 <- evaluate_test_preds(model_1, test_surv, test_data, c(covars_list, "cluster_bin"))
  } else {
    concordance_model_1 <- evaluate_test_preds(model_1, test_surv, test_data, c(covars_list, "cluster"))
  }

  cat("Concordance model 0:\n") 
  print(concordance_model_0)
  cat("\n========\n")
  cat("Concordance model 1:\n") 
  print(concordance_model_1)

  sink(NULL)
}

# transform variables for PH assumptions
### binary cluster membership
train_df$cluster_bin <- as.factor(ifelse(train_df$cluster == 1, 1, 0))
test_df$cluster_bin <- as.factor(ifelse(test_df$cluster == 1, 1, 0))
### log transform blood creatinine
train_df$log_creatinine <- log10(train_df$creatinine_mg.dl)
test_df$log_creatinine <- log10(test_df$creatinine_mg.dl)
### 5-year age bands
train_df$age_bin <- cut(train_df$age_recruit, breaks = c(40, 45, 50, 55, 60, 65, 70), right = FALSE)
test_df$age_bin <- cut(test_df$age_recruit, breaks = c(40, 45, 50, 55, 60, 65, 70), right = FALSE)

# plot KM curves
### time from diagnosis to ESRD date
km_esrd <- cluster_km_curves(
  rbind(train_df, test_df), 
  start_date = "ckd_date", 
  end_date = "esrd_date", 
  censored_date = HES_CENSORING_DATE, 
  label = "ESRD diagnosis",
  res_dir = "TTE_esrd",
  cluster_bin=FALSE
)
### time from diagnosis to death 
km_death <- cluster_km_curves(
  rbind(train_df, test_df), 
  start_date = "ckd_date", 
  end_date = "death_date", 
  censored_date = DEATH_CENSORING_DATE, 
  label = "Mortality",
  res_dir = "TTE_death",
  cluster_bin=FALSE
)
### combined plot 
merged_plot <- ((km_esrd$plot | km_death$plot) / (km_esrd$table | km_death$table)) + plot_layout(widths = c(6, 6), heights = c(6, 1.5), axes="collect", axis_titles="collect") + plot_annotation(tag_levels = 'A')
ggsave(filename = "km_merged.png", plot = merged_plot, device = "png", width = 12, height = 7.5)

# model 0 with log10(creatinine)
compare_cox_reg(
  train_df, test_df,
  start_date = "ckd_date", 
  end_date = "esrd_date", 
  censored_date = HES_CENSORING_DATE, 
  covars_list = c("age_bin", "sex", "log_creatinine"), 
  res_dir = "TTE_esrd/log10creatinine"
)

# time from diagnosis to death 
# model 0 with log10(creatinine)
compare_cox_reg(
  train_df, test_df,
  start_date = "ckd_date", 
  end_date = "death_date",
  censored_date = DEATH_CENSORING_DATE, 
  covars_list = c("age_bin", "sex", "log_creatinine"), 
  res_dir = "TTE_death/log10creatinine"
)
