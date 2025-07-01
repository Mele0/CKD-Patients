save_calibration_plot <- function(sharp_var_selection, clusterID, fig_dir) {
  sharp_stability_calibration_plot_filename <- file.path(fig_dir, paste0("sharp_stability_calibration_cluster_", clusterID, ".png"))
  png(filename = sharp_stability_calibration_plot_filename, width = 1000, height = 800)
  par(mar = c(8, 8, 10, 8) + 0.1)  # Bottom, Left, Top, Right margins
  CalibrationPlot(sharp_var_selection)
  mtext(paste("Sharp Stability Calibration - cluster", clusterID), side = 3, line = 6.5, cex = 2)
  dev.off()
}

save_selprop_plot <- function(selprop, calibrated_params, clusterID, fig_dir, top=NULL) {
  selprop_nonzero <- selprop[selprop != 0]
  selprop_sorted <- sort(selprop_nonzero, decreasing = TRUE)  # Sorting in ascending order
  if (!is.null(top)) {
    selprop_sorted = selprop_sorted[1:top]
    sharp_selection_prop_plot_filename <- file.path(fig_dir, paste0("sharp_selection_prop_top_", top, "cluster_",clusterID,".png"))
  } else {
    sharp_selection_prop_plot_filename <- file.path(fig_dir, paste0("sharp_selection_prop_cluster_",clusterID,".png"))
  }

  png(filename = sharp_selection_prop_plot_filename, width = 1000, height = 800)
  par(mar = c(10, 5, 1, 1), oma = c(0, 0, 6, 0))  # Increase the top outer margin
  
  
  plot(selprop_sorted, type = "h", lwd = 3, las = 1, ylim=c(0, 1), xlab = "",
       ylab = "Selection Proportion", xaxt = "n", col = ifelse(selprop_sorted >= calibrated_params[2], 
                                                               yes = "red", no = "grey"), cex.lab = 1.5)
  
  mtext(paste("Stability Lasso Selected Proteins Cluster", clusterID, "Membership"), 
        side = 3, line = 2, cex = 2, outer = TRUE) 
  mtext(paste("λ =", round(calibrated_params[1], 5), "; π =", round(calibrated_params[2], 2)), 
        side = 3, line = 0, cex = 1, outer = TRUE) 
  
  abline(h = calibrated_params[2], lty = 2, col = "darkred")
  
  for (i in 1:length(selprop_sorted)) {
    axis(side = 1, at = i, labels = names(selprop_sorted)[i],
         las = 2, col = ifelse(selprop_sorted[i] >= calibrated_params[2],
                               yes = "red", no = "grey"), col.axis = ifelse(selprop_sorted[i] >=
                                                                              calibrated_params[2], yes = "red", no = "grey"))
  }
  
  dev.off()
}

save_stable_proteins <- function(selprop, calibrated_params, clusterID, fig_dir) {
  stable_mask <- selprop >= calibrated_params[2]
  selprop_stable <- selprop[stable_mask]

  sink(file.path(fig_dir, paste0("stable_proteins_cluster_", clusterID, ".txt")))
  print(paste("Non-zero selection proportion proteins for cluster", clusterID))
  print(table(selprop != 0))
  cat("\n")
  print(paste("Stably selected proteins for cluster", clusterID))
  print(table(stable_mask))
  cat("\n")
  print(paste(names(selprop_stable), collapse = ", "))
  cat("\n")
  print("Calibrated params:")
  print(calibrated_params)
  cat("\n")
  print(selprop)
  sink(NULL)
}

save_stable_lasso_model <- function(model, confusion_matrix, clusterID, fig_dir) {
  saveRDS(model, paste0("stable_lasso_model_cluster_", clusterID, ".rds"))

  sink(file.path(fig_dir, paste0("model_summary_cluster_", clusterID, ".txt")))
  print(paste("Logistic regression for cluster", clusterID, "membership"))
  print(summary(model))
  cat("\n")
  print("Odds ratios and 95% CIs:")
  odds_ratios <- tryCatch(exp(coef(model)), error = function(e) {
    cat("Error getting ORs:", conditionMessage(e), "\n")
    NA
  })
  conf_ints <- tryCatch(exp(confint(model)), error = function(e) {
    cat("Error getting CIs:", conditionMessage(e), "\n")
    NA
  })
  tryCatch(print(cbind(odds_ratios, conf_ints)), error = function(e) {
    cat("Error combining OR and CIs:", conditionMessage(e), "\n")
    print(odds_ratios)
    print(conf_ints)
  })
  cat("\n")
  print("Test data confusion matrix:")
  print(confusion_matrix)
  sink(NULL)
}

save_betas_sorted_plot <- function(model, clusterID, fig_dir) {
  betas <- summary(model)$coefficients[-1, 1] # remove intercept, get est coeffs

  if (length(betas) > 0) {
    beta_plot_filename <- file.path(fig_dir, paste0("betas_sorted_cluster_", clusterID, ".png"))
    betas_sort_idx <- sort(betas, decreasing = TRUE, index.return = TRUE)$ix  # Sorting in ascending order
    png(filename = beta_plot_filename, width = 800, height = 600)
    par(mar = c(5, 5, 4, 2))
    plot(betas[betas_sort_idx], type = 'h', col='navy', lwd=3, xaxt='n', xlab='', ylab=expression(beta))
    axis(side = 1, at = 1:length(betas), labels = names(betas)[betas_sort_idx], las=2)
    title(main=paste("Cluster", clusterID, "Lasso Coefficients"))
    abline(h=0, lty=2)
    dev.off()
  } else {
    print(paste("No proteins selected, skipping beta plot for cluster", clusterID))
  }
}

save_roc_curve <- function(roc_curve, auc, clusterID, fig_dir) {
  roc_cluster_plot_filename <- file.path(fig_dir, paste0("roc_cluster_", clusterID, ".png"))
  png(filename = roc_cluster_plot_filename, width = 800, height = 800)
  par(pty="s", oma = c(1, 1, 8, 1))
  plot(roc_curve, asp = 1)
  mtext(paste("ROC for Proteins Predicting Cluster", clusterID, "Membership"), side = 3, line = 2, cex = 2, outer = TRUE)
  mtext(sprintf("AUC = %.3f", auc), side = 3, line = 0.5, cex = 1.5, outer = TRUE)
  dev.off()
}
