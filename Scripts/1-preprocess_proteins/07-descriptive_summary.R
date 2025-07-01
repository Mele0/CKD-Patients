# OLINK: Descriptive summary


rm(list=ls())

## Set up data and output directories
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  data_dir <- args[1]
  print(paste("Data directory:", data_dir))
  output_dir <- args[2]
  print(paste("Output directory:", output_dir))
} else {
  stop(paste0("Invalid number of arguments provided (", length(args),
              ", expected 2). Please specify input and output directories."))
}

setwd(data_dir)
figures_dir <- file.path(output_dir, "Figures")
tables_dir <- file.path(output_dir, "Tables")
dir.create(tables_dir)

## Spearman's correlation matrices ----
library(pheatmap)

proteins = readRDS("Proteins_nd_imputed.rds")

annot = readRDS("Panel_annotation.rds")
annotsub = annot[colnames(proteins)]

panel.colour = c("#f4a9be","#61bf1a","#00afac","#9364cc")

mat_col = data.frame(Panel = annotsub)
rownames(mat_col) = colnames(proteins)
mat_colors = list(Panel = panel.colour)
names(mat_colors$Panel) = levels(annotsub)

cor = cor(proteins, method = "spearman")
height = 1/144*nrow(cor) + 1
width = height + 2
mybreaks = seq(-1, 1, length = 100)
ifelse(dir.exists(figures_dir),"",dir.create(figures_dir))
{png(file.path(figures_dir, "Spearmans_correlation_matrix_hclust.png"), width = width, height = height, unit='in', res=500, pointsize = 12)
  pheatmap(cor,
           cellwidth = 0.5, cellheight = 0.5,
           border_color = NA, breaks = mybreaks,
           treeheight_row = 0, treeheight_col = 0,
           show_rownames = FALSE, show_colnames = FALSE,
           annotation_row = mat_col,annotation_col = mat_col, #annotation_colors = mat_colors,
           legend = TRUE, annotation_legend = TRUE,
           annotation_names_row = FALSE, annotation_names_col = FALSE)
  dev.off()
}

{png(file.path(figures_dir, "Spearmans_correlation_matrix.png"), width = width, height = height, unit='in', res=500, pointsize = 12)
  pheatmap(cor[order(annotsub),order(annotsub)],
           cellwidth = 0.5, cellheight = 0.5,
           border_color = NA, breaks = mybreaks,
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE,
           annotation_row = mat_col, annotation_col = mat_col, #annotation_colors = mat_colors,
           legend = TRUE, annotation_legend = TRUE,
           annotation_names_row = FALSE, annotation_names_col = FALSE)
  dev.off()
}

tmp = cor
tmp[lower.tri(tmp, diag = T)] = NA
high_corr = cbind(rownames(tmp)[which(abs(tmp) > 0.7, arr.ind = TRUE)[,1]],
                  colnames(tmp)[which(abs(tmp) > 0.7, arr.ind = TRUE)[,2]],
                  tmp[which(abs(tmp) > 0.7, arr.ind = TRUE)])

write.csv(high_corr, file.path(tables_dir, "Spearmans_correlation_High.csv"))

tmp = cor
tmp[lower.tri(tmp, diag = T)] = NA
high_corr = cbind(rownames(tmp)[which(abs(tmp) > 0.90, arr.ind = TRUE)[,1]],
                  colnames(tmp)[which(abs(tmp) > 0.90, arr.ind = TRUE)[,2]],
                  tmp[which(abs(tmp) > 0.90, arr.ind = TRUE)])

write.csv(high_corr, file.path(tables_dir, "Spearmans_correlation_Very_high.csv"))