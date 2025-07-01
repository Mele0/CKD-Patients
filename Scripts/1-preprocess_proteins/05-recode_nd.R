# OLINK: Recoding values below LOD

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

proteins = readRDS("Proteins_na_imputed.rds")

dict_qc = readRDS(file.path(output_dir, "Dictionaries/QC_information_UKB.rds"))
dict_sample = readRDS(file.path(output_dir, "Dictionaries/UKB_sample_information.rds"))

## Recode samples with values below LOD ----
X = proteins

all(rownames(X)==dict_sample$eid)

pb = utils::txtProgressBar(style = 3)
LOD_table = NULL
for(k in 1:ncol(X)){
  qc = dict_qc[make.names(dict_qc$Assay)==colnames(X)[k],]
  LOD_list = NULL
  for(i in 1:nrow(X)){
    LOD = qc$LOD[qc$PlateID==dict_sample$PlateID[i]]
    LOD_list = c(LOD_list, LOD) 
  }
  LOD_table = cbind(LOD_table, LOD_list)
  utils::setTxtProgressBar(pb, k/ncol(proteins))
}

X[X < LOD_table] = NA

## Filter assays not detected (<LOD) >=50% of samples ----
ndprop = apply(X, 2, function(x) sum(is.na(x), na.rm = T)/nrow(X))
prop = sort(ndprop, decreasing = TRUE)

pdf(file.path(output_dir, "Figures/Proportion_low_detection.pdf"), width = 7, height = 5)
par(mar=c(1,5,2,1), pty = "m")
plot(prop, type = "h", ylab="Proportion of non-detects", xlab="", xaxt = "n", lwd=1, ylim = c(0,1),
     col=ifelse(prop>=0.5, yes="tomato", no="black"), cex.lab=1.5)
abline(h=0.5, lty=2, col="darkred", lwd = 1)
for (i in 1:length(prop)){
  axis(side=1, at=i, labels=names(prop)[i], las=2,
       col=ifelse(prop[i]>=0.5, yes="tomato", no="black"),
       col.axis=ifelse(prop[i]>=0.5, yes="tomato", no="black"))
}
dev.off()

proteins_sub = X[,ndprop<0.5]

saveRDS(proteins_sub, "Proteins_nd_filtered.rds")
