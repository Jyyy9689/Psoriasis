
setwd("/home/datahup/pjy/PSE/Analysis/Rdata")
rm(list = ls())
gc()
# load(file = 'Step2_output.Rdata')
# library(Seurat)
# library(dplyr)
# Idents(object = scRNA_harmony) <- "RNA_snn_res.0.3"
# table(Idents(scRNA_harmony))
# markers = FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
# top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# table(top10$cluster)
# write.csv(markers,file = 'First_markers.csv')
# write.csv(top10,file = 'First_markers_top10.csv')



source("/home/datahup/pjy/BRCA/NK/Rpro/Cibersort.R")
library(preprocessCore)
library(e1071)
library(parallel)
Normal_result <- CIBERSORT("Step5_Healthy_Cell_top50Marker.txt", "Step5_Normal_expr.txt", perm = 1000, QN = T)
PSE_result <- CIBERSORT("Step5_PSO_Cell_top50Marker.txt", "Step5_Pse_expr.txt", perm = 1000, QN = T)

save(Normal_result,PSE_result,
     file = "Step5_CibersortResults.Rdata")





