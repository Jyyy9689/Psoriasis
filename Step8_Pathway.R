
setwd('/home/datahup/pjy/PSE/Analysis/Rdata/')
rm(list = ls())
gc()

library(Seurat)
SeuratWork <- function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  library(harmony)
  table(x@meta.data$patien)
  x = RunHarmony(x, group.by.vars = 'patien')
  x = RunUMAP(x, reduction = "harmony", dims = 1:25)
  x = RunTSNE(x, reduction = "harmony", dims = 1:25)
  x = FindNeighbors(x, reduction = "harmony", dims = 1:25) 
  x = FindClusters(object = x, resolution = c(seq(0,1,by = 0.1))) #根据不同分辨率对细胞群聚类
}

singlecell_gene_test <- function(SerautObj, 
                                 genes.use, 
                                 group.by=NULL, 
                                 assay = "RNA", 
                                 comp = NULL, 
                                 alpha_start = .05, 
                                 Bonferroni = T,
                                 only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      #t_out = t.test(group1_exp, group2_exp)
      t_out = wilcox.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
  else{
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname]
      group1_exp <- group1_exp[which(group1_exp>0)] 
      
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      group2_exp <- group2_exp[which(group2_exp>0)] 
      
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
    
  }
  
  if (Bonferroni == T){
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
    
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
    
  }
  
  else{
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
  }
  
  return(dfOUT)
}

# 用于寻找 GO ID
findGO <- function(pattern, method = "key"){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method == "key"){
    pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
  } else if(method == "gene"){
    pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  
  colnames(pathways) = "pathway"
  
  if(length(pathways) == 0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
} 


# 获取 GO geneSet
getGO <- function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)     
  } else{
    cat("No results!\n")
  }
} 

source("getGoTerm.R")
#GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
#save(GO_DATA, file = "GO_DATA.RData")
load("GO_DATA.RData") # 载入数据 GO_DATA

# limma 差异分析函数
library(limma)
deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

####分组####
load(file = 'Step5_output.Rdata')
table(Idents(PSO))
Basal <- PSO[, Idents(PSO) %in% 'Basal cell']
Basal
MLgene <- c("DEFB4A", "GJB2", "SERPINB3", "SERPINB13")
gene_dat <- FetchData(Basal,MLgene) # 使用FetchData提取出gene的表达量

# DEFB4A分组
quantile(gene_dat$DEFB4A)
gene_dat1 <- subset(gene_dat, gene_dat$DEFB4A != 0)
gene_dat1$DEFB4A_Group <- ifelse(gene_dat1$DEFB4A >= median(gene_dat1$DEFB4A), 'High', 'Low')
table(gene_dat1$DEFB4A_Group)
gene_dat1 <- data.frame(row.names = rownames(gene_dat1),
                        Cell = rownames(gene_dat1),
                        DEFB4A_Group = gene_dat1$DEFB4A_Group)
# GJB2
quantile(gene_dat$GJB2)
gene_dat2 <- subset(gene_dat, gene_dat$GJB2 != 0)
quantile(gene_dat2$GJB2)
gene_dat2$GJB2_Group <- ifelse(gene_dat2$GJB2 >= median(gene_dat2$GJB2), 'High', 'Low')
table(gene_dat2$GJB2_Group)
gene_dat2 <- data.frame(row.names = rownames(gene_dat2),
                        Cell = rownames(gene_dat2),
                        GJB2_Group = gene_dat2$GJB2_Group)
# SERPINB3
quantile(gene_dat$SERPINB3)
gene_dat3 <- subset(gene_dat, gene_dat$SERPINB3 != 0)
quantile(gene_dat3$SERPINB3)
gene_dat3$SERPINB3_Group <- ifelse(gene_dat3$SERPINB3 >= median(gene_dat3$SERPINB3), 'High', 'Low')
table(gene_dat3$SERPINB3_Group)
gene_dat3 <- data.frame(row.names = rownames(gene_dat3),
                        Cell = rownames(gene_dat3),
                        SERPINB3_Group = gene_dat3$SERPINB3_Group)
# SERPINB13
quantile(gene_dat$SERPINB13)
gene_dat4 <- subset(gene_dat, gene_dat$SERPINB13 != 0)
quantile(gene_dat4$SERPINB13)
gene_dat4$SERPINB13_Group <- ifelse(gene_dat4$SERPINB13 >= median(gene_dat4$SERPINB13), 'High', 'Low')
table(gene_dat4$SERPINB13_Group)
gene_dat4 <- data.frame(row.names = rownames(gene_dat4),
                        Cell = rownames(gene_dat4),
                        SERPINB13_Group = gene_dat4$SERPINB13_Group)
# merge
gene_dat0 <- subset(gene_dat, gene_dat$DEFB4A == 0 & gene_dat$GJB2 == 0 & gene_dat$SERPINB3 == 0 & gene_dat$SERPINB13 == 0)
gene_dat0 <- data.frame(row.names = rownames(gene_dat0),
                        Cell = rownames(gene_dat0),
                        Group = rep('None', nrow(gene_dat0)))
gene_dat_merge <- Reduce(function(x,y) merge(x,y,by='Cell',all=T), list(gene_dat0, gene_dat1, gene_dat2, gene_dat3, gene_dat4), accumulate = F)
gene_dat_merge[is.na(gene_dat_merge)] <- 'None'
table(gene_dat_merge$DEFB4A_Group)
table(gene_dat_merge$GJB2_Group)
table(gene_dat_merge$SERPINB3_Group)
table(gene_dat_merge$SERPINB13_Group)
rownames(gene_dat_merge) <- gene_dat_merge$Cell
gene_dat_merge <- gene_dat_merge[,-c(1:2)]

Basal <- AddMetaData(Basal, gene_dat_merge)


####DEFB4A####
table(Basal@meta.data[["DEFB4A_Group"]])
DEFB4A <- Basal[, Basal@meta.data[["DEFB4A_Group"]] != 'None']
DEFB4A
DEFB4A <- SeuratWork(DEFB4A)
Idents(DEFB4A) <- DEFB4A@meta.data$DEFB4A_Group
table(Idents(DEFB4A))

library(ggplot2)
library(ggpubr)
library(ggsignif)
# UMAP
plot1 = DimPlot(DEFB4A, reduction = "umap",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot1
ggsave(filename = "DEFB4A_UMAP.pdf", plot = plot1, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# TSNE
plot2 = DimPlot(DEFB4A, reduction = "tsne",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot2
ggsave(filename = "DEFB4A_TSNE.pdf", plot = plot2, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# Vln
A <- singlecell_gene_test(DEFB4A, 
                          genes.use = 'DEFB4A',
                          group.by = 'DEFB4A_Group', 
                          comp = c("High", "Low"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(DEFB4A, 
                         cols = c("#D62728FF", "#17BECFFF"),
                         pt.size = 0,
                         group.by = "DEFB4A_Group",
                         features = 'DEFB4A', 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)
data <- plots_violins[[1]]$data
colnames(data)[1] <- 'gene'
plots_violins[[1]]+ theme_classic() + 
  theme(axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        axis.title.y= element_text(size=12,color="black"),
        axis.title.x = element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  scale_x_discrete(labels = c("DEFB4A_Group_High","DEFB4A_Group_Low"))+
  geom_signif(annotations = anno_sig,
              y_position = max(data$gene)+0.5,
              xmin = 1,
              xmax = 2,
              tip_length = 0)
ggsave(filename = "DEFB4A_VlnPlot.pdf", width = 6,height = 5,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# FeatupePlot
library(viridis)
FeaturePlot(DEFB4A, features = 'DEFB4A',reduction = 'tsne') + scale_color_viridis(option="C") # 更改颜色
ggsave(filename = "DEFB4A_FeaturePlot.pdf", width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# DEG
DEG_DEFB4A <- FindMarkers(DEFB4A, ident.1 = 'High', ident.2 = 'Low', logfc.threshold = 0.3, min.pct = 0.2, only.pos = F)
table(DEG_DEFB4A$p_val<0.05)
DEG_DEFB4A_p <- subset(DEG_DEFB4A, DEG_DEFB4A$p_val < 0.05)
write.csv(DEG_DEFB4A_p, file = 'Step8_DEG_DEFB4A_p.csv')

# Pathway
# library(clusterProfiler)
# library(org.Hs.eg.db)
# # deg_entr <- bitr(rownames(DEG_DEFB4A_p), fromType = "SYMBOL",
# #                  toType = c("ENTREZID"),
# #                  OrgDb = org.Hs.eg.db)
# 
# DEG_DEFB4A_p_change <- read.csv(file = 'Step8_DEG_DEFB4A_p_change.csv')
# table(duplicated(DEG_DEFB4A_p_change$external_gene_name))
# DEG_DEFB4A_p_change <- DEG_DEFB4A_p_change[!duplicated(DEG_DEFB4A_p_change$external_gene_name),]
# table(is.na(DEG_DEFB4A_p_change$entrezgene_id))
# gene_list <- na.omit(DEG_DEFB4A_p_change$entrezgene_id)
# # KEGG
# DEFB4A_kk_deg <- enrichKEGG(gene = gene_list,
#                      organism = 'hsa',
#                      pvalueCutoff = 0.05)
# # GO
# DEFB4A_go_deg <- enrichGO(gene    = gene_list, 
#                    OrgDb          = org.Hs.eg.db,
#                    ont            = 'ALL', 
#                    pAdjustMethod  = "BH",
#                    pvalueCutoff   = 0.05,
#                    readable       = TRUE)

load(file = 'Step8_DEFB4A_Pathway.Rdata')


##提取p<0.05通路list
library(KEGGREST) 
library(stringr)
# 提取KEGG通路
keggdf = DEFB4A_kk_deg@result %>% filter(pvalue<0.05)
keggpathway = keggdf$ID
kegg_genelist <- list()
for (j in keggpathway) {
  gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路的信息，并缓存
  genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
  genelist <- genes[1:length(genes)%%3 ==2] 
  gs[[1]]$NAME # 通路名称
  kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
}

# 提取GO通路
godf = DEFB4A_go_deg@result %>% filter(pvalue < 0.05)
gopathway <- godf$ID
go_genelist <- list()
for (x in gopathway) {
  a <- getGO(x)
  go_genelist[[x]] = a[[1]]
}


## AUCell
library(AUCell)
cells_rankings <- AUCell_buildRankings(DEFB4A@assays$RNA@data, splitByBlocks=TRUE)
cells_AUC_kegg <- AUCell_calcAUC(kegg_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_kegg)
cells_AUC_go <- AUCell_calcAUC(go_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_go)
#提取打分
DEFB4A_AUC_df_kegg <- getAUC(cells_AUC_kegg)
rownames(DEFB4A_AUC_df_kegg) = str_sub(rownames(DEFB4A_AUC_df_kegg), start = 1, end = -24)
DEFB4A_AUC_df_go <- getAUC(cells_AUC_go)
rownames(DEFB4A_AUC_df_go) = ifelse(rownames(DEFB4A_AUC_df_go) %in% godf$ID, godf$Description, 'unknow')


save(Basal, 
     DEFB4A, DEG_DEFB4A, DEG_DEFB4A_p, DEFB4A_kk_deg, DEFB4A_go_deg, DEFB4A_AUC_df_kegg, DEFB4A_AUC_df_go,
     file = 'Step8_Basal.Rdata')







####GJB2####
table(Basal@meta.data[["GJB2_Group"]])
GJB2 <- Basal[, Basal@meta.data[["GJB2_Group"]] != 'None']
GJB2
GJB2 <- SeuratWork(GJB2)
Idents(GJB2) <- GJB2@meta.data$GJB2_Group
table(Idents(GJB2))

library(ggplot2)
library(ggpubr)
library(ggsignif)
# UMAP
plot1 = DimPlot(GJB2, reduction = "umap",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot1
ggsave(filename = "GJB2_UMAP.pdf", plot = plot1, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# TSNE
plot2 = DimPlot(GJB2, reduction = "tsne",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot2
ggsave(filename = "GJB2_TSNE.pdf", plot = plot2, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# Vln
A <- singlecell_gene_test(GJB2, 
                          genes.use = 'GJB2',
                          group.by = 'GJB2_Group', 
                          comp = c("High", "Low"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(GJB2, 
                         cols = c("#D62728FF", "#17BECFFF"),
                         pt.size = 0,
                         group.by = "GJB2_Group",
                         features = 'GJB2', 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)
data <- plots_violins[[1]]$data
colnames(data)[1] <- 'gene'
plots_violins[[1]]+ theme_classic() + 
  theme(axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        axis.title.y= element_text(size=12,color="black"),
        axis.title.x = element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  scale_x_discrete(labels = c("GJB2_Group_High","GJB2_Group_Low"))+
  geom_signif(annotations = anno_sig,
              y_position = max(data$gene)+0.5,
              xmin = 1,
              xmax = 2,
              tip_length = 0)
ggsave(filename = "GJB2_VlnPlot.pdf", width = 6,height = 5,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# FeatupePlot
library(viridis)
FeaturePlot(GJB2, features = 'GJB2',reduction = 'tsne') + scale_color_viridis(option="C") # 更改颜色
ggsave(filename = "GJB2_FeaturePlot.pdf", width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# DEG
DEG_GJB2 <- FindMarkers(GJB2, ident.1 = 'High', ident.2 = 'Low', logfc.threshold = 0.3, min.pct = 0.2, only.pos = F)
table(DEG_GJB2$p_val<0.05)
DEG_GJB2_p <- subset(DEG_GJB2, DEG_GJB2$p_val < 0.05)
write.csv(DEG_GJB2_p, file = 'Step8_DEG_GJB2_p.csv')

# Pathway
# library(clusterProfiler)
# library(org.Hs.eg.db)
# # deg_entr <- bitr(rownames(DEG_GJB2_p), fromType = "SYMBOL",
# #                  toType = c("ENTREZID"),
# #                  OrgDb = org.Hs.eg.db)
# 
# DEG_GJB2_p_change <- read.csv(file = 'Step8_DEG_GJB2_p_change.csv')
# table(duplicated(DEG_GJB2_p_change$external_gene_name))
# DEG_GJB2_p_change <- DEG_GJB2_p_change[!duplicated(DEG_GJB2_p_change$external_gene_name),]
# table(is.na(DEG_GJB2_p_change$entrezgene_id))
# gene_list <- na.omit(DEG_GJB2_p_change$entrezgene_id)
# # KEGG
# GJB2_kk_deg <- enrichKEGG(gene = gene_list,
#                             organism = 'hsa',
#                             pvalueCutoff = 0.05)
# # GO
# GJB2_go_deg <- enrichGO(gene    = gene_list, 
#                           OrgDb          = org.Hs.eg.db,
#                           ont            = 'ALL', 
#                           pAdjustMethod  = "BH",
#                           pvalueCutoff   = 0.05,
#                           readable       = TRUE)

load(file = 'Step8_GJB2_Pathway.Rdata')


##提取p<0.05通路list
library(KEGGREST) 
library(stringr)
# 提取KEGG通路
keggdf = GJB2_kk_deg@result %>% filter(pvalue<0.05)
keggpathway = keggdf$ID
kegg_genelist <- list()
for (j in keggpathway) {
  gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路的信息，并缓存
  genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
  genelist <- genes[1:length(genes)%%3 ==2] 
  gs[[1]]$NAME # 通路名称
  kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
}

# 提取GO通路
godf = GJB2_go_deg@result %>% filter(pvalue < 0.05)
gopathway <- godf$ID
go_genelist <- list()
for (x in gopathway) {
  a <- getGO(x)
  go_genelist[[x]] = a[[1]]
}


## AUCell
library(AUCell)
cells_rankings <- AUCell_buildRankings(GJB2@assays$RNA@data, splitByBlocks=TRUE)
cells_AUC_kegg <- AUCell_calcAUC(kegg_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_kegg)
cells_AUC_go <- AUCell_calcAUC(go_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_go)
#提取打分
GJB2_AUC_df_kegg <- getAUC(cells_AUC_kegg)
rownames(GJB2_AUC_df_kegg) = str_sub(rownames(GJB2_AUC_df_kegg), start = 1, end = -24)
GJB2_AUC_df_go <- getAUC(cells_AUC_go)
rownames(GJB2_AUC_df_go) = ifelse(rownames(GJB2_AUC_df_go) %in% godf$ID, godf$Description, 'unknow')


save(Basal, 
     DEFB4A, DEG_DEFB4A, DEG_DEFB4A_p, DEFB4A_kk_deg, DEFB4A_go_deg, DEFB4A_AUC_df_kegg, DEFB4A_AUC_df_go,
     GJB2, DEG_GJB2, DEG_GJB2_p, GJB2_kk_deg, GJB2_go_deg, GJB2_AUC_df_kegg, GJB2_AUC_df_go,
     file = 'Step8_Basal.Rdata')





####SERPINB3####
table(Basal@meta.data[["SERPINB3_Group"]])
SERPINB3 <- Basal[, Basal@meta.data[["SERPINB3_Group"]] != 'None']
SERPINB3
SERPINB3 <- SeuratWork(SERPINB3)
Idents(SERPINB3) <- SERPINB3@meta.data$SERPINB3_Group
table(Idents(SERPINB3))

library(ggplot2)
library(ggpubr)
library(ggsignif)
# UMAP
plot1 = DimPlot(SERPINB3, reduction = "umap",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot1
ggsave(filename = "SERPINB3_UMAP.pdf", plot = plot1, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# TSNE
plot2 = DimPlot(SERPINB3, reduction = "tsne",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot2
ggsave(filename = "SERPINB3_TSNE.pdf", plot = plot2, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# Vln
A <- singlecell_gene_test(SERPINB3, 
                          genes.use = 'SERPINB3',
                          group.by = 'SERPINB3_Group', 
                          comp = c("High", "Low"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(SERPINB3, 
                         cols = c("#D62728FF", "#17BECFFF"),
                         pt.size = 0,
                         group.by = "SERPINB3_Group",
                         features = 'SERPINB3', 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)
data <- plots_violins[[1]]$data
colnames(data)[1] <- 'gene'
plots_violins[[1]]+ theme_classic() + 
  theme(axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        axis.title.y= element_text(size=12,color="black"),
        axis.title.x = element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  scale_x_discrete(labels = c("SERPINB3_Group_High","SERPINB3_Group_Low"))+
  geom_signif(annotations = anno_sig,
              y_position = max(data$gene)+0.5,
              xmin = 1,
              xmax = 2,
              tip_length = 0)
ggsave(filename = "SERPINB3_VlnPlot.pdf", width = 6,height = 5,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# FeatupePlot
library(viridis)
FeaturePlot(SERPINB3, features = 'SERPINB3',reduction = 'tsne') + scale_color_viridis(option="C") # 更改颜色
ggsave(filename = "SERPINB3_FeaturePlot.pdf", width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# DEG
DEG_SERPINB3 <- FindMarkers(SERPINB3, ident.1 = 'High', ident.2 = 'Low', logfc.threshold = 0.3, min.pct = 0.2, only.pos = F)
table(DEG_SERPINB3$p_val<0.05)
DEG_SERPINB3_p <- subset(DEG_SERPINB3, DEG_SERPINB3$p_val < 0.05)
write.csv(DEG_SERPINB3_p, file = 'Step8_DEG_SERPINB3_p.csv')

# Pathway
# library(clusterProfiler)
# library(org.Hs.eg.db)
# # deg_entr <- bitr(rownames(DEG_SERPINB3_p), fromType = "SYMBOL",
# #                  toType = c("ENTREZID"),
# #                  OrgDb = org.Hs.eg.db)
# 
# DEG_SERPINB3_p_change <- read.csv(file = 'Step8_DEG_SERPINB3_p_change.csv')
# table(duplicated(DEG_SERPINB3_p_change$external_gene_name))
# DEG_SERPINB3_p_change <- DEG_SERPINB3_p_change[!duplicated(DEG_SERPINB3_p_change$external_gene_name),]
# table(is.na(DEG_SERPINB3_p_change$entrezgene_id))
# gene_list <- na.omit(DEG_SERPINB3_p_change$entrezgene_id)
# # KEGG
# SERPINB3_kk_deg <- enrichKEGG(gene = gene_list,
#                             organism = 'hsa',
#                             pvalueCutoff = 0.05)
# # GO
# SERPINB3_go_deg <- enrichGO(gene    = gene_list, 
#                           OrgDb          = org.Hs.eg.db,
#                           ont            = 'ALL', 
#                           pAdjustMethod  = "BH",
#                           pvalueCutoff   = 0.05,
#                           readable       = TRUE)

load(file = 'Step8_SERPINB3_Pathway.Rdata')


##提取p<0.05通路list
library(KEGGREST) 
library(stringr)
# 提取KEGG通路
keggdf = SERPINB3_kk_deg@result %>% filter(pvalue<0.05)
keggpathway = keggdf$ID
kegg_genelist <- list()
for (j in keggpathway) {
  gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路的信息，并缓存
  genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
  genelist <- genes[1:length(genes)%%3 ==2] 
  gs[[1]]$NAME # 通路名称
  kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
}

# 提取GO通路
godf = SERPINB3_go_deg@result %>% filter(pvalue < 0.05)
gopathway <- godf$ID
go_genelist <- list()
for (x in gopathway) {
  a <- getGO(x)
  go_genelist[[x]] = a[[1]]
}


## AUCell
library(AUCell)
cells_rankings <- AUCell_buildRankings(SERPINB3@assays$RNA@data, splitByBlocks=TRUE)
cells_AUC_kegg <- AUCell_calcAUC(kegg_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_kegg)
cells_AUC_go <- AUCell_calcAUC(go_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_go)
#提取打分
SERPINB3_AUC_df_kegg <- getAUC(cells_AUC_kegg)
rownames(SERPINB3_AUC_df_kegg) = str_sub(rownames(SERPINB3_AUC_df_kegg), start = 1, end = -24)
SERPINB3_AUC_df_go <- getAUC(cells_AUC_go)
rownames(SERPINB3_AUC_df_go) = ifelse(rownames(SERPINB3_AUC_df_go) %in% godf$ID, godf$Description, 'unknow')


save(Basal, 
     DEFB4A, DEG_DEFB4A, DEG_DEFB4A_p, DEFB4A_kk_deg, DEFB4A_go_deg, DEFB4A_AUC_df_kegg, DEFB4A_AUC_df_go,
     GJB2, DEG_GJB2, DEG_GJB2_p, GJB2_kk_deg, GJB2_go_deg, GJB2_AUC_df_kegg, GJB2_AUC_df_go,
     SERPINB3, DEG_SERPINB3, DEG_SERPINB3_p, SERPINB3_kk_deg, SERPINB3_go_deg, SERPINB3_AUC_df_kegg, SERPINB3_AUC_df_go,
     file = 'Step8_Basal.Rdata')




####SERPINB13####
table(Basal@meta.data[["SERPINB13_Group"]])
SERPINB13 <- Basal[, Basal@meta.data[["SERPINB13_Group"]] != 'None']
SERPINB13
SERPINB13 <- SeuratWork(SERPINB13)
Idents(SERPINB13) <- SERPINB13@meta.data$SERPINB13_Group
table(Idents(SERPINB13))

library(ggplot2)
library(ggpubr)
library(ggsignif)
# UMAP
plot1 = DimPlot(SERPINB13, reduction = "umap",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot1
ggsave(filename = "SERPINB13_UMAP.pdf", plot = plot1, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# TSNE
plot2 = DimPlot(SERPINB13, reduction = "tsne",label= T, label.box = T, raster = F, cols = c("#17BECFFF","#D62728FF")) + NoLegend()
plot2
ggsave(filename = "SERPINB13_TSNE.pdf", plot = plot2, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# Vln
A <- singlecell_gene_test(SERPINB13, 
                          genes.use = 'SERPINB13',
                          group.by = 'SERPINB13_Group', 
                          comp = c("High", "Low"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(SERPINB13, 
                         cols = c("#D62728FF", "#17BECFFF"),
                         pt.size = 0,
                         group.by = "SERPINB13_Group",
                         features = 'SERPINB13', 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)
data <- plots_violins[[1]]$data
colnames(data)[1] <- 'gene'
plots_violins[[1]]+ theme_classic() + 
  theme(axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        axis.title.y= element_text(size=12,color="black"),
        axis.title.x = element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  scale_x_discrete(labels = c("SERPINB13_Group_High","SERPINB13_Group_Low"))+
  geom_signif(annotations = anno_sig,
              y_position = max(data$gene)+0.5,
              xmin = 1,
              xmax = 2,
              tip_length = 0)
ggsave(filename = "SERPINB13_VlnPlot.pdf", width = 6,height = 5,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# FeatupePlot
library(viridis)
FeaturePlot(SERPINB13, features = 'SERPINB13',reduction = 'tsne') + scale_color_viridis(option="C") # 更改颜色
ggsave(filename = "SERPINB13_FeaturePlot.pdf", width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step8/")

# DEG
DEG_SERPINB13 <- FindMarkers(SERPINB13, ident.1 = 'High', ident.2 = 'Low', logfc.threshold = 0.25, min.pct = 0.1, only.pos = F)
table(DEG_SERPINB13$p_val<0.05)
DEG_SERPINB13_p <- subset(DEG_SERPINB13, DEG_SERPINB13$p_val < 0.05)
write.csv(DEG_SERPINB13_p, file = 'Step8_DEG_SERPINB13_p.csv')

# Pathway
# library(clusterProfiler)
# library(org.Hs.eg.db)
# # deg_entr <- bitr(rownames(DEG_SERPINB13_p), fromType = "SYMBOL",
# #                  toType = c("ENTREZID"),
# #                  OrgDb = org.Hs.eg.db)
# 
# DEG_SERPINB13_p_change <- read.csv(file = 'Step8_DEG_SERPINB13_p_change.csv')
# table(duplicated(DEG_SERPINB13_p_change$external_gene_name))
# DEG_SERPINB13_p_change <- DEG_SERPINB13_p_change[!duplicated(DEG_SERPINB13_p_change$external_gene_name),]
# table(is.na(DEG_SERPINB13_p_change$entrezgene_id))
# gene_list <- na.omit(DEG_SERPINB13_p_change$entrezgene_id)
# # KEGG
# SERPINB13_kk_deg <- enrichKEGG(gene = gene_list,
#                             organism = 'hsa',
#                             pvalueCutoff = 0.05)
# # GO
# SERPINB13_go_deg <- enrichGO(gene    = gene_list, 
#                           OrgDb          = org.Hs.eg.db,
#                           ont            = 'ALL', 
#                           pAdjustMethod  = "BH",
#                           pvalueCutoff   = 0.05,
#                           readable       = TRUE)

load(file = 'Step8_SERPINB13_Pathway.Rdata')


##提取p<0.05通路list
library(KEGGREST) 
library(stringr)
# 提取KEGG通路
keggdf = SERPINB13_kk_deg@result %>% filter(pvalue<0.05)
keggpathway = keggdf$ID
kegg_genelist <- list()
for (j in keggpathway) {
  gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路的信息，并缓存
  genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
  genelist <- genes[1:length(genes)%%3 ==2] 
  gs[[1]]$NAME # 通路名称
  kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
}

# 提取GO通路
godf = SERPINB13_go_deg@result %>% filter(pvalue < 0.05)
gopathway <- godf$ID
go_genelist <- list()
for (x in gopathway) {
  a <- getGO(x)
  go_genelist[[x]] = a[[1]]
}


## AUCell
library(AUCell)
cells_rankings <- AUCell_buildRankings(SERPINB13@assays$RNA@data, splitByBlocks=TRUE)
cells_AUC_kegg <- AUCell_calcAUC(kegg_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_kegg)
cells_AUC_go <- AUCell_calcAUC(go_genelist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC_go)
#提取打分
SERPINB13_AUC_df_kegg <- getAUC(cells_AUC_kegg)
rownames(SERPINB13_AUC_df_kegg) = str_sub(rownames(SERPINB13_AUC_df_kegg), start = 1, end = -24)
SERPINB13_AUC_df_go <- getAUC(cells_AUC_go)
rownames(SERPINB13_AUC_df_go) = ifelse(rownames(SERPINB13_AUC_df_go) %in% godf$ID, godf$Description, 'unknow')


save(Basal, 
     DEFB4A, DEG_DEFB4A, DEG_DEFB4A_p, DEFB4A_kk_deg, DEFB4A_go_deg, DEFB4A_AUC_df_kegg, DEFB4A_AUC_df_go,
     GJB2, DEG_GJB2, DEG_GJB2_p, GJB2_kk_deg, GJB2_go_deg, GJB2_AUC_df_kegg, GJB2_AUC_df_go,
     SERPINB3, DEG_SERPINB3, DEG_SERPINB3_p, SERPINB3_kk_deg, SERPINB3_go_deg, SERPINB3_AUC_df_kegg, SERPINB3_AUC_df_go,
     SERPINB13, DEG_SERPINB13, DEG_SERPINB13_p, SERPINB13_kk_deg, SERPINB13_go_deg, SERPINB13_AUC_df_kegg, SERPINB13_AUC_df_go,
     file = 'Step8_Basal.Rdata')




####高低组通路活性差异显著的交集####
library(dplyr)
load(file = 'Step8_Basal.Rdata')
MLgene <- c("DEFB4A", "GJB2", "SERPINB3", "SERPINB13")
kegg_deg_list <- list()
go_deg_list <- list()
for (i in MLgene) {
  # 表达矩阵
  dat <- get(paste0(i,'_AUC_df_kegg'))
  keggdf = as.data.frame(dat)
  #rownames(keggdf) = paste0('KEGG_',rownames(keggdf))
  dat1 <- get(paste0(i,'_AUC_df_go'))
  godf = as.data.frame(dat1)
  #rownames(godf) = paste0('GO_',rownames(godf))
  # 分组矩阵
  seu <- get(i)
  Group = as.data.frame(Idents(seu))
  names(Group)[1] <- 'group'
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # analyse
  kegg_deg = deg(keggdf,design,contrast.matrix)
  go_deg = deg(godf,design,contrast.matrix)
  
  kegg_deg1 = kegg_deg %>% filter(P.Value<0.05) 
  kegg_deg1 = kegg_deg1[order(kegg_deg1$logFC,decreasing = T),]
  kegg_deg_list[[i]] = kegg_deg1
  
  go_deg1 = go_deg %>% filter(P.Value<0.05) 
  go_deg1 = go_deg1[order(go_deg1$logFC,decreasing = T),]
  go_deg_list[[i]] = go_deg1
}

inter_kegg <- Reduce(function(x,y) intersect(x, y), list(rownames(kegg_deg_list[["DEFB4A"]]), rownames(kegg_deg_list[["GJB2"]]),
                                                         rownames(kegg_deg_list[["SERPINB3"]]), rownames(kegg_deg_list[["SERPINB13"]])))

inter_go <- Reduce(function(x,y) intersect(x, y), list(rownames(go_deg_list[["DEFB4A"]]), rownames(go_deg_list[["GJB2"]]),
                                                         rownames(go_deg_list[["SERPINB3"]]), rownames(go_deg_list[["SERPINB13"]])))

library(venn)
Venn_list <- list(DEFB4A_kegg = rownames(kegg_deg_list[["DEFB4A"]]), GJB2_kegg = rownames(kegg_deg_list[["GJB2"]]), 
                  SERPINB13_kegg = rownames(kegg_deg_list[["SERPINB3"]]), SERPINB3_kegg= rownames(kegg_deg_list[["SERPINB13"]]))
pdf(file = './Fig/Step8/Venn_KEGG.pdf',width = 10,height = 8)
venn(Venn_list,
     zcolor=c('green','purple','blue','red'), # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1,     # 数字大小
     sncs = 1        # 组名字体大小
)
dev.off()

Venn_list1 <- list(DEFB4A_GO = rownames(go_deg_list[["DEFB4A"]]), GJB2_GO = rownames(go_deg_list[["GJB2"]]), 
                   SERPINB13_GO = rownames(go_deg_list[["SERPINB3"]]), SERPINB3_GO= rownames(go_deg_list[["SERPINB13"]]))
pdf(file = './Fig/Step8/Venn_GO.pdf',width = 10,height = 8)
venn(Venn_list1,
     zcolor=c('green','purple','blue','red'), # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1,     # 数字大小
     sncs = 1        # 组名字体大小
)
dev.off()


# Plot
for (i in MLgene) {
  dat <- get(paste0(i,'_AUC_df_kegg'))
  df1 <- as.data.frame(t(dat[inter_kegg,]))
  df1$Cell <- rownames(df1)
  seu <- get(i)
  CellGroup <- as.data.frame(Idents(seu))
  CellGroup$Cell <- rownames(CellGroup)
  names(CellGroup)[1] <- 'Group'
  df2 <- merge(df1, CellGroup, by = 'Cell')
  library(reshape2)
  df3 <- melt(df2, id.vars = c('Cell', 'Group'))
  df3$variable <- paste0('KEGG_',df3$variable)
  
  dat1 <- get(paste0(i,'_AUC_df_go'))
  df4 <- as.data.frame(t(dat1[inter_go,]))
  df4$Cell <- rownames(df4)
  df5 <- merge(df4, CellGroup, by = 'Cell')
  df6 <- melt(df5, id.vars = c('Cell', 'Group'))
  df6$variable <- paste0('GO_',df6$variable)
  
  df7 <- rbind(df3, df6)
  library(ggsci)
  library(ggpubr)
  library(introdataviz)
  colnames(df7)
  mypalette <- pal_jco()(10)
  #show_col(mypalette)
  plot <- ggplot(df7,aes(x = reorder(variable, value), y = value, fill = Group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.5) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = "Pathway", y = "AUC Score", title = paste0('Pathway Activity with ',i, ' Group')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette)+
    stat_compare_means(aes(group = Group),method = "t.test",label.y = 1.1,label = 'p.signif') + # wilcox.test
    coord_flip() 
  ggsave(filename = paste0('Pathway_Score_',i,'.pdf'), plot, width = 8,height = 10, path = './Fig/Step8')
  
}

save(Basal, 
     DEFB4A, DEG_DEFB4A, DEG_DEFB4A_p, DEFB4A_kk_deg, DEFB4A_go_deg, DEFB4A_AUC_df_kegg, DEFB4A_AUC_df_go,
     GJB2, DEG_GJB2, DEG_GJB2_p, GJB2_kk_deg, GJB2_go_deg, GJB2_AUC_df_kegg, GJB2_AUC_df_go,
     SERPINB3, DEG_SERPINB3, DEG_SERPINB3_p, SERPINB3_kk_deg, SERPINB3_go_deg, SERPINB3_AUC_df_kegg, SERPINB3_AUC_df_go,
     SERPINB13, DEG_SERPINB13, DEG_SERPINB13_p, SERPINB13_kk_deg, SERPINB13_go_deg, SERPINB13_AUC_df_kegg, SERPINB13_AUC_df_go,
     kegg_deg_list, go_deg_list, inter_kegg, inter_go,
     file = 'Step8_Basal.Rdata')




####Bulk RNA中验证高低分组后的通路活性####
keggid <- subset(DEFB4A_kk_deg, DEFB4A_kk_deg@result$Description %in% inter_kegg)
keggid <- keggid$ID
goid <- subset(DEFB4A_go_deg, DEFB4A_go_deg@result$Description %in% inter_go)
goid <- goid$ID

library(KEGGREST) 
library(stringr)
# 提取KEGG通路
keggid_list <- list()
for (j in keggid) {
  gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路的信息，并缓存
  genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
  genelist <- genes[1:length(genes)%%3 ==2] 
  gs[[1]]$NAME # 通路名称
  keggid_list[[gs[[1]]$NAME]] <- genelist # 保存为list
}

# 提取GO通路
goid_list <- list()
for (x in goid) {
  a <- getGO(x)
  goid_list[[x]] = a[[1]]
}

## ssGSEA评分
load(file = 'Step4_PseMerge.Rdata')
table(PSE_phe$diagnosis)
pso_cli <- subset(PSE_phe, PSE_phe$diagnosis == 'Psoriasis')
pso_expr <- as.data.frame(t(PSE_expr_limma[,rownames(pso_cli)]))
pso_expr[1:4,1:4]

Group_exp <- data.frame(row.names = rownames(pso_expr),
                        DEFB4A = pso_expr$DEFB4A,
                        GJB2 = pso_expr$GJB2,
                        SERPINB3 = pso_expr$SERPINB3,
                        SERPINB13 = pso_expr$SERPINB13 )
Group_exp$DEFB4A_Group <- ifelse(Group_exp$DEFB4A >= median(Group_exp$DEFB4A), 'High', 'Low')
table(Group_exp$DEFB4A_Group)
Group_exp$GJB2_Group <- ifelse(Group_exp$GJB2 >= median(Group_exp$GJB2), 'High', 'Low')
table(Group_exp$GJB2_Group)
Group_exp$SERPINB3_Group <- ifelse(Group_exp$SERPINB3 >= median(Group_exp$SERPINB3), 'High', 'Low')
table(Group_exp$SERPINB3_Group)
Group_exp$SERPINB13_Group <- ifelse(Group_exp$SERPINB13 >= median(Group_exp$SERPINB13), 'High', 'Low')
table(Group_exp$SERPINB13_Group)
Group_exp1 <- Group_exp[,5:8]
Group_exp1$sample <- rownames(Group_exp1)
names(Group_exp1) <- c("DEFB4A","GJB2", "SERPINB3","SERPINB13","sample" )


library(GSVA)
pso_expr1 <- as.matrix(t(pso_expr))
ssgesa_kegg <- gsva(pso_expr1, keggid_list, kcdf="Gaussian",method = "ssgsea",parallel.sz=14) #parallel.sz为线程数
rownames(ssgesa_kegg) = str_sub(rownames(ssgesa_kegg), start = 1, end = -24)
ssgesa_go <- gsva(pso_expr1, goid_list, kcdf="Gaussian",method = "ssgsea",parallel.sz=14) #parallel.sz为线程
rownames(ssgesa_go) = ifelse(rownames(ssgesa_go) %in% DEFB4A_go_deg@result$ID, DEFB4A_go_deg@result$Description, 'unknow')

ssgesa_kegg1 <- as.data.frame(t(ssgesa_kegg))
ssgesa_kegg1$sample <- rownames(ssgesa_kegg1)
ssgesa_go1 <- as.data.frame(t(ssgesa_go))
ssgesa_go1$sample <- rownames(ssgesa_go1)

for (i in MLgene) {
  Group_exp2 <- data.frame(row.names = rownames(Group_exp1),
                           Group = Group_exp1[,i],
                           sample = rownames(Group_exp1))
  ssgesa_kegg2 <- merge(ssgesa_kegg1, Group_exp2, by = 'sample')
  df <- melt(ssgesa_kegg2, id.vars = c('sample', 'Group'))
  df$variable <- paste0('KEGG_',df$variable)
  
  ssgesa_go2 <- merge(ssgesa_go1, Group_exp2, by = 'sample')
  df1 <- melt(ssgesa_go2, id.vars = c('sample', 'Group'))
  df1$variable <- paste0('GO_',df1$variable)
  
  df2 <- rbind(df, df1)
  mypalette <- pal_aaas()(10)
  #show_col(mypalette)
  plot <- ggplot(df2,aes(x = reorder(variable, value), y = value, fill = Group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.5) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = "Pathway", y = "ssGSEA Score", title = paste0('Pathway Activity with ',i, ' Group')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[2:3])+
    stat_compare_means(aes(group = Group),method = "wilcox.test",label.y = 1.1,label = 'p.signif') + # wilcox.test
    coord_flip() 
  ggsave(filename = paste0('Pathway_Score_',i,'_Bulk.pdf'), plot, width = 8,height = 10, path = './Fig/Step8')
  
}



####交集通路AUC热图####
load(file = 'Step8_Basal.Rdata')
library(stringr)
aucMat_list <- list()
for (i in MLgene) {
  seu <- get(i)
  aucMat0 <- as.data.frame(Idents(seu))
  aucMat0$Cell <- rownames(aucMat0)
  names(aucMat0)[1] <- 'Group'
  aucMat0$Group <- paste0(i,'_',aucMat0$Group)
  
  df_kegg <- get(paste0(i,"_AUC_df_kegg"))
  aucMat_kegg <- as.data.frame(t(df_kegg[inter_kegg,]))
  colnames(aucMat_kegg) <- paste0('KEGG_',colnames(aucMat_kegg))
  aucMat_kegg$Cell <- rownames(aucMat_kegg)
  
  df_go <- get(paste0(i,'_AUC_df_go'))
  aucMat_go <- as.data.frame(t(df_go[inter_go,]))
  colnames(aucMat_go) <- paste0('GO_',colnames(aucMat_go))
  aucMat_go$Cell <- rownames(aucMat_go)
  
  aucMat1 <- merge(aucMat_kegg, aucMat_go, by = 'Cell')
  aucMat2 <- merge(aucMat1, aucMat0, by = 'Cell')
  rownames(aucMat2) <- aucMat2$Cell
  aucMat2 <- aucMat2[,-1]
  
  aucMat3 <- apply(aucMat2[,1:(ncol(aucMat2)-1)], 2 , function(x){aggregate(x, by=list(type=aucMat2$Group),mean)})# aggregate 分组求均值(mean)/和(sum)/计数(length)
  aucMat4 <- do.call(cbind, aucMat3)
  
  rownames(aucMat4) <- aucMat4[,1]
  aucMat4 <- aucMat4[,-seq(1, ncol(aucMat4), by = 2)] #删除单数列
  colnames(aucMat4) <- str_sub(colnames(aucMat4), 1, -3)
  
  aucMat_list[[i]] <- aucMat4
}

aucMat5 <- do.call(rbind, aucMat_list)
rownames(aucMat5) <- str_split(rownames(aucMat5), '[.]', simplify = T)[,2]
aucMat6 <- as.matrix(t(aucMat5))
aucMat6[1:4,1:4]
# DEFB4A_High DEFB4A_Low GJB2_High  GJB2_Low
# KEGG_Staphylococcus aureus infection   0.1576396  0.1508788 0.1535812 0.1406447
# KEGG_Estrogen signaling pathway        0.2055750  0.2317163 0.2128228 0.2167632
# GO_keratinization                      0.2035761  0.1736243 0.1569944 0.1321636
# GO_peptide cross-linking               0.2491576  0.2256340 0.1975440 0.1771622
library(pheatmap)
library(viridis)
p_heat <- pheatmap(aucMat6, show_colnames = T,show_rownames = T,
                   scale = 'row',
                   #annotation_col = anno, 
                   cluster_rows = F, cluster_cols = F,
                   gaps_row = 2,
                   gaps_col = c(2,4,6),
                   color=colorRampPalette(c("#17BECFFF", 'white', "#D62728FF"))(100),
                   angle_col = 90,
                   main = 'Pathway AUCell Score'
)
ggsave(filename = 'AUCell_pheatmap.pdf', p_heat, width = 8, height = 6, path = "./Fig/Step8")






