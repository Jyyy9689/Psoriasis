

####Merge celldata####
rm(list = ls())
gc()
library(Seurat)
load(file = 'Step2_output.Rdata')
scRNA_harmony
table(Idents(scRNA_harmony))
First = as.data.frame(Idents(scRNA_harmony))
First$cell = rownames(First)
names(First) = c('FirstAnnotation', 'Cell')
str(First)
First$FirstAnnotation = as.character(First$FirstAnnotation)
table(First$FirstAnnotation)

# Lymph
load(file = 'Step2_Lymph.Rdata')
table(Idents(Lymph))
Second_Lymph = as.data.frame(Idents(Lymph))
Second_Lymph$cell = rownames(Second_Lymph)
names(Second_Lymph) = c('SecondAnnotation', 'Cell')
str(Second_Lymph)
Second_Lymph$SecondAnnotation = as.character(Second_Lymph$SecondAnnotation)
table(Second_Lymph$SecondAnnotation)

# Myeloid
load(file = 'Step2_Myeloid.Rdata')
table(Idents(Myeloid))
Second_Mye = as.data.frame(Idents(Myeloid))
Second_Mye$cell = rownames(Second_Mye)
names(Second_Mye) = c('SecondAnnotation', 'Cell')
str(Second_Mye)
Second_Mye$SecondAnnotation = as.character(Second_Mye$SecondAnnotation)
table(Second_Mye$SecondAnnotation)


# Merge
library(dplyr)
Second = rbind(Second_Mye,Second_Lymph)
table(Second$SecondAnnotation)
Second1 = full_join(First,Second)
Second2 = subset(Second1, is.na(Second1$SecondAnnotation))
table(Second2$FirstAnnotation)
Second2$SecondAnnotation = Second2$FirstAnnotation
rownames(Second2) = Second2$Cell
Second3 = subset(Second1, !is.na(Second1$SecondAnnotation))
rownames(Second3) = Second3$Cell
Second4 = rbind(Second3,Second2)
table(Second4$SecondAnnotation)
Second5 = data.frame(row.names = rownames(Second4), SecondAnnotation = as.factor(Second4$SecondAnnotation))

scRNA_harmony = AddMetaData(scRNA_harmony, metadata = Second5)
table(scRNA_harmony@meta.data[["SecondAnnotation"]])
Idents(scRNA_harmony) = scRNA_harmony@meta.data[["SecondAnnotation"]]
save(scRNA_harmony,file = 'Step2_output.Rdata')


####Plot####
## 1.堆积直方图
## 病人细胞丰度 
rm(list = ls())
gc()
load(file = 'Step2_output.Rdata')
b = as.data.frame(table(scRNA_harmony@meta.data$patien))
rownames(b) = b$Var1
v = as.character(b$Var1)
pat = list()
for (i in v) {
  patien = scRNA_harmony[,scRNA_harmony@meta.data$patien %in% i]
  df1 = as.data.frame(table(Idents(patien)))
  df1$id = rep(i,length(nrow(df1)))
  df1$pre = df1$Freq / b[i,2]
  pat[[i]] = df1
}
# liat -> dataframe
pat2 = do.call(rbind,pat)
names(pat2) = c("Cell Type" ,"Freq" ,"Patient"  , "Percentage" )
unique(pat2$Patient)
rownames(pat2) <- NULL
pat2[pat2=='Ctrl1'] <- 'GSM4946161_Control1'
pat2[pat2=='Ctrl2'] <- 'GSM4946162_Control2'
pat2[pat2=='Ctrl3'] <- 'GSM4946163_Control3'
pat2[pat2=='Psor1'] <- 'GSM4946164_Psoriasis1'
pat2[pat2=='Psor2'] <- 'GSM4946165_Psoriasis2'
pat2[pat2=='Psor3'] <- 'GSM4946166_Psoriasis3'




library(ggsci)
pal1 = pal_ucscgb()(10)
pal2 = pal_d3('category20', alpha = 0.8)(20)
pal3 = c(pal1,pal2)
p1 = ggplot(pat2) +
  geom_bar(aes(x=Patient, y=Percentage, fill=`Cell Type`), stat="identity") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold.italic", vjust = 0.5, hjust = 1, angle = 0),
        axis.title.y = element_text(size = 14, color = "black",face = "bold"),
        axis.title.x = element_text(size = 14, color = "black",face = "bold"),
        legend.title = element_text(color="black",size=14,face="bold"),
        legend.text = element_text(color="black", size = 12, face = "bold"),
        legend.position = 'top'
  )+
  scale_fill_manual(values=pal3) 
  #coord_flip()
p1
ggsave(filename = "All_Patient.pdf",p1,width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step3")


##不同组别间细胞丰度
b = as.data.frame(table(scRNA_harmony@meta.data$status))
rownames(b) = b$Var1
v = as.character(b$Var1)
pat = list()
for (i in v) {
  status = scRNA_harmony[,scRNA_harmony@meta.data$status %in% i]
  df1 = as.data.frame(table(Idents(status)))
  df1$id = rep(i,length(nrow(df1)))
  df1$pre = df1$Freq / b[i,2]
  pat[[i]] = df1
}
# liat -> dataframe
pat2 = do.call(rbind,pat)
names(pat2) = c("Cell Type" ,"Freq" ,"Group"  , "Percentage" )

# 直方图
pat3 <- subset(pat2, pat2$Group == 'Healthy')
pat3$Percentage <- pat3$Percentage * 100
library(ggpubr)
ggplot(pat3,aes(x = reorder(`Cell Type`,`Percentage`) , y = Percentage, fill = `Cell Type`)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 16) +
  labs(x = "Cell Type", y = 'Percentage(Sub_cell/Healthy_TotalCell,%)', title = 'Healthy') +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold", vjust = 0.5, hjust = 1, angle = 0),
        axis.title.y = element_text(size = 14, color = "black",face = "bold"),
        axis.title.x = element_text(size = 14, color = "black",face = "bold"),
        plot.title  = element_text(size = 14, color = "black",face = "bold",vjust = 0.5, hjust = 0.5),
        legend.title = element_text(color="black",size=14,face="bold"),
        legend.text = element_text(color="black", size = 12, face = "bold"),
        legend.position = 'none'
  )+
  scale_fill_manual(values = pal3)+
  coord_flip()

ggsave(filename = 'Healthy_bar.pdf',width = 8,height = 10,path = '/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step3')


pat3 <- subset(pat2, pat2$Group == 'PSO')
pat3$Percentage <- pat3$Percentage * 100
library(ggpubr)
ggplot(pat3,aes(x = reorder(`Cell Type`,`Percentage`) , y = Percentage, fill = `Cell Type`)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 16) +
  labs(x = "Cell Type", y = 'Percentage(Sub_cell/Psoriasis_TotalCell,%)', title = 'Psoriasis') +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size = 12, color = "black",face = "bold"),
        axis.text.y = element_text(size = 12, color = "black",face = "bold", vjust = 0.5, hjust = 1, angle = 0),
        axis.title.y = element_text(size = 14, color = "black",face = "bold"),
        axis.title.x = element_text(size = 14, color = "black",face = "bold"),
        plot.title  = element_text(size = 14, color = "black",face = "bold",vjust = 0.5, hjust = 0.5),
        legend.title = element_text(color="black",size=14,face="bold"),
        legend.text = element_text(color="black", size = 12, face = "bold"),
        legend.position = 'none'
  )+
  scale_fill_manual(values = pal3)+
  coord_flip()

ggsave(filename = 'Psoriasis_bar.pdf',width = 8,height = 10,path = '/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step3')




### 2.箱图
table(scRNA_harmony@meta.data$status)
# Healthy
H <- scRNA_harmony[,scRNA_harmony@meta.data$status %in% 'Healthy']
b = as.data.frame(table(H@meta.data$patien))
rownames(b) = b$Var1
v = as.character(b$Var1)
pat = list()
for (i in v) {
  patien = H[,H@meta.data$patien %in% i]
  df1 = as.data.frame(table(Idents(patien)))
  df1$id = rep(i,length(nrow(df1)))
  df1$pre = df1$Freq / b[i,2]
  pat[[i]] = df1
}
# liat -> dataframe
pat2 = do.call(rbind,pat)
names(pat2) = c("Cell Type" ,"Freq" ,"Patient"  , "Percentage" )
pat2$Group <- rep('Healthy', nrow(pat2))

# PSO
P <- scRNA_harmony[,scRNA_harmony@meta.data$status %in% 'PSO']
b = as.data.frame(table(P@meta.data$patien))
rownames(b) = b$Var1
v = as.character(b$Var1)
pat = list()
for (i in v) {
  patien = P[,P@meta.data$patien %in% i]
  df1 = as.data.frame(table(Idents(patien)))
  df1$id = rep(i,length(nrow(df1)))
  df1$pre = df1$Freq / b[i,2]
  pat[[i]] = df1
}
# liat -> dataframe
pat3 = do.call(rbind,pat)
names(pat3) = c("Cell Type" ,"Freq" ,"Patient"  , "Percentage" )
pat3$Group <- rep('PSO', nrow(pat3))

pat4 <- rbind(pat2, pat3)
library(ggplot2)
library(tidyverse)
library(ggpubr)
# plot
colnames(pat4)
ggplot(pat4,aes(x = reorder(`Cell Type`,-`Percentage`) , y = Percentage, fill = Group)) +
  # 小提琴图层
  geom_violin(alpha = 1.2,
              #width = 1.8,
              trim = T,
              color = NA) +
  # 箱线图图层
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  theme_bw(base_size = 16) +
  labs(x = "Cell Type", y = 'Cell Count') +
  theme(axis.text.x = element_text(angle = 65,hjust = 1,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  scale_fill_manual(values = pal3)+ 
  # 添加显著性标记
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")   #kruskal.test  p.signif
ggsave(filename = 'All_Group_box.pdf',width = 16,height = 10,path = '/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step3')





# 2.Marker 表达量展示
# library(Seurat)
# table(Idents(scRNA_harmony))
# markers = FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
# top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# table(top10$cluster)
# 
# Gene <- subset(top10, top10$cluster=='Regulatory T (Treg) cell') # 依次检查每个簇的top10marker基因
# p_all_markers=DotPlot(scRNA_harmony,
#                       features = mar$gene,
#                       cols = c("lightgrey", "purple"),
#                       scale = T,assay='RNA' )+
#   theme(axis.text.x=element_text(angle=45,hjust = 1))
# p_all_markers
# Gene$gene
# 
# 
# mar <- FindMarkers(scRNA_harmony, only.pos = T, ident.1 = 'CD8+ T cell', min.pct = 0.6, logfc.threshold = 0.5)
# 
# Gene <- c("S100A12","S100A9","S100A8") # 检查感兴趣的基因
# p_all_markers1=DotPlot(scRNA_harmony,
#                        features = rownames(mar),
#                        cols = c("lightgrey", "purple"),
#                        scale = T,assay='RNA' )+
#   theme(axis.text.x=element_text(angle=45,hjust = 1))
# p_all_markers1

gene = c("CD79A", "MS4A1", "BANK1",   # B cell
         "KRT14", "KRT1",  "KRT5",    # Basal cell
         "IL7R",  "LTB",   "TRAC",    # CD4+ T cell
         'TXN', "HLA-DRA", "CD1C",    # Dendritic cell
         'VWF', 'RAMP2', 'GNG11',     # Endothelial cell
         "GBP1", "MX1",  "IFI6" ,     # Exhausted CD8+ T cell
         "DCN", "COL1A1", "COL3A1",   # Fibroblast
         "LCE3D", "CRCT1", "CNFN",    # keratinocyte
         "C1QA", "C1QB", "C1QC",      # Macrophages
         "TPSB2", "HPGD", "TPSAB1",   # Mast cell
         "PPBP", "PF4", 'NOP53',      # Megakaryocyte
         "MKI67", "TOP2A", "CENPF",   # MKI67+ progenitor cell
         "LYZ", "FCN1", 'CTSS',       # Monocyte:CD16-
         "FCGR3A", "LST1", "AIF1",    # Monocyte:CD16+
         "ACTA2", "C11orf96", "TAGLN",# Myofibroblast 
         'TCF7', "FHIT", "MAL",       # Naive CD4+ T cell
         "CD8A",  "CD8B", 'OXNAD1',   # Naive CD8+ T cell
         'GNLY', 'PRF1', 'CD247',     # Natural killer (NK) cell
         "CCL5", 'TRGC2', 'DUSP2',    # Natural killer T (NKT) cell
         "S100A12","S100A9","S100A8", # Neutrophil
         "IGHA1","IGKC","TNFRSF17",   # Plasma cell
         "CCDC50", "IRF8", "PPP1R14B",# Plasmacytoid dendritic cell
         "RTKN2", "IL32", "TRBC2"     # Regulatory T (Treg) cell
         ) 


#library(devtools) # 加载失败的话需要更新rlang包 : install.packages("rlang")
#install_github("jokergoo/ComplexHeatmap")
#devtools::install_github("sajuukLyu/ggunchull", type = "source")
#devtools::install_github('junjunlab/scRNAtoolVis')
library(Seurat)
library(ComplexHeatmap)
library(ggunchull)
library(scRNAtoolVis)
table(Idents(scRNA_harmony))
# gene = c('HLA-DRA','MS4A1', 'CD79A', 'BANK1', 'CD79B',  # B cell
#          "WFDC2", "FDCSP", "KRT17","KRT14", # Cancer Cell
#          'IL7R', 'ANXA1', 'LTB',  # CD4+ T cell
#          'GZMK', 'GZMH', 'DNAJB1',  # CD8+ T cell
#          'CLEC10A', 'CPVL', 'LGALS2', 'CD1A', 'CD1E',  # Dendritic cell
#          'VWF', 'RAMP2', 'GNG11',  # Endothelial Cell
#          "KRT19", "KRT8", "EPCAM","CD24", # Epithelial Cell
#          'CXCL13','NR3C1','ITM2A','DUSP4','RBPJ', # Exhausted T(Tex) cell
#          'IGFBP7', 'ACTA2', 'RGS5', 'SPARCL1','BGN','COL1A2',  # Fibroblast
#          'PLTP','C1QC','C1QB','C1QA','SPP1','CTSD', # Macrophage
#          'TPSAB1', 'TPSB2', 'CPA3', 'HPGDS', 'MS4A2',  # Mast cell
#          'S100A9', 'S100A8', 'LYZ', "S100P","PI3",  # Neutrophil
#          'GNLY', 'TRDC', 'XCL1', 'XCL2', 'KLRD1', 'NKG7', 'KLRC1', 'GZMB',  # NK cell
#          'HIST1H4C', 'TUBA1B', 'TUBB', 'TYMS', 'HMGN2', 'HMGB2',  # NKT cell
#          'MMP9', 'CTSK', 'ACP5',  # Osteoclast
#          'IGKV3-15', 'IGKV3-20', 'IGHV6-1', 'IGKV1-8', 'IGKV1-9',  # Plasma cell
#          'PTGDS', 'PLAC8', 'TCF4', 'LILRA4', 'IRF7', 'TCL1A', 'IRF8', 'TSPAN13', 'PLD4',  # Plasmacytoid dendritic cell
#          'TNFRSF4', 'TNFRSF18', 'BATF', 'TIGIT', 'SPOCK2', 'ARID5B' # Regulatory T(Treg) cell	
# )
library(scales) 
show_col(pal_d3('category20')(20))
pal = pal_d3('category20')(20)
table(duplicated(gene))
pdf(file = '/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step3/All_marker.pdf',width = 10,height = 12)
AverageHeatmap(object = scRNA_harmony,
               markerGene = gene,
               column_names_rot = 90, 
               row_title = 'Marker Genes',
               htCol = c("#D62728FF", "white", "#9467BDFF")
)
dev.off()



