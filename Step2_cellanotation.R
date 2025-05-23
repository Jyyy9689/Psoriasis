rm(list = ls())
gc()

####First_all####
#Seurat Workfolw
library(Seurat)
library(dplyr)
load(file = 'Step1_Finaldata.Rdata')
scRNA_harmony <- pse_merge_list2

# QC
# 计算线粒体基因比例 线粒体是独立遗传的，不是染色体上基因控制的 默认<10%
scRNA_harmony[["percent.mt"]] <- PercentageFeatureSet(scRNA_harmony, pattern = "^MT-")
table(scRNA_harmony[["percent.mt"]] < 10)
scRNA_harmony = subset(scRNA_harmony, subset = percent.mt < 10)
scRNA_harmony

#计算红血细胞基因比例 红细胞没有细胞核，没有转录组 默认<3%
rownames(scRNA_harmony)[grep("^HB[^(p)]", rownames(scRNA_harmony))]
scRNA_harmony <- PercentageFeatureSet(scRNA_harmony, "^HB[^(p)]", col.name = "percent_hb")
table(scRNA_harmony[["percent_hb"]] < 1)
scRNA_harmony = subset(scRNA_harmony, subset = percent_hb < 1)
scRNA_harmony

# 标准工作流
length(table(scRNA_harmony@meta.data[["patien"]]))
scRNA_harmony <- NormalizeData(scRNA_harmony)
scRNA_harmony <- FindVariableFeatures(scRNA_harmony, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA_harmony)
scRNA_harmony <- ScaleData(scRNA_harmony, features = all.genes)
scRNA_harmony <- RunPCA(scRNA_harmony, features = VariableFeatures(object = scRNA_harmony))

# Harmony整合
library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "patien")

#降维聚类
scRNA_harmony = RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:25)
scRNA_harmony = FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:25) 
scRNA_harmony = FindClusters(object = scRNA_harmony, resolution = c(seq(0,1,by = 0.1))) #根据不同分辨率对细胞群聚类
library(clustree)
clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "scRNA_harmony_resolution(0-1).pdf",width = 20,height = 14,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

save(scRNA_harmony,file = 'Step2_output.Rdata')

# 作图
rm(list = ls())
gc()
load(file = 'Step2_output.Rdata')
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.4"
head(Idents(scRNA_harmony), 5)#查看前5个细胞的分类ID
table(Idents(scRNA_harmony))
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T, raster=T, label.box = T) 
plot1
ggsave(filename = "scRNA_harmony_res0.4.pdf", plot = plot1, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")


plot2 = DimPlot(scRNA_harmony, reduction = "umap", label=F, group.by = 'patien', raster=T) 
plot2
ggsave(filename = "scRNA_harmony_patient.pdf", plot = plot2, width = 14,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")


library(ggsci)
Idents(object = scRNA_harmony) <- "status"
head(Idents(scRNA_harmony), 5)#查看前5个细胞的分类ID
table(Idents(scRNA_harmony))
pal = pal_ucscgb(alpha = 0.8)(26)
plot3 = DimPlot(scRNA_harmony, reduction = "umap", label=F, raster=T, cols = pal) 
plot3
ggsave(filename = "scRNA_harmony_group.pdf", plot = plot3, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")


# Find marker
rm(list = ls())
gc()
load(file = 'Step2_output.Rdata')
library(Seurat)
library(dplyr)
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.4"
table(Idents(scRNA_harmony))
markers = FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
table(top10$cluster)
write.csv(markers,file = 'Step2_First_markers.csv')
write.csv(top10,file = 'Step2_First_markers_top10.csv')


# SingleR
sce = scRNA_harmony
library(Seurat)
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.4
clusters=sce@meta.data$seurat_clusters

Blue.ref <- celldex::BlueprintEncodeData()
pred.Blue.ref <- SingleR(test = sce_for_SingleR, ref = Blue.ref, labels = Blue.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

DICE.ref <- celldex::DatabaseImmuneCellExpressionData()
pred.DICE.ref <- SingleR(test = sce_for_SingleR, ref = DICE.ref, labels = DICE.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

HPCA.ref <- celldex::HumanPrimaryCellAtlasData()
pred.HPCA.ref <- SingleR(test = sce_for_SingleR, ref = HPCA.ref, labels = HPCA.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Mona.ref <- celldex::MonacoImmuneData()
pred.Mona.ref <- SingleR(test = sce_for_SingleR, ref = Mona.ref, labels = Mona.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Nover.ref <- celldex::NovershternHematopoieticData()
pred.Nover.ref <- SingleR(test = sce_for_SingleR, ref = Nover.ref, labels = Nover.ref$label.fine ,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    Blue=pred.Blue.ref$labels,
                    DICE=pred.DICE.ref$labels,
                    HPCA=pred.HPCA.ref$labels,
                    Mona=pred.Mona.ref$labels,
                    Nover=pred.Nover.ref$labels )
head(cellType)
sce@meta.data$singleR_Blue=cellType[match(clusters,cellType$ClusterID),'Blue']
sce@meta.data$singleR_DICE=cellType[match(clusters,cellType$ClusterID),'DICE']
sce@meta.data$singleR_HPCA=cellType[match(clusters,cellType$ClusterID),'HPCA']
sce@meta.data$singleR_Nover=cellType[match(clusters,cellType$ClusterID),'Nover']
sce@meta.data$singleR_Mona=cellType[match(clusters,cellType$ClusterID),'Mona']


pro='All_SingleR_anno_res0.4'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

scRNA_harmony = sce
save(scRNA_harmony, file = "Step2_output.Rdata")
write.csv(cellType,file = 'Step2_First_celltype.csv')


# check marker
library(Seurat)
library(ggplot2)
load(file = 'Step2_output.Rdata')
table(Idents(scRNA_harmony))
top10 <- read.csv('Step2_First_markers_top10.csv',row.names = 1)
Gene <- subset(top10, top10$cluster==22) # 依次从0开始检查每个簇的top10marker基因
p_all_markers=DotPlot(scRNA_harmony,
                      features = Gene$gene,
                      cols = c("lightgrey", "purple"),
                      scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
Gene$gene

Gene <- c('FAP','THY1') # 检查感兴趣的基因
p_all_markers1=DotPlot(scRNA_harmony,
                      features = Gene,
                      cols = c("lightgrey", "purple"),
                      scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers1


# 结合Cellmarker数据库和SingleR结果，半监督注释细胞
celltype <- read.csv(file = 'Step2_First_celltype2.csv',header = T)
unique(celltype$Final)
new.cluster.ids = celltype$Final
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
table(Idents(scRNA_harmony))
scRNA_harmony@meta.data$FirstAnnotation = Idents(scRNA_harmony)

#plot
library(ggsci)
pal = pal_ucscgb(alpha = 0.8)(26)
pal = pal[-c(7,8)]
plot = DimPlot(scRNA_harmony, reduction = "umap", label=T, label.box = T, cols = pal, raster = T) #+ NoLegend()
plot
ggsave(filename = "scRNA_harmony_First.pdf", plot = plot, width = 12,height = 8,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

save(scRNA_harmony, file = "Step2_output.Rdata")





####Second_Lymph-sub####
table(Idents(scRNA_harmony))
Lymph = scRNA_harmony[ , Idents(scRNA_harmony) %in% c('T cell', 'NK cell', 'B cell')]
Lymph
table(Idents(Lymph))
rm(scRNA_harmony)
gc()

Lymph <- NormalizeData(Lymph)
Lymph <- FindVariableFeatures(Lymph, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Lymph)
Lymph <- ScaleData(Lymph, features = all.genes)
Lymph <- RunPCA(Lymph, features = VariableFeatures(object = Lymph))
library(harmony)
table(Lymph@meta.data$patien)
Lymph = RunHarmony(Lymph, group.by.vars = 'patien')
Lymph = RunUMAP(Lymph, reduction = "harmony", dims = 1:25)
Lymph = FindNeighbors(Lymph, reduction = "harmony", dims = 1:25) 
Lymph = FindClusters(object = Lymph, resolution = c(seq(0,1,by = 0.1))) #根据不同分辨率对细胞群聚类
library(clustree)
clustree(Lymph@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "Lymph_resolution(0-1).pdf",width = 20,height = 14,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")
save(Lymph,file = 'Step2_Lymph.Rdata')

load(file = 'Step2_Lymph.Rdata')
Idents(object = Lymph) <- "RNA_snn_res.0.5"
Lymph@meta.data$seurat_clusters = Lymph@meta.data$RNA_snn_res.0.5
head(Idents(Lymph), 5)#查看前5个细胞的分类ID

library(ggplot2)
plot1 = DimPlot(Lymph, reduction = "umap", label= T, label.box = T, raster = T) + NoLegend()
plot1
ggsave(filename = "Lymph_res0.5.pdf", plot = plot1, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")


library(dplyr)
Idents(object = Lymph) <- "RNA_snn_res.0.5"
table(Idents(Lymph))
markers = FindAllMarkers(Lymph, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

table(top10$cluster)
write.csv(markers,file = 'Step2_Lymph_markers.csv')
write.csv(top10,file = 'Step2_Lymph_markers_top10.csv')


# SingleR
sce = Lymph
library(Seurat)
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.5
clusters=sce@meta.data$seurat_clusters

Blue.ref <- celldex::BlueprintEncodeData()
pred.Blue.ref <- SingleR(test = sce_for_SingleR, ref = Blue.ref, labels = Blue.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

DICE.ref <- celldex::DatabaseImmuneCellExpressionData()
pred.DICE.ref <- SingleR(test = sce_for_SingleR, ref = DICE.ref, labels = DICE.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

HPCA.ref <- celldex::HumanPrimaryCellAtlasData()
pred.HPCA.ref <- SingleR(test = sce_for_SingleR, ref = HPCA.ref, labels = HPCA.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Mona.ref <- celldex::MonacoImmuneData()
pred.Mona.ref <- SingleR(test = sce_for_SingleR, ref = Mona.ref, labels = Mona.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Nover.ref <- celldex::NovershternHematopoieticData()
pred.Nover.ref <- SingleR(test = sce_for_SingleR, ref = Nover.ref, labels = Nover.ref$label.fine ,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    Blue=pred.Blue.ref$labels,
                    DICE=pred.DICE.ref$labels,
                    HPCA=pred.HPCA.ref$labels,
                    Mona=pred.Mona.ref$labels,
                    Nover=pred.Nover.ref$labels )
head(cellType)
sce@meta.data$singleR_Blue=cellType[match(clusters,cellType$ClusterID),'Blue']
sce@meta.data$singleR_DICE=cellType[match(clusters,cellType$ClusterID),'DICE']
sce@meta.data$singleR_HPCA=cellType[match(clusters,cellType$ClusterID),'HPCA']
sce@meta.data$singleR_Nover=cellType[match(clusters,cellType$ClusterID),'Nover']
sce@meta.data$singleR_Mona=cellType[match(clusters,cellType$ClusterID),'Mona']


pro='Lymph_SingleR_anno_res0.5'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

Lymph = sce
save(Lymph, file = "Step2_Lymph.Rdata")
write.csv(cellType,file = 'Step2_Lymph_celltype.csv')


# check marker
library(Seurat)
library(ggplot2)
load(file = 'Step2_Lymph.Rdata')
table(Idents(Lymph))
top10 <- read.csv('Step2_Lymph_markers_top10.csv',row.names = 1)
Gene <- subset(top10, top10$cluster==2) # 依次从0开始检查每个簇的top10marker基因
p_all_markers=DotPlot(Lymph,
                      features = Gene$gene,
                      cols = c("lightgrey", "purple"),
                      scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
Gene$gene

Gene <- c('CD25','CD4','FOXP3','CD127') # 检查感兴趣的基因
p_all_markers1=DotPlot(Lymph,
                       features = Gene,
                       cols = c("lightgrey", "purple"),
                       scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers1


# 结合Cellmarker数据库和SingleR结果，半监督注释细胞
celltype <- read.csv(file = 'Step2_Lymph_celltype2.csv',header = T)
unique(celltype$Final)
new.cluster.ids = celltype$Final
names(new.cluster.ids) <- levels(Lymph)
Lymph <- RenameIdents(Lymph, new.cluster.ids)
table(Idents(Lymph))
Lymph@meta.data$SecondAnnotation = Idents(Lymph)

#plot
library(ggsci)
pal = pal_d3(alpha = 0.8)(26)
#pal = pal[-c(7,8)]
plot = DimPlot(Lymph, reduction = "umap", label=T, label.box = T, cols = pal, raster = T) #+ NoLegend()
plot
ggsave(filename = "Lymph_Final.pdf", plot = plot, width = 12,height = 8,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

save(Lymph, file = "Step2_Lymph.Rdata")






####Second_Myeloid-sub####
rm(list = ls())
gc()
library(Seurat)
load(file = 'Step2_output.Rdata')
table(Idents(scRNA_harmony))
Myeloid = scRNA_harmony[ , Idents(scRNA_harmony) %in% c('Monocytes', 'Dendritic cell', 'Macrophage', 'Plasmacytoid dendritic cell', 'Mast cell')]
Myeloid
table(Idents(Myeloid))
rm(scRNA_harmony)
gc()

Myeloid <- NormalizeData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Myeloid)
Myeloid <- ScaleData(Myeloid, features = all.genes)
Myeloid <- RunPCA(Myeloid, features = VariableFeatures(object = Myeloid))

library(harmony)
table(Myeloid@meta.data$patien)
Myeloid = RunHarmony(Myeloid, group.by.vars = 'patien')
Myeloid = RunUMAP(Myeloid, reduction = "harmony", dims = 1:25)
Myeloid = FindNeighbors(Myeloid, reduction = "harmony", dims = 1:25) 
Myeloid = FindClusters(object = Myeloid, resolution = c(seq(0,1,by = 0.1))) #根据不同分辨率对细胞群聚类
library(clustree)
clustree(Myeloid@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "Myeloid_resolution(0-1).pdf",width = 20,height = 14,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")
save(Myeloid,file = 'Step2_Myeloid.Rdata')

load(file = 'Step2_Myeloid.Rdata')
Idents(object = Myeloid) <- "RNA_snn_res.0.6"
Myeloid@meta.data$seurat_clusters = Myeloid@meta.data$RNA_snn_res.0.6
head(Idents(Myeloid), 5)#查看前5个细胞的分类ID

library(ggplot2)
plot1 = DimPlot(Myeloid, reduction = "umap", label= T, label.box = T, raster = T) + NoLegend()
plot1
ggsave(filename = "Myeloid_res0.6.pdf", plot = plot1, width = 8,height = 6,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")


library(dplyr)
Idents(object = Myeloid) <- "RNA_snn_res.0.6"
table(Idents(Myeloid))
markers = FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
table(top10$cluster)
write.csv(markers,file = 'Step2_Myeloid_markers.csv')
write.csv(top10,file = 'Step2_Myeloid_markers_top10.csv')


# SingleR
sce = Myeloid
library(Seurat)
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.6
clusters=sce@meta.data$seurat_clusters

Blue.ref <- celldex::BlueprintEncodeData()
pred.Blue.ref <- SingleR(test = sce_for_SingleR, ref = Blue.ref, labels = Blue.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

DICE.ref <- celldex::DatabaseImmuneCellExpressionData()
pred.DICE.ref <- SingleR(test = sce_for_SingleR, ref = DICE.ref, labels = DICE.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

HPCA.ref <- celldex::HumanPrimaryCellAtlasData()
pred.HPCA.ref <- SingleR(test = sce_for_SingleR, ref = HPCA.ref, labels = HPCA.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Mona.ref <- celldex::MonacoImmuneData()
pred.Mona.ref <- SingleR(test = sce_for_SingleR, ref = Mona.ref, labels = Mona.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Nover.ref <- celldex::NovershternHematopoieticData()
pred.Nover.ref <- SingleR(test = sce_for_SingleR, ref = Nover.ref, labels = Nover.ref$label.fine ,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    Blue=pred.Blue.ref$labels,
                    DICE=pred.DICE.ref$labels,
                    HPCA=pred.HPCA.ref$labels,
                    Mona=pred.Mona.ref$labels,
                    Nover=pred.Nover.ref$labels )
head(cellType)
sce@meta.data$singleR_Blue=cellType[match(clusters,cellType$ClusterID),'Blue']
sce@meta.data$singleR_DICE=cellType[match(clusters,cellType$ClusterID),'DICE']
sce@meta.data$singleR_HPCA=cellType[match(clusters,cellType$ClusterID),'HPCA']
sce@meta.data$singleR_Nover=cellType[match(clusters,cellType$ClusterID),'Nover']
sce@meta.data$singleR_Mona=cellType[match(clusters,cellType$ClusterID),'Mona']


pro='Myeloid_SingleR_anno_res0.6'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=T)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 12,height = 10,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

Myeloid = sce
save(Myeloid, file = "Step2_Myeloid.Rdata")
write.csv(cellType,file = 'Step2_Myeloid_celltype.csv')


# check marker
library(Seurat)
library(ggplot2)
load(file = 'Step2_Myeloid.Rdata')
table(Idents(Myeloid))
top10 <- read.csv('Step2_Myeloid_markers_top10.csv',row.names = 1)
Gene <- subset(top10, top10$cluster==20) # 依次从0开始检查每个簇的top10marker基因
p_all_markers=DotPlot(Myeloid,
                      features = Gene$gene,
                      cols = c("lightgrey", "purple"),
                      scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
Gene$gene

Gene <- c('FCGR3A','FCGR3B','CD14') # 检查感兴趣的基因
p_all_markers1=DotPlot(Myeloid,
                       features = Gene,
                       cols = c("lightgrey", "purple"),
                       scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers1


# 结合Cellmarker数据库和SingleR结果，半监督注释细胞
celltype <- read.csv(file = 'Step2_Myeloid_celltype2.csv',header = T)
unique(celltype$Final)
new.cluster.ids = celltype$Final
names(new.cluster.ids) <- levels(Myeloid)
Myeloid <- RenameIdents(Myeloid, new.cluster.ids)
table(Idents(Myeloid))
Myeloid@meta.data$SecondAnnotation = Idents(Myeloid)

#plot
library(ggsci)
pal = pal_d3(alpha = 0.8)(26)
#pal = pal[-c(7,8)]
plot = DimPlot(Myeloid, reduction = "umap", label=T, label.box = T, cols = pal, raster = T) #+ NoLegend()
plot
ggsave(filename = "Myeloid_Final.pdf", plot = plot, width = 12,height = 8,path = "/home/datahup/pjy/PSE/Analysis/Rdata/Fig/Step2")

save(Myeloid, file = "Step2_Myeloid.Rdata")


