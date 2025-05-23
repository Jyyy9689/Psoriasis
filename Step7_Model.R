
setwd('/home/datahup/pjy/PSE/Analysis/Rdata/')
rm(list = ls())
gc()

load(file = 'Step4_PseMerge.Rdata')
load(file = 'Step5_output.Rdata')
load(file = 'Step6_WGCNAOutput.Rdata')

# 筛选特征
library(Seurat)
table(Idents(PSO))
Basal <- FindMarkers(PSO, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, ident.1 = 'Basal cell')
table(Basal$p_val < 0.05)
# table(PSO_marker$cluster)
# Basal <- subset(PSO_marker,PSO_marker$cluster == 'Basal cell')
library(venn)
venn_list <- list(`WGCNA Hub Gene` = rownames(hub1),
                  `Basal Cell DEG` = rownames(Basal))
#作图
pdf(file = './Fig/Step7/Venn_Hub_Basal.pdf',width = 10,height = 8)
venn(venn_list,
     zcolor=c('green','purple'), # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1,     # 数字大小
     sncs = 1        # 组名字体大小
)
dev.off()


# 构建训练队列
inter <- intersect(venn_list[["WGCNA Hub Gene"]], venn_list[["Basal Cell DEG"]])
table(PSE_phe$Type)
expr <- as.data.frame(t(PSE_expr_limma[inter,]))
expr$ID <- rownames(expr)
phe <- data.frame(row.names = rownames(PSE_phe),
                  ID = rownames(PSE_phe),
                  Group = PSE_phe$Type)
expr1 <- merge(expr, phe, by = 'ID')
rownames(expr1) <- expr1$ID
expr1 <- expr1[,-1]
table(expr1$Group)

pse_nor <- subset(expr1, expr1$Group != 'Psoriasis-no-lesion')
pse_nor[pse_nor=='Psoriasis-lesion'] <- 1
pse_nor[pse_nor=='normal-normal'] <- 0
table(pse_nor$Group)
write.csv(pse_nor, file = 'Step7_pse_nor.csv')

pse_lesion <- subset(expr1, expr1$Group != 'normal-normal')
pse_lesion[pse_lesion=='Psoriasis-lesion'] <- 1
pse_lesion[pse_lesion=='Psoriasis-no-lesion'] <- 0
table(pse_lesion$Group)
write.csv(pse_lesion, file = 'Step7_pse_lesion.csv')


# 训练脚本为python
####整理模型训练结果####

# pse_nor
svm_psenor <- read.csv(file = 'Step7_Resluts_SVM_pse_nor.csv')
svm_psenor <- as.data.frame(t(svm_psenor))
svm_psenor$gene <- rownames(svm_psenor)
svm_psenor <- svm_psenor[order(abs(svm_psenor$V1), decreasing = T),]
names(svm_psenor) <- c('weight','gene')

library(plyr)
library(ggplot2)
library(ggsci)
pal_igv()(10)
svm_psenor1 <- svm_psenor[1:30,]
p1 <- ggplot(svm_psenor1,aes(reorder(gene, abs(weight)),abs(weight),fill=gene))+
  geom_bar(stat="identity",position = position_stack(),fill="#5DB1DDFF")+
  scale_fill_grey()+
  #scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Genes",y="Weight", fill="", title = 'Support Vector Machine Top30 Feature Weights')+
  geom_text(aes(label=abs(round(weight,4))),vjust=0.5,size=4,color="black")+
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) + #也就加上这一行
  coord_flip() 
p1


rfc_psenor <- read.csv(file = 'Step7_Resluts_RFC_pse_nor.csv')
rfc_psenor <- as.data.frame(t(rfc_psenor))
rfc_psenor$gene <- rownames(rfc_psenor)
rfc_psenor <- rfc_psenor[order(abs(rfc_psenor$V1), decreasing = T),]
names(rfc_psenor) <- c('Feature Importance','gene')

rfc_psenor1 <- rfc_psenor[1:30,]
p2 <- ggplot(rfc_psenor1,aes(reorder(gene, `Feature Importance`),`Feature Importance`,fill=gene))+
  geom_bar(stat="identity",position = position_stack(),fill="#802268FF")+
  scale_fill_grey()+
  #scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Genes",y="Feature Importance", fill="", title = 'Random Forest Top30 Features Importance')+
  geom_text(aes(label=abs(round(`Feature Importance`,4))),vjust=0.5,size=4,color="black")+
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip() 
p2

p1+p2
ggsave(filename = './Fig/Step7/Pse_Nor_bar.pdf',width = 14,height = 8)


# pse_lesion
svm_pselesion <- read.csv(file = 'Step7_Resluts_SVM_pse_lesion.csv')
svm_pselesion <- as.data.frame(t(svm_pselesion))
svm_pselesion$gene <- rownames(svm_pselesion)
svm_pselesion <- svm_pselesion[order(abs(svm_pselesion$V1), decreasing = T),]
names(svm_pselesion) <- c('weight', 'gene')

svm_pselesion1 <- svm_pselesion[1:30,]
p3 <- ggplot(svm_pselesion1,aes(reorder(gene, abs(weight)),abs(weight),fill=gene))+
  geom_bar(stat="identity",position = position_stack(),fill="#5DB1DDFF")+
  scale_fill_grey()+
  #scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Genes",y="Weight", fill="", title = 'Support Vector Machine Top30 Feature Weights')+
  geom_text(aes(label=abs(round(weight,4))),vjust=0.5,size=4,color="black")+
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) + #也就加上这一行
  coord_flip() 
p3


rfc_pselesion <- read.csv(file = 'Step7_Resluts_RFC_pse_lesion.csv')
rfc_pselesion <- as.data.frame(t(rfc_pselesion))
rfc_pselesion$gene <- rownames(rfc_pselesion)
rfc_pselesion <- rfc_pselesion[order(abs(rfc_pselesion$V1), decreasing = T),]
names(rfc_pselesion) <- c('Feature Importance','gene')

rfc_pselesion1 <- rfc_pselesion[1:30,]
p4 <- ggplot(rfc_pselesion1,aes(reorder(gene, `Feature Importance`),`Feature Importance`,fill=gene))+
  geom_bar(stat="identity",position = position_stack(),fill="#802268FF")+
  scale_fill_grey()+
  #scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Genes",y="Feature Importance", fill="", title = 'Random Forest Top30 Features Importance')+
  geom_text(aes(label=abs(round(`Feature Importance`,4))),vjust=0.5,size=4,color="black")+
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) + #也就加上这一行
  coord_flip() 
p4

p3+p4
ggsave(filename = './Fig/Step7/Pse_lesion_bar.pdf',width = 14,height = 8)


# 取交集
pse_nor <- intersect(svm_psenor1$gene, rfc_psenor1$gene)
pse_lesion <- intersect(svm_pselesion1$gene, rfc_pselesion1$gene)
inter <- intersect(pse_nor, pse_lesion)
inter
#  "DEFB4A"    "GJB2"      "SERPINB3"  "SERPINB13"

library(venn)
venn_list <- list(`SVM (Lesion-Normal)` = svm_psenor1$gene,
                  `RFC (Lesion-Normal)` = rfc_psenor1$gene,
                  `SVM (Lesion-NoLesion)` = svm_pselesion1$gene,
                  `RFC (Lesion-NoLesion)` = rfc_pselesion1$gene
                  )
#作图
pdf(file = './Fig/Step7/Venn_Gene.pdf',width = 10,height = 8)
venn(venn_list,
     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1,     # 数字大小
     sncs = 1        # 组名字体大小
)
dev.off()

# 查看交集详情,并导出结果
library(VennDiagram)
inter1 <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter1)) inter1[i,'values'] <- paste(inter1[[i,'..values..']], collapse = '|')
inter1 <- subset(inter1, select = -..values.. )



####结果可视化####
# Bulk RNA-seq
load(file = 'Step4_PseMerge.Rdata')
PSE_expr_limma[1:4,1:4]
df0 <- as.data.frame(t(PSE_expr_limma))
df <- data.frame(DEFB4A = df0$DEFB4A,
                 GJB2 = df0$GJB2,
                 SERPINB3 = df0$SERPINB3,
                 SERPINB13 = df0$SERPINB13,
                 row.names = rownames(df0)
                 )
df$Group <- ifelse(rownames(df) %in% rownames(PSE_phe) ,PSE_phe$Type, 'none')
table(df$Group)

# box
library(reshape2)
# 使用melt函数将宽格式转换为长格式
melted_df <- melt(df, id.vars = "Group")
library(ggsci)
library(ggpubr)
mypalette = pal_d3("category20")(20)
mypalette1 = pal_ucscgb()(10)
ggplot(melted_df,aes(variable,value,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Genes", y = "Expression", title = 'Psoriasis Key Driver Genes Expression in Bulk RNA-seq') +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = c("#6699FFFF", "#FF0000FF",  "#CC33FFFF")) + 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test") # wilcox.test
ggsave(filename = './Fig/Step7/Bulk_box.pdf',width = 10,height = 6)


# cor 
load(file = 'Step5_CibersortResults.Rdata')
Ciber <- rbind(PSE_result,Normal_result)
Ciber <- as.data.frame(Ciber[,1:23])
Ciber$sample <- rownames(Ciber)
df$sample <- rownames(df)
Ciber1 <- merge(df, Ciber, by = 'sample')
rownames(Ciber1) <- Ciber1$sample
Ciber1 <- Ciber1[,-c(1,6)]

Cor1 = cor(Ciber1, method = "pearson")
head(Cor1)
Cor2 <- Cor1[1:4,]
library(corrplot)
pdf(file = './Fig/Step7/Bulk_Cor_pearson.pdf',width = 14,height = 12)
corrplot(Cor2, type="upper",
         # addCoef.col = "black", #添加相关系数
         #title = 'Pearson Correlation',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="pie",
         col = colorRampPalette(c("purple","white","firebrick3"))(100),
         # addCoefasPercent = T,
         diag = F)
dev.off()




# scRNA-seq
library(Seurat)
load(file = 'Step5_output.Rdata')
MLgene <- c("DEFB4A", "GJB2", "SERPINB3","SERPINB13")
# table(Idents(PSO))
# PSO@meta.data[["seurat_clusters"]] <- PSO@meta.data[["SecondAnnotation"]]
# vln.dat=FetchData(PSO,c(MLgene,"seurat_clusters")) # 使用FetchData提取出marker gene的表达量
# vln.dat$Cell <- rownames(vln.dat)
# #宽数据转长数据
# vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell","seurat_clusters"), 
#                                measure.vars = MLgene,
#                                variable.name = "gene", 
#                                value.name = "Expr") %>%
#   group_by(seurat_clusters,gene) %>% #分组
#   mutate(fillcolor=mean(Expr)) #计算均值
# 
# library(cowplot)
# #mycolor <- c(pal_d3('category20')(20))
# p1 <- ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
#   geom_violin(scale = "width", adjust = 1, trim = TRUE) +
#   scale_y_continuous(expand = c(0, 0), position="right", labels = function(x) {
#     if (length(x) >= 2) {
#       return(c(rep(x = "", times = length(x) - 2), x[length(x) - 1], ""))
#     } else {
#       return(x)
#     }
#   }) +
#   facet_grid(rows = vars(seurat_clusters), scales = "free", switch = "y") +
#   #scale_fill_manual(values = mycolor) + 
#   theme_cowplot(font_size = 12) +
#   theme(legend.position = "none", panel.spacing = unit(0, "lines"),
#         plot.title = element_text(hjust = 0.5),
#         panel.background = element_rect(fill = NA, color = "black"),
#         plot.margin = margin(7, 7, 0, 7, "pt"),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold"),
#         strip.text.y.left = element_text(angle = 0),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
#   ) +
#   ggtitle("Marker Gene") + ylab("Expression Level")
# p1

p1 = VlnPlot(PSO, features = MLgene[1], pt.size = 0, raster = T) + NoLegend() + labs(x='')
p1

p2 = VlnPlot(PSO, features = MLgene[2], pt.size = 0, raster = T) + NoLegend() + labs(x='')
p2

p3 = VlnPlot(PSO, features = MLgene[3], pt.size = 0, raster = T) + NoLegend() + labs(x='')
p3

p4 = VlnPlot(PSO, features = MLgene[4], pt.size = 0, raster = T) + NoLegend() + labs(x='')
p4

library(patchwork)
pt <- (p1 + p2 + p3 + p4) + plot_layout(ncol = 2)
ggsave(filename = './Fig/Step7/scRNA_Vln_AllCell.pdf',pt, width = 10,height = 10)
pt2 <- (p1 + p2 + p3 + p4) + plot_layout(ncol = 4)
ggsave(filename = './Fig/Step7/scRNA_Vln_AllCell1.pdf',pt2, width = 20,height = 6)




load(file = 'Step2_output.Rdata')
table(Idents(scRNA_harmony))
Basal <- scRNA_harmony[, Idents(scRNA_harmony) %in% 'Basal cell']
table(Basal@meta.data$status)
# VlnPlot(Basal, features =MLgene,pt.size=0,group.by = 'status')+ 
#   NoLegend() +
#   stat_compare_means(group = status,method = "t.test",label = 'p.signif')  # wilcox.test


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
      # t_out = t.test(group1_exp, group2_exp)
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
A <- singlecell_gene_test(Basal, 
                          genes.use = MLgene,
                          group.by = 'status', 
                          comp = c("PSO", "Healthy"))
# A1 <- singlecell_gene_test(mouse_data,
#                            genes.use = c('S100a8','Ltf','Ncf1','Ly6g','Anxa1','Il1b'),
#                            group.by = 'orig.ident', 
#                            comp = c("10X_ntph_F", "10X_ntph_M"),
#                            only_postive = T)

anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(Basal, 
                         cols = c("limegreen", "purple"),
                         pt.size = 0,
                         group.by = "status",
                         features = MLgene, 
                         ncol = 3, 
                         log = FALSE,
                         combine = FALSE)

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none',
          plot.title = element_text(hjust = 0.5, face = "bold")
          )+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Healthy","Psoriasis"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

pt3 <- (plots_violins[[1]] + plots_violins[[2]] + plots_violins[[3]] + plots_violins[[4]]) + plot_layout(ncol = 4)
#CombinePlots(plots_violins)
ggsave(filename = './Fig/Step7/scRNA_Vln_Group.pdf', pt3, width = 16,height = 4)



# Basal拟时序
# Step1 Seurat -> monocle
library(monocle)
Basal_PSO <- PSO[, Idents(PSO) %in% 'Basal cell']
data <- as(as.matrix(Basal_PSO@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Basal_PSO@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

# Step2 估计size factor和离散度
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
# 过滤低质量的细胞
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))


# Step3 细胞分类
# Clustering cells without marker genes  无监督方法
HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)

# 可视化 tSNE
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 12,
                        reduction_method = 'tSNE', verbose = T,
                        check_duplicates = FALSE)
HSMM <- clusterCells(HSMM, num_clusters = 4)
plot_cell_clusters(HSMM)


# Step4 构建轨迹
# Step4.1: 选择定义过程的基因 三种方法均为无监督

#使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(NK)
# diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
# HSMM <- setOrderingFilter(HSMM, diff.genes)
# plot_ordering_genes(HSMM)

##使用seurat选择的高变基因
# var.seurat <- VariableFeatures(NK)
# # var.seurat <- NK@assays[["RNA"]]@var.features
# HSMM <- setOrderingFilter(HSMM, var.seurat)
# plot_ordering_genes(HSMM)

##使用monocle选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


# Step4.2 降维
HSMM <- reduceDimension(HSMM, max_components = 3,
                        num_dim = 5,
                        method = 'DDRTree')

# Step4.3 按照轨迹排序细胞
HSMM <- orderCells(HSMM) 
# 每次重新加载monocle包都会报以下错误：
# Error in if (class(projection) != "matrix") projection <- as.matrix(projection) : 
#   the condition has length > 
# 解决办法：
#trace('project2MST', edit = T, where = asNamespace("monocle")) # 在弹出窗口中把if (elass(projection) != "matrix")注释，并保存，重新运行即可

# 使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
# cds <- orderCells(cds, root_state = 3) #把State3设成拟时间轴的起始点
## !!!若时序得分和预期得分颠倒，可将Pseudotime的值进行反向转换，将Pseudotime列的值更新为最大值减去原始Pseudotime值。即原始值越大，转换后的值越小，原始值越小，转换后的值越大。
# pData(cds)$Pseudotime <- max(pData(cds)$Pseudotime) - pData(cds)$Pseudotime

library(ggplot2)
library(scales)
library(ggsci)
library(viridis)
pal = pal_d3('category20')(20)

plot_cell_trajectory(HSMM, color_by = "Pseudotime",
                     show_branch_points = F # 去除根标记
) + 
  scale_color_viridis(option="C") # 更改颜色
ggsave(filename = 'Monocle_pseudotime.pdf',width = 8,height = 6,path = './Fig/Step7/')


# plot_cell_trajectory(HSMM, color_by = "seurat_clusters",
#                      show_branch_points = F) + 
#   scale_color_discrete(names('')) + # 图注标题隐藏
#   scale_color_manual(values = pal) 
# ggsave(filename = 'Monocle_NKthird.pdf',width = 8,height = 6,path = '/home/datahup/pjy/BRCA/NK/Fig/PlanC/2.1/')


# plot_cell_trajectory(HSMM, color_by = "seurat_clusters", 
#                      show_branch_points = F)  + 
#   scale_color_manual(values = pal) +
#   scale_color_discrete(names('')) +
#   facet_wrap(~seurat_clusters, nrow = 1) # 逐个查看
# ggsave(filename = 'Monocle_NKthird_split.pdf',width = 20,height = 6,path = '/home/datahup/pjy/BRCA/NK/Fig/PlanC/2.1/')


Time_diff <- differentialGeneTest(HSMM[MLgene,], cores = 12, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
#gene_short_name <- fData$gene_short_name
#Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# num_clusters为人为设置的聚类
p=plot_pseudotime_heatmap(HSMM[MLgene,],
                          num_clusters=1,
                          show_rownames= T,
                          return_heatmap= T,
                          #cluster_rows = F,
                          #add_annotation_row =annotation_row,
                          #hmcols = c("#2CA02CFF", "#9467BDFF")
)
p
ggsave("Monocle_MLgene.pdf", p, width = 6, height = 4,path = './Fig/Step7/')


# 美化拟时序热图
Time_genes <- MLgene
{
  newdata <- data.frame(Pseudotime = seq(min(pData(HSMM)$Pseudotime),
                                         max(pData(HSMM)$Pseudotime),
                                         length.out = 100)) # 默认值: 把时序切割成100个连续变量
  m <- genSmoothCurves(HSMM[Time_genes,],
                       trend_formula = '~sm.ns(Pseudotime, df=3)',
                       relative_expr = T, new_data = newdata)
  # remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  m=log10(m+1)
  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>3] = 3
  m[m<-3] = -3
  heatmap_matrix <- m
  
  # 这样我们就可以画个基本的heatmap了
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
}

# 添加行列注释
library(arm)
annotation_col = data.frame(pseudotime = newdata$Pseudotime,
                            row.names = colnames(m))
# annotation_row = data.frame(row.names = Time_genes,
#                             Cluster = ifelse(rownames(m) %in% unique(top10$gene), as.character(top10$cluster), 'unkonw'))
#table(annotation_row$Cluster)
# 作图：拟时序热图
library(scales) 
library(ggsci)
#show_col(pal_d3('category20')(20))
pal = pal_d3('category20')(20)
pal
#rowcolor <- pal[1:6]
#names(rowcolor) <- unique(annotation_row$Cluster) #类型颜色
library(viridis)
#ann_colors <- list(pseudotime=viridis(100),Cluster=rowcolor) #颜色设置
library(pheatmap)
pdf(file = "./Fig/Step7/Monocle_MLgene2.pdf", width = 8, height = 6)
pheatmap(m,
         #useRaster = T,
         cluster_cols=FALSE,
         cluster_rows=F, # 行聚类
         show_rownames=T,
         show_colnames=F,
         #clustering_method = "ward.D2",
         #clustering_distance_rows=row_dist,
         #cutree_rows=5,
         border_color = NA,
         filename=NA,
         color=colorRampPalette(c("#1F77B4FF","white","#FF7F0EFF"))(100),
         annotation_col = annotation_col
         #annotation_colors = ann_colors,
         #annotation_row = annotation_row
)
dev.off()




