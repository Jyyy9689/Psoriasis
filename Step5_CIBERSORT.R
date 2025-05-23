
rm(list = ls())
gc()

####细胞信息处理####
library(dplyr)
library(Seurat)
load(file = 'Step2_output.Rdata')
table(Idents(scRNA_harmony))
table(scRNA_harmony@meta.data[["status"]])
Healthy <- scRNA_harmony[ , scRNA_harmony@meta.data[["status"]] %in% 'Healthy']
Healthy
table(Idents(Healthy))

Healthy_marker = FindAllMarkers(Healthy, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
Healthy_top50 = Healthy_marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
table(duplicated(Healthy_top50$gene))
marker0 <- unique(Healthy_top50$gene)


options(max.print=1000000)
# 将稀疏矩阵数据重新录入matrix
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x

  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

count_list = list()
for (i in as.character(unique(Idents(Healthy)))) {
  state = Healthy[, Idents(Healthy) == i]
  state
  state_dat = state@assays[["RNA"]]@counts
  state_dat1 <- as_matrix(state_dat)
  state_dat1 = state_dat1[unique(marker0),]
  # state_dat1$ExpMean =  rowMeans(state_dat1)
  state_dat2 = data.frame(row.names = rownames(state_dat1),
                          EXPSum = rowMeans(state_dat1))
  names(state_dat2) = i
  range(state_dat2)
  count_list[[i]] = state_dat2
}

Healthy_exp = do.call(cbind,count_list)
write.table(Healthy_exp,file="Step5_Healthy_Cell_top50Marker.txt",sep = "\t")

# PSO
table(Idents(scRNA_harmony))
table(scRNA_harmony@meta.data[["status"]])
PSO <- scRNA_harmony[ , scRNA_harmony@meta.data[["status"]] %in% 'PSO']
PSO
table(Idents(PSO))

PSO_marker = FindAllMarkers(PSO, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
PSO_top50 = PSO_marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
table(duplicated(PSO_top50$gene))
marker0 <- unique(PSO_top50$gene)

count_list = list()
for (i in as.character(unique(Idents(PSO)))) {
  state = PSO[, Idents(PSO) == i]
  state
  state_dat = state@assays[["RNA"]]@counts
  state_dat1 <- as_matrix(state_dat)
  state_dat1 = state_dat1[unique(marker0),]
  # state_dat1$ExpMean =  rowMeans(state_dat1)
  state_dat2 = data.frame(row.names = rownames(state_dat1),
                          EXPSum = rowMeans(state_dat1))
  names(state_dat2) = i
  range(state_dat2)
  count_list[[i]] = state_dat2
}

PSO_exp = do.call(cbind,count_list)
write.table(PSO_exp,file="Step5_PSO_Cell_top50Marker.txt",sep = "\t")

save(Healthy, Healthy_marker, Healthy_top50, 
     PSO, PSO_marker, PSO_top50,
     file = 'Step5_output.Rdata')



####Bluk数据处理####
load(file = 'Step4_PseMerge.Rdata')
library(limma)
exp <- PSE_expr_limma
exp[1:4,1:4]

# Normal
table(PSE_phe$diagnosis)
nor <- PSE_phe[PSE_phe$diagnosis == 'normal',]
exp1 <- exp[,rownames(nor)]
dim(exp1) # [1] 14046   155
# 整理一下行名，列名，删除表达特别低的基因
data=as.matrix(exp1)
dimnames=list(rownames(exp1),colnames(exp1))
data=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
# 去除负值
data[data<0] = NA
data = na.omit(data)
is.recursive(data)
is.atomic(data)
data[1:4,1:4]
dim(data)
write.table(data,file="Step5_Normal_expr.txt",sep = "\t")

# PSE
table(PSE_phe$diagnosis)
pse <- PSE_phe[PSE_phe$diagnosis == 'Psoriasis',]
exp1 <- exp[,rownames(pse)]
dim(exp1) # [1] 14046   155
# 整理一下行名，列名，删除表达特别低的基因
data=as.matrix(exp1)
dimnames=list(rownames(exp1),colnames(exp1))
data=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
# 去除负值
data[data<0] = NA
data = na.omit(data)
is.recursive(data)
is.atomic(data)
data[1:4,1:4]
dim(data)
write.table(data,file="Step5_Pse_expr.txt",sep = "\t")



####Cibersort####
source("/home/datahup/pjy/BRCA/NK/Rpro/Cibersort.R")
library(preprocessCore)
library(e1071)
library(parallel)
Normal_result <- CIBERSORT("Step5_Healthy_Cell_top50Marker.txt", "Step5_Normal_expr.txt", perm = 1000, QN = T)
PSE_result <- CIBERSORT("Step5_PSO_Cell_top50Marker.txt", "Step5_Pse_expr.txt", perm = 1000, QN = T)

save(Normal_result,PSE_result,
     file = "Step5_CibersortResults.Rdata")


# 可视化
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
load(file = 'Step5_CibersortResults.Rdata')

# PSE
Cib_df1 <- as.data.frame(PSE_result)
Cib_df1 <- Cib_df1[,c(1:(ncol(Cib_df1)-3))]
Cib_df2 = Cib_df1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

# bar
library(ggsci)
mypalette = pal_d3("category20")(20)
mypalette1 = pal_ucscgb()(10)
ggplot(Cib_df2,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Cibersort Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = c(mypalette,mypalette1))
ggsave(filename = './Fig/Step5/PSE_Bar.pdf',width = 14,height = 8)

# box
ggplot(Cib_df2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Cibersort Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = c(mypalette,mypalette1))
ggsave(filename = './Fig/Step5/PSE_Box.pdf',width = 10,height = 8)


# normal
Cib_df3 <- as.data.frame(Normal_result)
Cib_df3 <- Cib_df3[,c(1:(ncol(Cib_df3)-3))]
Cib_df4 = Cib_df3 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

library(ggsci)
mypalette = pal_d3("category20")(20)
mypalette1 = pal_ucscgb()(10)
ggplot(Cib_df4,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Cibersort Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = c(mypalette,mypalette1))
ggsave(filename = './Fig/Step5/Nor_Bar.pdf',width = 14,height = 8)

# box
ggplot(Cib_df4,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Cibersort Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = c(mypalette,mypalette1))
ggsave(filename = './Fig/Step5/Nor_Box.pdf',width = 10,height = 8)


# PSE VS nor 
PSE_phe1 <- subset(PSE_phe, PSE_phe$Type != 'normal-normal')
table(PSE_phe1$Type)
Cib_df2$Type <- ifelse(Cib_df2$Sample %in% rownames(PSE_phe1), PSE_phe1$Type, 'No')
table(Cib_df2$Type)
Cib_df4$Type <- rep('Normal',nrow(Cib_df4))

Cib_df5 <- rbind(Cib_df2, Cib_df4)
table(Cib_df5$Type)
table(Cib_df5$Cell_type)
Cib_df5 <- subset(Cib_df5, Cib_df5$Cell_type != 'Megakaryocyte')

library(ggpubr)
mypalette2 <- pal_lancet()(9)
ggplot(Cib_df5,aes(Cell_type,Proportion,fill = Type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Cibersort  Proportion") +
  theme(legend.position = "right") + 
  theme(axis.title.y = element_text(),
        axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0),
        axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0),
        legend.text = element_text(color="black", size = 10, face = "bold")
        )+
  scale_fill_manual(values = mypalette2)+
  coord_flip() +
  stat_compare_means(aes(group = Type,label = ..p.signif..),method = "kruskal.test") # wilcox.test
ggsave(filename = './Fig/Step5/Group_Box.pdf',width = 8,height = 8)





