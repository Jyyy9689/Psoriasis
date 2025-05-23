# 2023.9.25-2123.9.27
# step1:scRNA数据清洗
rm(list = ls())
gc()

library(Seurat)
library(stringr)


# GSE194315
phe <- read.delim(gzfile('/home/datahup/pjy/PSE/IuputData/scRNA/GSE194315/GSE194315_CellMetadata-PSA_TotalCiteseq_20220103.tsv.gz'))
table(phe$Status)
phe_healthy <- subset(phe, phe$Status == 'Healthy')
table(phe_healthy$Sample)
phe_pso <- subset(phe, phe$Status == 'PSO')
table(phe_pso$Sample)

sample1 <- unique(c(unique(phe_healthy$Sample), unique(phe_pso$Sample)))
pse_list <- list()
for (i in sample1) {
  expr <- Read10X(paste0('/home/datahup/pjy/PSE/IuputData/scRNA/GSE194315/sample/',i))
  if (is.list(expr) == F) {
    pse1 <- CreateSeuratObject(counts = expr, project = i,min.cells = 10, min.features = 300)
    pse_list[[i]] <- pse1
    
  } else{
    pse1 <- CreateSeuratObject(counts = expr[["Gene Expression"]], project = i,min.cells = 10, min.features = 300)
    pse_list[[i]] <- pse1
    
  }
}

pse_merge <- merge(pse_list[['PBMC-01-3']], y = c(pse_list[['PBMC-01-4']], pse_list[['PBMC-02-1']], pse_list[['PBMC-02-2']], 
                                                  pse_list[['PBMC-02-3']], pse_list[['PBMC-02-4']], pse_list[['PBMC-03-1']],
                                                  pse_list[['PBMC-03-2']], pse_list[['PBMC-03-3']], pse_list[['PBMC-03-4']],
                                                  pse_list[['PBMC-04-1']], pse_list[['PBMC-04-2']], pse_list[['PBMC-04-3']],
                                                  pse_list[['PBMC-04-4']], pse_list[['PBMC-05-1']], pse_list[['PBMC-05-2']],
                                                  pse_list[['PBMC-05-3']], pse_list[['PBMC-05-4']], pse_list[['PBMC-06-1']],
                                                  pse_list[['PBMC-06-2']], pse_list[['PBMC-06-3']], pse_list[['PBMC-06-4']],
                                                  pse_list[['PBMC-07-1']], pse_list[['PBMC-07-2']], pse_list[['PBMC-07-3']],
                                                  pse_list[['PBMC-07-4']]),
                   add.cell.ids = sample1, 
                   project = 'GSE194315')

pse_cell <- as.data.frame(pse_merge@assays[["RNA"]]@counts@Dimnames[[2]])
names(pse_cell) <- 'Firstname'
pse_cell$Sendname <- str_sub(pse_cell$Firstname,1,-3)
#pse_phe <- rbind(phe_pso,phe_healthy)
#names(pse_phe)[1] <- 'Sendname'
names(phe)[1] <- 'Sendname'
pse_cell1 <- merge(pse_cell, phe, by = 'Sendname')
table(pse_cell1$Status)
rownames(pse_cell1) <- pse_cell1$Firstname
pse_cell1 <- pse_cell1[,-c(1,2)]

pse_merge <- AddMetaData(pse_merge, pse_cell1)
table(pse_merge@meta.data[["Status"]])
pse_merge1 <- pse_merge[ , pse_merge@meta.data[["Status"]] %in% c('Healthy', 'PSO')]
pse_merge1

save(pse_merge1, pse_cell1, file = 'Step1_GSE194315.Rdata')



# GSE151177
rm(list = ls())
gc()
sample <- read.delim(file = '/home/datahup/pjy/PSE/IuputData/scRNA/GSE151177/sample.txt',header = F)
sample <- sample$V1

pse_list2 <- list()
for (i in sample) {
  expr0 <- Read10X(paste0('/home/datahup/pjy/PSE/IuputData/scRNA/GSE151177/sample/',i))
  pse0 <- CreateSeuratObject(counts = expr0, project = i,min.cells = 10, min.features = 300)
  pse_list2[[i]] <- pse0
}

pse_merge2 <- merge(pse_list2[["GSM4567877_Control01"]], y = c(pse_list2[["GSM4567878_Control02"]], pse_list2[["GSM4567879_Control03"]],
                                                               pse_list2[["GSM4567880_Control04"]], pse_list2[["GSM4567881_Control05"]],
                                                               pse_list2[["GSM4567882_Control05F"]], pse_list2[["GSM4567883_Psoriasis01"]],
                                                               pse_list2[["GSM4567883_Psoriasis01"]], pse_list2[["GSM4567884_Psoriasis02"]],
                                                               pse_list2[["GSM4567885_Psoriasis02F"]], pse_list2[["GSM4567886_Psoriasis03"]],
                                                               pse_list2[["GSM4567887_Psoriasis03F"]], pse_list2[["GSM4567888_Psoriasis04"]],
                                                               pse_list2[["GSM4567889_Psoriasis04F"]], pse_list2[["GSM4567890_Psoriasis05"]],
                                                               pse_list2[["GSM4567891_Psoriasis06"]], pse_list2[["GSM4567892_Psoriasis06F"]],
                                                               pse_list2[["GSM4567893_Psoriasis07"]], pse_list2[["GSM4567894_Psoriasis08"]],
                                                               pse_list2[["GSM4567895_Psoriasis09"]], pse_list2[["GSM4567896_Psoriasis10"]],
                                                               pse_list2[["GSM4567897_Psoriasis11"]], pse_list2[["GSM4567898_Psoriasis12"]],
                                                               pse_list2[["GSM4567899_Psoriasis13"]]), 
                    #add.cell.ids = sample,
                    project = 'GSE151177'
)

pse2_cell <- as.data.frame(Idents(pse_merge2))
pse2_cell$patien <- str_sub(pse2_cell$`Idents(pse_merge2)`, 12, -3)
table(pse2_cell$patien)
pse2_cell[pse2_cell == 'Control0'] <- 'Control'
pse2_cell[pse2_cell == 'Psoriasis0'] <- 'Psoriasis'
table(pse2_cell$patien)
names(pse2_cell)[1] <- 'sample'
pse_merge2 <- AddMetaData(pse_merge2, pse2_cell)

save(pse_merge2, file = 'Step1_GSE151177.Rdata')  




# GSE221648
rm(list = ls())
gc()
sample <- read.delim(file = '/home/datahup/pjy/PSE/IuputData/scRNA/GSE221648/sample.txt',header = F)
sample <- sample$V1

pse_list3 <- list()
for (i in sample) {
  expr0 <- Read10X(paste0('/home/datahup/pjy/PSE/IuputData/scRNA/GSE221648/sample/',i))
  pse0 <- CreateSeuratObject(counts = expr0, project = i,min.cells = 10, min.features = 300)
  pse_list3[[i]] <- pse0
}

pse_merge3 <- merge(pse_list3[['GSM6892082']], y = c(pse_list3[['GSM6892083']], pse_list3[['GSM6892084']],
                                                     pse_list3[['GSM7528393']], pse_list3[['GSM7528394']],
                                                     pse_list3[['GSM7528395']], pse_list3[['GSM7528396']],
                                                     pse_list3[['GSM7528397']]),
                    add.cell.ids = sample,
                    project = 'GSE221648')
pse_merge3

save(pse_merge3, file = 'Step1_GSE221648.Rdata')


# GSE162183
rm(list = ls())
gc()
expr0 <- read.delim(file = gzfile('/home/datahup/pjy/PSE/IuputData/scRNA/GSE162183_Raw_gene_counts_matrix.tab.gz'))
expr0[1:4,1:4]
rownames(expr0) <- expr0[,1]
expr0 <- expr0[,-1]
pse0 <- CreateSeuratObject(counts = expr0, project = 'GSE162183',min.cells = 10, min.features = 300)
pse0
pse0@meta.data$patien <- Idents(pse0)

save(pse0,file = 'Step1_GSE162183.Rdata')


## 整合四个数据集
rm(list = ls())
gc()
load(file = 'Step1_GSE151177.Rdata')
load(file = 'Step1_GSE162183.Rdata')
load(file = 'Step1_GSE194315.Rdata')
load(file = 'Step1_GSE221648.Rdata')
pse_merge_list <- list('GSE194315' = pse_merge1,
                       'GSE151177' = pse_merge2,
                       'GSE221648' = pse_merge3,
                       'GSE162183' = pse0)
pse_merge_list2 <- merge(pse_merge_list[['GSE194315']], y = c(pse_merge_list[['GSE151177']], 
                                                              pse_merge_list[['GSE221648']], 
                                                              pse_merge_list[['GSE162183']]),
                         project = 'PSE')
pse_merge_list2

cellgroup <- as.data.frame(pse_merge_list2@active.ident)
names(cellgroup)[1] <- 'patien'
cellgroup$status <- pse_merge_list2@meta.data[["Status"]]
table(cellgroup$status)

cellgroup1 <- subset(cellgroup, is.na(cellgroup$status))
cellgroup1$patien <- as.character(cellgroup1$patien)
str(cellgroup1)
table(cellgroup1$patien)
cellgroup1$status <- ifelse(cellgroup1$patien == 'Ctrl1' | cellgroup1$patien == 'Ctrl2' | cellgroup1$patien == 'Ctrl3' | 
                              cellgroup1$patien == 'GSM4567877_Control01' | cellgroup1$patien == 'GSM4567878_Control02' |
                              cellgroup1$patien == 'GSM4567879_Control03' | cellgroup1$patien == 'GSM4567880_Control04' | 
                              cellgroup1$patien == 'GSM4567881_Control05' | cellgroup1$patien == 'GSM4567882_Control05F', 
                            'Healthy', 'PSO')
table(cellgroup1$status)

cellgroup2 <- subset(cellgroup, !(is.na(cellgroup$status)))
table(cellgroup2$status)
cellgroup3 <- rbind(cellgroup1, cellgroup2)
table(cellgroup3$status)
table(cellgroup3$patien)
pse_merge_list2 <- AddMetaData(pse_merge_list2, cellgroup3)


save(pse_merge_list2, file = 'Step1_Finaldata.Rdata')



