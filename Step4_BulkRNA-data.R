
rm(list = ls())
gc()


####GEO1:GSE186063####
library(GEOquery)
eSet <- getGEO('GSE186063', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
table(phe$`diagnosis:ch1`)
phe1 <- subset(phe, phe$`diagnosis:ch1` == 'Psoriasis')
GEO1_phe <- data.frame(row.names = rownames(phe1),
                       source = phe1$source_name_ch1,
                       age = phe1$`age:ch1`,
                       sex = phe1$`Sex:ch1`,
                       diagnosis = phe1$`diagnosis:ch1`)

# 表达信息处理
expr <- read.table(gzfile('../../IuputData/BulkRNA/GSE186063_Raw_gene_counts_matrix.txt.gz'),sep = ',',header = T, row.names = 1)
colnames(expr) <- str_split(colnames(expr), 'X', simplify = T)[,2]
colnames(expr) <- gsub('[.]','-',colnames(expr))
expr1 <- as.data.frame(t(expr))
expr1$sample <- rownames(expr1)

sam <- read.csv(file = '../../IuputData/BulkRNA/GSE186063_sample.csv',header = F) # 样本与检索号对应
names(sam) <- c('ID','sample')

expr2 <- merge(sam, expr1, by = 'sample')
rownames(expr2) <- expr2$ID
expr2[1:4,1:4]
expr2 <- expr2[,-c(1,2)]

GEO1_expr <- as.data.frame(t(expr2[rownames(GEO1_phe),]))
range(GEO1_expr)
GEO1_expr <- log2(GEO1_expr+1) # 数据缩放
range(GEO1_expr)
gene <- rownames(GEO1_expr)
write.table(gene, file = 'Step4_geneid.txt') # 服务器无法在线查找，替换笔记本

# library(biomaRt)
# library(stringr)
# gene <- read.table(file = 'Step4_geneid.txt')
# mart <- useMart("ensembl",'hsapiens_gene_ensembl')
# listFilters(mart)#查看可以输入ID的类型
# listAttributes(mart)
# # gene <- str_split(rownames(expr), '[.]', simplify = T)[,1] #需要转换的ENSN id不能小数点
# gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
#                    filters = "ensembl_gene_id",
#                    values = gene$x, mart = mart)
# write.csv(gene_name, file = 'Step4_gene_name.csv')

ids <- read.csv(file = 'Step4_gene_name.csv', row.names = 1)
names(ids) <- c('ENSG','symbol','ENTREZ')
GEO1_expr1 <- GEO1_expr
GEO1_expr1$ENSG <- rownames(GEO1_expr1)
GEO1_expr1 <- merge(ids, GEO1_expr1, by = 'ENSG')
GEO1_expr1[1:4,1:4]
GEO1_expr1 <- GEO1_expr1[!duplicated(GEO1_expr1$symbol),]
rownames(GEO1_expr1) <- GEO1_expr1$symbol
GEO1_expr1 <- GEO1_expr1[,-c(1:3)]

save(GEO1_phe, GEO1_expr1, file = 'Step4_GEOdata.Rdata')



####GEO2:GSE54456####
library(GEOquery)
eSet <- getGEO('GSE54456', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
GEO2_phe <- data.frame(row.names = rownames(phe),
                       soure = phe$`tissue type:ch1`,
                       diagnosis = phe$source_name_ch1)
table(GEO2_phe$soure)
GEO2_phe$soure <- ifelse(GEO2_phe$soure == 'lesional psoriatic skin', 'lesion', 'normal')
table(GEO2_phe$diagnosis)
GEO2_phe$diagnosis <- ifelse(GEO2_phe$diagnosis == 'Psoriasis_skin', 'Psoriasis', 'normal')

# 表达信息处理
expr <- read.table(file = '../../IuputData/BulkRNA/GSE54456_RPKM_samples.txt.gz')
range(expr)
str(expr)
expr <- log2(expr+1)
range(expr)

sam <- read.csv(file = '../../IuputData/BulkRNA/GSE54456_sample.csv',header = F)
table(sam$V2 %in% colnames(expr))
names(sam) <- c('sample','ID')
expr1 <- as.data.frame(t(expr))
expr1$ID <- rownames(expr1)
expr1 <- merge(sam,expr1,by = 'ID')
expr1[1:4,1:4]
rownames(expr1) <- expr1$sample
expr1 <- expr1[,-c(1,2)]
GEO2_expr <- as.data.frame(t(expr1))
GEO2_expr[1:4,1:4]

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     file = 'Step4_GEOdata.Rdata')



#GEO3:GSE63741
# library(GEOquery)
# eSet <- getGEO('GSE63741', destdir=".",
#                AnnotGPL = F,
#                getGPL = F)
# b = eSet[[1]]
# raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载
# 
# # 临床信息处理
# phe=pData(b)
# GEO3_phe <- data.frame(row.names = rownames(phe),
#                        sample = phe$title,
#                        soure = phe$`sample type:ch2`)
# table(GEO3_phe$soure)
# GEO3_phe <- subset(GEO3_phe, GEO3_phe$soure == 'Control' | GEO3_phe$soure == 'Psoriasis Vulgaris')
# GEO3_phe$soure <- ifelse(GEO3_phe$soure == 'Psoriasis Vulgaris', 'lesion', 'normal')
# GEO3_phe$diagnosis <- ifelse(GEO3_phe$soure == 'lesion', 'Psoriasis', 'normal')
# table(GEO3_phe$diagnosis)
# 
# # 表达信息处理
# dat <- read.delim(file = '../../IuputData/BulkRNA/GSE63741/GSM1556392_000404_070_1b.txt.gz')



####GEO4:GSE121212####
library(GEOquery)
eSet <- getGEO('GSE121212', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
GEO4_phe <- data.frame(row.names = rownames(phe),
                       sample = phe$title,
                       soure = phe$`skin type:ch1`,
                       diagnosis = phe$`patient's condition:ch1`)
table(GEO4_phe$diagnosis)
GEO4_phe <- subset(GEO4_phe, GEO4_phe$diagnosis == 'PSO' | GEO4_phe$diagnosis == 'CTRL')
GEO4_phe$diagnosis <- ifelse(GEO4_phe$diagnosis == 'PSO', 'Psoriasis', 'normal')
table(GEO4_phe$soure)
GEO4_phe$soure <- ifelse(GEO4_phe$soure == 'healthy' , 'normal', GEO4_phe$soure)


# 表达数据处理
expr <- read.delim(file = '../../IuputData/BulkRNA/GSE121212_readcount.txt.gz')
table(duplicated(expr$X))
expr <- expr[!duplicated(expr$X),]
rownames(expr) <- expr$X
colnames(expr) <- gsub('[.]', '-', colnames(expr))
expr1 <- as.data.frame(t(expr[,GEO4_phe$sample]))
expr1$ID <- rownames(expr1)

sam <- data.frame(row.names = rownames(GEO4_phe),
                  ID = GEO4_phe$sample,
                  sample = rownames(GEO4_phe))
GEO4_expr <- merge(sam, expr1, by = 'ID')
GEO4_expr[1:4,1:4]
rownames(GEO4_expr) <- GEO4_expr$sample
GEO4_expr <- as.data.frame(t(GEO4_expr[,-c(1,2)]))

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     file = 'Step4_GEOdata.Rdata')




####GEO5:GSE66511####
library(GEOquery)
eSet <- getGEO('GSE66511', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
GEO5_phe <- data.frame(row.names = phe$geo_accession,
                       sample = phe$title,
                       age = phe$`age:ch1`,
                       source = ifelse(phe$`disease:ch1` == 'LP', 'lesion',
                                       ifelse(phe$`disease:ch1` == 'NLP', 'no-lesion', 'normal')),
                       sex = phe$`gender:ch1`,
                       diagnosis = ifelse(phe$`disease:ch1` != 'C', 'Psoriasis','normal'))

# 表达信息处理
expr <- read.table(file = '../../IuputData/BulkRNA/GSE66511_Psoriasis_counts.txt.gz',header = T)
table(duplicated(expr$Symbol))
expr1 <- expr[!duplicated(expr$Symbol),]
rownames(expr1) <- expr1$Symbol
expr1 <- expr1[,-1]
range(expr1)
expr1 <- log2(expr1+1)
range(expr1)

table(colnames(expr1) %in% GEO5_phe$sample)
expr2 <- as.data.frame(t(expr1))
expr2$ID <- rownames(expr2)
sam <- data.frame(row.names = rownames(GEO5_phe),
                  sample = rownames(GEO5_phe),
                  ID = GEO5_phe$sample)
GEO5_expr <- merge(sam, expr2, by = 'ID')
rownames(GEO5_expr) <- GEO5_expr$sample
GEO5_expr <- GEO5_expr[,-c(1,2)]
GEO5_expr <- as.data.frame(t(GEO5_expr))
GEO5_expr[1:4,1:4]

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     file = 'Step4_GEOdata.Rdata')



####GEO6:GSE117405####
library(GEOquery)
eSet <- getGEO('GSE117405', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
table(phe$source_name_ch1)
table(phe$`disease state:ch1`)
GEO6_phe <- data.frame(row.names = phe$geo_accession,
                       ID = phe$title,
                       source = ifelse(phe$source_name_ch1 == 'Control skin', 'normal', 'lesion'),
                       diagnosis = ifelse(phe$`disease state:ch1` == 'healthy control', 'normal', 'Psoriasis'))

# 表达信息处理
expr <- read.table(file = '../../IuputData/BulkRNA/GSE117405_SCALP_HF_P_KT_FPKM_Table_05_18_2017.txt.gz', header = T)
table(duplicated(expr$Gene))
rownames(expr) <- expr$Gene
expr <- expr[,-1]
range(expr)
expr <- log2(expr+1)
range(expr)

table(colnames(expr) %in% GEO6_phe$ID) # 临床信息不全，前往GSE117405补全
setdiff(colnames(expr) , GEO6_phe$ID) # "HF_3" "SP_5" "SP_8" "SP_6" "SP_7" "CP_8" "CP_5" "CP_6" "CP_7"
GEO6_phe1 <- data.frame(row.names = c('GSM3293899','GSM3293921','GSM3293924','GSM3293922', 'GSM3293923', 'GSM3293916','GSM3293913', 'GSM3293914','GSM3293915'),
                        ID = c('HF_3', 'SP_5', 'SP_8','SP_6','SP_7','CP_8','CP_5','CP_6','CP_7'),
                        source = c('lesion', 'lesion', 'lesion', 'lesion','lesion', 'lesion', 'lesion','lesion','lesion'),
                        diagnosis = c('Psoriasis', 'Psoriasis', 'Psoriasis', 'Psoriasis', 'Psoriasis', 'Psoriasis','Psoriasis','Psoriasis','Psoriasis'))
GEO6_phe2 <- rbind(GEO6_phe, GEO6_phe1)

table(colnames(expr) %in% GEO6_phe2$ID)
sam <- data.frame(row.names = rownames(GEO6_phe2),
                  ID = GEO6_phe2$ID,
                  sample = rownames(GEO6_phe2))
expr1 <- as.data.frame(t(expr))
expr1$ID <- rownames(expr1)
expr1 <- merge(sam, expr1, by = 'ID')
expr1[1:4,1:4]
rownames(expr1) <- expr1$sample
expr1 <- expr1[,-c(1,2)]
GEO6_expr <- as.data.frame(t(expr1))
GEO6_expr[1:4,1:4]

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     GEO6_phe2,GEO6_expr,
     file = 'Step4_GEOdata.Rdata')



####GEO7:GSE205748####
library(GEOquery)
eSet <- getGEO('GSE205748', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
table(phe$`tissue type:ch1`)
GEO7_phe <- data.frame(row.names = phe$geo_accession,
                       ID = phe$title,
                       souce = ifelse(phe$`tissue type:ch1` == 'Healthy control skin', 'normal',
                                      ifelse(phe$`tissue type:ch1` == 'Psoriatic arthritis skin lesion', 'lesion', 'no-lesion')),
                       diagnosis = ifelse(phe$`tissue type:ch1` == 'Healthy control skin', 'normal', 'Psoriasis'))

# 表达信息处理
expr <- read.delim(file = '../../IuputData/BulkRNA/GSE205748_read_counts.csv.gz')
table(duplicated(expr$ID))
rownames(expr) <- expr$ID
expr <- expr[,-1]
range(expr)
expr <- log2(expr+1)
range(expr)

library(stringr)
sam <- data.frame(row.names = rownames(GEO7_phe),
                  ID = str_split(GEO7_phe$ID, '\\[', simplify = T)[,2],
                  sample = rownames(GEO7_phe))
sam$ID <- str_split(sam$ID, '\\]', simplify = T)[,1]
table(sam$ID %in% colnames(expr))

expr1 <- as.data.frame(t(expr))
expr1$ID <- rownames(expr1)
expr1 <- merge(sam, expr1, by = 'ID')
expr1[1:4,1:4]
rownames(expr1) <- expr1$sample
expr1 <- expr1[,-c(1,2)]
GEO7_expr <- as.data.frame(t(expr1))
GEO7_expr[1:4,1:4]

gene <- rownames(GEO7_expr)
write.table(gene, file = 'Step4_geneid2.txt') # 服务器无法转换ID，换台式机
# library(biomaRt)
# gene <- read.table(file = 'Step4_geneid2.txt')
# mart <- useMart("ensembl",'hsapiens_gene_ensembl')
# listFilters(mart)#查看可以输入ID的类型
# listAttributes(mart)
# # gene <- str_split(rownames(expr), '[.]', simplify = T)[,1] #需要转换的ENSN id不能小数点
# gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
#                    filters = "ensembl_gene_id",
#                    values = gene$x, mart = mart)
# write.csv(gene_name, file = 'Step4_gene_name2.csv')
ids <- read.csv(file = './Step4_gene_name2.csv', row.names = 1)
names(ids) <- c('ENSG', 'symble', 'entrez')
GEO7_expr$ENSG <- rownames(GEO7_expr)
GEO7_expr1 <- merge(ids, GEO7_expr, by = 'ENSG')
GEO7_expr1[1:4,1:4]
GEO7_expr1 <- GEO7_expr1[!duplicated(GEO7_expr1$symble),]
rownames(GEO7_expr1) <- GEO7_expr1$symble
GEO7_expr1 <- GEO7_expr1[,-c(1:3)]

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     GEO6_phe2,GEO6_expr,
     GEO7_phe, GEO7_expr1,
     file = 'Step4_GEOdata.Rdata')



####GEO8:GSE114286####
library(GEOquery)
eSet <- getGEO('GSE114286', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
table(phe$`participant condition:ch1`)
GEO8_phe <- data.frame(row.names = phe$geo_accession,
                       ID = phe$title,
                       source = ifelse(phe$`participant condition:ch1` == 'healthy volunteer', 'normal', 'lesion'),
                       diagnosis = ifelse(phe$`participant condition:ch1` == 'healthy volunteer', 'normal', 'Psoriasis'),
                       age = phe$`age (years):ch1`,
                       sex = phe$`gender:ch1`)

# 表达数据处理
expr <- read.table(file = '../../IuputData/BulkRNA/GSE114286_All_Counts.txt.gz', header = T)
table(duplicated(expr$GeneID))
rownames(expr) <- expr$GeneID
expr <- expr[,-1]
range(expr)
expr <- log2(expr+1)
range(expr)

table(colnames(expr) %in% GEO8_phe$ID)
expr1 <- as.data.frame(t(expr))
expr1$ID <- rownames(expr1)
sam <- data.frame(row.names = rownames(GEO8_phe),
                  ID = GEO8_phe$ID,
                  sample = rownames(GEO8_phe))
expr2 <- merge(sam, expr1, by = 'ID')
expr2[1:4,1:4]
rownames(expr2) <- expr2$sample
expr2 <- expr2[,-c(1:2)]
GEO8_expr <- as.data.frame(t(expr2))
GEO8_expr[1:4,1:4]

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     GEO6_phe2,GEO6_expr,
     GEO7_phe, GEO7_expr1,
     GEO8_phe, GEO8_expr,
     file = 'Step4_GEOdata.Rdata')



####GEO9:GSE83645####
library(GEOquery)
eSet <- getGEO('GSE83645', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
table(phe$`tissue type:ch1`)
GEO9_phe <- data.frame(row.names = phe$geo_accession,
                       ID = phe$title,
                       source = ifelse(phe$`tissue type:ch1` == 'psoriasis', 'lesion', 'no-lesion'),
                       diagnosis = rep('Psoriasis', nrow(phe)))

# 表达信息处理
expr <- read.table(file = '../../IuputData/BulkRNA/GSE83645_rawdatacount.txt.gz')
range(expr)
expr <- log2(expr+1)
range(expr)

colnames(expr) <- str_split(colnames(expr), '_',simplify = T)[,2]
table(colnames(expr) %in% GEO9_phe$ID)
sam <- data.frame(row.names = rownames(GEO9_phe),
                  ID = GEO9_phe$ID,
                  sample = rownames(GEO9_phe))
expr1 <- as.data.frame(t(expr))
expr1$ID <- rownames(expr1)
expr2 <- merge(sam, expr1, by = 'ID')
expr2[1:4,1:4]
rownames(expr2) <- expr2$sample
expr2 <- expr2[,-c(1,2)]
GEO9_expr <- as.data.frame(t(expr2))
GEO9_expr[1:4,1:4]

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     GEO6_phe2,GEO6_expr,
     GEO7_phe, GEO7_expr1,
     GEO8_phe, GEO8_expr,
     GEO9_phe, GEO9_expr,
     file = 'Step4_GEOdata.Rdata')



####GEO10:GSE142582####
library(GEOquery)
eSet <- getGEO('GSE142582', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
table(phe$`tissue subtype:ch1`)
table(phe$`subject status:ch1`)
GEO10_phe <- data.frame(row.names = phe$geo_accession,
                        ID = phe$title,
                        age = phe$`age:ch1`,
                        source = ifelse(phe$`tissue subtype:ch1` == 'Psoriatic plaque tissue (PS)', 'lesion',
                                        ifelse(phe$`tissue subtype:ch1` == 'adjacent normal skin tissue (PN)', 'no-lesion', 'normal')),
                        diagnosis = ifelse(phe$`subject status:ch1` == 'healthy control', 'normal', 'Psoriasis'))

# 表达信息处理
expr <- read.delim(file = '../../IuputData/BulkRNA/GSE142582_mRNA_counts.txt.gz')
expr <- expr[,1:17]
table(duplicated(expr$gene))
expr1 <- expr[!duplicated(expr$gene),]
rownames(expr1) <- expr1$gene
expr1 <- expr1[,-c(1,17)]
range(expr1)
expr1 <- log2(expr1+1)
range(expr1)

sam <- data.frame(row.names = rownames(GEO10_phe),
                  ID = str_split(GEO10_phe$ID, '\\[', simplify = T)[,2],
                  sample = rownames(GEO10_phe))
sam$ID <- str_split(sam$ID, '\\]', simplify = T)[,1]
table(colnames(expr1) %in% sam$ID)
expr2 <- as.data.frame(t(expr1))
expr2$ID <- rownames(expr2)
expr2 <- merge(sam, expr2, by = 'ID')
expr2[1:4,1:4]
rownames(expr2) <- expr2$sample
expr2 <- expr2[,-c(1,2)]
GEO10_expr <- as.data.frame(t(expr2))
GEO10_expr[1:4,1:4]



save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     GEO6_phe2,GEO6_expr,
     GEO7_phe, GEO7_expr1,
     GEO8_phe, GEO8_expr,
     GEO9_phe, GEO9_expr,
     GEO10_phe, GEO10_expr,
     file = 'Step4_GEOdata.Rdata')



####GEO11:GSE41745####
library(GEOquery)
eSet <- getGEO('GSE41745', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
table(phe$`group:ch1`)
GEO11_phe <- data.frame(row.names = phe$geo_accession,
                        ID = phe$title,
                        source = ifelse(phe$`group:ch1` == 'lesional', 'lesion', 'no-lesion'),
                        diagnosis = rep('Psoriasis', nrow(phe)))

# 表达信息处理
expr_list <- list()
for (i in paste0(rownames(GEO11_phe), '_',GEO11_phe$ID)) {
  expr <- read.delim(file = paste0('../../IuputData/BulkRNA/GSE41745/',i,'.counts.v2.txt.gz'), header = F)
  expr <- expr[-c(33656:33659),]
  names(expr) <- c('gene', i)
  expr_list[[i]] <- expr
}
expr1 <- do.call(cbind, expr_list)
expr1 <- expr1[,c(1,2,4,6,8,10,12)]
table(duplicated(expr1$GSM1023434_NL251.gene))
colnames(expr1) <- str_split(colnames(expr1),'[.]',simplify = T)[,2]
colnames(expr1)[2:7] <- str_split(colnames(expr1)[2:7],'_',simplify = T)[,1]

write.table(gene, file = 'Step4_geneid3.txt') # 去笔记本上做ID转换
# library(biomaRt)
# gene <- read.table(file = 'Step4_geneid3.txt')
# mart <- useMart("ensembl",'hsapiens_gene_ensembl')
# listFilters(mart)#查看可以输入ID的类型
# listAttributes(mart)
# # gene <- str_split(rownames(expr), '[.]', simplify = T)[,1] #需要转换的ENSN id不能小数点
# gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "refseq_mrna"),
#                    filters = "refseq_mrna",
#                    values = gene$x, mart = mart)
# write.csv(gene_name, file = 'Step4_gene_name3.csv')
ids <- read.csv(file = 'Step4_gene_name3.csv')
table(expr1$gene %in% ids$refseq_mrna)
expr2 <- expr1[expr1$gene %in% ids$refseq_mrna,]
ids1 <- ids[,3:4]
names(expr2)[1] <- 'refseq_mrna'
expr3 <- merge(ids1, expr2, by = 'refseq_mrna')
table(duplicated(expr3$external_gene_name))
expr3 <- expr3[!duplicated(expr3$external_gene_name),]
rownames(expr3) <- expr3$external_gene_name
expr3 <- expr3[,-c(1:2)]
range(expr3)
expr3 <- log2(expr3+1)
range(expr3)
GEO11_expr <- expr3

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     GEO6_phe2,GEO6_expr,
     GEO7_phe, GEO7_expr1,
     GEO8_phe, GEO8_expr,
     GEO9_phe, GEO9_expr,
     GEO10_phe, GEO10_expr,
     GEO11_phe, GEO11_expr,
     file = 'Step4_GEOdata.Rdata')



####GEO12:GSE186117####
library(GEOquery)
eSet <- getGEO('GSE186117', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
raw_exprSet= exprs(b)  # 表达谱为空，需要单独下载

# 临床信息处理
phe=pData(b)
GEO12_phe <- data.frame(row.names = phe$geo_accession,
                    ID = str_split(phe$title, '\\[', simplify = T)[,2],
                    source = ifelse(str_split(phe$`volunteer:ch1`, ' ', simplify = T)[,4] == 'lesion', 'lesion', 'no-lesion'),
                    diagnosis = rep('Psoriasis', nrow(phe)))
GEO12_phe$ID <- str_split(GEO12_phe$ID, '\\]', simplify = T)[,1]

# 表达信息处理
library(openxlsx)
expr <- read.xlsx('../../IuputData/BulkRNA/GSE186117_Lesional_skin_Expression_Profile.hg19.transcript.xlsx')
expr <- expr[,c(3,7:12)]
table(duplicated(expr$Gene_Symbol))
expr <- expr[!duplicated(expr$Gene_Symbol),]
rownames(expr) <- expr$Gene_Symbol
expr <- expr[,-1]

expr1 <- read.xlsx('../../IuputData/BulkRNA/GSE186117_Uninvovled_Expression_Profile.hg19.transcript.xlsx')
expr1 <- expr1[,c(3,7:10)]
table(duplicated(expr1$Gene_Symbol))
expr1 <- expr1[!duplicated(expr1$Gene_Symbol),]
rownames(expr1) <- expr1$Gene_Symbol
expr1 <- expr1[,-1]

expr2 <- cbind(expr[,1:3], expr1[,1:2])
colnames(expr2) <- str_split(colnames(expr2), '_', simplify = T)[,1]
colnames(expr2) <- ifelse(colnames(expr2) %in% GEO12_phe$ID, rownames(GEO12_phe), 'No')

range(expr2)
expr2 <- log2(expr2+1)
range(expr2)
GEO12_expr <- expr2

save(GEO1_phe, GEO1_expr1,
     GEO2_phe, GEO2_expr,
     GEO4_phe, GEO4_expr,
     GEO5_phe, GEO5_expr,
     GEO6_phe2,GEO6_expr,
     GEO7_phe, GEO7_expr1,
     GEO8_phe, GEO8_expr,
     GEO9_phe, GEO9_expr,
     GEO10_phe, GEO10_expr,
     GEO11_phe, GEO11_expr,
     GEO12_phe, GEO12_expr,
     file = 'Step4_GEOdata.Rdata')


####Merge####
rm(list = ls())
gc()
load(file = 'Step4_GEOdata.Rdata')

# 临床信息合并
GEO1_phe1 <- GEO1_phe[,c(1,4)]
GEO1_phe1$ID <- rownames(GEO1_phe1)
GEO1_phe1$GSE <- rep('GSE186063',nrow(GEO1_phe1))
names(GEO1_phe1)

GEO2_phe$ID <- rownames(GEO2_phe)
GEO2_phe$GSE <- rep('GSE54456',nrow(GEO2_phe))
names(GEO2_phe)[1] <- 'source'
names(GEO2_phe)

names(GEO4_phe)[1] <- 'ID'
GEO4_phe$GSE <- rep('GSE121212',nrow(GEO4_phe))
GEO4_phe <- GEO4_phe[,c(2,3,1,4)]
names(GEO4_phe)[1] <- 'source'
names(GEO4_phe)

GEO5_phe1 <- GEO5_phe[,c(1,3,5)]
names(GEO5_phe1)[1] <- 'ID'
GEO5_phe1$GSE <- rep('GSE66511',nrow(GEO5_phe1))
GEO5_phe1 <- GEO5_phe1[,c(2,3,1,4)]
names(GEO5_phe1)

GEO6_phe2$GSE <- rep('GSE117405',nrow(GEO6_phe2))
GEO6_phe2 <- GEO6_phe2[,c(2,3,1,4)]
names(GEO6_phe2)

GEO7_phe$GSE <- rep('GSE205748',nrow(GEO7_phe))
GEO7_phe <- GEO7_phe[,c(2,3,1,4)]
names(GEO7_phe)[1] <- 'source'
names(GEO7_phe)

GEO8_phe1 <- GEO8_phe[,c(1,2,3)]
GEO8_phe1$GSE <- rep('GSE114286',nrow(GEO8_phe1))
GEO8_phe1 <- GEO8_phe1[,c(2,3,1,4)]
names(GEO8_phe1)

GEO9_phe$GSE <- rep('GSE83645',nrow(GEO9_phe))
GEO9_phe <- GEO9_phe[,c(2,3,1,4)]
names(GEO9_phe)

GEO10_phe1 <- GEO10_phe[,c(1,3,4)]
GEO10_phe1$GSE <- rep('GSE142582',nrow(GEO10_phe1))
GEO10_phe1 <- GEO10_phe1[,c(2,3,1,4)]
names(GEO10_phe1)

GEO11_phe$GSE <- rep('GSE41745',nrow(GEO11_phe))
GEO11_phe <- GEO11_phe[,c(2,3,1,4)]
names(GEO11_phe)

GEO12_phe$GSE <- rep('GSE186117',nrow(GEO12_phe))
GEO12_phe <- GEO12_phe[,c(2,3,1,4)]
names(GEO12_phe)

PSE_phe <- rbind(GEO1_phe1, GEO2_phe, GEO4_phe, GEO5_phe1, GEO6_phe2,
                 GEO7_phe, GEO8_phe1, GEO9_phe, GEO10_phe1, GEO11_phe, GEO12_phe)
table(PSE_phe$source)
PSE_phe[PSE_phe=='lesional'] <- 'lesion'
PSE_phe[PSE_phe=='non-lesion'] <- 'no-lesion'
PSE_phe[PSE_phe=='non-lesional'] <- 'no-lesion'
table(PSE_phe$source)
table(PSE_phe$diagnosis)
table(PSE_phe$GSE)

# 表达信息合并
gene_int <- Reduce(function(x,y) intersect(x,y), list(rownames(GEO1_expr1), rownames(GEO2_expr), rownames(GEO4_expr),
                                                      rownames(GEO5_expr), rownames(GEO6_expr), rownames(GEO7_expr1),
                                                      rownames(GEO8_expr), rownames(GEO9_expr), rownames(GEO10_expr),
                                                      rownames(GEO11_expr), rownames(GEO12_expr)),
                   accumulate =FALSE)
length(unique(gene_int))

GEO1_expr1[1:4,1:4]
range(GEO1_expr1)
GEO1_expr2 <- GEO1_expr1[gene_int,]
table(is.na(GEO1_expr2))
GEO1_expr2[1:4,1:4]

GEO2_expr[1:4,1:4]
range(GEO2_expr)
GEO2_expr1 <- GEO2_expr[gene_int,]
table(is.na(GEO2_expr1))
GEO2_expr1[1:4,1:4]

GEO4_expr[1:4,1:4]
range(GEO4_expr)
GEO4_expr1 <- log2(GEO4_expr+1)
range(GEO4_expr1)
GEO4_expr1[1:4,1:4]
GEO4_expr1 <- GEO4_expr1[gene_int,]
table(is.na(GEO4_expr1))

GEO5_expr[1:4,1:4]
range(GEO5_expr)
GEO5_expr1 <- GEO5_expr[gene_int,]
GEO5_expr1[1:4,1:4]
table(is.na(GEO5_expr1))

GEO6_expr[1:4,1:4]
range(GEO6_expr)
GEO6_expr1 <- GEO6_expr[gene_int,]
GEO6_expr1[1:4,1:4]
table(is.na(GEO6_expr1))

GEO7_expr1[1:4,1:4]
range(GEO7_expr1)
GEO7_expr2 <- GEO7_expr1[gene_int,]
GEO7_expr2[1:4,1:4]
table(is.na(GEO7_expr2))

GEO8_expr[1:4,1:4]
range(GEO8_expr)
GEO8_expr1 <- GEO8_expr[gene_int,]
GEO8_expr1[1:4,1:4]
table(is.na(GEO8_expr1))

GEO9_expr[1:4,1:4]
range(GEO9_expr)
GEO9_expr1 <- GEO9_expr[gene_int,]
GEO9_expr1[1:4,1:4]
table(is.na(GEO9_expr1))

GEO10_expr[1:4,1:4]
range(GEO10_expr)
GEO10_expr1 <- GEO10_expr[gene_int,]
GEO10_expr1[1:4,1:4]
table(is.na(GEO10_expr1))

GEO11_expr[1:4,1:4]
range(GEO11_expr)
GEO11_expr1 <- GEO11_expr[gene_int,]
GEO11_expr1[1:4,1:4]
table(is.na(GEO11_expr1))

GEO12_expr[1:4,1:4]
range(GEO12_expr)
GEO12_expr1 <- GEO12_expr[gene_int,]
GEO12_expr1[1:4,1:4]
table(is.na(GEO12_expr1))

#合并
PSE_expr <- cbind(GEO1_expr2, GEO2_expr1, GEO4_expr1, GEO5_expr1,
                  GEO6_expr1, GEO7_expr2, GEO8_expr1, GEO9_expr1,
                  GEO10_expr1, GEO11_expr1, GEO12_expr1)
dim(PSE_expr)
PSE_expr[1:4,1:4] # 行为基因，列为样本

# 批次校正,两种方法尝试后选最佳
library(bladderbatch)
PSE_phe$Type <- paste0(PSE_phe$diagnosis, '-', PSE_phe$source)
table(PSE_phe$Type)
model <- model.matrix(~as.factor(PSE_phe$Type)) #建立批次效应的模型，data$type表示的是数据中除了有不同的批次，还有生物学上的差异。

library(sva)
PSE_expr_combat <- ComBat(dat = PSE_expr, batch = PSE_phe$GSE, mod = model)

library(limma)
PSE_expr_limma <- removeBatchEffect(PSE_expr, batch = PSE_phe$GSE, design=model)

# PCA
library(factoextra) 
library(FactoMineR)
shape_level <- length(unique(PSE_phe$GSE)) # type 需要改成自己映射到形状的列名
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}

# 批次校正前
PSE_expr_pca <- PCA(t(PSE_expr), graph = FALSE)  # 行为样本 列为基因
fviz_pca_ind(PSE_expr_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = PSE_phe$GSE, # color by groups
             #palette = mypalette[c(1,3)],
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
) + labs(title = 'Batch Effect') + scale_shape_manual(values=shapes)
ggsave(filename = './Fig/Step4/BatchEffect_Before.pdf', width = 8, height = 6)

# 批次校正后
PSE_expr_combat_pca <- PCA(t(PSE_expr_combat), graph = FALSE)
fviz_pca_ind(PSE_expr_combat_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = PSE_phe$GSE, # color by groups
             #palette = mypalette[c(1,3)],
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
) + labs(title = 'Remove Batch Effect') + scale_shape_manual(values=shapes)
ggsave(filename = './Fig/Step4/BatchEffect_After_combat.pdf', width = 8, height = 6)

PSE_expr_limma_pca <- PCA(t(PSE_expr_limma), graph = FALSE)
fviz_pca_ind(PSE_expr_limma_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = PSE_phe$GSE, # color by groups
             #palette = mypalette[c(1,3)],
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
) + labs(title = 'Remove Batch Effect') + scale_shape_manual(values=shapes)
ggsave(filename = './Fig/Step4/BatchEffect_After_limma.pdf', width = 8, height = 6)


# 保存最佳结果
PSE_expr_limma[1:4,1:4]
PSE_expr1 <- as.data.frame(PSE_expr_limma)
save(PSE_phe, PSE_expr_limma, file = 'Step4_PseMerge.Rdata')







