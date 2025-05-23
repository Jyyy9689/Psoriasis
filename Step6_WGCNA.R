
setwd('/home/datahup/pjy/PSE/Analysis/Rdata/')
rm(list = ls())
gc()

load(file = 'Step4_PseMerge.Rdata')
load(file = 'Step5_CibersortResults.Rdata')
library(WGCNA)
table(PSE_phe$diagnosis)
PSE_phe1 <- subset(PSE_phe, PSE_phe$diagnosis == 'Psoriasis')
expr <- as.data.frame(PSE_expr_limma[,rownames(PSE_phe1)])
expr[1:4,1:4]

## 检查，要求表达值都是数值类型
{
  datExpr0 = as.data.frame(t(expr))
  datExpr0 = apply(datExpr0,2,as.numeric)
  datExpr0 <- as.data.frame(datExpr0)
  rownames(datExpr0) <- colnames(expr)
  
  gsg = goodSamplesGenes(datExpr0, verbose = 3)
  gsg$allOK
}
##  去除不好的gene
dataExpr = datExpr0
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
dim(datExpr0)
head(dataExpr)[,1:8]


## 查看是否有离群样品
datExpr0 = dataExpr
sampleTree = hclust(dist(datExpr0), method = "average")
#窗口大小
sizeGrWindow(12,9)
#聚类树
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# cut线
abline(h = 200, col = "red")
#选择cut线以下的样本
clust = cutreeStatic(sampleTree, cutHeight = (200), minSize = 10)
table(clust)
##clust 1 是我们想留下的
keepSamples = (clust== 1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##移除不需要的临床数据

allTraits = PSE_phe1
##形成一个类似于表达数据的数据框，该数据框将保留临床特征。
femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, rownames(allTraits))
datTraits = allTraits[traitRows, ]
collectGarbage()

##Re-cluster samples
# sampleTree2 = hclust(dist(datExpr), method = "average")
# ##将特征转换为颜色表示：白色表示低，红色表示高，灰色表示缺少输入
# traitColors = numbers2colors(datTraits[,1], signed = FALSE);
# ##样本分类
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits),
#                     main = "Sample dendrogram and trait heatmap")

##选择一个soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#第power选取
##调用网络拓扑分析功能
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
##对sft结果进行可视化
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
##无标度拓扑拟合指数与软阈值功率的关系
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
##这条线对应于使用h的R ^ 2截止
abline(h=0.9,col="red")
##平均连通性与软阈值功率的关系
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sft$powerEstimate 
##自动选择power,调整参数，设置模块多少
# 报错：
# Error in (new("standardGeneric", .Data = function (x, y = NULL, use = "everything",  : 
#                                                   参数没有用(weights.x = NULL, weights.y = NULL, cosine = FALSE)
# 解决办法：函数冲突
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, deepSplit = 4,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)
table(net$colors)
cor<-stats::cor # 使用完后改回来
sizeGrWindow(12, 9)
##将标签转换为颜色以进行绘图
mergedColors = labels2colors(net$colors)
##在下方绘制树状图和模块颜色
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)

##量化module-trait 之间的关系
##定义基因和样品的数量
##用颜色标签重新计算ME
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
##MEs可以被认为是基因表达谱的代表
MEs = orderMEs(MEs0)
##Cibersort信息
load(file = 'Step5_CibersortResults.Rdata')
cib_df <- as.data.frame(PSE_result[rownames(MEs),-c((ncol(PSE_result)-2):ncol(PSE_result))])
cib_df$Lesion <- ifelse(rownames(cib_df) %in% rownames(datTraits), datTraits$source,'No')
cib_df$Lesion <- ifelse(cib_df$Lesion == 'lesion', 1, 0)
cib_df$No_Lesion <- ifelse(rownames(cib_df) %in% rownames(datTraits), datTraits$source,'No')
cib_df$No_Lesion <- ifelse(cib_df$No_Lesion == 'no-lesion', 1, 0)
  
datTraits1 <- cib_df[,c(24,25,1:23)]
datTraits1 <- datTraits1[,-c(12,15)] # 这两列没有表达量
dim(MEs)
dim(datTraits1)
moduleTraitCor = cor(MEs, datTraits1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
#＃将显示相关性及其p值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
##在热图图中显示相关值
pdf(file = './Fig/Step6/WGCNA_heatmap_group.pdf',width = 13,height = 6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colorRampPalette(c("#2CA02CFF","white","firebrick3"))(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,zlim = c(-1,1),
               xLabelsAngle = 45,
               main = paste("Relationships of Cell Infiltration (Cibersort) and Module Trait")) 
dev.off()


# 对p/brown(条件/模块)具体分析,模块内基因与表型数据关联:
# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达值算出相关系数
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要

# Firstly, calculate the correlation matrix between module and genes.(首先计算模块与基因的相关性矩)
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
# Calculate the Pearson correlation coefficient matrix of each module and its genes(算出每个模块跟基因的皮尔森相关系数矩)
# MEs is the value of each module in each sample.(MEs是每个模块在每个样本里面的)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Secondly, calculate the correlation matrix between conditions and genes. (再计算性状与基因的相关性矩)
# Only continuous properties can be computed. If the variables are discrete, the matrix is converted to 0-1 when the sample table is constructed.(只有连续型性状才能只有计算,如果是离散变量，在构建样品表时就转为0-1矩阵)
# Here, the variable whether or not it belongs to the P condition is numeralized with 0 and 1.(这里把是否属于 P 实验条件这个变量进行数值化，０和１表示)
group = 'Basal cell'
P = as.data.frame(datTraits1[,group]) # choose an interested condition!!
geneTraitSignificance = as.data.frame(cor(datExpr, P, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(P), sep="")
names(GSPvalue) = paste("p.GS.", names(P), sep="")

## Then, combine aboved two correlation matrixes(最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析)
colnames(geneModuleMembership)
module = 'yellow' # choose interested module
column = match(module, modNames)
moduleGenes = moduleColors == module ## get the genes in the interested module(获取模块内的基因)

library(ggsci)
library(scales)
pal = pal_d3()(10)
show_col(pal)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in Yellow Module"),
                   ylab = "Gene Significance for Basal cell Infiltration",
                   main = paste("Module membership vs Gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = "#FF7F0EFF"
)

## 优化散点图
# 筛选Hub基因
c <- data.frame(MM = abs(geneModuleMembership[moduleGenes, column]), 
                GS = abs(geneTraitSignificance[moduleGenes, 1]),
                row.names = rownames(geneModuleMembership[moduleGenes, ]))
head(c)
hub <- abs(c$MM)>0.8&abs(c$GS)>0.2
table(hub)

# 对基因进行分组
c$group<-ifelse(hub == TRUE, TRUE, FALSE)
head(c)
table(c$group)
hub1 <- subset(c, c$group == TRUE)

# 利用ggplot2作图
library(ggplot2)
ggplot(data=c,
       aes(x=MM, y=GS, color=group))+
  geom_point(size=1.5)+
  scale_colour_manual(values=c("grey60", "#FF7F0EFF"))+ 
  theme_bw()+  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+  
  labs(x="Module Membership in Yellow Module", 
       y="Gene Significance for Basal cell Infiltration",
       title = "Module membership vs Gene significance",
       subtitle = "Cor = 0.51, P < 0.001")+
  theme(axis.title.x =element_text(size=14), 
        axis.title.y=element_text(size=14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 12,face = "bold"),
        plot.margin = unit(rep(2,4),'lines')) +
  theme(legend.position = 'none')+
  geom_hline(aes(yintercept=0.2),colour="darkred",lwd=1,linetype=5)+
  geom_vline(aes(xintercept=0.8),colour="darkred",lwd=1,linetype=5)+
  scale_x_continuous(breaks=seq(0, 1, 0.2))+ ## X 轴每隔 5 个单位显示一个刻度
  scale_y_continuous(breaks=seq(0, 1, 0.2))
ggsave(filename = "./Fig/Step6/MMGS_ggplot2.pdf.pdf",width = 8,height = 6)


save(c, hub1, geneModuleMembership, geneTraitSignificance, datExpr, MEs,
     file = 'Step6_WGCNAOutput.Rdata')
write.csv(hub1,file = 'Step6_Hubgene.csv')


# Hub基因通路富集
library(clusterProfiler)
library(org.Hs.eg.db)
#keytypes(org.Mm.eg.db)
deg_entr <- bitr(rownames(hub1), fromType = "SYMBOL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db)
# KEGG
kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                     organism = 'hsa',
                     pvalueCutoff = 0.05)
# GO
go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                   OrgDb          = org.Hs.eg.db,
                   ont            = 'ALL', 
                   pAdjustMethod  = "BH",
                   pvalueCutoff   = 0.05, 
                   # qvalueCutoff   = 0.2, 
                   readable       = TRUE)

library(ggplot2)
kk_df1 <- as.data.frame(kk_deg@result)
kk_df2 <- kk_df1[order(kk_df1$Count,decreasing = T),][1:15,]
ggplot(data = kk_df2, mapping = aes(x= reorder(Description,Count),y=Count)) + 
  geom_point(aes(size=Count,color=log10(pvalue))) + #size 形状大小 color 形状颜色  shape 形状类型
  coord_flip() + #将X轴与Y轴倒过来
  labs(title= 'KEGG Enrichment: Hub Gene',x="Term",y = "Count") + # 添加标题
  scale_color_gradient(low="green",high ="red") +  # high，low是颜色变化的上下限
  theme(legend.title = element_text(size = 6, face = 2))+
  theme(legend.key.size=unit(0.6,'cm'))+
  theme(legend.text = element_text(size = 6,face = 'bold'))+
  theme(legend.position = "right")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5)) 
ggsave(filename = './Fig/Step6/Hub_Gene_KEGG_Top15.pdf', width = 10,height = 7)

library(stringr)
p3 = dotplot(go_deg, split="ONTOLOGY",showCategory = 15, title = 'GO Enrichment: Hub Gene')+ facet_grid(ONTOLOGY~., scale="free") + 
  scale_color_continuous(low='green', high='red')+ scale_y_discrete(labels = function(x) str_wrap(x, width = 70))
p3
ggsave(filename = './Fig/Step6/Hub_Gene_GO_Top15.pdf', p3, width = 10,height = 12)















