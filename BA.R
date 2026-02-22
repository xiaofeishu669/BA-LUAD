setwd("C:/R/R4.4.1")
# 加载包----
library(DESeq2)
library(edgeR)
library(limma)
library(TCGAbiolinks)
library(rjson)
library(tidyverse)
# 1 读取mRNA数据----
#获取样本名称及文件名称
json <- jsonlite::fromJSON("mRNA-metadata.cart.2025-04-03.json")
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)
count_file <-list.files('gdc-mRNA/',
                        pattern = '*.tsv',recursive=TRUE)
count_file_name<-strsplit(count_file,split = '/')
count_file_name<- sapply(count_file_name,function(x){x[2]})
#构建空数据框
matrix = data.frame(matrix(nrow=60660,ncol = 0))
#逐个读取并合并
for (i in 1:length(count_file)){
  path = paste0('gdc-mRNA//',count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data<-data[-c(1:6),]
  #3:unstranded,counts; 4:stranded_first;5:stranded_second;6:tpm_unstranded;7:fpkm_unstranded;8:fpkm_uq_unstranded
  data<-data[3]
  colnames(data) <-file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <-cbind(matrix,data)
}
#转化为gene_symbol
path = paste0('gdc-mRNA//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
matrix0 <-cbind(gene_name,matrix)
gene_type <-data[-c(1:6),2]
matrix0 <-cbind(gene_type,matrix0)
matrix0 <-aggregate(.~gene_name,data=matrix0,max)
table(gene_name)
#保留mRNA
matrix0<- subset(x =matrix0,gene_type=="protein_coding")
table(gene_name)
#将gene_name设为行名并转化为导出格式
rownames(matrix0) <- matrix0[,1]
matrix0 <-matrix0[,-c(1,2)]
matrix1 = data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1) = gsub('[.]','-',colnames(matrix1))
write.table(matrix1,'TCGA_LUAD_count.txt',sep="\t",quote = F,row.names = F)
# 2 edgeR----
library(edgeR)
# 1. 数据读取与预处理
data <- read.table("TCGA_LUAD_count.txt", header=T, sep="\t", 
                   row.names=1, check.names=F)
# 2. 精确提取样本类型（关键修改点）
sample_ids <- colnames(data)
sample_parts <- strsplit(sample_ids, "-")
sample_codes <- sapply(sample_parts, function(x) x[4])  # 提取第四段
sample_types <- substr(sample_codes, 1, 2)  # 取前两位字符（01/11）
# 3. 严格筛选01A肿瘤和11A正常样本
keep_samples <- sample_types %in% c("01", "11")
data <- data[, keep_samples]
sample_types <- sample_types[keep_samples]
# 4. 创建正确分组因子（设置Normal为参考组）
group <- factor(ifelse(sample_types == "01", "Tumor", "Normal"),
                levels = c("Normal", "Tumor"))  # 注意level顺序
# 查看样本分布
cat("样本数量统计:\n")
table(group)   #Normal  Tumor 59    539
# 创建样本信息表（在原代码基础上修改）
sample_info <- data.frame(
  Sample_ID = colnames(data),
  # 生成Normal和Tumor的二进制列
  Normal = ifelse(group == "Normal", 1, 0),
  Tumor = ifelse(group == "Tumor", 1, 0)
)
write.csv(sample_info, "sample_info.csv", row.names = FALSE)
# 5. 创建DGEList对象
DGElist <- DGEList(counts = data, group = group)
# 6. 智能过滤低表达基因（优于原代码的rowMeans方法）
keep_genes <- filterByExpr(DGElist, group = group)
DGElist <- DGElist[keep_genes, , keep.lib.sizes = FALSE]
# 7. 标准化与离散度估计
DGElist <- calcNormFactors(DGElist, method = "TMM")
design <- model.matrix(~ group)  # 自动生成对比矩阵
DGElist <- estimateDisp(DGElist, design)
# 8. 拟合模型与差异分析
fit <- glmQLFit(DGElist, design)
qlf <- glmQLFTest(fit, coef = 2)  # coef=2对应Tumor vs Normal
# 9. 结果提取与过滤
nrDEG_edgeR <- topTags(qlf, n = Inf)$table
sig_genes <- nrDEG_edgeR[
  nrDEG_edgeR$FDR < 0.05 & abs(nrDEG_edgeR$logFC) > 2,
]
# 查看Top差异基因
head(sig_genes[order(sig_genes$FDR), ])
# 检查已知标志物的表达方向
#sig_genes[grep("CD274", rownames(sig_genes)), ]
# 为所有结果添加基因名列
nrDEG_edgeR_withID <- nrDEG_edgeR %>%
  tibble::rownames_to_column("Gene_Symbol")
# 火山图----
#为显著结果添加基因名列
sig_genes_withID <- sig_genes %>%
  tibble::rownames_to_column("Gene_Symbol")
library(writexl)
write_xlsx(nrDEG_edgeR_withID, path = "DEG_results_all.xlsx")
write_xlsx(sig_genes_withID, path = "DEG_results_significant_2FOLD.xlsx")

# 添加差异表达标签
nrDEG_edgeR_withID$change <- ifelse(
  nrDEG_edgeR_withID$FDR < 0.05 & abs(nrDEG_edgeR_withID$logFC) > 2,
  ifelse(nrDEG_edgeR_withID$logFC > 2, "Up", "Down"),
  "Stable"
)
# 查看上调/下调/不变基因数量
table(nrDEG_edgeR_withID$change)

# 需要的对象：nrDEG_edgeR_withID（含 logFC, FDR）
stopifnot(exists("nrDEG_edgeR_withID"))
stopifnot(all(c("logFC","FDR") %in% colnames(nrDEG_edgeR_withID)))

## 1) 设定差异分类（与你原文一致：FDR<0.05 且 |logFC|>2）
nrDEG_edgeR_withID$change <- ifelse(
  nrDEG_edgeR_withID$FDR < 0.05 & abs(nrDEG_edgeR_withID$logFC) > 2,
  ifelse(nrDEG_edgeR_withID$logFC > 2, "Up", "Down"),
  "Stable"
)

# 查看数量
print(table(nrDEG_edgeR_withID$change))

## 2) 计算绘图所需的 max_neglogFDR / max_logFC（避免找不到对象）
nrDEG_edgeR_withID$neglogFDR <- -log10(nrDEG_edgeR_withID$FDR)

max_neglogFDR_raw <- max(nrDEG_edgeR_withID$neglogFDR, na.rm = TRUE)
max_logFC_raw <- max(abs(nrDEG_edgeR_withID$logFC), na.rm = TRUE)

# 先取上整，保证范围略大于数据
max_neglogFDR <- ceiling(max_neglogFDR_raw)
max_logFC <- ceiling(max_logFC_raw)

## 3) 坐标轴刻度（按你原代码：y每40，x每4）
y_breaks <- seq(0, max_neglogFDR, by = 40)
x_breaks <- seq(-max_logFC, max_logFC, by = 4)

# 把范围强制对齐到刻度端点
if(length(y_breaks) == 0) y_breaks <- c(0, max_neglogFDR)
if(length(x_breaks) == 0) x_breaks <- c(-max_logFC, 0, max_logFC)

max_logFC <- 12
max_neglogFDR <- 160
x_breaks <- seq(-12, 12, by = 4)
y_breaks <- seq(0, 160, by = 40)

## 4) 颜色：深边框 + 浅填充（更接近你原图视觉）
nrDEG_edgeR_withID$color_border <- ifelse(
  nrDEG_edgeR_withID$change == "Up",   "#C85A40",
  ifelse(nrDEG_edgeR_withID$change == "Down", "#315E94", "grey60")
)
nrDEG_edgeR_withID$color_inner <- ifelse(
  nrDEG_edgeR_withID$change == "Up",   "#FA8260",
  ifelse(nrDEG_edgeR_withID$change == "Down", "#4D8FD1", "grey75")
)

## 5) 百分比（与你原代码一致：Up/Down 占总基因的比例）
total_genes <- nrow(nrDEG_edgeR_withID)
up_percent <- sum(nrDEG_edgeR_withID$change == "Up") / total_genes * 100
down_percent <- sum(nrDEG_edgeR_withID$change == "Down") / total_genes * 100


png("volcano_plot_edgeR.png", width = 1200, height = 1000, res = 160)

par(bty = "n", mgp = c(1.5, 0.33, 0),
    mar = c(4.5, 3.5, 3, 3) + 0.1,
    las = 1, tcl = -0.3)

plot(nrDEG_edgeR_withID$logFC,
     nrDEG_edgeR_withID$neglogFDR,
     col = nrDEG_edgeR_withID$color_border,
     bg  = nrDEG_edgeR_withID$color_inner,
     pch = 21, cex = 1.5, lwd = 2.2,
     xlim = c(-max_logFC, max_logFC),
     ylim = c(0, max_neglogFDR),
     xlab = "", ylab = "-log10(FDR)",
     yaxt = "n", xaxt = "n")

axis(2, at = y_breaks, labels = y_breaks, lwd = 2.5)
axis(1, at = x_breaks, labels = x_breaks, lwd = 2.5)

abline(h = -log10(0.05), col = "grey60", lwd = 2.5, lty = 3)
abline(v = c(-2, 2), col = "grey60", lwd = 2.5, lty = 3)
abline(v = 0, col = "black", lwd = 2.5)

mtext(expression(paste(Log[2], " Fold Change")),
      side = 1, line = 2.2, cex = 1.2)

# 箭头位置（与你原逻辑一样）
arrow_y <- max_neglogFDR * 0.95
arrows(0.5,  arrow_y,  max_logFC*0.9,  arrow_y,
       length = 0.1, col = "#C85A40", lwd = 4)
arrows(-0.5, arrow_y, -max_logFC*0.9, arrow_y,
       length = 0.1, col = "#315E94", lwd = 4)

text(max_logFC*0.7,  arrow_y, "UP",   pos = 3, col = "#FA8260", cex = 1.1)
text(-max_logFC*0.7, arrow_y, "DOWN", pos = 3, col = "#4D8FD1", cex = 1.1)

text(max_logFC*0.7,  arrow_y*0.9, paste0(sprintf("%.1f", up_percent), "%"),
     col = "#FA8260", cex = 1.1)
text(-max_logFC*0.7, arrow_y*0.9, paste0(sprintf("%.1f", down_percent), "%"),
     col = "#4D8FD1", cex = 1.1)

title(main = "Volcano Plot")

dev.off()

pdf("volcano_plot_edgeR.pdf", width = 8, height = 6.7)

par(bty = "n", mgp = c(1.5, 0.33, 0),
    mar = c(4.5, 3.5, 3, 3) + 0.1,
    las = 1, tcl = -0.3)

plot(nrDEG_edgeR_withID$logFC,
     nrDEG_edgeR_withID$neglogFDR,
     col = nrDEG_edgeR_withID$color_border,
     bg  = nrDEG_edgeR_withID$color_inner,
     pch = 21, cex = 1.2, lwd = 1.8,
     xlim = c(-max_logFC, max_logFC),
     ylim = c(0, max_neglogFDR),
     xlab = "", ylab = "-log10(FDR)",
     yaxt = "n", xaxt = "n")

axis(2, at = y_breaks, labels = y_breaks, lwd = 2)
axis(1, at = x_breaks, labels = x_breaks, lwd = 2)

abline(h = -log10(0.05), col = "grey60", lwd = 2.5, lty = 3)
abline(v = c(-2, 2), col = "grey60", lwd = 2.5, lty = 3)
abline(v = 0, col = "black", lwd = 2.5)

mtext(expression(paste(Log[2], " Fold Change")),
      side = 1, line = 2.2, cex = 1.2)
arrow_y <- max_neglogFDR - 2

arrows(0.5,  arrow_y,  max_logFC*0.9,  arrow_y,
       length = 0.1, col = "#C85A40", lwd = 4)
arrows(-0.5, arrow_y, -max_logFC*0.9, arrow_y,
       length = 0.1, col = "#315E94", lwd = 4)

text(max_logFC*0.7,  arrow_y, "UP",   pos = 3, col = "#FA8260", cex = 1.1)
text(-max_logFC*0.7, arrow_y, "DOWN", pos = 3, col = "#4D8FD1", cex = 1.1)

text(max_logFC*0.7,  arrow_y*0.9, paste0(sprintf("%.1f", up_percent), "%"),
     col = "#FA8260", cex = 1.1)
text(-max_logFC*0.7, arrow_y*0.9, paste0(sprintf("%.1f", down_percent), "%"),
     col = "#4D8FD1", cex = 1.1)

title(main = "Volcano Plot")

dev.off()

#WGCNA-----
library(limma)
library(WGCNA)

expFile <- "TPM.txt"
notesFile <- "notes.txt"

expression_data <- read.table(expFile,sep = "\t",header = TRUE,check.names = FALSE)
expression_data_matrix<- as.matrix(expression_data) #转换为矩阵格式
rownames(expression_data_matrix) <- expression_data_matrix[,1] #第一列作为行名（基因ID）
expression_data<- expression_data_matrix[,2:ncol(expression_data_matrix)] #去掉第一列（基因ID列）
# 转换为数值矩阵并保留行列名
data_matrix <- matrix(as.numeric(as.matrix(expression_data)),
                      nrow=nrow(expression_data),
                      dimnames=dimnames(expression_data))
data_matrix <- avereps(data_matrix) #对重复基因取均值
# 降低数据偏态  如果表达量数据过大,可以取消下面这行的注释以进行1og2转换
data_matrix <- log2(data_matrix + 1) 
#对表达矩阵进行数组间的标准化（消批次效应）
#data_matrix <- normalizeBetweenArrays(data_matrix)
# 过滤低变异基因（标准差<0.5的基因被移除）
data_matrix <- data_matrix[apply(data_matrix, 1, sd) > 0.5, ]
# 读取样本信息
sample_info <- read.table(notesFile, sep = ",", header = TRUE) 
colnames(sample_info) <- c("Sample_ID", "group")  # 标准化列名
# 统一分组标签（Normal/Tumor）
sample_info$group <- ifelse(grepl("normal", sample_info$group, ignore.case = TRUE), "Normal", "Tumor")
# 提取正常和肿瘤样本ID
normal_samples <- sample_info[sample_info$group == "Normal", "Sample_ID"]
Tumor_samples <- sample_info[sample_info$group == "Tumor", "Sample_ID"]
# 按分组拆分表达矩阵
control_data <- data_matrix[, as.vector(normal_samples)]
treat_data <- data_matrix[, as.vector(Tumor_samples)]
# 合并为新的矩阵（正常样本在前，肿瘤样本在后）
data_matrix <- cbind(control_data, treat_data)
# 记录各组样本数
control_count <- ncol(control_data)
treat_count <- ncol(treat_data)
#检查处理缺失值
dat_expr <- t(data_matrix)   # 转置矩阵（WGCNA要求行=样本，列=基因）
gsg <- goodSamplesGenes(dat_expr, verbose = 3)      # 检查基因和样本质量
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {     # 移除低质量基因
    cat("Removing genes:", paste(names(dat_expr)[!gsg$goodGenes], collapse = ","))
  }
  if (sum(!gsg$goodsamples) > 0) {   # 移除低质量样本
    cat("Removing samples:", paste(rownames(dat_expr)[!gsg$goodSamples], collapse = ","))
  }
  dat_expr <- dat_expr[gsg$goodsamples, gsg$goodGenes]   # 保留合格数据
}
# 1 样品聚类可视化----
sample_tree <- hclust(dist(dat_expr),method = "average")  # 层次聚类（使用平均距离法）
pdf("1.样品聚类.pdf",width = 15,height = 10)
par(cex =0.6,mar =c(0,4,2,0))   # 全局文字缩放比例`/ `下左上右的边距（单位：行）
plot(sample_tree,main = "Sample clustering to detect outliers",xlab="",cex.lab = 1.5,cex.axis =1.5,cex.main = 1.2)  #不显示X坐标，坐标轴标签大小，坐标轴刻度文字大小，主标题大小
abline(h = 20000,col = "red")  # 添加红色参考线（高度=20000）
dev.off()
# 删除聚类中高度较低的样品
#clust <- cutreeStatic(sample_tree, cutHeight = 20000, minSize = 10)
#table(clust)
#keep_samples <- clust == 1
#dat_expr <- dat_expr[keep_samples, ]

# 准备临床数据并匹配样本
clinical_data <- data.frame(Normal = rep(1, control_count),   # 正常样本标记为1，肿瘤为0
                            Tumor = rep(0, control_count),
                            stringsAsFactors = FALSE)
clinical_data <- rbind(clinical_data,
                       data.frame(Normal = rep(0, treat_count),     # 合并肿瘤样本数据（正常标记0，肿瘤标记1）
                                  Tumor = rep(1, treat_count),
                                  stringsAsFactors = FALSE))
rownames(clinical_data) <- colnames(data_matrix)         # 设置行名为样本ID（需与表达矩阵列名一致）

# 样本匹配
fpkm_samples <- rownames(dat_expr)                        # 表达矩阵样本名
trait_samples <- rownames(clinical_data)                  # 临床数据样本名
same_samples <- intersect(fpkm_samples, trait_samples)    # 取交集
dat_expr <- dat_expr[same_samples, ]                    # 筛选共有样本
clinical_data <- clinical_data[same_samples, ]

# 2 样品聚类热图----
sample_tree2 <- hclust(dist(dat_expr), method = "average")        # 重新聚类（确保使用匹配后的数据）
trait_colors <- numbers2colors(clinical_data, signed = FALSE)     # 将临床数据转换为颜色编码

pdf("2.样品聚类热图.pdf", width = 15, height = 10)
plotDendroAndColors(sample_tree2, trait_colors,                   # 颜色矩阵
                    groupLabels = names(clinical_data),           # 分组标签
                    main = "Sample dendrogram and trait heatmap") # 主标题
dev.off()
# 3 power值散点图----
library(WGCNA)
enableWGCNAThreads()  #多线程工作
powers <- 1:20  #测试幂指数范围1:20
sft <- pickSoftThreshold(dat_expr, powerVector = powers, verbose = 5)  #计算软阈值指标
pdf("3. power值散点图.pdf", width = 15, height = 10)
par (mfrow = c(1, 2))   # 1行2列布局
cex1 <- 0.9
# 左图：无标度拓扑拟合
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")    # 推荐R²阈值线
# 右图：平均连通性
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
#左图寻找首个R²≥0.9的power值（红线上方），右图确保平均连接度>50
dev.off()
# 获取自动推荐的power值
soft_power <- sft$powerEstimate  #最优power值
# 计算邻接矩阵
adjacency <- adjacency(dat_expr, power = soft_power)
soft_power
# 手动验证power=9的连接度分布
k <- softConnectivity(dat_expr, power = 9)
scaleFreePlot(k, main = "Power=9 的连接度分布")

#TOM矩阵
TOM <- TOMsimilarity(adjacency)   # 基于邻接矩阵计算拓扑重叠矩阵(TOM)
dissTOM <-1 -TOM                  # 将TOM转换为相异度矩阵（范围0-1，0=高度相似，1=完全不相似）
# 4 基因聚类----
gene_tree <- hclust(as.dist(dissTOM),method = "average")   # 基于TOM相异度进行层次聚类（平均连接法）
pdf("4.基因聚类.pdf",width = 15,height = 10)
plot(gene_tree,xlab = "",sub = "",main = "Gene clustering on TOM-based dissmilarity",   # 不显示坐标轴标签
     labels = FALSE,hang = 0.04)   # 不显示基因名标签  # 所有分支底部对齐
dev.off()

min_module_size <- 50    # 设置最小模块基因数为50
# 动态剪裁树状图识别模块
mic_mods <- cutreeDynamic(dendro = gene_tree, distM = dissTOM,
                          deepSplit = 2, pamRespectsDendro = FALSE,    # 分割深度（0-4，越大模块越多）  # 是否优先保持树结构
                          minClusterSize = min_module_size)            # 模块最小基因数
dynamic_colors <- labels2colors(mic_mods)               # 将模块编号转换为颜色标签      
# 5 基因模块对应关系----
pdf("5.基因模块对应关系.pdf",width=15,height=10)
plotDendroAndColors(gene_tree, dynamic_colors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,    # 不显示基因标签  # 分支悬挂长度
                    addGuide = TRUE, guideHang = 0.05,    #添加辅助虚线     # 辅助线悬挂长度
                    main = "Gene dendrogram and module colors")
dev.off()

# 6 模块聚类----
# 计算模块特征基因（第一主成分）
MEList <- moduleEigengenes(dat_expr, colors = dynamic_colors)
MEs <- MEList$eigengenes    # 提取特征基因矩阵
MEDiss <- 1 - cor(MEs)      # 计算模块间相关性（1 - Pearson相关系数）
METree <- hclust(as.dist(MEDiss), method = "average")   # 基于模块相关性进行层次聚类

pdf("6.模块聚类.pdf", width = 7, height = 6)
plot(METree, main = "clustering of module eigengenes",
     xlab = "", sub ="")
MEDissThres <- 0.25   # 设置合并阈值（相关系数≥0.75的模块合并） 1 - 0.75 = 0.25
abline(h = MEDissThres, col = "red")     # 添加参考线
dev.off()
# 7 合并高度相关的模块（基于预设的MEDissThres阈值）----
merge <- mergeCloseModules(dat_expr, 
                           dynamic_colors,  # 动态切割得到的模块颜色标签
                           cutHeight = MEDissThres,   # 合并阈值（当模块特征基因间的相异度（1 - cor）MEDiss < 0.25时合并）
                           verbose = 3)            # 输出详细信息
# 提取合并后的模块颜色标签和特征基因
merged_colors <- merge$colors  # 新模块颜色
merged_MEs <- merge$newMEs     # 新模块特征基因（第一主成分）
pdf("7.合并后模块.pdf",width=15,height=10)
plotDendroAndColors(gene_tree, merged_colors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# 8 模块与性状热图----
# 计算模块特征基因与临床性状的相关性及p值
module_trait_cor <- cor(merged_MEs, clinical_data, use = "p") # Pearson相关
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(dat_expr)) # 计算p值
pdf("8.模块与性状热图.pdf",width=10, height=10)
# 生成带统计值的文本标签（格式：相关系数\n(p值)）
text_matrix <- paste(signif(module_trait_cor, 2), "\n(",                    # 保留2位有效数字的相关系数
                     signif(module_trait_pvalue, 1), ")", sep ="")          # 保留1位有效数字的p值
dim(text_matrix) <- dim(module_trait_cor)          # 保持矩阵维度一致
#绘制热图
par(mar = c(5, 10, 3, 3))   # 调整边距（下,左,上,右）
labeledHeatmap(Matrix = module_trait_cor,            # 输入相关系数矩阵
               xLabels = names(clinical_data),      # X轴标签（如Normal/Tumor）
               yLabels = names (merged_MEs),               # Y轴标签（模块名称）
               ySymbols = names (merged_MEs),              # Y轴符号（同模块名称）
               colorLabels = FALSE,                 # 不使用颜色标签
               colors = blueWhiteRed(50),           # 颜色梯度（蓝-白-红，50级）
               textMatrix = text_matrix,            # 添加文本标签
               setStdMargins = FALSE,               # 禁用自动边距
               cex.text = 0.5,                      # 文本大小
               zlim = c(-1, 1),                     # 颜色映射范围（-1到1）
               main = "Module-trait relationships")
dev.off()

# 查看各模块基因数量分布
table(merged_colors)
modules <- unique(merged_colors)
for (mod in modules) {
  genes_in_module <- colnames(dat_expr)[merged_colors == mod]
  write.table(genes_in_module,
              file = paste0("Merged_Module_", mod, "_genes.txt"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}










# 9 Hub基因筛选----
library(WGCNA)
#可视化GS和MM,GS所有基因表达谱与这个模块的eigengene相关性cor，每个值代表这个基因与模块间关系
#如果这个值绝对值接近0，这个基因就不是这个模块中一部分，绝对值接近1则与模块高度相关
#MM基因和表型性状如肿瘤正常的相关性的绝对值

# 选MEsalmon
module <- "salmon"                         # 目标模块颜色
Selectedclinical <- "Tumor"                # 目标性状名称

# 提取目标性状数据
datTraits <- clinical_data                 # 确保临床数据变量名与模板一致
Selectedclinical_df <- as.data.frame(datTraits[, Selectedclinical])
names(Selectedclinical_df) <- Selectedclinical

# 准备模块特征基因数据
MEs <- merged_MEs                          # 使用合并后的模块特征基因
modNames <- substring(names(MEs), 3)       # 去除特征基因名前缀"ME"

# 获取表达矩阵(样本x基因)
datExpr0 <- t(data_matrix)                 # 原始表达矩阵需要转置为样本x基因格式
datExpr0 <- datExpr0[rownames(dat_expr), ] # 确保样本顺序与过滤后一致
datExpr1 <- datExpr0                       # 创建副本用于分析

# 计算基因模块成员关系(MM)
nSamples <- nrow(dat_expr)                 # 样本数量
geneModuleMembership <- as.data.frame(cor(datExpr1, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# 命名列
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

# 计算基因性状相关性(GS)
geneTraitSignificance <- as.data.frame(cor(datExpr1, Selectedclinical_df, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

# 命名列
names(geneTraitSignificance) <- paste0("GS.", Selectedclinical)
names(GSPvalue) <- paste0("p.GS.", Selectedclinical)

# 筛选目标模块基因
module_column <- match(module, modNames)   # 定位目标模块
moduleGenes <- merged_colors == module     # 获取模块内基因

# 绘制相关性散点图
# 绘制相关性散点图（移除回归线，添加阈值线）
pdf(paste0("0.", Selectedclinical, "_", module, "_MMvsGS.pdf"), width = 7, height = 7)
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, module_column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = paste("Gene significance"),
  main = paste("Module membership vs. Gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
  col = module, 
  abline = FALSE  # 禁用默认回归线
)
# 添加自定义阈值线
abline(h = 0.5, v = 0.8, col = "black", lty = 2)
dev.off()

library(dplyr)
# 保存Hub基因列表(MM > 0.8 & GS > 0.5)
hub_genes <- data.frame(
  Gene = colnames(datExpr1)[moduleGenes],
  MM = geneModuleMembership[moduleGenes, module_column],
  GS = geneTraitSignificance[moduleGenes, 1]
) %>% 
  filter(abs(MM) > 0.8 & abs(GS) > 0.5) %>% 
  arrange(desc(abs(GS)))

write.csv(hub_genes, 
          file = paste0("HubGenes_", module, "_", Selectedclinical, ".csv"),
          row.names = FALSE)

# 合并所有基因的MM、GS及其p值
all_genes_results <- data.frame(
  Gene = colnames(datExpr1),                            # 所有基因名称
  Module = merged_colors,                               # 基因所属模块
  MM = geneModuleMembership[, paste0("MM", module)],   # 目标模块的MM值
  MM_pval = MMPvalue[, paste0("p.MM", module)],         # MM的p值
  GS = geneTraitSignificance[, paste0("GS.", Selectedclinical)],  # GS值
  GS_pval = GSPvalue[, paste0("p.GS.", Selectedclinical)]        # GS的p值
)

# 按MM绝对值降序排列（查看模块内共表达强度）
all_genes_sorted_by_MM <- all_genes_results %>% 
  arrange(desc(abs(MM)))

# 按GS绝对值降序排列（查看与性状关联强度）
all_genes_sorted_by_GS <- all_genes_results %>% 
  arrange(desc(abs(GS)))

# 输出为CSV文件
write.csv(all_genes_results, 
          file = paste0("AllGenes_MM_GS_", module, "_", Selectedclinical, ".csv"), 
          row.names = FALSE)



#GO----
#气泡一
library(clusterProfiler)
library(org.Hs.eg.db)
# 读取final.csv文件
final_genes <- read.csv("final.csv", stringsAsFactors = FALSE)
# 确保final_genes存在且包含Gene列
if (!exists("final_genes") || !"Gene" %in% colnames(final_genes)) {
  stop("final_genes数据框或Gene列不存在，请检查数据")
}

# 提取Gene列中的基因符号（需确保是字符向量）
gene_list <- final_genes$Gene

# GO富集分析（关键修改：gene参数替换为gene_list）
enrich.go <- enrichGO(
  gene = gene_list,            # 替换为final_genes$Gene
  OrgDb = org.Hs.eg.db,        # 人类基因数据库
  keyType = "SYMBOL",          # 输入基因类型为Symbol
  ont = "ALL",                 # 包括BP/MF/CC
  pAdjustMethod = "BH",        # 校正方法
  pvalueCutoff = 0.05,         # p值阈值
  qvalueCutoff = 0.05,          # q值阈值
  minGSSize = 10,              # 最小基因集大小
  readable = FALSE             # 不自动转换ID
)

# 保存结果
go <- as.data.frame(enrich.go)
write.table(go, "enrich.go.txt", sep = "\t", row.names = FALSE, quote = FALSE)
library("ggplot2")
pdf(file = 'GO_ALL_DOT.pdf', width = 10, height = 13) 
dotplot(enrich.go,split="ONTOLOGY",showCategory=10)+
  facet_grid(ONTOLOGY~.,scales = "free") +
  theme(axis.text.y = element_text(size = 12)) # 增加y轴文字大小
dev.off()
#KEGG-----
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
# 1. 基因名转换（增加NA值处理）
entrezid_all <- mapIds(
  org.Hs.eg.db,
  keys = unique(gene_list), 
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

entrezid_all <- unique(na.omit(as.character(entrezid_all)))  # ← 去 NA + 去重

cat("成功转换的基因数量:", length(entrezid_all), "/", length(gene_list), "\n")
# 3. KEGG富集分析
enrich.kegg <- enrichKEGG(
  gene = entrezid_all,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
# 4. 转换为数据框（修正错误）
if (!is.null(enrich.kegg)) {
  kegg_df <- as.data.frame(enrich.kegg)
  write.csv(kegg_df, "kegg_enrichment_results.csv", row.names = FALSE)
} else {  warning("KEGG富集无显著结果，请放宽阈值或检查输入基因")}

#将ENTREZID转换为基因符号
kegg_df$geneID <- sapply(kegg_df$geneID, function(entrez_ids) {
  # 分割字符串获取单个ENTREZID
  ids <- strsplit(entrez_ids, "/")[[1]]
  # 转换ENTREZID为基因符号
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ids,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  # 合并符号并处理NA值
  paste(na.omit(symbols), collapse = "/")
})

# 保存结果
write.csv(kegg_df, "kegg_enrichment_results.csv", row.names = FALSE)

#气泡图
pdf (file = "KEGG_dot2.pdf",width=8,height = 9)
enrich.kegg <- setReadable(enrich.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
dotplot(enrich.kegg, showCategory = 20, title = "KEGG Enrichment")
dev.off ()


#CluserKEGG-----
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readxl)

d1 <- read_excel("1.xlsx")
gene_list1 <- d1$Gene

# 1) SYMBOL -> ENTREZ（去NA + 去重，保证稳定）
entrezid_all <- mapIds(
  org.Hs.eg.db,
  keys = unique(gene_list1),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)
entrezid_all <- unique(na.omit(as.character(entrezid_all)))

cat("成功转换的基因数量:", length(entrezid_all), "/", length(unique(gene_list1)), "\n")

# 2) KEGG enrichment
enrich.kegg <- enrichKEGG(
  gene = entrezid_all,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

if (is.null(enrich.kegg) || nrow(as.data.frame(enrich.kegg)) == 0) {
  stop("KEGG富集无显著结果：请放宽阈值或检查输入基因。")
}

# 3) ✅ 把 enrich.kegg 对象里的 geneID 从 ENTREZ 直接转成 SYMBOL（图和表都会变）
enrich.kegg <- setReadable(enrich.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# 4) 导出（geneID 已经是 SYMBOL 了）
kegg_df <- as.data.frame(enrich.kegg)
write.csv(kegg_df, "kegg_results1.csv", row.names = FALSE)

# 5) 画图
pdf("KEGG_dotcluser1.pdf", width = 8, height = 9)
print(dotplot(enrich.kegg, showCategory = 20, title = "KEGG Enrichment"))
dev.off()


# 免疫浸润分析-----
library(devtools)
library(e1071)
library (preprocessCore)
library(parallel)
library(ggplot2)
library(CIBERSORT)
library(corrplot)
library(vioplot)

data=read.table("TPM.txt", header=T, sep="\t", check.names=F, row.names = 1)
dimnames=list(rownames (data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
results <- cibersort(sig_matrix = LM22, mixture_file = data)
results=as.matrix(results[,1:(ncol(results)-3)])
results=rbind(id=colnames(results), results)
write.table(results, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)

immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1, check.names=F)
immune=as.matrix(immune)
immune_data=t (immune)
#柱状图1-----
col=rainbow(nrow(immune_data),s=0.7,v=0.7)
pdf("CIBERSORT1.pdf",height=10,width=22)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1 = barplot(immune_data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n", cex.lab=1.8)
a2 = axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(immune_data),adj=1,cex=0.6);par(srt=0)
ytick2 =cumsum(immune_data[,ncol(immune_data)])
ytick1=c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(immune_data),col=col,pch=15,bty="n",cex=1.3)
dev.off()
#柱状图2-----
cellnum <- read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
cell.prop <- apply(cellnum, 1, function(x){x/sum(x)})  #将绝对值转为比例（每行样本细胞比例总和为1）
my36colors <- c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3',
                '#476D87','#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C',
                '#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4c755',
                '#F7F398','#AA9A59','#E63863','#E39A35','#C1E6F3','#6778AE','#91D0BE',
                '#B53E2B','#712820','#DCC1DD','#CCE0F5','#CCC9E6','#624D9E','#68A180',
                '#3A6963','#968175')
data4plot <- data.frame()  #创建新数据框转化数据为长格式（绘图）
for(i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],rownames(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}
colnames (data4plot) <- c('proportion', 'celltype' , 'sample') #设列名
#转换比例为数值类型
data4plot$proportion <- as.numeric(data4plot$proportion)
pdf(file="CIBERSORT2.pdf",height=10,width=22)
ggplot(data4plot,aes(sample,proportion,fill=celltype))+
  geom_bar (stat="identity",position="fill")+
  scale_fill_manual(values=my36colors)+ #自定义fill的颜色
  labs(x = "ID", title = "Cell Proportion") +  # 修改x轴标签为"ID"
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.5, 'cm'),
    axis.text.x = element_blank(),        # 隐藏x轴样本标签
    axis.ticks.x = element_blank(),       # 隐藏x轴刻度线
    axis.title.x = element_text(size=14)  # 设置x轴标题字体
  ) +
  guides(fill=guide_legend(title=NULL))
dev.off()

#相关性图-----
pdf("cor.pdf",height=13,width=13)
par(oma=c(0.5,1,1,1.2))
immune=immune [,colMeans(immune)>0]  #过滤低丰度细胞类型（保留在至少部分样本中存在的细胞类型）
M=cor (immune)  #计算细胞类型间的Pearson相关系数矩阵
corrplot(M,
         order="hclust",  # 使用层次聚类重组矩阵行列
         method = "color",  # 颜色方块表示相关强度
         addCoef.col ='black',  #在色块上显示黑色相关系数
         diag = TRUE,  #显示对角线（自身相关=1）
         tl.col="black",  #设置标签颜色为黑色
         col=colorRampPalette(c("blue","white","red"))(50))  #蓝-白-红渐变色
dev.off()


#制作小提琴图-----
rt=read.table("CIBERSORT-Results.txt", header=T, sep="\t",check.names=F, row.names=1)
# 提取样本类型字段（第4部分）
sample_type <- sapply(strsplit(rownames(rt), "-"), "[", 4)
# 提取前两位字符进行筛选
sample_type_part <- substr(sample_type, 1, 2)  # 关键修改点
keep <- sample_type_part %in% c("01", "11")    # 筛选01xx和11xx样本
# 创建分组变量（0=肿瘤，1=正常）
group <- ifelse(sample_type_part[keep] == "01", 0, 1)  # 0=Tumor, 1=Normal
# 保持原代码变量名不变
rt <- rt[keep, ]  # 应用筛选
conNum=length(group[group == "1"])  # 正常样本数
treatNum=length(group[group == "0"]) # 肿瘤样本数
# 保持原排序逻辑
rt1 = rt[group == "1", ]
rt2 = rt[group == "0", ]
rt = rbind(rt1, rt2)
#绘制小提琴图
outTab=data.frame()
pdf(file="vioplot.pdf", width=13, height=8)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63),ylim=c(min(rt), max(rt)+0.05),
     main="", xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
rt2_mat <- rt2   # 备份肿瘤样本矩阵
#免疫细行循环
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  #提取数据
  rt1=rt[1:conNum, i]
  rt2=rt[(conNum+1):(conNum+treatNum),i]
  #绘制小提琴图
  vioplot(rt1,at=3*(i-1),lty=1,add= T,col ='blue')
  vioplot(rt2,at=3*(i-1)+1,Ity=1,add = T,col = 'red')
  #统计检验
  wilcoxTest=wilcox.test(rt1,rt2)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }        
  mx=max(c(rt1,rt2))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  #加上pvalue
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",sprintf("%.03f",p))),cex=0.8)
}   
legend("topright",
       c("Normal", "Tumor"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
text(seq(1,64,3),-0.04,xpd=NA,labels=colnames(rt),cex=1,srt=45,pos=2)
dev.off()
#输出免疫细胞和p
write.table(outTab,file="Diff.xls",sep="\t",row.names=F,quote=F)
#热图-----
# 重新还原肿瘤/正常矩阵
rt_raw <- read.table("CIBERSORT-Results.txt", header=TRUE, sep="\t",
                     row.names=1, check.names=FALSE)

# 提取样本类型（第4段）
sample_type <- sapply(strsplit(rownames(rt_raw), "-"), "[", 4)
sample_type_part <- substr(sample_type, 1, 2)

# 只保留 01(肿瘤) 与 11(正常)
keep <- sample_type_part %in% c("01", "11")
rt_raw <- rt_raw[keep, , drop=FALSE]
sample_type_part <- sample_type_part[keep]

# 还原矩阵：rt2_mat=肿瘤(01), rt1_mat=正常(11)
rt2_mat <- rt_raw[sample_type_part == "01", , drop=FALSE]  # Tumor
rt1_mat <- rt_raw[sample_type_part == "11", , drop=FALSE]  # Normal

cat("Tumor matrix dim:", dim(rt2_mat), "\n")
cat("Normal matrix dim:", dim(rt1_mat), "\n")
cat("Tumor rownames NULL?", is.null(rownames(rt2_mat)), "\n")
# 1. 提取肿瘤样本ID
tumor_samples <- rownames(rt2_mat)

# 2. 从TPM里提取这些样本（关键：先做交集，否则仍可能空）
tumor_samples <- intersect(tumor_samples, colnames(data))
data_tumor <- data[, tumor_samples, drop=FALSE]

# 3. 转置
expr_data <- as.data.frame(t(data_tumor))

# 4. 免疫丰度同样按相同样本与顺序提取
tumor_samples <- intersect(tumor_samples, rownames(immune))
expr_data <- expr_data[tumor_samples, , drop=FALSE]
immu_data <- as.data.frame(immune[tumor_samples, , drop=FALSE])

cat("expr_data dim:", dim(expr_data), "\n")
cat("immu_data dim:", dim(immu_data), "\n")
# 加载必要的包
library(pheatmap)
library(reshape2)

# 初始化存储所有基因结果的列表
all_cor <- list()
all_pval <- list()

# 修改原有循环，收集所有结果
for (i in genes) {
  y <- as.numeric(expr_data[, i])
  cor_data <- do.call(rbind, lapply(colnames(immu_data), function(x) {
    dd <- cor.test(as.numeric(immu_data[, x]), y, method = "spearman", exact = FALSE)
    data.frame(gene = i, cell = x, 
               cor = dd$estimate, 
               p.value = dd$p.value)
  }))
  
  all_cor[[i]] <- cor_data$cor
  all_pval[[i]] <- cor_data$p.value
}

cor_matrix <- acast(do.call(rbind, lapply(genes, function(g) {
  data.frame(gene = g, cell = colnames(immu_data), cor = all_cor[[g]])
})), gene ~ cell, value.var = "cor")

pval_matrix <- acast(do.call(rbind, lapply(genes, function(g) {
  data.frame(gene = g, cell = colnames(immu_data), p.value = all_pval[[g]])
})), gene ~ cell, value.var = "p.value")

sig_matrix <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
sig_matrix[pval_matrix < 0.05] <- "*"
sig_matrix[pval_matrix < 0.01] <- "**"
sig_matrix[pval_matrix < 0.001] <- "***"

print(range(cor_matrix, na.rm = TRUE))  # 输出数据实际的最小值和最大值
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-0.4, 0.4, length.out = 101)

pdf("Gene_Immune_Correlation_Heatmap.pdf", family = "Times", 
    width = 12, height = 10)  # 增大画布尺寸

pheatmap(cor_matrix,
         color = color_palette,
         breaks = breaks,
         display_numbers = sig_matrix,
         number_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 10,  # 减小行字体
         fontsize_col = 8,   # 减小列字体
         angle_col = 45,     # 垂直显示列标签
         margins = c(5, 2), # 调整边距（底部，右侧）
         main = "Gene-Immune Cell Correlation",
         cellwidth = 20,     # 固定单元格宽度
         cellheight = 50     # 固定单元格高度
)

dev.off()

