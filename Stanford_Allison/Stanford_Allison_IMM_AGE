library(GSVA)
library(pheatmap)
library(data.table)
library(GSEABase)
expr_data <- fread("F:/Framingham Heart Study_2013_15/Z_GSE123698.txt")
expr_data <- na.omit(expr_data )
expr_data <- as.matrix(expr_data)
rownames(expr_data) <- expr_data[,1]
expr_data <- expr_data[,-1]
storage.mode(expr_data) <- "numeric"

gene_sets<-getGmt("D:/HuaweiMoveData/Users/Zheng/Desktop/PIK3IP1_immue/FHS_IMM_AGE/IMM_AGE_gene_set.txt")

gsva_data <-gsva(expr = expr_data,gset.idx.list= gene_sets, method ="ssgsea",kcdf="Gaussian")
gsva_data <- t(gsva_data)
gsva_data <- as.data.frame(gsva_data)
gsva_data$ID <- rownames(gsva_data)
rownames(gsva_data) <- NULL
fwrite(gsva_data,"IMM_AGE_ID_2015.txt",quote = F,sep = "\t",row.names = F)
