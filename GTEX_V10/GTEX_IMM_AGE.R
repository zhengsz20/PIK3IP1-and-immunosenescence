library(GSVA)
library(pheatmap)
library(data.table)
library(GSEABase)
expr_data <- fread("E:/GTEX_v10/Z_blood.txt")
expr_data <- as.matrix(expr_data)
rownames(expr_data) <- expr_data[,1]
expr_data <- expr_data[,-1]
expr_data <- t(expr_data)
storage.mode(expr_data) <- "numeric"

##The matrix format with rows representing genes and columns representing samples.
gene_sets<-getGmt("D:/HuaweiMoveData/Users/Zheng/Desktop/PIK3IP1_immue/GTEX_v10/IMM_AGE_gene_set.txt")

gsva_data <-gsva(expr = expr_data,gset.idx.list= gene_sets, method ="ssgsea",kcdf="Gaussian")
gsva_data <- t(gsva_data)
gsva_data <- as.data.frame(gsva_data)
gsva_data$ID <- rownames(gsva_data)
fwrite(gsva_data,"gsva_data.txt",quote = F,sep = "\t",row.names = F)
