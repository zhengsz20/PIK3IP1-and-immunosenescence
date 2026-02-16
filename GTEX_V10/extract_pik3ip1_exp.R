library(data.table)

data <- fread("F:/GTEX_v10/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz")

data <- subset(data,data$Description=="PIK3IP1")
data <- t(data)
data <- as.data.frame(data[-c(1:2),])
data$ID <- rownames(data)
rownames(data) <- NULL
data$Gene <- "PIK3IP1"
data$data <- as.numeric(data$data)
data$LOG10 <- log10(data$data+1)
data <- data[,c(2,3,1,4)]
colnames(data)[3] <- "TPM"

data2 <- fread("F:/GTEX_v10/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
data <- merge(data,data2,by.x="ID",by.y = "SAMPID")
data <- data[,c(1:4,10)]

library(stringr)
data$Sample <- str_extract(data$ID, "^[^-]+-[^-]+")

data3 <- fread("F:/GTEX_v10/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
data <- merge(data,data3,by.x="Sample",by.y = "SUBJID")
data$Z_score<- scale(data$TPM,center = T,scale = T)

fwrite(data,"pik3ip1_exp.txt",quote = F,sep = "\t",row.names = F)
