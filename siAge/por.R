library(data.table)
data <- fread("D:/HuaweiMoveData/Users/Zheng/Desktop/PIK3IP1_immue/SIAGE/pcor.txt")
library(ppcor)
library(fastDummies)

data$BMI <- as.numeric(data$BMI)
pcor_result <- pcor.test(
  x = data$CD8_TCM_HAVCR2, 
  y = data$siage, 
  z = data$BMI, 
  method = "spearman"
)
result <- as.data.frame(print(pcor_result))
aresult <- result[,c(1,2,6)]

