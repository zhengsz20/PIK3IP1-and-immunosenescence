library(data.table)
a <- read.delim("D:/HuaweiMoveData/Users/Zheng/Desktop/PIK3IP1_immue/distribution_test.txt")


ks_result <- ks.test(
  x = a$Group1,     
  y = a$Group2    
)

print(ks_result)
