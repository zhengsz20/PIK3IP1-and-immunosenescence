data <- read.delim("~/data.txt")
plot(data$A, data$B)
cor(data, method = "spearman")
cor.test(data$A, data$B, method = "spearman")

