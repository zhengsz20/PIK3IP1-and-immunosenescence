library(data.table)
data <- fread("~/SIAGE/explained_variance.txt")
data$Sex <- as.factor(data$Sex)
data$blood_type <- as.factor(data$blood_type)
data$BMI <- as.numeric(data$BMI)
# 构建模型
model <- lm(siage ~ Sex+BMI, data = data)
model <- lm(siage ~ blood_type, data = data)
model <- lm(siage ~ B_Atypical_Memory , data = data)
model <- lm(siage ~ CD8_TCM_HAVCR2, data = data)
model <- lm(siage ~ CD8_TEM_ZNF683, data = data)

# 查看模型摘要，获取调整后的R²
summary(model)
