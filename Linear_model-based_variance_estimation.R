library(lme4)
library(lmerTest) 
library(data.table)
Data <- fread("~/LMM_all_IHM.txt")

Data$Time <- factor(Data$Time)
Data$Sex <- factor(Data$Sex)

# fitting model
model <- lmer(
  IHM ~ Pik3ip1 + Time + Age  +Sex + (1 | ID),
  data = Data
)

# Check the results.
summary(model)    # Fixed effect coefficients and p-values
ranef(model)      # Individual random intercept deviation

Data <- fread("~/LMM_female_IHM.txt")
Data$Time <- factor(Data$Time)
Data$Sex <- factor(Data$Sex)

# fitting model
model <- lmer(
  IHM ~ Pik3ip1 + Time + Age + (1 | ID),
  data = Data
)
# Check the results.
summary(model)    
ranef(model)      
