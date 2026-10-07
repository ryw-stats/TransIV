source("../../R/model_functions.R", local = TRUE)

library(glmnet)
library(MASS)
library(data.table)

DATA <- as.data.frame(fread("BMI_SBP_EUR_SAS_SEL.csv"))
pval <- as.data.frame(fread("pval_BMI_SAS.csv"))
sel <- as.numeric(which(c(rep(T, 10), pval <= 1e-2)))
DATA <- DATA[, sel]

DATA$sex_male <- as.numeric(DATA$sex_male == "Male")

D <- DATA$bmi
D[is.na(D)] <- mean(D, na.rm = T)
XZ <- as.matrix(DATA[, 4:ncol(DATA)])
XZ[is.na(XZ)] <- 0

soucelasso_model_info <- train_lasso_model(XZ, D)
model_main <- soucelasso_model_info$model_residual
lambda_main <- soucelasso_model_info$lambda_residual

save(model_main, file = "source_lasso_EUR_SAS.Rdata")

DATA <- as.data.frame(fread("BMI_SBP_AFR_SAS_SEL.csv"))
sel <- as.numeric(which(c(rep(T, 10), pval <= 1e-2)))
DATA <- DATA[, sel]

DATA$sex_male <- as.numeric(DATA$sex_male == "Male")

D <- DATA$bmi
D[is.na(D)] <- mean(D, na.rm = T)
XZ <- as.matrix(DATA[, 4:ncol(DATA)])
XZ[is.na(XZ)] <- 0

soucelasso_model_info <- train_lasso_model(XZ, D)
model_main <- soucelasso_model_info$model_residual
lambda_main <- soucelasso_model_info$lambda_residual

save(model_main, file = "source_lasso_AFR_SAS.Rdata")
