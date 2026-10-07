iv_functions_file <- normalizePath("../../R/iv_functions.R", mustWork = TRUE)
source(iv_functions_file, local = TRUE)

source("../../R/model_functions.R", local = TRUE)

library(glmnet)
library(MASS)
library(data.table)

DATA <- as.data.frame(fread("BMI_SBP_AFR_SEL.csv"))
pval <- as.data.frame(fread("pval_BMI_AFR.csv"))
beta_GWAS <- as.data.frame(fread("beta_BMI_AFR.csv"))

sel <- as.numeric(which(c(rep(T, 10), pval <= 1e-2)))
DATA <- DATA[, sel]

DATA$sex_male <- as.numeric(DATA$sex_male == "Male")

D <- DATA$bmi
D[is.na(D)] <- mean(D, na.rm = T)
XZ <- as.matrix(DATA[, 4:ncol(DATA)])
XZ[is.na(XZ)] <- 0

load("source_lasso_EUR_AFR.Rdata")
D_source1 <- predict(model_main, XZ)

load("source_lasso_SAS_AFR.Rdata")
D_source2 <- predict(model_main, XZ)

residual1 <- D - D_source1

residual2 <- D - D_source2

n <- length(D)
K_split <- 3

set.seed(-1)
split_ind <- cut(sample(1:n, replace = F), K_split, labels = F)

D_Target <- rep(NA, n)
D_FT1 <- rep(NA, n)
D_FT2 <- rep(NA, n)

for (k in 1:K_split) {
  targetlasso_model_info <- train_lasso_model(XZ[split_ind != k, ], D[split_ind != k])
  model_target <- targetlasso_model_info$model_residual
  D_Target[split_ind == k] <- predict(model_target, XZ[split_ind == k, ])

  FTlasso_model_info <- train_lasso_model(XZ[split_ind != k, ], residual1[split_ind != k])
  model_FT <- FTlasso_model_info$model_residual
  D_FT1[split_ind == k] <- predict(model_FT, XZ[split_ind == k, ]) + D_source1[split_ind == k]

  FTlasso_model_info <- train_lasso_model(XZ[split_ind != k, ], residual2[split_ind != k])
  model_FT <- FTlasso_model_info$model_residual
  D_FT2[split_ind == k] <- predict(model_FT, XZ[split_ind == k, ]) + D_source2[split_ind == k]
}

TargetLasso_tmp <- c()
FTLasso1_tmp <- c()
FTLasso2_tmp <- c()
CF_tmp <- c()
KPI_tmp <- c()
v_TargetLasso_tmp <- c()
v_FTLasso1_tmp <- c()
v_FTLasso2_tmp <- c()
v_KPI_tmp <- c()
Strength_TargetLasso_tmp <- c()
Strength_FTLasso1_tmp <- c()
Strength_FTLasso2_tmp <- c()
Strength_KPI_tmp <- c()

IV_all <- cbind(D_Target, D_source1, D_FT1, D_source2, D_FT2)

Y_all <- DATA$sbp
Y_all[is.na(Y_all)] <- mean(Y_all, na.rm = T)
D_all <- D

X <- cbind(1, XZ[, 1:7])
MX <- diag(1, n) - X %*% solve(t(X) %*% X) %*% t(X)
IV <- MX %*% IV_all[, 2]
Y <- c(MX %*% Y_all)
D <- c(MX %*% D_all)
Source1 <- sum(IV * Y) / sum(IV * D)

Strength_Source1 <- n * cov(D, IV)^{2} / var(IV) / (var(D) - cov(D, IV)^{2} / var(IV))

v_Source1 <- sum((Y - Source1 * D)^2 * IV^2) / sum(IV * D)^2

IV <- MX %*% IV_all[, 4]
Y <- c(MX %*% Y_all)
D <- c(MX %*% D_all)
Source2 <- sum(IV * Y) / sum(IV * D)

Strength_Source2 <- n * cov(D, IV)^{2} / var(IV) / (var(D) - cov(D, IV)^{2} / var(IV))

v_Source2 <- sum((Y - Source2 * D)^2 * IV^2) / sum(IV * D)^2

MX <- list()
for (k in 1:K_split) {
  MX[[k]] <- diag(1, sum(split_ind == k)) - X[split_ind == k, ] %*% solve(t(X[split_ind == k, ]) %*% X[split_ind == k, ]) %*% t(X[split_ind == k, ])
}

for (k in 1:K_split) {
  IV <- MX[[k]] %*% IV_all[split_ind == k, ]
  Y <- c(MX[[k]] %*% Y_all[split_ind == k])
  D <- c(MX[[k]] %*% D_all[split_ind == k])

  TargetLasso_tmp[k] <- sum(IV[, 1] * Y) / sum(IV[, 1] * D)
  FTLasso1_tmp[k] <- sum(IV[, 3] * Y) / sum(IV[, 3] * D)
  FTLasso2_tmp[k] <- sum(IV[, 5] * Y) / sum(IV[, 5] * D)

  v_TargetLasso_tmp[k] <- sum((Y - TargetLasso_tmp[k] * D)^2 * IV[, 1]^2) / sum(IV[, 1] * D)^2
  v_FTLasso1_tmp[k] <- sum((Y - FTLasso1_tmp[k] * D)^2 * IV[, 3]^2) / sum(IV[, 3] * D)^2
  v_FTLasso2_tmp[k] <- sum((Y - FTLasso2_tmp[k] * D)^2 * IV[, 5]^2) / sum(IV[, 5] * D)^2

  Strength_TargetLasso_tmp[k] <- n * cov(D, IV[, 1])^{2} / var(IV[, 1]) / (var(D) - cov(D, IV[, 1])^{2} / var(IV[, 1]))
  Strength_FTLasso1_tmp[k] <- n * cov(D, IV[, 3])^{2} / var(IV[, 3]) / (var(D) - cov(D, IV[, 3])^{2} / var(IV[, 3]))
  Strength_FTLasso2_tmp[k] <- n * cov(D, IV[, 5])^{2} / var(IV[, 5]) / (var(D) - cov(D, IV[, 5])^{2} / var(IV[, 5]))
}

TargetLasso <- mean(TargetLasso_tmp)
FTLasso1 <- mean(FTLasso1_tmp)
FTLasso2 <- mean(FTLasso2_tmp)

v_TargetLasso <- mean(v_TargetLasso_tmp) / 3
v_FTLasso1 <- mean(v_FTLasso1_tmp) / 3
v_FTLasso2 <- mean(v_FTLasso2_tmp) / 3

Strength_TargetLasso <- mean(Strength_TargetLasso_tmp)
Strength_FTLasso1 <- mean(Strength_FTLasso1_tmp)
Strength_FTLasso2 <- mean(Strength_FTLasso2_tmp)

d <- ncol(IV_all)

w_cf <- c()
w_cf_tmp <- c()

for (k in 1:K_split) {
  for (j in 1:(K_split - 1)) {
    ind <- (1:K_split)[-k]

    IV <- MX[[ind[j]]] %*% IV_all[split_ind == ind[j], ]
    Y <- c(MX[[ind[j]]] %*% Y_all[split_ind == ind[j]])
    D <- c(MX[[ind[j]]] %*% D_all[split_ind == ind[j]])

    w_cf_tmp <- cbind(w_cf_tmp, nlminb(rep(0, d), iv_lasso_objective, D = D, IV = IV, penalty = 0.5 * sqrt(n))$par)
  }
  w_cf <- cbind(w_cf, rowMeans(w_cf_tmp))
}

for (k in 1:K_split) {
  IV <- MX[[k]] %*% IV_all[split_ind == k, ]
  Y <- c(MX[[k]] %*% Y_all[split_ind == k])
  D <- c(MX[[k]] %*% D_all[split_ind == k])

  CF_tmp[k] <- sum(IV %*% w_cf[, k] * Y) / sum(IV %*% w_cf[, k] * D)
}

CF <- mean(CF_tmp)

INIT <- max(min(10, CF), -10)

for (k in 1:K_split) {
  X_tmp <- X[split_ind == k, ]
  eps2 <- c(MX[[k]] %*%
              (Y_all[split_ind == k] - D_all[split_ind == k] * INIT))^2

  IV <- MX[[k]] %*% IV_all[split_ind == k, ]
  Y <- c(MX[[k]] %*% Y_all[split_ind == k])
  D <- c(MX[[k]] %*% D_all[split_ind == k])

  eta_h <- colSums(IV * D)
  V_h <-  t(IV) %*% (c(eps2) * IV)
  kappa0 <- (t(eta_h) %*% ginv(V_h) %*% eta_h / n) / ncol(IV)

  ev_min <- c()
  for (j in 1:50) {
    ev_min <- c(ev_min, min(Re(eigen(iv_variance_hessian(runif(d, - 0.25, 0.25), IV_raw = IV_all[split_ind == k, ], X = X_tmp, D_raw = D_all[split_ind == k], eps2 = eps2))$values)))
  }
  lambda <- max(0, min(0.25, 0.25 / (n * kappa0)) -  min(ev_min) / 2)

  b_opt <- nlminb(rep(0, d), iv_variance_objective,
    IV_raw = IV_all[split_ind == k, ], X = X_tmp,
    D_raw = D_all[split_ind == k], eps2 = eps2, lambda = lambda,
    projection_order = "right", lower = -0.25, upper = 0.25)$par

  q_opt <- c((exp(IV_all[split_ind == k, ] %*% b_opt) < 20) * exp(IV_all[split_ind == k, ] %*% b_opt) +
               (exp(IV_all[split_ind == k, ] %*% b_opt) >= 20) * 20)

  MX_opt <- diag(1, sum(split_ind == k)) - X_tmp %*% solve(t(X_tmp) %*% (q_opt * X_tmp)) %*% t(X_tmp) %*% diag(q_opt)
  IV <- MX_opt %*% IV_all[split_ind == k, ]
  Y <- c(MX_opt %*% Y_all[split_ind == k])
  D <- c(MX_opt %*% D_all[split_ind == k])

  V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

  H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)
  eta_h <- colSums(IV * D * q_opt)
  V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

  gamma_h <- colMeans(D * q_opt * IV)

  Strength_KPI_tmp[k] <- t(gamma_h) %*% ginv(sum(split_ind == k)^{-1} * t(IV) %*% (q_opt * IV)) %*% gamma_h

  Strength_KPI_tmp[k] <- n * Strength_KPI_tmp[k] / (mean(D^{2} * q_opt) - Strength_KPI_tmp[k])

  v_KPI_tmp[k] <- (t(eta_h) %*% ginv(V_h) %*% eta_h)^{-1}

  H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt) %*%
    (diag(1, sum(split_ind == k)) - X[split_ind == k, ] %*% solve(t(X[split_ind == k, ]) %*% (q_opt * X[split_ind == k, ])) %*% t(X[split_ind == k, ]) %*% diag(q_opt))

  beta_D_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, lower = - 5, upper = 5)$par

  KPI_tmp[k] <- beta_D_h
}

KPI <- mean(KPI_tmp)

v_KPI <- mean(v_KPI_tmp) / 3

Strength_KPI <- mean(Strength_KPI_tmp)

Res <- cbind(c(TargetLasso, Source1, FTLasso1, Source2, FTLasso2, KPI),
             sqrt(c(v_TargetLasso, v_Source1, v_FTLasso1, v_Source2, v_FTLasso2, v_KPI)))

Res <- cbind(Res, 2 * (1 - pnorm(abs(Res[, 1] / Res[, 2]))), Res[, 1] - qnorm(0.975) * Res[, 2],
                  Res[, 1] + qnorm(0.975) * Res[, 2], c(Strength_TargetLasso, Strength_Source1, Strength_FTLasso1,
                                                        Strength_Source2, Strength_FTLasso2, Strength_KPI))

rownames(Res) <- c("Target", "Source1", "FT1", "Source2", "FT2", "KPI")
colnames(Res) <- c("EST", "SD", "pval", "lb", "ub", "IV Strength")

save(Res, file = "Res_BMI_AFR.Rdata")
round(Res, 2)
sum(pval <= 1e-2)
