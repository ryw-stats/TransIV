iv_functions_file <- normalizePath("../../R/iv_functions.R", mustWork = TRUE)
source(iv_functions_file, local = TRUE)

load("nonlinear_new.Rdata")

Results_new <- Results_list_best_new_med

library(MASS)

n <- 350

split <- cbind(1:350, 351:700, 701:1050)
gamma <- c(rep(0.4, 5), rep(-0.4, 5), rep(0, 9))

#' Run one replicate of the Nonlinear-sim_Target_Transfer analysis.
#'
#' Uses the scenario data and constants configured below the imports.
#' @param sim Replicate index, also used as the random seed.
#' @return Numeric vector of causal estimates in the order assembled at the end.
simone <- function(sim) {
  set.seed(sim)
  Target_Lasso_tmp <- c()
  Target_RF_tmp <- c()
  Target_GBM_tmp <- c()
  TransLasso_tmp <- c()
  TransRF_tmp <- c()
  TransGBM_tmp <- c()
  TransNN_wide_tmp <- c()
  TransNN_deep_tmp <- c()
  CF_tmp <- c()
  KNIVES_opt_tmp <- c()
  v_KNIVES_tmp <- c()

  name <- paste("Results", arg, sep = "_")
  DATA <- eval(parse(text = name))

  IV_all <- cbind(DATA[[sim]]$InitialPredictions,
                  DATA[[sim]]$InitialPredictions_rf,
                  DATA[[sim]]$InitialPredictions_gbm,
                  DATA[[sim]]$InitialPredictions_lasso,
                  DATA[[sim]]$TargetPredictions,
                  DATA[[sim]]$RFPredictions,
                  DATA[[sim]]$GBMPredictions,
                  DATA[[sim]]$CombinedPredictions1,
                  DATA[[sim]]$CombinedPredictions2,
                  DATA[[sim]]$CombinedPredictions3,
                  DATA[[sim]]$CombinedPredictions4,
                  DATA[[sim]]$CombinedPredictions5
  )

  D_ZX <- DATA[[sim]]$D_given_ZX

  ZX <- matrix(unlist(DATA[[sim]][17:35]), 3 * n)
  Y_all <- DATA[[sim]]$Y + exp(ZX %*% gamma) * rnorm(3 * n, 0, 1)
  D_all <- DATA[[sim]]$D

  X <- cbind(1, ZX[, 15:19])
  MX <- diag(1, 3 * n) - X %*% solve(t(X) %*% X) %*% t(X)
  Source_NN <- sum(IV_all[, 1] * MX %*% Y_all) / sum(IV_all[, 1] * MX %*% D_all)
  Source_RF <- sum(IV_all[, 2] * MX %*% Y_all) / sum(IV_all[, 2] * MX %*% D_all)
  Source_GBM <- sum(IV_all[, 3] * MX %*% Y_all) / sum(IV_all[, 3] * MX %*% D_all)
  Source_Lasso <- sum(IV_all[, 4] * MX %*% Y_all) / sum(IV_all[, 4] * MX %*% D_all)

  MX <- list()
  for (k in 1:3) {
    X <- cbind(1, ZX[split[, k], 15:19])
    MX[[k]] <- diag(1, n) - X %*% solve(t(X) %*% X) %*% t(X)
  }

  for (k in 1:3) {
    IV <- MX[[k]] %*% IV_all[split[, k], ]
    Y <- c(MX[[k]] %*% Y_all[split[, k]])
    D <- c(MX[[k]] %*% D_all[split[, k]])

    Target_Lasso_tmp[k] <- sum(IV[, 5] * Y) / sum(IV[, 5] * D)
    Target_RF_tmp[k] <- sum(IV[, 6] * Y) / sum(IV[, 6] * D)
    Target_GBM_tmp[k] <- sum(IV[, 7] * Y) / sum(IV[, 7] * D)
    TransLasso_tmp[k] <- sum(IV[, 8] * Y) / sum(IV[, 8] * D)
    TransRF_tmp[k] <- sum(IV[, 9] * Y) / sum(IV[, 9] * D)
    TransGBM_tmp[k] <- sum(IV[, 10] * Y) / sum(IV[, 10] * D)
    TransNN_wide_tmp[k] <- sum(IV[, 11] * Y) / sum(IV[, 11] * D)
    TransNN_deep_tmp[k] <- sum(IV[, 12] * Y) / sum(IV[, 12] * D)
  }
  IV_all <- IV_all[, -11]
  d <- ncol(IV_all)

  w_cf <- c()
  w_cf_tmp <- c()

  for (k in 1:3) {
    for (j in 1:2) {
      ind <- (1:3)[-k]

      IV <- MX[[ind[j]]] %*% IV_all[split[, ind[j]], ]
      Y <- c(MX[[ind[j]]] %*% Y_all[split[, ind[j]]])
      D <- c(MX[[ind[j]]] %*% D_all[split[, ind[j]]])

      w_cf_tmp <- cbind(w_cf_tmp, nlminb(rep(0, d), iv_lasso_objective, D = D, IV = IV, penalty = 0.5 * sqrt(n))$par)
    }
    w_cf <- cbind(w_cf, rowMeans(w_cf_tmp))
  }

  for (k in 1:3) {
    IV <- MX[[k]] %*% IV_all[split[, k], ]
    Y <- c(MX[[k]] %*% Y_all[split[, k]])
    D <- c(MX[[k]] %*% D_all[split[, k]])

    CF_tmp[k] <- sum(IV %*% w_cf[, k] * Y) / sum(IV %*% w_cf[, k] * D)
  }

  CF <- mean(CF_tmp)

  INIT <- max(min(10, CF), -10)

  for (k in 1:3) {
    X <- cbind(1, ZX[split[, k], 15:19])
    eps2 <- c(MX[[k]] %*%
                (Y_all[split[, k]] - D_all[split[, k]] * INIT))^2

    IV <- MX[[k]] %*% IV_all[split[, k], ]
    Y <- c(MX[[k]] %*% Y_all[split[, k]])
    D <- c(MX[[k]] %*% D_all[split[, k]])

    eta_h <- colSums(IV * D)
    V_h <-  t(IV) %*% (c(eps2) * IV)
    kappa0 <- (t(eta_h) %*% ginv(V_h) %*% eta_h / n) / ncol(IV)

    ev_min <- c()
    for (j in 1:50) {
      ev_min <- c(ev_min, min(Re(eigen(iv_variance_hessian(runif(d, - 0.25, 0.25), IV_raw = IV_all[split[, k], ], X = X, D_raw = D_all[split[, k]], eps2 = eps2))$values)))
    }
    lambda <- max(0, min(0.25, 0.25 / (n * kappa0)) -  min(ev_min) / 2)

    b_opt <- nlminb(rep(0, d), iv_variance_objective, IV_raw = IV_all[split[, k], ], X = X, D_raw = D_all[split[, k]], eps2 = eps2, lambda = lambda, lower = -0.25, upper = 0.25)$par

    q_opt <- c((exp(IV_all[split[, k], ] %*% b_opt) < 20) * exp(IV_all[split[, k], ] %*% b_opt) +
                 (exp(IV_all[split[, k], ] %*% b_opt) >= 20) * 20)

    MX_opt <- diag(1, n) - X %*% solve(t(X) %*% (q_opt * X)) %*% t(X) %*% diag(q_opt)
    IV <- MX_opt %*% IV_all[split[, k], ]
    Y <- c(MX_opt %*% Y_all[split[, k]])
    D <- c(MX_opt %*% D_all[split[, k]])

    eta_h <- colSums(IV * D * q_opt)
    V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

    H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)

    v_KNIVES_tmp[k] <- (t(eta_h) %*% ginv(V_h) %*% eta_h)^{-1}

    beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, lower = - 10, upper = 10)$par

    KNIVES_opt_tmp[k] <- beta_h
  }

  Target_Lasso <- mean(Target_Lasso_tmp)
  Target_RF <- mean(Target_RF_tmp)
  Target_GBM <- mean(Target_GBM_tmp)
  TransLasso <- mean(TransLasso_tmp)
  TransRF <- mean(TransRF_tmp)
  TransGBM <- mean(TransGBM_tmp)
  TransNN_wide <- mean(TransNN_wide_tmp)
  TransNN_deep <- mean(TransNN_deep_tmp)
  KNIVES_opt <- mean(KNIVES_opt_tmp)
  v_KNIVES <- mean(v_KNIVES_tmp) / 3

  c(Target_Lasso, Target_RF, Target_GBM, Source_Lasso, Source_NN, Source_RF, Source_GBM, TransLasso, TransNN_deep, TransRF, TransGBM, KNIVES_opt, CF)
}

library(parallel)

arg = "new"

cl <- makeCluster(20)
clusterExport(cl, c("n", "Results_new", "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res_new <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res_new <- t(Res_new)[, -13]

Res_new <- pmin(pmax(Res_new, -10), 10)

colnames(Res_new) <-  c("Target Lasso", "Target RF", "Target GBM", "Source Lasso", "Source NN", "Source RF", "Source GBM", "Fine-Tuned Lasso", "Fine-Tuned NN", "Fine-Tuned RF", "Fine-Tuned GBM",  "KNIVES")

beta <- 2

mean_all <- c(apply(Res_new, 2, abr_mean))
se_all <- c(apply(Res_new, 2, abr_sd))

library(ggplot2)
library(latex2exp)
library(viridis)
Res_all <- data.frame(c(1:12), c(mean_all), c(se_all), c(mean_all) - c(se_all), c(mean_all) + c(se_all),
                       factor(rep(c("Target Lasso", "Target RF", "Target GBM", "Source Lasso", "Source NN", "Source RF", "Source GBM", "Fine-Tuned Lasso", "Fine-Tuned NN", "Fine-Tuned RF", "Fine-Tuned GBM",  "KNIVES"), 1),
                              levels = c("Target Lasso", "Target RF", "Target GBM", "Source Lasso", "Source NN", "Source RF", "Source GBM", "Fine-Tuned Lasso", "Fine-Tuned NN", "Fine-Tuned RF", "Fine-Tuned GBM",  "KNIVES")))
colnames(Res_all) <- c("X", "Mean", "SD", "LB", "UB", "Method")

ggplot(Res_all) +
  geom_errorbar(aes(x = X, ymin = LB, ymax = UB), linetype = "longdash", cex = 2.5, width = 0.3) +
  geom_point(aes(x= X, y = Mean, color = Method), shape = 19, size = 12) +
  scale_colour_viridis_d(
    alpha = 1,
    begin = 0.95,
    end = 0.05,
    direction = 1,
    option = "H",
    aesthetics = "colour"
  ) +
  geom_hline(yintercept = beta, linetype = 2, color = "darkgrey", size = 2.5) +
  ylim(c(1, 3)) +
  ylab("Estimate") +
  xlab("Method") +
  scale_x_continuous(breaks = 1:12,
                     labels = c("Target Lasso", "Target RF", "Target GBM", "Source Lasso", "Source NN", "Source RF", "Source GBM", "Fine-Tuned Lasso", "Fine-Tuned NN", "Fine-Tuned RF", "Fine-Tuned GBM",  "KNIVES")) +
  theme_light() +
  guides(color = "none") +
  theme(axis.text.x = element_text(size = 35, angle = 60, vjust = 0.55),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 40, face = "bold"),
        axis.title.x = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 35, face = "bold"),
        legend.text = element_text(size = 35, face = "bold"),
        legend.key.size=unit(2,'cm'),
        legend.position = "top")

ggsave(paste(paste("Mean-SD", "SG", sep = "_"), ".pdf", sep = ""),
       path = "./fig", device = "pdf", width = 20, height = 10, units = "in")
