iv_functions_file <- normalizePath("../../R/iv_functions.R", mustWork = TRUE)
source(iv_functions_file, local = TRUE)

load("nonlinear_new.Rdata")

Results_new <- Results_list_best_new_med

library(MASS)

n <- 350

split <- cbind(1:350, 351:700, 701:1050)
gamma <- c(rep(0.4, 5), - rep(0.4, 5), rep(0, 9))

mod <- "M11"

#' Run one replicate of the Nonlinear-sim_MM_CF analysis.
#'
#' Uses the scenario data and constants configured below the imports.
#' @param sim Replicate index, also used as the random seed.
#' @return Numeric vector of causal estimates in the order assembled at the end.
simone <- function(sim) {
  set.seed(sim)
  TSLS_tmp <- c()
  JIV_tmp <- c()
  rJIV_tmp <- c()
  GEL_tmp <- c()
  HFUL_tmp <- c()
  orDEEM_tmp <- c()
  KNIVES_opt0_tmp <- c()
  KNIVES_opt1_tmp <- c()
  KNIVES_opt2_tmp <- c()
  KNIVES_opt3_tmp <- c()
  KNIVES_opt4_tmp <- c()
  v_tmp <- c()
  CF_tmp <- c()

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
                  DATA[[sim]]$CombinedPredictions5
  )

  D_ZX <- DATA[[sim]]$D_given_ZX

  ZX <- matrix(unlist(DATA[[sim]][17:35]), 3 * n)
  Y_all <- DATA[[sim]]$Y + exp(ZX %*% gamma) * rnorm(3 * n, 0, 1)
  D_all <- DATA[[sim]]$D

  S0 <- c(1, 7, 8)
  S_all <- 1:11

  if (mod == "M3") {
    IV_all <- IV_all[, S0]
  } else {
    IV_all <- IV_all[, S_all]
  }

  d <- ncol(IV_all)
  MX <- list()
  for (k in 1:3) {
    X <- cbind(1, ZX[split[, k], 15:19])
    MX[[k]] <- diag(1, n) - X %*% solve(t(X) %*% X) %*% t(X)
  }

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
    IV <- MX[[k]] %*% IV_all[split[, k], ]
    Y <- c(MX[[k]] %*% Y_all[split[, k]])
    D <- c(MX[[k]] %*% D_all[split[, k]])

    gamma_h <- c(t(IV) %*% D / colSums(IV^{2}))
    Gamma_h <- c(t(IV) %*% Y / colSums(IV^{2}))

    Sig_d <- diag(colMeans(IV^{2})^{-1}) %*% var(t((D - t(IV) * gamma_h) * t(IV))) %*% diag(colMeans(IV^{2})^{-1}) / n
    Sig_y <- diag(colMeans(IV^{2})^{-1}) %*% var(t((Y - t(IV) * Gamma_h) * t(IV))) %*% diag(colMeans(IV^{2})^{-1}) / n

    V_DEEM <- ginv(Sig_y)

    D_d <- diag(diag(V_DEEM %*% Sig_d) / diag(V_DEEM))
    D_y <- diag(diag(V_DEEM %*% Sig_y) / diag(V_DEEM))

    rho_h <- c((t(gamma_h) %*% V_DEEM %*% (Gamma_h - INIT * gamma_h) + INIT * sum(diag(V_DEEM %*% D_d))) /
                 sum(diag(V_DEEM %*% D_d)))

    orDEEM_tmp[k] <- nlminb(INIT, iv_deem_objective, gamma_h = gamma_h, Gamma_h = Gamma_h,
      rho_h = rho_h, D_d = D_d, D_y = D_y, V_DEEM = V_DEEM, lower = - 10, upper = 10)$par
  }

  for (k in 1:3) {
    IV <- MX[[k]] %*% IV_all[split[, k], ]
    Y <- c(MX[[k]] %*% Y_all[split[, k]])
    D <- c(MX[[k]] %*% D_all[split[, k]])

    H_uni <- IV %*% ginv(t(IV) %*% IV) %*% t(IV)

    beta_TSLS <- c(t(D) %*% H_uni %*% Y / t(D) %*% H_uni %*% D)
    TSLS_tmp[k] <- beta_TSLS

    h_IV <- diag(H_uni)

    D_hat <- H_uni %*% D
    e_IV <- D - D_hat
    D_JIV <- D_hat - h_IV / (1 - h_IV) * e_IV

    JIV_tmp[k] <- sum(D_JIV * Y) / sum(D_JIV * D)

    DE <- c(D, rep(0, d))
    IVE <- rbind(IV, diag(rep(sqrt(n / 10), d)))
    H_IVE <- IVE %*% ginv(t(IVE) %*% IVE) %*% t(IVE)

    h_IVE <- diag(H_IVE)

    D_hatE <- H_IVE %*% DE
    e_IVE <- DE - D_hatE
    D_rJIV <- D_hatE[-((n + 1):(n + d))] - h_IVE[-((n + 1):(n + d))] / (1 - h_IVE[-((n + 1):(n + d))]) * e_IVE[-((n + 1):(n + d))]

    rJIV_tmp[k] <- sum(D_rJIV * Y) / sum(D_rJIV * D)

    beta_h <- nlminb(INIT, iv_gel_objective, D = D, Y = Y, IV = IV, lower = -10, upper = 10)$par

    GEL_tmp[k] <- beta_h

    YD <- cbind(Y, D)

    alpha <- eigen(solve(t(YD) %*% YD) %*% (t(YD) %*% H_uni %*% YD - t(YD) %*% (diag(H_uni) * YD)))$value[2]
    alpha <- (alpha - (1 - alpha) / n) / (1 - (1 - alpha) / n)
    HFUL_tmp[k] <- solve(t(D) %*% H_uni %*% D - t(D) %*% (diag(H_uni) * D) - alpha * t(D) %*% D) %*%
      (t(D) %*% H_uni %*% Y - t(diag(H_uni) * D) %*% Y - alpha * t(D) %*% Y)

  }

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

    b_opt <- nlminb(rep(0, d), iv_variance_objective, IV_raw = IV_all[split[, k], ], X = X, D_raw = D_all[split[, k]], eps2 = eps2, lambda = lambda, lower = - 0.25, upper = 0.25)$par

    q_opt <- c((exp(IV_all[split[, k], ] %*% b_opt) < 20) * exp(IV_all[split[, k], ] %*% b_opt) +
                 (exp(IV_all[split[, k], ] %*% b_opt) >= 20) * 20)

    MX_opt <- diag(1, n) - X %*% solve(t(X) %*% (q_opt * X)) %*% t(X) %*% diag(q_opt)
    IV <- MX_opt %*% IV_all[split[, k], ]
    Y <- c(MX_opt %*% Y_all[split[, k]])
    D <- c(MX_opt %*% D_all[split[, k]])

    V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

    H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)

    H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)

    KNIVES_opt0_tmp[k] <- c(t(D) %*% H_opt %*% Y / (t(D) %*% H_opt %*% D))

    beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, lower = - 10, upper = 10)$par

    KNIVES_opt1_tmp[k] <- beta_h

    rho_h <- sum(diag(H_opt) * D * (Y - INIT * D))

    beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, rho_h = rho_h, lower = - 10, upper = 10)$par

    KNIVES_opt2_tmp[k] <- beta_h

    beta_h <- ((t(D) %*% H_opt %*% Y) - sum(diag(H_opt) * D * Y)) /
      ((t(D) %*% H_opt %*% D) - sum(diag(H_opt) * D^{2}))

    KNIVES_opt3_tmp[k] <- beta_h

    rho_h <- sum(diag(H_opt) * D * (Y - INIT * D))

    beta_h <- (t(D) %*% H_opt %*% Y - rho_h) / (t(D) %*% H_opt %*% D)

    KNIVES_opt4_tmp[k] <- beta_h

    gamma_h <- c((t(D) %*% diag(q_opt) %*% IV) / length(Y))

    Omega_h <- ginv(V_h / length(Y))

    kappa_h <- max(1e-5, t(gamma_h) %*% Omega_h %*% gamma_h - sum(diag(cov(D * q_opt * IV) %*% Omega_h)) / length(Y))

    e_h <- Y - KNIVES_opt1_tmp[k] * D

    c_h <- sum(diag(H_opt) * D * e_h) / sum(diag(H_opt) * e_h^2)

    M_tmp <- ((t(IV) %*% (q_opt^2 * e_h * (D - c_h * e_h) * IV)) / length(Y)) %*% Omega_h
    ev <- Re(eigen(Omega_h)$values)
    v_tmp[k] <- (kappa_h + max(0, sum(diag(Omega_h %*% var(q_opt * D * IV))) - c_h^2 * sum(ev > 0.05 * max(ev)) +
    sum(diag(M_tmp %*% M_tmp))) / length(Y)) / kappa_h^2
  }

  TSLS <- mean(TSLS_tmp)

  JIV <- mean(JIV_tmp)

  rJIV <- mean(rJIV_tmp)

  GEL <- mean(GEL_tmp)

  HFUL <- mean(HFUL_tmp)

  orDEEM <- mean(orDEEM_tmp)

  KNIVES_opt0 <- mean(KNIVES_opt0_tmp)

  KNIVES_opt1 <- mean(KNIVES_opt1_tmp)

  KNIVES_opt2 <- mean(KNIVES_opt2_tmp)

  KNIVES_opt3 <- mean(KNIVES_opt3_tmp)

  KNIVES_opt4 <- mean(KNIVES_opt4_tmp)

  v_KNIVES <- mean(v_tmp) / (3 * n)

  c(TSLS, JIV, rJIV, GEL, HFUL, orDEEM,  CF, KNIVES_opt0, KNIVES_opt3, KNIVES_opt4, KNIVES_opt1, KNIVES_opt2, v_KNIVES)
}

library(parallel)

arg = "new"

cl <- makeCluster(20)
clusterExport(cl, c("n", "mod", "Results_new", "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res_new <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res_new <- t(Res_new)[, c(1, 2, 4, 11, 13)]

Res_new <- pmin(pmax(Res_new, -10), 10)

colnames(Res_new) <-  c("TSLS", "JIV", "GEL", "KNIVES", "v_KNIVES")

beta <- 2

mean_all <- c(apply(Res_new, 2, abr_mean))[-5]
se_all <- c(apply(Res_new, 2, abr_sd))[-5]

library(ggplot2)
library(latex2exp)
library(viridis)
Res_all <- data.frame(c(4:1), c(mean_all), c(se_all), c(mean_all) - c(se_all), c(mean_all) + c(se_all),
                      factor(rep(c("TSLS", "JIV", "GEL", "KNIVES"), 1),
                             levels = c("TSLS", "JIV", "GEL", "KNIVES")))
colnames(Res_all) <- c("X", "Mean", "SD", "LB", "UB", "Method")

ggplot(Res_all) +
  geom_errorbar(aes(y = X, xmin = LB, xmax = UB), linetype = "longdash", cex = 2.5, width = 0.2) +
  geom_point(aes(x=Mean, y = X, color = Method), shape = 19, size = 12) +
  scale_colour_viridis_d(
    alpha = 1,
    begin = 0.9,
    end = 0.1,
    direction = 1,
    option = "H",
    aesthetics = "colour"
  ) +
  geom_vline(xintercept = beta, linetype = 2, color = "darkgrey", size = 2.5) +
  xlim(c(1, 3)) +
  ylab("Method") +
  xlab("Estimate") +
  scale_y_continuous(breaks = 4:1,
                     labels = c("TSLS", "JIV", "GEL", "KNIVES")) +
  theme_light() +
  guides(color = "none") +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 40, face = "bold"),
        axis.title.x = element_text(size = 40, face = "bold"),
        legend.title = element_text(size = 35, face = "bold"),
        legend.text = element_text(size = 35, face = "bold"),
        legend.key.size=unit(2,'cm'),
        legend.position = "top")

ggsave(paste(paste("Mean-SD", "Comb", mod, sep = "_"), ".pdf", sep = ""),
       path = "./fig", device = "pdf", width = 10, height = 6, units = "in")
