iv_functions_file <- normalizePath("../../R/iv_functions.R", mustWork = TRUE)
source(iv_functions_file, local = TRUE)

load("linear_1.Rdata")
load("linear_2.Rdata")
load("linear_3.Rdata")
load("linear_4.Rdata")
load("linear_5.Rdata")

library(MASS)

n <- 350

split <- cbind(1:350, 351:700, 701:1050)
gamma <- c(rep(0.5, 10), - rep(0.5, 10))

mod <- "M7"

#' Run one replicate of the Linear-sim_MM_CF analysis.
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
  KPI_opt0_tmp <- c()
  KPI_opt1_tmp <- c()
  KPI_opt2_tmp <- c()
  KPI_opt3_tmp <- c()
  KPI_opt4_tmp <- c()
  CF_tmp <- c()

  name <- paste("Results", arg, sep = "")
  DATA <- eval(parse(text = name))

  IV_all <- cbind(DATA[[sim]]$Source1,
                  DATA[[sim]]$Target1,
                  DATA[[sim]]$Trans1,
                  DATA[[sim]]$Source2,
                  DATA[[sim]]$Trans2,
                  DATA[[sim]]$Source3,
                  DATA[[sim]]$Trans3
  )

  D_ZX <- DATA[[sim]]$D_given_ZX

  ZX <- matrix(unlist(DATA[[sim]][13:117]), 3 * n)
  D_all <- DATA[[sim]]$D
  Y_all <- DATA[[sim]]$Y + exp(ZX %*% c(gamma, rep(0, 85))) * rnorm(3 * n, 0, 1)

  X <- cbind(1, ZX[, 101:105])
  MX <- diag(1, 3 * n) - X %*% solve(t(X) %*% X) %*% t(X)

  S0 <- c(2, 3, 5, 7)
  S_all <- 1:7

  if (mod == "M4") {
    IV_all <- IV_all[, S0]
  } else {
    IV_all <- IV_all[, S_all]
  }

  d <- ncol(IV_all)
  MX <- list()
  for (k in 1:3) {
    X <- cbind(1, ZX[split[, k], 101:105])
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

      w_cf_tmp <- cbind(w_cf_tmp, nlminb(rep(0, d), iv_lasso_objective, D = D, IV = IV, penalty = 0.5 * sqrt(n), lower = 0.01)$par)
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

    beta_h <- nlminb(INIT, iv_gel_objective, D = D, Y = Y, IV = IV, lower = -10, upper =10)$par

    GEL_tmp[k] <- beta_h

    YD <- cbind(Y, D)

    alpha <- eigen(solve(t(YD) %*% YD) %*% (t(YD) %*% H_uni %*% YD - t(YD) %*% (diag(H_uni) * YD)))$value[2]
    alpha <- (alpha - (1 - alpha) / n) / (1 - (1 - alpha) / n)
    HFUL_tmp[k] <- solve(t(D) %*% H_uni %*% D - t(D) %*% (diag(H_uni) * D) - alpha * t(D) %*% D) %*%
      (t(D) %*% H_uni %*% Y - t(diag(H_uni) * D) %*% Y - alpha * t(D) %*% Y)

  }

  for (k in 1:3) {
    X <- cbind(1, ZX[split[, k], 101:105])
    eps2 <- c(MX[[k]] %*%
                (Y_all[split[, k]] - D_all[split[, k]] * INIT))^2

    ev_min <- c()
    for (j in 1:50) {
      ev_min <- c(ev_min, min(Re(eigen(iv_variance_hessian(runif(d, - 0.25, 0.25), IV_raw = IV_all[split[, k], ], X = X, D_raw = D_all[split[, k]], eps2 = eps2))$values)))
    }
    lambda <- max(0,  0.25 - min(ev_min) / 2)

    b_opt <- nlminb(rep(0, d), iv_variance_objective, IV_raw = IV_all[split[, k], ], X = X, D_raw = D_all[split[, k]], eps2 = eps2, lambda = lambda, lower = - 0.25, upper = 0.25)$par

    q_opt <- c((exp(IV_all[split[, k], ] %*% b_opt) < 20) * exp(IV_all[split[, k], ] %*% b_opt) +
                 (exp(IV_all[split[, k], ] %*% b_opt) >= 20) * 20)

    MX_opt <- diag(1, n) - X %*% solve(t(X) %*% (q_opt * X)) %*% t(X) %*% diag(q_opt)
    IV <- MX_opt %*% IV_all[split[, k], ]
    Y <- c(MX_opt %*% Y_all[split[, k]])
    D <- c(MX_opt %*% D_all[split[, k]])

    V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

    H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)

    KPI_opt0_tmp[k] <- c(t(D) %*% H_opt %*% Y / (t(D) %*% H_opt %*% D))

    beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, lower = - 10, upper = 10)$par

    KPI_opt1_tmp[k] <- beta_h

    rho_h <- sum(diag(H_opt) * D * (Y - INIT * D))

    beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, rho_h = rho_h, lower = - 10, upper = 10)$par

    KPI_opt2_tmp[k] <- beta_h

    beta_h <- ((t(D) %*% H_opt %*% Y) - sum(diag(H_opt) * D * Y)) /
      ((t(D) %*% H_opt %*% D) - sum(diag(H_opt) * D^{2}))

    KPI_opt3_tmp[k] <- beta_h

    rho_h <- sum(diag(H_opt) * D * (Y - INIT * D))

    beta_h <- (t(D) %*% H_opt %*% Y - rho_h) / (t(D) %*% H_opt %*% D)

    KPI_opt4_tmp[k] <- beta_h
  }

  TSLS <- mean(TSLS_tmp)

  JIV <- mean(JIV_tmp)

  rJIV <- mean(rJIV_tmp)

  GEL <- mean(GEL_tmp)

  HFUL <- mean(HFUL_tmp)

  orDEEM <- mean(orDEEM_tmp)

  KPI_opt0 <- mean(KPI_opt0_tmp)

  KPI_opt1 <- mean(KPI_opt1_tmp)

  KPI_opt2 <- mean(KPI_opt2_tmp)

  KPI_opt3 <- mean(KPI_opt3_tmp)

  KPI_opt4 <- mean(KPI_opt4_tmp)

  c(TSLS, JIV, rJIV, GEL, HFUL, orDEEM, CF, KPI_opt0, KPI_opt3, KPI_opt4, KPI_opt1, KPI_opt2)
}

library(parallel)

arg = "1"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "mod", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res1 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res1 <- t(Res1)[, -c(3, 5, 6, 9, 10)]

Res1 <- pmin(pmax(Res1, -10), 10)

colnames(Res1) <-  c("TSLS", "JIV", "GEL", "Initial", "Plug-in", "KPI", "mKPI")

arg = "2"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "mod", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res2 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res2 <- t(Res2)[, -c(3, 5, 6, 9, 10)]

Res2 <- pmin(pmax(Res2, -10), 10)

colnames(Res1) <-  c("TSLS", "JIV", "GEL", "Initial", "Plug-in", "KPI", "mKPI")

arg = "3"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "mod", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res3 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res3 <- t(Res3)[, -c(3, 5, 6, 9, 10)]

Res3 <- pmin(pmax(Res3, -10), 10)

colnames(Res3) <-  c("TSLS", "JIV", "GEL", "Initial", "Plug-in", "KPI", "mKPI")

arg = "4"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "mod", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res4 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res4 <- t(Res4)[, -c(3, 5, 6, 9, 10)]

Res4 <- pmin(pmax(Res4, -10), 10)

colnames(Res4) <-  c("TSLS", "JIV", "GEL", "Initial", "Plug-in", "KPI", "mKPI")

arg = "5"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "mod", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res5 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res5 <- t(Res5)[, -c(3, 5, 6, 9, 10)]

Res5 <- pmin(pmax(Res5, -10), 10)

colnames(Res5) <-  c("TSLS", "JIV", "GEL", "Initial", "Plug-in", "KPI", "mKPI")

beta <- 0.5

mean_all <- c(apply(Res1, 2, abr_mean), apply(Res2, 2, abr_mean), apply(Res3, 2, abr_mean), apply(Res4, 2, abr_mean), apply(Res5, 2, abr_mean))
se_all <- c(apply(Res1, 2, abr_sd), apply(Res2, 2, abr_sd), apply(Res3, 2, abr_sd), apply(Res4, 2, abr_sd), apply(Res5, 2, abr_sd))
rmse_all <- c(apply(Res1, 2, abr_rmse, beta = beta), apply(Res2, 2, abr_rmse, beta = beta), apply(Res3, 2, abr_rmse, beta = beta), apply(Res4, 2, abr_rmse, beta = beta), apply(Res5, 2, abr_rmse, beta = beta))

mean_all <- matrix(mean_all, 7)
se_all <- matrix(se_all, 7)
rmse_all <- matrix(rmse_all, 7)

library(ggplot2)
library(latex2exp)
library(reshape2)
library(viridis)

S <- rep(c(0, 0.25, 0.5, 0.75, 1), 3)

dat = data.frame(S, rbind(abs(t(mean_all - beta)), t(se_all), t(rmse_all)))
names(dat) = c("Shift","TSLS", "JIV", "GEL", "Initial", "Plug-in", "KPI", "mKPI")
dat$cri = factor(rep(c("Bias", "SE", "RMSE"), rep(5, 3)), levels = c("Bias", "SE", "RMSE"))
dat = melt(dat,variable.name="Method",value.name = "Res", id = c("Shift", "cri"))
head(dat)

dat1 <- dat[dat$cri == "Bias", ]
a1 = ggplot(dat1)+
  geom_bar(
    mapping = aes(x = Shift,y = Res, fill = Method),
    position = 'dodge',
    stat = 'identity') +

  scale_colour_viridis_d(
    alpha = 1,
    begin = 1,
    end = 0.0,
    direction = 1,
    option = "D",
    aesthetics = "fill"
  ) +
  coord_cartesian(ylim=c(0, 0.4)) +

  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "h",y = "", title = "")+
  ggtitle("Bias") +

  guides(fill = guide_legend(reverse = F, nrow = 1))+
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "bottom",
        legend.key.size=unit(1,'cm'),
        strip.text = element_text(size=20))

dat2 <- dat[dat$cri == "SE", ]
a2 = ggplot(dat2)+
  geom_bar(
    mapping = aes(x = Shift,y = Res, fill = Method),
    position = 'dodge',
    stat = 'identity') +

  scale_colour_viridis_d(
    alpha = 1,
    begin = 1,
    end = 0.0,
    direction = 1,
    option = "D",
    aesthetics = "fill"
  ) +
  coord_cartesian(ylim=c(0, 0.8)) +

  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "h",y = "", title = "")+
  ggtitle("SE") +

  guides(fill = guide_legend(reverse = F, nrow = 1))+
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "bottom",
        legend.key.size=unit(1,'cm'),
        strip.text = element_text(size=20))

dat3 <- dat[dat$cri == "RMSE", ]
a3 = ggplot(dat3)+
  geom_bar(
    mapping = aes(x = Shift,y = Res, fill = Method),
    position = 'dodge',
    stat = 'identity') +

  scale_colour_viridis_d(
    alpha = 1,
    begin = 1,
    end = 0.0,
    direction = 1,
    option = "D",
    aesthetics = "fill"
  ) +
  coord_cartesian(ylim=c(0, 0.8)) +

  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "h",y = "", title = "")+
  ggtitle("RMSE") +

  guides(fill = guide_legend(reverse = F, nrow = 1))+
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "bottom",
        legend.key.size=unit(1,'cm'),
        strip.text = element_text(size=20))

library(patchwork)

layoutplot <- "
bbbbccccdddd
"

plotlist <- list(b= a1, c=a2, d=a3)

wrap_plots(plotlist, guides = 'collect', nrow = 1, design = layoutplot) & theme(legend.position = 'bottom',  legend.key.size=unit(1.5,'cm'),
                                                                                strip.text = element_text(size=25))
ggsave(paste(paste("Trend", "Comb", mod, sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 20, height = 10, units = "in")
