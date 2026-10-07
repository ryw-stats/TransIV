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

#' Run one replicate of the Linear-sim_Target_Transfer analysis.
#'
#' Uses the scenario data and constants configured below the imports.
#' @param sim Replicate index, also used as the random seed.
#' @return Numeric vector of causal estimates in the order assembled at the end.
simone <- function(sim) {
  set.seed(sim)
  Target_tmp <- c()
  Trans1_tmp <- c()
  Target2_tmp <- c()
  Trans2_tmp <- c()
  Target3_tmp <- c()
  Trans3_tmp <- c()
  CF_tmp <- c()
  KPI_opt_tmp <- c()

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

  Source1 <- sum(IV_all[, 1] * MX %*% Y_all) / sum(IV_all[, 1] * MX %*% D_all)
  Source2 <- sum(IV_all[, 4] * MX %*% Y_all) / sum(IV_all[, 4] * MX %*% D_all)
  Source3 <- sum(IV_all[, 6] * MX %*% Y_all) / sum(IV_all[, 6] * MX %*% D_all)

  MX <- list()
  for (k in 1:3) {
    X <- cbind(1, ZX[split[, k], 101:105])
    MX[[k]] <- diag(1, n) - X %*% solve(t(X) %*% X) %*% t(X)
  }

  for (k in 1:3) {
    IV <- MX[[k]] %*% IV_all[split[, k], ]
    Y <- c(MX[[k]] %*% Y_all[split[, k]])
    D <- c(MX[[k]] %*% D_all[split[, k]])

    Target_tmp[k] <- sum(IV[, 2] * Y) / sum(IV[, 2] * D)
    Trans1_tmp[k] <- sum(IV[, 3] * Y) / sum(IV[, 3] * D)
    Trans2_tmp[k] <- sum(IV[, 5] * Y) / sum(IV[, 5] * D)
    Trans3_tmp[k] <- sum(IV[, 7] * Y) / sum(IV[, 7] * D)
  }

  d <- ncol(IV_all)

  w_cf <- c()
  w_cf_tmp <- c()

  for (k in 1:3) {
    for (j in 1:2) {
      ind <- (1:3)[-k]

      IV <- MX[[ind[j]]] %*% IV_all[split[, ind[j]], ]
      Y <- c(MX[[ind[j]]] %*% Y_all[split[, ind[j]]])
      D <- c(MX[[ind[j]]] %*% D_all[split[, ind[j]]])

      w_cf_tmp <- cbind(w_cf_tmp, nlminb(rep(0, d), iv_lasso_objective, D = D, IV = IV, penalty = sqrt(n), lower = 0.01)$par)
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
    X <- cbind(1, ZX[split[, k], 101:105])
    eps2 <- c(MX[[k]] %*%
                (Y_all[split[, k]] - D_all[split[, k]] * INIT))^2

    ev_min <- c()
    for (j in 1:50) {
      ev_min <- c(ev_min, min(Re(eigen(iv_variance_hessian(runif(d, -0.25, 0.25), IV_raw = IV_all[split[, k], ], X = X, D_raw = D_all[split[, k]], eps2 = eps2))$values)))
    }
    lambda <- max(0, 0.25 - min(ev_min) / 2)

    b_opt <- nlminb(rep(0, d), iv_variance_objective, IV_raw = IV_all[split[, k], ], X = X, D_raw = D_all[split[, k]], eps2 = eps2, lambda = lambda, lower = -0.25, upper = 0.25)$par

    q_opt <- c((exp(IV_all[split[, k], ] %*% b_opt) < 20) * exp(IV_all[split[, k], ] %*% b_opt) +
                 (exp(IV_all[split[, k], ] %*% b_opt) >= 20) * 20)

    MX_opt <- diag(1, n) - X %*% solve(t(X) %*% (q_opt * X)) %*% t(X) %*% diag(q_opt)
    IV <- MX_opt %*% IV_all[split[, k], ]
    Y <- c(MX_opt %*% Y_all[split[, k]])
    D <- c(MX_opt %*% D_all[split[, k]])

    V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

    H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)

    rho_h <- sum(diag(H_opt) * D * (Y - INIT * D))

    beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, rho_h = rho_h, lower = - 10, upper = 10)$par

    KPI_opt_tmp[k] <- beta_h
  }

  Target <- mean(Target_tmp)
  Trans1 <- mean(Trans1_tmp)
  Trans2 <- mean(Trans2_tmp)
  Trans3 <- mean(Trans3_tmp)
  KPI_opt <- mean(KPI_opt_tmp)

  c(Target, Source1, Source2, Source3, Trans1, Trans2, Trans3, KPI_opt, CF)
}

library(parallel)

arg = "1"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res1 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res1 <- t(Res1)[, -9]

Res1 <- pmin(pmax(Res1, -10), 10)

colnames(Res1) <-  c("Target", "Source 1", "Source 2", "Source 3", "Refit 1", "Refit 2", "Refit 3", "mKPI")

arg = "2"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res2 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res2 <- t(Res2)[, -9]

Res2 <- pmin(pmax(Res2, -10), 10)

colnames(Res2) <-  c("Target", "Source 1", "Source 2", "Source 3", "Refit 1", "Refit 2", "Refit 3", "mKPI")

arg = "3"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res3 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res3 <- t(Res3)[, -9]

Res3 <- pmin(pmax(Res3, -10), 10)

colnames(Res3) <-  c("Target", "Source 1", "Source 2", "Source 3", "Refit 1", "Refit 2", "Refit 3", "mKPI")

arg = "4"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res4 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res4 <- t(Res4)[, -9]

Res4 <- pmin(pmax(Res4, -10), 10)

colnames(Res4) <-  c("Target", "Source 1", "Source 2", "Source 3", "Refit 1", "Refit 2", "Refit 3", "mKPI")

arg = "5"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res5 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res5 <- t(Res5)[, -9]

Res5 <- pmin(pmax(Res5, -10), 10)

colnames(Res5) <-  c("Target", "Source 1", "Source 2", "Source 3", "Refit 1", "Refit 2", "Refit 3", "mKPI")

beta <- 0.5

mean_all <- c(apply(Res1, 2, abr_mean), apply(Res2, 2, abr_mean), apply(Res3, 2, abr_mean), apply(Res4, 2, abr_mean), apply(Res5, 2, abr_mean))
se_all <- c(apply(Res1, 2, abr_sd), apply(Res2, 2, abr_sd), apply(Res3, 2, abr_sd), apply(Res4, 2, abr_sd), apply(Res5, 2, abr_sd))
rmse_all <- c(apply(Res1, 2, abr_rmse, beta = beta), apply(Res2, 2, abr_rmse, beta = beta), apply(Res3, 2, abr_rmse, beta = beta), apply(Res4, 2, abr_rmse, beta = beta), apply(Res5, 2, abr_rmse, beta = beta))

mean_all <- matrix(mean_all, 8)
se_all <- matrix(se_all, 8)
rmse_all <- matrix(rmse_all, 8)

library(ggplot2)
library(latex2exp)
library(reshape2)
library(viridis)

S <- rep(c(0, 0.25, 0.5, 0.75, 1), 3)

dat = data.frame(S, rbind(abs(t(mean_all - beta)), t(se_all), t(rmse_all)))
names(dat) = c("Shift", "Target", "Source 1", "Source 2", "Source 3", "Refit 1", "Refit 2", "Refit 3", "mKPI")
dat$cri = factor(rep(c("Bias", "SE", "RMSE"), rep(5, 3)), levels = c("Bias", "SE", "RMSE"))
dat = melt(dat,variable.name="Method",value.name = "Res", id = c("Shift", "cri"))
head(dat)

shape <- c(0:2, 5:9)

dat1 <- dat[dat$cri == "Bias", ]
a1 = ggplot(dat1, aes(x = Shift,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = Shift,y = Res, color = Method), linewidth = 1.5, linetype = 2)+
  geom_point(aes(x = Shift,y = Res, color = Method, shape = Method), size = 4)+
  scale_colour_viridis_d(
    alpha = 1,
    begin = 0.95,
    end = 0.05,
    direction = 1,
    option = "H",
    aesthetics = "colour"
  ) +
  coord_cartesian(ylim=c(0, 0.5)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "h",y = "", title = "")+
  ggtitle("Bias") +

  guides(colour = guide_legend(reverse = F, nrow = 1))+
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "bold"),
        legend.position = "bottom",
        legend.key.size=unit(1.2,'cm'),
        strip.text = element_text(size=20))

dat2 <- dat[dat$cri == "SE", ]
a2 = ggplot(dat2, aes(x = Shift,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = Shift,y = Res, color = Method), linewidth = 1.5, linetype = 2)+
  geom_point(aes(x = Shift,y = Res, color = Method, shape = Method), size = 4)+
  scale_colour_viridis_d(
    alpha = 1,
    begin = 0.95,
    end = 0.05,
    direction = 1,
    option = "H",
    aesthetics = "colour"
  ) +
  coord_cartesian(ylim=c(0, 2)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "h",y = "", title = "")+
  ggtitle("SE") +

  guides(colour = guide_legend(reverse = F, nrow = 1))+
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "bold"),
        legend.position = "bottom",
        legend.key.size=unit(1.2,'cm'),
        strip.text = element_text(size=20))

dat3 <- dat[dat$cri == "RMSE", ]
a3 = ggplot(dat3, aes(x = Shift,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = Shift,y = Res, color = Method), linewidth = 1.5, linetype = 2)+
  geom_point(aes(x = Shift,y = Res, color = Method, shape = Method), size = 4)+
  scale_colour_viridis_d(
    alpha = 1,
    begin = 0.95,
    end = 0.05,
    direction = 1,
    option = "H",
    aesthetics = "colour"
  ) +
  coord_cartesian(ylim=c(0, 2)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "h",y = "", title = "")+
  ggtitle("RMSE") +

  guides(colour = guide_legend(reverse = F, nrow = 1))+
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "bold"),
        legend.position = "bottom",
        legend.key.size=unit(1.2,'cm'),
        strip.text = element_text(size=20))

library(patchwork)

layoutplot <- "
bbbbccccdddd
"

plotlist <- list(b= a1, c=a2, d=a3)

wrap_plots(plotlist, guides = 'collect', nrow = 1, design = layoutplot) & theme(legend.position = 'bottom',  legend.key.size=unit(1.5,'cm'),
                                                                                strip.text = element_text(size=25))
ggsave(paste(paste("Trend", "sg", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 26, height = 12, units = "in")
