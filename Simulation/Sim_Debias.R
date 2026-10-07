iv_functions_file <- normalizePath("../R/iv_functions.R", mustWork = TRUE)
source(iv_functions_file, local = TRUE)

library(MASS)

d <- 10

n <- 500

split <- cbind(1:n, (n + 1):(2 * n), (2 * n + 1):(3 * n))

gamma <- rep(0.1, d) / sqrt(d)

#' Run one replicate of the Sim_Debias analysis.
#'
#' Uses the scenario data and constants configured below the imports.
#' @param sim Replicate index, also used as the random seed.
#' @return Numeric vector of causal estimates in the order assembled at the end.
simone <- function(sim) {
  set.seed(sim)
  CF_tmp <- c()

  IV_all <- mvrnorm(3 * n, rep(0, d), 0.1 * diag(1, d) + 0.9 * matrix(1, d, d))
  U_all <- rnorm(3 * n)
  D_all <- IV_all %*% gamma - sqrt(2 + rho) * U_all + sqrt(3 - rho) * rnorm(3 * n, 0, 1)
  Y_all <- D_all + sqrt(2 + rho) * U_all + sqrt(3 - rho) * rnorm(3 * n, 0, 1)

  d <- ncol(IV_all)
  MX <- list()
  for (k in 1:3) {
    X <- matrix(1, n, 1)
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

  X <- matrix(1, 3 * n, 1)
  MX <- diag(1, 3 * n) - X %*% solve(t(X) %*% X) %*% t(X)

  eps2 <- c(MX %*%
              (Y_all - D_all * INIT))^2

  IV <- MX %*% IV_all[, ]
  Y <- c(MX %*% Y_all[])
  D <- c(MX %*% D_all[])

  eta_h <- colSums(IV * D)
  V_h <-  t(IV) %*% (c(eps2) * IV)
  kappa0 <- (t(eta_h) %*% ginv(V_h) %*% eta_h / n) / ncol(IV)

  ev_min <- c()
  for (j in 1:50) {
    ev_min <- c(ev_min, min(Re(eigen(iv_variance_hessian(runif(d, - 0.25, 0.25), IV_raw = IV_all, X = X, D_raw = D_all, eps2 = eps2))$values)))
  }
  lambda <- max(0, min(0.25, 0.25 / (n * kappa0)) -  min(ev_min) / 2)

  b_opt <- nlminb(rep(0, d), iv_variance_objective, IV_raw = IV_all, X = X, D_raw = D_all, eps2 = eps2, lambda = lambda, lower = -0.25, upper = 0.25)$par

  q_opt <- c((exp(IV_all[, ] %*% b_opt) < 20) * exp(IV_all[, ] %*% b_opt) +
               (exp(IV_all[, ] %*% b_opt) >= 20) * 20)

  MX_opt <- diag(1, 3 * n) - X %*% solve(t(X) %*% (q_opt * X)) %*% t(X) %*% diag(q_opt)
  IV <- MX_opt %*% IV_all[, ]
  Y <- c(MX_opt %*% Y_all[])
  D <- c(MX_opt %*% D_all[])

  V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

  H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)

  KNIVES_opt0 <- c(t(D) %*% H_opt %*% Y / (t(D) %*% H_opt %*% D))

  beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, lower = - 10, upper = 10)$par

  KNIVES_opt1 <- beta_h

  r <- (t(D) %*% H_opt %*% D - sum(diag(H_opt) * D^{2})) / (t(D) %*% H_opt %*% D)

  KNIVES_opt2 <- r * KNIVES_opt1 + (1 - r) * INIT

  beta_h <- (t(D) %*% H_opt %*% Y - sum(diag(H_opt) * D * Y)) /
    (t(D) %*% H_opt %*% D - sum(diag(H_opt) * D^{2}))

  KNIVES_opt3 <- beta_h

  rho_h <- sum(diag(H_opt) * D * (Y - INIT * D))

  beta_h <- (t(D) %*% H_opt %*% Y - rho_h) / (t(D) %*% H_opt %*% D)

  KNIVES_opt4 <- beta_h

  c(INIT, KNIVES_opt0, KNIVES_opt3, KNIVES_opt4, KNIVES_opt1, KNIVES_opt2)
}

library(parallel)

beta <- 1

Res_debias <- array(NA, c(6, 1000, 3))

rho <- 1

cl <- makeCluster(20)
clusterExport(cl, c("n", "split", "gamma", "d", "rho"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res_debias[ , , 1] <- parSapply(cl, 1:1000, simone, simplify = "array")
stopCluster(cl)

rho <- 2

cl <- makeCluster(20)
clusterExport(cl, c("n", "split", "gamma", "d", "rho"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res_debias[ , , 2] <- parSapply(cl, 1:1000, simone, simplify = "array")
stopCluster(cl)

rho <- 3

cl <- makeCluster(20)
clusterExport(cl, c("n", "split", "gamma", "d", "rho"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res_debias[ , , 3] <- parSapply(cl, 1:1000, simone, simplify = "array")
stopCluster(cl)

rownames(Res_debias) <-  c("Init", "Plug-in", "DRE", "SUB", "KNIVES", "mKNIVES")

save(Res_debias,file = "Res_debiased.Rdata")

mean_all <- apply(Res_debias[c(2, 3, 5), ,], c(1, 3), abr_mean)

se_all <- apply(Res_debias[c(2, 3, 5), ,], c(1, 3), abr_sd)

rmse_all <- apply(Res_debias[c(2, 3, 5), ,], c(1, 3), abr_rmse, beta = beta)

library(ggplot2)
library(latex2exp)
library(reshape2)
library(viridis)

S <- rep(c(1, 2, 3), 3)

dat = data.frame(S, rbind(abs(t(mean_all - beta)), t(se_all), t(rmse_all)))
names(dat) = c("b", "Plug-in", "DRE", "KNIVES")
dat$cri = factor(rep(c("Bias", "SE", "RMSE"), rep(3, 3)), levels = c("Bias", "SE", "RMSE"))
dat = melt(dat,variable.name="Method",value.name = "Res", id = c("b", "cri"))
head(dat)

color <- c("#CC0000", "#006633", "#000066")

shape <- c(5, 6, 17)

dat1 <- dat[dat$cri == "Bias", ]
a1 = ggplot(dat1, aes(x = b,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = b,y = Res, color = Method), linewidth = 3, linetype = 2)+
  geom_point(aes(x = b,y = Res, color = Method, shape = Method), size = 7.5)+
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_colour_manual(values = color) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  labs(x = "b",y = "", title = "") +
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
        legend.key.size=unit(1.5,'cm'),
        strip.text = element_text(size=20))

dat2 <- dat[dat$cri == "SE", ]
a2 = ggplot(dat2, aes(x = b,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = b,y = Res, color = Method), linewidth = 3, linetype = 2)+
  geom_point(aes(x = b,y = Res, color = Method, shape = Method), size = 7.5)+
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_colour_manual(values = color) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  labs(x = "b",y = "", title = "")+
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
        legend.key.size=unit(1.5,'cm'),
        strip.text = element_text(size=20))

dat3 <- dat[dat$cri == "RMSE", ]
a3 = ggplot(dat3, aes(x = b,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = b,y = Res, color = Method), linewidth = 3, linetype = 2)+
  geom_point(aes(x = b,y = Res, color = Method, shape = Method), size = 7.5)+
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_colour_manual(values = color) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  labs(x = "b",y = "", title = "")+
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
        legend.key.size=unit(1.5,'cm'),
        strip.text = element_text(size=20))

library(patchwork)

layoutplot <- "
bbbbccccdddd
"

plotlist <- list(b= a1, c=a2, d=a3)

wrap_plots(plotlist, guides = 'collect', nrow = 1, design = layoutplot) & theme(legend.position = 'bottom',  legend.key.size=unit(3,'cm'),
                                                                                strip.text = element_text(size=25))
ggsave(paste(paste("Trend", "debias", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 20, height = 10, units = "in")

mean_all <- apply(Res_debias[ - 6, ,], c(1, 3), abr_mean)

se_all <- apply(Res_debias[ - 6, ,], c(1, 3), abr_sd)

rmse_all <- apply(Res_debias[ - 6, ,], c(1, 3), abr_rmse, beta = beta)

S <- rep(c(1, 2, 3), 3)

dat = data.frame(S, rbind(abs(t(mean_all - beta)), t(se_all), t(rmse_all)))
names(dat) = c("b", "Init", "Plug-in", "DRE", "SUB", "KNIVES")
dat$cri = factor(rep(c("Bias", "SE", "RMSE"), rep(3, 3)), levels = c("Bias", "SE", "RMSE"))
dat = melt(dat,variable.name="Method",value.name = "Res", id = c("b", "cri"))
head(dat)

color <- viridis(5, begin = 0.95,
                 end = 0.05, option = "H")
shape <- c(2, 5, 6, 15, 17)

dat1 <- dat[dat$cri == "Bias", ]
a1 = ggplot(dat1, aes(x = b,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = b,y = Res, color = Method), linewidth = 3, linetype = 2)+
  geom_point(aes(x = b,y = Res, color = Method, shape = Method), size = 7.5)+
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_colour_manual(values = color) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  labs(x = "b",y = "", title = "") +
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
        legend.key.size=unit(1.5,'cm'),
        strip.text = element_text(size=20))

dat2 <- dat[dat$cri == "SE", ]
a2 = ggplot(dat2, aes(x = b,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = b,y = Res, color = Method), linewidth = 3, linetype = 2)+
  geom_point(aes(x = b,y = Res, color = Method, shape = Method), size = 7.5)+
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_colour_manual(values = color) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  labs(x = "b",y = "", title = "")+
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
        legend.key.size=unit(1.5,'cm'),
        strip.text = element_text(size=20))

dat3 <- dat[dat$cri == "RMSE", ]
a3 = ggplot(dat3, aes(x = b,y = Res, color = Method, shape = Method))+
  geom_line(aes(x = b,y = Res, color = Method), linewidth = 3, linetype = 2)+
  geom_point(aes(x = b,y = Res, color = Method, shape = Method), size = 7.5)+
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = 3) +
  scale_colour_manual(values = color) +
  scale_shape_manual(values=shape)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  labs(x = "b",y = "", title = "")+
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
        legend.key.size=unit(1.5,'cm'),
        strip.text = element_text(size=20))

layoutplot <- "
bbbbccccdddd
"

plotlist <- list(b= a1, c=a2, d=a3)

wrap_plots(plotlist, guides = 'collect', nrow = 1, design = layoutplot) & theme(legend.position = 'bottom',  legend.key.size=unit(3,'cm'),
                                                                                strip.text = element_text(size=25))
ggsave(paste(paste("Trend", "debias", "all", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 20, height = 10, units = "in")
