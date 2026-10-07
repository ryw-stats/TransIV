iv_functions_file <- normalizePath("../../R/iv_functions.R", mustWork = TRUE)
source(iv_functions_file, local = TRUE)

load("linear_3.Rdata")

library(MASS)

n <- 350

split <- cbind(1:350, 351:700, 701:1050)
gamma <- c(rep(0.5, 10), - rep(0.5, 10))

#' Run one replicate of the Illustrative Example analysis.
#'
#' Uses the scenario data and constants configured below the imports.
#' @param sim Replicate index, also used as the random seed.
#' @return Numeric vector of causal estimates in the order assembled at the end.
simone <- function(sim) {
  set.seed(sim)
  Target_tmp <- c()
  CF_tmp <- c()
  PlugIn_tmp <- c()
  KNIVES_opt_tmp <- c()

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

    V_h <-  t(IV) %*% (c(eps2 * q_opt^{2}) * IV)

    H_opt <- diag(q_opt) %*% IV %*% ginv(V_h) %*% t(IV) %*% diag(q_opt)

    PlugIn_tmp[k] <- (t(D) %*% H_opt %*% Y) / (t(D) %*% H_opt %*% D)

    beta_h <- nlminb(INIT, iv_kpi_objective, D = D, Y = Y, H = H_opt, lower = - 10, upper = 10)$par

    KNIVES_opt_tmp[k] <- beta_h
  }

  Target <- mean(Target_tmp)
  PlugIn <- mean(PlugIn_tmp)
  KNIVES_opt <- mean(KNIVES_opt_tmp)

  c(Target, PlugIn, KNIVES_opt)
}

library(parallel)

arg = "3"

cl <- makeCluster(20)
clusterExport(cl, c("n", paste("Results", arg, sep = ""), "split", "arg", "gamma"))
clusterEvalQ(cl, library("MASS"))
clusterCall(cl, source, iv_functions_file)
Res1 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

Res1 <- data.frame(t(Res1))

beta <- 0.5

Res1 <- Res1[apply(abs(Res1 - beta), 1, max) < 4, ]

colnames(Res1) <-  c("Target", "PlugIn", "KNIVES")

library(ggplot2)

ggplot(Res1, aes(x=Target)) + geom_density(position = "identity", alpha = 0.2, linewidth = 1.5, bw = 0.3, color = "#000066", fill = "#000066") +
  labs(x = "Estimate",y = "Density", title = "") +
  geom_vline(xintercept=beta, color = "#333333",
             linetype="dashed", linewidth = 1.5) +
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35, face = "bold"),
        axis.title.y = element_text(size = 35, face = "bold"),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "bold"),
        legend.position = "bottom",
        legend.key.size=unit(1.2,'cm'),
        strip.text = element_text(size=20))

ggsave("Illustrate_target.pdf",
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")

Res2 <- data.frame(Est = c(Res1[, 1], Res1[, 2], Res1[, 3]), Method = factor(rep(c("Target-based", "Plug-in", "KNIVES"), each = nrow(Res1)),
                                                         levels = c("Target-based", "Plug-in", "KNIVES")))

ggplot(Res2, aes(x=Est, color = Method, fill = Method)) +
  geom_density(position = "identity", alpha = 0.2, linewidth = 1.5, bw = 0.3) +
  scale_color_manual(values = c("#000066", "#006600", "#CC0000")) +
  scale_fill_manual(values = c("#000066", "#006600", "#CC0000")) +
  labs(x = "Estimate",y = "Density", title = "") +
  geom_vline(xintercept=beta, color = "#333333",
             linetype="dashed", linewidth = 1.5) +
  theme_light() +
  theme(plot.title = element_text(size = 35,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 35, face = "bold"),
        axis.title.y = element_text(size = 35, face = "bold"),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30, face = "bold"),
        legend.position = "top",
        legend.key.size=unit(1.2,'cm'),
        legend.spacing.x = unit(1.0, 'cm'),
        strip.text = element_text(size=20))

ggsave("Illustrate_plug-in.pdf",
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")
