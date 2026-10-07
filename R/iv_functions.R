# Shared IV objectives and simulation summaries. No analysis runs on source().

#' Evaluate the lasso objective for an instrument combination.
#'
#' @param w Numeric instrument-weight vector.
#' @param D Residualized exposure vector.
#' @param IV Residualized instrument matrix, with observations in rows.
#' @param penalty Multiplier of the L1 penalty, set by the calling analysis.
#' @return Scalar penalized residual sum of squares.
iv_lasso_objective <- function(w, D, IV, penalty) {
  sum((D - IV %*% w)^2) + penalty * sum(abs(w))
}

#' Evaluate the variance Hessian used to choose the ridge penalty.
#'
#' @param b Log-weight coefficients, one per instrument.
#' @param IV_raw Instrument matrix before covariate adjustment.
#' @param X Covariate design matrix, including the intercept.
#' @param D_raw Exposure vector before covariate adjustment.
#' @param eps2 Squared residuals computed using the initial causal estimate.
#' @return Square Hessian matrix with one row and column per instrument.
iv_variance_hessian <- function(b, IV_raw, X, D_raw, eps2) {
  q <- c(exp(IV_raw %*% b))
  cX_q <- solve(t(X) %*% (q * X)) %*% t(X) %*% diag(q)
  IV <- IV_raw - X %*% cX_q %*% IV_raw
  D <- c(D_raw - X %*% cX_q %*% D_raw)

  eta_h <- colSums(IV * D * q)
  V_h <- t(IV) %*% (c(eps2 * q^2) * IV)
  h0 <- c(t(eta_h) %*% MASS::ginv(V_h) %*% eta_h)
  h1 <- c(IV %*% MASS::ginv(V_h) %*% eta_h)
  vec_tmp1 <- colSums((D * h1 - q * eps2 * h1^2) * q * IV)
  vec_tmp2 <- q * (2 * q * eps2 * h1 - D)
  Mat_tmp <- t(IV) %*% (vec_tmp2 * IV)
  8 / h0^3 * vec_tmp1 %*% t(vec_tmp1) - 1 / h0^2 * (
    2 * t(IV) %*% (q * (D * h1 - 2 * q * eps2 * h1^2) * IV) +
      2 * Mat_tmp %*% MASS::ginv(V_h) %*% Mat_tmp)
}

#' Evaluate the penalized variance objective for instrument weighting.
#'
#' @param b Log-weight coefficients, one per instrument.
#' @param IV_raw Instrument matrix before covariate adjustment.
#' @param X Covariate design matrix, including the intercept.
#' @param D_raw Exposure vector before covariate adjustment.
#' @param eps2 Squared residuals computed using the initial causal estimate.
#' @param lambda Ridge penalty selected by the calling analysis.
#' @param projection_order Matrix multiplication order: "left" for simulations,
#'   "right" for the original real-data calculations. Preserves rounding behavior.
#' @return One-by-one matrix containing the penalized variance criterion.
iv_variance_objective <- function(b, IV_raw, X, D_raw, eps2, lambda,
                                  projection_order = "left") {
  projection_order <- match.arg(projection_order, c("left", "right"))
  q <- c(exp(IV_raw %*% b))
  cX_q <- solve(t(X) %*% (q * X)) %*% t(X) %*% diag(q)
  if (projection_order == "right") {
    IV <- IV_raw - X %*% (cX_q %*% IV_raw)
    D <- c(D_raw - X %*% (cX_q %*% D_raw))
  } else {
    IV <- IV_raw - X %*% cX_q %*% IV_raw
    D <- c(D_raw - X %*% cX_q %*% D_raw)
  }
  eta_h <- colSums(IV * D * q)
  V_h <- t(IV) %*% (c(eps2 * q^2) * IV)
  (t(eta_h) %*% MASS::ginv(V_h) %*% eta_h)^(-1) + lambda * sum(b^2)
}

#' Evaluate the squared bias-corrected IV estimating equation.
#'
#' @param b Candidate causal effect.
#' @param D Residualized exposure vector.
#' @param Y Residualized outcome vector.
#' @param H Weighted instrument projection matrix.
#' @param rho_h Fixed bias correction; NULL recomputes it at b.
#' @return One-by-one matrix containing the squared estimating equation.
iv_kpi_objective <- function(b, D, Y, H, rho_h = NULL) {
  if (is.null(rho_h)) {
    rho_h <- sum(diag(H) * D * (Y - b * D))
  }
  (t(D) %*% H %*% (Y - D * b) -
     rho_h / sum(diag(H) * (Y - b * D)^2) *
     t(Y - D * b) %*% H %*% (Y - D * b))^2
}

#' Evaluate the generalized empirical likelihood criterion.
#'
#' @param b Candidate causal effect.
#' @param D Residualized exposure vector.
#' @param Y Residualized outcome vector.
#' @param IV Residualized instrument matrix.
#' @return One-by-one matrix containing the moment-based criterion.
iv_gel_objective <- function(b, D, Y, IV) {
  g <- (Y - D * b) * IV
  t(colSums(g)) %*% MASS::ginv(t(g) %*% g) %*% colSums(g)
}

#' Evaluate the squared oracle DEEM estimating equation.
#'
#' @param beta Candidate causal effect.
#' @param gamma_h Estimated instrument-exposure coefficients.
#' @param Gamma_h Estimated instrument-outcome coefficients.
#' @param rho_h Estimated exposure-outcome error covariance adjustment.
#' @param D_d Diagonal exposure covariance adjustment matrix.
#' @param D_y Diagonal outcome covariance adjustment matrix.
#' @param V_DEEM Instrument weighting matrix.
#' @return One-by-one matrix containing the squared estimating equation.
iv_deem_objective <- function(beta, gamma_h, Gamma_h, rho_h, D_d, D_y, V_DEEM) {
  g <- t(gamma_h - (rho_h - beta) * D_d %*%
           solve(D_y + (beta^2 - 2 * beta * rho_h) * D_d) %*%
           (Gamma_h - beta * gamma_h)) %*% V_DEEM %*% (Gamma_h - beta * gamma_h)
  g^2
}

#' Retain estimates strictly within the original 1.5-IQR fences.
#'
#' @param x Numeric simulation estimates without missing values.
#' @return Subvector inside the fences; constant vectors yield an empty vector.
trim_iqr <- function(x) {
  IQR <- quantile(x, 0.75) - quantile(x, 0.25)
  x[(x < (quantile(x, 0.75) + 1.5 * IQR)) &
      (x > (quantile(x, 0.25) - 1.5 * IQR))]
}

#' Calculate the mean after the original IQR trimming.
#' @param x Numeric simulation estimates without missing values.
#' @return Scalar trimmed mean; NaN when no estimates remain.
abr_mean <- function(x) {
  mean(trim_iqr(x))
}

#' Calculate the standard deviation after the original IQR trimming.
#' @param x Numeric simulation estimates without missing values.
#' @return Scalar sample standard deviation; NA with fewer than two retained estimates.
abr_sd <- function(x) {
  sd(trim_iqr(x))
}

#' Calculate root mean squared error after the original IQR trimming.
#' @param x Numeric simulation estimates without missing values.
#' @param beta True causal effect used as the error reference.
#' @return Scalar trimmed RMSE; NaN when no estimates remain.
abr_rmse <- function(x, beta) {
  sqrt(mean((trim_iqr(x) - beta)^2))
}
