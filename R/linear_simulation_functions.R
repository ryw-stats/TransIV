# Shared functions; source this file from a workflow before calling them.
# Source model_functions.R first. Load the packages listed in Readme.txt.

#' Prepare raw linear simulation outputs for the IV analysis scripts.
#'
#' @param results List of replicate outputs from run_linear_simulation.
#' @return List of target data and predictions with the analysis field names.
prepare_linear_results <- function(results) {
  lapply(seq_along(results), function(i) {
    replicate <- results[[i]]
    if (!is.list(replicate) || length(replicate) != 2L ||
        !is.list(replicate[[1]]) || length(replicate[[1]]) != 117L) {
      stop("Replicate ", i, " is incomplete or has an unexpected format.")
    }
    data <- replicate[[1]]
    names(data)[4:12] <- c("Source1", "Trans1", "Target1", "Source2", "Trans2",
                          "Target2", "Source3", "Trans3", "Target3")
    data
  })
}

#' Generate linear source data and three target folds.
#'
#' @param n Number of target observations per held-out fold; total target size is 3 * n.
#' @param n_training Number of observations in each source training sample.
#' @param beta1 True exposure effect in the outcome model.
#' @param seed Random seed for this replicate.
#' @param random_vector Length-100 vector of target genotype coefficients.
#' @param selected_data_eu Source genotype matrix with 2,000 matched variant columns.
#' @param selected_data_af Target genotype matrix with the same 2,000 variant columns.
#' @param near_weight Weight on near-source coefficients; the far-source weight is 1 - near_weight.
#' @return List of three source samples, target training/held-out folds, and conditional exposure means.
simulate_linear_data <- function(n, n_training, beta1, seed, random_vector,
                                 selected_data_eu, selected_data_af, near_weight) {
  set.seed(seed)

  n_test=n*2
  total_n <- n + n_test
  nx <- 5

  # Sample source and target genotypes independently.
  sample_indices_source1 <- sample(1:nrow(selected_data_eu), n_training, replace = TRUE)
  sample_indices_source2 <- sample(1:nrow(selected_data_eu), n_training, replace = TRUE)
  sample_indices_source3 <- sample(1:nrow(selected_data_eu), n_training, replace = TRUE)

  sample_indices_target <- sample(1:nrow(selected_data_af), total_n, replace = TRUE)

  Z <- selected_data_af[sample_indices_target, ]
  Z_training1 <- selected_data_eu[sample_indices_source1, ]
  Z_training2 <- selected_data_eu[sample_indices_source2, ]
  Z_training3 <- selected_data_eu[sample_indices_source3, ]

  colnames(Z) <- NULL
  rownames(Z) <- NULL
  colnames(Z_training1) <- NULL
  rownames(Z_training1) <- NULL
  colnames(Z_training2) <- NULL
  rownames(Z_training2) <- NULL
  colnames(Z_training3) <- NULL
  rownames(Z_training3) <- NULL

  Z <- as.matrix(Z)
  Z_training1 <- as.matrix(Z_training1)
  Z_training2 <- as.matrix(Z_training2)
  Z_training3 <- as.matrix(Z_training3)

  # Generate source and target covariates.
  X <- mvrnorm(total_n, mu = rep(0, nx), Sigma = diag(nx))
  X_training1 <- mvrnorm(n_training, mu = rep(0, nx), Sigma = diag(nx))
  X_training2 <- mvrnorm(n_training, mu = rep(0, nx), Sigma = diag(nx))
  X_training3 <- mvrnorm(n_training, mu = rep(0, nx), Sigma = diag(nx))

  beta_Z <- random_vector

  beta_X <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  beta_X_Y <- c(0.1, 0.1, 0.1, 0.1, 0.1)

  # Mix near and far source coefficients using the notebook's original weight.

  C <- rnorm(total_n, mean = 0, sd = 1)
  D <- Z[, 1] * beta_Z[1] +
       Z[, 2] * beta_Z[2] +
       Z[, 3] * beta_Z[3] +
       Z[, 4] * beta_Z[4] +
       Z[, 5]  * beta_Z[5] +
       Z[, 6] * beta_Z[6] +
       Z[, 7] * beta_Z[7] +
       Z[, 8] * beta_Z[8] +
       Z[, 9] * beta_Z[9] +
       Z[, 10] * beta_Z[10] +
       Z[, 11] * beta_Z[11] +
       Z[, 12] * beta_Z[12] +
       Z[, 13] * beta_Z[13] +
       Z[, 14] * beta_Z[14] +
       Z[, 15]  * beta_Z[15] +
       Z[, 16] * beta_Z[16] +
       Z[, 17] * beta_Z[17] +
       Z[, 18] * beta_Z[18] +
       Z[, 19] * beta_Z[19] +
       Z[, 20] * beta_Z[20] +
       Z[, 21] * beta_Z[21] +
       Z[, 22] * beta_Z[22] +
       Z[, 23] * beta_Z[23] +
       Z[, 24] * beta_Z[24] +
       Z[, 25]  * beta_Z[25] +
       Z[, 26] * beta_Z[26] +
       Z[, 27] * beta_Z[27] +
       Z[, 28] * beta_Z[28] +
       Z[, 29] * beta_Z[29] +
       Z[, 30] * beta_Z[30] +
       Z[, 31] * beta_Z[31] +
       Z[, 32] * beta_Z[32] +
       Z[, 33] * beta_Z[33] +
       Z[, 34] * beta_Z[34] +
       Z[, 35]  * beta_Z[35] +
       Z[, 36] * beta_Z[36] +
       Z[, 37] * beta_Z[37] +
       Z[, 38] * beta_Z[38] +
       Z[, 39] * beta_Z[39] +
       Z[, 40] * beta_Z[40] +
       Z[, 41] * beta_Z[41] +
       Z[, 42] * beta_Z[42] +
       Z[, 43] * beta_Z[43] +
       Z[, 44] * beta_Z[44] +
       Z[, 45]  * beta_Z[45] +
       Z[, 46] * beta_Z[46] +
       Z[, 47] * beta_Z[47] +
       Z[, 48] * beta_Z[48] +
       Z[, 49] * beta_Z[49] +
       Z[, 50] * beta_Z[50] +
       Z[, 51] * beta_Z[51] +
       Z[, 52] * beta_Z[52] +
       Z[, 53] * beta_Z[53] +
       Z[, 54] * beta_Z[54] +
       Z[, 55]  * beta_Z[55] +
       Z[, 56] * beta_Z[56] +
       Z[, 57] * beta_Z[57] +
       Z[, 58] * beta_Z[58] +
       Z[, 59] * beta_Z[59] +
       Z[, 60] * beta_Z[60] +
       Z[, 61] * beta_Z[61] +
       Z[, 62] * beta_Z[62] +
       Z[, 63] * beta_Z[63] +
       Z[, 64] * beta_Z[64] +
       Z[, 65]  * beta_Z[65] +
       Z[, 66] * beta_Z[66] +
       Z[, 67] * beta_Z[67] +
       Z[, 68] * beta_Z[68] +
       Z[, 69] * beta_Z[69] +
       Z[, 70] * beta_Z[70] +
       Z[, 71] * beta_Z[71] +
       Z[, 72] * beta_Z[72] +
       Z[, 73] * beta_Z[73] +
       Z[, 74] * beta_Z[74] +
       Z[, 75]  * beta_Z[75] +
       Z[, 76] * beta_Z[76] +
       Z[, 77] * beta_Z[77] +
       Z[, 78] * beta_Z[78] +
       Z[, 79] * beta_Z[79] +
       Z[, 80] * beta_Z[80] +
       Z[, 81] * beta_Z[81] +
       Z[, 82] * beta_Z[82] +
       Z[, 83] * beta_Z[83] +
       Z[, 84] * beta_Z[84] +
       Z[, 85]  * beta_Z[85] +
       Z[, 86] * beta_Z[86] +
       Z[, 87] * beta_Z[87] +
       Z[, 88] * beta_Z[88] +
       Z[, 89] * beta_Z[89] +
       Z[, 90] * beta_Z[90] +
       Z[, 91] * beta_Z[91] +
       Z[, 92] * beta_Z[92] +
       Z[, 93] * beta_Z[93] +
       Z[, 94] * beta_Z[94] +
       Z[, 95]  * beta_Z[95] +
       Z[, 96] * beta_Z[96] +
       Z[, 97] * beta_Z[97] +
       Z[, 98] * beta_Z[98] +
       Z[, 99] * beta_Z[99] +
       Z[, 100] * beta_Z[100] +
       X[, 1] * beta_X[1] +
       X[, 2] * beta_X[2] +
       X[, 3] * beta_X[3] +
       X[, 4] * beta_X[4] +
       X[, 5] * beta_X[5] +
       2*C+
       rnorm(total_n, mean = 0, sd = 2)

  D_given_ZX <- Z[, 1] * beta_Z[1] +
       Z[, 2] * beta_Z[2] +
       Z[, 3] * beta_Z[3] +
       Z[, 4] * beta_Z[4] +
       Z[, 5]  * beta_Z[5] +
       Z[, 6] * beta_Z[6] +
       Z[, 7] * beta_Z[7] +
       Z[, 8] * beta_Z[8] +
       Z[, 9] * beta_Z[9] +
       Z[, 10] * beta_Z[10] +
       Z[, 11] * beta_Z[11] +
       Z[, 12] * beta_Z[12] +
       Z[, 13] * beta_Z[13] +
       Z[, 14] * beta_Z[14] +
       Z[, 15]  * beta_Z[15] +
       Z[, 16] * beta_Z[16] +
       Z[, 17] * beta_Z[17] +
       Z[, 18] * beta_Z[18] +
       Z[, 19] * beta_Z[19] +
       Z[, 20] * beta_Z[20] +
       Z[, 21] * beta_Z[21] +
       Z[, 22] * beta_Z[22] +
       Z[, 23] * beta_Z[23] +
       Z[, 24] * beta_Z[24] +
       Z[, 25]  * beta_Z[25] +
       Z[, 26] * beta_Z[26] +
       Z[, 27] * beta_Z[27] +
       Z[, 28] * beta_Z[28] +
       Z[, 29] * beta_Z[29] +
       Z[, 30] * beta_Z[30] +
       Z[, 31] * beta_Z[31] +
       Z[, 32] * beta_Z[32] +
       Z[, 33] * beta_Z[33] +
       Z[, 34] * beta_Z[34] +
       Z[, 35]  * beta_Z[35] +
       Z[, 36] * beta_Z[36] +
       Z[, 37] * beta_Z[37] +
       Z[, 38] * beta_Z[38] +
       Z[, 39] * beta_Z[39] +
       Z[, 40] * beta_Z[40] +
       Z[, 41] * beta_Z[41] +
       Z[, 42] * beta_Z[42] +
       Z[, 43] * beta_Z[43] +
       Z[, 44] * beta_Z[44] +
       Z[, 45]  * beta_Z[45] +
       Z[, 46] * beta_Z[46] +
       Z[, 47] * beta_Z[47] +
       Z[, 48] * beta_Z[48] +
       Z[, 49] * beta_Z[49] +
       Z[, 50] * beta_Z[50] +
       Z[, 51] * beta_Z[51] +
       Z[, 52] * beta_Z[52] +
       Z[, 53] * beta_Z[53] +
       Z[, 54] * beta_Z[54] +
       Z[, 55]  * beta_Z[55] +
       Z[, 56] * beta_Z[56] +
       Z[, 57] * beta_Z[57] +
       Z[, 58] * beta_Z[58] +
       Z[, 59] * beta_Z[59] +
       Z[, 60] * beta_Z[60] +
       Z[, 61] * beta_Z[61] +
       Z[, 62] * beta_Z[62] +
       Z[, 63] * beta_Z[63] +
       Z[, 64] * beta_Z[64] +
       Z[, 65]  * beta_Z[65] +
       Z[, 66] * beta_Z[66] +
       Z[, 67] * beta_Z[67] +
       Z[, 68] * beta_Z[68] +
       Z[, 69] * beta_Z[69] +
       Z[, 70] * beta_Z[70] +
       Z[, 71] * beta_Z[71] +
       Z[, 72] * beta_Z[72] +
       Z[, 73] * beta_Z[73] +
       Z[, 74] * beta_Z[74] +
       Z[, 75]  * beta_Z[75] +
       Z[, 76] * beta_Z[76] +
       Z[, 77] * beta_Z[77] +
       Z[, 78] * beta_Z[78] +
       Z[, 79] * beta_Z[79] +
       Z[, 80] * beta_Z[80] +
       Z[, 81] * beta_Z[81] +
       Z[, 82] * beta_Z[82] +
       Z[, 83] * beta_Z[83] +
       Z[, 84] * beta_Z[84] +
       Z[, 85]  * beta_Z[85] +
       Z[, 86] * beta_Z[86] +
       Z[, 87] * beta_Z[87] +
       Z[, 88] * beta_Z[88] +
       Z[, 89] * beta_Z[89] +
       Z[, 90] * beta_Z[90] +
       Z[, 91] * beta_Z[91] +
       Z[, 92] * beta_Z[92] +
       Z[, 93] * beta_Z[93] +
       Z[, 94] * beta_Z[94] +
       Z[, 95]  * beta_Z[95] +
       Z[, 96] * beta_Z[96] +
       Z[, 97] * beta_Z[97] +
       Z[, 98] * beta_Z[98] +
       Z[, 99] * beta_Z[99] +
       Z[, 100] * beta_Z[100] +
       X[, 1] * beta_X[1] +
       X[, 2] * beta_X[2] +
       X[, 3] * beta_X[3] +
       X[, 4] * beta_X[4] +
       X[, 5] * beta_X[5]

  Y <- beta1 * D +
       X[, 1] * beta_X_Y[1] +
       X[, 2] * beta_X_Y[2] +
       X[, 3] * beta_X_Y[3] +
       X[, 4] * beta_X_Y[4] +
       X[, 5] * beta_X_Y[5] +
       5*C +
       rnorm(total_n, mean = 0, sd = 0.5)

  shuffled_indices <-sample(1:total_n)

  split_size <- floor(total_n / 3)

  fold1 <- shuffled_indices[1:split_size]
  fold2 <- shuffled_indices[(split_size + 1):(2 * split_size)]
  fold3 <- shuffled_indices[(2 * split_size + 1):total_n]

  test_indices <- fold1
  train_indices <- setdiff(1:total_n, test_indices)

  Z1 <- Z[test_indices, ]
  X1 <- X[test_indices, ]
  D1 <- D[test_indices]
  Y1  <- Y[test_indices]
  D_given_ZX1 <- D_given_ZX[test_indices]

  Z_test <- Z[train_indices, ]
  X_test <- X[train_indices, ]
  D_test <- D[train_indices]
  Y_test  <- Y[train_indices]

  # Construct the second held-out fold.
  test_indices2 <- fold2
  train_indices2 <- setdiff(1:total_n, test_indices2)
  Z2 <- Z[test_indices2, ]
  X2 <- X[test_indices2, ]
  D2 <- D[test_indices2]
  Y2 <- Y[test_indices2]
  D_given_ZX2 <- D_given_ZX[test_indices2]

  Z_test_2 <- Z[train_indices2, ]
  X_test_2 <- X[train_indices2, ]
  D_test_2 <- D[train_indices2]
  Y_test_2 <- Y[train_indices2]

  # Construct the third held-out fold.
  test_indices3 <- fold3
  train_indices3 <- setdiff(1:total_n, test_indices3)
  Z3 <- Z[test_indices3, ]
  X3 <- X[test_indices3, ]
  D3 <- D[test_indices3]
  Y3 <- Y[test_indices3]
  D_given_ZX3 <- D_given_ZX[test_indices3]

  Z_test_3 <- Z[train_indices3, ]
  X_test_3 <- X[train_indices3, ]
  D_test_3 <- D[train_indices3]
  Y_test_3 <- Y[train_indices3]

  beta_Z_far1 <- beta_Z
  beta_Z_far1[c(5:20)] <- 0.1 * beta_Z_far1[c(5:20)]

  beta_Z_near1 <- beta_Z
  beta_Z_near1[17:20] <- 0.1 * beta_Z_near1[17:20]

  C_training1 <- rnorm(n_training, mean = 0, sd = 1)
  D_training1 <-
       Z_training1[, -(101:2000)] %*% (near_weight * beta_Z_near1 + (1 - near_weight) * beta_Z_far1) +
       X_training1 %*% beta_X +
       2*C_training1 +
       rnorm(n_training, mean = 0, sd = 2)

  beta_Z_far2 <- beta_Z
  beta_Z_far2[c(1:12, 17:20)] <- 0.1 * beta_Z_far2[c(1:12, 17:20)]

  beta_Z_near2 <- beta_Z
  beta_Z_near2[6:9] <- 0.1 * beta_Z_near2[6:9]

  C_training2 <- rnorm(n_training, mean = 0, sd = 1)
  D_training2 <-
       Z_training2[, -(101:2000)] %*% (near_weight * beta_Z_near2 + (1 - near_weight) * beta_Z_far2) +
       X_training2 %*% beta_X +
       2*C_training2 +
       rnorm(n_training, mean = 0, sd = 2)

  beta_Z_far3 <- beta_Z
  beta_Z_far3[c(1:5, 10:20)] <- 0.1 * beta_Z_far3[c(1:5, 10:20)]

  beta_Z_near3 <- beta_Z
  beta_Z_near3[12:15] <- 0.1 * beta_Z_near3[12:15]

  C_training3 <- rnorm(n_training, mean = 0, sd = 1)
  D_training3 <-
       Z_training3[, -(101:2000)] %*% (near_weight * beta_Z_near3 + (1 - near_weight) * beta_Z_far3) +
       X_training3 %*% beta_X +
       2*C_training3 +
       rnorm(n_training, mean = 0, sd = 2)

  return(list(
  Z1=Z1, X1=X1, D1=D1, Y1=Y1, D_given_ZX1=D_given_ZX1,
  Z_test_1=Z_test, X_test_1=X_test, D_test_1=D_test, Y_test_1=Y_test,

  Z2=Z2, X2=X2, D2=D2, Y2=Y2,D_given_ZX2=D_given_ZX2,
  Z_test_2=Z_test_2, X_test_2=X_test_2, D_test_2=D_test_2, Y_test_2=Y_test_2,

  Z3=Z3, X3=X3, D3=D3, Y3=Y3,D_given_ZX3=D_given_ZX3,
  Z_test_3=Z_test_3, X_test_3=X_test_3, D_test_3=D_test_3, Y_test_3=Y_test_3,

  Z_training1=Z_training1, Z_training2=Z_training2, Z_training3=Z_training3, X_training1=X_training1, X_training2=X_training2, X_training3=X_training3, D_training1=D_training1, D_training2=D_training2,D_training3=D_training3
))
}

#' Compute source-model residuals observation by observation.
#'
#' @param model_main Fitted source prediction model.
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @param lambda Lasso penalty used for prediction.
#' @return Numeric residual vector in observation order.
compute_residuals <- function(model_main, x0, y0, lambda) {

  residuals <- numeric(length(y0))

  for (i in 1:nrow(x0)) {

    X_i_0 <- matrix(x0[i, ], nrow = 1, byrow = TRUE)
    Y_i_0 <- y0[i]

    prediction <- predict(model_main, X_i_0, s = lambda)

    residual <- Y_i_0 - prediction
    residuals[i] <- residual
  }

  return(residuals)
}

#' Select and refit a random forest by out-of-bag error.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param residuals Response or residual vector aligned with x0.
#' @return List containing the fitted residual forest and selected hyperparameters.
train_rf_model <- function(x0, residuals) {
  x0_matrix <- as.matrix(x0)
  residuals_matrix <- as.vector(residuals)

  colnames(x0_matrix) <- paste0("V", 1:ncol(x0_matrix))

  ntree_range <- seq(200)
  mtry_range <- seq(2)
  nodesize_range <- c(10)

  param_grid <- expand.grid(ntree = ntree_range,
                            mtry = mtry_range,
                            nodesize = nodesize_range)

  error_rates <- numeric(nrow(param_grid))

  for (i in seq_len(nrow(param_grid))) {
    rf_model <- randomForest(x = x0_matrix,
                             y = residuals_matrix,
                             ntree = param_grid$ntree[i],
                             mtry = param_grid$mtry[i],
                             nodesize = param_grid$nodesize[i])

    error_rates[i] <- rf_model$mse[which.min(rf_model$mse)]
  }

  best_index <- which.min(error_rates)
  best_params <- param_grid[best_index, ]

  best_rf_model <- randomForest(x = x0_matrix,
                                y = residuals_matrix,
                                ntree = best_params$ntree,
                                mtry = best_params$mtry,
                                nodesize = best_params$nodesize)

  cat("Best hyperparameters:\n")
  cat("ntree:", best_params$ntree, "\n")
  cat("mtry:", best_params$mtry, "\n")
  cat("nodesize:", best_params$nodesize, "\n")

  return(list(model_residual_rf = best_rf_model, best_params = best_params))
}

#' Fit a gradient-boosted residual model.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param residuals Response or residual vector aligned with x0.
#' @return List containing the fitted residual GBM and its tuning information.
train_gbm_model <- function(x0, residuals) {
  x0_matrix <- as.matrix(x0)
  data_for_gbm <- as.data.frame(cbind(residuals, x0_matrix))

  colnames(data_for_gbm) <- c("residuals", paste0("V", 1:ncol(x0_matrix)))

  variable_names <- setdiff(colnames(data_for_gbm), "residuals")
  formula_gbm <- as.formula(paste("residuals ~", paste(variable_names, collapse = " + ")))

  model_gbm <- gbm(
    formula_gbm,
    data = data_for_gbm,
    distribution = "gaussian",
    n.trees = 500,
    interaction.depth = 5,
    shrinkage = 0.01,
    cv.folds = 5,
    verbose = FALSE
  )

  return(list(model_residual_gbm = model_gbm, best_parameters = model_gbm$bestTune))
}

#' Select exposure predictions by mean squared error.
#'
#' @param initial_predictions Candidate exposure predictions aligned with actuals.
#' @param combined_predictions Candidate exposure predictions aligned with actuals.
#' @param combined_predictions2 Candidate exposure predictions aligned with actuals.
#' @param combined_predictions3 Candidate exposure predictions aligned with actuals.
#' @param combined_predictions4 Candidate exposure predictions aligned with actuals.
#' @param combined_predictions5 Candidate exposure predictions aligned with actuals.
#' @param target_predictions Candidate exposure predictions aligned with actuals.
#' @param target_predictions2 Candidate exposure predictions aligned with actuals.
#' @param actuals Observed exposure values aligned with the predictions.
#' @return Integer index of the candidate with the smallest mean squared error.
select_optimal_model <- function(initial_predictions, combined_predictions,combined_predictions2,combined_predictions3,combined_predictions4,combined_predictions5, target_predictions,target_predictions2, actuals) {
  mse_initial <- mean((initial_predictions - actuals) ^ 2)
  mse_combined <- mean((combined_predictions - actuals) ^ 2)
  mse_combined2 <- mean((combined_predictions2 - actuals) ^ 2)
  mse_combined3 <- mean((combined_predictions3 - actuals) ^ 2)
  mse_combined4 <- mean((combined_predictions4 - actuals) ^ 2)
  mse_combined5 <- mean((combined_predictions5 - actuals) ^ 2)
  mse_target <- mean((target_predictions - actuals) ^ 2)
  mse_target2 <- mean((target_predictions2 - actuals) ^ 2)

  mse_values <- c(mse_initial, mse_combined, mse_combined2, mse_combined3,mse_combined4,mse_combined5, mse_target,mse_target2)
  optimal_model_index <- which.min(mse_values)

  return(optimal_model_index)
}

#' Predict exposure with source, transfer, and target models.
#'
#' @param model_main First source fit list returned by train_lasso_model.
#' @param model_main2 Second source fit list returned by train_lasso_model.
#' @param model_main3 Third source fit list returned by train_lasso_model.
#' @param model_transfer Fitted residual model for the first source.
#' @param model_transfer2 Fitted residual model for the second source.
#' @param model_transfer3 Fitted residual model for the third source.
#' @param model_target Fitted target-only model.
#' @param test_data Predictor matrix for held-out observations.
#' @param actuals Observed exposure values aligned with the predictions.
#' @return Named list of prediction vectors for each fitted model.
evaluate_and_get_predictions <- function(model_main,model_main2,model_main3, model_transfer,model_transfer2,model_transfer3,model_target,test_data, actuals) {
  colnames(test_data) <- paste0("V", 1:ncol(test_data))

  initial_predictions <- predict(model_main$model_residual, test_data, s = model_main$lambda_residual)
  initial_predictions_three <- predict(model_main2$model_residual, test_data, s = model_main2$lambda_residual)
  initial_predictions_third <- predict(model_main3$model_residual, test_data, s = model_main3$lambda_residual)

  residuals_predictions <- predict(model_transfer, newx = test_data, s = model_transfer$lambda)
  residuals_predictions_three <- predict(model_transfer2, newx = test_data, s = model_transfer2$lambda)
  residuals_predictions_third <- predict(model_transfer3, newx = test_data, s = model_transfer3$lambda)

  combined_predictions <- initial_predictions + residuals_predictions
  combined_predictions_three <- initial_predictions_three + residuals_predictions_three
  combined_predictions_third <-initial_predictions_third + residuals_predictions_third

  target_predictions <- predict(model_target$glmnet.fit, newx = as.matrix(test_data), s = model_target$lambda.min[1])
  target_predictions_three <- predict(model_target$glmnet.fit, newx = as.matrix(test_data), s = model_target$lambda.min[1]*3)
  target_predictions_third <- predict(model_target$glmnet.fit, newx = as.matrix(test_data), s = model_target$lambda.min[1]/3)

  predictions_list <- list(
    initial_predictions = initial_predictions,
    combined_predictions = combined_predictions,
    target_predictions = target_predictions,
    initial_predictions_three = initial_predictions_three,
    combined_predictions_three = combined_predictions_three,
    target_predictions_three = target_predictions_three,
    initial_predictions_third = initial_predictions_third,
    combined_predictions_third = combined_predictions_third,
    target_predictions_third = target_predictions_third
  )

  return(predictions_list)
}

#' Summarize exposure prediction errors.
#'
#' @param predictions_list Named list of exposure predictions in the evaluation order.
#' @param actuals Observed exposure values aligned with the predictions.
#' @return List with mse_values, mae_values, and the mean-only baseline mse_mean.
calculate_mse_mae <- function(predictions_list, actuals) {
  mse_values <- sapply(predictions_list, function(predictions) mean((predictions - actuals) ^ 2))
  mae_values <- sapply(predictions_list, function(predictions) mean(abs(predictions - actuals)))

  names(mse_values) <- c("initial", "combined", "target","initial_three", "combined_three", "target_three","initial_third", "combined_third", "target_third")
  names(mae_values) <- c("initial", "combined","target","initial_three", "combined_three", "target_three","initial_third", "combined_third", "target_third")

  D_mean <- mean(actuals)
  mse_mean <- mean((actuals - D_mean)^2)

  list(
    mse_values = mse_values,
    mae_values = mae_values,
    mse_mean = mse_mean
  )
}

#' Evaluate outcome regressions using predicted exposures.
#'
#' @param predictions Named list of source, transfer, and target exposure predictions.
#' @param actuals_D Observed exposure vector.
#' @param actuals_Y Observed outcome vector.
#' @param true_coefficient True causal coefficient for squared-error calculations.
#' @return List with beta_coefficients and squared errors in mse_values.
calculate_regression_coefficients_and_mse <- function(predictions, actuals_D, actuals_Y, true_coefficient = 2) {

  weighted_initial_predictions <- predictions$initial_predictions
  weighted_combined_predictions1 <- predictions$combined_predictions
  weighted_target_predictions3 <- predictions$target_predictions
  weighted_initial_predictions_three <- predictions$initial_predictions_three
  weighted_combined_predictions1_three <- predictions$combined_predictions_three
  weighted_target_predictions3_three <- predictions$target_predictions_three
  weighted_initial_predictions_third <- predictions$initial_predictions_third
  weighted_combined_predictions1_third <- predictions$combined_predictions_third
  weighted_target_predictions3_third <- predictions$target_predictions_third

  lm_Y_ols <- lm(actuals_Y ~ actuals_D)
  lm_Y_initial <- lm(actuals_Y ~ weighted_initial_predictions)
  lm_Y_combined1 <- lm(actuals_Y ~ weighted_combined_predictions1)
  lm_Y_target3 <- lm(actuals_Y ~ weighted_target_predictions3)
  lm_Y_initial_three <- lm(actuals_Y ~ weighted_initial_predictions_three)
  lm_Y_combined1_three <- lm(actuals_Y ~ weighted_combined_predictions1_three)
  lm_Y_target3_three <- lm(actuals_Y ~ weighted_target_predictions3_three)
  lm_Y_initial_third <- lm(actuals_Y ~ weighted_initial_predictions_third)
  lm_Y_combined1_third <- lm(actuals_Y ~ weighted_combined_predictions1_third)
  lm_Y_target3_third <- lm(actuals_Y ~ weighted_target_predictions3_third)

  beta_values <-c(
  beta_ols <- coef(lm_Y_ols)[2],
  beta_initial <- coef(lm_Y_initial)[2],
  beta_combined1 <- coef(lm_Y_combined1)[2],
  beta_target3 <- coef(lm_Y_target3)[2],
  beta_initial_three <- coef(lm_Y_initial_three)[2],
  beta_combined1_three <- coef(lm_Y_combined1_three)[2],
  beta_target3_three <- coef(lm_Y_target3_three)[2],
  beta_initial_third <- coef(lm_Y_initial_third)[2],
  beta_combined1_third <- coef(lm_Y_combined1_third)[2],
  beta_target3_third <- coef(lm_Y_target3_third)[2]
  )

  mse_values <- c(
    mse_ols = (beta_ols - true_coefficient)^2,
    mse_initial = (beta_initial - true_coefficient)^2,
    mse_combined1 = (beta_combined1 - true_coefficient)^2,
    mse_target3 = (beta_target3 - true_coefficient)^2,
    mse_initial_three = (beta_initial_three - true_coefficient)^2,
    mse_combined1_three = (beta_combined1_three - true_coefficient)^2,
    mse_target3_three = (beta_target3_three - true_coefficient)^2,
    mse_initial_third = (beta_initial_third - true_coefficient)^2,
    mse_combined1_third = (beta_combined1_third - true_coefficient)^2,
    mse_target3_third = (beta_target3_third - true_coefficient)^2
  )

  return(list(
    beta_coefficients =  beta_values, mse_values = mse_values
  ))
}

#' Train and evaluate all lasso models for one linear replicate.
#'
#' Source and transfer lasso helpers come from model_functions.R. The three folds and output
#' ordering match the original notebooks.
#' The historical coefficient-error diagnostics use reference effect 2; beta1 controls
#' data generation. This convention is retained independently of the downstream IV summaries.
#'
#' @param i Replicate index, also used as the random seed.
#' @param selected_data_eu Source genotype matrix with 2,000 matched variant columns.
#' @param selected_data_af Target genotype matrix with the same variant columns.
#' @param random_vector Length-100 target genotype coefficient vector.
#' @param near_weight Mixture weight passed to simulate_linear_data.
#' @param results_template Empty evaluation data frame initialized by the notebook.
#' @param n Number of target observations per held-out fold; total target size is 3 * n.
#' @param n_training Number of observations in each source training sample.
#' @param beta1 True exposure effect in the outcome model.
#' @return Two-element list: stacked target data/predictions, then replicate evaluation metrics.
run_linear_simulation <- function(i, selected_data_eu, selected_data_af,
                                  random_vector, near_weight, results_template,
                                  n = 350, n_training = 20000, beta1 = 0.5) {

  results <- simulate_linear_data(n = n, n_training = n_training, beta1 = beta1, seed = i,
                                  random_vector = random_vector, selected_data_eu = selected_data_eu,
                                  selected_data_af = selected_data_af, near_weight = near_weight)

  combined_training_data <- cbind(results$Z_training1, results$X_training1)
  model_main <- train_lasso_model(combined_training_data, as.matrix(results$D_training1))

  combined_training_data <- cbind(results$Z_training2, results$X_training2)
  model_main2 <- train_lasso_model(combined_training_data, as.matrix(results$D_training2))

  combined_training_data <- cbind(results$Z_training3, results$X_training3)
  model_main3 <- train_lasso_model(combined_training_data, as.matrix(results$D_training3))

  # Fit nuisance models on the training part of fold 1.
  x0 <- cbind(results$Z_test_1,results$X_test_1)
  y0 <- results$D_test_1

  residuals <-compute_residuals(model_main$model_residual,x0, y0, lambda=model_main$lambda_residual)
  residuals2 <-compute_residuals(model_main2$model_residual,x0, y0, lambda=model_main2$lambda_residual)
  residuals3 <-compute_residuals(model_main3$model_residual,x0, y0, lambda=model_main3$lambda_residual)

  lasso_model_info <- train_lasso_model(x0, residuals)
  model_transfer <- lasso_model_info$model_residual

  lasso_model_info <- train_lasso_model(x0, residuals2)
  model_transfer2 <- lasso_model_info$model_residual

  lasso_model_info <- train_lasso_model(x0, residuals3)
  model_transfer3 <- lasso_model_info$model_residual

  combined_testing_data <- cbind(results$Z_test_1, results$X_test_1)

  model_target <- train_target_lasso_model(combined_testing_data, as.matrix(results$D_test_1))

  model_target <- model_target$model_residual

  final_testing_data <- cbind(results$Z1,results$X1)
  actuals1 <- results$D1
  predictions_result <- evaluate_and_get_predictions(model_main,model_main2,model_main3, model_transfer,model_transfer2,model_transfer3,model_target,final_testing_data,actuals1)

  initial_predictions <- predictions_result$initial_predictions
  combined_predictions <- predictions_result$combined_predictions
  target_predictions <- predictions_result$target_predictions
  initial_predictions_three <- predictions_result$initial_predictions_three
  combined_predictions_three <- predictions_result$combined_predictions_three
  target_predictions_three <- predictions_result$target_predictions_three
  initial_predictions_third <- predictions_result$initial_predictions_third
  combined_predictions_third <- predictions_result$combined_predictions_third
  target_predictions_third <- predictions_result$target_predictions_third

  predictions_list1 <- list(
    initial_predictions = initial_predictions,
    combined_predictions = combined_predictions,
    target_predictions = target_predictions,
    initial_predictions_three = initial_predictions_three,
    combined_predictions_three = combined_predictions_three,
    target_predictions_three = target_predictions_three,
    initial_predictions_third = initial_predictions_third,
    combined_predictions_third = combined_predictions_third,
    target_predictions_third = target_predictions_third
  )

  # Fit nuisance models on the training part of fold 2.
  x0 <- cbind(results$Z_test_2,results$X_test_2)
  y0 <- results$D_test_2
  residuals <-compute_residuals(model_main$model_residual,x0, y0, lambda=model_main$lambda_residual)
  residuals2 <-compute_residuals(model_main2$model_residual,x0, y0, lambda=model_main2$lambda_residual)
  residuals3 <-compute_residuals(model_main3$model_residual,x0, y0, lambda=model_main3$lambda_residual)

  lasso_model_info <- train_lasso_model(x0, residuals)
  model_transfer <- lasso_model_info$model_residual

  lasso_model_info <- train_lasso_model(x0, residuals2)
  model_transfer2 <- lasso_model_info$model_residual

  lasso_model_info <- train_lasso_model(x0, residuals3)
  model_transfer3 <- lasso_model_info$model_residual

  combined_testing_data <- cbind(results$Z_test_2, results$X_test_2)

  model_target <- train_target_lasso_model(combined_testing_data, as.matrix(results$D_test_2))
  model_target <- model_target$model_residual

  final_testing_data <- cbind(results$Z2,results$X2)
  actuals2 <- results$D2
  predictions_result <- evaluate_and_get_predictions(model_main,model_main2,model_main3, model_transfer,model_transfer2,model_transfer3,model_target,final_testing_data,actuals2)

  initial_predictions <- predictions_result$initial_predictions
  combined_predictions <- predictions_result$combined_predictions
  target_predictions <- predictions_result$target_predictions
  initial_predictions_three <- predictions_result$initial_predictions_three
  combined_predictions_three <- predictions_result$combined_predictions_three
  target_predictions_three <- predictions_result$target_predictions_three
  initial_predictions_third <- predictions_result$initial_predictions_third
  combined_predictions_third <- predictions_result$combined_predictions_third
  target_predictions_third <- predictions_result$target_predictions_third

  predictions_list2 <- list(
    initial_predictions = initial_predictions,
    combined_predictions = combined_predictions,
    target_predictions = target_predictions,
    initial_predictions_three = initial_predictions_three,
    combined_predictions_three = combined_predictions_three,
    target_predictions_three = target_predictions_three,
    initial_predictions_third = initial_predictions_third,
    combined_predictions_third = combined_predictions_third,
    target_predictions_third = target_predictions_third
  )

    # Fit nuisance models on the training part of fold 3.
  x0 <- cbind(results$Z_test_3,results$X_test_3)
  y0 <- results$D_test_3
  residuals <-compute_residuals(model_main$model_residual,x0, y0, lambda=model_main$lambda_residual)
  residuals2 <-compute_residuals(model_main2$model_residual,x0, y0, lambda=model_main2$lambda_residual)
  residuals3 <-compute_residuals(model_main3$model_residual,x0, y0, lambda=model_main3$lambda_residual)

  lasso_model_info <- train_lasso_model(x0, residuals)
  model_transfer <- lasso_model_info$model_residual

  lasso_model_info <- train_lasso_model(x0, residuals2)
  model_transfer2 <- lasso_model_info$model_residual

  lasso_model_info <- train_lasso_model(x0, residuals3)
  model_transfer3 <- lasso_model_info$model_residual

  combined_testing_data <- cbind(results$Z_test_3, results$X_test_3)

  model_target <- train_target_lasso_model(combined_testing_data, as.matrix(results$D_test_3))

  model_target <- model_target$model_residual

  final_testing_data <- cbind(results$Z3,results$X3)
  actuals3 <- results$D3
  predictions_result <- evaluate_and_get_predictions(model_main,model_main2,model_main3, model_transfer,model_transfer2,model_transfer3,model_target,final_testing_data,actuals3)

  initial_predictions <- predictions_result$initial_predictions
  combined_predictions <- predictions_result$combined_predictions
  target_predictions <- predictions_result$target_predictions
  initial_predictions_three <- predictions_result$initial_predictions_three
  combined_predictions_three <- predictions_result$combined_predictions_three
  target_predictions_three <- predictions_result$target_predictions_three
  initial_predictions_third <- predictions_result$initial_predictions_third
  combined_predictions_third <- predictions_result$combined_predictions_third
  target_predictions_third <- predictions_result$target_predictions_third

  predictions_list3 <- list(
    initial_predictions = initial_predictions,
    combined_predictions = combined_predictions,
    target_predictions = target_predictions,
    initial_predictions_three = initial_predictions_three,
    combined_predictions_three = combined_predictions_three,
    target_predictions_three = target_predictions_three,
    initial_predictions_third = initial_predictions_third,
    combined_predictions_third = combined_predictions_third,
    target_predictions_third = target_predictions_third
  )

  result1 <- calculate_mse_mae(predictions_list1, results$D1)
  result2 <- calculate_mse_mae(predictions_list2, results$D2)
  result3 <- calculate_mse_mae(predictions_list3, results$D3)

  mse_main_model1 <- result1$mse_values["initial"]
  mse_trans_combined11 <- result1$mse_values["combined"]
  mse_target_model1 <- result1$mse_values["target"]
  mse_main_model1_three <- result1$mse_values["initial_three"]
  mse_trans_combined11_three <- result1$mse_values["combined_three"]
  mse_target_model1_three <- result1$mse_values["target_three"]
  mse_main_model1_third <- result1$mse_values["initial_third"]
  mse_trans_combined11_third <- result1$mse_values["combined_third"]
  mse_target_model1_third <- result1$mse_values["target_third"]
  mse_mean1 <- result1$mse_mean

  mse_main_model2 <- result2$mse_values["initial"]
  mse_trans_combined21 <- result2$mse_values["combined"]
  mse_target_model2 <- result2$mse_values["target"]
  mse_main_model2_three <- result2$mse_values["initial_three"]
  mse_trans_combined21_three <- result2$mse_values["combined_three"]
  mse_target_model2_three <- result2$mse_values["target_three"]
  mse_main_model2_third <- result2$mse_values["initial_third"]
  mse_trans_combined21_third <- result2$mse_values["combined_third"]
  mse_target_model2_third <- result2$mse_values["target_third"]
  mse_mean2 <- result2$mse_mean

  mse_main_model3 <- result3$mse_values["initial"]
  mse_trans_combined31 <- result3$mse_values["combined"]
  mse_target_model3 <- result3$mse_values["target"]
  mse_main_model3_three <- result3$mse_values["initial_three"]
  mse_trans_combined31_three <- result3$mse_values["combined_three"]
  mse_target_model3_three <- result3$mse_values["target_three"]
  mse_main_model3_third <- result3$mse_values["initial_third"]
  mse_trans_combined31_third <- result3$mse_values["combined_third"]
  mse_target_model3_third<- result3$mse_values["target_third"]
  mse_mean3 <- result3$mse_mean

  mse_main_model <- mean(c(mse_main_model1, mse_main_model2, mse_main_model3))
  mse_trans_combined1 <- mean(c(mse_trans_combined11, mse_trans_combined21, mse_trans_combined31))
  mse_target_model <- mean(c(mse_target_model1, mse_target_model2, mse_target_model3))
  mse_main_model_three <- mean(c(mse_main_model1_three, mse_main_model2_three, mse_main_model3_three))
  mse_trans_combined1_three <- mean(c(mse_trans_combined11_three, mse_trans_combined21_three, mse_trans_combined31_three))
  mse_target_model_three <- mean(c(mse_target_model1_three, mse_target_model2_three, mse_target_model3_three))
  mse_main_model_third <- mean(c(mse_main_model1_third, mse_main_model2_third, mse_main_model3_third))
  mse_trans_combined1_third <- mean(c(mse_trans_combined11_third, mse_trans_combined21_third, mse_trans_combined31_third))
  mse_target_model_third <- mean(c(mse_target_model1_third, mse_target_model2_third, mse_target_model3_third))
  mse_mean <- mean(c(mse_mean1, mse_mean2, mse_mean3))

  results_final1 <- calculate_regression_coefficients_and_mse(predictions_list1, results$D1, results$Y1, true_coefficient = 2)
  results_final2 <- calculate_regression_coefficients_and_mse(predictions_list2, results$D2, results$Y2, true_coefficient = 2)
  results_final3 <- calculate_regression_coefficients_and_mse(predictions_list3, results$D3, results$Y3, true_coefficient = 2)

  beta_ols1<-results_final1$beta_coefficients["actuals_D"]
  beta_initial_predictions1 <- results_final1$beta_coefficients["weighted_initial_predictions"]
  beta_combined_predictions11 <- results_final1$beta_coefficients["weighted_combined_predictions1"]
  beta_target_predictions13 <- results_final1$beta_coefficients["weighted_target_predictions3"]
  beta_initial_predictions1_three <- results_final1$beta_coefficients["weighted_initial_predictions_three"]
  beta_combined_predictions11_three <- results_final1$beta_coefficients["weighted_combined_predictions1_three"]
  beta_target_predictions13_three <- results_final1$beta_coefficients["weighted_target_predictions3_three"]
  beta_initial_predictions1_third <- results_final1$beta_coefficients["weighted_initial_predictions_third"]
  beta_combined_predictions11_third <- results_final1$beta_coefficients["weighted_combined_predictions1_third"]
  beta_target_predictions13_third <- results_final1$beta_coefficients["weighted_target_predictions3_third"]

  mse_ols1<-  results_final1$mse_values["mse_ols.actuals_D"]
  mse_initial1 <- results_final1$mse_values["mse_initial.weighted_initial_predictions"]
  mse_combined11 <- results_final1$mse_values["mse_combined1.weighted_combined_predictions1"]
  mse_target13 <- results_final1$mse_values["mse_target3.weighted_target_predictions3"]
  mse_initial1_three <- results_final1$mse_values["mse_initial_three.weighted_initial_predictions_three"]
  mse_combined11_three <- results_final1$mse_values["mse_combined1_three.weighted_combined_predictions1_three"]
  mse_target13_three <- results_final1$mse_values["mse_target3_three.weighted_target_predictions3_three"]
  mse_initial1_third <- results_final1$mse_values["mse_initial_third.weighted_initial_predictions_third"]
  mse_combined11_third <- results_final1$mse_values["mse_combined1_third.weighted_combined_predictions1_third"]
  mse_target13_third <- results_final1$mse_values["mse_target3_third.weighted_target_predictions3_third"]

  beta_ols2<-results_final2$beta_coefficients["actuals_D"]
  beta_initial_predictions2 <- results_final2$beta_coefficients["weighted_initial_predictions"]
  beta_combined_predictions21 <- results_final2$beta_coefficients["weighted_combined_predictions1"]
  beta_target_predictions23 <- results_final2$beta_coefficients["weighted_target_predictions3"]
  beta_initial_predictions2_three <- results_final2$beta_coefficients["weighted_initial_predictions_three"]
  beta_combined_predictions21_three <- results_final2$beta_coefficients["weighted_combined_predictions1_three"]
  beta_target_predictions23_three <- results_final2$beta_coefficients["weighted_target_predictions3_three"]
  beta_initial_predictions2_third <- results_final2$beta_coefficients["weighted_initial_predictions_third"]
  beta_combined_predictions21_third <- results_final2$beta_coefficients["weighted_combined_predictions1_third"]
  beta_target_predictions23_third <- results_final2$beta_coefficients["weighted_target_predictions3_third"]

  mse_ols2<-  results_final2$mse_values["mse_ols.actuals_D"]
  mse_initial2 <- results_final2$mse_values["mse_initial.weighted_initial_predictions"]
  mse_combined21 <- results_final2$mse_values["mse_combined1.weighted_combined_predictions1"]
  mse_target23 <- results_final2$mse_values["mse_target3.weighted_target_predictions3"]
  mse_initial2_three <- results_final2$mse_values["mse_initial_three.weighted_initial_predictions_three"]
  mse_combined21_three <- results_final2$mse_values["mse_combined1_three.weighted_combined_predictions1_three"]
  mse_target23_three <- results_final2$mse_values["mse_target3_three.weighted_target_predictions3_three"]
  mse_initial2_third <- results_final2$mse_values["mse_initial_third.weighted_initial_predictions_third"]
  mse_combined21_third <- results_final2$mse_values["mse_combined1_third.weighted_combined_predictions1_third"]
  mse_target23_third <- results_final2$mse_values["mse_target3_third.weighted_target_predictions3_third"]

  beta_ols3<-results_final3$beta_coefficients["actuals_D"]
  beta_initial_predictions3 <- results_final3$beta_coefficients["weighted_initial_predictions"]
  beta_combined_predictions31 <- results_final3$beta_coefficients["weighted_combined_predictions1"]
  beta_target_predictions33 <- results_final3$beta_coefficients["weighted_target_predictions3"]
  beta_initial_predictions3_three <- results_final3$beta_coefficients["weighted_initial_predictions_three"]
  beta_combined_predictions31_three <- results_final3$beta_coefficients["weighted_combined_predictions1_three"]
  beta_target_predictions33_three <- results_final3$beta_coefficients["weighted_target_predictions3_three"]
  beta_initial_predictions3_third <- results_final3$beta_coefficients["weighted_initial_predictions_third"]
  beta_combined_predictions31_third <- results_final3$beta_coefficients["weighted_combined_predictions1_third"]
  beta_target_predictions33_third <- results_final3$beta_coefficients["weighted_target_predictions3_third"]

  mse_ols3<-  results_final3$mse_values["mse_ols.actuals_D"]
  mse_initial3 <- results_final3$mse_values["mse_initial.weighted_initial_predictions"]
  mse_combined31 <- results_final3$mse_values["mse_combined1.weighted_combined_predictions1"]
  mse_target33 <- results_final3$mse_values["mse_target3.weighted_target_predictions3"]
  mse_initial3_three <- results_final3$mse_values["mse_initial_three.weighted_initial_predictions_three"]
  mse_combined31_three <- results_final3$mse_values["mse_combined1_three.weighted_combined_predictions1_three"]
  mse_target33_three <- results_final3$mse_values["mse_target3_three.weighted_target_predictions3_three"]
  mse_initial3_third <- results_final3$mse_values["mse_initial_third.weighted_initial_predictions_third"]
  mse_combined31_third <- results_final3$mse_values["mse_combined1_third.weighted_combined_predictions1_third"]
  mse_target33_third <- results_final3$mse_values["mse_target3_third.weighted_target_predictions3_third"]

  beta_ols_mean <-mean(c(beta_ols1, beta_ols2, beta_ols3))
  beta_initial_predictions <- mean(c(beta_initial_predictions1, beta_initial_predictions2, beta_initial_predictions3))
  beta_combined_predictions1 <- mean(c(beta_combined_predictions11, beta_combined_predictions21, beta_combined_predictions31))
  beta_target_predictions3 <- mean(c(beta_target_predictions13, beta_target_predictions23, beta_target_predictions33))
  beta_initial_predictions_three <- mean(c(beta_initial_predictions1_three, beta_initial_predictions2_three, beta_initial_predictions3_three))
  beta_combined_predictions1_three <- mean(c(beta_combined_predictions11_three, beta_combined_predictions21_three, beta_combined_predictions31_three))
  beta_target_predictions3_three <- mean(c(beta_target_predictions13_three, beta_target_predictions23_three, beta_target_predictions33_three))
  beta_initial_predictions_third <- mean(c(beta_initial_predictions1_third, beta_initial_predictions2_third, beta_initial_predictions3_third))
  beta_combined_predictions1_third <- mean(c(beta_combined_predictions11_third, beta_combined_predictions21_third, beta_combined_predictions31_third))
  beta_target_predictions3_third <- mean(c(beta_target_predictions13_third, beta_target_predictions23_third, beta_target_predictions33_third))

  mse_ols <- mean(c(mse_ols1, mse_ols2, mse_ols3))
  mse_initial <- mean(c(mse_initial1, mse_initial2, mse_initial3))
  mse_combined1 <- mean(c(mse_combined11, mse_combined21, mse_combined31))
  mse_target3 <- mean(c(mse_target13, mse_target23, mse_target33))
  mse_initial_three <- mean(c(mse_initial1_three, mse_initial2_three, mse_initial3_three))
  mse_combined1_three <- mean(c(mse_combined11_three, mse_combined21_three, mse_combined31_three))
  mse_target3_three <- mean(c(mse_target13_three, mse_target23_three, mse_target33_three))
  mse_initial_third <- mean(c(mse_initial1_third, mse_initial2_third, mse_initial3_third))
  mse_combined1_third <- mean(c(mse_combined11_third, mse_combined21_third, mse_combined31_third))
  mse_target3_third <- mean(c(mse_target13_third, mse_target23_third, mse_target33_third))

  results_template <- rbind(results_template, data.frame(

    mse_main_model = mse_main_model,
    mse_trans_combined1 = mse_trans_combined1,
    target_predictions = mse_target_model,
    mse_main_model_three = mse_main_model_three,
    mse_trans_combined1_three = mse_trans_combined1_three,
    target_predictions_three = mse_target_model_three,
    mse_main_model_third = mse_main_model_third,
    mse_trans_combined1_third = mse_trans_combined1_third,
    target_predictions_third = mse_target_model_third,
    mse_mean = mse_mean,
    mse_ols = mse_ols,
    mse_initial = mse_initial,
    mse_combined1 = mse_combined1,
    mse_target3 = mse_target3,
    mse_initial_three = mse_initial_three,
    mse_combined1_three = mse_combined1_three,
    mse_target3_three = mse_target3_three,
    mse_initial_third = mse_initial_third,
    mse_combined1_third = mse_combined1_third,
    mse_target3_third = mse_target3_third,
    beta_ols_mean = beta_ols_mean,
    beta_initial_predictions = beta_initial_predictions,
    beta_combined_predictions1 = beta_combined_predictions1,
    beta_target_predictions3 = beta_target_predictions3,
    beta_initial_predictions_three = beta_initial_predictions_three,
    beta_combined_predictions1_three = beta_combined_predictions1_three,
    beta_target_predictions3_three= beta_target_predictions3_three,
    beta_initial_predictions_third = beta_initial_predictions_third,
    beta_combined_predictions1_third = beta_combined_predictions1_third,
    beta_target_predictions3_third = beta_target_predictions3_third
  ))

  iteration_results_best_near <- list(
  D = c(unlist(results$D1), unlist(results$D2), unlist(results$D3)),
  D_given_ZX = c(unlist(results$D_given_ZX1), unlist(results$D_given_ZX2), unlist(results$D_given_ZX3)),
  Y = c(unlist(results$Y1), unlist(results$Y2), unlist(results$Y3)),
  InitialPredictions = c(unlist(predictions_list1$initial_predictions), unlist(predictions_list2$initial_predictions), unlist(predictions_list3$initial_predictions)),
  CombinedPredictions1 = c(unlist(predictions_list1$combined_predictions), unlist(predictions_list2$combined_predictions), unlist(predictions_list3$combined_predictions)),
  TargetPredictions = c(unlist(predictions_list1$target_predictions), unlist(predictions_list2$target_predictions), unlist(predictions_list3$target_predictions)),
    InitialPredictions_three = c(unlist(predictions_list1$initial_predictions_three), unlist(predictions_list2$initial_predictions_three), unlist(predictions_list3$initial_predictions_three)),
  CombinedPredictions1_three = c(unlist(predictions_list1$combined_predictions_three), unlist(predictions_list2$combined_predictions_three), unlist(predictions_list3$combined_predictions_three)),
  TargetPredictions_three = c(unlist(predictions_list1$target_predictions_three), unlist(predictions_list2$target_predictions_three), unlist(predictions_list3$target_predictions_three)),
    InitialPredictions_third = c(unlist(predictions_list1$initial_predictions_third), unlist(predictions_list2$initial_predictions_third), unlist(predictions_list3$initial_predictions_third)),
  CombinedPredictions1_third = c(unlist(predictions_list1$combined_predictions_third), unlist(predictions_list2$combined_predictions_third), unlist(predictions_list3$combined_predictions_third)),
  TargetPredictions_third = c(unlist(predictions_list1$target_predictions_third), unlist(predictions_list2$target_predictions_third), unlist(predictions_list3$target_predictions_third)),
  Z1 = c(unlist(results$Z1[, 1]), unlist(results$Z2[, 1]), unlist(results$Z3[, 1])),
  Z2 = c(unlist(results$Z1[, 2]), unlist(results$Z2[, 2]), unlist(results$Z3[, 2])),
  Z3 = c(unlist(results$Z1[, 3]), unlist(results$Z2[, 3]), unlist(results$Z3[, 3])),
  Z4 = c(unlist(results$Z1[, 4]), unlist(results$Z2[, 4]), unlist(results$Z3[, 4])),
  Z5 = c(unlist(results$Z1[, 5]), unlist(results$Z2[, 5]), unlist(results$Z3[, 5])),
  Z6 = c(unlist(results$Z1[, 6]), unlist(results$Z2[, 6]), unlist(results$Z3[, 6])),
  Z7 = c(unlist(results$Z1[, 7]), unlist(results$Z2[, 7]), unlist(results$Z3[, 7])),
  Z8 = c(unlist(results$Z1[, 8]), unlist(results$Z2[, 8]), unlist(results$Z3[, 8])),
  Z9 = c(unlist(results$Z1[, 9]), unlist(results$Z2[, 9]), unlist(results$Z3[, 9])),
  Z10 = c(unlist(results$Z1[, 10]), unlist(results$Z2[, 10]), unlist(results$Z3[, 10])),
  Z11 = c(unlist(results$Z1[, 11]), unlist(results$Z2[, 11]), unlist(results$Z3[, 11])),
  Z12 = c(unlist(results$Z1[, 12]), unlist(results$Z2[, 12]), unlist(results$Z3[, 12])),
  Z13 = c(unlist(results$Z1[, 13]), unlist(results$Z2[, 13]), unlist(results$Z3[, 13])),
  Z14 = c(unlist(results$Z1[, 14]), unlist(results$Z2[, 14]), unlist(results$Z3[, 14])),
  Z15 = c(unlist(results$Z1[, 15]), unlist(results$Z2[, 15]), unlist(results$Z3[, 15])),
  Z16 = c(unlist(results$Z1[, 16]), unlist(results$Z2[, 16]), unlist(results$Z3[, 16])),
  Z17 = c(unlist(results$Z1[, 17]), unlist(results$Z2[, 17]), unlist(results$Z3[, 17])),
  Z18 = c(unlist(results$Z1[, 18]), unlist(results$Z2[, 18]), unlist(results$Z3[, 18])),
  Z19 = c(unlist(results$Z1[, 19]), unlist(results$Z2[, 19]), unlist(results$Z3[, 19])),
  Z20 = c(unlist(results$Z1[, 20]), unlist(results$Z2[, 20]), unlist(results$Z3[, 20])),
  Z21 = c(unlist(results$Z1[, 21]), unlist(results$Z2[, 21]), unlist(results$Z3[, 21])),
  Z22 = c(unlist(results$Z1[, 22]), unlist(results$Z2[, 22]), unlist(results$Z3[, 22])),
  Z23 = c(unlist(results$Z1[, 23]), unlist(results$Z2[, 23]), unlist(results$Z3[, 23])),
  Z24 = c(unlist(results$Z1[, 24]), unlist(results$Z2[, 24]), unlist(results$Z3[, 24])),
  Z25 = c(unlist(results$Z1[, 25]), unlist(results$Z2[, 25]), unlist(results$Z3[, 25])),
  Z26 = c(unlist(results$Z1[, 26]), unlist(results$Z2[, 26]), unlist(results$Z3[, 26])),
  Z27 = c(unlist(results$Z1[, 27]), unlist(results$Z2[, 27]), unlist(results$Z3[, 27])),
  Z28 = c(unlist(results$Z1[, 28]), unlist(results$Z2[, 28]), unlist(results$Z3[, 28])),
  Z29 = c(unlist(results$Z1[, 29]), unlist(results$Z2[, 29]), unlist(results$Z3[, 29])),
  Z30 = c(unlist(results$Z1[, 30]), unlist(results$Z2[, 30]), unlist(results$Z3[, 30])),
  Z31 = c(unlist(results$Z1[, 31]), unlist(results$Z2[, 31]), unlist(results$Z3[, 31])),
  Z32 = c(unlist(results$Z1[, 32]), unlist(results$Z2[, 32]), unlist(results$Z3[, 32])),
  Z33 = c(unlist(results$Z1[, 33]), unlist(results$Z2[, 33]), unlist(results$Z3[, 33])),
  Z34 = c(unlist(results$Z1[, 34]), unlist(results$Z2[, 34]), unlist(results$Z3[, 34])),
  Z35 = c(unlist(results$Z1[, 35]), unlist(results$Z2[, 35]), unlist(results$Z3[, 35])),
  Z36 = c(unlist(results$Z1[, 36]), unlist(results$Z2[, 36]), unlist(results$Z3[, 36])),
  Z37 = c(unlist(results$Z1[, 37]), unlist(results$Z2[, 37]), unlist(results$Z3[, 37])),
  Z38 = c(unlist(results$Z1[, 38]), unlist(results$Z2[, 38]), unlist(results$Z3[, 38])),
  Z39 = c(unlist(results$Z1[, 39]), unlist(results$Z2[, 39]), unlist(results$Z3[, 39])),
  Z40 = c(unlist(results$Z1[, 40]), unlist(results$Z2[, 40]), unlist(results$Z3[, 40])),
  Z41 = c(unlist(results$Z1[, 41]), unlist(results$Z2[, 41]), unlist(results$Z3[, 41])),
  Z42 = c(unlist(results$Z1[, 42]), unlist(results$Z2[, 42]), unlist(results$Z3[, 42])),
  Z43 = c(unlist(results$Z1[, 43]), unlist(results$Z2[, 43]), unlist(results$Z3[, 43])),
  Z44 = c(unlist(results$Z1[, 44]), unlist(results$Z2[, 44]), unlist(results$Z3[, 44])),
  Z45= c(unlist(results$Z1[, 45]), unlist(results$Z2[, 45]), unlist(results$Z3[, 45])),
  Z46 = c(unlist(results$Z1[, 46]), unlist(results$Z2[, 46]), unlist(results$Z3[, 46])),
  Z47 = c(unlist(results$Z1[, 47]), unlist(results$Z2[, 47]), unlist(results$Z3[, 47])),
  Z48 = c(unlist(results$Z1[, 48]), unlist(results$Z2[, 48]), unlist(results$Z3[, 48])),
  Z49 = c(unlist(results$Z1[, 49]), unlist(results$Z2[, 49]), unlist(results$Z3[, 49])),
  Z50 = c(unlist(results$Z1[, 50]), unlist(results$Z2[, 50]), unlist(results$Z3[, 50])),
  Z51 = c(unlist(results$Z1[, 51]), unlist(results$Z2[, 51]), unlist(results$Z3[, 51])),
  Z52 = c(unlist(results$Z1[, 52]), unlist(results$Z2[, 52]), unlist(results$Z3[, 52])),
  Z53 = c(unlist(results$Z1[, 53]), unlist(results$Z2[, 53]), unlist(results$Z3[, 53])),
  Z54 = c(unlist(results$Z1[, 54]), unlist(results$Z2[, 54]), unlist(results$Z3[, 54])),
  Z55 = c(unlist(results$Z1[, 55]), unlist(results$Z2[, 55]), unlist(results$Z3[, 55])),
  Z56 = c(unlist(results$Z1[, 56]), unlist(results$Z2[, 56]), unlist(results$Z3[, 56])),
  Z57 = c(unlist(results$Z1[, 57]), unlist(results$Z2[, 57]), unlist(results$Z3[, 57])),
  Z58 = c(unlist(results$Z1[, 58]), unlist(results$Z2[, 58]), unlist(results$Z3[, 58])),
  Z59 = c(unlist(results$Z1[, 59]), unlist(results$Z2[, 59]), unlist(results$Z3[, 59])),
  Z60 = c(unlist(results$Z1[, 60]), unlist(results$Z2[, 60]), unlist(results$Z3[, 60])),
  Z61 = c(unlist(results$Z1[, 61]), unlist(results$Z2[, 61]), unlist(results$Z3[, 61])),
  Z62 = c(unlist(results$Z1[, 62]), unlist(results$Z2[, 62]), unlist(results$Z3[, 62])),
  Z63 = c(unlist(results$Z1[, 63]), unlist(results$Z2[, 63]), unlist(results$Z3[, 63])),
  Z64 = c(unlist(results$Z1[, 64]), unlist(results$Z2[, 64]), unlist(results$Z3[, 64])),
  Z65 = c(unlist(results$Z1[, 65]), unlist(results$Z2[, 65]), unlist(results$Z3[, 65])),
  Z66 = c(unlist(results$Z1[, 66]), unlist(results$Z2[, 66]), unlist(results$Z3[, 66])),
  Z67 = c(unlist(results$Z1[, 67]), unlist(results$Z2[, 67]), unlist(results$Z3[, 67])),
  Z68 = c(unlist(results$Z1[, 68]), unlist(results$Z2[, 68]), unlist(results$Z3[, 68])),
  Z69 = c(unlist(results$Z1[, 69]), unlist(results$Z2[, 69]), unlist(results$Z3[, 69])),
  Z70 = c(unlist(results$Z1[, 70]), unlist(results$Z2[, 70]), unlist(results$Z3[, 70])),
  Z71 = c(unlist(results$Z1[, 71]), unlist(results$Z2[, 71]), unlist(results$Z3[, 71])),
  Z72 = c(unlist(results$Z1[, 72]), unlist(results$Z2[, 72]), unlist(results$Z3[, 72])),
  Z73 = c(unlist(results$Z1[, 73]), unlist(results$Z2[, 73]), unlist(results$Z3[, 73])),
  Z74 = c(unlist(results$Z1[, 74]), unlist(results$Z2[, 74]), unlist(results$Z3[, 74])),
  Z75 = c(unlist(results$Z1[, 75]), unlist(results$Z2[, 75]), unlist(results$Z3[, 75])),
  Z76 = c(unlist(results$Z1[, 76]), unlist(results$Z2[, 76]), unlist(results$Z3[, 76])),
  Z77 = c(unlist(results$Z1[, 77]), unlist(results$Z2[, 77]), unlist(results$Z3[, 77])),
  Z78 = c(unlist(results$Z1[, 78]), unlist(results$Z2[, 78]), unlist(results$Z3[, 78])),
  Z79 = c(unlist(results$Z1[, 79]), unlist(results$Z2[, 79]), unlist(results$Z3[, 79])),
  Z80 = c(unlist(results$Z1[, 80]), unlist(results$Z2[, 80]), unlist(results$Z3[, 80])),
  Z81 = c(unlist(results$Z1[, 81]), unlist(results$Z2[, 81]), unlist(results$Z3[, 81])),
  Z82 = c(unlist(results$Z1[, 82]), unlist(results$Z2[, 82]), unlist(results$Z3[, 82])),
  Z83 = c(unlist(results$Z1[, 83]), unlist(results$Z2[, 83]), unlist(results$Z3[, 83])),
  Z84 = c(unlist(results$Z1[, 84]), unlist(results$Z2[, 84]), unlist(results$Z3[, 84])),
  Z85 = c(unlist(results$Z1[, 85]), unlist(results$Z2[, 85]), unlist(results$Z3[, 85])),
  Z86 = c(unlist(results$Z1[, 86]), unlist(results$Z2[, 86]), unlist(results$Z3[, 86])),
  Z87 = c(unlist(results$Z1[, 87]), unlist(results$Z2[, 87]), unlist(results$Z3[, 87])),
  Z88 = c(unlist(results$Z1[, 88]), unlist(results$Z2[, 88]), unlist(results$Z3[, 88])),
  Z89 = c(unlist(results$Z1[, 89]), unlist(results$Z2[, 89]), unlist(results$Z3[, 89])),
  Z90 = c(unlist(results$Z1[, 90]), unlist(results$Z2[, 90]), unlist(results$Z3[, 90])),
  Z91 = c(unlist(results$Z1[, 91]), unlist(results$Z2[, 91]), unlist(results$Z3[, 91])),
  Z92 = c(unlist(results$Z1[, 92]), unlist(results$Z2[, 92]), unlist(results$Z3[, 92])),
  Z93 = c(unlist(results$Z1[, 93]), unlist(results$Z2[, 93]), unlist(results$Z3[, 93])),
  Z94 = c(unlist(results$Z1[, 94]), unlist(results$Z2[, 94]), unlist(results$Z3[, 94])),
  Z95= c(unlist(results$Z1[, 95]), unlist(results$Z2[, 95]), unlist(results$Z3[, 95])),
  Z96 = c(unlist(results$Z1[, 96]), unlist(results$Z2[, 96]), unlist(results$Z3[, 96])),
  Z97 = c(unlist(results$Z1[, 97]), unlist(results$Z2[, 97]), unlist(results$Z3[, 97])),
  Z98 = c(unlist(results$Z1[, 98]), unlist(results$Z2[, 98]), unlist(results$Z3[, 98])),
  Z99 = c(unlist(results$Z1[, 99]), unlist(results$Z2[, 99]), unlist(results$Z3[, 99])),
  Z100 = c(unlist(results$Z1[, 100]), unlist(results$Z2[, 100]), unlist(results$Z3[, 100])),
  X1 = c(unlist(results$X1[, 1]), unlist(results$X2[, 1]), unlist(results$X3[, 1])),
  X2 = c(unlist(results$X1[, 2]), unlist(results$X2[, 2]), unlist(results$X3[, 2])),
  X3 = c(unlist(results$X1[, 3]), unlist(results$X2[, 3]), unlist(results$X3[, 3])),
  X4 = c(unlist(results$X1[, 4]), unlist(results$X2[, 4]), unlist(results$X3[, 4])),
  X5 = c(unlist(results$X1[, 5]), unlist(results$X2[, 5]), unlist(results$X3[, 5]))
)

  Results_list_best_new_near <- iteration_results_best_near
  list(Results_list_best_new_near, results_template)
}
