# Shared functions; source this file from a workflow before calling them.
# Source model_functions.R first. Load the packages listed in Readme.txt.

#' Generate nonlinear source data and three target folds.
#'
#' @param n Number of target observations per held-out fold; total target size is 3 * n.
#' @param n_training Number of observations in each source training sample.
#' @param beta1 True exposure effect in the outcome model.
#' @return List of source data, target folds, and conditional exposure means.
simulate_and_calculate_variances <- function(n=350, n_training=20000, beta1=2) {
  nc <- 15
  n_test=n*2
  total_n <- n + n_test
  nx <- 5

  Z <- mvrnorm(total_n, mu = rep(0, nc), Sigma = diag(nc))
  Z_training <- mvrnorm(n_training, mu = rep(0, nc), Sigma = diag(nc))

  X <- mvrnorm(total_n, mu = rep(0, nx), Sigma = diag(nx))
  X_training <- mvrnorm(n_training, mu = rep(0, nx), Sigma = diag(nx))

  beta_Z <- c(0.06878480, -0.02157827, 0.21446026, 0.02387094, -0.01881773, 0.29892867, -0.07078311, -0.09823684, 0.06865451, -0.23106484, -0.59164275, 0.36938943, -0.21947054, -0.02041420, -0.02733346) *2
  beta_X<-c(0.2, 0.2, 0.2, 0.2, 0.2)

  beta_X_Y<-c(0.1, 0.1, 0.1, 0.1, 0.1)

  # Source and target exposure models differ in their functional terms.
  beta_Z_training <- beta_Z

  C <- rnorm(total_n, mean = 0, sd = 1)
  D <- log(1 + exp(Z[, 6] * beta_Z[6] + Z[, 7] * beta_Z[7] +
                  Z[, 8] * beta_Z[8] + Z[, 9] * beta_Z[9] +
                  Z[, 10] * beta_Z[10]) +
            (Z[, 11] * beta_Z[11] + Z[, 12] * beta_Z[12] +
             Z[, 13] * beta_Z[13] + Z[, 14] * beta_Z[14] +
             Z[, 15] * beta_Z[15])^2) +
       Z[, 1] * beta_Z[1] +
       Z[, 2] * beta_Z[2] +
       Z[, 3] * beta_Z[3] +
       Z[, 4] * beta_Z[4] +
       Z[, 5] * beta_Z[5] +
       X[, 1] * beta_X[1] +
       X[, 2] * beta_X[2] +
       X[, 3] * beta_X[3] +
       log(abs(X[, 4]) + 1) * beta_X[4] +
       (X[, 5] * X[, 5]) * beta_X[5] +
       2*C+
       rnorm(total_n, mean = 0, sd = 2)

  D_given_ZX <- log(1 + exp(Z[, 6] * beta_Z[6] + Z[, 7] * beta_Z[7] +
                  Z[, 8] * beta_Z[8] + Z[, 9] * beta_Z[9] +
                  Z[, 10] * beta_Z[10]) +
            (Z[, 11] * beta_Z[11] + Z[, 12] * beta_Z[12] +
             Z[, 13] * beta_Z[13] + Z[, 14] * beta_Z[14] +
             Z[, 15] * beta_Z[15])^2) +
       Z[, 1] * beta_Z[1] +
       Z[, 2] * beta_Z[2] +
       Z[, 3] * beta_Z[3] +
       Z[, 4] * beta_Z[4] +
       Z[, 5] * beta_Z[5] +
       X[, 1] * beta_X[1] +
       X[, 2] * beta_X[2] +
       X[, 3] * beta_X[3] +
       log(abs(X[, 4]) + 1) * beta_X[4] +
       (X[, 5] * X[, 5]) * beta_X[5]

  Y <- beta1 * D +
       X[, 1] * beta_X_Y[1] +
       X[, 2] * beta_X_Y[2] +
       X[, 3] * beta_X_Y[3] +
       X[, 4] * beta_X_Y[4] +
       X[, 5] * beta_X_Y[5] +
       4*C +
       rnorm(total_n, mean = 0, sd = 0.5)

  shuffled_indices <- sample(1:total_n)

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

  C_training <- rnorm(n_training, mean = 0, sd = 1)
  D_training <- log(1 + exp(Z_training[, 6] * beta_Z_training[6] +
                  Z_training[, 7] * beta_Z_training[7] +
                  Z_training[, 8] * beta_Z_training[8] +
                  Z_training[, 9] * beta_Z_training[9] +
                  Z_training[, 10] * beta_Z_training[10]) +
            (Z_training[, 11] * beta_Z_training[11] +
             Z_training[, 12] * beta_Z_training[12] +
             Z_training[, 13] * beta_Z_training[13] +
             Z_training[, 14] * beta_Z_training[14] +
             Z_training[, 15] * beta_Z_training[15])^2) +
       X_training[, 1] * beta_X[1] +
       X_training[, 2] * beta_X[2] +
       X_training[, 3] * beta_X[3] +
       log(abs(X_training[, 4]) + 1) * beta_X[4] +
       (X_training[, 5] * X_training[, 5]) * beta_X[5] +
       2*C_training +
       2*rnorm(n_training, mean = 0, sd = 2)

  D_training_givenZX <- log(1 + exp(Z_training[, 6] * beta_Z_training[6] +
                  Z_training[, 7] * beta_Z_training[7] +
                  Z_training[, 8] * beta_Z_training[8] +
                  Z_training[, 9] * beta_Z_training[9] +
                  Z_training[, 10] * beta_Z_training[10]) +
            (Z_training[, 11] * beta_Z_training[11] +
             Z_training[, 12] * beta_Z_training[12] +
             Z_training[, 13] * beta_Z_training[13] +
             Z_training[, 14] * beta_Z_training[14] +
             Z_training[, 15] * beta_Z_training[15])^2) +
       X_training[, 1] * beta_X[1] +
       X_training[, 2] * beta_X[2] +
       X_training[, 3] * beta_X[3] +
       log(abs(X_training[, 4]) + 1) * beta_X[4] +
       (X_training[, 5] * X_training[, 5]) * beta_X[5]

  Y_training <-  beta1 * D_training +
                 X_training[, 1] * beta_X_Y[1] +
                 X_training[, 2] * beta_X_Y[2] +
                 X_training[, 3] * beta_X_Y[3] +
                 X_training[, 4] * beta_X_Y[4] +
                 X_training[, 5] * beta_X_Y[5] +
                 4*C_training +
                 rnorm(n_training, mean = 0, sd = 0.5)

  return(list(
  Z1=Z1, X1=X1, D1=D1, Y1=Y1, D_given_ZX1=D_given_ZX1,
  Z_test_1=Z_test, X_test_1=X_test, D_test_1=D_test, Y_test_1=Y_test,

  Z2=Z2, X2=X2, D2=D2, Y2=Y2,D_given_ZX2=D_given_ZX2,
  Z_test_2=Z_test_2, X_test_2=X_test_2, D_test_2=D_test_2, Y_test_2=Y_test_2,

  Z3=Z3, X3=X3, D3=D3, Y3=Y3,D_given_ZX3=D_given_ZX3,
  Z_test_3=Z_test_3, X_test_3=X_test_3, D_test_3=D_test_3, Y_test_3=Y_test_3,

  Z_training=Z_training, X_training=X_training, D_training=D_training, Y_training=Y_training,
  D_training_givenZX=D_training_givenZX
))
}

#' Compute source-model residuals observation by observation.
#'
#' @param model_main Fitted source prediction model.
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @return Numeric residual vector in observation order.
compute_residuals <- function(model_main, x0, y0) {

  residuals <- numeric(length(y0))

  for (i in 1:nrow(x0)) {

    X_i_0 <- matrix(x0[i, ], nrow = 1, byrow = TRUE)
    Y_i_0 <- y0[i]

    prediction <- model_main$predict(X_i_0,verbose=0)

    residual <- Y_i_0 - prediction
    residuals[i] <- residual
  }

  return(residuals)
}

#' Select and refit a random forest by out-of-bag error.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param residuals Response or residual vector aligned with x0.
#' @param param_grid Data frame with one candidate hyperparameter combination per row.
#' @return List containing the fitted residual forest and selected hyperparameters.
train_rf_model <- function(x0, residuals, param_grid) {
  x0_matrix <- as.matrix(x0)
  residuals_matrix <- as.vector(residuals)

  colnames(x0_matrix) <- paste0("V", 1:ncol(x0_matrix))

  error_rates <- numeric(nrow(param_grid))

  for (i in seq_len(nrow(param_grid))) {
    params <- param_grid[i, ]

    rf_model <- randomForest(x0_matrix, residuals_matrix,
                             ntree = params$ntree,
                             mtry = params$mtry,
                             nodesize = params$nodesize,
                             maxnodes = ifelse(is.finite(params$maxnodes), params$maxnodes, NULL))

    error_rates[i] <- rf_model$mse[which.min(rf_model$mse)]}

  best_index <- which.min(error_rates)
  best_params <- param_grid[best_index, ]

  best_rf_model <- randomForest(x0_matrix, residuals_matrix,
                                ntree = best_params$ntree,
                                mtry = best_params$mtry,
                                nodesize = best_params$nodesize,
                                maxnodes = ifelse(is.finite(best_params$maxnodes), best_params$maxnodes, NULL))

  return(list(model_residual_rf = best_rf_model, best_residual_params = best_params))
}

#' Fit a gradient-boosted residual model.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @param param_grid Data frame with one candidate hyperparameter combination per row.
#' @return List containing the fitted residual GBM and its tuning information.
train_gbm_model <- function(x0, y0, param_grid) {
  x0_matrix <- as.matrix(x0)
  y0_vector <- y0
  data_for_gbm <- as.data.frame(cbind(y0_vector, x0_matrix))

  colnames(data_for_gbm) <- c("y0", paste0("V", 1:ncol(x0_matrix)))

  variable_names <- setdiff(colnames(data_for_gbm), "y0")
  formula_gbm <- as.formula(paste("y0 ~", paste(variable_names, collapse = " + ")))

  best_model <- NULL
  best_params <- NULL
  best_mse_train <- Inf

  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]
    model_gbm <- gbm(
      formula_gbm,
      data = data_for_gbm,
      distribution = "gaussian",
      n.trees = params$n.trees,
      interaction.depth = params$interaction.depth,
      shrinkage = params$shrinkage,
      verbose = FALSE
    )

    predicted_values <- predict(model_gbm, newdata = data_for_gbm, n.trees = params$n.trees)

    mse_train <- mean((y0_vector - predicted_values)^2)

    if (mse_train < best_mse_train) {
      best_mse_train <- mse_train
      best_model <- model_gbm
      best_params <- params
    }
  }

  variance_y0 <- var(y0_vector)

  return(list(
    model_transfer_gbm = best_model,
    best_transfer_parameters = best_params,
    mse_train = best_mse_train,
    variance_y0 = variance_y0
  ))
}

#' Tune a target-only gradient-boosted model.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @param param_grid Data frame with one candidate hyperparameter combination per row.
#' @return List containing the target GBM and its selected hyperparameters.
train_target_gbm_model <- function(x0, y0, param_grid) {
  x0_matrix <- as.matrix(x0)
  y0_vector <- y0
  data_for_gbm <- as.data.frame(cbind(y0_vector, x0_matrix))

  colnames(data_for_gbm) <- c("y0", paste0("V", 1:ncol(x0_matrix)))

  variable_names <- setdiff(colnames(data_for_gbm), "y0")
  formula_gbm <- as.formula(paste("y0 ~", paste(variable_names, collapse = " + ")))

  best_model <- NULL
  best_params <- NULL
  best_mse_train <- Inf

  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]
    model_gbm <- gbm(
      formula_gbm,
      data = data_for_gbm,
      distribution = "gaussian",
      n.trees = params$n.trees,
      interaction.depth = params$interaction.depth,
      shrinkage = params$shrinkage,
      verbose = FALSE
    )

    predicted_values <- predict(model_gbm, newdata = data_for_gbm, n.trees = params$n.trees)

    mse_train <- mean((y0_vector - predicted_values)^2)

    if (mse_train < best_mse_train) {
      best_mse_train <- mse_train
      best_model <- model_gbm
      best_params <- params
    }
  }

  variance_y0 <- var(y0_vector)

  return(list(
    model_target_gbm = best_model,
    best_target_parameters = best_params,
    mse_train = best_mse_train,
    variance_y0 = variance_y0
  ))
}

#' Predict exposure with source, transfer, and target models.
#'
#' @param best_model_source_nn Fitted model source nn object.
#' @param best_model_source_rf Fitted model source rf object.
#' @param best_model_source_lasso Fitted model source lasso object.
#' @param best_model_source_gbm Fitted model source gbm object.
#' @param param_source_gbm Selected GBM hyperparameters, including n.trees.
#' @param param_transfer_gmb Selected GBM hyperparameters, including n.trees.
#' @param param_target_gmb Selected GBM hyperparameters, including n.trees.
#' @param lambda_source_lasso Selected lasso prediction penalty for the corresponding model.
#' @param model_residual Fitted model residual object.
#' @param model_residual_rf Fitted model residual rf object.
#' @param model_residual_gbm Fitted model residual gbm object.
#' @param model_tnn1 Fitted model tnn1 object.
#' @param model_tnn2 Fitted model tnn2 object.
#' @param model_target Fitted target-only model.
#' @param rf_model Fitted rf model object.
#' @param model_target_gbm Fitted model target gbm object.
#' @param test_data Predictor matrix for held-out observations.
#' @param actuals Observed exposure values aligned with the predictions.
#' @param source_best_flag Selected source model family: nn, rf, gbm, or lasso.
#' @return Named list of prediction vectors for each fitted model.
evaluate_and_get_predictions <- function(best_model_source_nn,
                                         best_model_source_rf,
                                         best_model_source_lasso,
                                         best_model_source_gbm,
                                         param_source_gbm,
                                         param_transfer_gmb,
                                         param_target_gmb,
                                         lambda_source_lasso,
                                         model_residual, model_residual_rf,
                                         model_residual_gbm, model_tnn1, model_tnn2,
                                         model_target, rf_model, model_target_gbm,
                                         test_data, actuals, source_best_flag) {

  colnames(test_data) <- paste0("V", 1:ncol(test_data))
  test_data_matrix <- as.matrix(test_data)
  test_data_frame <- as.data.frame(test_data)

  initial_predictions_nn <- best_model_source_nn$predict(test_data_matrix, verbose = 0)
  initial_predictions_rf <- predict(best_model_source_rf, newdata = test_data_matrix)
  initial_predictions_lasso <- predict(best_model_source_lasso, newx = test_data_matrix, s = lambda_source_lasso)
  initial_predictions_gbm <- predict(best_model_source_gbm, newdata = test_data_frame, n.trees = param_source_gbm$n.trees)

  residuals_predictions <- predict(model_residual, newx = test_data_matrix, s = model_residual$lambda.min)
  residual_predictions_rf <- predict(model_residual_rf, newdata = test_data_matrix)
  residual_predictions_gbm <- predict(model_residual_gbm, newdata = test_data_frame, n.trees = param_transfer_gmb$n.trees)
  residual_predictions_tnn1 <- model_tnn1$predict(test_data_matrix)
  residual_predictions_tnn2 <- model_tnn2$predict(test_data_matrix)

  if (source_best_flag == "nn") {
      selected_initial_predictions <- initial_predictions_nn
      combined_predictions <- initial_predictions_nn + residuals_predictions
      combined_predictions2 <- initial_predictions_nn + residual_predictions_rf
      combined_predictions3 <- initial_predictions_nn + residual_predictions_gbm
      combined_predictions4 <- initial_predictions_nn + residual_predictions_tnn1
      combined_predictions5 <- initial_predictions_nn + residual_predictions_tnn2
  } else if (source_best_flag == "rf") {
      selected_initial_predictions <- initial_predictions_rf
      combined_predictions <- initial_predictions_rf + residuals_predictions
      combined_predictions2 <- initial_predictions_rf + residual_predictions_rf
      combined_predictions3 <- initial_predictions_rf + residual_predictions_gbm
      combined_predictions4 <- initial_predictions_rf + residual_predictions_tnn1
      combined_predictions5 <- initial_predictions_rf + residual_predictions_tnn2
  } else if (source_best_flag == "lasso") {
      selected_initial_predictions <- initial_predictions_lasso
      combined_predictions <- initial_predictions_lasso + residuals_predictions
      combined_predictions2 <- initial_predictions_lasso + residual_predictions_rf
      combined_predictions3 <- initial_predictions_lasso + residual_predictions_gbm
      combined_predictions4 <- initial_predictions_lasso + residual_predictions_tnn1
      combined_predictions5 <- initial_predictions_lasso + residual_predictions_tnn2
  } else if (source_best_flag == "gbm") {
      selected_initial_predictions <- initial_predictions_gbm
      combined_predictions <- initial_predictions_gbm + residuals_predictions
      combined_predictions2 <- initial_predictions_gbm + residual_predictions_rf
      combined_predictions3 <- initial_predictions_gbm + residual_predictions_gbm
      combined_predictions4 <- initial_predictions_gbm + residual_predictions_tnn1
      combined_predictions5 <- initial_predictions_gbm + residual_predictions_tnn2
  } else {
      stop("Error: Invalid source_best_flag. Choose from 'nn', 'rf', 'lasso', or 'gbm'.")
  }

  target_predictions <- model_target$predict(test_data_matrix)

  rf_predictions_vector <- predict(rf_model, newdata = test_data_matrix)
  rf_predictions_matrix <- matrix(rf_predictions_vector, ncol = 1)

  gbm_predictions_vector <- predict(model_target_gbm, newdata = test_data_frame, n.trees = param_target_gmb$n.trees)
  gbm_predictions_matrix <- matrix(gbm_predictions_vector, ncol = 1)

  predictions_list <- list(
    initial_predictions = selected_initial_predictions,
    initial_predictions_lasso=initial_predictions_lasso,
    initial_predictions_rf=initial_predictions_rf,
    initial_predictions_gbm=initial_predictions_gbm,
    initial_predictions_nn =initial_predictions_nn,
    combined_predictions = combined_predictions,
    combined_predictions2 = combined_predictions2,
    combined_predictions3 = combined_predictions3,
    combined_predictions4 = combined_predictions4,
    combined_predictions5 = combined_predictions5,
    target_predictions = target_predictions,
    rf_predictions = rf_predictions_matrix,
    gbm_predictions = gbm_predictions_matrix
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

  names(mse_values) <- c("initial","initial_nn","initial_gbm", "initial_lasso","initial_rf", "combined","combined2","combined3","combined4","combined5", "rf","gbm","target")
  names(mae_values) <- c("initial","initial_nn","initial_gbm", "initial_lasso","initial_rf", "combined","combined2","combined3","combined4","combined5", "rf","gbm","target")

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
  weighted_initial_predictions_nn <- predictions$initial_predictions_nn
  weighted_initial_predictions_rf <- predictions$initial_predictions_rf
  weighted_initial_predictions_gbm <- predictions$initial_predictions_gbm
  weighted_initial_predictions_lasso <- predictions$initial_predictions_lasso
  weighted_combined_predictions1 <- predictions$combined_predictions
  weighted_combined_predictions2 <- predictions$combined_predictions2
  weighted_combined_predictions3 <- predictions$combined_predictions3
  weighted_combined_predictions4 <- predictions$combined_predictions4
  weighted_combined_predictions5 <- predictions$combined_predictions4
  weighted_target_predictions <- predictions$rf_predictions
  weighted_target_predictions2 <- predictions$gbm_predictions
  weighted_target_predictions3 <- predictions$target_predictions

  lm_Y_ols <- lm(actuals_Y ~ actuals_D)
  lm_Y_initial <- lm(actuals_Y ~ weighted_initial_predictions)
  lm_Y_initial_nn <- lm(actuals_Y ~ weighted_initial_predictions_nn)
  lm_Y_initial_rf <- lm(actuals_Y ~ weighted_initial_predictions_rf)
  lm_Y_initial_gbm <- lm(actuals_Y ~ weighted_initial_predictions_gbm)
  lm_Y_initial_lasso <- lm(actuals_Y ~ weighted_initial_predictions_lasso)
  lm_Y_combined1 <- lm(actuals_Y ~ weighted_combined_predictions1)
  lm_Y_combined2 <- lm(actuals_Y ~ weighted_combined_predictions2)
  lm_Y_combined3 <- lm(actuals_Y ~ weighted_combined_predictions3)
  lm_Y_combined4 <- lm(actuals_Y ~ weighted_combined_predictions4)
  lm_Y_combined5 <- lm(actuals_Y ~ weighted_combined_predictions5)
  lm_Y_target <- lm(actuals_Y ~ weighted_target_predictions)
  lm_Y_target2 <- lm(actuals_Y ~ weighted_target_predictions2)
  lm_Y_target3 <- lm(actuals_Y ~ weighted_target_predictions3)

  beta_values <- c(
  beta_ols <- coef(lm_Y_ols)[2],
  beta_initial <- coef(lm_Y_initial)[2],
  beta_initial_nn <- coef(lm_Y_initial_nn)[2],
  beta_initial_rf <- coef(lm_Y_initial_rf)[2],
  beta_initial_gbm <- coef(lm_Y_initial_gbm)[2],
  beta_initial_lasso <- coef(lm_Y_initial_lasso)[2],
  beta_combined1 <- coef(lm_Y_combined1)[2],
  beta_combined2 <- coef(lm_Y_combined2)[2],
  beta_combined3 <- coef(lm_Y_combined3)[2],
  beta_combined4 <- coef(lm_Y_combined4)[2],
  beta_combined5 <- coef(lm_Y_combined5)[2],
  beta_target <- coef(lm_Y_target)[2],
  beta_target2 <- coef(lm_Y_target2)[2],
  beta_target3 <- coef(lm_Y_target3)[2]
  )

  mse_values <- c(
    ols = (beta_ols - true_coefficient)^2,
    initial = (beta_initial - true_coefficient)^2,
    initial_nn = (beta_initial_nn - true_coefficient)^2,
    initial_rf = (beta_initial_rf - true_coefficient)^2,
    initial_gbm = (beta_initial_gbm - true_coefficient)^2,
    initial_lasso = (beta_initial_lasso - true_coefficient)^2,
    combined1 = (beta_combined1 - true_coefficient)^2,
    combined2 = (beta_combined2 - true_coefficient)^2,
    combined3 = (beta_combined3 - true_coefficient)^2,
    combined4 = (beta_combined4 - true_coefficient)^2,
    combined5 = (beta_combined5 - true_coefficient)^2,
    target = (beta_target - true_coefficient)^2,
    target2 = (beta_target2 - true_coefficient)^2,
    target3 = (beta_target3 - true_coefficient)^2
  )

  return(list(
    beta_coefficients = beta_values,
    mse_values = mse_values
  ))
}

#' Tune the source gradient-boosted exposure model.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @param param_grid Data frame with one candidate hyperparameter combination per row.
#' @return Source GBM, best_parameters, validation MSE in mse_train, and variance_y0.
train_source_gbm_model <- function(x0, y0, param_grid) {
  x0_matrix <- as.matrix(x0)

  y0_vector <- y0
  data_for_gbm <- as.data.frame(cbind(y0_vector, x0_matrix))
  colnames(data_for_gbm) <- c("y0", paste0("V", 1:ncol(x0_matrix)))

  train_indices <- sample(1:nrow(data_for_gbm), size = floor(0.7 * nrow(data_for_gbm)))
  train_data <- data_for_gbm[train_indices, ]
  valid_data <- data_for_gbm[-train_indices, ]

  variable_names <- setdiff(colnames(data_for_gbm), "y0")
  formula_gbm <- as.formula(paste("y0 ~", paste(variable_names, collapse = " + ")))

  best_model <- NULL
  best_params <- NULL
  best_mse_valid <- Inf

  for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
    model_gbm <- gbm(
      formula_gbm,
      data = train_data,
      distribution = "gaussian",
      n.trees = params$n.trees,
      interaction.depth = params$interaction.depth,
      shrinkage = params$shrinkage,
      verbose = FALSE
    )

    predicted_values <- predict(model_gbm, newdata = valid_data, n.trees = params$n.trees)

    mse_valid <- mean((valid_data$y0 - predicted_values)^2)

    if (mse_valid < best_mse_valid) {
      best_mse_valid <- mse_valid
      best_model <- model_gbm
      best_params <- params
    }
  }

  variance_y0 <- var(y0_vector)

  return(list(
    model_source_gbm = best_model,
    best_parameters = best_params,
    mse_train = best_mse_valid,
    variance_y0 = variance_y0
  ))
}

#' Tune the source random-forest exposure model.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @param param_grid Data frame with one candidate hyperparameter combination per row.
#' @return Source forest, best_source_params, out-of-bag MSE in mse_train, and variance_y0.
train_source_rf_model <- function(x0, y0, param_grid) {
  x0_matrix <- as.matrix(x0)
  y0_vector <- as.vector(y0)

  colnames(x0_matrix) <- paste0("V", 1:ncol(x0_matrix))

  error_rates <- numeric(nrow(param_grid))

  for (i in seq_len(nrow(param_grid))) {
    params <- param_grid[i, ]

    rf_model <- randomForest(x0_matrix, y0_vector,
                             ntree = params$ntree,
                             mtry = params$mtry,
                             nodesize = params$nodesize,
                             maxnodes = ifelse(is.finite(params$maxnodes), params$maxnodes, NULL))
    error_rates[i] <- rf_model$mse[which.min(rf_model$mse)]
  }

  best_index <- which.min(error_rates)
  best_params <- param_grid[best_index, ]

  best_rf_model <- randomForest(x0_matrix, y0_vector,
                                ntree = best_params$ntree,
                                mtry = best_params$mtry,
                                nodesize = best_params$nodesize,
                                maxnodes = ifelse(is.finite(best_params$maxnodes), best_params$maxnodes, NULL))

  mse_train <- best_rf_model$mse[which.min(best_rf_model$mse)]

  variance_y0 <- var(y0_vector)

  return(list(
    model_source_rf = best_rf_model,
    best_source_params = best_params,
    mse_train = mse_train,
    variance_y0 = variance_y0
  ))
}

#' Fit and evaluate the source lasso exposure model.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param residuals Response or residual vector aligned with x0.
#' @return List with source fit, penalties, coefficients, training MSE, and response variance.
train_source_lasso_model <- function(x0, residuals) {
  x0_matrix <- as.matrix(x0)
  residuals_matrix <- as.matrix(residuals)

  colnames(x0_matrix) <- paste0("V", 1:ncol(x0_matrix))

  cv_lasso <- cv.glmnet(x0_matrix, residuals_matrix, alpha = 1)

  best_lasso_model <- glmnet(x0_matrix, residuals_matrix, alpha = 1, lambda = cv_lasso$lambda.min)

  lambda_1se <- cv_lasso$lambda.1se
  lambda_min <- cv_lasso$lambda.min

  coefficients <- coef(best_lasso_model, s = cv_lasso$lambda.min)
  coefficients2 <- coef(best_lasso_model, s = cv_lasso$lambda.1se)

  predicted_y_train <- predict(best_lasso_model, newx = x0_matrix, s = lambda_min)

  mse_train <- mean((residuals_matrix - predicted_y_train)^2)

  variance_y_train <- var(residuals_matrix)

  return(list(
    model_source_lasso = best_lasso_model,
    lambda_source_lasso = lambda_min,
    lambda_1se = lambda_1se,
    coefficients = coefficients,
    coefficients2 = coefficients2,
    mse_train = mse_train,
    variance_y_train = variance_y_train
  ))
}

#' Select a source neural network by validation loss.
#'
#' @param grid_search_params Data frame of neural-network hyperparameter candidates.
#' @param combined_training_data Source predictor matrix used to train the neural network.
#' @param D_training Source exposure vector aligned with combined_training_data.
#' @param input_shape Input dimensions accepted by the Keras model.
#' @param epochs Maximum number of neural-network training epochs.
#' @param validation_split Proportion reserved for Keras validation.
#' @return List with best_model and its minimum validation loss, best_val_loss.
grid_search_best_model <- function(grid_search_params, combined_training_data, D_training, input_shape = c(20), epochs = 200, validation_split = 0.2) {

  best_model <- NULL
  best_val_loss <- Inf

  for (i in 1:nrow(grid_search_params)) {
    params <- grid_search_params[i, ]

    model_info <- define_compile_model(
      input_shape = input_shape,
      units = unlist(params$units),
      dropout_rate = params$dropout_rate,
      l2_rate = params$l2_rate,
      patience = params$patience
    )

    history <- model_info$model$fit(
      x = combined_training_data,
      y = as.matrix(D_training),
      epochs = as.integer(epochs),
      validation_split = validation_split,
      callbacks = list(model_info$early_stopping),
      verbose = 0
    )

    val_loss <- min(history$history$val_loss)

    if (val_loss < best_val_loss) {
      best_val_loss <- val_loss
      best_model <- model_info$model
    }
  }

  return(list(
    best_model = best_model,
    best_val_loss = best_val_loss
  ))
}

#' Select a source model using target-sample prediction error.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @param best_model_lasso Fitted model lasso object.
#' @param lambda_lasso Selected source-lasso prediction penalty.
#' @param best_model_nn Fitted model nn object.
#' @param best_model_rf Fitted model rf object.
#' @param best_model_gbm Fitted model gbm object.
#' @param best_gbm_trees Number of trees used for source GBM predictions.
#' @return List with best_model_name and the selected model residuals.
evaluate_and_select_best_model <- function(x0, y0,
                                           best_model_lasso, lambda_lasso,
                                           best_model_nn,
                                           best_model_rf,
                                           best_model_gbm, best_gbm_trees) {
  x0_matrix <- as.matrix(x0)
  x0_df <- as.data.frame(x0)
  y0_vector <- as.matrix(y0)

  predicted_y_lasso <- predict(best_model_lasso, newx = x0_matrix, s = lambda_lasso)
  mse_lasso <- mean((y0_vector - predicted_y_lasso)^2)

  predicted_y_nn <- best_model_nn$predict(x0_matrix, verbose = 0)
  mse_nn <- mean((y0_vector - predicted_y_nn)^2)

  predicted_y_rf <- predict(best_model_rf, x0_matrix)
  mse_rf <- mean((y0_vector - predicted_y_rf)^2)

  predicted_y_gbm <- predict(best_model_gbm, newdata = x0_df, n.trees = best_gbm_trees)
  mse_gbm <- mean((y0_vector - predicted_y_gbm)^2)

  mse_values <- c(lasso = mse_lasso, nn = mse_nn, rf = mse_rf, gbm = mse_gbm)

  best_model_index <- which.min(mse_values)
  best_model_name <- names(mse_values)[best_model_index]

  best_model <- switch(best_model_name,
                       lasso = best_model_lasso,
                       nn = best_model_nn,
    rf = best_model_rf,
                       gbm = best_model_gbm)

  residuals <- numeric(length(y0_vector))

  for (i in 1:nrow(x0_matrix)) {
    X_i_0 <- matrix(x0_matrix[i, ], nrow = 1, byrow = TRUE)
    Y_i_0 <- y0_vector[i]

    if (best_model_name == "nn") {
      prediction <- best_model$predict(X_i_0, verbose = 0)
    } else if (best_model_name == "lasso") {
      prediction <- predict(best_model, newx = X_i_0, s = lambda_lasso)
    } else if (best_model_name == "gbm") {
      prediction <- predict(best_model, newdata = as.data.frame(X_i_0), n.trees = best_gbm_trees)
    }else {
      prediction <- predict(best_model, X_i_0)
    }

    residual <- Y_i_0 - prediction
    residuals[i] <- residual
  }

  return(list(
    best_model_name = best_model_name,
    residuals = residuals
  ))
}

#' Predict and evaluate the selected source and transferred models.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param y0 Observed response vector aligned with x0.
#' @param combined_testing_data Predictor matrix for fitting target-only models.
#' @param combined_testing_y Target-only training exposure vector.
#' @param final_testing_data Predictor matrix for final held-out evaluation.
#' @param actuals Observed exposure values aligned with the predictions.
#' @param best_model_source_lasso Fitted model source lasso object.
#' @param lambda_source_lasso Selected lasso prediction penalty for the corresponding model.
#' @param best_model_source_nn Fitted model source nn object.
#' @param best_model_source_rf Fitted model source rf object.
#' @param best_model_source_gbm Fitted model source gbm object.
#' @param model_result_source_gbm Source GBM fit list containing best_parameters.
#' @param i Replicate index used by the training workflow.
#' @param param_grid_rf_residual Candidate hyperparameter grid for the corresponding model.
#' @param param_grid_gbm_transfer Candidate hyperparameter grid for the corresponding model.
#' @param param_grid_gbm_target Candidate hyperparameter grid for the corresponding model.
#' @param param_grid_rf_target Candidate hyperparameter grid for the corresponding model.
#' @return List with predictions_list and the selected source family in source_best_flagreturn.
process_vs_return_predict_list <- function(x0, y0,
                                           combined_testing_data, combined_testing_y,
                                           final_testing_data,
                                           actuals,
                                           best_model_source_lasso, lambda_source_lasso,
                                           best_model_source_nn, best_model_source_rf,
                                           best_model_source_gbm, model_result_source_gbm,
                                            i,
                                           param_grid_rf_residual,
                                           param_grid_gbm_transfer,
                                           param_grid_gbm_target,
                                           param_grid_rf_target) {

  evaluate_source_best_model <- evaluate_and_select_best_model(
    x0 = x0,
    y0 = y0,
    best_model_lasso = best_model_source_lasso,
    lambda_lasso = lambda_source_lasso,
    best_model_nn = best_model_source_nn,
    best_model_rf = best_model_source_rf,
    best_model_gbm = best_model_source_gbm,
    best_gbm_trees = model_result_source_gbm$best_parameters$n.trees
  )

  residuals <- evaluate_source_best_model$residuals
  source_best_flagreturn <- evaluate_source_best_model$best_model_name

  lasso_model_info <- train_lasso_model(x0, residuals, include_1se = TRUE)
  model_residual <- lasso_model_info$model_residual

  rf_model_info <- train_rf_model(x0, residuals, param_grid_rf_residual)
  model_residual_rf <- rf_model_info$model_residual_rf

  gbm_model_info <- train_gbm_model(x0, residuals, param_grid_gbm_transfer)
  model_residual_gbm <- gbm_model_info$model_transfer_gbm
  param_transfer_gmb <- gbm_model_info$best_transfer_parameters

  tnn_model_info1 <- define_compile_target_model(input_shape = c(ncol(x0)), units_list = c(64))
  model_tnn1 <- tnn_model_info1$model

  model_tnn1$fit(
    x = as.matrix(x0),
    y = as.matrix(residuals),
    epochs =as.integer(100) ,
    validation_split = 0.2,
    callbacks = list(tnn_model_info1$early_stopping),
    verbose = 0
  )

  tnn_model_info2 <- define_compile_target_model(input_shape = c(ncol(x0)), units_list = c(16))
  model_tnn2 <- tnn_model_info2$model

  model_tnn2$fit(
    x = as.matrix(x0),
    y = as.matrix(residuals),
    epochs = as.integer(100) ,
    validation_split = 0.2,
    callbacks = list(tnn_model_info2$early_stopping),
    verbose = 0
  )

  model_target_info <- define_compile_target_model(input_shape = c(ncol(combined_testing_data)))
  model_target <- model_target_info$model

  model_target$fit(
    x = as.matrix(combined_testing_data),
    y = as.matrix(combined_testing_y),
    epochs =as.integer(100),
    validation_split = 0.2,
    callbacks = list(model_target_info$early_stopping),
    verbose = 0
  )

  gbm_target_model_info <- train_target_gbm_model(combined_testing_data, as.matrix(combined_testing_y), param_grid_gbm_target)
  model_target_gbm <- gbm_target_model_info$model_target_gbm
  param_target_gmb <- gbm_target_model_info$best_target_parameters

  rf_target_model_info <- train_rf_model(combined_testing_data, combined_testing_y, param_grid_rf_target)
  rf_model <- rf_target_model_info$model_residual_rf

  predictions_result <- evaluate_and_get_predictions(
    best_model_source_nn = best_model_source_nn,
    best_model_source_rf = best_model_source_rf,
    best_model_source_lasso = best_model_source_lasso,
    best_model_source_gbm = best_model_source_gbm,
    param_source_gbm = model_result_source_gbm$best_parameters,
    param_transfer_gmb = param_transfer_gmb,
    param_target_gmb = param_target_gmb,
    lambda_source_lasso = lambda_source_lasso,
    model_residual = model_residual,
    model_residual_rf = model_residual_rf,
    model_residual_gbm = model_residual_gbm,
    model_tnn1 = model_tnn1,
    model_tnn2 = model_tnn2,
    model_target = model_target,
    rf_model = rf_model,
    model_target_gbm = model_target_gbm,
    test_data = final_testing_data,
    actuals = actuals,
    source_best_flag = source_best_flagreturn
  )

  return(list(predictions_list = predictions_result,
              source_best_flagreturn = source_best_flagreturn))
}
