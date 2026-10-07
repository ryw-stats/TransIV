# Shared functions; source this file from a workflow before calling them.

#' Construct the correlation matrix used by the original generators.
#'
#' This preserves the active second definition from the notebooks: it is an equicorrelation
#' matrix, not an AR(1) matrix.
#'
#' @param nc Number of rows and columns.
#' @param rc Common off-diagonal correlation.
#' @return Square matrix with ones on the diagonal and rc elsewhere.
create_toeplitz_matrix <- function(nc, rc) {
  toeplitz(c(1, rep(rc, nc - 1)))
}

#' Calculate residual sum of squares.
#'
#' @param y_true Observed response vector.
#' @param y_pred Predictions aligned with y_true.
#' @return Scalar sum of squared prediction errors.
calculate_rss <- function(y_true, y_pred) {
  sum((y_true - y_pred) ^ 2)
}

#' Fit a cross-validated lasso and refit at lambda.min.
#'
#' The optional lambda_min field retains its historical name but contains lambda.1se.
#' coefficients2 uses the same single-lambda refit as the original nonlinear notebook. Cross-
#' validation uses the caller's RNG state.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param residuals Response or residual vector aligned with x0.
#' @param penalty.factor Per-predictor penalty multipliers; zero leaves a predictor unpenalized.
#' @param include_1se Include the historical nonlinear-workflow fields for lambda.1se.
#' @return List with model_residual, lambda_residual, and coefficients; optionally lambda_min and coefficients2.
train_lasso_model <- function(x0, residuals, penalty.factor = rep(1, ncol(x0)), include_1se = FALSE) {
  x0_matrix <- as.matrix(x0)
  residuals_matrix <- as.matrix(residuals)

  colnames(x0_matrix) <- paste0("V", 1:ncol(x0_matrix))

  cv_lasso <- glmnet::cv.glmnet(x0_matrix, residuals_matrix, alpha = 1, penalty.factor = penalty.factor)

  best_lasso_model <- glmnet::glmnet(x0_matrix, residuals_matrix, alpha = 1, lambda = cv_lasso$lambda.min, penalty.factor = penalty.factor)

  lambda_residual <- cv_lasso$lambda.min

  cat("Best lambda:", lambda_residual, "\n")

  coefficients <- coef(best_lasso_model, s = cv_lasso$lambda.min)

  result <- list(model_residual = best_lasso_model, lambda_residual = lambda_residual, coefficients = coefficients)
  if (include_1se) {
    result <- list(model_residual = best_lasso_model, lambda_residual = lambda_residual,
                   lambda_min = cv_lasso$lambda.1se, coefficients = coefficients,
                   coefficients2 = coef(best_lasso_model, s = cv_lasso$lambda.1se))
  }
  return(result)
}

#' Fit a cross-validated target-only lasso.
#'
#' @param x0 Numeric predictor matrix or data frame; rows are observations.
#' @param residuals Response or residual vector aligned with x0.
#' @return List with model_residual (cv.glmnet fit) and lambda_residual.
train_target_lasso_model <- function(x0, residuals) {
  x0_matrix <- as.matrix(x0)
  residuals_matrix <- as.matrix(residuals)

  colnames(x0_matrix) <- paste0("V", 1:ncol(x0_matrix))

  cv_lasso <- glmnet::cv.glmnet(x0_matrix, residuals_matrix, alpha = 1)

  return(list(model_residual = cv_lasso, lambda_residual = cv_lasso$lambda.min))
}

#' Build and compile the source neural network.
#'
#' @param input_shape Input dimensions accepted by the Keras model.
#' @param units Number of hidden units in each source network layer.
#' @param dropout_rate Dropout proportion for hidden layers.
#' @param l2_rate L2 regularization strength for dense-layer weights.
#' @param patience Epochs without validation improvement before early stopping.
#' @return List with the compiled model and early_stopping callback.
define_compile_model <- function(input_shape, units = c(256, 256, 128, 128, 64), dropout_rate = 0.1, l2_rate = 0.001, patience = 20) {

  l2_regularizer <- tf$keras$regularizers$L2(l2_rate)

  model <- keras_model_sequential()
  model$add(layer_dense(units = units[1], activation = 'relu', input_shape = input_shape, kernel_regularizer = l2_regularizer))
  model$add(layer_dropout(rate = dropout_rate))

  if (length(units) > 1) {
    for (unit in units[-1]) {
      model$add(layer_dense(units = unit, activation = 'relu', kernel_regularizer = l2_regularizer))
      model$add(layer_dropout(rate = dropout_rate))
    }
  }

  model$add(layer_dense(units = 1))

  model$compile(
    optimizer = 'adam',
    loss = 'mse',
    metrics = list('mean_absolute_error')
  )

  early_stopping <- callback_early_stopping(
    monitor = 'val_loss',
    patience = patience,
    restore_best_weights = TRUE
  )

  return(list(model = model, early_stopping = early_stopping))
}

#' Build and compile the target neural network.
#'
#' @param input_shape Input dimensions accepted by the Keras model.
#' @param units_list Number of hidden units in each target network layer.
#' @param dropout_rate Dropout proportion for hidden layers.
#' @param l2_rate L2 regularization strength for dense-layer weights.
#' @return List with the compiled model and early_stopping callback.
define_compile_target_model <- function(input_shape, units_list = c(64, 32), dropout_rate = 0.1, l2_rate = 0.001) {
  l2_regularizer <- tf$keras$regularizers$L2(l2_rate)

  model_target <- keras_model_sequential()
  model_target$add(layer_dense(units = units_list[1], activation = 'relu', input_shape = input_shape, kernel_regularizer = l2_regularizer))

  for (unit in units_list[-1]) {
    model_target$add(layer_dense(units = unit, activation = 'relu', kernel_regularizer = l2_regularizer))
    model_target$add(layer_dropout(rate = dropout_rate))
  }

  model_target$add(layer_dense(units = 1))

  model_target$compile(optimizer = 'adam', loss = 'mse', metrics = list('mean_absolute_error'))

  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 20, restore_best_weights = TRUE)

  return(list(model = model_target, early_stopping = early_stopping))
}
