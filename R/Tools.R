#' Compute Variance-Covariance Matrix
#'
#' @description
#' Computes the variance-covariance matrix for a given matrix of estimated coefficients.
#'
#' @param coef_matrix A numeric matrix where each row corresponds to individual estimates of parameters.
#' @param method A character string indicating the method to compute the covariance matrix. Options are "standard" (default) or "robust".
#' @return A variance-covariance matrix.
#' @export
compute_vcov <- function(coef_matrix, method = "standard") {

  if (!is.matrix(coef_matrix)) {
    stop("coef_matrix must be a numeric matrix.")
  }

  n <- nrow(coef_matrix)
  mean_coefs <- colMeans(coef_matrix, na.rm = TRUE)

  if (method == "standard") {
    vcov_matrix <- cov(coef_matrix, use = "complete.obs") / n
  } else if (method == "robust") {
    deviations <- sweep(coef_matrix, 2, mean_coefs, "-")
    vcov_matrix <- t(deviations) %*% deviations / (n - 1)
  } else {
    stop("Invalid method. Choose either 'standard' or 'robust'.")
  }

  return(vcov_matrix)
}

#' Predict Y using estimated coefficients
#'
#' @description
#' Predicts the dependent variable using estimated coefficients and input features.
#'
#' @param coef_vector A numeric vector of estimated coefficients.
#' @param X_matrix A numeric matrix of independent variables (including intercept if applicable).
#' @return A numeric vector of predicted values.
#' @export
#'
predict_y <- function(coef_vector, X_matrix) {

  coef_vector <- coef_vector[!is.na(coef_vector)]
  X_matrix <- X_matrix[, names(coef_vector)]

  if (!is.vector(coef_vector) || !is.numeric(coef_vector)) {
    stop("coef_vector must be a numeric vector.")
  }

  if (!is.matrix(X_matrix) || !is.numeric(X_matrix)) {
    stop("X_matrix must be a numeric matrix.")
  }

  if (ncol(X_matrix) != length(coef_vector)) {
    stop("Number of columns in X_matrix must match the length of coef_vector.")
  }

  return(X_matrix %*% coef_vector)
}

