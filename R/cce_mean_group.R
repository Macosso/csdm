#' Common Correlated Effects Mean Group (CCE-MG) Estimator
#'
#' @description
#' The CCE Mean Group (CCE-MG) estimator, introduced by Pesaran (2006), accounts
#' for cross-sectional dependence in panel data by augmenting regressions with
#' cross-sectional averages of dependent and independent variables.
#'
#' @details
#' Standard panel data estimators assume cross-sectional independence, which is often violated
#' in empirical applications due to unobserved common factors. The CCE-MG estimator
#' corrects for this by including cross-sectional averages as additional regressors.
#'
#' The estimation follows three key steps:
#' 1. Compute cross-sectional averages of the dependent and independent variables.
#' 2. Estimate individual id-level regressions with cross-sectional averages included.
#' 3. Compute mean group estimates by averaging individual coefficients.
#'
#' This estimator is robust to the presence of unobserved common factors, even if their
#' factor loadings vary across cross-sectional ids.
#'
#' Reference:
#' Pesaran, M. H. (2006). "Estimation and Inference in Large Heterogeneous Panels with a
#' Multifactor Error Structure," *Econometrica*, 74(4), 967–1012.
#' DOI: [10.1111/j.1468-0262.2006.00692.x](https://doi.org/10.1111/j.1468-0262.2006.00692.x)
#'
#' @import plm stats utils
#' @importFrom plm pdata.frame
#' @importFrom stats model.frame model.response model.matrix lm coef sd
#' @importFrom utils head tail
#' @param formula Model formula (e.g., `y ~ x1 + x2`).
#' @param data A `data.frame` or `pdata.frame` containing the panel data. If `pdata.frame` is provided,
#' id and time variables are inferred from `index`.
#' @param id Name of the cross-sectional identifier variable (required if `data` is a `data.frame`).
#' @param time Name of the time variable (required if `data` is a `data.frame`).
#' @param ... Additional options (e.g., robust standard errors).
#' @return An object of class `"cce_mean_group"` containing:
#' \item{mean_coefs}{Mean group estimated coefficients.}
#' \item{std_errors}{Standard errors of the mean group estimates.}
#' \item{individual_coefs}{Matrix of individual id-specific estimates.}
#' @export
cce_mean_group <- function(formula, data, id = NULL, time = NULL,
                           na.action = na.omit, ...) {

  if (!requireNamespace("plm", quietly = TRUE)) stop("Package 'plm' is required.")
  if (!requireNamespace("stats", quietly = TRUE)) stop("Package 'stats' is required.")

  if (!inherits(data, "pdata.frame")) {
    if (is.null(id) || is.null(time)) {
      stop("If data is not a pdata.frame, both 'id' and 'time' must be provided.")
    }
    data <- plm::pdata.frame(data, index = c(id, time))
  }

  clean_data <- na.action(data)

  model_frame <- model.frame(formula, clean_data, na.action = na.action)
  y <- model.response(model_frame)
  X <- model.matrix(formula, clean_data)

  ids <- unique(index(clean_data)[[1]])
  time_periods <- unique(index(clean_data)[[2]])

  avg_y <- tapply(y, index(clean_data)[[2]], mean)
  avg_X <- apply(X[, setdiff(colnames(X), "(Intercept)")], 2,
                 function(col) tapply(col, index(clean_data)[[2]], mean))

  coef_list <- list()
  rsquared <- c()
  s2mg <- c()
  s2 <- c()
  N <- length(ids)
  K <-  ncol(X)

  for (i in ids) {
    sub_data <- clean_data[index(clean_data)[[1]] == i, ]
    sub_data$avg_y <- avg_y[match(index(sub_data)[[2]], names(avg_y))]
    sub_data$avg_X <- avg_X[match(index(sub_data)[[2]], rownames(avg_X)), ]
    reg_formula <- update(formula, . ~ . + avg_y + avg_X)
    model <- lm(reg_formula, data = sub_data)
    coef_list[[as.character(i)]] <- coef(model)

    model_summary <- summary(model)
    rsquared[i] <- model_summary$r.squared

    r <- model$residuals
    TT <- nrow(sub_data)
    s2mg[i] <- (1/(N*(TT -1))*(t(r) %*% r))/((TT - 2*K -2)/N)
    s2[i] <- sum((sub_data$log_rgdpo - mean(sub_data$log_rgdpo))^2)
  }


  1 - sum(s2mg)/sum(s2)

  coef_matrix <- do.call(rbind, coef_list)

  mean_coefs <- colMeans(coef_matrix, na.rm = TRUE)
  vcvm <- var(coef_matrix) # Variance covariance matrix

  std_errors <- sqrt(diag(vcvm))/sqrt(nrow(coef_matrix)) # standard errors

  result <- list(
    call = match.call(),
    coefficients = mean_coefs,
    vcvm = vcvm,
    std_errors = std_errors,
    individual_coefs = coef_matrix,
    z_stat = mean_coefs/std_errors,
    p_value = (1-pnorm(abs(mean_coefs/std_errors)))*2,
    data = data,
    formula = formula,
    rsquared = mean(rsquared)
  )

  class(result) <- "cce_mean_group"
  return(result)
}


#' Print method for CCE-MG estimator
#' @param x An object of class 'cce_mean_group'.
#' @export
print.cce_mean_group <- function(x) {
  cat("\nCommon Correlated Effects Mean Group (CCE-MG) Estimator\n")
  print(data.frame(Estimate = x$coefficients, Std.Error = x$std_errors))
}



#' Summary method for CCE-MG estimator
#' @param object An object of class 'cce_mean_group'.
#' @importFrom dplyr mutate case_when
#' @export
summary.cce_mean_group <- function(object) {
  res <- data.frame(Estimate = object$coefficients,
             Std.Error = object$std_errors,
             z.value = object$z_stat,
             p.value = object$p_value)

  res <- res |>
    mutate(signi = case_when(p.value < 0.001 ~ "***",
                          p.value < 0.01 ~ "**",
                          p.value < 0.05 ~ "*",
                          p.value < 0.1 ~ ".",
                          .default = ""))

  cat("\nSummary of CCE-MG Estimator\n")
  cat("Panel Model Formula: ", deparse(object$formula), "\n")
  print(res)
  cat("------------------\n")
  cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
}



#' Coefficients method for CCE-MG estimator
#' @param object An object of class 'cce_mean_group'.
#' @export
coef.cce_mean_group <- function(object) {
  return(object$mean_coefs)
}


