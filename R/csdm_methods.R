# csdm_methods.R

#' @export
print.csdm_fit <- function(x, digits = 4, ...) {
  cat("csdm fit (", x$model, ")\n", sep = "")
  if (!is.null(x$formula)) cat("Formula: ", deparse(x$formula), "\n", sep = "")
  if (!is.null(x$meta$N) && !is.null(x$meta$T)) {
    cat("N = ", x$meta$N, ", T = ", x$meta$T, "\n", sep = "")
  }

  if (length(x$coef_mg)) {
    est <- x$coef_mg
    se  <- x$se_mg
    tab <- data.frame(
      Estimate = as.numeric(est),
      Std.Error = as.numeric(se[names(est)]),
      check.names = FALSE
    )
    rownames(tab) <- names(est)
    print(round(tab, digits))
  }
  invisible(x)
}


#' @export
summary.csdm_fit <- function(object, digits = 4, ...) {
  est <- object$coef_mg
  se  <- object$se_mg
  z   <- est / se[names(est)]
  p   <- 2 * (1 - stats::pnorm(abs(z)))

  coef_table <- data.frame(
    term = names(est),
    estimate = as.numeric(est),
    se = as.numeric(se[names(est)]),
    z = as.numeric(z),
    p_value = as.numeric(p),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  out <- list(
    call = object$call,
    formula = object$formula,
    model = object$model,
    N = object$meta$N,
    T = object$meta$T,
    coef_table = coef_table
  )
  class(out) <- "summary.csdm_fit"
  out
}


#' @export
print.summary.csdm_fit <- function(x, digits = 4, ...) {
  cat("csdm summary (", x$model, ")\n", sep = "")
  if (!is.null(x$formula)) cat("Formula: ", deparse(x$formula), "\n", sep = "")
  if (!is.null(x$N) && !is.null(x$T)) {
    cat("N = ", x$N, ", T = ", x$T, "\n\n", sep = "")
  }
  tab <- x$coef_table
  num_cols <- intersect(c("estimate", "se", "z", "p_value"), names(tab))
  tab[num_cols] <- lapply(tab[num_cols], function(col) round(col, digits))
  print(tab, row.names = FALSE)
  invisible(x)
}


#' @export
coef.csdm_fit <- function(object, ...) {
  object$coef_mg
}


#' @export
vcov.csdm_fit <- function(object, ...) {
  object$vcov_mg
}


#' @export
residuals.csdm_fit <- function(object, type = c("e", "u"), ...) {
  type <- match.arg(type)
  if (type == "u") stop("residuals(type='u') not implemented yet")
  object$residuals_e
}


#' @export
predict.csdm_fit <- function(object, newdata = NULL, type = c("xb", "residuals"), ...) {
  type <- match.arg(type)
  if (!is.null(newdata)) stop("predict(newdata=...) not implemented yet")
  if (type == "residuals") return(stats::residuals(object, type = "e"))
  if (!is.null(object$fitted_xb)) return(object$fitted_xb)
  stop("predict(type='xb') not implemented yet")
}
