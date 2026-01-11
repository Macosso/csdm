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

  tab <- data.frame(
    Estimate = as.numeric(est),
    Std.Error = as.numeric(se[names(est)]),
    z.value = as.numeric(z),
    p.value = as.numeric(p),
    check.names = FALSE
  )
  rownames(tab) <- names(est)

  out <- list(
    call = object$call,
    formula = object$formula,
    model = object$model,
    coefficients = tab,
    vcov = object$vcov_mg,
    meta = object$meta
  )
  class(out) <- "summary.csdm_fit"
  out
}


#' @export
print.summary.csdm_fit <- function(x, digits = 4, ...) {
  cat("csdm summary (", x$model, ")\n", sep = "")
  if (!is.null(x$formula)) cat("Formula: ", deparse(x$formula), "\n", sep = "")
  if (!is.null(x$meta$N) && !is.null(x$meta$T)) {
    cat("N = ", x$meta$N, ", T = ", x$meta$T, "\n\n", sep = "")
  }
  print(round(x$coefficients, digits))
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
predict.csdm_fit <- function(object, type = c("xb", "residuals"), ...) {
  type <- match.arg(type)
  if (type == "residuals") return(object$residuals_e)
  stop("predict(type='xb') not implemented yet (fit does not store model frame/newdata support).")
}
