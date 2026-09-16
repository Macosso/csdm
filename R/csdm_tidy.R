#' Tidy model coefficients, fit statistics, and observations
#' @param x A csdm_fit object.
#' @param component Parameter component, as in coef().
#' @param conf.int Include normal-approximation confidence intervals.
#' @param conf.level Confidence level between zero and one.
#' @param data Original fitting data, in its original order.
#' @param newdata Not supported.
#' @param ... Further arguments; currently unused.
#' @importFrom generics tidy glance augment
#' @name csdm_tidy
NULL

#' @rdname csdm_tidy
#' @export
tidy.csdm_fit <- function(x, component = c("all", "levels", "adjustment", "long_run"),
                         conf.int = FALSE, conf.level = 0.95, ...) {
  component <- match.arg(component)
  .csdm_flag(conf.int, "conf.int")
  if (!is.numeric(conf.level) || length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) stop("'conf.level' must be between zero and one.")
  b <- stats::coef(x, component = component)
  se <- sqrt(diag(stats::vcov(x, component = component)))
  statistic <- b / se
  statistic[!is.finite(se) | se <= 0] <- NA_real_
  out <- tibble::tibble(term = names(b), estimate = unname(b), std.error = unname(se),
    statistic = unname(statistic), p.value = unname(2 * stats::pnorm(abs(statistic), lower.tail = FALSE)))
  out$n_used <- if (is.null(x$components)) x$meta$N_used else x$components[[component]]$n_used
  if (conf.int) {
    critical <- stats::qnorm((1 + conf.level) / 2)
    out$conf.low <- unname(out$estimate - critical * se)
    out$conf.high <- unname(out$estimate + critical * se)
  }
  out
}

#' @rdname csdm_tidy
#' @export
glance.csdm_fit <- function(x, ...) {
  tibble::tibble(nobs = stats::nobs(x), n_units = x$meta$N_used,
    n_units_observed = x$meta$N_observed, n_periods = x$meta$T,
    r.squared = x$stats$R2_ols_mg, adj.r.squared = x$stats$R2_mg,
    cd.statistic = x$stats$CD_stat, cd.p.value = x$stats$CD_p)
}

#' @rdname csdm_tidy
#' @export
augment.csdm_fit <- function(x, data = x$data, newdata = NULL, ...) {
  if (!is.null(newdata)) stop("Augmentation of new data is not implemented.")
  if (!identical(as.data.frame(data), x$data)) stop("'data' must equal the original fitting data in its original order.")
  out <- tibble::as_tibble(data)
  if (any(c(".row", ".fitted", ".resid", ".used") %in% names(out))) stop("Augmentation column names already exist in data.")
  out$.row <- seq_len(nrow(out))
  out$.fitted <- unname(stats::fitted(x, format = "vector"))
  out$.resid <- unname(stats::residuals(x, format = "vector"))
  out$.used <- is.finite(out$.resid)
  out
}
