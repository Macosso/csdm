# csdm.R

#' Unified front door for csdm estimators
#'
#' @param formula A model formula like `y ~ x1 + x2`.
#' @param data A `data.frame` (or `plm::pdata.frame`) containing the variables in `formula`.
#' @param id,time Column names (strings) for the unit and time indexes.
#'   If `data` is a `pdata.frame`, these are taken from its index and the provided values are ignored.
#' @param model Which estimator to fit. Currently implemented: `"mg"`, `"cce"`, `"dcce"`.
#' @param csa Cross-sectional-average specification, created by [csdm_csa()].
#' @param lr Long-run specification (currently a stub), created by [csdm_lr()].
#' @param pooled Pooled specification (currently a stub), created by [csdm_pooled()].
#' @param trend One of `"none"` or `"unit"` (a linear trend per unit is added as a regressor).
#'   `"pooled"` is reserved and not implemented.
#' @param fullsample Logical; reserved for future extensions.
#' @param mgmissing Logical; reserved for future extensions.
#' @param vcov Variance-covariance specification, created by [csdm_vcov()].
#' @param ... Reserved for future extensions.
#'
#' @return An object of class `csdm_fit`.
#' @export
csdm <- function(
  formula, data, id, time,
  model = c("mg", "cce", "dcce", "cs_ardl", "cs_ecm", "cs_dl"),
  csa = csdm_csa(),
  lr = csdm_lr(),
  pooled = csdm_pooled(),
  trend = c("none", "unit", "pooled"),
  fullsample = FALSE,
  mgmissing = FALSE,
  vcov = csdm_vcov(),
  ...
) {
  model <- match.arg(model)
  trend <- match.arg(trend)

  if (inherits(data, "pdata.frame")) {
    idx <- attr(data, "index")
    if (!is.null(idx)) {
      id <- names(idx)[1L]
      time <- names(idx)[2L]
    }
  }

  panel_df <- .csdm_prepare_panel_df(data = data, id = id, time = time)

  if (trend == "pooled") {
    stop("trend='pooled' is not implemented yet")
  }
  if (trend == "unit") {
    panel_df$.csdm_trend__ <- .csdm_time_index(panel_df[[time]])
    formula <- stats::update(formula, . ~ . + .csdm_trend__)
  }

  fit <- switch(
    model,
    mg = .csdm_fit_mg(panel_df = panel_df, formula = formula, id = id, time = time, lr = lr, vcov = vcov, ...),
    cce = .csdm_fit_cce(panel_df = panel_df, formula = formula, id = id, time = time, csa = csa, lr = lr, vcov = vcov, ...),
    dcce = .csdm_fit_dcce(panel_df = panel_df, formula = formula, id = id, time = time, csa = csa, lr = lr, vcov = vcov, ...),
    cs_ardl = stop("Not implemented yet"),
    cs_ecm  = stop("Not implemented yet"),
    cs_dl   = stop("Not implemented yet")
  )

  fit$call <- match.call()
  fit$formula <- formula
  fit$model <- model
  fit$id <- id
  fit$time <- time
  fit$meta$trend <- trend
  fit$meta$fullsample <- isTRUE(fullsample)
  fit$meta$mgmissing <- isTRUE(mgmissing)
  fit$meta$lr <- lr
  fit$meta$pooled <- pooled
  fit$meta$vcov <- vcov

  class(fit) <- "csdm_fit"
  fit
}
