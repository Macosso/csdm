# csdm.R

#' Panel Model Estimation with Cross-Sectional Dependence
#'
#' Estimate heterogeneous panel data models with optional cross-sectional
#' augmentation and dynamic structure. The interface supports Mean Group (MG),
#' Common Correlated Effects (CCE), Dynamic CCE (DCCE), and
#' Cross-Sectionally Augmented ARDL (CS-ARDL) estimators with a consistent
#' specification workflow for cross-sectional averages, lag structure, and
#' variance-covariance estimation.
#'
#' @param formula Model formula of the form \code{y ~ x1 + x2}.
#' @param data A \code{data.frame} (or \code{plm::pdata.frame}) containing the
#'   variables in \code{formula}.
#' @param id,time Column names (strings) for the unit and time indexes. If
#'   \code{data} is a \code{pdata.frame}, these are taken from its index and the
#'   provided values are ignored.
#' @param model Estimator to fit. One of \code{"mg"}, \code{"cce"},
#'   \code{"dcce"}, or \code{"cs_ardl"}.
#' @param csa Cross-sectional-average specification, created by [csdm_csa()].
#' @param lr Long-run or dynamic specification, created by [csdm_lr()].
#' @param pooled Pooled specification (reserved for future use), created by
#'   [csdm_pooled()].
#' @param trend One of \code{"none"} or \code{"unit"} (adds a linear unit trend).
#'   \code{"pooled"} is reserved and not implemented.
#' @param fullsample Logical; reserved for future extensions.
#' @param mgmissing Logical; reserved for future extensions.
#' @param vcov Variance-covariance specification, created by [csdm_vcov()].
#' @param ... Reserved for future extensions.
#'
#' @return An object of class \code{csdm_fit} containing estimated coefficients,
#'   residuals, variance-covariance estimates, model metadata, and diagnostics.
#'   Use [summary()], [coef()], [residuals()], [vcov()], and
#'   [cd_test()] to access standard outputs.
#'
#' @details
#' Let \eqn{i = 1, \ldots, N} index cross-sectional units and
#' \eqn{t = 1, \ldots, T} index time. A baseline heterogeneous panel model is
#'
#' \deqn{y_{it} = \alpha_i + \beta_i^T x_{it} + u_{it}.}
#'
#' Here \eqn{\alpha_i} is a unit-specific intercept, \eqn{x_{it}} is a vector
#' of regressors, \eqn{\beta_i} is a vector of unit-specific slopes, and
#' \eqn{u_{it}} is an error term that may exhibit cross-sectional dependence.
#'
#' Cross-sectional averages are specified through [csdm_csa()] and dynamic or
#' long-run structure is specified through [csdm_lr()]. This keeps the model
#' interface consistent across estimators while allowing the degree of
#' cross-sectional augmentation and lag structure to vary by application.
#'
#' \strong{Implemented estimators}
#'
#' \strong{MG (Pesaran and Smith, 1995)}
#'
#' The Mean Group estimator fits separate regressions for each unit and averages
#' the resulting coefficients:
#'
#' \deqn{\hat{\beta}_{MG} = \frac{1}{N}\sum_{i=1}^N \hat{\beta}_i.}
#'
#' This estimator accommodates slope heterogeneity but does not explicitly model
#' cross-sectional dependence.
#'
#' \strong{CCE (Pesaran, 2006)}
#'
#' Regressions are augmented with cross-sectional averages to proxy unobserved
#' common factors:
#'
#' \deqn{y_{it} = \alpha_i + \beta_i^T x_{it} + \gamma_i^T \bar{z}_{t} + v_{it}.}
#'
#' A common choice is
#'
#' \deqn{\bar{z}_t = (\bar{y}_t, \bar{x}_t),}
#'
#' with
#'
#' \deqn{\bar{x}_t = \frac{1}{N}\sum_{i=1}^N x_{it}, \qquad
#' \bar{y}_t = \frac{1}{N}\sum_{i=1}^N y_{it}.}
#'
#' More generally, \eqn{\bar{z}_t} collects the cross-sectional averages
#' specified in \code{csa}.
#'
#' \strong{DCCE (Chudik and Pesaran, 2015)}
#'
#' Dynamic CCE extends CCE by allowing lagged dependent variables and lagged
#' cross-sectional averages:
#'
#' \deqn{y_{it} = \alpha_i + \sum_{p=1}^{P} \phi_{ip} y_{i,t-p}
#' + \beta_i^T x_{it}
#' + \sum_{q=0}^{Q} \delta_{iq}^T \bar{z}_{t-q}
#' + e_{it}.}
#'
#' In the package implementation, lagged dependent variables and distributed
#' lags of regressors are controlled through \code{lr}, while contemporaneous
#' and lagged cross-sectional averages are controlled through \code{csa}.
#'
#' \strong{CS-ARDL (Chudik and Pesaran, 2015)}
#'
#' In the package implementation, \code{model = "cs_ardl"} is obtained by first
#' estimating a cross-sectionally augmented ARDL-style regression in levels,
#' using the same dynamic specification as \code{model = "dcce"}, and then
#' transforming the unit-specific coefficients into adjustment and long-run
#' parameters.
#'
#' The underlying unit-level regression is of the form
#'
#' \deqn{y_{it} = \alpha_i + \sum_{p=1}^{P} \phi_{ip} y_{i,t-p}
#' + \sum_{q=0}^{Q} \beta_{iq}^T x_{i,t-q}
#' + \sum_{s=0}^{S} \omega_{is}^T \bar{z}_{t-s}
#' + e_{it}.}
#'
#' From this dynamic specification, the package recovers the implied
#' error-correction form
#'
#' \deqn{\Delta y_{it} =
#' \alpha_i +
#' \varphi_i \left(y_{i,t-1} - \theta_i^T x_{i,t-1}\right)
#' + \sum_{j=1}^{P-1} \lambda_{ij} \Delta y_{i,t-j}
#' + \sum_{j=0}^{Q-1} \psi_{ij}^T \Delta x_{i,t-j}
#' + \sum_{s=0}^{S} \tilde{\omega}_{is}^T \bar{z}_{t-s}
#' + e_{it},}
#'
#' where \eqn{\varphi_i} is the adjustment coefficient and \eqn{\theta_i} is
#' the implied long-run relationship. In the current implementation, these
#' quantities are computed from the estimated lag polynomials rather than from a
#' direct ECM regression.
#'
#' \strong{Identification and assumptions}
#'
#' MG requires sufficient time-series variation within each unit.
#'
#' CCE relies on cross-sectional averages acting as proxies for latent common
#' factors, together with adequate cross-sectional and time dimensions.
#'
#' DCCE additionally requires enough time periods to support lagged dependent
#' variables, distributed lags, and lagged cross-sectional averages.
#'
#' CS-ARDL requires sufficient time length for the distributed-lag structure and
#' is intended for applications where both short-run dynamics and long-run
#' relationships are of interest in the presence of common factors.
#'
#' @references
#' \insertRef{PesaranSmith1995}{csdm}
#'
#' \insertRef{Pesaran2006}{csdm}
#'
#' \insertRef{ChudikPesaran2015a}{csdm}
#'
#' @examples
#' library(csdm)
#' data(PWT_60_07, package = "csdm")
#' df <- PWT_60_07
#'
#' # Keep examples fast but fully runnable
#' keep_ids <- unique(df$id)[1:10]
#' df_small <- df[df$id %in% keep_ids & df$year >= 1970, ]
#'
#' # Mean Group (MG)
#' mg <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small, id = "id", time = "year", model = "mg"
#' )
#' summary(mg)
#'
#' # Common Correlated Effects (CCE)
#' cce <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small, id = "id", time = "year", model = "cce",
#'   csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"))
#' )
#' summary(cce)
#'
#' # Dynamic CCE (DCCE)
#' dcce <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small, id = "id", time = "year", model = "dcce",
#'   csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), lags = 3),
#'   lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 0)
#' )
#' summary(dcce)
#'
#' # CS-ARDL
#' cs_ardl <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small, id = "id", time = "year", model = "cs_ardl",
#'   csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), lags = 3),
#'   lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 0)
#' )
#' summary(cs_ardl)
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
    cs_ardl = .csdm_fit_cs_ardl(panel_df = panel_df, formula = formula, id = id, time = time, csa = csa, lr = lr, vcov = vcov, ...),
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
