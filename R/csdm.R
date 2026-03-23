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
#' \deqn{
#' y_{it} = \alpha_i + x_{it}^{\top}\beta_i + u_{it},
#' }
#'
#' where \eqn{\alpha_i} is a unit-specific intercept, \eqn{x_{it}} is a vector
#' of regressors, \eqn{\beta_i} is a vector of unit-specific slopes, and
#' \eqn{u_{it}} is an error term that may exhibit cross-sectional dependence.
#' The inner product \eqn{x_{it}^{\top}\beta_i} is scalar-valued. The estimators
#' implemented in \code{csdm()} differ in how they handle slope heterogeneity,
#' common factors, and dynamic adjustment.
#'
#' \strong{Cross-sectional averages and dynamic structure}
#'
#' Cross-sectional averages are specified through [csdm_csa()] and dynamic or
#' long-run structure is specified through [csdm_lr()]. This keeps the model
#' interface consistent across estimators while allowing the degree of
#' cross-sectional augmentation and lag structure to vary by application.
#'
#' \strong{Implemented estimators}
#'
#' \describe{
#'   \item{MG (Pesaran and Smith, 1995)}{
#'     Unit-by-unit estimation with heterogeneous slopes:
#'
#'     \deqn{
#'     y_{it} = \alpha_i + x_{it}^{\top}\beta_i + u_{it}.
#'     }
#'
#'     The reported coefficients are cross-sectional averages of the
#'     unit-specific estimates:
#'
#'     \deqn{
#'     \hat{\beta}_{MG} = \frac{1}{N}\sum_{i=1}^N \hat{\beta}_i.
#'     }
#'
#'     This estimator accommodates slope heterogeneity but does not explicitly
#'     model cross-sectional dependence.
#'   }
#'   \item{CCE (Pesaran, 2006)}{
#'     Regressions are augmented with cross-sectional averages to proxy
#'     unobserved common factors:
#'
#'     \deqn{
#'     y_{it} = \alpha_i + x_{it}^{\top}\beta_i + \bar{z}_{t}^{\top}\gamma_i + v_{it},
#'     }
#'
#'     where \eqn{\bar{z}_t} collects the cross-sectional averages specified in
#'     \code{csa}, for example averages of the dependent variable and
#'     regressors. This estimator is suitable when cross-sectional dependence is
#'     driven by latent common factors.
#'   }
#'   \item{DCCE (Chudik and Pesaran, 2015)}{
#'     Dynamic CCE extends CCE by allowing lagged dependent variables and lagged
#'     cross-sectional averages:
#'
#'     \deqn{
#'     y_{it} =
#'     \alpha_i
#'     + \sum_{p=1}^{P} \phi_{ip} y_{i,t-p}
#'     + x_{it}^{\top}\beta_i
#'     + \sum_{q=0}^{Q} \bar{z}_{t-q}^{\top}\delta_{iq}
#'     + e_{it}.
#'     }
#'
#'     The lag structure is controlled through \code{lr} and \code{csa}. This
#'     estimator is designed for dynamic panels with persistent outcomes and
#'     common factor dependence.
#'   }
#'   \item{CS-ARDL (Chudik and Pesaran, 2015)}{
#'     Cross-sectionally augmented ARDL combines distributed-lag dynamics with
#'     cross-sectional augmentation. A convenient error-correction
#'     representation is
#'
#'     \deqn{
#'     \Delta y_{it}
#'     =
#'     \alpha_i
#'     + \phi_i \left( y_{i,t-1} - \theta_i^{\top} x_{i,t-1} \right)
#'     + \sum_{j=1}^{P-1} \lambda_{ij} \Delta y_{i,t-j}
#'     + \sum_{j=0}^{Q-1} \psi_{ij}^{\top} \Delta x_{i,t-j}
#'     + \sum_{s=0}^{S} \bar{z}_{t-s}^{\top} \omega_{is}
#'     + e_{it}.
#'     }
#'
#'     Here \eqn{\theta_i} denotes the unit-specific long-run relationship,
#'     \eqn{\phi_i} is the speed of adjustment, and the remaining terms capture
#'     short-run dynamics and common cross-sectional components. The lag
#'     structure is controlled through \code{lr} and \code{csa}.
#'   }
#' }
#'
#' \strong{Identification and assumptions}
#'
#' \describe{
#'   \item{MG}{
#'   Requires sufficient time-series variation within each unit. Consistency is
#'   based on unit-by-unit estimation and averaging across units.}
#'   \item{CCE}{
#'   Identification relies on cross-sectional averages acting as proxies for
#'   latent common factors, together with adequate cross-sectional and time
#'   dimensions and weak dependence in the remaining idiosyncratic component.}
#'   \item{DCCE}{
#'   In addition to the CCE assumptions, dynamic identification requires enough
#'   time periods to support lagged dependent variables and lagged
#'   cross-sectional averages.}
#'   \item{CS-ARDL}{
#'   Requires sufficient time length for the distributed-lag structure and is
#'   intended for applications where both short-run dynamics and long-run
#'   relationships are of interest in the presence of common factors.}
#' }
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
