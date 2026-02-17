# csdm.R

#' Panel Model Estimation with Cross Section Dependence
#'
#' Estimate panel data models that allow for cross-sectional dependence and
#' heterogeneous slopes. The interface supports Mean Group (MG), Common
#' Correlated Effects (CCE), Dynamic CCE (DCCE), and Cross-Sectionally
#' Augmented ARDL (CS-ARDL) estimators with consistent handling of
#' cross-sectional averages, dynamic structure, and robust inference.
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
#'   Use \code{summary()}, \code{coef()}, \code{residuals()}, \code{vcov()}, and
#'   \code{cd_test()} to access standard outputs.
#'
#' @details
#' ## Model equations
#'
#' \describe{
#'   \item{MG (Pesaran and Smith, 1995)}{
#'     \deqn{y_{it} = x_{it}^\top \beta_i + u_{it}}
#'   }
#'   \item{CCE (Pesaran, 2006)}{
#'     \deqn{y_{it} = x_{it}^\top \beta_i + \lambda_i^\top F_t + u_{it}}
#'   }
#'   \item{DCCE (Chudik and Pesaran, 2015)}{
#'     \deqn{\Delta y_{it} = \Delta x_{it}^\top \beta_i + \lambda_i^\top \Delta F_t + u_{it}}
#'   }
#'   \item{CS-ARDL (Chudik and Pesaran, 2015)}{
#'     \deqn{y_{it} = \phi_i y_{it-1} + x_{it}^\top \theta_i + \lambda_i^\top F_t + u_{it}}
#'   }
#' }
#'
#' ## Estimation, identification, and assumptions
#'
#' \describe{
#'   \item{MG}{Unit-by-unit estimation with heterogeneous slopes. The reported
#'   coefficients are cross-sectional averages of unit estimates. Requires
#'   sufficient time series per unit and weak serial dependence in errors.}
#'   \item{CCE}{Augments regressions with cross-sectional averages (CSA) to proxy
#'   unobserved common factors. Identification relies on large N and T, weak
#'   dependence in idiosyncratic errors after CSA, and weak exogeneity of
#'   regressors.}
#'   \item{DCCE}{Extends CCE to dynamic settings with lagged dependent variables
#'   and CSA lags. Identification relies on weak exogeneity, adequate time length
#'   for dynamic lags, and a stable factor structure.}
#'   \item{CS-ARDL}{Specifies dynamic distributed lags with CSA terms. Estimation
#'   follows ARDL-style dynamics in each unit and aggregates to panel averages.
#'   Assumes weak exogeneity and sufficient time length for lag structure.}
#' }
#'
#' @references
#' Pesaran, M.H. and Smith, R. (1995). "Estimating long-run relationships from
#' dynamic heterogeneous panels." Journal of Econometrics, 68(1), 79-113.
#'
#' Pesaran, M.H. (2006). "Estimation and inference in large heterogeneous panels
#' with multifactor error structure." Econometrica, 74(4), 967-1012.
#'
#' Chudik, A. and Pesaran, M.H. (2015). "Common correlated effects estimation of
#' heterogeneous dynamic panel data models with weakly exogenous regressors."
#' Journal of Econometrics, 188(2), 393-420.
#'
#' Chudik, A. and Pesaran, M.H. (2015). "Large panel data models with
#' cross-sectional dependence: A survey." Annals of Economics and Finance, 16(1),
#' 53-78.
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
