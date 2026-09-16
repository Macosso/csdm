# csdm_specs.R

#' Specification: Cross-sectional averages (CSA)
#'
#' @param vars Character. One of "_all", "_none", or a character vector of variable names.
#' @param lags Integer. Either a scalar integer >= 0 applied to all CSA variables,
#'   or a named integer vector giving per-variable maximum lags.
#'   Named specifications apply only to named variables; other CSA lags are zero.
#' @param scope CSA sample scope. Must be `"estimation"`.
#' @param cluster Must be `NULL`; clustered CSA construction is not implemented.
#'
#' @return A spec object (list) used by csdm().
#' @export
#' @examples
#' # Cross-sectional averages (CSA) configuration for DCCE
#' csa <- csdm_csa(
#'   vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"),
#'   lags = 3
#' )
#' csa
csdm_csa <- function(
  vars = "_all",
  lags = 0,
  scope = "estimation",
  cluster = NULL
) {
  # vars
  if (is.null(vars) || length(vars) == 0L) {
    stop("'vars' must be '_all', '_none', or a character vector.")
  }
  if (length(vars) == 1L && vars %in% c("_all", "_none")) {
    vars <- as.character(vars)
  } else {
    if (!is.character(vars)) stop("'vars' must be a character vector.")
    if (anyNA(vars) || any(!nzchar(vars))) stop("'vars' contains missing or empty strings.")
    vars <- unique(as.character(vars))
  }

  lags <- .csdm_integer(lags, "lags", scalar = FALSE)
  if (length(lags) > 1L && is.null(names(lags))) stop("Multiple CSA lags must be named.")
  if (!is.null(names(lags)) && (anyNA(names(lags)) || any(!nzchar(names(lags))) || anyDuplicated(names(lags)))) {
    stop("CSA lag names must be nonempty, unique, and nonmissing.")
  }

  if (!is.character(scope) || length(scope) != 1L || is.na(scope) || scope != "estimation") {
    stop("Only scope='estimation' is implemented.", call. = FALSE)
  }
  if (!is.null(cluster)) {
    stop("'cluster' is not implemented and must be NULL.", call. = FALSE)
  }

  spec <- list(
    vars = vars,
    lags = lags,
    scope = scope,
    cluster = cluster
  )
  class(spec) <- "csdm_csa_spec"
  spec
}


#' Specification: Long-run configuration
#'
#' @param vars Must be `NULL`; variable-specific long-run restrictions are not
#'   implemented.
#' @param type Either `"none"` or `"ardl"`.
#' @param ylags Integer >= 0. Within-unit lags of the dependent variable to include
#'   when supported by the chosen model/type.
#' @param xdlags Integer >= 0. Scalar distributed lags to apply to each RHS regressor
#'   when supported by the chosen model/type.
#' @param options Must be an empty list; additional long-run options are not
#'   implemented.
#'
#' @return A spec object (list) used by csdm().
#' @export
#' @examples
#' # Long-run / dynamic configuration (ARDL-style lags)
#' lr <- csdm_lr(type = "ardl", ylags = 1)
#' lr
#'
#' # Minimal end-to-end DCCE example (kept small for speed)
#' data(PWT_60_07, package = "csdm")
#' df <- PWT_60_07
#' keep_ids <- unique(df$id)[1:10]
#' df_small <- df[df$id %in% keep_ids & df$year >= 1970, ]
#' fit <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small,
#'   id = "id",
#'   time = "year",
#'   model = "dcce",
#'   csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), lags = 3),
#'   lr = csdm_lr(type = "ardl", ylags = 1)
#' )
#' summary(fit)
csdm_lr <- function(vars = NULL,
                    type = c("none", "ardl"),
                    ylags = 0,
                    xdlags = 0,
                    options = list()) {
  type <- match.arg(type)

  if (!is.null(vars)) {
    stop("'vars' is not implemented and must be NULL.", call. = FALSE)
  }
  if (!identical(options, list())) {
    stop("'options' is not implemented and must be an empty list.", call. = FALSE)
  }

  ylags <- .csdm_integer(ylags, "ylags")
  xdlags <- .csdm_integer(xdlags, "xdlags")

  spec <- list(vars = vars, type = type, ylags = ylags, xdlags = xdlags, options = options)
  class(spec) <- "csdm_lr_spec"
  spec
}


#' Deprecated pooled-constraint specification
#'
#' `csdm_pooled()` is deprecated because pooled restrictions are not implemented.
#' Explicit pooled specifications continue to be rejected by [csdm()].
#'
#' @param vars Deprecated; formerly reserved for pooled variables.
#' @param constant Deprecated logical pooled-constant indicator.
#' @param trend Deprecated logical pooled-trend indicator.
#'
#' @return A deprecated specification object retained for compatibility.
#' @export
csdm_pooled <- function(vars = NULL, constant = FALSE, trend = FALSE) {
  .Deprecated(
    package = "csdm",
    old = "csdm_pooled",
    msg = "'csdm_pooled()' is deprecated because pooled restrictions are not implemented."
  )
  .csdm_pooled_spec(vars, constant, trend)
}

.csdm_pooled_spec <- function(vars = NULL, constant = FALSE, trend = FALSE) {
  spec <- list(vars = vars, constant = .csdm_flag(constant, "constant"), trend = .csdm_flag(trend, "trend"))
  class(spec) <- "csdm_pooled_spec"
  spec
}


#' Specification: Mean-group variance-covariance estimator
#'
#' @param type Must be `"mg"`, the implemented mean-group variance estimator.
#' @param ... Must be empty; additional variance estimators are not implemented.
#'
#' @return A spec object (list) used by csdm().
#' @export
csdm_vcov <- function(type = "mg", ...) {
  if (!is.character(type) || length(type) != 1L || is.na(type) || type != "mg") {
    stop("Only type='mg' is implemented.", call. = FALSE)
  }
  if (...length()) {
    stop("Additional variance-covariance options are not implemented.", call. = FALSE)
  }
  spec <- list(type = type, options = list(...))
  class(spec) <- "csdm_vcov_spec"
  spec
}
