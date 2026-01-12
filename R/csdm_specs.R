# csdm_specs.R

#' Specification: Cross-sectional averages (CSA)
#'
#' @param vars Character. One of "_all", "_none", or a character vector of variable names.
#' @param lags Integer. Either a scalar integer >= 0 applied to all CSA variables,
#'   or a named integer vector giving per-variable maximum lags.
#' @param scope Character vector. One or more of c("estimation","global","cluster").
#' @param cluster Reserved for future use.
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
  scope = c("estimation", "global", "cluster"),
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
    if (any(!nzchar(vars))) stop("'vars' contains empty strings.")
    vars <- unique(as.character(vars))
  }

  # lags
  if (length(lags) == 1L) {
    if (!is.numeric(lags) || !is.finite(lags)) stop("'lags' must be a finite integer >= 0.")
    lags <- as.integer(lags)
    if (lags < 0L) stop("'lags' must be >= 0.")
  } else {
    if (is.null(names(lags)) || any(!nzchar(names(lags)))) {
      stop("If 'lags' is not scalar, it must be a *named* integer vector.")
    }
    if (!is.numeric(lags) || any(!is.finite(lags))) stop("'lags' must be finite integers.")
    lags <- as.integer(lags)
    if (any(lags < 0L)) stop("All entries of 'lags' must be >= 0.")
    lags <- lags[unique(names(lags))]
  }

  # scope
  allowed <- c("estimation", "global", "cluster")
  if (!is.character(scope) || length(scope) == 0L) stop("'scope' must be a character vector.")
  if (any(!scope %in% allowed)) {
    bad <- setdiff(unique(scope), allowed)
    stop("Invalid 'scope': ", paste(bad, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "))
  }
  scope <- unique(scope)

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
#' @param vars Reserved for future use.
#' @param type One of c("none","ecm","ardl","csdl").
#' @param ylags Integer >= 0. Within-unit lags of the dependent variable to include
#'   when supported by the chosen model/type.
#' @param xdlags Integer >= 0. Scalar distributed lags to apply to each RHS regressor
#'   when supported by the chosen model/type.
#' @param options Reserved for future use.
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
                    type = c("none", "ecm", "ardl", "csdl"),
                    ylags = 0,
                    xdlags = 0,
                    options = list()) {
  type <- match.arg(type)

  if (!is.numeric(ylags) || length(ylags) != 1L || !is.finite(ylags)) {
    stop("'ylags' must be a finite integer >= 0.")
  }
  ylags <- as.integer(ylags)
  if (ylags < 0L) stop("'ylags' must be >= 0.")

  if (!is.numeric(xdlags) || length(xdlags) != 1L || !is.finite(xdlags)) {
    stop("'xdlags' must be a finite integer >= 0.")
  }
  xdlags <- as.integer(xdlags)
  if (xdlags < 0L) stop("'xdlags' must be >= 0.")

  spec <- list(vars = vars, type = type, ylags = ylags, xdlags = xdlags, options = options)
  class(spec) <- "csdm_lr_spec"
  spec
}


#' Specification: Pooled constraints (stub)
#'
#' @param vars Reserved for future use.
#' @param constant Logical; pooled constant.
#' @param trend Logical; pooled trend.
#'
#' @return A spec object (list) used by csdm().
#' @export
csdm_pooled <- function(vars = NULL, constant = FALSE, trend = FALSE) {
  spec <- list(vars = vars, constant = isTRUE(constant), trend = isTRUE(trend))
  class(spec) <- "csdm_pooled_spec"
  spec
}


#' Specification: Variance-covariance for MG output (stub)
#'
#' @param type One of c("mg","np","nw","wpn","ols").
#' @param ... Reserved for future use.
#'
#' @return A spec object (list) used by csdm().
#' @export
csdm_vcov <- function(type = c("mg", "np", "nw", "wpn", "ols"), ...) {
  type <- match.arg(type)
  spec <- list(type = type, options = list(...))
  class(spec) <- "csdm_vcov_spec"
  spec
}
