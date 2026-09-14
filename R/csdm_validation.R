.csdm_integer <- function(x, name, scalar = TRUE) {
  if (!is.numeric(x) || !length(x) || (scalar && length(x) != 1L) ||
      any(!is.finite(x)) || any(x < 0 | x != floor(x) | x > .Machine$integer.max)) {
    stop("'", name, "' must contain nonnegative finite integers.", call. = FALSE)
  }
  stats::setNames(as.integer(x), names(x))
}

.csdm_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.csdm_validate_specs <- function(model, csa, lr, pooled, vcov, fullsample, mgmissing) {
  specs <- list(csa = csa, lr = lr, pooled = pooled, vcov = vcov)
  for (nm in names(specs)) {
    if (!inherits(specs[[nm]], paste0("csdm_", nm, "_spec"))) {
      stop("Use csdm_", nm, "() to construct '", nm, "'.", call. = FALSE)
    }
    expected <- names(do.call(get(paste0("csdm_", nm)), list()))
    if (!identical(sort(names(specs[[nm]])), sort(expected))) {
      stop("Invalid '", nm, "' specification fields.", call. = FALSE)
    }
    do.call(get(paste0("csdm_", nm)), if (nm == "vcov") list(type = vcov$type) else specs[[nm]])
  }
  if (.csdm_flag(fullsample, "fullsample") || .csdm_flag(mgmissing, "mgmissing")) {
    stop("fullsample=TRUE and mgmissing=TRUE are not implemented.", call. = FALSE)
  }
  if (!identical(csa$scope, "estimation") || !is.null(csa$cluster)) {
    stop("Only csa scope='estimation' without cluster is implemented.", call. = FALSE)
  }
  if (!is.null(pooled$vars) || pooled$constant || pooled$trend) {
    stop("Pooled restrictions are not implemented.", call. = FALSE)
  }
  if (!identical(vcov$type, "mg") || length(vcov$options)) {
    stop("Only vcov=csdm_vcov('mg') without options is implemented.", call. = FALSE)
  }
  if (!is.null(lr$vars) || length(lr$options)) stop("lr vars/options are not implemented.", call. = FALSE)
  if (model %in% c("mg", "cce") &&
      (lr$type != "none" || lr$ylags != 0L || lr$xdlags != 0L)) {
    stop("Use model='dcce' or 'cs_ardl' for an ARDL specification.", call. = FALSE)
  }
  if (lr$type == "none" && (lr$ylags != 0L || lr$xdlags != 0L)) {
    stop("Positive model lags require lr type='ardl'.", call. = FALSE)
  }
}

.csdm_csa_lags <- function(lags, vars) {
  if (is.null(names(lags))) return(stats::setNames(rep(lags, length(vars)), vars))
  unknown <- setdiff(names(lags), vars)
  if (length(unknown)) stop("Unknown CSA lag variable(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  out <- stats::setNames(integer(length(vars)), vars)
  out[names(lags)] <- lags
  out
}
