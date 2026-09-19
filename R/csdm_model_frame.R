.csdm_design <- function(formula, data) {
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || is.matrix(y)) stop("The response must evaluate to a numeric vector.")
  if (length(attr(stats::terms(mf), "offset"))) stop("Offsets are not implemented.")
  X <- stats::model.matrix(stats::terms(mf), mf)
  list(frame = mf, terms = stats::terms(mf), X = X, y = y)
}

.csdm_augment <- function(data, formula, econ_formula, id, time, csa, fullsample = FALSE) {
  base <- .csdm_design(formula, data)
  if (identical(csa$vars, "_none")) return(list(data = data, formula = econ_formula, csa = csa))
  if (identical(csa$vars, "_all")) {
    keep <- !colnames(base$X) %in% c("(Intercept)", ".csdm_trend__")
    values <- cbind(base$y, base$X[, keep, drop = FALSE])
    colnames(values)[1L] <- paste(deparse(formula[[2L]]), collapse = "")
  } else {
    if (!all(csa$vars %in% names(data))) stop("Unknown CSA variable.")
    if (!all(vapply(data[csa$vars], is.numeric, logical(1)))) stop("Explicit CSA variables must be numeric.")
    values <- as.matrix(data[csa$vars])
  }
  labels <- colnames(values)
  if (anyDuplicated(labels)) stop("Default CSA labels are ambiguous; specify CSA variables explicitly.")
  lags <- .csdm_csa_lags(csa$lags, labels)
  # Estimation scope uses complete base model rows before lag trimming.
  eligible <- is.finite(base$y) & rowSums(!is.finite(base$X)) == 0
  times <- sort(unique(data[[time]]))
  terms <- character()
  source_n <- integer(length(labels))
  for (j in seq_along(labels)) {
    source <- if (fullsample) is.finite(values[, j]) else eligible
    source_n[j] <- sum(source)
    means <- vapply(times, function(t) {
      v <- values[source & data[[time]] == t, j]
      v <- v[is.finite(v)]
      if (length(v)) mean(v) else NA_real_
    }, numeric(1))
    for (k in 0:lags[j]) {
      nm <- paste0(".csdm_csa", j, "_lag", k)
      if (nm %in% names(data)) stop("CSA column collision: ", nm)
      v <- if (k == 0) means else .csdm_lag(means, times, k, attr(data, "csdm_time_step"))
      data[[nm]] <- v[match(data[[time]], times)]
      terms <- c(terms, nm)
    }
  }
  csa$resolved_vars <- labels
  csa$resolved_lags <- lags
  source_n <- stats::setNames(source_n, labels)
  csa$source_n <- if (fullsample) source_n else sum(eligible)
  csa$source_n_by_variable <- source_n
  csa$fullsample <- fullsample
  list(data = data, formula = stats::update(econ_formula,
    paste(". ~ . +", paste(terms, collapse = " + "))), csa = csa)
}
