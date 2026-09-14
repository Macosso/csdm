.csdm_fit_units <- function(panel_df, formula, id, time, econ_formula = formula, model = "mg", csa = NULL) {
  design <- .csdm_design(formula, panel_df)
  economic <- .csdm_design(econ_formula, panel_df)
  X <- design$X
  y <- design$y
  econ_names <- colnames(economic$X)
  ids <- unique(panel_df[[id]])
  B <- matrix(NA_real_, length(ids), length(econ_names), dimnames = list(ids, econ_names))
  e <- xb <- rep(NA_real_, nrow(panel_df))
  units <- vector("list", length(ids)); names(units) <- ids
  complete <- is.finite(y) & rowSums(!is.finite(X)) == 0
  na_fun <- attr(panel_df, "csdm_na_action")
  if (is.null(na_fun)) na_fun <- stats::na.omit
  if (!all(complete) && identical(na_fun, stats::na.fail)) stop("Missing values in the estimation design.")
  if (!any(complete)) stop("No complete observations in the estimation design.")
  df_e <- stats::setNames(rep(NA_real_, length(ids)), ids)
  for (uid in ids) {
    rows <- which(panel_df[[id]] == uid & complete)
    if (!length(rows)) {
      units[[uid]] <- list(rows = integer(), nobs = 0L, rank = 0L, df.residual = 0L,
        reason = "no complete observations", aliased = econ_names)
      next
    }
    fit <- .csdm_unit_regression(X[rows, , drop = FALSE], y[rows], econ_names)
    if (!is.na(fit$reason)) {
      units[[uid]] <- list(rows = panel_df$.csdm_rowid__[rows], nobs = length(rows),
        rank = fit$rank, df.residual = fit$df.residual, reason = fit$reason, aliased = fit$aliased)
      next
    }
    B[uid, ] <- fit$coefficients[econ_names]
    e[rows] <- fit$residuals
    xb[rows] <- fit$fitted.values
    df_e[uid] <- fit$df.residual
    units[[uid]] <- list(rows = panel_df$.csdm_rowid__[rows], nobs = length(rows),
      rank = fit$rank, df.residual = fit$df.residual, coefficients = fit$coefficients,
      reason = NA_character_, aliased = character())
  }
  dropped <- names(units)[vapply(units, function(u) !is.na(u$reason), logical(1))]
  if (length(dropped)) warning("Excluded unit(s): ", paste(dropped, collapse = ", "),
    ". See fit$units for identification/sample reasons.", call. = FALSE)
  if (length(dropped) == length(ids)) stop("No identified units with positive residual degrees of freedom.")
  res <- panel_df[c(id, time)]; res$residual <- e
  fv <- panel_df[c(id, time)]; fv$xb <- xb
  E <- .csdm_residual_matrix(panel_df, id, time, res)
  fitted <- .csdm_fitted_matrix(panel_df, id, time, fv)
  V <- .csdm_mg_vcov(B)
  panel_df$.csdm_response__ <- y
  stats <- .csdm_compute_fit_stats(panel_df, id, time, ".csdm_response__", E)
  r2 <- .csdm_residual_matrix_r2(E, panel_df, id, time, ".csdm_response__", df_e,
    intercept = attr(economic$terms, "intercept") == 1L)
  stats$R2_i <- r2$R2_i
  stats$R2_mg <- r2$R2_mg
  stats$R2_ols_mg <- r2$R2_ols_mg
  used <- is.finite(e)
  list(model = model, id = id, time = time, coef_mg = colMeans(B, na.rm = TRUE),
    se_mg = sqrt(diag(V)), vcov_mg = V, coef_i = B, residuals_e = E, fitted_xb = fitted,
    stats = stats, units = units, model_frame = economic$frame[used, , drop = FALSE],
    model_matrix = economic$X[used, , drop = FALSE], augmented_matrix = X[used, , drop = FALSE],
    terms = economic$terms, fitted_formula = econ_formula,
    contrasts = attr(economic$X, "contrasts"),
    sample = data.frame(row = panel_df$.csdm_rowid__, used = used, residual = e, fitted = xb),
    meta = list(N = length(ids), T = ncol(E), csa = csa,
      dropped_units = dropped))
}
