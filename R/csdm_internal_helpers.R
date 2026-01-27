# R-style significance codes for p-values (see ?summary.lm)
.csdm_significance_stars <- function(p) {
  cut(p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", ".", ""),
      right = FALSE)
}

# csdm_internal_helpers.R

#' Compute R-squared (mg) from the final wide residual matrix and strictly aligned Y
#'
#' This function computes R² (mg) for mean-group models using only the cells present in the final
#' wide residual matrix (residuals_e), with Y strictly aligned by id/time. This ensures R² reflects
#' the model's actual predictions and is robust to data alignment and NA-handling issues. Optionally
#' also returns the per-unit OLS R² for debugging.
#'
#' @param residuals_e Wide matrix of residuals (units x time), as from .csdm_residual_matrix().
#' @param panel_df Original panel data.frame.
#' @param id Name of id column.
#' @param time Name of time column.
#' @param yname Name of outcome variable.
#' @return List with R2_i (by unit), R2_mg (mean), Y (aligned), and optionally R2_ols_i, R2_ols_mg.
.csdm_residual_matrix_r2 <- function(residuals_e, panel_df, id, time, yname) {
  # Strictly align Y to residuals_e: same dimnames, types, and only cells present in E
  E <- residuals_e
  ids_levels <- rownames(E)
  time_levels <- colnames(E)
  Y <- matrix(NA_real_, nrow = length(ids_levels), ncol = length(time_levels),
              dimnames = list(ids_levels, time_levels))
  ii <- match(as.character(panel_df[[id]]), ids_levels)
  tt <- match(as.character(panel_df[[time]]), time_levels)
  keep <- is.finite(ii) & is.finite(tt)
  if (any(keep)) {
    Y[cbind(ii[keep], tt[keep])] <- as.numeric(panel_df[[yname]][keep])
  }

  # Compute R2 for each unit using only cells present in E and Y
  R2_i <- rep(NA_real_, nrow(E)); names(R2_i) <- ids_levels
  for (r in seq_len(nrow(E))) {
    er <- E[r, ]; yr <- Y[r, ]
    ok <- is.finite(er) & is.finite(yr)
    if (sum(ok) < 2L) next
    sse <- sum((er[ok])^2)
    yc <- yr[ok]
    sst <- sum((yc - mean(yc))^2)
    if (!is.finite(sst) || sst <= 0) next
    R2_i[[r]] <- 1 - sse / sst
  }
  R2_mg <- if (all(is.na(R2_i))) NA_real_ else mean(R2_i, na.rm = TRUE)

  # Optionally: also compute OLS-based R2 for each unit (using lm on regression sample)
  R2_ols_i <- rep(NA_real_, nrow(E)); names(R2_ols_i) <- ids_levels
  for (r in seq_len(nrow(E))) {
    uid <- ids_levels[[r]]
    sub <- panel_df[as.character(panel_df[[id]]) == uid, , drop = FALSE]
    if (!nrow(sub)) next
    if (!all(c(yname) %in% names(sub))) next
    # Try to fit OLS on available data for this unit
    yvals <- sub[[yname]]
    oky <- is.finite(yvals)
    if (sum(oky) < 2L) next
    yc <- yvals[oky]
    sst <- sum((yc - mean(yc))^2)
    if (!is.finite(sst) || sst <= 0) next
    # OLS R2: need to know the regression sample (not just y)
    # For debugging, use all available y for this unit
    # (If you want exact OLS R2, use the regression sample as in the fit loop)
    # Here, just for reference
    # R2_ols_i[[r]] <- ...
    # For now, leave as NA or implement if needed
  }
  R2_ols_mg <- if (all(is.na(R2_ols_i))) NA_real_ else mean(R2_ols_i, na.rm = TRUE)

  list(R2_i = R2_i, R2_mg = R2_mg, Y = Y, R2_ols_i = R2_ols_i, R2_ols_mg = R2_ols_mg)
}

.csdm_mg_vcov <- function(B) {
  B <- as.matrix(B)
  ok <- stats::complete.cases(B)
  B <- B[ok, , drop = FALSE]
  N <- nrow(B)
  if (N < 2L) {
    V <- matrix(NA_real_, ncol(B), ncol(B))
    colnames(V) <- rownames(V) <- colnames(B)
    return(V)
  }
  bbar <- colMeans(B)
  dev <- sweep(B, 2L, bbar, "-")
  V <- crossprod(dev) / (N * (N - 1))
  dimnames(V) <- list(colnames(B), colnames(B))
  V
}


.csdm_econ_names <- function(formula, panel_df, na.action = stats::na.omit) {
  mf <- stats::model.frame(formula, panel_df, na.action = na.action)
  X <- stats::model.matrix(formula, mf)
  colnames(X)
}


.csdm_residual_matrix <- function(panel_df, id, time, res_long) {
  .or <- function(x, y) if (is.null(x)) y else x

  ids <- .or(attr(panel_df, "csdm_id_levels"), sort(unique(panel_df[[id]])))
  times <- .or(attr(panel_df, "csdm_time_levels"), sort(unique(panel_df[[time]])))

  E <- matrix(NA_real_, nrow = length(ids), ncol = length(times))
  rownames(E) <- as.character(ids)
  colnames(E) <- as.character(times)

  if (nrow(res_long) == 0L) return(E)

  ii <- match(res_long[[id]], ids)
  tt <- match(res_long[[time]], times)
  keep <- is.finite(ii) & is.finite(tt)
  if (any(keep)) E[cbind(ii[keep], tt[keep])] <- res_long$residual[keep]
  E
}


.csdm_fitted_matrix <- function(panel_df, id, time, fit_long) {
  .or <- function(x, y) if (is.null(x)) y else x

  ids <- .or(attr(panel_df, "csdm_id_levels"), sort(unique(panel_df[[id]])))
  times <- .or(attr(panel_df, "csdm_time_levels"), sort(unique(panel_df[[time]])))

  Xb <- matrix(NA_real_, nrow = length(ids), ncol = length(times))
  rownames(Xb) <- as.character(ids)
  colnames(Xb) <- as.character(times)

  if (nrow(fit_long) == 0L) return(Xb)

  ii <- match(fit_long[[id]], ids)
  tt <- match(fit_long[[time]], times)
  keep <- is.finite(ii) & is.finite(tt)
  if (any(keep)) Xb[cbind(ii[keep], tt[keep])] <- fit_long$xb[keep]
  Xb
}


.csdm_make_coef_table <- function(estimate, se, n_used = NULL) {
  estimate <- as.numeric(estimate)
  se <- as.numeric(se)
  z <- estimate / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  crit <- stats::qnorm(0.975)
  ci_lo <- estimate - crit * se
  ci_hi <- estimate + crit * se

  bad <- !is.finite(se) | se <= 0
  if (any(bad)) {
    z[bad] <- NA_real_
    p[bad] <- NA_real_
    ci_lo[bad] <- NA_real_
    ci_hi[bad] <- NA_real_
  }

  # Add R-style significance codes as in summary.lm
  tab <- data.frame(
    `Coef.` = estimate,
    `Std. Err.` = se,
    z = z,
    `P>|z|` = p,
    check.names = FALSE
  )
  tab$Signif. <- .csdm_significance_stars(tab$`P>|z|`)
  tab$`CI 2.5%` <- ci_lo
  tab$`CI 97.5%` <- ci_hi
  # Reorder columns to put Signif. after P>|z|
  tab <- tab[, c("Coef.", "Std. Err.", "z", "P>|z|", "Signif.", "CI 2.5%", "CI 97.5%")]
  if (!is.null(n_used)) {
    tab$n_used <- as.integer(n_used)
  }
  tab
}


.csdm_format_csa_footer <- function(csa_spec) {
  if (is.null(csa_spec) || identical(csa_spec$vars, "_none")) {
    return(list(csa_vars = "none", csa_lags = "0"))
  }

  vars_txt <- if (identical(csa_spec$vars, "_all")) {
    "_all"
  } else {
    paste(as.character(csa_spec$vars), collapse = ", ")
  }

  lags_txt <- if (is.null(csa_spec$lags)) {
    "0"
  } else if (length(csa_spec$lags) == 1L) {
    as.character(as.integer(csa_spec$lags))
  } else {
    "(per-variable)"
  }

  list(csa_vars = vars_txt, csa_lags = lags_txt)
}



# Compute fit statistics for csdm models, including robust R² (mg)
#
# R² (mg) is computed for the *model-predicted values on exactly the sample present in the final residual matrix*,
# ensuring robustness to index/panel alignment and NA-handling issues. This matches the definition used in xtdcce2.
# See .csdm_residual_matrix_r2 for the main logic.
.csdm_compute_fit_stats <- function(panel_df, id, time, yname, residuals_e, cd_min_overlap = 2L) {
  E <- residuals_e
  nobs <- if (is.matrix(E)) sum(is.finite(E)) else NA_integer_

  # R2 (mg) and R2_i are now computed in .csdm_residual_matrix_r2 and attached by the fit engine.
  # This fallback block is retained for legacy or edge cases only.
  R2_i <- NULL
  R2_mg <- NA_real_
  if (is.matrix(E) && is.character(yname) && length(yname) == 1L && yname %in% names(panel_df)) {
    ids_levels <- rownames(E)
    time_levels <- colnames(E)
    Y <- matrix(NA_real_, nrow = length(ids_levels), ncol = length(time_levels),
                dimnames = list(ids_levels, time_levels))
    ii <- match(as.character(panel_df[[id]]), ids_levels)
    tt <- match(as.character(panel_df[[time]]), time_levels)
    keep <- is.finite(ii) & is.finite(tt)
    if (any(keep)) {
      Y[cbind(ii[keep], tt[keep])] <- as.numeric(panel_df[[yname]][keep])
    }
    R2_i <- rep(NA_real_, nrow(E)); names(R2_i) <- ids_levels
    for (r in seq_len(nrow(E))) {
      er <- E[r, ]; yr <- Y[r, ]
      ok <- is.finite(er) & is.finite(yr)
      if (sum(ok) < 2L) next
      sse <- sum((er[ok])^2)
      yc <- yr[ok]
      sst <- sum((yc - mean(yc))^2)
      if (!is.finite(sst) || sst <= 0) next
      R2_i[[r]] <- 1 - sse / sst
    }
    R2_mg <- if (all(is.na(R2_i))) NA_real_ else mean(R2_i, na.rm = TRUE)
  }


  # CD test (unchanged)
  cd_stat <- NA_real_
  cd_p_value <- NA_real_
  if (is.matrix(E) && nrow(E) >= 2L && ncol(E) >= 2L) {
    cd <- tryCatch(
      cd_test(E, min_overlap = as.integer(cd_min_overlap)),
      error = function(e) NULL
    )
    if (!is.null(cd)) {
      cd_stat <- as.numeric(cd$statistic)
      cd_p_value <- as.numeric(cd$p.value)
    }
  }

  list(
    nobs = as.integer(nobs),
    R2_i = R2_i,
    R2_mg = as.numeric(R2_mg),
    cd_stat = as.numeric(cd_stat),
    cd_p_value = as.numeric(cd_p_value)
  )
}
