# R-style significance codes for p-values (see ?summary.lm)
.csdm_significance_stars <- function(p) {
  cut(p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", ".", ""),
      right = FALSE)
}

# csdm_internal_helpers.R

# Compute R-squared (mg) from the final wide residual matrix and strictly aligned Y
#
# This function computes R² (mg) for mean-group models using only the cells present in the final
# wide residual matrix (residuals_e), with Y strictly aligned by id/time. This ensures R² reflects
# the model's actual predictions and is robust to data alignment and NA-handling issues. Optionally
# also returns the per-unit OLS R² for debugging.
.csdm_residual_matrix_r2 <- function(residuals_e, panel_df, id, time, yname, df_e = NULL) {

  E <- residuals_e
  ids_levels  <- rownames(E)
  time_levels <- colnames(E)

  # --- Align Y to E (IMPORTANT: yname must be the ACTUAL LHS used in estimation, e.g. dy for CS-ARDL/ECM) ---
  Y <- matrix(NA_real_, nrow = length(ids_levels), ncol = length(time_levels),
              dimnames = list(ids_levels, time_levels))

  ii <- match(as.character(panel_df[[id]]), ids_levels)
  tt <- match(as.character(panel_df[[time]]), time_levels)
  keep <- is.finite(ii) & is.finite(tt)

  if (any(keep)) {
    Y[cbind(ii[keep], tt[keep])] <- as.numeric(panel_df[[yname]][keep])
  }

  # --- Per unit statistics ---
  n <- nrow(E)
  R2_i <- rep(NA_real_, n); names(R2_i) <- ids_levels
  Ti   <- rep(NA_integer_, n); names(Ti) <- ids_levels
  SSE  <- rep(NA_real_, n); names(SSE)  <- ids_levels
  SST  <- rep(NA_real_, n); names(SST)  <- ids_levels
  dfy  <- rep(NA_integer_, n); names(dfy) <- ids_levels

  for (r in seq_len(n)) {
    er <- E[r, ]
    yr <- Y[r, ]
    ok <- is.finite(er) & is.finite(yr)
    Ti[r] <- sum(ok)

    if (Ti[r] < 2L) next

    SSE[r] <- sum(er[ok]^2)

    yc <- yr[ok]
    SST[r] <- sum((yc - mean(yc))^2)
    dfy[r] <- Ti[r] - 1L

    if (is.finite(SST[r]) && SST[r] > 0) {
      R2_i[r] <- 1 - SSE[r] / SST[r]   # plain per-unit R²
    }
  }

  # --- Plain MG R² (mean of unit R²) ---
  R2_ols_mg <- if (all(is.na(R2_i))) NA_real_ else mean(R2_i, na.rm = TRUE)

  # --- xtdcce2-style R²(CCEMG): 1 - s_mg^2 / s^2 ---
  # s^2 = pooled variance around unit means using the SAME trimmed sample as E
  ok_s2 <- is.finite(SST) & is.finite(dfy) & dfy > 0
  s2 <- if (!any(ok_s2)) NA_real_ else sum(SST[ok_s2]) / sum(dfy[ok_s2])

  # s_mg^2 = mean of unit sigma^2_i = SSE_i / df_ei
  # df_e must be provided for exact replication (varies with lags, rank, collinearity)
  if (is.null(df_e)) {
    # fallback: use Ti-1 (not exact, but safe)
    df_e <- Ti - 1L
  } else {
    df_e <- as.numeric(df_e)
    if (length(df_e) == 1L) df_e <- rep(df_e, n)
    names(df_e) <- ids_levels
  }

  ok_sig <- is.finite(SSE) & is.finite(df_e) & df_e > 0
  s_mg2 <- if (!any(ok_sig)) NA_real_ else mean(SSE[ok_sig] / df_e[ok_sig])

  R2_mg <- if (!is.finite(s2) || s2 <= 0 || !is.finite(s_mg2)) NA_real_ else 1 - s_mg2 / s2

  list(
    R2_i = R2_i,
    R2_ols_mg = R2_ols_mg,
    R2_mg = R2_mg,
    Ti = Ti,
    SSE = SSE,
    SST = SST,
    s2 = s2,
    s_mg2 = s_mg2
  )
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
  }
  # Classic CD test only (for model summary)
  CD_stat <- NA_real_
  CD_p <- NA_real_
  # Advanced tests (not shown in summary by default)
  CDw_stat <- NA_real_
  CDw_p <- NA_real_
  CDw_plus_stat <- NA_real_
  CDw_plus_p <- NA_real_
  CDstar_stat <- NA_real_
  CDstar_p <- NA_real_
  CDstar_n_pc <- NA_integer_

  if (is.matrix(E) && nrow(E) >= 2L && ncol(E) >= 2L) {
    # Compute classic CD
    cd_classic <- tryCatch(
      cd_test(E, type = "CD"),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(cd_classic) && !is.null(cd_classic$tests$CD)) {
      CD_stat <- as.numeric(cd_classic$tests$CD$statistic)
      CD_p <- as.numeric(cd_classic$tests$CD$p.value)
    }

    # Compute all tests for storage (user can access via cd_test later)
    cd_all <- tryCatch(
      suppressWarnings(cd_test(E, type = "all", n_pc = 4L)),
      error = function(e) NULL
    )
    if (!is.null(cd_all) && !is.null(cd_all$tests)) {
      if (!is.null(cd_all$tests$CDw)) {
        CDw_stat <- as.numeric(cd_all$tests$CDw$statistic)
        CDw_p <- as.numeric(cd_all$tests$CDw$p.value)
      }
      if (!is.null(cd_all$tests$CDw_plus)) {
        CDw_plus_stat <- as.numeric(cd_all$tests$CDw_plus$statistic)
        CDw_plus_p <- as.numeric(cd_all$tests$CDw_plus$p.value)
      }
      if (!is.null(cd_all$tests$CDstar)) {
        CDstar_stat <- as.numeric(cd_all$tests$CDstar$statistic)
        CDstar_p <- as.numeric(cd_all$tests$CDstar$p.value)
        if (!is.null(cd_all$tests$CDstar$n_pc)) {
          CDstar_n_pc <- as.integer(cd_all$tests$CDstar$n_pc)
        }
      }
    }
  }

  list(
    nobs = as.integer(nobs),
    CD_stat = as.numeric(CD_stat),
    CD_p = as.numeric(CD_p),
    CDw_stat = as.numeric(CDw_stat),
    CDw_p = as.numeric(CDw_p),
    CDw_plus_stat = as.numeric(CDw_plus_stat),
    CDw_plus_p = as.numeric(CDw_plus_p),
    CDstar_stat = as.numeric(CDstar_stat),
    CDstar_p = as.numeric(CDstar_p),
    CDstar_n_pc = as.integer(CDstar_n_pc)  )
}
