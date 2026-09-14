.csdm_fit_units <- function(panel_df, formula, id, time, econ_formula = formula, model = "mg", csa = NULL) {
  ids <- unique(panel_df[[id]])

  econ_names <- .csdm_econ_names(econ_formula, panel_df)

  coef_i <- matrix(NA_real_, nrow = length(ids), ncol = length(econ_names),
                   dimnames = list(as.character(ids), econ_names))

  res_long <- data.frame(stringsAsFactors = FALSE, check.names = FALSE)
  res_long[[id]] <- character(0)
  res_long[[time]] <- numeric(0)
  res_long$residual <- numeric(0)

  fit_long <- data.frame(stringsAsFactors = FALSE, check.names = FALSE)
  fit_long[[id]] <- character(0)
  fit_long[[time]] <- numeric(0)
  fit_long$xb <- numeric(0)

  dropped <- character(0)
  r2_i <- stats::setNames(rep(NA_real_, length(ids)), as.character(ids))

  # residual degrees of freedom for R2_i calculation
  df_e <- stats::setNames(rep(NA_real_, length(ids)), as.character(ids))

  for (uid in ids) {
    sub <- panel_df[panel_df[[id]] == uid, , drop = FALSE]

    mf <- tryCatch(
      stats::model.frame(stats::update(formula, . ~ . + .csdm_rowid__), sub, na.action = stats::na.omit),
      error = function(e) NULL
    )
    if (is.null(mf) || nrow(mf) == 0L) {
      dropped <- c(dropped, as.character(uid))
      next
    }

    fit <- tryCatch(stats::lm(formula, data = sub, na.action = stats::na.omit), error = function(e) NULL)
    if (is.null(fit)) {
      dropped <- c(dropped, as.character(uid))
      next
    }

    df_e[[as.character(uid)]] <- stats::df.residual(fit)

    # extract residuals
    # Unit-level R2 on the exact estimation sample
    y_used <- tryCatch(stats::model.response(stats::model.frame(fit)), error = function(e) NULL)
    e_used <- tryCatch(as.numeric(fit$residuals), error = function(e) NULL)
    if (!is.null(y_used) && !is.null(e_used) && length(y_used) == length(e_used) && length(y_used) >= 2L) {
      ok <- is.finite(y_used) & is.finite(e_used)
      if (sum(ok) >= 2L) {
        sse <- sum((e_used[ok])^2)
        yc <- y_used[ok]
        sst <- sum((yc - mean(yc))^2)
        if (is.finite(sse) && is.finite(sst) && sst > 0) {
          r2_i[[as.character(uid)]] <- 1 - sse / sst
        }
      }
    }

    cf <- stats::coef(fit)
    keep <- intersect(econ_names, names(cf))
    if (length(keep)) coef_i[as.character(uid), keep] <- cf[keep]

    used_rowid <- mf[[".csdm_rowid__"]]
    if (is.null(used_rowid)) next
    idx_used <- match(used_rowid, sub$.csdm_rowid__)
    if (anyNA(idx_used)) next
    sub_used <- sub[idx_used, , drop = FALSE]
    if (nrow(sub_used) != length(fit$residuals)) next

    chunk <- data.frame(
      sub_used[[id]],
      sub_used[[time]],
      residual = as.numeric(fit$residuals),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    names(chunk)[1:2] <- c(id, time)
    res_long <- rbind(res_long, chunk)

    fchunk <- data.frame(
      sub_used[[id]],
      sub_used[[time]],
      xb = as.numeric(fit$fitted.values),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    names(fchunk)[1:2] <- c(id, time)
    fit_long <- rbind(fit_long, fchunk)
  }

  coef_mg <- colMeans(coef_i, na.rm = TRUE)
  V <- .csdm_mg_vcov(coef_i)
  se <- sqrt(diag(V))

  E <- .csdm_residual_matrix(panel_df, id, time, res_long)
  Xb <- .csdm_fitted_matrix(panel_df, id, time, fit_long)

  yname <- NA_character_
  lhs <- formula[[2]]
  if (is.name(lhs)) yname <- as.character(lhs)

  fit <- list(
    call = NULL,
    formula = NULL,
    model = model,
    id = id,
    time = time,
    coef_mg = coef_mg,
    se_mg = se,
    vcov_mg = V,
    coef_i = coef_i,
    residuals_e = E,
    fitted_xb = Xb,
    meta = list(
      N = length(unique(panel_df[[id]])),
      T = length(unique(panel_df[[time]])),
      csa = csa,
      dropped_units = dropped
    )
  )


  # Compute R2 (mg) from residual matrix and strictly aligned Y
  r2_resid <- .csdm_residual_matrix_r2(
    residuals_e = fit$residuals_e,
    panel_df = panel_df,
    id = id,
    time = time,
    yname = yname,
    df_e = df_e
  )
  fit$stats <- .csdm_compute_fit_stats(
    panel_df = panel_df,
    id = id,
    time = time,
    yname = yname,
    residuals_e = fit$residuals_e
  )
  fit$stats$R2_i <- r2_resid$R2_i
  fit$stats$R2_mg <- r2_resid$R2_mg
  fit$stats$R2_ols_mg <- r2_resid$R2_ols_mg

  fit
}
