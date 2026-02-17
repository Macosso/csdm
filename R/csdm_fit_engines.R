# csdm_fit_engines.R

.csdm_fit_mg <- function(panel_df, formula, id, time, lr = NULL, vcov, ...) {
  ids <- unique(panel_df[[id]])

  econ_names <- .csdm_econ_names(formula, panel_df)

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
    model = "mg",
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
      csa = NULL,
      dropped_units = dropped
    )
  )


  # Compute R2 (mg) from residual matrix and strictly aligned Y
  r2_resid <- .csdm_residual_matrix_r2(
    residuals_e = fit$residuals_e,
    panel_df = panel_df,
    id = id,
    time = time,
    yname = yname
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
  fit$stats$R2_ols_i <- r2_resid$R2_ols_i
  fit$stats$R2_ols_mg <- r2_resid$R2_ols_mg

  fit
}


.csdm_fit_cce <- function(panel_df, formula, id, time, csa, lr = NULL, vcov, ...) {
  if (length(csa$lags) == 1L && as.integer(csa$lags) != 0L) {
    stop("For model='cce', csa$lags must be 0")
  }
  if (length(csa$lags) > 1L && any(as.integer(csa$lags) != 0L)) {
    stop("For model='cce', all csa$lags entries must be 0")
  }

  if (identical(csa$vars, "_none")) {
    csa_vars <- character(0)
  } else if (identical(csa$vars, "_all")) {
    csa_vars <- setdiff(unique(all.vars(formula)), ".csdm_trend__")
  } else {
    csa_vars <- as.character(csa$vars)
  }

  if (length(csa_vars)) {
    csa_attached <- cross_sectional_avg(
      data = panel_df,
      id = id,
      time = time,
      vars = csa_vars,
      leave_out = FALSE,
      suffix = "csa",
      return_mode = "attach",
      na.rm = TRUE
    )
  } else {
    csa_attached <- panel_df
  }
  attr(csa_attached, "csdm_time_levels") <- attr(panel_df, "csdm_time_levels")
  attr(csa_attached, "csdm_id_levels") <- attr(panel_df, "csdm_id_levels")

  econ_names <- .csdm_econ_names(formula, panel_df)

  csa_terms <- if (length(csa_vars)) paste0("csa_", csa_vars) else character(0)
  fml <- if (length(csa_terms)) stats::update(formula, paste0(". ~ . + ", paste(csa_terms, collapse = " + "))) else formula

  ids <- unique(csa_attached[[id]])
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

  for (uid in ids) {
    sub <- csa_attached[csa_attached[[id]] == uid, , drop = FALSE]

    mf <- tryCatch(
      stats::model.frame(stats::update(fml, . ~ . + .csdm_rowid__), sub, na.action = stats::na.omit),
      error = function(e) NULL
    )
    if (is.null(mf) || nrow(mf) == 0L) {
      dropped <- c(dropped, as.character(uid))
      next
    }

    fit <- tryCatch(stats::lm(fml, data = sub, na.action = stats::na.omit), error = function(e) NULL)
    if (is.null(fit)) {
      dropped <- c(dropped, as.character(uid))
      next
    }

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
  E <- .csdm_residual_matrix(csa_attached, id, time, res_long)
  Xb <- .csdm_fitted_matrix(csa_attached, id, time, fit_long)

  yname <- NA_character_
  lhs <- formula[[2]]
  if (is.name(lhs)) yname <- as.character(lhs)

  fit <- list(
    call = NULL,
    formula = NULL,
    model = "cce",
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


  r2_resid <- .csdm_residual_matrix_r2(
    residuals_e = fit$residuals_e,
    panel_df = panel_df,
    id = id,
    time = time,
    yname = yname
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
  fit$stats$R2_ols_i <- r2_resid$R2_ols_i
  fit$stats$R2_ols_mg <- r2_resid$R2_ols_mg

  fit
}


.csdm_fit_dcce <- function(panel_df, formula, id, time, csa, lr, vcov, ...) {
  lr_type <- if (!is.null(lr) && !is.null(lr$type)) as.character(lr$type) else "none"
  lr_ylags <- if (!is.null(lr) && !is.null(lr$ylags)) as.integer(lr$ylags) else 0L
  lr_xdlags <- if (!is.null(lr) && !is.null(lr$xdlags)) as.integer(lr$xdlags) else 0L

  if (!identical(lr_type, "none") && !identical(lr_type, "ardl")) {
    stop("Not implemented yet")
  }

  panel_work <- panel_df
  econ_formula <- formula

  if (identical(lr_type, "ardl") && (isTRUE(lr_ylags > 0L) || isTRUE(lr_xdlags > 0L))) {
    lhs <- formula[[2]]
    if (!is.name(lhs)) {
      stop("For lr(type='ardl') with ylags>0 and/or xdlags>0, the dependent variable in 'formula' must be a simple column name.")
    }
    y_name <- as.character(lhs)
    if (!y_name %in% names(panel_work)) stop("Dependent variable not found in data: ", y_name)

    # Collect all lag term names first, and stop on any collision
    lag_terms_y <- character(0)
    if (isTRUE(lr_ylags > 0L)) {
      lag_terms_y <- paste0("lag", seq_len(lr_ylags), "_", y_name)
    }

    lag_terms_x <- character(0)
    xnames <- character(0)
    if (isTRUE(lr_xdlags > 0L)) {
      tt <- stats::terms(formula)
      rhs_terms <- attr(tt, "term.labels")
      if (length(rhs_terms)) {
        for (term in rhs_terms) {
          is_simple <- grepl("^[.A-Za-z][.A-Za-z0-9._]*$", term) && (term %in% names(panel_work))
          if (!is_simple) {
            stop("xdlags currently supports only simple RHS variable names (no transformations/interactions); offending term: ", term)
          }
        }
        xnames <- rhs_terms
      }

      if (length(xnames)) {
        lag_terms_x <- unlist(
          lapply(xnames, function(xn) paste0("lag", seq_len(lr_xdlags), "_", xn)),
          use.names = FALSE
        )
      }
    }

    lag_terms_all <- c(lag_terms_y, lag_terms_x)
    collisions <- intersect(lag_terms_all, names(panel_work))
    if (length(collisions)) {
      stop("Lag column(s) already exist in data: ", paste(collisions, collapse = ", "))
    }

    # Create columns
    if (length(lag_terms_all)) {
      for (nm in lag_terms_all) panel_work[[nm]] <- NA_real_
    }

    # Fill within-unit y lags
    if (length(lag_terms_y)) {
      ids <- unique(panel_work[[id]])
      for (uid in ids) {
        idx <- which(panel_work[[id]] == uid)
        if (length(idx) == 0L) next
        o <- order(panel_work[[time]][idx])
        idxo <- idx[o]
        yv <- panel_work[[y_name]][idxo]
        for (k in seq_len(lr_ylags)) {
          if (length(yv) <= k) {
            lagv <- rep(NA_real_, length(yv))
          } else {
            lagv <- c(rep(NA_real_, k), yv[seq_len(length(yv) - k)])
          }
          panel_work[[paste0("lag", k, "_", y_name)]][idxo] <- lagv
        }
      }
    }

    # Fill within-unit x distributed lags
    if (length(xnames) && isTRUE(lr_xdlags > 0L)) {
      ids <- unique(panel_work[[id]])
      for (uid in ids) {
        idx <- which(panel_work[[id]] == uid)
        if (length(idx) == 0L) next
        o <- order(panel_work[[time]][idx])
        idxo <- idx[o]
        for (xn in xnames) {
          xv <- panel_work[[xn]][idxo]
          for (k in seq_len(lr_xdlags)) {
            if (length(xv) <= k) {
              lagv <- rep(NA_real_, length(xv))
            } else {
              lagv <- c(rep(NA_real_, k), xv[seq_len(length(xv) - k)])
            }
            panel_work[[paste0("lag", k, "_", xn)]][idxo] <- lagv
          }
        }
      }
    }

    if (length(lag_terms_all)) {
      econ_formula <- stats::update(formula, paste0(". ~ . + ", paste(lag_terms_all, collapse = " + ")))
    }
  }

  if (identical(csa$vars, "_none")) {
    csa_vars <- character(0)
  } else if (identical(csa$vars, "_all")) {
    csa_vars <- setdiff(unique(all.vars(formula)), ".csdm_trend__")
  } else {
    csa_vars <- as.character(csa$vars)
  }

  if (length(csa_vars)) {
    csa_time <- cross_sectional_avg(
      data = panel_work,
      id = id,
      time = time,
      vars = csa_vars,
      leave_out = FALSE,
      suffix = "csa",
      return_mode = "time",
      na.rm = TRUE
    )
  } else {
    csa_time <- unique(panel_work[, time, drop = FALSE])
  }

  o <- order(csa_time[[time]])
  csa_time <- csa_time[o, , drop = FALSE]

  lag_spec <- csa$lags
  if (length(csa_vars) && length(lag_spec) > 1L) {
    lag_spec <- lag_spec[csa_vars]
    lag_spec[is.na(lag_spec)] <- 0L
  }

  add_lag <- function(x, L) {
    if (L <= 0L) return(NULL)
    out <- vector("list", L)
    for (l in seq_len(L)) {
      if (length(x) <= l) {
        out[[l]] <- rep(NA_real_, length(x))
      } else {
        out[[l]] <- c(rep(NA_real_, l), x[seq_len(length(x) - l)])
      }
    }
    out
  }

  csa_term_names <- character(0)
  if (length(csa_vars)) {
    for (v in csa_vars) {
      base_col <- paste0("csa_", v)
      if (!base_col %in% names(csa_time)) next

      maxL <- if (length(lag_spec) == 1L) as.integer(lag_spec) else as.integer(lag_spec[[v]])
      maxL <- ifelse(is.na(maxL), 0L, maxL)
      if (maxL < 0L) stop("csa$lags must be >= 0")

      csa_term_names <- c(csa_term_names, base_col)

      if (maxL > 0L) {
        lag_list <- add_lag(csa_time[[base_col]], maxL)
        for (l in seq_len(maxL)) {
          nm <- paste0(base_col, "_lag", l)
          csa_time[[nm]] <- lag_list[[l]]
          csa_term_names <- c(csa_term_names, nm)
        }
      }
    }
  }

  key <- match(panel_work[[time]], csa_time[[time]])
  csa_attached <- panel_work
  if (nrow(csa_time) && length(setdiff(names(csa_time), time))) {
    for (nm in setdiff(names(csa_time), time)) {
      csa_attached[[nm]] <- csa_time[[nm]][key]
    }
  }

  econ_names <- .csdm_econ_names(econ_formula, panel_work)

  fml <- if (length(csa_term_names)) {
    stats::update(econ_formula, paste0(". ~ . + ", paste(unique(csa_term_names), collapse = " + ")))
  } else {
    econ_formula
  }

  ids <- unique(csa_attached[[id]])
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

  for (uid in ids) {
    sub <- csa_attached[csa_attached[[id]] == uid, , drop = FALSE]

    mf <- tryCatch(
      stats::model.frame(stats::update(fml, . ~ . + .csdm_rowid__), sub, na.action = stats::na.omit),
      error = function(e) NULL
    )
    if (is.null(mf) || nrow(mf) == 0L) {
      dropped <- c(dropped, as.character(uid))
      next
    }

    fit <- tryCatch(stats::lm(fml, data = sub, na.action = stats::na.omit), error = function(e) NULL)
    if (is.null(fit)) {
      dropped <- c(dropped, as.character(uid))
      next
    }

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
  E <- .csdm_residual_matrix(csa_attached, id, time, res_long)
  Xb <- .csdm_fitted_matrix(csa_attached, id, time, fit_long)

  yname <- NA_character_
  lhs <- formula[[2]]
  if (is.name(lhs)) yname <- as.character(lhs)

  fit <- list(
    call = NULL,
    formula = NULL,
    model = "dcce",
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


  r2_resid <- .csdm_residual_matrix_r2(
    residuals_e = fit$residuals_e,
    panel_df = panel_df,
    id = id,
    time = time,
    yname = yname
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
  fit$stats$R2_ols_i <- r2_resid$R2_ols_i
  fit$stats$R2_ols_mg <- r2_resid$R2_ols_mg

  fit
}


.csdm_fit_cs_ardl <- function(panel_df, formula, id, time, csa, lr, vcov, ...) {
  lr_type <- if (!is.null(lr) && !is.null(lr$type)) as.character(lr$type) else "none"
  lr_ylags <- if (!is.null(lr) && !is.null(lr$ylags)) as.integer(lr$ylags) else 0L
  lr_xdlags <- if (!is.null(lr) && !is.null(lr$xdlags)) as.integer(lr$xdlags) else 0L

  if (!identical(lr_type, "ardl")) {
    stop("model='cs_ardl' currently requires lr = csdm_lr(type='ardl', ...) ")
  }
  if (!isTRUE(lr_ylags >= 1L)) {
    stop("model='cs_ardl' requires lr(type='ardl', ylags >= 1) to compute adjustment and long-run effects")
  }

  fit <- .csdm_fit_dcce(
    panel_df = panel_df,
    formula = formula,
    id = id,
    time = time,
    csa = csa,
    lr = lr,
    vcov = vcov,
    ...
  )
  fit$model <- "cs_ardl"

  lhs <- formula[[2]]
  if (!is.name(lhs)) {
    stop("model='cs_ardl' requires a simple dependent variable name in 'formula'")
  }
  yname <- as.character(lhs)

  tt <- stats::terms(formula)
  rhs_terms <- attr(tt, "term.labels")
  xnames <- character(0)
  if (length(rhs_terms)) {
    for (term in rhs_terms) {
      is_simple <- grepl("^[.A-Za-z][.A-Za-z0-9._]*$", term) && (term %in% names(panel_df))
      if (!is_simple) {
        stop(
          "xdlags currently supports only simple RHS variable names (no transformations/interactions); offending term: ",
          term
        )
      }
    }
    xnames <- rhs_terms
  }

  coef_i <- as.matrix(fit$coef_i)
  unit_ids <- rownames(coef_i)

  alpha_terms <- paste0("lag", seq_len(lr_ylags), "_", yname)
  alpha_mat <- coef_i[, intersect(alpha_terms, colnames(coef_i)), drop = FALSE]
  alpha_sum <- if (ncol(alpha_mat)) rowSums(alpha_mat) else rep(0, nrow(coef_i))
  denom_i <- 1 - alpha_sum
  adj_i <- -denom_i

  eps <- 1e-8
  ok_denom <- is.finite(denom_i) & (abs(denom_i) > eps)

  lr_i <- matrix(NA_real_, nrow = nrow(coef_i), ncol = length(xnames))
  colnames(lr_i) <- paste0("lr_", xnames)
  rownames(lr_i) <- unit_ids

  if (length(xnames)) {
    for (j in seq_along(xnames)) {
      xn <- xnames[[j]]
      beta_terms <- c(xn, paste0("lag", seq_len(lr_xdlags), "_", xn))
      beta_mat <- coef_i[, intersect(beta_terms, colnames(coef_i)), drop = FALSE]
      beta_sum <- if (ncol(beta_mat)) rowSums(beta_mat) else rep(NA_real_, nrow(coef_i))
      lr_i[, j] <- beta_sum / denom_i
    }
  }

  if (any(!ok_denom)) {
    warning("cs_ardl: excluding ", sum(!ok_denom), " unit(s) from adjustment/long-run reporting due to near-zero or missing (1 - sum(alpha_y_lags))")
  }

  used_adj <- ok_denom & is.finite(adj_i)
  n_used_adj <- sum(used_adj)
  adj_est <- if (n_used_adj) mean(adj_i[used_adj]) else NA_real_
  adj_se <- if (n_used_adj >= 2L) stats::sd(adj_i[used_adj]) / sqrt(n_used_adj) else NA_real_

  lr_terms <- paste0("lr_", xnames)
  lr_est <- rep(NA_real_, length(xnames))
  lr_se <- rep(NA_real_, length(xnames))
  lr_n  <- rep(0L, length(xnames))
  if (length(xnames)) {
    for (j in seq_along(xnames)) {
      v <- lr_i[, j]
      used <- ok_denom & is.finite(v)
      lr_n[[j]] <- sum(used)
      lr_est[[j]] <- if (lr_n[[j]]) mean(v[used]) else NA_real_
      lr_se[[j]] <- if (lr_n[[j]] >= 2L) stats::sd(v[used]) / sqrt(lr_n[[j]]) else NA_real_
    }
  }

  fit$cs_ardl <- list(
    y = yname,
    x = xnames,
    ylags = lr_ylags,
    xdlags = lr_xdlags,
    unit_ids = unit_ids,
    denom_i = stats::setNames(as.numeric(denom_i), unit_ids),
    adj_i = stats::setNames(as.numeric(adj_i), unit_ids),
    lr_i = lr_i,
    mg = list(
      adj = c(estimate = adj_est, se = adj_se),
      lr = data.frame(
        term = lr_terms,
        estimate = lr_est,
        se = lr_se,
        n_used = lr_n,
        stringsAsFactors = FALSE
      )
    )
  )

  return(fit)
}
