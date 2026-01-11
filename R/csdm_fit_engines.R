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

  list(
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

  list(
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
}


.csdm_fit_dcce <- function(panel_df, formula, id, time, csa, lr, vcov, ...) {
  lr_type <- if (!is.null(lr) && !is.null(lr$type)) as.character(lr$type) else "none"
  lr_ylags <- if (!is.null(lr) && !is.null(lr$ylags)) as.integer(lr$ylags) else 0L

  if (!identical(lr_type, "none") && !identical(lr_type, "ardl")) {
    stop("Not implemented yet")
  }

  panel_work <- panel_df
  econ_formula <- formula

  if (identical(lr_type, "ardl") && isTRUE(lr_ylags > 0L)) {
    lhs <- formula[[2]]
    if (!is.name(lhs)) {
      stop("For lr(type='ardl', ylags>0), the dependent variable in 'formula' must be a simple column name.")
    }
    y_name <- as.character(lhs)
    if (!y_name %in% names(panel_work)) stop("Dependent variable not found in data: ", y_name)

    lag_terms <- character(lr_ylags)
    for (k in seq_len(lr_ylags)) {
      lag_terms[[k]] <- paste0("lag", k, "_", y_name)
      if (!lag_terms[[k]] %in% names(panel_work)) panel_work[[lag_terms[[k]]]] <- NA_real_
    }

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

    econ_formula <- stats::update(formula, paste0(". ~ . + ", paste(lag_terms, collapse = " + ")))
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

  list(
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
}
