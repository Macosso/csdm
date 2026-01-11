# csdm.R

#' Unified front door for csdm estimators
#'
#' @export
csdm <- function(
  formula, data, id, time,
  model = c("mg", "cce", "dcce", "cs_ardl", "cs_ecm", "cs_dl"),
  csa = csdm_csa(),
  lr = csdm_lr(),
  pooled = csdm_pooled(),
  trend = c("none", "unit", "pooled"),
  fullsample = FALSE,
  mgmissing = FALSE,
  vcov = csdm_vcov(),
  ...
) {
  model <- match.arg(model)
  trend <- match.arg(trend)

  if (inherits(data, "pdata.frame")) {
    idx <- attr(data, "index")
    if (!is.null(idx)) {
      id <- names(idx)[1L]
      time <- names(idx)[2L]
    }
  }

  panel_df <- .csdm_prepare_panel_df(data = data, id = id, time = time)

  if (trend == "pooled") {
    stop("trend='pooled' is not implemented yet")
  }
  if (trend == "unit") {
    panel_df$.csdm_trend__ <- .csdm_time_index(panel_df[[time]])
    formula <- stats::update(formula, . ~ . + .csdm_trend__)
  }

  fit <- switch(
    model,
    mg = .csdm_fit_mg(panel_df = panel_df, formula = formula, id = id, time = time, vcov = vcov, ...),
    cce = .csdm_fit_cce(panel_df = panel_df, formula = formula, id = id, time = time, csa = csa, vcov = vcov, ...),
    dcce = .csdm_fit_dcce(panel_df = panel_df, formula = formula, id = id, time = time, csa = csa, vcov = vcov, ...),
    cs_ardl = stop("Not implemented yet"),
    cs_ecm  = stop("Not implemented yet"),
    cs_dl   = stop("Not implemented yet")
  )

  fit$call <- match.call()
  fit$formula <- formula
  fit$model <- model
  fit$id <- id
  fit$time <- time
  fit$meta$trend <- trend
  fit$meta$fullsample <- isTRUE(fullsample)
  fit$meta$mgmissing <- isTRUE(mgmissing)
  fit$meta$lr <- lr
  fit$meta$pooled <- pooled
  fit$meta$vcov <- vcov

  class(fit) <- "csdm_fit"
  fit
}


.csdm_prepare_panel_df <- function(data, id, time) {
  if (inherits(data, "pdata.frame")) {
    df <- as.data.frame(data)
    idx <- attr(data, "index")
    if (!is.null(idx)) {
      id <- names(idx)[1L]
      time <- names(idx)[2L]
    }
  } else {
    df <- as.data.frame(data)
  }

  if (is.null(id) || is.null(time)) stop("Both 'id' and 'time' must be provided.")
  if (!all(c(id, time) %in% names(df))) stop("'data' must contain 'id' and 'time' columns.")

  df[[id]] <- as.character(df[[id]])

  o <- order(df[[id]], df[[time]])
  df <- df[o, , drop = FALSE]

  df$.csdm_rowname__ <- seq_len(nrow(df))
  rownames(df) <- as.character(df$.csdm_rowname__)

  df
}


.csdm_time_index <- function(time_vec) {
  tt <- sort(unique(time_vec))
  match(time_vec, tt)
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


.csdm_residual_matrix <- function(panel_df, id, time, res_long) {
  ids <- sort(unique(panel_df[[id]]))
  times <- sort(unique(panel_df[[time]]))

  E <- matrix(NA_real_, nrow = length(ids), ncol = length(times))
  rownames(E) <- as.character(ids)
  colnames(E) <- as.character(times)

  if (nrow(res_long) == 0L) return(E)

  ii <- match(as.character(res_long[[id]]), rownames(E))
  tt <- match(as.character(res_long[[time]]), colnames(E))
  keep <- is.finite(ii) & is.finite(tt)
  if (any(keep)) E[cbind(ii[keep], tt[keep])] <- res_long$residual[keep]
  E
}


.csdm_fit_mg <- function(panel_df, formula, id, time, vcov, ...) {
  ids <- unique(panel_df[[id]])

  # determine economic coefficient names from pooled model matrix (on complete rows)
  econ_names <- colnames(stats::model.matrix(formula, stats::model.frame(formula, panel_df, na.action = stats::na.omit)))

  coef_i <- matrix(NA_real_, nrow = length(ids), ncol = length(econ_names),
                   dimnames = list(as.character(ids), econ_names))

  res_long <- data.frame(
    id = character(0),
    time = character(0),
    residual = numeric(0),
    stringsAsFactors = FALSE
  )
  names(res_long)[1:2] <- c(id, time)

  dropped <- character(0)

  for (uid in ids) {
    sub <- panel_df[panel_df[[id]] == uid, , drop = FALSE]

    mf <- tryCatch(stats::model.frame(formula, sub, na.action = stats::na.omit), error = function(e) NULL)
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

    used_rows <- rownames(mf)
    sub_used <- sub[used_rows, , drop = FALSE]
    if (nrow(sub_used) != length(fit$residuals)) {
      # conservative fallback: skip residuals if mapping fails
      next
    }

    res_long <- rbind(
      res_long,
      data.frame(
        sub_used[[id]],
        sub_used[[time]],
        as.numeric(fit$residuals),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
    names(res_long)[1:3] <- c(id, time, "residual")
  }

  coef_mg <- colMeans(coef_i, na.rm = TRUE)

  V <- .csdm_mg_vcov(coef_i)
  se <- sqrt(diag(V))

  E <- .csdm_residual_matrix(panel_df, id, time, res_long)

  out <- list(
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
    meta = list(
      N = length(unique(panel_df[[id]])),
      T = length(unique(panel_df[[time]])),
      csa = NULL,
      dropped_units = dropped
    )
  )
  out
}


.csdm_fit_cce <- function(panel_df, formula, id, time, csa, vcov, ...) {
  # Validate CSA lags for static CCE
  if (length(csa$lags) == 1L && as.integer(csa$lags) != 0L) {
    stop("For model='cce', csa$lags must be 0")
  }
  if (length(csa$lags) > 1L && any(as.integer(csa$lags) != 0L)) {
    stop("For model='cce', all csa$lags entries must be 0")
  }

  # CSA vars
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

  econ_names <- colnames(stats::model.matrix(formula, stats::model.frame(formula, panel_df, na.action = stats::na.omit)))

  csa_terms <- if (length(csa_vars)) paste0("csa_", csa_vars) else character(0)
  fml <- if (length(csa_terms)) stats::update(formula, paste0(". ~ . + ", paste(csa_terms, collapse = " + "))) else formula

  ids <- unique(csa_attached[[id]])
  coef_i <- matrix(NA_real_, nrow = length(ids), ncol = length(econ_names),
                   dimnames = list(as.character(ids), econ_names))

  res_long <- data.frame(
    id = character(0),
    time = character(0),
    residual = numeric(0),
    stringsAsFactors = FALSE
  )
  names(res_long)[1:2] <- c(id, time)

  dropped <- character(0)

  for (uid in ids) {
    sub <- csa_attached[csa_attached[[id]] == uid, , drop = FALSE]

    mf <- tryCatch(stats::model.frame(fml, sub, na.action = stats::na.omit), error = function(e) NULL)
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

    used_rows <- rownames(mf)
    sub_used <- sub[used_rows, , drop = FALSE]
    if (nrow(sub_used) != length(fit$residuals)) next

    res_long <- rbind(
      res_long,
      data.frame(
        sub_used[[id]],
        sub_used[[time]],
        as.numeric(fit$residuals),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
    names(res_long)[1:3] <- c(id, time, "residual")
  }

  coef_mg <- colMeans(coef_i, na.rm = TRUE)
  V <- .csdm_mg_vcov(coef_i)
  se <- sqrt(diag(V))
  E <- .csdm_residual_matrix(csa_attached, id, time, res_long)

  out <- list(
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
    meta = list(
      N = length(unique(panel_df[[id]])),
      T = length(unique(panel_df[[time]])),
      csa = csa,
      dropped_units = dropped
    )
  )
  out
}


.csdm_fit_dcce <- function(panel_df, formula, id, time, csa, vcov, ...) {
  # CSA vars
  if (identical(csa$vars, "_none")) {
    csa_vars <- character(0)
  } else if (identical(csa$vars, "_all")) {
    csa_vars <- setdiff(unique(all.vars(formula)), ".csdm_trend__")
  } else {
    csa_vars <- as.character(csa$vars)
  }

  # build time-level CSA table first
  if (length(csa_vars)) {
    csa_time <- cross_sectional_avg(
      data = panel_df,
      id = id,
      time = time,
      vars = csa_vars,
      leave_out = FALSE,
      suffix = "csa",
      return_mode = "time",
      na.rm = TRUE
    )
  } else {
    csa_time <- unique(panel_df[, time, drop = FALSE])
  }

  # sort time table
  o <- order(csa_time[[time]])
  csa_time <- csa_time[o, , drop = FALSE]

  # add CSA lags by time only
  lag_spec <- csa$lags
  if (length(csa_vars) && length(lag_spec) > 1L) {
    # per-variable
    lag_spec <- lag_spec[csa_vars]
    lag_spec[is.na(lag_spec)] <- 0L
  }

  add_lag <- function(x, L) {
    if (L <= 0L) return(NULL)
    out <- vector("list", L)
    for (l in seq_len(L)) {
      out[[l]] <- c(rep(NA_real_, l), x[seq_len(length(x) - l)])
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

      # include contemporaneous
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

  # merge back to panel by time
  key <- match(as.character(panel_df[[time]]), as.character(csa_time[[time]]))
  csa_attached <- panel_df
  if (nrow(csa_time) && length(setdiff(names(csa_time), time))) {
    for (nm in setdiff(names(csa_time), time)) {
      csa_attached[[nm]] <- csa_time[[nm]][key]
    }
  }

  econ_names <- colnames(stats::model.matrix(formula, stats::model.frame(formula, panel_df, na.action = stats::na.omit)))

  fml <- if (length(csa_term_names)) stats::update(formula, paste0(". ~ . + ", paste(unique(csa_term_names), collapse = " + "))) else formula

  ids <- unique(csa_attached[[id]])
  coef_i <- matrix(NA_real_, nrow = length(ids), ncol = length(econ_names),
                   dimnames = list(as.character(ids), econ_names))

  res_long <- data.frame(
    id = character(0),
    time = character(0),
    residual = numeric(0),
    stringsAsFactors = FALSE
  )
  names(res_long)[1:2] <- c(id, time)

  dropped <- character(0)

  for (uid in ids) {
    sub <- csa_attached[csa_attached[[id]] == uid, , drop = FALSE]

    mf <- tryCatch(stats::model.frame(fml, sub, na.action = stats::na.omit), error = function(e) NULL)
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

    used_rows <- rownames(mf)
    sub_used <- sub[used_rows, , drop = FALSE]
    if (nrow(sub_used) != length(fit$residuals)) next

    res_long <- rbind(
      res_long,
      data.frame(
        sub_used[[id]],
        sub_used[[time]],
        as.numeric(fit$residuals),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
    names(res_long)[1:3] <- c(id, time, "residual")
  }

  coef_mg <- colMeans(coef_i, na.rm = TRUE)
  V <- .csdm_mg_vcov(coef_i)
  se <- sqrt(diag(V))
  E <- .csdm_residual_matrix(csa_attached, id, time, res_long)

  out <- list(
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
    meta = list(
      N = length(unique(panel_df[[id]])),
      T = length(unique(panel_df[[time]])),
      csa = csa,
      dropped_units = dropped
    )
  )
  out
}
