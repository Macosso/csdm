# csdm_fit_engines.R

.csdm_fit_mg <- function(panel_df, formula, id, time, lr = NULL, vcov, ...) {
  .csdm_fit_units(panel_df, formula, id, time)
}


.csdm_fit_cce <- function(panel_df, formula, id, time, csa, lr = NULL, vcov,
                          fullsample = FALSE, ...) {
  if (length(csa$lags) == 1L && as.integer(csa$lags) != 0L) {
    stop("For model='cce', csa$lags must be 0")
  }
  if (length(csa$lags) > 1L && any(as.integer(csa$lags) != 0L)) {
    stop("For model='cce', all csa$lags entries must be 0")
  }

  augmented <- .csdm_augment(panel_df, formula, formula, id, time, csa, fullsample)
  .csdm_fit_units(augmented$data, augmented$formula, id, time, formula, "cce", augmented$csa)
}


.csdm_fit_dcce <- function(panel_df, formula, id, time, csa, lr, vcov,
                           fullsample = FALSE, ...) {
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
          lagv <- .csdm_lag(yv, panel_work[[time]][idxo], k, attr(panel_df, "csdm_time_step"))
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
            lagv <- .csdm_lag(xv, panel_work[[time]][idxo], k, attr(panel_df, "csdm_time_step"))
            panel_work[[paste0("lag", k, "_", xn)]][idxo] <- lagv
          }
        }
      }
    }

    if (length(lag_terms_all)) {
      econ_formula <- stats::update(formula, paste0(". ~ . + ", paste(lag_terms_all, collapse = " + ")))
    }
  }

  augmented <- .csdm_augment(
    panel_work, formula, econ_formula, id, time, csa, fullsample
  )
  .csdm_fit_units(augmented$data, augmented$formula, id, time, econ_formula, "dcce", augmented$csa)
}


.csdm_fit_cs_ardl <- function(panel_df, formula, id, time, csa, lr, vcov,
                              fullsample = FALSE, ...) {
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
    fullsample = fullsample,
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

  ok_denom <- is.finite(denom_i) & abs(denom_i) > 1e-8
  lr_i <- matrix(NA_real_, nrow(coef_i), length(xnames),
    dimnames = list(unit_ids, if (length(xnames)) paste0("lr_", xnames) else character()))
  for (j in seq_along(xnames)) {
    beta_terms <- c(xnames[j], if (lr_xdlags > 0L) paste0("lag", seq_len(lr_xdlags), "_", xnames[j]))
    if (!all(beta_terms %in% colnames(coef_i))) stop("Long-run effects require numeric scalar regressors.")
    lr_i[ok_denom, j] <- rowSums(coef_i[ok_denom, beta_terms, drop = FALSE]) / denom_i[ok_denom]
  }
  adjustment <- matrix(adj_i, ncol = 1L, dimnames = list(unit_ids, paste0("lr_", yname)))
  combined <- cbind(coef_i, adjustment, lr_i)
  if (anyDuplicated(colnames(combined))) stop("Long-run parameter names collide with economic terms.")
  components <- lapply(list(all = combined, levels = coef_i, adjustment = adjustment, long_run = lr_i),
    .csdm_parameter_component)
  excluded_lr <- setdiff(fit$meta$included_units, components$all$units)
  if (length(excluded_lr)) warning("Undefined long-run ratios for unit(s): ", paste(excluded_lr, collapse = ", "),
    ". Component samples are recorded in fit$components.", call. = FALSE)
  roots <- lapply(seq_len(nrow(alpha_mat)), function(i) {
    a <- alpha_mat[i, ]
    if (any(!is.finite(a))) return(NA_complex_)
    polyroot(c(1, -a))
  })
  names(roots) <- unit_ids
  stable <- vapply(roots, function(r) if (anyNA(r)) NA else all(Mod(r) > 1), logical(1))
  fit$components <- components
  default <- components$all
  lr_terms <- colnames(lr_i)
  fit$cs_ardl <- list(y = yname, x = xnames, ylags = lr_ylags, xdlags = lr_xdlags,
    unit_ids = unit_ids, denom_i = stats::setNames(denom_i, unit_ids),
    adj_i = stats::setNames(adj_i, unit_ids), lr_i = lr_i, ar_roots = roots, stable = stable,
    mg = list(adj = c(estimate = default$coefficients[colnames(adjustment)],
      se = sqrt(diag(default$vcov))[colnames(adjustment)]),
      lr = data.frame(term = lr_terms, estimate = unname(default$coefficients[lr_terms]),
        se = unname(sqrt(diag(default$vcov))[lr_terms]), n_used = rep(default$n_used, length(lr_terms)))))
  names(fit$cs_ardl$mg$adj) <- c("estimate", "se")
  fit
}

.csdm_parameter_component <- function(B) {
  keep <- rowSums(!is.finite(B)) == 0L
  used <- B[keep, , drop = FALSE]
  estimate <- if (nrow(used)) colMeans(used) else stats::setNames(rep(NA_real_, ncol(B)), colnames(B))
  list(coefficients = estimate, vcov = .csdm_mg_vcov(used),
    n_used = nrow(used), units = rownames(used), unit_coefficients = used)
}
