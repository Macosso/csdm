## Legacy implementation preserved for reference.
## Stored under inst/attic and NOT loaded by the package.

#' Cross-Sectionally Augmented ARDL (CS-ARDL) Mean Group Estimator
#'
#' @description
#' Estimates unit-specific ARDL(p, q) regressions augmented with cross-sectional
#' averages (CSAs) of the dependent variable and regressors - plus their lags up to `P` -
#' then averages the coefficients across units (Mean Group). The CSA terms proxy
#' latent common factors and mitigate cross-sectional dependence.
#'
#' @noRd
cs_ardl_mg <- function(formula, data, id, time,
                       p = 1, q = 1, P = 1, loo = TRUE,
                       trend = FALSE,
                       robust_se = c("none","HC1","HC3"),
                       min_T = 10) {

  robust_se <- match.arg(robust_se)
  stopifnot(p >= 1, P >= 0)

  tf <- stats::terms(formula)
  yname  <- as.character(attr(tf, "variables"))[2L]
  xnames <- attr(tf, "term.labels")
  if (length(xnames) == 0L) stop("Formula must have regressors on the RHS.")

  if (length(q) == 1L) {
    qvec <- rep.int(as.integer(q), length(xnames)); names(qvec) <- xnames
  } else {
    if (is.null(names(q))) stop("If `q` has length>1, it must be a *named* vector matching RHS vars.")
    qvec <- as.integer(q[xnames])
    if (any(is.na(qvec))) stop("`q` must provide lags for all RHS variables in the same order/names.")
  }

  df <- as.data.frame(data)
  if (!all(c(id, time, yname, xnames) %in% names(df))) {
    stop("`data` must contain id, time, y, and all RHS variables.")
  }

  o  <- order(df[[id]], df[[time]])
  df <- df[o, , drop = FALSE]

  csa_frame <- stats::aggregate(df[, c(yname, xnames)],
                         list(tt = df[[time]]),
                         function(v) mean(v, na.rm = TRUE))
  names(csa_frame)[1L] <- time

  lagv <- function(v, L) {
    if (L == 0L) return(cbind(v))
    out <- matrix(NA_real_, nrow = length(v), ncol = L + 1L)
    for (j in 0:L) out[(j+1L):length(v), j+1L] <- v[seq_len(length(v)-j)]
    out
  }

  unit_ids <- unique(df[[id]])
  res_list <- vector("list", length(unit_ids))
  names(res_list) <- unit_ids

  cd_store <- list(e = numeric(0), id = vector("list", 0), tt = numeric(0))

  for (ii in seq_along(unit_ids)) {
    uid <- unit_ids[ii]
    dfi <- df[df[[id]] == uid, , drop = FALSE]
    dfi <- merge(dfi, csa_frame, by = time, suffixes = c("", ".csa"), sort = FALSE)

    if (loo) {
      tt <- dfi[[time]]
      idx_t <- match(tt, csa_frame[[time]])
      agg_sum <- stats::aggregate(df[, c(yname, xnames)], list(tt = df[[time]]), function(v) sum(v, na.rm = TRUE))
      agg_cnt <- stats::aggregate(df[, c(yname, xnames)], list(tt = df[[time]]), function(v) sum(!is.na(v)))
      names(agg_sum)[1L] <- names(agg_cnt)[1L] <- time
      sums <- agg_sum[idx_t, c(yname, xnames), drop = FALSE]
      cnts <- agg_cnt[idx_t, c(yname, xnames), drop = FALSE]
      for (vn in c(yname, xnames)) {
        num <- sums[[vn]] - dfi[[vn]]
        den <- pmax(cnts[[vn]] - 1L, 1L)
        dfi[[paste0(vn, ".csa")]] <- ifelse(den > 0, num / den, NA_real_)
      }
    }

    T_i <- nrow(dfi)
    if (T_i < (max(p, max(qvec)) + P + 2L)) next

    ylagM <- lagv(dfi[[yname]], p)
    colnames(ylagM) <- c(yname, paste0("L", 1:p, ".", yname))

    xlag_list <- list()
    for (k in seq_along(xnames)) {
      XK  <- xnames[k]; qk <- qvec[k]
      xlagM <- lagv(dfi[[XK]], qk)
      colnames(xlagM) <- paste0(XK, "L", 0:qk)
      xlag_list[[XK]] <- xlagM
    }

    ycsaM <- lagv(dfi[[paste0(yname, ".csa")]], P)
    colnames(ycsaM) <- paste0("ycsaL", 0:P)

    xcsa_list <- list()
    for (XK in xnames) {
      xcsaM <- lagv(dfi[[paste0(XK, ".csa")]], P)
      colnames(xcsaM) <- paste0(XK, "csaL", 0:P)
      xcsa_list[[XK]] <- xcsaM
    }

    Xreg <- cbind(
      ylagM[, -1, drop = FALSE],
      do.call(cbind, xlag_list),
      ycsaM,
      do.call(cbind, xcsa_list)
    )
    yreg <- ylagM[, 1]
    keep <- rowSums(is.na(cbind(yreg, Xreg))) == 0
    yreg <- yreg[keep]
    Xreg <- Xreg[keep, , drop = FALSE]
    Ti_eff <- length(yreg)
    if (Ti_eff < min_T) next
    if (trend) Xreg <- cbind(Xreg, trend = seq_len(nrow(Xreg)))

    df_i <- data.frame(y = yreg, Xreg, check.names = FALSE)
    fit  <- lm(y ~ ., data = df_i)

    r2_i   <- summary(fit)$r.squared
    coef_i <- coef(fit)

    lag_y_names <- paste0("L", seq_len(p), ".", yname)
    alpha_y_lags <- coef_i[lag_y_names]
    phi_sum_i <- sum(alpha_y_lags, na.rm = TRUE)

    bx_names <- unlist(lapply(xnames, function(XK) paste0(XK, "L", 0:qvec[XK])))
    bx_hat <- coef_i[bx_names]

    theta_i <- vapply(xnames, function(XK) {
      num <- sum(coef_i[paste0(XK, "L", 0:qvec[XK])], na.rm = TRUE)
      den <- 1 - phi_sum_i
      if (is.finite(num) && is.finite(den) && abs(den) > 1e-8) num / den else NA_real_
    }, numeric(1))

    kept_times <- dfi[[time]][keep]
    cd_store$e  <- c(cd_store$e, stats::residuals(fit))
    cd_store$id <- c(cd_store$id, rep(list(uid), length(stats::residuals(fit))))
    cd_store$tt <- c(cd_store$tt, kept_times)

    res_list[[ii]] <- list(
      id    = uid,
      coef  = coef_i,
      phi_sum = phi_sum_i,
      bx      = bx_hat,
      theta   = theta_i,
      r2    = r2_i,
      T_eff = Ti_eff
    )
  }

  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0L) stop("No units passed the lagging/min_T filters.")

  r2_mg <- mean(vapply(res_list, function(z) z$r2, numeric(1)), na.rm = TRUE)

  df_e <- data.frame(id = unlist(cd_store$id),
                     time = cd_store$tt,
                     e = cd_store$e, check.names = FALSE)
  ids_used   <- sort(unique(df_e$id))
  times_used <- sort(unique(df_e$time))
  Ew <- matrix(NA_real_, nrow = length(times_used), ncol = length(ids_used))
  rownames(Ew) <- times_used; colnames(Ew) <- ids_used
  for (jj in seq_len(nrow(df_e))) {
    rr <- df_e$time[jj]; cc <- as.character(df_e$id[jj])
    Ew[as.character(rr), as.character(cc)] <- df_e$e[jj]
  }
  N <- ncol(Ew); Tbar <- nrow(Ew)
  C  <- suppressWarnings(stats::cor(Ew, use = "pairwise.complete.obs"))
  iu <- upper.tri(C, diag = FALSE)
  rho_bar <- mean(C[iu], na.rm = TRUE)
  cd_stat <- sqrt(2 * N * (N - 1) / Tbar) * rho_bar
  cd_pval <- 2 * (1 - pnorm(abs(cd_stat)))

  all_bx_names <- sort(unique(unlist(lapply(res_list, function(z) names(z$bx)))))
  BX <- do.call(rbind, lapply(res_list, function(z) {
    out <- rep(NA_real_, length(all_bx_names)); names(out) <- all_bx_names
    out[names(z$bx)] <- z$bx
    out
  }))
  bx_mg <- colMeans(BX, na.rm = TRUE)
  bx_se <- apply(BX, 2, function(v) stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))

  phi_vec <- vapply(res_list, function(z) z$phi_sum, numeric(1))
  phi_mg  <- c(
    mean = mean(phi_vec, na.rm = TRUE),
    se   = stats::sd(phi_vec,  na.rm = TRUE) / sqrt(sum(!is.na(phi_vec)))
  )

  TH <- do.call(rbind, lapply(res_list, function(z) z$theta))
  colnames(TH) <- xnames
  theta_mg <- colMeans(TH, na.rm = TRUE)
  theta_se <- apply(TH, 2, function(v) stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))

  units_tbl <- do.call(rbind, lapply(res_list, function(z) {
    cbind(
      data.frame(id = z$id, r2 = z$r2, phi_sum = z$phi_sum, T_eff = z$T_eff,
                 stringsAsFactors = FALSE, check.names = FALSE),
      t(z$theta)
    )
  }))
  rownames(units_tbl) <- NULL

  out <- list(
    call = match.call(),
    coef_mg_short = data.frame(term = names(bx_mg), estimate = unname(bx_mg), se = unname(bx_se)),
    coef_mg_long  = data.frame(variable = names(theta_mg), estimate = unname(theta_mg), se = unname(theta_se)),
    phi_mg  = phi_mg,
    r2_mg   = r2_mg,
    cd_stat = cd_stat,
    cd_pval = cd_pval,
    units   = units_tbl
  )
  class(out) <- "cs_ardl_mg"
  return(out)
}


print.cs_ardl_mg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  stopifnot(inherits(x, "cs_ardl_mg"))
  N_eff <- if (!is.null(x$units)) nrow(x$units) else NA_integer_
  cat("CS-ARDL Mean Group fit\n")
  cat("Call: "); print(x$call)
  cat(sprintf("Units (N): %s\n", N_eff))
  if (!is.null(x$r2_mg)) cat(sprintf("Mean R-squared: %.*f\n", digits, x$r2_mg))
  if (!is.null(x$phi_mg)) {
    cat(sprintf("Mean sum(phi_lags): %.*f  (SE: %.*f)\n",
                digits, x$phi_mg["mean"], digits, x$phi_mg["se"]))
  }
  if (!is.null(x$cd_stat)) {
    cat(sprintf("CD test: stat = %.*f, p = %.*g\n",
                digits, x$cd_stat, digits, x$cd_pval))
  }

  if (!is.null(x$coef_mg_long)) {
    cat("\nLong-run MG coefficients (preview):\n")
    print(utils::head(x$coef_mg_long, 10L), row.names = FALSE)
  }
  if (!is.null(x$coef_mg_short)) {
    cat("\nShort-run MG coefficients (preview):\n")
    print(utils::head(x$coef_mg_short, 10L), row.names = FALSE)
  }
  invisible(x)
}


summary.cs_ardl_mg <- function(object, z_dist = TRUE, ...) {
  stopifnot(inherits(object, "cs_ardl_mg"))
  N_eff <- if (!is.null(object$units)) nrow(object$units) else NA_integer_
  df <- if (is.na(N_eff)) NA_integer_ else max(N_eff - 1L, 1L)

  add_stats <- function(tab) {
    if (is.null(tab) || nrow(tab) == 0L) return(tab)
    se <- tab$se
    zt <- tab$estimate / se
    if (z_dist) {
      p <- 2 * stats::pnorm(abs(zt), lower.tail = FALSE)
    } else {
      p <- 2 * stats::pt(abs(zt), df = df, lower.tail = FALSE)
    }
    cbind(tab, statistic = zt, p.value = p, row.names = NULL)
  }

  long_tab  <- add_stats(object$coef_mg_long)
  short_tab <- add_stats(object$coef_mg_short)

  out <- list(
    call   = object$call,
    N      = N_eff,
    r2_mg  = object$r2_mg,
    phi_mg = object$phi_mg,
    cd_stat = object$cd_stat,
    cd_pval = object$cd_pval,
    long_run  = long_tab,
    short_run = short_tab,
    z_dist = z_dist,
    df     = if (z_dist) Inf else df
  )
  class(out) <- "summary.cs_ardl_mg"
  out
}


print.summary.cs_ardl_mg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     signif.stars = TRUE, ...) {
  stopifnot(inherits(x, "summary.cs_ardl_mg"))
  cat("CS-ARDL Mean Group - Summary\n")
  cat("Call: "); print(x$call)
  cat(sprintf("Units (N): %s\n", x$N))
  if (!is.null(x$r2_mg)) cat(sprintf("Mean R-squared: %.*f\n", digits, x$r2_mg))
  if (!is.null(x$phi_mg)) {
    cat(sprintf("Mean sum(phi_lags): %.*f  (SE: %.*f)\n",
                digits, x$phi_mg["mean"], digits, x$phi_mg["se"]))
  }
  if (!is.null(x$cd_stat)) {
    cat(sprintf("CD test: stat = %.*f, p = %.*g\n\n",
                digits, x$cd_stat, digits, x$cd_pval))
  }

  pfmt <- function(tab, title) {
    if (is.null(tab) || nrow(tab) == 0L) return(invisible())
    cat(title, "\n", sep = "")
    if (signif.stars && "p.value" %in% names(tab)) {
      stars <- stats::symnum(tab$p.value, corr = FALSE,
                      cutpoints = c(0, .001, .01, .05, .1, 1),
                      symbols = c("***","**","*","."," "))
      tab$` ` <- stars
    }
    print(tab, digits = digits, row.names = FALSE)
    cat("\n")
  }

  pfmt(x$long_run,  "Long-run MG coefficients:")
  pfmt(x$short_run, "Short-run MG coefficients:")

  dist_str <- if (isTRUE(is.infinite(x$df))) "z" else sprintf("t(df=%d)", x$df)
  cat(sprintf("Test statistics use %s distribution for p-values.\n", dist_str))
  invisible(x)
}


coef.cs_ardl_mg <- function(object, type = c("long_run","short_run","phi"), ...) {
  stopifnot(inherits(object, "cs_ardl_mg"))
  type <- match.arg(type)
  if (type == "long_run") {
    with(object$coef_mg_long, stats::setNames(estimate, variable))
  } else if (type == "short_run") {
    with(object$coef_mg_short, stats::setNames(estimate, term))
  } else {
    stats::setNames(object$phi_mg["mean"], "phi_sum")
  }
}


vcov.cs_ardl_mg <- function(object, type = c("long_run","short_run"), ...) {
  stopifnot(inherits(object, "cs_ardl_mg"))
  type <- match.arg(type)
  if (type == "long_run") {
    nm <- object$coef_mg_long$variable
    se <- object$coef_mg_long$se
  } else {
    nm <- object$coef_mg_short$term
    se <- object$coef_mg_short$se
  }
  V <- diag(se^2, nrow = length(se), ncol = length(se))
  dimnames(V) <- list(nm, nm)
  V
}
