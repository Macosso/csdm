#' Cross-Sectionally Augmented ARDL (CS-ARDL) Mean Group Estimator
#'
#' @description
#' Estimates unit-specific ARDL(p, q) regressions augmented with cross-sectional
#' averages (CSAs) of the dependent variable and regressors—plus their lags up to `P`—
#' then averages the coefficients across units (Mean Group). The CSA terms proxy
#' latent common factors and mitigate cross-sectional dependence.
#'
#' @details
#' **Model.** For unit \(i=1,\dots,N\) and time \(t\), let \(y_{it}\) be the dependent
#' variable and \(x_{it}\) a K-vector of regressors. A CS-ARDL(p, q, P) for unit \(i\) is
#' \deqn{
#' y_{it} = \alpha_i + \sum_{j=1}^{p} \phi_{ij} y_{i,t-j}
#'          + \sum_{k=1}^{K} \sum_{s=0}^{q_k} \beta_{ik,s} x_{k,i,t-s}
#'          + \sum_{\ell=0}^{P} \delta_{i\ell}^{(y)} \bar{y}_{t-\ell}
#'          + \sum_{k=1}^{K} \sum_{\ell=0}^{P} \delta_{ik,\ell}^{(x)} \bar{x}_{k,t-\ell}
#'          + u_{it},
#' }
#' where \(\bar{y}_t\) and \(\bar{x}_{k,t}\) are cross-sectional averages (CSAs) that can be
#' computed either **including** or **excluding** unit \(i\) (leave-one-out, `loo=TRUE`).
#'
#' **Long-run effects.** For each unit \(i\), the long-run coefficient on regressor \(k\) is
#' \(\theta_{ik} = \big(\sum_{s=0}^{q_k} \beta_{ik,s}\big) / \big(1 - \sum_{j=1}^{p} \phi_{ij}\big)\),
#' provided the AR polynomial is stable. The function reports MG means of \(\theta_{ik}\) and
#' their MG standard errors (cross-sectional dispersion / \(\sqrt{N}\)).
#'
#' **Unbalanced panels.** Estimation drops rows with insufficient lags; CSAs at time \(t\)
#' are computed using all units with non-missing values at \(t\) (and excluding \(i\) if `loo=TRUE`).
#'
#' **Returns.**
#' - `coef_mg_short`: MG averages of the unit-level **short-run** coefficients (on the lagged \(y\) and lagged \(x\)'s).
#' - `coef_mg_long`: MG averages of **long-run** coefficients \(\theta_k\) with MG SEs.
#' - `phi_mg`: MG average of the sum of lagged-\(y\) coefficients \(\sum_j \phi_{ij}\) and its SE.
#' - `r2_mg`: average of unit-level R² from each OLS.
#' - `cd_stat`, `cd_pval`: Pesaran-style CD test on pooled residuals.
#' - `units`: tibble with unit-level results (coefficients, long-run, R², residuals index).
#'
#' **Related literature in your PDFs.**
#' - MG/PMG background and ARDL motivation: Pesaran & Smith (1995); Pesaran, Shin & Smith (1999).
#' - CS augmentation to address unobserved common factors: see the dynamic CCE/CS-ARDL approach
#'   developed by Chudik & Pesaran (2015).
#' - Error-correction and panel cointegration testing context: Westerlund (2007).
#'
#' This implementation is aligned with those ideas but written from scratch for the package.
#'
#' @param formula An object like `y ~ x1 + x2 + ...`. All regressors on RHS are treated as levels
#'   for ARDL construction; lags are generated internally.
#' @param data A data.frame (or tibble) with columns for `id`, `time`, `y`, `x`'s.
#' @param id,time Column names (strings) identifying unit and time indices.
#' @param p Integer \(\ge 1\). Number of lags of \(y\).
#' @param q Integer (recycled) or named integer vector for lags of each regressor
#'   (names must match RHS variable names). `q` counts max lag; lag 0 is included by default.
#' @param P Integer \(\ge 0\). Number of lags of cross-sectional averages added (lag 0..P).
#' @param loo Logical. If `TRUE`, compute CSAs in leave-one-out fashion.
#' @param trend Logical. If `TRUE`, adds a deterministic linear trend (common slope) to each unit regression.
#' @param robust_se One of `c("none","HC1","HC3")` for unit OLS variance; MG SEs are based on cross-sectional dispersion of unit estimates.
#' @param min_T Integer. Minimum time observations required per unit after lagging; units with fewer are dropped.
#' @return A list with MG estimates, their SEs, CD test, and unit-level results (see Details).
#' @export
cs_ardl_mg <- function(formula, data, id, time,
                       p = 1, q = 1, P = 1, loo = TRUE,
                       trend = FALSE,
                       robust_se = c("none","HC1","HC3"),
                       min_T = 10) {

  # --- Helpers ---------------------------------------------------------------
  robust_se <- match.arg(robust_se)
  stopifnot(p >= 1, P >= 0)

  # Parse formula
  tf <- terms(formula)
  yname <- as.character(attr(tf, "variables"))[2L]
  xnames <- attr(tf, "term.labels")
  if (length(xnames) == 0L) stop("Formula must have regressors on the RHS.")

  # Normalize q to named integer vector
  if (length(q) == 1L) {
    qvec <- rep.int(as.integer(q), length(xnames))
    names(qvec) <- xnames
  } else {
    if (is.null(names(q))) stop("If `q` has length>1, it must be a *named* vector matching RHS vars.")
    qvec <- as.integer(q[xnames])
    if (any(is.na(qvec))) stop("`q` must provide lags for all RHS variables in the same order/names.")
  }

  df <- as.data.frame(data)
  if (!all(c(id, time, yname, xnames) %in% names(df))) {
    stop("`data` must contain id, time, y, and all RHS variables.")
  }

  # Order by id, time
  o <- order(df[[id]], df[[time]])
  df <- df[o, , drop = FALSE]

  # Build per-time CSAs (including or excluding i later)
  # For efficiency, precompute inclusive averages; for loo, adjust per unit.
  csa_frame <- aggregate(df[, c(yname, xnames)], list(tt = df[[time]]), function(v) mean(v, na.rm = TRUE))
  names(csa_frame)[1L] <- time

  # Function to create lags for a vector with NA padding
  lagv <- function(v, L) {
    if (L == 0L) return(cbind(v))
    out <- matrix(NA_real_, nrow = length(v), ncol = L + 1L)
    for (j in 0:L) out[(j + 1L):length(v), j + 1L] <- v[seq_len(length(v) - j)]
    out
  }

  # Storage
  unit_ids <- unique(df[[id]])
  res_list <- vector("list", length(unit_ids))
  names(res_list) <- unit_ids

  # For CD test, collect residuals with (id,time)
  cd_store <- list(e = numeric(0), id = vector("list", 0), tt = numeric(0))

  # Iterate units
  for (ii in seq_along(unit_ids)) {
    uid <- unit_ids[ii]
    dfi <- df[df[[id]] == uid, , drop = FALSE]

    # Merge with CSA (inclusive), keep only times present in unit
    dfi <- merge(dfi, csa_frame, by = time, suffixes = c("", ".csa"), sort = FALSE)

    # If leave-one-out needed, re-compute CSAs excluding i for those times
    if (loo) {
      # Compute LOO means time-by-time using all available units at each time
      # Efficient approach: for each time, (sum_all - value_i) / (count_all - 1)
      # Precompute sums & counts over all units by time for y and each x
      # We do it on-the-fly for observed rows of this unit
      tt <- dfi[[time]]
      idx_t <- match(tt, csa_frame[[time]])
      # sums & counts
      # Build aggregate sums and counts once:
      # (We reconstruct here for clarity; could be cached for speed.)
      agg_sum <- aggregate(df[, c(yname, xnames)], list(tt = df[[time]]), function(v) sum(v, na.rm = TRUE))
      agg_cnt <- aggregate(df[, c(yname, xnames)], list(tt = df[[time]]), function(v) sum(!is.na(v)))
      names(agg_sum)[1L] <- names(agg_cnt)[1L] <- time

      # Map for this unit's rows
      sums <- agg_sum[idx_t, c(yname, xnames), drop = FALSE]
      cnts <- agg_cnt[idx_t, c(yname, xnames), drop = FALSE]

      # LOO CSAs
      for (vn in c(yname, xnames)) {
        num <- sums[[vn]] - dfi[[vn]]
        den <- pmax(cnts[[vn]] - 1L, 1L)  # avoid zero division
        dfi[[paste0(vn, ".csa")]] <- ifelse(den > 0, num / den, NA_real_)
      }
    }

    # Build design: y_t, lags of y (1..p), lags of each x (0..q_k),
    # CSA(y) lags (0..P), CSA(x_k) lags (0..P), optional trend
    T_i <- nrow(dfi)
    if (T_i < (max(p, max(qvec)) + P + 2L)) next

    # Generate lag blocks
    ylagM <- lagv(dfi[[yname]], p)
    colnames(ylagM) <- c("y0", paste0("yL", 1:p))  # y0 is current y

    xlag_list <- list()
    for (k in seq_along(xnames)) {
      XK <- xnames[k]
      qk <- qvec[k]
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

    # Assemble regression frame, drop initial rows with NAs due to lagging
    Xreg <- cbind(
      ylagM[, -1, drop = FALSE],                     # y lags (no current y)
      do.call(cbind, xlag_list),                     # x lags incl. L0
      ycsaM,                                          # CSA y lags
      do.call(cbind, xcsa_list)                      # CSA x lags
    )
    yreg <- ylagM[, 1]                                # current y

    keep <- rowSums(is.na(cbind(yreg, Xreg))) == 0
    yreg <- yreg[keep]
    Xreg <- Xreg[keep, , drop = FALSE]
    Ti_eff <- length(yreg)
    if (Ti_eff < min_T) next

    if (trend) {
      Xreg <- cbind(Xreg, trend = seq_len(nrow(Xreg)))
    }

    # OLS with optional robust SE at unit level (for reporting only; MG SEs are cross-sectional)
    df_i <- data.frame(y = yreg, Xreg, check.names = FALSE)
    fit <- lm(y ~ ., data = df_i)

    # Unit-level R2
    r2_i <- summary(fit)$r.squared

    coef_i <- coef(fit)
    # per unit i, after estimating the unit regression
    alpha_y_levels <- coef_i[grepl("^L[1-9]+\\.", names(coef_i)) & names(coef_i) == paste0("L1.", yname)]

    # If LR block can have more L.y levels, adapt the grepl to pick only those in LR.
    phi_i <- -(1 - sum(alpha_y_levels, na.rm = TRUE))


    # Short-run x-coefs (each lag separately)

    bx_names <- unlist(lapply(xnames, function(XK) paste0(XK, "L", 0:qvec[XK])))
    bx_hat <- coef(fit)[bx_names]
    # # Long-run theta_i per regressor
    # theta_i <- vapply(xnames, function(XK) {
    #   num <- sum(coef(fit)[paste0(XK, "L", 0:qvec[XK])], na.rm = TRUE)
    #   den <- (1 - phi_i)
    #   if (is.finite(num) && is.finite(den) && abs(den) > 1e-8) num / den else NA_real_
    # }, numeric(1))
    #

    # x-level coefficients that belong to the LR block (e.g., contemporaneous levels of x)
    bx_levels <- coef_i[xnames]                 # e.g., c("log_hc","log_ck","log_ngd")
    theta_i   <- bx_levels / (1 - sum(alpha_y_levels, na.rm = TRUE))


    # Collect residuals for CD test with (id,time)
    # Map kept rows back to original times
    kept_times <- dfi[[time]][keep]
    cd_store$e <- c(cd_store$e, residuals(fit))
    cd_store$id <- c(cd_store$id, rep(list(uid), length(residuals(fit))))
    cd_store$tt <- c(cd_store$tt, kept_times)

    res_list[[ii]] <- list(
      id = uid,
      coef = coef(fit),
      phi_sum = phi_i,
      bx = bx_hat,
      theta = theta_i,
      r2 = r2_i,
      T_eff = Ti_eff
    )
  }

  # Drop NULLs (units that didn't pass min_T)
  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0L) stop("No units passed the lagging/min_T filters.")

  # Average R2
  r2_mg <- mean(vapply(res_list, function(z) z$r2, numeric(1)), na.rm = TRUE)

  # CD test (Pesaran-type) on pooled residuals
  # Build residual matrix E_{it} aligned by (id,time)
  # 1) Identify common time points
  df_e <- data.frame(id = unlist(cd_store$id), time = cd_store$tt, e = cd_store$e, check.names = FALSE)
  # Wide by id
  ids_used <- sort(unique(df_e$id))
  times_used <- sort(unique(df_e$time))
  Ew <- matrix(NA_real_, nrow = length(times_used), ncol = length(ids_used))
  rownames(Ew) <- times_used; colnames(Ew) <- ids_used
  # Fill
  for (jj in seq_len(nrow(df_e))) {
    rr <- df_e$time[jj]; cc <- as.character(df_e$id[jj])
    Ew[as.character(rr), as.character(cc)] <- df_e$e[jj]
  }
  N <- ncol(Ew); Tbar <- nrow(Ew)

  # De-mean by time if desired (classic CD uses raw residuals); we use raw here
  # Compute pairwise correlations at each time alignment (listwise over time)
  # Implementation: compute correlation matrix of columns, omitting rows with any NA pairwise.
  # We'll use pairwise-complete obs correlation.
  C <- suppressWarnings(stats::cor(Ew, use = "pairwise.complete.obs"))
  # Upper triangular mean correlation
  iu <- upper.tri(C, diag = FALSE)
  rho_bar <- mean(C[iu], na.rm = TRUE)
  # CD statistic ~ N(0,1) under weak dependence

  cd_stat <- sqrt(2 * N * (N - 1) / Tbar) * rho_bar
  cd_pval <- 2 * (1 - pnorm(abs(cd_stat)))

  # MG: short-run (align names)
  all_bx_names <- sort(unique(unlist(lapply(res_list, function(z) names(z$bx)))))
  BX <- do.call(rbind, lapply(res_list, function(z) {
    out <- rep(NA_real_, length(all_bx_names)); names(out) <- all_bx_names
    out[names(z$bx)] <- z$bx
    out
  }))
  bx_mg <- colMeans(BX, na.rm = TRUE)
  bx_se <- apply(BX, 2, function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))

  # MG: phi sums
  phi_vec <- vapply(res_list, function(z) z$phi_sum, numeric(1))
  phi_mg   <- mean(phi_vec,  na.rm = TRUE)
  phi_se   <- sqrt( sum((phi_vec - phi_mg)^2, na.rm = TRUE) / (N*(N-1)) )

  # MG: long-run
  TH <- do.call(rbind, lapply(res_list, function(z) z$theta))
  colnames(TH) <- xnames
  # theta_mg <- colMeans(TH, na.rm = TRUE)
  # theta_se <- apply(TH, 2, function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))


  theta_mg <- colMeans(TH, na.rm = TRUE)
  theta_se <- sqrt( colSums( sweep(TH, 2, theta_mg, "-")^2, na.rm = TRUE ) / (N*(N-1)) )


  # Unit table
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
    phi_mg = phi_mg,
    r2_mg = r2_mg,
    cd_stat = cd_stat,
    cd_pval = cd_pval,
    units = units_tbl
  )
  class(out) <- "cs_ardl_mg"

  return(out)
}



#' @title Print method for cs_ardl_mg objects
#' @description Compact display of key results.
#' @param x A \code{cs_ardl_mg} object.
#' @param digits Number of digits to print.
#' @param ... Not used.
#' @return \code{invisible(x)}.
#' @export
#' @method print cs_ardl_mg
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

  # Show a small preview of coefficients
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


#' @title Summary for cs_ardl_mg objects
#' @description Adds test statistics and p-values to MG coefficients.
#' @param object A \code{cs_ardl_mg} object.
#' @param z_dist Logical; if \code{TRUE}, use standard normal for p-values (default).
#'   If \code{FALSE}, uses a t-approximation with df = N - 1.
#' @param ... Not used.
#' @return An object of class \code{summary.cs_ardl_mg}.
#' @export
#' @method summary cs_ardl_mg
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

#' @title Print method for summary.cs_ardl_mg
#' @description Nicely formats the summary output.
#' @param x A \code{summary.cs_ardl_mg} object.
#' @param digits Number of digits.
#' @param signif.stars Logical; print significance stars.
#' @param ... Not used.
#' @return \code{invisible(x)}.
#' @export
#' @method print summary.cs_ardl_mg
print.summary.cs_ardl_mg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     signif.stars = TRUE, ...) {
  stopifnot(inherits(x, "summary.cs_ardl_mg"))
  cat("CS-ARDL Mean Group — Summary\n")
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
    # basic stars
    if (signif.stars && "p.value" %in% names(tab)) {
      stars <- symnum(tab$p.value, corr = FALSE,
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


#' @title Coefficients for cs_ardl_mg
#' @description Extracts MG coefficients.
#' @param object A \code{cs_ardl_mg} object.
#' @param type One of \code{"long_run"}, \code{"short_run"}, or \code{"phi"}.
#'   \code{"long_run"} returns a named numeric vector of MG long-run estimates;
#'   \code{"short_run"} returns a named numeric vector of MG short-run estimates;
#'   \code{"phi"} returns a named scalar \code{c(phi_sum = ...)}.
#' @param ... Not used.
#' @return A named numeric vector.
#' @export
#' @method coef cs_ardl_mg
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


#' @title Variance-covariance for cs_ardl_mg
#' @description Diagonal VCOV based on MG standard errors.
#' @param object A \code{cs_ardl_mg} object.
#' @param type One of \code{"long_run"} (default) or \code{"short_run"}.
#' @param ... Not used.
#' @return A variance-covariance matrix with appropriate row/column names.
#' @export
#' @method vcov cs_ardl_mg
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
