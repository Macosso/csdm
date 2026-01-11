#' Common Correlated Effects Mean Group (CCE-MG) Estimator (static)
#'
#' @description
#' CCE-MG estimates unit-specific regressions augmented with cross-sectional averages
#' (CSAs) of the dependent variable and each regressor, then averages the slopes on the
#' original regressors across units (Mean Group). CSAs proxy latent common factors.
#'
#' @details
#' For each unit i, estimate:
#'   y_it = a_i + x_it' * beta_i + \eqn{\bar{y}_t} * delta_{y,i} + \eqn{\bar{x}_t'} * delta_{x,i} + e_it,
#' then report beta_MG = \eqn{N^{-1} \sum_i beta_i} and its MG standard errors based on the
#' cross-sectional dispersion of {beta_i}.
#'
#' This implementation supports unbalanced panels and an optional leave-one-out CSA.
#' The CSA coefficients (delta_y, delta_x) are treated as nuisances and not averaged.
#'
#' @param formula Model formula (e.g., y ~ x1 + x2). Intercept is included by default.
#' @param data data.frame or plm::pdata.frame.
#' @param id,time Column names for unit and time (required if data is not a pdata.frame).
#' @param leave_out logical; if TRUE, compute CSAs excluding the unit itself. Default FALSE.
#' @param na.action function; default stats::na.omit.
#' @return An object of class "cce_mean_group" with components:
#'   - beta_mg, beta_se, beta_z, beta_p: MG point estimates, SEs, and tests
#'   - beta_i: matrix of unit-specific slopes for original regressors (rows=ids)
#'   - gamma_i: matrix of unit-specific nuisance CSA slopes (csa_y and csa_xk)
#'   - residuals_cce, residuals_pca, residuals_pca_std: N x T residual matrices
#'   - sigma_i: unit std dev used for standardization
#'   - vcov: full covariance matrix for beta_mg
#'   - call, formula, ids, times, N, K, leave_out
#'
#' @importFrom plm pdata.frame index
#' @importFrom stats model.frame model.response model.matrix lm coef pnorm var
#' @export
cce_mean_group <- function(formula, data, id = NULL, time = NULL,
                           leave_out = FALSE, na.action = stats::na.omit) {

  .Deprecated("csdm", package = "csdm")

  # indices
  if (!inherits(data, "pdata.frame")) {
    if (is.null(id) || is.null(time)) {
      stop("If data is not a pdata.frame, both 'id' and 'time' must be provided.")
    }
    data <- plm::pdata.frame(data, index = c(id, time))
  } else {
    idx <- attr(data, "index")
    id   <- names(idx)[1]
    time <- names(idx)[2]
  }

  data[[id]] <- as.character(data[[id]])
  df <- na.action(as.data.frame(data))

  # model parts
  mf <- stats::model.frame(formula, df, na.action = na.action)
  y  <- stats::model.response(mf)
  X  <- stats::model.matrix(formula, mf)
  Xn <- if ("(Intercept)" %in% colnames(X)) X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE] else X

  # base panel table
  tab <- cbind.data.frame(
    .id   = df[[id]],
    .time = df[[time]],
    .y    = as.numeric(y)
  )
  if (ncol(Xn) > 0) tab <- cbind(tab, as.data.frame(Xn))

  # --- CSAs via utils_avg::cross_sectional_avg (always attached) ---
  vars <- c(".y", colnames(Xn))
  # avoid leading-dot name inside util
  tab_tmp <- stats::setNames(tab, sub("^\\.y$", "y__tmp__", names(tab)))
  ca <- cross_sectional_avg(
    data = tab_tmp, id = ".id", time = ".time",
    vars = sub("^\\.y$", "y__tmp__", vars),
    leave_out = leave_out, suffix = "csa", return_mode = "attach", na.rm = TRUE
  )
  # map back to expected names
  ca <- stats::setNames(ca, sub("^csa_y__tmp__$", "csa_.y", names(ca)))
  tab$csa_y <- ca[["csa_.y"]]
  xnames <- colnames(Xn)
  if (length(xnames)) {
    for (xn in xnames) tab[[paste0("csa_", xn)]] <- ca[[paste0("csa_", xn)]]
  }

  # --- per-id CCE OLS ---
  ids  <- unique(tab$.id)
  K    <- ncol(Xn)
  beta_i  <- matrix(NA_real_, nrow = length(ids), ncol = K,
                    dimnames = list(as.character(ids), xnames))
  gam_cols <- c("csa_y", paste0("csa_", xnames))
  gamma_i <- matrix(NA_real_, nrow = length(ids), ncol = length(gam_cols),
                    dimnames = list(as.character(ids), gam_cols))

  # ... keep everything up to the creation of beta_i and gamma_i

  ids  <- unique(tab$.id)
  K    <- ncol(Xn)
  beta_i  <- matrix(NA_real_, nrow = length(ids), ncol = K,
                    dimnames = list(as.character(ids), colnames(Xn)))
  gam_cols <- c("csa_y", paste0("csa_", colnames(Xn)))
  gamma_i <- matrix(NA_real_, nrow = length(ids), ncol = length(gam_cols),
                    dimnames = list(as.character(ids), gam_cols))

  # NEW: containers for R2
  r2_i      <- rep(NA_real_, length(ids))
  adj_r2_i  <- rep(NA_real_, length(ids))
  ssr_i     <- rep(NA_real_, length(ids))
  tss_i     <- rep(NA_real_, length(ids))
  dof_i     <- rep(NA_integer_, length(ids))

  for (j in seq_along(ids)) {
    ii <- ids[j]
    sub <- tab[tab$.id == ii, , drop = FALSE]

    rhs_terms <- c(colnames(Xn), "csa_y", paste0("csa_", colnames(Xn)))
    rhs_terms <- rhs_terms[rhs_terms %in% colnames(sub)]
    fml <- stats::as.formula(paste(".y ~", paste(rhs_terms, collapse = " + ")))

    fit <- stats::lm(fml, data = sub)
    cf  <- stats::coef(fit)

    # save original-regressor slopes
    if (K > 0) {
      keep_b <- intersect(colnames(Xn), names(cf))
      beta_i[j, keep_b] <- cf[keep_b]
    }
    # save CSA slopes
    keep_g <- intersect(gam_cols, names(cf))
    if (length(keep_g)) gamma_i[j, keep_g] <- cf[keep_g]

    # --- NEW: R-squared bits (computed on the actual estimation sample)
    mfi   <- stats::model.frame(fit)
    yi    <- stats::model.response(mfi)
    ei    <- fit$residuals
    Ti    <- length(yi)
    pi    <- sum(!is.na(cf))                       # #estimated params incl. intercept
    ssr   <- sum(ei^2)
    tss   <- sum((yi - mean(yi))^2)

    ssr_i[j]    <- ssr
    tss_i[j]    <- tss
    dof_i[j]    <- Ti - pi
    r2_i[j]     <- if (tss > 0) 1 - ssr/tss else NA_real_
    adj_r2_i[j] <- if (tss > 0 && (Ti - 1) > 0 && dof_i[j] > 0)
      1 - (ssr/dof_i[j]) / (tss/(Ti - 1)) else NA_real_
  }

  # MG aggregation for original regressors (unchanged) ...
  # ...

  # --- NEW: Aggregate R2 measures
  mg_r2_mean      <- mean(r2_i, na.rm = TRUE)                       # "R-squared (mg)" like xtdcce2
  mg_r2_overall   <- 1 - sum(ssr_i, na.rm = TRUE)/sum(tss_i, na.rm = TRUE)
  mg_adj_r2_mean  <- mean(adj_r2_i, na.rm = TRUE)

  out <- list(
    call      = match.call(),
    formula   = formula,
    ids       = ids,
    times     = sort(unique(tab$.time)),
    N         = length(ids),
    K         = K,
    beta_mg   = beta_bar,
    beta_se   = beta_se,
    beta_z    = beta_z,
    beta_p    = beta_p,
    beta_i    = beta_i,
    gamma_i   = gamma_i,
    leave_out = leave_out,

    # NEW: return R2 info
    r2 = list(
      per_unit     = stats::setNames(r2_i,  ids),
      adj_per_unit = stats::setNames(adj_r2_i, ids),
      mg           = mg_r2_mean,      # this is what we'll print by default
      mg_overall   = mg_r2_overall,
      mg_adj       = mg_adj_r2_mean,
      ssr          = stats::setNames(ssr_i, ids),
      tss          = stats::setNames(tss_i, ids),
      dof          = stats::setNames(dof_i, ids)
    )
  )
  class(out) <- "cce_mean_group"
  out
}

#' @export
print.cce_mean_group <- function(x, ...) {
  cat("\nCommon Correlated Effects - Mean Group (CCE-MG)\n")
  if (length(x$beta_mg)) {
    res <- data.frame(Estimate = x$beta_mg, Std.Error = x$beta_se,
                      z.value = x$beta_z, p.value = x$beta_p)
    print(res)
  } else {
    cat("No regressors to report (K = 0).\n")
  }
}

#' @export
summary.cce_mean_group <- function(object, ...) {
  print(object)
  cat("\nNotes:\n",
      "- Slopes shown are MG averages of unit-specific coefficients on ORIGINAL regressors only.\n",
      "- CSA coefficients are nuisances; see $gamma_i.\n",
      "- Residual matrices for CD/CD*: $residuals_cce, $residuals_pca, $residuals_pca_std.\n", sep = "")
}

#' @export
coef.cce_mean_group <- function(object, ...) {
  object$beta_mg
}

#' @export
vcov.cce_mean_group <- function(object, ...) {
  object$vcov
}

#' @export
residuals.cce_mean_group <- function(object, type = c("auto","pca_std","pca","cce"), ...) {
  # always use utils_residuals accessor (attached)
  get_residuals(object, type = match.arg(type), strict = TRUE)
}
