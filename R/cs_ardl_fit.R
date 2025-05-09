#' Internal: Prepare Data for Cross‐Sectionally Augmented ARDL Using plm
#'
#' This function converts input to a `pdata.frame`, computes cross‐sectional means,
#' creates panel lags using `plm::lag`, and drops initial rows lacking full lags.
#'
#' @param panel_data A `pdata.frame` or a data.frame with panel structure.
#'   Must contain the dependent variable `y` and regressors.
#' @param p Integer lag order (non-negative).
#' @return A `data.frame` (not `pdata.frame`) with:
#'   - original index columns (`<idx>`)
#'   - `y`, regressors
#'   - `cs_<var>`: cross‐sectional means at each time
#'   - `<var>.lag<k>`: panel lags for `y`, regs, and `cs_` vars for k = 1..p
#' Only complete cases after lagging are returned.
#' @keywords internal
.prepare_cs_ardl_data <- function(panel_data, p) {
  # Require plm
  if (!inherits(panel_data, "pdata.frame")) {
    stopifnot(is.data.frame(panel_data))
    # user must have .index_names attribute or specify via attr
    idx <- attr(panel_data, "index")
    if (is.null(idx) || length(idx) != 2) {
      stop("panel_data must be a pdata.frame or have 'index' attribute of length 2")
    }
    panel_data <- plm::pdata.frame(panel_data, index = idx)
  }

  # Extract index names
  idx_names <- names(attr(panel_data, "index"))
  group_idx <- idx_names[1]
  time_idx  <- idx_names[2]

  # Variables to be used
  all_vars <- setdiff(names(panel_data), idx_names)

  # Compute cross‐sectional means
  cs_list <- lapply(all_vars, function(v) {
    cs <- ave(panel_data[[v]], panel_data[, time_idx], FUN = mean, na.rm = TRUE)
    cs
  })
  names(cs_list) <- paste0("cs_", all_vars)

  # Combine original and cs vars
  df <- cbind(as.data.frame(panel_data), as.data.frame(cs_list))

  # Generate panel lags for each var (including cs_ vars)
  lag_vars <- c(all_vars, names(cs_list))
  for (k in seq_len(p)) {
    lagged <- plm::lag(as.data.frame(df[lag_vars]), k)
    names(lagged) <- paste0(names(lagged), ".lag", k)
    df <- cbind(df, lagged)
  }

  # Drop rows with any NA in lagged columns
  if (p > 0) {
    lag_cols <- grep("\\.lag[0-9]+$", names(df), value = TRUE)
    df <- df[stats::complete.cases(df[, lag_cols]), ]
  }

  # Return a plain data.frame for estimation
  rownames(df) <- NULL
  df
}

#' Internal: Fit Cross‐Sectionally Augmented ARDL Models
#'
#' Assumes data prepared by `.prepare_cs_ardl_data()`, i.e. a plain `data.frame` with
#' correct lagged variables. Uses the same pooling logic as before.
#' @inherit .fit_cs_ardl params
#' @keywords internal
.fit_cs_ardl <- function(dat, p, method) {
  stopifnot(is.data.frame(dat), is.numeric(p), p >= 0,
            method %in% c("mg", "pmg"))

  # Identify index columns
  idx_cols <- grep("^(?!cs_|y$|\\.|[[:alnum:]_]+\\.lag)$(?!Intercept)", names(dat), perl=TRUE, value=TRUE)
  # Better: pick y and all columns ending in .lag or starting with cs_
  y_col <- "y"
  lag_cols <- grep("\\.lag[0-9]+$", names(dat), value=TRUE)
  reg_cols <- setdiff(names(dat), c(idx_cols, y_col, lag_cols))

  groups <- unique(dat[[idx_cols[1]]])
  nG     <- length(groups)

  coefs_list <- vector("list", nG)
  names(coefs_list) <- groups
  theta_i    <- setNames(numeric(nG), groups)
  resid_var  <- numeric(nG)
  nobs       <- numeric(nG)

  for (i in seq_along(groups)) {
    group_val <- groups[i]
    dat_i     <- dat[dat[[idx_cols[1]]] == group_val, , drop = FALSE]

    Y <- dat_i[[y_col]]
    X <- as.matrix(dat_i[, setdiff(names(dat_i), c(idx_cols, y_col)), drop = FALSE])
    X <- cbind(Intercept = 1, X)

    fit_i <- stats::lm(Y ~ X - 1)
    b_i   <- stats::coef(fit_i)
    coefs_list[[i]] <- b_i

    alpha_names <- grep("^y\\.lag", names(b_i), value=TRUE)
    beta_names  <- setdiff(names(b_i), c("Intercept", alpha_names))

    sum_alpha <- sum(b_i[alpha_names])
    sum_beta  <- sum(b_i[beta_names])
    theta_i[i] <- -sum_beta / sum_alpha

    e_i        <- stats::residuals(fit_i)
    resid_var[i] <- sum(e_i^2) / length(e_i)
    nobs[i]    <- length(e_i)
  }

  if (method == "mg") {
    theta <- mean(theta_i)
  } else {
    w     <- nobs / resid_var
    theta <- sum(w * theta_i) / sum(w)
  }

  list(coefs         = coefs_list,
       theta_i       = theta_i,
       theta         = theta,
       method        = method,
       data_prepared = dat)
}
