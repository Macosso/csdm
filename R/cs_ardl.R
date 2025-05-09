#' Fit a Cross‐Sectionally Augmented ARDL Model
#'
#' Estimates an ARDL(p,p) with cross‐sectional averages (Pesaran 2006),
#' using Mean‐Group (MG) or Pooled Mean‐Group (PMG) pooling.
#'
#' @param formula A model formula `y ~ x1 + x2` for the dependent and regressors.
#' @param data A `data.frame` or `pdata.frame` containing the panel data.
#' @param p Integer lag order (non-negative) for both y and regressors.
#' @param method One of `"mg"` or `"pmg"`.
#' @param index Optional panel index (as in `plm`), e.g. `c("id","time")`.
#' @return An object of class `cs_ardl` with components:
#'   - `coefs`: list of coefficient vectors by group
#'   - `theta_i`: long‐run slopes per group
#'   - `theta`: pooled long‐run slope
#'   - `method`, `formula`, `call`, and `data_prepared`.
#' @export
cs_ardl <- function(formula, data, p = 1,
                    method = c("mg","pmg"), index = NULL) {
  method <- match.arg(method)
  # Convert to pdata.frame if necessary
  if (!inherits(data, "pdata.frame")) {
    if (is.null(index)) stop("Provide `index` or supply a pdata.frame")
    data <- plm::pdata.frame(data, index = index)
  }

  # Build model frame and extract response & regressors
  mf <- stats::model.frame(formula, data, na.action = na.fail)
  y <- stats::model.response(mf)
  x <- stats::model.matrix(stats::terms(mf), mf)[, -1, drop = FALSE]
  vars_needed <- c(colnames(x), "y")

  # Construct pdata.frame with only needed vars plus index
  idx <- names(attr(data, "index"))
  df_sub <- data[, idx, drop = FALSE]
  df_sub$y <- y
  df_sub <- cbind(df_sub, x)
  df_sub <- plm::pdata.frame(df_sub, index = idx)

  # Prepare lags and cross‐section means
  dat_prep <- .prepare_cs_ardl_data(df_sub, p)

  # Estimate and pool
  fit <- .fit_cs_ardl(dat_prep, p, method)
  fit$formula <- formula
  fit$call    <- match.call()
  class(fit)  <- c("cs_ardl","list")
  fit
}

#' Internal: Prepare Data for Cross‐Sectionally Augmented ARDL Using plm
#'
#' @inheritParams cs_ardl
#' @keywords internal
.prepare_cs_ardl_data <- function(panel_data, p) {
  if (!inherits(panel_data, "pdata.frame")) stop("panel_data must be pdata.frame")
  idx <- names(attr(panel_data, "index"))
  # Identify only y and regressors
  all_vars <- setdiff(names(panel_data), idx)

  # Compute cross‐sectional means by time
  cs_list <- lapply(all_vars, function(v) {
    ave(panel_data[[v]], panel_data[ , idx[2]], FUN = function(z) mean(z, na.rm = TRUE))
  })
  names(cs_list) <- paste0("cs_", all_vars)

  # Combine base data and cs means
  df <- cbind(as.data.frame(panel_data), as.data.frame(cs_list))
  df <- plm::pdata.frame(df)

  # Generate panel lags for each relevant var & its cs mean
  lag_vars <- c(all_vars, names(cs_list))
  if (p > 0) {
    for (k in seq_len(p)) {

      df <- df |>
        group_by(!!rlang::sym(idx[2])) |>
        mutate(across(all_of(lag_vars), ~ dplyr::lag(., k), .names = "{.col}.lag{k}")) |>
        ungroup()
    }
    # Retain only complete cases of lagged cols
    lag_cols <- grep("\\.lag[0-9]+$", names(df), value = TRUE)
    df <- df[stats::complete.cases(df[, lag_cols]), , drop = FALSE]
  }

  rownames(df) <- NULL
  df <- plm::pdata.frame(df, index = idx)
  df
}

#' Internal: Fit Cross‐Sectionally Augmented ARDL Models
#'
#' @inheritParams cs_ardl
#' @keywords internal
.fit_cs_ardl <- function(dat, p, method) {
  # Extract index and variable columns
  idx      <- names(attr(dat, "index"))
  y_col    <- "y"
  lag_cols <- grep("\\.lag[0-9]+$", names(dat), value = TRUE)
  reg_cols <- setdiff(names(dat), c(idx, y_col, lag_cols))

  groups <- as.character(unique(dat[[idx[1]]]))
  nG     <- length(groups)

  coefs_list <- setNames(vector("list", nG), groups)
  theta_i    <- setNames(numeric(nG), groups)
  resid_var  <- numeric(nG)
  nobs       <- numeric(nG)

  for (i in seq_along(groups)) {
    dat_i <- dat[dat[[idx[1]]] == groups[i], , drop = FALSE]
    Y     <- dat_i[[y_col]]
    X     <- as.matrix(dat_i[, c(reg_cols, lag_cols), drop = FALSE])
    X     <- cbind(Intercept = 1, X)

    fit_i <- stats::lm(Y ~ X - 1)
    b_i   <- stats::coef(fit_i)
    names(b_i) <- gsub("^X", "", names(b_i))

    coefs_list[[groups[i]]] <- b_i

    alpha_names <- grep(paste0("^", y_col, "\\.lag[1-9]"), names(b_i), value = TRUE)
    beta_names  <- setdiff(names(b_i), c("Intercept", alpha_names))
    theta_i[i]  <- -sum(b_i[beta_names]) / sum(b_i[alpha_names])

    e_i          <- stats::residuals(fit_i)
    resid_var[i] <- sum(e_i^2) / length(e_i)
    nobs[i]      <- length(e_i)
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


#' Print method for cs_ardl objects
#'
#' @param x A `cs_ardl` object.
#' @param ... Additional arguments (unused).
#' @export
print.cs_ardl <- function(x, ...) {
  cat("Cross‐Sectionally Augmented ARDL Model\n")
  cat("Method: ", x$method, "\n", sep = "")
  cat("Formula: ", deparse(x$formula), "\n", sep = "")
  cat("Long‐run coefficient (theta): ", round(x$theta, 4), "\n", sep = "")
  invisible(x)
}

#' Summary method for cs_ardl objects
#'
#' @param object A `cs_ardl` object.
#' @param ... Additional arguments (unused).
#' @return A summary list invisibly and prints details.
#' @export
summary.cs_ardl <- function(object, ...) {
  cat("Summary of CS‐ARDL Model\n")
  cat("-------------------------\n")
  cat("Method: ", object$method, "\n", sep = "")
  cat("Formula: ", deparse(object$formula), "\n\n", sep = "")

  cat("Pooled long‐run slope (theta): ", round(object$theta, 4), "\n\n", sep = "")
  cat("Individual long‐run slopes (theta_i):\n")
  print(round(object$theta_i, 4))
  cat("\nCoefficients by group: \n")
  lapply(names(object$coefs), function(g) {
    cat("Group ", g, ":\n", sep = "")
    print(round(object$coefs[[g]], 4))
    cat("\n")
  })

  invisible(list(theta = object$theta,
                 theta_i = object$theta_i,
                 coefs = object$coefs))
}
