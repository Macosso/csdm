# csdm_methods.R

#' @export
print.csdm_fit <- function(x, digits = 4, ...) {
  cat("csdm fit (", x$model, ")\n", sep = "")
  if (!is.null(x$formula)) cat("Formula: ", deparse(x$formula), "\n", sep = "")
  if (!is.null(x$meta$N) && !is.null(x$meta$T)) {
    cat("N = ", x$meta$N, ", T = ", x$meta$T, "\n", sep = "")
  }

  if (length(x$coef_mg)) {
    est <- x$coef_mg
    se  <- x$se_mg
    tab <- data.frame(
      Estimate = as.numeric(est),
      Std.Error = as.numeric(se[names(est)]),
      check.names = FALSE
    )
    rownames(tab) <- names(est)
    print(round(tab, digits))
  }
  invisible(x)
}


#' @export
summary.csdm_fit <- function(object, digits = 4, ...) {
  make_table <- function(estimate, se, n_used = NULL) {
    estimate <- as.numeric(estimate)
    se <- as.numeric(se)
    z <- estimate / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    crit <- stats::qnorm(0.975)
    ci_lo <- estimate - crit * se
    ci_hi <- estimate + crit * se

    # handle missing/invalid SE
    bad <- !is.finite(se) | se <= 0
    if (any(bad)) {
      z[bad] <- NA_real_
      p[bad] <- NA_real_
      ci_lo[bad] <- NA_real_
      ci_hi[bad] <- NA_real_
    }

    tab <- data.frame(
      `Coef.` = estimate,
      `Std. Err.` = se,
      z = z,
      `P>|z|` = p,
      `CI 2.5%` = ci_lo,
      `CI 97.5%` = ci_hi,
      check.names = FALSE
    )
    if (!is.null(n_used)) {
      tab$n_used <- as.integer(n_used)
    }
    tab
  }

  est <- object$coef_mg
  se  <- object$se_mg
  z   <- est / se[names(est)]
  p   <- 2 * (1 - stats::pnorm(abs(z)))

  coef_table <- data.frame(
    term = names(est),
    estimate = as.numeric(est),
    se = as.numeric(se[names(est)]),
    z = as.numeric(z),
    p_value = as.numeric(p),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  out <- list(
    call = object$call,
    formula = object$formula,
    model = object$model,
    N = object$meta$N,
    T = object$meta$T,
    coef_table = coef_table
  )

  if (identical(object$model, "cs_ardl") && !is.null(object$cs_ardl) && !is.null(object$cs_ardl$mg)) {
    nobs <- NA_integer_
    R2_mg <- NA_real_
    if (!is.null(object$stats)) {
      if (!is.null(object$stats$nobs)) nobs <- as.integer(object$stats$nobs)
      if (!is.null(object$stats$R2_mg)) R2_mg <- as.numeric(object$stats$R2_mg)
    }
    out$nobs <- nobs
    out$stats <- list(R2_mg = R2_mg)

    short_run_tab <- make_table(
      estimate = as.numeric(est),
      se = as.numeric(se[names(est)])
    )
    rownames(short_run_tab) <- names(est)

    adj_term <- paste0("lr_", object$cs_ardl$y)
    adj_est <- as.numeric(object$cs_ardl$mg$adj[["estimate"]])
    adj_se <- as.numeric(object$cs_ardl$mg$adj[["se"]])
    adjust_tab <- make_table(estimate = adj_est, se = adj_se)
    rownames(adjust_tab) <- adj_term

    lr_tab_raw <- object$cs_ardl$mg$lr
    if (is.null(lr_tab_raw) || !nrow(lr_tab_raw)) {
      long_run_tab <- data.frame(
        `Coef.` = numeric(0),
        `Std. Err.` = numeric(0),
        z = numeric(0),
        `P>|z|` = numeric(0),
        `CI 2.5%` = numeric(0),
        `CI 97.5%` = numeric(0),
        n_used = integer(0),
        check.names = FALSE
      )
    } else {
      long_run_tab <- make_table(
        estimate = lr_tab_raw$estimate,
        se = lr_tab_raw$se,
        n_used = lr_tab_raw$n_used
      )
      rownames(long_run_tab) <- lr_tab_raw$term
    }

    # Variable lists (for printing + testing)
    yname <- object$cs_ardl$y
    xnames <- object$cs_ardl$x
    p_ylags <- as.integer(object$cs_ardl$ylags)
    q_xdlags <- as.integer(object$cs_ardl$xdlags)

    mg_vars <- c(
      if (is.finite(p_ylags) && p_ylags > 0L) paste0("lag", seq_len(p_ylags), "_", yname) else character(0),
      xnames,
      if (length(xnames) && is.finite(q_xdlags) && q_xdlags > 0L)
        unlist(lapply(xnames, function(xn) paste0("lag", seq_len(q_xdlags), "_", xn)), use.names = FALSE)
      else character(0)
    )

    csa_spec <- object$meta$csa
    csa_vars_txt <- "_none"
    csa_lags_txt <- "0"
    if (!is.null(csa_spec)) {
      csa_vars_txt <- if (identical(csa_spec$vars, "_all") || identical(csa_spec$vars, "_none")) {
        as.character(csa_spec$vars)
      } else {
        paste(as.character(csa_spec$vars), collapse = ", ")
      }
      csa_lags_txt <- if (length(csa_spec$lags) == 1L) {
        as.character(as.integer(csa_spec$lags))
      } else {
        "(per-variable)"
      }
    }

    out$tables <- list(
      short_run = short_run_tab,
      adjust_term = adjust_tab,
      long_run = long_run_tab
    )
    out$lists <- list(
      mean_group_variables = mg_vars,
      csa_vars = csa_vars_txt,
      csa_lags = csa_lags_txt,
      long_run_variables = xnames,
      cointegration_variables = yname
    )
  }

  class(out) <- "summary.csdm_fit"
  out
}


#' @export
print.summary.csdm_fit <- function(x, digits = 4, ...) {
  cat("csdm summary (", x$model, ")\n", sep = "")
  if (!is.null(x$formula)) cat("Formula: ", deparse(x$formula), "\n", sep = "")

  if (identical(x$model, "cs_ardl") && !is.null(x$tables) && !is.null(x$stats) && !is.null(x$lists)) {
    if (!is.null(x$N) && !is.null(x$T)) {
      cat("N = ", x$N, ", T = ", x$T, "\n", sep = "")
    }
    if (!is.null(x$nobs)) {
      cat("Number of obs = ", x$nobs, "\n", sep = "")
    }
    if (!is.null(x$stats$R2_mg)) {
      cat("R-squared (mg) = ", round(x$stats$R2_mg, digits), "\n\n", sep = "")
    } else {
      cat("\n")
    }

    cat("Short Run Est.\n")
    tab <- x$tables$short_run
    tab[] <- lapply(tab, function(col) if (is.numeric(col)) round(col, digits) else col)
    print(tab)

    cat("\nAdjust. Term\n")
    atab <- x$tables$adjust_term
    atab[] <- lapply(atab, function(col) if (is.numeric(col)) round(col, digits) else col)
    print(atab)

    cat("\nLong Run Est.\n")
    lrtab <- x$tables$long_run
    lrtab[] <- lapply(lrtab, function(col) if (is.numeric(col)) round(col, digits) else col)
    print(lrtab)

    cat("\n")
    cat("Mean Group Variables: ", paste(x$lists$mean_group_variables, collapse = ", "), "\n", sep = "")
    cat("Cross Sectional Averaged Variables: ", x$lists$csa_vars, " (lags=", x$lists$csa_lags, ")\n", sep = "")
    cat("Long Run Variables: ", paste(x$lists$long_run_variables, collapse = ", "), "\n", sep = "")
    cat("Cointegration variable(s): ", x$lists$cointegration_variables, "\n", sep = "")
  } else {
    if (!is.null(x$N) && !is.null(x$T)) {
      cat("N = ", x$N, ", T = ", x$T, "\n\n", sep = "")
    }
    tab <- x$coef_table
    num_cols <- intersect(c("estimate", "se", "z", "p_value"), names(tab))
    tab[num_cols] <- lapply(tab[num_cols], function(col) round(col, digits))
    print(tab, row.names = FALSE)
  }

  invisible(x)
}


#' @export
coef.csdm_fit <- function(object, ...) {
  if (identical(object$model, "cs_ardl") && !is.null(object$cs_ardl) && !is.null(object$cs_ardl$mg)) {
    lr_y <- setNames(
      as.numeric(object$cs_ardl$mg$adj[["estimate"]]),
      paste0("lr_", object$cs_ardl$y)
    )
    lr_x <- numeric(0)
    if (!is.null(object$cs_ardl$mg$lr) && nrow(object$cs_ardl$mg$lr)) {
      lr_x <- stats::setNames(
        as.numeric(object$cs_ardl$mg$lr$estimate),
        as.character(object$cs_ardl$mg$lr$term)
      )
    }
    return(c(object$coef_mg, lr_y, lr_x))
  }
  object$coef_mg
}


#' @export
vcov.csdm_fit <- function(object, ...) {
  object$vcov_mg
}


#' @export
residuals.csdm_fit <- function(object, type = c("e", "u"), ...) {
  type <- match.arg(type)
  if (type == "u") stop("residuals(type='u') not implemented yet")
  object$residuals_e
}


#' @export
predict.csdm_fit <- function(object, newdata = NULL, type = c("xb", "residuals"), ...) {
  type <- match.arg(type)
  if (!is.null(newdata)) stop("predict(newdata=...) not implemented yet")
  if (type == "residuals") return(stats::residuals(object, type = "e"))
  if (!is.null(object$fitted_xb)) return(object$fitted_xb)
  stop("predict(type='xb') not implemented yet")
}
