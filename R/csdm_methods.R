# csdm_methods.R

#' Compact print method for fitted csdm models
#'
#' Prints a concise overview of a fitted \code{csdm_fit} object, including the
#' model type, formula, panel dimensions, and a coefficient table with standard
#' errors when available.
#'
#' @param x A fitted object of class \code{csdm_fit}.
#' @param digits Number of printed digits.
#' @param ... Currently unused.
#'
#' @return Invisibly returns \code{x}.
#' @seealso [summary.csdm_fit()], [coef.csdm_fit()], [residuals.csdm_fit()]
#' @export
print.csdm_fit <- function(x, digits = 4, ...) {
  cat("csdm fit (", x$model, ")\n", sep = "")
  if (!is.null(x$formula)) cat("Formula: ", deparse(x$formula), "\n", sep = "")
  if (!is.null(x$meta$N) && !is.null(x$meta$T)) {
    cat("N: ", x$meta$N, ", T: ", x$meta$T, "\n", sep = "")
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

#' Summarize csdm model estimation results
#'
#' Computes post-estimation summaries for \code{csdm_fit} objects, including
#' mean-group coefficient inference, model-level diagnostics, and model-specific
#' summary tables (for example, short-run and long-run blocks for CS-ARDL).
#'
#' @param object A fitted model object of class \code{csdm_fit}.
#' @param digits Number of digits to print.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' ## Reported inference
#'
#' For each coefficient \eqn{\hat\beta_k}, the summary reports standard errors,
#' \eqn{z}-statistics, and two-sided normal-approximation p-values:
#' \deqn{z_k = \frac{\hat\beta_k}{\operatorname{se}(\hat\beta_k)}, \qquad
#' p_k = 2\{1-\Phi(|z_k|)\}.}
#'
#' ## Diagnostics
#'
#' The printed summary shows the classic Pesaran CD diagnostic by default. Extended
#' diagnostics (CDw, CDw+, CD*) are available through [cd_test()].
#'
#' @return An object of class \code{summary.csdm_fit} with core metadata
#'   (call/formula/model/N/T), coefficient tables, fit statistics, and
#'   model-specific components for printing and downstream inspection.
#' @seealso [print.summary.csdm_fit()], [cd_test()], [coef.csdm_fit()], [vcov.csdm_fit()]
#'
#' @examples
#' data(PWT_60_07, package = "csdm")
#' df <- PWT_60_07
#' ids <- unique(df$id)[1:10]
#' df_small <- df[df$id %in% ids & df$year >= 1970, ]
#' fit <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small,
#'   id = "id",
#'   time = "year",
#'   model = "cce",
#'   csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"))
#' )
#' s <- summary(fit)
#' s
#' @export
summary.csdm_fit <- function(object, digits = 4, ...) {
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
    nobs <- if (!is.null(object$stats$nobs)) as.integer(object$stats$nobs) else NA_integer_
    out$nobs <- nobs
    out$stats <- list(
      R2_mg = if (!is.null(object$stats$R2_mg)) as.numeric(object$stats$R2_mg) else NA_real_,
      R2_ols_mg = if (!is.null(object$stats$R2_ols_mg)) as.numeric(object$stats$R2_ols_mg) else NA_real_,
      CD_stat = if (!is.null(object$stats$CD_stat)) as.numeric(object$stats$CD_stat) else NA_real_,
      CD_p = if (!is.null(object$stats$CD_p)) as.numeric(object$stats$CD_p) else NA_real_,
      CDw_stat = if (!is.null(object$stats$CDw_stat)) as.numeric(object$stats$CDw_stat) else NA_real_,
      CDw_p = if (!is.null(object$stats$CDw_p)) as.numeric(object$stats$CDw_p) else NA_real_,
      CDw_plus_stat = if (!is.null(object$stats$CDw_plus_stat)) as.numeric(object$stats$CDw_plus_stat) else NA_real_,
      CDw_plus_p = if (!is.null(object$stats$CDw_plus_p)) as.numeric(object$stats$CDw_plus_p) else NA_real_,
      CDstar_stat = if (!is.null(object$stats$CDstar_stat)) as.numeric(object$stats$CDstar_stat) else NA_real_,
      CDstar_p = if (!is.null(object$stats$CDstar_p)) as.numeric(object$stats$CDstar_p) else NA_real_,
      CDstar_n_pc = if (!is.null(object$stats$CDstar_n_pc)) as.integer(object$stats$CDstar_n_pc) else NA_integer_
    )

    short_run_tab <- .csdm_make_coef_table(
      estimate = as.numeric(est),
      se = as.numeric(se[names(est)])
    )
    rownames(short_run_tab) <- names(est)

    adj_term <- paste0("lr_", object$cs_ardl$y)
    adj_est <- as.numeric(object$cs_ardl$mg$adj[["estimate"]])
    adj_se <- as.numeric(object$cs_ardl$mg$adj[["se"]])
    adjust_tab <- .csdm_make_coef_table(estimate = adj_est, se = adj_se)
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
      long_run_tab <- .csdm_make_coef_table(
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
  } else {
    # Common summary structure for MG / CCE / DCCE
    nobs <- if (!is.null(object$stats$nobs)) as.integer(object$stats$nobs) else NA_integer_
    out$nobs <- nobs
    out$stats <- list(
      R2_mg = if (!is.null(object$stats$R2_mg)) as.numeric(object$stats$R2_mg) else NA_real_,
      CD_stat = if (!is.null(object$stats$CD_stat)) as.numeric(object$stats$CD_stat) else NA_real_,
      CD_p = if (!is.null(object$stats$CD_p)) as.numeric(object$stats$CD_p) else NA_real_,
      CDw_stat = if (!is.null(object$stats$CDw_stat)) as.numeric(object$stats$CDw_stat) else NA_real_,
      CDw_p = if (!is.null(object$stats$CDw_p)) as.numeric(object$stats$CDw_p) else NA_real_,
      CDw_plus_stat = if (!is.null(object$stats$CDw_plus_stat)) as.numeric(object$stats$CDw_plus_stat) else NA_real_,
      CDw_plus_p = if (!is.null(object$stats$CDw_plus_p)) as.numeric(object$stats$CDw_plus_p) else NA_real_,
      CDstar_stat = if (!is.null(object$stats$CDstar_stat)) as.numeric(object$stats$CDstar_stat) else NA_real_,
      CDstar_p = if (!is.null(object$stats$CDstar_p)) as.numeric(object$stats$CDstar_p) else NA_real_,
      CDstar_n_pc = if (!is.null(object$stats$CDstar_n_pc)) as.integer(object$stats$CDstar_n_pc) else NA_integer_
    )

    mg_tab <- .csdm_make_coef_table(
      estimate = as.numeric(est),
      se = as.numeric(se[names(est)])
    )
    rownames(mg_tab) <- names(est)
    out$tables <- list(mean_group = mg_tab)

    mg_vars <- setdiff(names(est), c("(Intercept)"))
    csa_footer <- .csdm_format_csa_footer(object$meta$csa)
    out$lists <- list(
      mean_group_variables = mg_vars,
      csa_vars = csa_footer$csa_vars,
      csa_lags = csa_footer$csa_lags
    )
  }

  class(out) <- "summary.csdm_fit"
  out
}


#' Print method for csdm summary objects
#'
#' Formats and prints a \code{summary.csdm_fit} object. Output adapts to model
#' type and includes coefficient tables, selected goodness-of-fit diagnostics,
#' and compact model metadata.
#'
#' @param x A \code{summary.csdm_fit} object.
#' @param digits Number of digits to print.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' The printout includes classic Pesaran CD diagnostics from the summary object.
#' For a full CD diagnostic panel (CD, CDw, CDw+, CD*), use [cd_test()] on the
#' fitted model.
#'
#' @return Invisibly returns \code{x}.
#' @seealso [summary.csdm_fit()], [cd_test()]
#'
#' @export
print.summary.csdm_fit <- function(x, digits = 4, ...) {

  model_amapping <- c(
    "mg" = "Mean Group Model (MG)",
    "cce" = "Static Common Correlated Error Model (CCE)",
    "dcce" = "Dynamic Common Correlated Error Model (DCCE)",
    "cs_ardl" = "Cross-Sectional ARDL (CS-ARDL)",
    "cs_ecm" = "Cross-Sectional ECM (CS-ECM)",
    "cs_dl" = "Cross-Sectional Distributed Lag (CS-DL)"
    )

  cat("csdm summary: ", model_amapping[x$model], "\n", sep = "")
  if (!is.null(x$formula)) cat("Formula: ", deparse(x$formula), "\n", sep = "")

  signif_footer_printed <- FALSE
  if (identical(x$model, "cs_ardl") && !is.null(x$tables) && !is.null(x$stats) && !is.null(x$lists)) {
    if (!is.null(x$N) && !is.null(x$T)) {
      cat("N: ", x$N, ", T: ", x$T, "\n", sep = "")
    }
    if (!is.null(x$nobs)) {
      cat("Number of obs: ", x$nobs, "\n", sep = "")
    }
    if (!is.null(x$stats$R2_mg)) {
      cat("R-squared (mg): ", round(x$stats$R2_mg, digits), "\n\n", sep = "")
    }
    # Print CD test (classic only)
    .csdm_print_cd_tests(x$stats, digits, classic_only = TRUE)

    cat("Short Run Est.\n")
    tab <- x$tables$short_run
    tab[] <- lapply(tab, function(col) if (is.numeric(col)) round(col, digits) else col)
    print(tab)
    if ("Signif." %in% colnames(tab)) signif_footer_printed <- TRUE

    cat("\nAdjust. Term\n")
    atab <- x$tables$adjust_term
    atab[] <- lapply(atab, function(col) if (is.numeric(col)) round(col, digits) else col)
    print(atab)
    if ("Signif." %in% colnames(atab)) signif_footer_printed <- TRUE

    cat("\nLong Run Est.\n")
    lrtab <- x$tables$long_run
    lrtab[] <- lapply(lrtab, function(col) if (is.numeric(col)) round(col, digits) else col)
    print(lrtab)
    if ("Signif." %in% colnames(lrtab)) signif_footer_printed <- TRUE

    cat("\n")
    cat("Mean Group Variables: ", paste(x$lists$mean_group_variables, collapse = ", "), "\n", sep = "")
    cat("Cross Sectional Averaged Variables: ", x$lists$csa_vars, " (lags=", x$lists$csa_lags, ")\n", sep = "")
    cat("Long Run Variables: ", paste(x$lists$long_run_variables, collapse = ", "), "\n", sep = "")
    cat("Cointegration variable(s): ", x$lists$cointegration_variables, "\n", sep = "")
  } else if (!is.null(x$tables) && !is.null(x$stats) && !is.null(x$lists) && !is.null(x$tables$mean_group)) {
    if (!is.null(x$N) && !is.null(x$T)) {
      cat("N: ", x$N, ", T: ", x$T, "\n", sep = "")
    }
    if (!is.null(x$nobs)) {
      cat("Number of obs: ", x$nobs, "\n", sep = "")
    }
    if (!is.null(x$stats$R2_mg)) {
      cat("R-squared (mg): ", round(x$stats$R2_mg, digits), "\n", sep = "")
      # Print CD test (classic only)
      .csdm_print_cd_tests(x$stats, digits, classic_only = TRUE)
    }

    cat("Mean Group:\n")
    tab <- x$tables$mean_group
    tab[] <- lapply(tab, function(col) if (is.numeric(col)) round(col, digits) else col)
    print(tab)
    if ("Signif." %in% colnames(tab)) signif_footer_printed <- TRUE

    cat("\n")
    cat("Mean Group Variables: ", paste(x$lists$mean_group_variables, collapse = ", "), "\n", sep = "")
    cat("Cross Sectional Averaged Variables: ", x$lists$csa_vars, " (lags=", x$lists$csa_lags, ")\n", sep = "")
  } else {
    if (!is.null(x$N) && !is.null(x$T)) {
      cat("N: ", x$N, ", T: ", x$T, "\n\n", sep = "")
    }
    tab <- x$coef_table
    num_cols <- intersect(c("estimate", "se", "z", "p_value"), names(tab))
    tab[num_cols] <- lapply(tab[num_cols], function(col) round(col, digits))
    print(tab, row.names = FALSE)
    if ("Signif." %in% colnames(tab)) signif_footer_printed <- TRUE
  }

  if (signif_footer_printed) {
    cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  invisible(x)
}

#' Extract model coefficients from a fitted csdm model
#'
#' Returns estimated mean-group coefficients from a \code{csdm_fit} object. For
#' \code{model = "cs_ardl"}, the returned vector includes short-run mean-group
#' coefficients, the adjustment coefficient (named \code{lr_<y>}), and long-run
#' coefficients when available.
#'
#' @param object A fitted object of class \code{csdm_fit}.
#' @param ... Currently unused.
#'
#' @return A named numeric vector of estimated coefficients.
#' @seealso [summary.csdm_fit()], [vcov.csdm_fit()]
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


#' Extract coefficient covariance matrix from a fitted csdm model
#'
#' @param object A fitted object of class \code{csdm_fit}.
#' @param ... Currently unused.
#'
#' @return A numeric variance-covariance matrix aligned with \code{coef(object)}
#'   for models where this is available.
#' @seealso [coef.csdm_fit()], [summary.csdm_fit()]
#' @export
vcov.csdm_fit <- function(object, ...) {
  object$vcov_mg
}


#' Extract residual matrix from a fitted csdm model
#'
#' Returns residuals as an \eqn{N x T} matrix (rows are units, columns are time).
#' This method is designed for panel diagnostics and downstream tools such as
#' [cd_test()].
#'
#' @param object A fitted object of class \code{csdm_fit}.
#' @param type Residual type. Currently only \code{"e"} is implemented.
#' @param ... Currently unused.
#'
#' @return A numeric matrix of residuals with dimensions \eqn{N x T}.
#' @seealso [get_residuals()], [cd_test()], [predict.csdm_fit()]
#' @export
residuals.csdm_fit <- function(object, type = c("e", "u"), ...) {
  type <- match.arg(type)
  if (type == "u") stop("residuals(type='u') not implemented yet")
  object$residuals_e
}


#' Predict method for csdm models
#'
#' Produces fitted values (index \code{"xb"}) when available, or returns model
#' residuals. Prediction on new data is not yet implemented.
#'
#' @param object A fitted object of class \code{csdm_fit}.
#' @param newdata Optional new data (not yet supported).
#' @param type One of \code{"xb"} for fitted values or \code{"residuals"}.
#' @param ... Currently unused.
#'
#' @return A numeric matrix of fitted values or residuals, depending on
#'   \code{type}.
#' @seealso [residuals.csdm_fit()], [summary.csdm_fit()]
#' @export
predict.csdm_fit <- function(object, newdata = NULL, type = c("xb", "residuals"), ...) {
  type <- match.arg(type)
  if (!is.null(newdata)) stop("predict(newdata=...) not implemented yet")
  if (type == "residuals") return(stats::residuals(object, type = "e"))
  if (!is.null(object$fitted_xb)) return(object$fitted_xb)
  stop("predict(type='xb') not implemented yet")
}
# Helper function to print CD tests
.csdm_print_cd_tests <- function(stats, digits = 4, classic_only = FALSE) {
  # Helper to print a single test
  print_one_test <- function(stat_val, p_val, label) {
    if (!is.null(stat_val) && !is.na(stat_val)) {
      cat(label, " = ", round(stat_val, digits), sep = "")
      if (!is.null(p_val) && !is.na(p_val)) {
        cat(", p = ", round(p_val, digits), sep = "")
      }
      cat("\n")
      return(TRUE)
    }
    return(FALSE)
  }

  printed_any <- FALSE

  # Always print classic CD
  if (print_one_test(stats$CD_stat, stats$CD_p, "CD")) {
    printed_any <- TRUE
  }

  # If not classic_only, print advanced tests
  if (!classic_only) {
    if (print_one_test(stats$CDw_stat, stats$CDw_p, "CDw")) printed_any <- TRUE
    if (print_one_test(stats$CDw_plus_stat, stats$CDw_plus_p, "CDw+")) printed_any <- TRUE

    # CD* includes n_pc in label
    if (!is.null(stats$CDstar_stat) && !is.na(stats$CDstar_stat)) {
      n_pc_label <- if (!is.null(stats$CDstar_n_pc) && !is.na(stats$CDstar_n_pc)) {
        paste0("CD* (n_pc=", stats$CDstar_n_pc, ")")
      } else {
        "CD* (n_pc=4)"
      }
      print_one_test(stats$CDstar_stat, stats$CDstar_p, n_pc_label)
      printed_any <- TRUE
    }
  }

  if (printed_any) {
    if (classic_only) {
      cat("(For additional CD diagnostics, use cd_test())\n")
    }
    cat("\n")
  }
}
