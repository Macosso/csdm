#' @export
summary.mg_estimator <- function(object, digits = 4, ...) {
  # Extract from object:
  mg_coef <- object$mg_coef
  mg_var  <- object$mg_var
  N       <- object$N

  # Calculate standard errors, t-statistics, p-values (using df = N - 1)
  mg_se <- sqrt(mg_var)
  t_val <- mg_coef / mg_se
  df    <- N - 1
  p_val <- 2 * pt(abs(t_val), df = df, lower.tail = FALSE)

  # Significance stars
  get_stars <- function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else if (p < 0.10) return(".")
    else return("")
  }
  stars <- vapply(p_val, get_stars, character(1))

  # Combine into a data.frame for neat display
  results_table <- data.frame(
    Estimate = mg_coef,
    `Std.Err` = mg_se,
    t.value   = t_val,
    `Pr(>|t|)`= p_val,
    sign.     = stars,
    row.names = names(mg_coef)
  )

  # Optionally, you might compute (or store) additional summary stats here:
  # Example placeholders:
  #   R-squared, residual variance, MG-based CD statistic, etc.
  # E.g., R-squared is not directly computed by mg_estimator but might be
  # average of the group-specific R-squareds or something else:
  # mg_rsq <- NA_real_  # placeholder
  # mg_cd  <- NA_real_  # placeholder for CD test statistic

  # Build a return list
  out <- list(
    call       = object$formula,
    coefficients = results_table,
    df        = df,
    N         = N
    # mg_rsq   = mg_rsq,
    # mg_cd    = mg_cd
  )
  class(out) <- "summary.mg_estimator"
  return(out)
}

# A corresponding print method for summary.mg_estimator:
#' @export
print.summary.mg_estimator <- function(x, digits = 4, ...) {
  cat("Mean Group Estimation (Pesaran & Smith, 1995)\n")
  cat("Formula: ", deparse(x$call), "\n")
  cat("Number of groups (N): ", x$N, "\n")
  cat("Degrees of freedom:   ", x$df, "\n\n")

  cat("MG Coefficients:\n")
  # Format numeric columns with the given precision
  fmt_df <- x$coefficients
  fmt_df[, 1:4] <- lapply(fmt_df[, 1:4], function(col) formatC(col, digits = digits, format = "f"))
  print(fmt_df, row.names = TRUE, right = TRUE)

  # If you computed or stored other summary statistics:
  # cat("\nAverage R-squared:", formatC(x$mg_rsq, digits = digits, format = "f"), "\n")
  # cat("CD Statistic:", formatC(x$mg_cd, digits = digits, format = "f"), "\n")

  cat("\nSignif. codes:  0 ***  0.001 **  0.01 *  0.05 .  0.1\n\n")
  invisible(x)
}


#' @export
predict.mg_estimator <- function(object, newdata = NULL, id = NULL, time = NULL, ...) {
  # If newdata is not provided, use the stored original data.
  if (is.null(newdata)) {
    if (is.null(object$model)) {
      stop("No newdata provided and original data was not stored in the mg_estimator object.")
    }
    newdata <- object$model
  } else {
    # If newdata is a regular data.frame, convert it using id and time names.
    if (!("pdata.frame" %in% class(newdata))) {
      if (is.null(id) || is.null(time)) {
        stop("For data.frame input in newdata, please supply both 'id' and 'time' variable names.")
      }
      if (!requireNamespace("plm", quietly = TRUE)) {
        stop("Package 'plm' is required. Please install it.")
      }
      newdata <- plm::pdata.frame(newdata, index = c(id, time))
    }
  }
  # Compute the model matrix based on the formula.
  X <- model.matrix.lm(object$formula, data = newdata, na.action = "na.pass")
  # Predict using the mean group (MG) coefficients.
  pred <- as.vector(X %*% object$mg_coef)
  return(pred)
}
