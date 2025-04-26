mg_estimator <- function(formula, data, id = NULL, time = NULL, na.action = na.omit) {
  # If the data is not already a pdata.frame, require that id and time variable names are provided.
  if (!("pdata.frame" %in% class(data))) {
    if (is.null(id) || is.null(time)) {
      stop("For data.frame input, please supply both 'id' and 'time' variable names.")
    }
    if (!requireNamespace("plm", quietly = TRUE)) {
      stop("Package 'plm' is required. Please install it.")
    }
    data <- plm::pdata.frame(data, index = c(id, time))
  }

  # Extract the cross-sectional (panel) index
  panel_index <- plm::index(data)[, 1]
  groups <- unique(panel_index)
  N <- length(groups)

  # Initialize a list to store the coefficients for each group
  coef_list <- vector("list", length = N)
  names(coef_list) <- as.character(groups)

  rsquared_vec <- c()
  # Loop over each group (cross-sectional unit)
  for (g in groups) {
    subset_data <- data[panel_index == g, ]
    fit <- lm(formula, data = subset_data)
    coef_list[[as.character(g)]] <- coef(fit)
    rsquared_vec[g] <- summary(fit)$r.squared
  }

  # Combine individual coefficients into a matrix (each row corresponds to one unit)
  coef_mat <- do.call(rbind, coef_list)

  # Compute the MG estimator as the mean across cross-sectional units
  mg_coef <- colMeans(coef_mat, na.rm = TRUE)

  # Compute an estimate of the variance of the MG estimator (variance across units divided by N)
  mg_var <- apply(coef_mat, 2, var, na.rm = TRUE) / N

  # Return the results as a list with class "mg_estimator"
  result <- list(
    mg_coef = mg_coef,
    mg_var = mg_var,
    individual_coefs = coef_mat,
    N = N,
    formula = formula,
    model = data,
    r.squared = mean(rsquared_vec)
  )
  class(result) <- "mg_estimator"
  return(result)
}

# Print method for objects returned by mg_estimator
print.mg_estimator <- function(x, ...) {
  cat("Mean Group Estimator (MG) results (Pesaran & Smith, 1995):\n")
  cat("Panel Model Formula: ", deparse(x$formula), "\n")
  cat("Number of obs:", x$N, "\n\n")
  cat("Number of groups:", x$N, "\n\n")
  cat("MG Coefficient Estimates:\n")
  print(x$mg_coef)
  cat("\nEstimated Variance of the MG Coefficients (across groups / N):\n")
  print(x$mg_var)
  invisible(x)
}



