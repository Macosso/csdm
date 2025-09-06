#' Dynamic Common Correlated Effects Mean Group (CCE-MG) Estimator
#'
#' @description
#' Implements the levels ARDL(1) + Common Correlated Effects Mean Group estimator
#' with user-specified number of lags on cross-section averages. Supports a formula
#' interface and works with data.frame, tibble, or `plm::pdata.frame` objects.
#'
#' @param formula A two-sided formula of the form `y ~ x1 + x2 + ...` specifying the dependent and independent variables.
#' @param data A data.frame, tibble, or `plm::pdata.frame` containing the panel data.
#' @param id Character. Name of the cross-sectional identifier column (required if `data` is not a pdata.frame).
#' @param time Character. Name of the time identifier column (required if `data` is not a pdata.frame).
#' @param p Integer. Number of lags for cross-sectional averages (default = 3).
#' @param na.action Function to handle missing values (default = na.omit).
#'
#' @return An S3 object of class 'dynamic_cce_mg' containing:
#'   \item{coef}{Mean-group short-run coefficients (named vector).}
#'   \item{se}{Standard errors for the mean-group coefficients (named vector).}
#'   \item{r.squared}{Per-unit R-squared values (named vector).}
#'   \item{residuals}{Data frame of unit, time, and residuals for each observation.}
#'   \item{call}{The final levels ARDL formula used.}
#'   \item{p}{Lag order on cross-sectional averages.}
#'   \item{id, time}{Index names.}
#'
#' @importFrom dplyr group_by ungroup arrange mutate filter across all_of sym summarise nest unnest unnest_wider
#' @importFrom purrr map map_dbl map2
#' @importFrom tidyr nest data unnest_wider unnest
#' @importFrom stats lm sd model.frame as.formula
#' @export
#'
#' @examples
#' pf <- plm::pdata.frame(my_df, index = c("id","year"))
#' fit <- dynamic_cce_mg(log_rgdpo ~ log_hc + log_ck + log_ngd, data = pf, p = 3)
#' summary(fit)

dynamic_cce_mg <- function(formula, data, id = NULL, time = NULL,
                           p = 3, na.action = na.omit) {
  # process panel index
  if (inherits(data, "pdata.frame")) {
    idx <- names(index(data)); id <- idx[1]; time <- idx[2]
    df <- as.data.frame(data)
  } else {
    if (is.null(id) || is.null(time)) stop(
      "When 'data' is not a pdata.frame, please provide 'id' and 'time'."
    )
    df <- as.data.frame(data)
  }

  # parse formula
  mf <- model.frame(formula, data = df)
  vars <- all.vars(formula); y <- vars[1]; x_vars <- vars[-1]

  # cross-section averages
  df <- df |>
    select(!!sym(id), !!sym(time), all_of(c(y, x_vars))) |>
    group_by(!!sym(time)) |>
    mutate(across(all_of(vars), ~ mean(.x, na.rm=TRUE), .names="cs_{.col}")) |>
    ungroup() |>
    na.omit()

  # lags
  df <- df |> arrange(!!sym(id), !!sym(time)) |> group_by(!!sym(id)) |>
    mutate("lag_{y}" := lag(!!sym(y),1))
  avg_names <- paste0("cs_", vars)
  for (l in seq_len(p)) {
    df <- df |> mutate(across(all_of(avg_names), ~ lag(.x, l),
                               .names="{.col}_lag{l}") )
  }
  df <- df |> ungroup()

  # prepare regressors and filter
  lag_labs <- as.vector(outer(avg_names, seq_len(p),
                              FUN=function(v,l) paste0(v, "_lag", l)))
  lag_labs <- c(paste0("lag_", y), lag_labs)

  regressors <- c(lag_labs, x_vars, avg_names)
  regressors <- unique(regressors)
  if (!all(regressors %in% names(df))) {
    stop("Some regressors are not present in the data.")
  }
  if (length(regressors) == 0) {
    stop("No regressors found in the data.")
  }

  # filter out rows with NAs in y or regressors
  df2 <- na.omit(df)

  nG <- length(unique(df2[[id]]))
  nT <- length(unique(df2[[time]]))

  # obtain the final model fram
  model_matrix <- model.matrix(
    as.formula(paste(y, "~", paste(regressors, collapse=" + "))),
                                  data = df2, na.action = na.omit)

  # nest by unit
  nested <- df2 |>
    # group_by(!!sym(id)) |>
    nest(data = everything(), .by = !!sym(id)) |>
    select(-!!sym(id))

  # fit models and extract
  nested <- nested |>
    mutate(
      model = map(data, ~ lm(
        as.formula(paste(y, "~", paste(regressors, collapse=" + "))),
        data = .x
      )),
      coef = map(model, ~ coef(.x)),
      r.squared = map_dbl(model, ~ summary(.x)$r.squared),
      residuals = map2(data, model, ~ tibble(
        !!id := .x[[id]],
        !!time := .x[[time]],
        residual = resid(.y)
      ))
    )

  # build output components
  coef_df <- nested |> unnest_wider(coef) |> select(-data, -model, -r.squared, -residuals)
  res_df <- nested |> select(residuals) |> unnest(residuals)
  r2 <- setNames(nested$r.squared, nested[[id]])

  # compute mean-group
  mg_stats <- coef_df |>
    # select(-!!sym(id)) |>
    summarise(across(everything(), list("--coef" = mean, "--se" = ~ sd(.x)/sqrt(n()))))
  mg_tidy <- mg_stats |>
    tidyr::pivot_longer(everything(), names_to = c("term","stat"),
                        names_sep="_--", values_to = "value") |>
    tidyr::pivot_wider(names_from = stat, values_from = value) |>
    mutate(z = coef / se)

  out <- list(
    coef = setNames(mg_tidy$coef, mg_tidy$term)[c("(Intercept)", paste0("lag_", y), x_vars)],
    se = setNames(mg_tidy$se, mg_tidy$term)[c("(Intercept)", paste0("lag_", y), x_vars)],
    coef_df = coef_df,
    r.squared = r2,
    residuals = plm::pdata.frame(res_df, index = c(id, time)),
    call = formula,
    p = p,
    id = id,
    time = time,
    nG = nG,
    nT = nT,
    data = df2,
    df = nT - length(vars) - 1,
    dfcs = nT - length(regressors) - 1,
    regressors = regressors,
    cs_vars = vars
  )
  class(out) <- "dynamic_cce_mg"
  return(out)
}



#' Print method for dynamic_cce_mg
#' @export
print.dynamic_cce_mg <- function(x, ...) {
  cat("Dynamic CCE-MG Mean-Group Estimation\n")
  cat("Call:", deparse(x$call), "\n")
  cat("Lag order (p):", x$p, "\n\n")
  cat("Mean-group coefficients:\n")
  print(round(x$coef, 4))
  cat("\nStandard errors:\n")
  print(round(x$se, 4))
  invisible(x)
}



#' Summary method for dynamic_cce_mg
#' @export
summary.dynamic_cce_mg <- function(object, ...) {

  mg <- data.frame(
    row.names = names(object$coef),
    estimate = object$coef,
    std.error = object$se,
    z = object$coef / object$se,
    stringsAsFactors = FALSE
  ) |>
    mutate(
      p.value = 2 * (1 - pnorm(abs(z))),
      significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      )
    ) |>
    select(estimate, std.error, z, p.value, significance)

  r.squared <- mean(object$r.squared)

  # CD test
  CD <- cd(
    object$residuals,
    var = "residual"
  )

  out <- list(
    call = object$call,
    coefficients = mg,
    r.squared = r.squared,
    residuals = object$residuals,
    CD = CD,
    id  = object$id,
    time = object$time,
    nG = object$nG,
    nT = object$nT,
    p = object$p,
    df = object$df,
    dfcs = object$dfcs,
    data = object$data,
    regressors = object$regressors,
    cs_vars = object$cs_vars
  )
  class(out) <- "summary.dynamic_cce_mg"
  return(out)
}


print.summary.dynamic_cce_mg  <- function(x, digits = 7, ...){
  cat("Dynamic Common Correlated Error - Mean Group (CCE-MG) Estimation\n")
  cat("Call:", deparse(x$call), "\n")

  cat("\nNumber of observations: N = ", x$nG * x$nT, "\tn = ", x$nG, "\tT = ", x$nT, "\n")

  cat("\nMean-group coefficients:\n")
  x$coefficients$p.value <- format(round(x$coefficients$p.value,
                                   digits = min(digits, 3)),
                                   nsmall =  min(digits, 3))
  x$coefficients$z <- format(round(x$coefficients$z,
                             digits = min(digits, 3)),
                             nsmall =  min(digits, 3))

  x$coefficients$estimate <- format(round(x$coefficients$estimate,
                             digits = digits),  nsmall =  digits)

  x$coefficients$std.error <- format(round(x$coefficients$std.error,
                                    digits = digits),  nsmall =  digits)
  print(x$coefficients)
  cat("---", "\nSignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")



  cat("\nR-squared (mg):", format(x$r.squared, digits = 2, nsmall = 2), "\n")

  cat("Cross-sectional Dependence test (CD):\tstatistic = ",
      format(x$CD$statistic, digits = 4, nsmall = 4),
      "\tp-value = ", format(x$CD$p.value, digits = 3, nsmall = 3), "\n")

  cat("Degrees of freedom per group: \twithout cs averages = ", x$df,
      "\twith cs averages = ", x$dfcs, "\n\n")

  cat("cs lags = ", x$p,
      "\t variables in mean group reg = ",
      (length(x$regressors[!startsWith(x$regressors, "cs_")]) + 1) * x$nG,
      "\n")

  cat("Mean Group Variables: ", c("(Intercept)", x$regressors[!startsWith(x$regressors, "cs_")]))
  cat("\nCross Sectional Averaged Variables: ", x$cs_vars)

  invisible(x)
}


vcov.dynamic_cce_mg  <- function(object, ...) {
  # Extract the residuals
  res <- object$residuals$residual

  # Compute the covariance matrix of the residuals
  cov_matrix <- cov(res)

  # Return the covariance matrix
  return(cov_matrix)
}


coef.dynamic_cce_mg <- function(object, ...) {
  return(object$coef)
}
