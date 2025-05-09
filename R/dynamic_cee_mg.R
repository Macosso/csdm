#' Dynamic Common Correlated Effects Mean Group (CCE-MG) Estimator
#'
#' @description
#' Implements the levels ARDL(1) + Common Correlated Effects Mean Group estimator with user-specified
#' number of lags on cross-section averages. Supports a formula interface and works with data.frame,
#' tibble, or `plm::pdata.frame` objects.
#'
#' @param formula A two-sided formula of the form `y ~ x1 + x2 + ...` specifying the dependent and independent variables.
#' @param data A data.frame, tibble, or `plm::pdata.frame` containing the panel data.
#' @param id Character. Name of the cross-sectional identifier column (required if `data` is not a pdata.frame).
#' @param time Character. Name of the time identifier column (required if `data` is not a pdata.frame).
#' @param p Integer. Number of lags for cross-sectional averages (default = 3).
#'
#' @return A tibble with the short-run mean-group coefficients, standard errors, and z-statistics
#'   for the lagged dependent variable and each regressor.
#'
#' @importFrom dplyr group_by mutate ungroup arrange filter across all_of sym if_any summarise
#' @importFrom stats lm sd model.frame as.formula terms
#' @export
#'
#' @examples
#' # Using formula with data.frame:
#' # dynamic_cce_mg(log_rgdpo ~ log_hc + log_ck + log_ngd, data = my_df, id = "id", time = "year", p = 3)
#' # Using pdata.frame (id/time inferred):
#' # pf <- plm::pdata.frame(my_df, index = c("id","year"))
#' # dynamic_cce_mg(log_rgdpo ~ log_hc + log_ck + log_ngd, data = pf, p = 3)
#'
dynamic_cce_mg <- function(formula, data, id = NULL, time = NULL, p = 3, na.action = na.omit) {

  data <- na.action(data)
  # Handle pdata.frame vs data.frame
  if (inherits(data, "pdata.frame")) {
    idx <- attr(data, "index")
    id <- idx[1]; time <- idx[2]
    df <- as.data.frame(data)
  } else {
    if (is.null(id) || is.null(time)) {
      stop("When 'data' is not a pdata.frame, please provide 'id' and 'time' arguments.")
    }
    df <- as.data.frame(data)
  }

  # Validate and parse formula
  # model.frame ensures variables exist and handles transformations
  mf <- model.frame(formula, data = df)
  vars <- all.vars(formula)
  y <- vars[1]
  x_vars <- vars[-1]

  # Compute cross-sectional averages of y and x_vars
  avg_names <- paste0("bar_", vars)
  df <- df %>%
    group_by(!!sym(time)) %>%
    mutate(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "bar_{.col}")) %>%
    ungroup()

  # Create lag of dependent variable and lags of averages
  df <- df %>% arrange(!!sym(id), !!sym(time)) %>%
    group_by(!!sym(id)) %>%
    mutate(
      lag_y = lag(!!sym(y), 1)
    )
  for (l in seq_len(p)) {
    df <- df %>%
      mutate(across(all_of(avg_names), ~ lag(.x, l), .names = "{.col}_lag{l}"))
  }
  df <- df %>% ungroup()

  # Define regressors
  lag_vars <- as.vector(outer(paste0("bar_", vars, "_lag"), seq_len(p), FUN = "paste0"))
  regressors <- c("lag_y", x_vars, paste0("bar_", vars), lag_vars)

  # Filter complete cases
  df <- df %>% filter(!if_any(all_of(c(y, regressors)), is.na))

  # Unit-specific ARDL + CCE and collect short-run coefs
  coef_df <- df %>%
    group_by(!!sym(id)) %>%
    group_modify(~ {
      fmla <- as.formula(paste(y, "~", paste(regressors, collapse = " + ")))
      fit <- lm(fmla, data = .x)
      # coefs <- coef(fit)[c("lag_y", x_vars)]
      coefs <- coef(fit)
      tibble::as_tibble_row(coefs)
    }) %>% ungroup()

  # Mean-group summary
  mg <- coef_df %>%
    summarise(across(all_of(c("lag_y", regressors)),
                     list("-coef" = mean, "-se" = ~ sd(.x)/sqrt(n())))) %>%
    tidyr::pivot_longer(everything(), names_to = c("term","stat"),
                        names_sep = "_-", values_to = "value") %>%
    tidyr::pivot_wider(names_from = stat, values_from = value) %>%
    mutate(z = coef / se)

  coefs <- setNames(mg$coef, mg$term)
  se <- setNames(mg$se, mg$term)

  return(
    list(
      coef = coefs,
      se = se,
      call = as.formula(paste(y, "~", paste(regressors, collapse = " + ")))

    )
  )
}
