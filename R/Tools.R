#' Predict cce_mg
#'
#' @param object an cee-mg model
#' @param data optional dataset
#' @param id id panel data id indicator
#' @param time time panel data time indicator
#'
#' @returns predictions
#' @export
#'
#' @importFrom dplyr select ungroup rename across group_by
predict.cce_mean_group <- function(object, data = NULL, id = NULL, time = NULL){

  if(is.null(data)){
    data <- object$data
  }

  if (!requireNamespace("plm", quietly = TRUE)) stop("Package 'plm' is required.")
  if (!requireNamespace("stats", quietly = TRUE)) stop("Package 'stats' is required.")

  if (!inherits(data, "pdata.frame")) {
    if (is.null(id) || is.null(time)) {
      stop("If data is not a pdata.frame, both 'id' and 'time' must be provided.")
    }
    data <- plm::pdata.frame(data, index = c(id, time))
  }

#  data <- na.action(data)

  data$id <- index(data)[[1]]
  data$time <- index(data)[[2]]

  model_cols <- all.vars(formula)


  coefs <- object$coefficients
  model_matrix <- data |>
    group_by(time) |>
    mutate(across(all_of(model_cols),
                  list(avg = function(x) mean(x, na.rm = TRUE)),
                  .names = "avg_X{.col}"),
           `(Intercept)` = 1) |>
    rename(avg_y = paste0("avg_X",as.character(formula)[2])) |>
    ungroup() |>
    select(names(coefs)) |>
    as.matrix()

  predictions <- as.vector(model_matrix %*% coefs)

  return(predictions)

}


