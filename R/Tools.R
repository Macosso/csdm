
predict.cce_mean_group <- function(object, data = NULL, unit = NULL, time = NULL){

  if(is.null(data)){
    data <- object$data
  }

  if (!requireNamespace("plm", quietly = TRUE)) stop("Package 'plm' is required.")
  if (!requireNamespace("stats", quietly = TRUE)) stop("Package 'stats' is required.")

  if (!inherits(data, "pdata.frame")) {
    if (is.null(unit) || is.null(time)) {
      stop("If data is not a pdata.frame, both 'unit' and 'time' must be provided.")
    }
    data <- plm::pdata.frame(data, index = c(unit, time))
  }

#  data <- na.action(data)

  data$unit <- index(data)[[1]]
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
#
#
# data$pred <- predict(object = object)
#
# data$residual <- data$log_rgdpo - data$pred
#
#
# tss <- sum((data$log_rgdpo - mean(data$log_rgdpo, na.rm = TRUE))^2)
# rss <- sum((data$residual)^2, na.rm = TRUE)
#
#
# 1 - rss/tss
#
#
#
# sp2 = sum(i=1,N) e(i)'e(i) / [N ( T - k - 2) - k]
#
#
#
#

