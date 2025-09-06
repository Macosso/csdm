#' Title
#'
#' @param data dataset
#' @param id panel data id indicator
#' @param time panel data time indicator
#' @param var variable to be tested
#' @param na.action na actions
#' @param ... other argunents
#'
#' @returns CD test statistic and p-value
#' @export
#'
#' @importFrom dplyr select arrange
#' @importFrom tidyr pivot_wider
#' @importFrom rlang enquo eval_tidy
cd <- function(data, var, id = NULL, time = NULL, na.action = na.omit, ...){

  data <- na.action(data)
  if (inherits(data, "pdata.frame")) {
    idx <- attr(data, "index")
    id <- idx[1]
    time <- idx[[2]]
    data <- data |>
      mutate(id = id,
             time = time)

    id <- "id"
    time <- "time"

  } else {
    if (is.null(id) || is.null(time)) {
      stop("If data is not a pdata.frame, both 'id' and 'time' must be provided.")
    }
  }
  panel_matrix <- data |>
    select({{id}}, {{time}}, {{var}}) |>
    pivot_wider(names_from = {{id}}, values_from = {{var}}) |>
    arrange({{time}}) |>
    select(-{{time}})

  rho <- cor(panel_matrix)

  N <- ncol(panel_matrix)
  TT <- nrow(panel_matrix)



  cd_stat <- sqrt((2*TT)/(N*(N - 1)))*sum(rho[upper.tri(rho)])
  cd_pval <- 2 * (1 - pnorm(abs(cd_stat)))
  CD <- list(
    statistic = cd_stat,
    p.value = cd_pval,
    method = "Pesaran's CD test",
    alternative = "cross-sectional independence"
  )

  return(CD)
}

# cd(data, id, time = year, var = residual)
