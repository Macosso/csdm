#' Title
#'
#' @param data dataset
#' @param id panel data id indicator
#' @param time panel data time indicator
#' @param var variable to be tested
#' @param na.action na actions
#' @param ... other argunents
#'
#' @returns CD test
#' @export
#'
#' @importFrom dplyr select arrange
#' @importFrom tidyr pivot_wider
#' @importFrom rlang enquo eval_tidy
cd <- function(data, id, time, var, na.action = na.omit, ...){

  data <- na.action(data)

  panel_matrix <- data |>
    select({{id}}, {{time}}, {{var}}) |>
    pivot_wider(names_from = {{id}}, values_from = {{var}}) |>
    arrange({{time}}) |>
    select(-{{time}})

  rho <- cor(panel_matrix)

  N <- ncol(panel_matrix)
  TT <- nrow(panel_matrix)



  CD <- sqrt((2*TT)/(N*(N - 1)))*sum(rho[upper.tri(rho)])
  return(CD)
}

# cd(data, id, time = year, var = residual)
