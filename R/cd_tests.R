cd <- function(data, id, time, var, na.action = na.omit, ...){

  # data |>
  #   select(id, year, !!sym(var) )
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

cd(data, id, time = year, var = residual)
