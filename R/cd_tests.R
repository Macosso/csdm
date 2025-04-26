# cd <- function(data, var, ...){
#
#   data |>
#     select(id, year, quote(var) )
#
#   # panel_matrix <- data %>%
#   #   select(id, year, !!sym(var) ) %>%
#   #   pivot_wider(names_from = id, values_from = !!sym(var)) %>%
#   #   arrange(year)
#
# }
#
# cd(data, var = log_rgdpo)
