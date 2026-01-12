#' Penn World Tables panel (93 countries, 1960–2007)
#'
#' A panel of 93 countries (unit id) observed annually over 1960–2007 (time/year),
#' with the log-transformed variables used in 
#' xtdcce2-style examples.
#'
#' @format A data frame with 4464 rows and 6 variables:
#' \describe{
#'   \item{id}{Unit identifier (country id).}
#'   \item{year}{Time identifier (year, 1960–2007).}
#'   \item{log_rgdpo}{Log real GDP (output).}
#'   \item{log_hc}{Log human capital index.}
#'   \item{log_ck}{Log capital stock.}
#'   \item{log_ngd}{Log (net) government debt (or similar), used as a covariate/control.}
#' }
#' @source Penn World Table (PWT). This dataset is included as a small, convenient
#'   panel for examples and tests.
"PWT_60_07"
