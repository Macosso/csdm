#' Penn World Tables panel (93 countries, 1960-2007)
#'
#' A panel of 93 countries (unit id) observed annually over 1960-2007 (time/year),
#' with the log-transformed variables used in 
#' xtdcce2-style examples.
#'
#' @format A data frame with 4464 rows and 6 variables:
#' \describe{
#'   \item{id}{Unit identifier (country id).}
#'   \item{year}{Time identifier (year, 1960-2007).}
#'   \item{log_rgdpo}{Log real output.}
#'   \item{log_hc}{Log human capital index.}
#'   \item{log_ck}{Log physical capital.}
#'   \item{log_ngd}{Log population growth plus a 5 percent break-even investment rate.}
#' }
#' @source Penn World Table 8 example data distributed with Stata's
#'   \code{xtdcce2}. The variable descriptions follow the accompanying
#'   \code{xtdcce2} documentation.
"PWT_60_07"
