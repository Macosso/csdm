# utils_avg.R

#' Cross-sectional averages by time (with optional leave-one-out)
#'
#' @description
#' Computes cross-sectional averages (CSAs) of specified variables for each time
#' period, optionally in a leave-one-out (LOO) fashion per observation. Supports
#' unbalanced panels and observation weights.
#'
#' @param data A \code{data.frame} or \code{plm::pdata.frame}.
#' @param id,time Character scalar names of unit and time columns when \code{data}
#'   is a plain \code{data.frame}. If \code{data} is a \code{pdata.frame}, these are
#'   inferred from its index and can be omitted.
#' @param vars Character vector of column names to average cross-sectionally.
#' @param leave_out Logical; if \code{TRUE}, computes LOO means for each row:
#'   \eqn{\bar{x}_{-i,t} = (\sum_{j \neq i} w_{jt} x_{jt}) / (\sum_{j \neq i} w_{jt})}.
#'   If \code{FALSE}, computes standard time means:
#'   \eqn{\bar{x}_{t} = (\sum_j w_{jt} x_{jt}) / (\sum_j w_{jt})}.
#' @param weights Optional. Either:
#'   \itemize{
#'     \item a numeric vector of length \code{nrow(data)}, or
#'     \item the name of a column in \code{data} with nonnegative weights.
#'   }
#'   If \code{NULL}, uses equal weights (1 for observed values).
#' @param suffix Character suffix to append to CSA columns (default \code{"csa"}).
#' @param return_mode One of \code{"attach"} or \code{"time"}.
#'   \itemize{
#'     \item \code{"attach"} returns the original \code{data} with CSA columns added.
#'     \item \code{"time"} returns a unique-time table \code{[time, csa_*]}.
#'   }
#' @param na.rm Logical; if \code{TRUE}, excludes \code{NA}s from sums and
#'   denominators. If \code{FALSE}, any \code{NA} in a time slice yields \code{NA}
#'   for that time's CSA for that variable.
#'
#' @returns A \code{data.frame}:
#'   \itemize{
#'     \item If \code{return_mode="attach"}: original data + CSA columns
#'       named \code{paste0(suffix, "_", vars)}.
#'     \item If \code{return_mode="time"}: unique time rows with CSA columns.
#'   }
#'
#' @details
#' This is a standalone data utility. It does not configure the averages used by
#' `csdm()`; use [csdm_csa()] for that purpose. Model fitting constructs averages
#' from the evaluated model terms and its documented source sample, which can
#' differ from averages of raw data columns produced here.
#'
#' Efficiently computes, for each \code{v in vars} and time \code{t},
#' \deqn{\bar v_t = \frac{\sum_i w_{it}\, 1_{\{v_{it}\text{ finite}\}}\, v_{it}}
#'                 {\sum_i w_{it}\, 1_{\{v_{it}\text{ finite}\}}}}
#' For \code{leave_out=TRUE}, each row's CSA excludes its own contribution; if the
#' denominator becomes \eqn{\le 0} (e.g., only one finite observation at that time),
#' the LOO mean is set to \code{NA} for that row/variable.
#'
#' @export
cross_sectional_avg <- function(data,
                                id = NULL,
                                time = NULL,
                                vars,
                                leave_out = FALSE,
                                weights = NULL,
                                suffix = "csa",
                                return_mode = c("attach", "time"),
                                na.rm = TRUE) {
  return_mode <- match.arg(return_mode)

  .csdm_flag(leave_out, "leave_out")
  .csdm_flag(na.rm, "na.rm")
  if (leave_out && return_mode == "time") stop("Leave-one-out averages require return_mode='attach'.")
  if (inherits(data, "pdata.frame")) {
    idx <- attr(data, "index")
    if (is.null(id)) id <- names(idx)[1L]
    if (is.null(time)) time <- names(idx)[2L]
  }
  df <- as.data.frame(data)
  if (!length(id) || !length(time) || !all(c(id, time) %in% names(df))) stop("Specify valid id/time columns.")
  if (!is.character(vars) || !length(vars) || anyNA(vars) || anyDuplicated(vars) ||
      !all(vars %in% names(df))) stop("'vars' must name unique numeric columns.")
  if (anyNA(df[[id]]) || anyNA(df[[time]]) || anyDuplicated(df[c(id, time)])) stop("Invalid or duplicate panel keys.")
  if (!all(vapply(df[vars], is.numeric, logical(1)))) stop("CSA variables must be numeric.")
  if (!is.character(suffix) || length(suffix) != 1L || is.na(suffix) || !nzchar(suffix)) stop("Invalid suffix.")
  out_names <- paste0(suffix, "_", vars)
  if (any(out_names %in% names(df))) stop("CSA column names already exist in data.")
  if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(df)) stop("Unknown weights column.")
    weights <- df[[weights]]
  }
  w <- if (is.null(weights)) rep(1, nrow(df)) else weights
  if (!is.numeric(w) || length(w) != nrow(df) || any(!is.finite(w) | w < 0)) {
    stop("Weights must be nonnegative, finite, and have length nrow(data).")
  }
  times <- sort(unique(df[[time]]))
  index <- match(df[[time]], times)
  time_result <- data.frame(times, check.names = FALSE)
  names(time_result) <- time
  for (j in seq_along(vars)) {
    x <- df[[vars[j]]]
    finite <- is.finite(x)
    contribution <- ifelse(finite, x, 0) * w
    denominator <- w * finite
    sums <- as.numeric(rowsum(contribution, index, reorder = TRUE))
    denoms <- as.numeric(rowsum(denominator, index, reorder = TRUE))
    missing <- as.numeric(rowsum(as.integer(!finite), index, reorder = TRUE))
    if (leave_out) {
      numerator <- sums[index] - contribution
      denom <- denoms[index] - denominator
      missing <- missing[index] - !finite
    } else {
      numerator <- sums
      denom <- denoms
    }
    value <- rep(NA_real_, length(denom))
    ok <- denom > 0 & (na.rm | missing == 0)
    value[ok] <- numerator[ok] / denom[ok]
    if (leave_out) df[[out_names[j]]] <- value else {
      time_result[[out_names[j]]] <- value
      df[[out_names[j]]] <- value[index]
    }
  }
  if (return_mode == "time") time_result else df
}
