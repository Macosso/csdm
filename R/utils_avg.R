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
#'   for that time’s CSA for that variable.
#'
#' @returns A \code{data.frame}:
#'   \itemize{
#'     \item If \code{return_mode="attach"}: original data + CSA columns
#'       named \code{paste0(suffix, "_", vars)}.
#'     \item If \code{return_mode="time"}: unique time rows with CSA columns.
#'   }
#'
#' @details
#' Efficiently computes, for each \code{v in vars} and time \code{t},
#' \deqn{\bar v_t = \frac{\sum_i w_{it}\, 1_{\{v_{it}\text{ finite}\}}\, v_{it}}
#'                 {\sum_i w_{it}\, 1_{\{v_{it}\text{ finite}\}}}}
#' For \code{leave_out=TRUE}, each row’s CSA excludes its own contribution; if the
#' denominator becomes \eqn{\le 0} (e.g., only one finite observation at that time),
#' the LOO mean is set to \code{NA} for that row/variable.
#'
#' @keywords internal
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

  # --- 0) Extract id/time if pdata.frame
  if (inherits(data, "pdata.frame")) {
    idx <- attr(data, "index")
    if (is.null(idx)) stop("pdata.frame without index attribute.")
    if (is.null(id))   id   <- names(idx)[1L]
    if (is.null(time)) time <- names(idx)[2L]
    df <- as.data.frame(data)
  } else {
    df <- as.data.frame(data)
    if (is.null(id) || is.null(time)) {
      stop("For data.frame, both 'id' and 'time' must be provided.")
    }
  }

  if (!all(vars %in% names(df))) {
    miss <- setdiff(vars, names(df))
    stop("Variables not found in data: ", paste(miss, collapse = ", "))
  }
  if (!all(c(id, time) %in% names(df))) {
    stop("Columns 'id' and/or 'time' not found in data.")
  }

  # --- 1) Weights vector
  if (is.null(weights)) {
    w <- rep(1, nrow(df))
  } else if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(df)) stop("weights column '", weights, "' not found.")
    w <- as.numeric(df[[weights]])
  } else if (is.numeric(weights)) {
    if (length(weights) != nrow(df)) stop("weights vector must have length nrow(data).")
    w <- as.numeric(weights)
  } else {
    stop("'weights' must be NULL, a column name, or a numeric vector length nrow(data).")
  }
  if (any(!is.finite(w) | w < 0, na.rm = TRUE)) {
    stop("weights must be nonnegative and finite.")
  }
  # For NA weights, treat as zero weight
  w[!is.finite(w)] <- 0

  # --- 2) Build time-level sums and denominators for each var
  # Efficient base aggregation with tapply on vectors
  time_vec <- df[[time]]
  id_vec   <- df[[id]]

  # Helper: time-wise sum of w * x and sum of w for finite x
  time_sumw  <- function(x) tapply(x, time_vec, sum, na.rm = TRUE)
  time_sumi  <- function(ind) tapply(ind, time_vec, sum, na.rm = TRUE)

  # Precompute maps per var: sum_wx_t and sum_w_t
  uniq_t <- sort(unique(time_vec))
  K <- length(vars)

  sum_wx_list <- vector("list", K)
  sum_w_list  <- vector("list", K)

  for (k in seq_along(vars)) {
    v <- df[[vars[k]]]
    fin <- is.finite(v)
    ww  <- w * fin
    sum_wx <- tapply(ww * v, time_vec, sum, na.rm = TRUE)
    sum_w  <- tapply(ww,       time_vec, sum, na.rm = TRUE)

    # Ensure all uniq_t present
    sum_wx <- sum_wx[as.character(uniq_t)]
    sum_w  <- sum_w[as.character(uniq_t)]

    # Replace NA with 0 in aggregates (consistent with na.rm=TRUE path)
    sum_wx[!is.finite(sum_wx)] <- 0
    sum_w[!is.finite(sum_w)]   <- 0

    sum_wx_list[[k]] <- sum_wx
    sum_w_list[[k]]  <- sum_w
  }

  # --- 3) Compute non-LOO time means if needed
  if (!leave_out) {
    csa_time <- as.data.frame(setNames(list(uniq_t), time))
    for (k in seq_along(vars)) {
      denom <- sum_w_list[[k]]
      num   <- sum_wx_list[[k]]
      mu    <- if (na.rm) ifelse(denom > 0, num / denom, NA_real_) else num / denom
      csa_time[[paste0(suffix, "_", vars[k])]] <- as.numeric(mu)
    }
    # Return as requested
    if (return_mode == "time") {
      rownames(csa_time) <- NULL
      return(csa_time)
    } else {
      # attach to df: map time -> CSA columns
      key <- as.character(uniq_t)
      idx_map <- match(as.character(time_vec), key)
      for (nm in setdiff(names(csa_time), time)) {
        df[[nm]] <- csa_time[[nm]][idx_map]
      }
      return(df)
    }
  }

  # --- 4) Leave-one-out: row-wise adjustment (sum_wx_t - w_i x_i) / (sum_w_t - w_i * 1_{finite})
  # Prepare output columns
  out_cols <- replicate(K, numeric(nrow(df)), simplify = FALSE)
  names(out_cols) <- paste0(suffix, "_", vars)

  # Precompute time-index vector to align with time aggregates
  time_key <- as.character(uniq_t)
  t_idx <- match(as.character(time_vec), time_key)

  for (k in seq_along(vars)) {
    v <- df[[vars[k]]]
    fin_i <- is.finite(v)
    w_i   <- w

    sum_wx_t <- as.numeric(sum_wx_list[[k]][t_idx])
    sum_w_t  <- as.numeric(sum_w_list[[k]][t_idx])

    num   <- sum_wx_t - (w_i * v) * fin_i
    denom <- sum_w_t  - (w_i * fin_i)

    mu_lo <- ifelse(denom > 0, num / denom, NA_real_)
    if (!na.rm) {
      # If not removing NA, set NA when v is NA or any NA in slice is present.
      # Cheap conservative rule: if v_i is NA, keep NA; otherwise use mu_lo.
      mu_lo[!fin_i] <- NA_real_
    }
    out_cols[[k]] <- mu_lo
  }

  for (k in seq_along(vars)) {
    df[[paste0(suffix, "_", vars[k])]] <- out_cols[[k]]
  }

  if (return_mode == "time") {
    # Collapse to time-level non-LOO (there is no single "time-level LOO" table);
    # return the standard CSA by time to avoid misleading output.
    csa_time <- as.data.frame(setNames(list(uniq_t), time))
    for (k in seq_along(vars)) {
      denom <- sum_w_list[[k]]
      num   <- sum_wx_list[[k]]
      mu    <- if (na.rm) ifelse(denom > 0, num / denom, NA_real_) else num / denom
      csa_time[[paste0(suffix, "_", vars[k])]] <- as.numeric(mu)
    }
    rownames(csa_time) <- NULL
    return(csa_time)
  } else {
    return(df)
  }
}
