# csdm_internal_panel.R

.csdm_prepare_panel_df <- function(data, id, time, time_step = 1) {
  if (inherits(data, "pdata.frame")) {
    df <- as.data.frame(data)
    idx <- attr(data, "index")
    if (!is.null(idx)) {
      id <- names(idx)[1L]
      time <- names(idx)[2L]
    }
  } else {
    df <- as.data.frame(data)
  }

  if (!is.character(id) || length(id) != 1L || is.na(id) ||
      !is.character(time) || length(time) != 1L || is.na(time) || id == time) {
    stop("'id' and 'time' must name distinct columns.", call. = FALSE)
  }
  if (!all(c(id, time) %in% names(df))) stop("'data' must contain 'id' and 'time' columns.")
  if (!nrow(df) || anyDuplicated(names(df))) stop("Supply nonempty data with unique column names.")
  if (any(startsWith(names(df), ".csdm_"))) stop("Column names beginning '.csdm_' are reserved.")

  if (inherits(data, "pdata.frame") && is.factor(df[[time]])) {
    df[[time]] <- suppressWarnings(as.numeric(as.character(df[[time]])))
  }

  # time is required and must be numeric
  if (!is.numeric(df[[time]])) {
    stop("'time' must be a numeric (integer/double) column.")
  }
  if (anyNA(df[[id]]) || any(!nzchar(as.character(df[[id]]))) ||
      (is.numeric(df[[id]]) && any(!is.finite(df[[id]]))) || any(!is.finite(df[[time]]))) {
    stop("Panel keys must be nonmissing and finite; unit IDs cannot be empty.")
  }
  if (anyDuplicated(df[c(id, time)])) stop("Duplicate id/time cells are not allowed.")

  if (!is.numeric(time_step) || length(time_step) != 1L ||
      !is.finite(time_step) || time_step <= 0) stop("'time_step' must be finite and positive.")
  grid <- (df[[time]] - min(df[[time]])) / time_step
  if (any(abs(grid - round(grid)) > 1e-7)) stop("Time values do not lie on the specified time_step grid.")
  df[[id]] <- as.character(df[[id]])
  df$.csdm_rowid__ <- seq_len(nrow(df))

  o <- order(df[[id]], df[[time]])
  df <- df[o, , drop = FALSE]

  # robust row mapping (do not rely on rownames)
  rownames(df) <- as.character(df$.csdm_rowid__)

  # stable levels for downstream matrix shaping
  attr(df, "csdm_time_step") <- time_step
  attr(df, "csdm_time_levels") <- sort(unique(df[[time]]))
  attr(df, "csdm_id_levels") <- sort(unique(df[[id]]))

  df
}


.csdm_time_index <- function(time_vec, time_step = 1) {
  (time_vec - min(time_vec)) / time_step + 1
}

.csdm_lag <- function(x, time, lag, time_step = 1) {
  grid <- round((time - min(time)) / time_step)
  x[match(grid - lag, grid)]
}
