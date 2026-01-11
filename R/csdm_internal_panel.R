# csdm_internal_panel.R

.csdm_prepare_panel_df <- function(data, id, time) {
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

  if (is.null(id) || is.null(time)) stop("Both 'id' and 'time' must be provided.")
  if (!all(c(id, time) %in% names(df))) stop("'data' must contain 'id' and 'time' columns.")

  # time is required and must be numeric
  if (!is.numeric(df[[time]])) {
    stop("'time' must be a numeric (integer/double) column.")
  }

  df[[id]] <- as.character(df[[id]])

  o <- order(df[[id]], df[[time]])
  df <- df[o, , drop = FALSE]

  # robust row mapping (do not rely on rownames)
  df$.csdm_rowid__ <- seq_len(nrow(df))
  rownames(df) <- as.character(df$.csdm_rowid__)

  # stable levels for downstream matrix shaping
  attr(df, "csdm_time_levels") <- sort(unique(df[[time]]))
  attr(df, "csdm_id_levels") <- sort(unique(df[[id]]))

  df
}


.csdm_time_index <- function(time_vec) {
  tt <- sort(unique(time_vec))
  match(time_vec, tt)
}
