#' @rdname cd_test
#' @param id,time Column names (strings) for the unit and time indexes of a plain
#'   data frame. For a \code{pdata.frame}, these are inferred from its stored
#'   indexes and supplied values are ignored. Time labels need not be numeric.
#' @export
#' @method cd_test data.frame
cd_test.data.frame <- function(object, ..., id = NULL, time = NULL,
                               type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                               n_pc = 4L, seed = NULL, min_overlap = 2L,
                               na.action = c("pairwise", "drop.incomplete.times"),
                               reps = 1L) {
  vars <- .cd_select_variables(as.list(substitute(list(...)))[-1L],
                               names(object), parent.frame())
  panel <- .cd_panel_input(object, id, time, vars)
  type <- match.arg(type)
  na.action <- match.arg(na.action)
  call <- match.call()
  results <- lapply(vars, function(variable) {
    E <- matrix(NA_real_, length(panel$units), length(panel$times),
                dimnames = list(panel$units, as.character(panel$times)))
    E[panel$cells] <- as.numeric(object[[variable]])
    result <- cd_test.default(E, type = type, n_pc = n_pc, seed = seed,
                              min_overlap = min_overlap, na.action = na.action, reps = reps)
    result$variable <- variable
    result$units <- panel$units[!seq_along(panel$units) %in% result$excluded_units]
    result$excluded_unit_ids <- panel$units[result$excluded_units]
    result$call <- call
    result
  })
  names(results) <- vars
  structure(results, class = c("cd_test_list", "list"), call = call,
            id = panel$id, time = panel$time)
}

# The data-frame method reads the pdata.frame indexes without requiring plm.
#' @rdname cd_test
#' @export
#' @method cd_test pdata.frame
cd_test.pdata.frame <- cd_test.data.frame

.cd_select_variables <- function(expressions, columns, envir) {
  if (!length(expressions)) {
    stop("Select at least one variable explicitly in '...'.", call. = FALSE)
  }
  if (!is.null(names(expressions)) && any(nzchar(names(expressions)))) {
    stop("Variable selections in '...' must be unnamed; name test controls explicitly.",
         call. = FALSE)
  }
  selections <- lapply(expressions, function(expr) {
    if (is.symbol(expr) && as.character(expr) %in% columns) {
      return(as.character(expr))
    }
    value <- tryCatch(eval(expr, envir), error = function(e) NULL)
    if (!is.character(value) || !length(value) || anyNA(value) || any(!nzchar(value))) {
      stop("Select variables using bare column names, quoted names, or character vectors.",
           call. = FALSE)
    }
    value
  })
  vars <- unlist(selections, use.names = FALSE)
  if (anyDuplicated(vars)) stop("Select each variable only once.", call. = FALSE)
  unknown <- setdiff(vars, columns)
  if (length(unknown)) stop("Unknown variable(s): ", paste(unknown, collapse = ", "),
                            call. = FALSE)
  vars
}

.cd_panel_input <- function(object, id, time, vars) {
  if (!nrow(object) || anyDuplicated(names(object))) {
    stop("Supply nonempty data with unique column names.", call. = FALSE)
  }
  if (inherits(object, "pdata.frame")) {
    index <- attr(object, "index")
    if (!is.data.frame(index) || ncol(index) < 2L || nrow(index) != nrow(object)) {
      stop("The pdata.frame must have valid stored unit and time indexes.", call. = FALSE)
    }
    id <- names(index)[1L]
    time <- names(index)[2L]
    keys <- index[1:2]
  } else {
    if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id) ||
        !is.character(time) || length(time) != 1L || is.na(time) || !nzchar(time) ||
        id == time || !all(c(id, time) %in% names(object))) {
      stop("'id' and 'time' must name distinct columns in the data frame.", call. = FALSE)
    }
    keys <- object[c(id, time)]
  }
  valid_key <- function(x) {
    is.atomic(x) && is.null(dim(x)) && !is.complex(x) &&
      !anyNA(x) && all(nzchar(as.character(x))) &&
      (!(is.numeric(x) || inherits(x, "Date") || inherits(x, "POSIXt")) ||
         all(is.finite(x)))
  }
  if (!all(vapply(keys, valid_key, logical(1)))) {
    stop("Panel keys must be nonmissing and finite, with no empty labels.", call. = FALSE)
  }
  if (anyDuplicated(keys)) stop("Duplicate id/time cells are not allowed.", call. = FALSE)
  if (any(vars %in% c(id, time))) {
    stop("Select test variables other than the unit and time indexes.", call. = FALSE)
  }
  numeric_column <- function(x) is.numeric(x) && is.null(dim(x))
  if (!all(vapply(object[vars], numeric_column, logical(1)))) {
    stop("Selected test variables must be numeric vectors.", call. = FALSE)
  }
  units <- sort(unique(as.character(keys[[1L]])))
  times <- sort(unique(keys[[2L]]))
  list(id = id, time = time, units = units, times = times,
       cells = cbind(match(as.character(keys[[1L]]), units), match(keys[[2L]], times)))
}

#' @rdname cd_test
#' @export
#' @method print cd_test_list
print.cd_test_list <- function(x, digits = 3, ...) {
  rows <- lapply(names(x), function(variable) {
    result <- x[[variable]]
    selected <- if (identical(result$type, "CDw+")) "CDw_plus" else result$type
    tests <- result$tests[selected]
    data.frame(variable = variable, test = gsub("CDw_plus", "CDw+", names(tests), fixed = TRUE),
               statistic = vapply(tests, function(test) test$statistic, numeric(1)),
               p.value = vapply(tests, function(test) test$p.value, numeric(1)),
               N = result$N, T = result$T, row.names = NULL)
  })
  cat("Cross-sectional dependence tests by variable\n\n")
  table <- do.call(rbind, rows)
  for (column in c("statistic", "p.value")) {
    table[[column]] <- ifelse(is.na(table[[column]]), "NA",
                             formatC(table[[column]], digits = digits, format = "f"))
  }
  print(table, row.names = FALSE, right = TRUE)
  invisible(x)
}
