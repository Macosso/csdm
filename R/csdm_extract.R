#' Extract fitted model data and sample information
#'
#' The model frame and economic design matrix contain only estimated observations,
#' in panel order. formula() and terms() include constructed economic lags.
#' nobs() counts estimated observations. df.residual() returns Inf for asymptotic
#' normal MG inference; individual regression degrees of freedom are in object$units.
#' @param object,x,formula A csdm_fit object.
#' @param ... Further arguments.
#' @name csdm_extract
NULL

#' @rdname csdm_extract
#' @export
nobs.csdm_fit <- function(object, ...) sum(object$sample$used)

#' @rdname csdm_extract
#' @export
model.frame.csdm_fit <- function(formula, ...) formula$model_frame

#' @rdname csdm_extract
#' @export
model.matrix.csdm_fit <- function(object, ...) object$model_matrix

#' @rdname csdm_extract
#' @export
terms.csdm_fit <- function(x, ...) x$terms

#' @rdname csdm_extract
#' @export
formula.csdm_fit <- function(x, ...) x$fitted_formula

#' @rdname csdm_extract
#' @export
df.residual.csdm_fit <- function(object, ...) Inf

#' Extract fitted panel values
#' @param object A csdm_fit object.
#' @param format Matrix (units by times), vector (original rows), or long data.
#'   Vector and long formats pad excluded/unestimated rows with NA.
#' @param ... Further arguments.
#' @export
fitted.csdm_fit <- function(object, format = c("matrix", "vector", "long"), ...) {
  .csdm_observation_output(object, "fitted", match.arg(format))
}

.csdm_observation_output <- function(object, value, format) {
  if (format == "matrix") return(if (value == "fitted") object$fitted_xb else object$residuals_e)
  out <- rep(NA_real_, nrow(object$data))
  out[object$sample$row] <- object$sample[[value]]
  if (format == "vector") {
    names(out) <- rownames(object$data)
    return(out)
  }
  data <- object$data[c(object$id, object$time)]
  data$.row <- seq_len(nrow(data))
  data[[value]] <- out
  data
}

#' Update a fitted panel model
#' @param object A csdm_fit object.
#' @param formula. Formula update.
#' @param ... Named arguments replacing the original call arguments.
#' @param evaluate Evaluate the updated call.
#' @export
update.csdm_fit <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  call[[1L]] <- quote(csdm::csdm)
  call$data <- object$data
  if (isTRUE(object$meta$pdata)) {
    call$data[[object$time]] <- as.numeric(as.character(call$data[[object$time]]))
  }
  call$subset <- seq_len(nrow(object$data)) %in% object$meta$selected_rows
  call$na.action <- object$meta$na_action
  call$time_step <- object$meta$time_step
  call$model <- object$model
  call$fullsample <- object$meta$fullsample
  call$mgmissing <- object$meta$mgmissing
  call$formula <- object$formula
  call$id <- object$id
  call$time <- object$time
  for (nm in c("lr", "pooled", "vcov")) call[[nm]] <- object$meta[[nm]]
  call$csa <- object$meta$requested_csa
  call$trend <- object$meta$trend
  if (!missing(formula.)) call$formula <- stats::update.formula(object$formula, formula.)
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    if (is.null(names(extras)) || any(!nzchar(names(extras)))) stop("Update arguments must be named.")
    for (nm in names(extras)) call[[nm]] <- extras[[nm]]
  }
  if (evaluate) eval(call, parent.frame()) else call
}
