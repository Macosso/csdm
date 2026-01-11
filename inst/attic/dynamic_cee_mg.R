## Legacy implementation preserved for reference.
## Stored under inst/attic and NOT loaded by the package.

#' Dynamic Common Correlated Effects Mean Group (CCE-MG) Estimator
#'
#' @description
#' Levels ARDL(1) with Common Correlated Effects (CCE) augmentation and
#' Mean Group (MG) averaging. Includes p lags of cross-sectional averages (CSAs)
#' of y and regressors. Designed for heterogeneous panels with latent common factors
#' (Chudik & Pesaran; Pesaran & Smith).
#'
#' @param formula y ~ x1 + x2 + ...
#' @param data data.frame, tibble, or `plm::pdata.frame`.
#' @param id,time Character names of unit and time identifiers (required if not a `pdata.frame`).
#' @param p Integer, number of CSA lags (default 3).
#' @param na.action NA handler (default `na.omit`).
#' @param vcov_unit One of c("OLS","HC"); per-unit covariance for reporting only.
#' @param cd_method One of c("CDstar","CD","none"); residual cross-section test to run in summary.
#' @details
#' (Documentation preserved from legacy version.)
#' @noRd
dynamic_cce_mg <- function(formula, data, id = NULL, time = NULL,
													 p = 3, na.action = na.omit,
													 vcov_unit = c("OLS","HC"),
													 cd_method = c("CDstar","CD","none")) {
	.Deprecated("csdm", package = "csdm")
	vcov_unit <- match.arg(vcov_unit)
	cd_method <- match.arg(cd_method)

	# ---- indices ----
	if (inherits(data, "pdata.frame")) {
		idx <- names(plm::index(data)); id <- idx[1]; time <- idx[2]
		df <- as.data.frame(data)
	} else {
		if (is.null(id) || is.null(time)) stop("Provide 'id' and 'time' when data is not a pdata.frame.")
		df <- as.data.frame(data)
	}

	# ---- variables ----
	mf <- model.frame(formula, data = df)
	vars <- all.vars(formula)
	y <- vars[1]; x_vars <- vars[-1]
	if (length(x_vars) == 0L) stop("At least one regressor is required on the RHS.")

	# keep only needed cols
	df <- df |>
		dplyr::select(!!sym(id), !!sym(time), dplyr::all_of(c(y, x_vars)))

	# ---- build CSAs at time-level (correct) ----
	# time-level table of CSAs
	csa_time <- df |>
		dplyr::group_by(!!sym(time)) |>
		dplyr::summarise(dplyr::across(dplyr::all_of(c(y, x_vars)),
																	 ~ mean(.x, na.rm = TRUE),
																	 .names = "cs_{.col}"),
										 .groups = "drop") |>
		dplyr::arrange(!!sym(time))

	# add p lags of CSAs by time only
	add_lags <- function(tab, vnames, L) {
		for (l in seq_len(L)) {
			for (v in vnames) {
				tab[[paste0(v, "_lag", l)]] <- dplyr::lag(tab[[v]], l)
			}
		}
		tab
	}
	csa_cols <- names(csa_time)[startsWith(names(csa_time), "cs_")]
	if (p > 0L) csa_time <- add_lags(csa_time, csa_cols, p)

	# ---- merge CSAs back to panel ----
	df <- dplyr::left_join(df, csa_time, by = time)

	# ---- add ARDL(1) lag of y by unit ----
	df <- df |>
		dplyr::arrange(!!sym(id), !!sym(time)) |>
		dplyr::group_by(!!sym(id)) |>
		dplyr::mutate(!!paste0("lag_", y) := dplyr::lag(!!sym(y), 1L)) |>
		dplyr::ungroup()

	# regressors = lag_y + own x + contemporaneous CSAs + CSA lags
	avg_names <- csa_cols
	lag_labs  <- if (p > 0L) {
		as.vector(outer(avg_names, seq_len(p), FUN = function(v, l) paste0(v, "_lag", l)))
	} else character(0)
	regressors <- unique(c(paste0("lag_", y), x_vars, avg_names, lag_labs))

	# filter rows with complete cases for y and regressors (use provided na.action)
	df2 <- na.action(df[, c(id, time, y, regressors)])

	# ---- nest by unit and fit ----
	nested <- df2 |>
		tidyr::nest(.by = !!sym(id)) |>
		dplyr::rename(.unit_id = !!sym(id)) # preserve id column

	fit_formula <- stats::as.formula(paste(y, "~", paste(regressors, collapse = " + ")))

	get_vcov <- function(fit) {
		if (vcov_unit == "HC") {
			if (!requireNamespace("sandwich", quietly = TRUE))
				stop("Install 'sandwich' for HC vcov.")
			sandwich::vcovHC(fit, type = "HC1")
		} else {
			stats::vcov(fit)
		}
	}

	nested <- df2 |>
		tidyr::nest(.by = !!sym(id)) |>
		dplyr::rename(.unit_id = !!sym(id)) |>
		dplyr::mutate(
			results = purrr::map2(data, .unit_id, ~ {
				di  <- .x
				uid <- .y

				cols   <- c(y, regressors)
				mask   <- stats::complete.cases(di[, cols, drop = FALSE])
				di_fit <- di[mask, , drop = FALSE]

				if (nrow(di_fit) == 0L) {
					return(list(
						model = NULL,
						coef  = stats::setNames(numeric(0), character(0)),
						vc    = NULL,
						r2    = NA_real_,
						residuals = stats::setNames(
							data.frame(numeric(0), numeric(0), numeric(0)),
							c(id, time, "residual")
						)
					))
				}

				fit <- stats::lm(fit_formula, data = di_fit)
				r   <- stats::residuals(fit)

				# build residual frame using the group id (uid), not di_fit[[id]]
				df_res <- data.frame(uid, di_fit[[time]], r, check.names = FALSE)
				names(df_res) <- c(id, time, "residual")

				list(
					model = fit,
					coef  = stats::coef(fit),
					vc    = get_vcov(fit),
					r2    = summary(fit)$r.squared,
					residuals = df_res
				)
			}),
			model     = purrr::map(results, "model"),
			coef      = purrr::map(results, "coef"),
			vc        = purrr::map(results, "vc"),
			r2        = purrr::map_dbl(results, "r2"),
			residuals = purrr::map(results, "residuals")
		) |>
		dplyr::select(-results)

	# ---- align coefficients across units ----
	all_terms <- unique(unlist(purrr::map(nested$coef, names)))
	coef_mat <- do.call(rbind, lapply(nested$coef, function(b) {
		out <- stats::setNames(rep(NA_real_, length(all_terms)), all_terms)
		out[names(b)] <- b
		out
	}))
	rownames(coef_mat) <- nested$.unit_id

	# MG average on the "MG terms": intercept, lag_y, and user regressors only (exclude CSA controls)
	mg_terms <- intersect(c("(Intercept)", paste0("lag_", y), x_vars), colnames(coef_mat))
	B <- coef_mat[, mg_terms, drop = FALSE]
	b_mg <- colMeans(B, na.rm = TRUE)

	# Cross-section dispersion VCOV (Pesaran-Smith MG)
	complete_rows <- stats::complete.cases(B)
	Bc <- B[complete_rows, , drop = FALSE]
	N_eff <- nrow(Bc)
	if (N_eff < 2L) stop("Not enough complete unit estimates to compute MG dispersion.")
	dev <- sweep(Bc, 2L, b_mg, "-")
	V_mg <- crossprod(as.matrix(dev)) / (N_eff * (N_eff - 1))
	se_mg <- sqrt(diag(V_mg))

	# residuals and R2
	res_df <- tidyr::unnest(nested, residuals)
	r2 <- stats::setNames(nested$r2, nested$.unit_id)

	out <- list(
		coef      = b_mg,
		se        = stats::setNames(se_mg, names(b_mg)),
		vcov_mg   = V_mg,
		coef_i    = stats::setNames(split(coef_mat, row(coef_mat)), NULL),
		vcov_i    = stats::setNames(nested$vc, nested$.unit_id),
		coef_df   = data.frame(unit = nested$.unit_id, B, check.names = FALSE),
		r.squared = r2,
		residuals = plm::pdata.frame(res_df, index = c(id, time)),
		call      = formula,
		p         = p,
		id        = id,
		time      = time,
		meta      = list(
			N = length(unique(df2[[id]])),
			T = length(unique(df2[[time]])),
			regressors = regressors,
			csa_vars = c(y, x_vars),
			csa_lags = p,
			cd_method = cd_method,
			vcov_unit = vcov_unit
		)
	)
	class(out) <- "dynamic_cce_mg"
	out
}

print.dynamic_cce_mg <- function(x, ...) {
	cat("Dynamic CCE-MG (ARDL(1) + CCE controls)\n")
	cat("Call:", deparse(x$call), "\n")
	cat("CSA lags (p):", x$p, "\n\n")
	cat("Mean-group coefficients:\n")
	print(round(x$coef, 6))
	cat("\nStd. errors (cross-section dispersion):\n")
	print(round(x$se, 6))
	invisible(x)
}

summary.dynamic_cce_mg <- function(object, ...) {
	mg <- data.frame(
		estimate = object$coef,
		std.error = object$se,
		z = object$coef / object$se,
		row.names = names(object$coef)
	)
	mg$p.value <- 2 * (1 - pnorm(abs(mg$z)))
	mg$significance <- cut(mg$p.value,
												 breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
												 labels = c("***","**","*","."," "),
												 right = FALSE)

	CD <- NULL
	if (object$meta$cd_method != "none") {
		if (object$meta$cd_method == "CDstar" && exists("cd_star")) {
			CD <- cd_star(object$residuals, var = "residual")
		} else if (object$meta$cd_method == "CD" && exists("cd")) {
			CD <- cd(object$residuals, var = "residual")
		}
	}

	out <- list(
		call = object$call,
		coefficients = mg,
		r.squared = mean(object$r.squared, na.rm = TRUE),
		CD = CD,
		meta = object$meta
	)
	class(out) <- "summary.dynamic_cce_mg"
	out
}

print.summary.dynamic_cce_mg <- function(x, digits = 6, ...) {
	cat("Dynamic CCE-MG (heterogeneous ARDL(1) with CCE)\n")
	cat("Call:", deparse(x$call), "\n\n")

	x$coefficients <- x$coefficients |>
		dplyr::mutate(
			estimate = formatC(estimate, format = "f", digits = digits),
			std.error = formatC(std.error, format = "f", digits = digits),
			z = formatC(z, format = "f", digits = digits),
			p.value = formatC(p.value, format = "f", digits = 4)
		)
	print(x$coefficients)
	print("")
	cat("\nMean R^2 across units:", format(round(x$r.squared, 3), nsmall = 3), "\n")
	if (!is.null(x$CD)) {
		cat("Cross-section dependence test:", x$CD$method, "\n",
				"stat =", round(x$CD$statistic, 2), " p-value =", round(x$CD$p.value, 4), "\n")
	}
	invisible(x)
}

vcov.dynamic_cce_mg <- function(object, ...) {
	object$vcov_mg
}

coef.dynamic_cce_mg <- function(object, ...) {
	object$coef
}

