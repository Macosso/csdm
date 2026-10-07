# utils_cd.R - Cross-sectional dependence tests for panel models

#' Cross-sectional dependence (CD) tests for panel data and residuals
#'
#' Computes Pesaran CD, CDw, CDw+, and CD* tests for cross-sectional dependence
#' in panel variables or residuals. The implementation supports data frames,
#' indexed panel data frames, numeric matrices, and fitted \code{csdm_fit} objects.
#'
#' @param object A \code{data.frame}, \code{pdata.frame}, \code{csdm_fit} model,
#'   or numeric matrix with units in rows and time periods in columns (N x T).
#' @param ... For data frames, explicitly selected numeric columns as bare names,
#'   quoted names, or character vectors of names. Selections must be unnamed.
#'   For matrices and fitted models, additional arguments passed to methods.
#'
#' @return An object of class \code{cd_test} with fields \code{tests}, \code{type},
#'   \code{N}, \code{T}, \code{na.action}, \code{excluded_units},
#'   \code{excluded_times}, \code{kept_times}, and \code{call}. The \code{tests}
#'   list contains one or more test results, each with \code{statistic} and
#'   \code{p.value}.
#'   Data-frame methods return a \code{cd_test_list}: a named list containing one
#'   \code{cd_test} result per selected variable, including when only one variable
#'   is selected. Each result also records \code{variable}, \code{units}, and
#'   \code{excluded_unit_ids}. The list has \code{call}, \code{id}, and \code{time}
#'   attributes and prints a combined table with each variable's sample dimensions.
#'   For fitted models, CD* additionally records \code{input = "partial_residuals"}
#'   and its own sample exclusions in \code{tests$CDstar}. With \code{type = "all"},
#'   top-level sample fields describe the full residual panel; the CD* sample
#'   dimensions are recorded in \code{tests$CDstar$N} and \code{tests$CDstar$T}.
#'
#' @details
#' ## Selecting panel variables
#'
#' For a plain data frame, supply distinct unit and time column names through
#' \code{id} and \code{time}. For a \code{pdata.frame}, the stored indexes are
#' used, even when the index columns have been removed from the data. Select at
#' least one numeric non-index column explicitly. Each variable is tested
#' separately using its own available sample and the same missing-data policy.
#' No regression is fitted. CDw demeans the observations within each unit, as it
#' does for matrix inputs. Test controls and panel indexes follow \code{...} in
#' the data-frame methods and must be named.
#'
#' ## Notation
#'
#' Let \eqn{E} contain the selected observations or residuals, with \eqn{N}
#' cross-sectional units and \eqn{T} time periods. For each unit pair
#' \eqn{(i,j)}, let \eqn{T_{ij}} be the number of
#' overlapping time periods and \eqn{\widehat\rho_{ij}} the correlation computed
#' after demeaning both series over that pair's overlapping observations.
#'
#' ## Test statistics
#'
#' \describe{
#'   \item{CD (Pesaran, 2015)}{
#'     \deqn{CD = \sqrt{\frac{2}{N(N-1)}} \sum_{i<j} \sqrt{T_{ij}} \, \widehat\rho_{ij}}
#'     The sum includes pairs with at least \code{min_overlap} observations and
#'     a finite correlation; \eqn{N} remains the retained number of units.
#'   }
#'   \item{CDw (Juodis and Reese, 2022)}{
#'     Independent Rademacher weights \eqn{w_i \in \{-1,1\}} are applied by
#'     unit, independently of the data. Write \eqn{u_{it}=e_{it}-\bar e_i} for
#'     observations demeaned within unit and define the pooled variance
#'     \deqn{\widehat s_W^2=\frac{1}{NT}\sum_{i=1}^N\sum_{t=1}^T w_i^2 u_{it}^2.}
#'     For a balanced panel, equation (30) of Juodis and Reese is implemented as
#'     \deqn{CD_W = (\widehat s_W^2)^{-1}
#'     \sqrt{\frac{2}{TN(N-1)}}
#'     \sum_{t=1}^T\sum_{i=2}^N\sum_{j=1}^{i-1}w_i u_{it}w_j u_{jt}.}
#'     Since \eqn{w_i^2=1}, the pooled variance is \eqn{\operatorname{mean}(u_{it}^2)}.
#'     This is a weighted covariance statistic, with pooled rather than
#'     pair-specific scale normalization. The paper's estimated residuals have
#'     zero unit means when unit intercepts are included. With
#'     \code{reps = G}, draw \eqn{G} independent sets of unit weights and combine
#'     the statistics as in equation (33) of Juodis and Reese:
#'     \deqn{\overline{CD}_W = \frac{1}{\sqrt{G}}\sum_{g=1}^G CD_W^{(g)}.}
#'     The default \code{reps = 1} preserves the single-draw statistic. Use
#'     \code{seed} for reproducibility. The authors suggest a modest number of
#'     draws (for example, 30), since large \eqn{G} can amplify lower-order terms.
#'   }
#'   \item{CDw+ (Juodis and Reese, 2022; Fan, Liao, and Yao, 2015)}{
#'     Equation (32) of Juodis and Reese applies the power-enhancement principle
#'     of Fan, Liao, and Yao. Its statistic is
#'     \deqn{CD_{W+} = CD_W + \sum_{i<j}|\widehat\rho_{ij}|
#'     1\left\{|\widehat\rho_{ij}| > 2\sqrt{\log(N)/T}\right\}.}
#'     Here \eqn{\widehat\rho_{ij}} is the ordinary residual correlation, without
#'     multiplication by \eqn{\sqrt{T}}. The nonnegative screening term is
#'     asymptotically zero under the conditions of the null hypothesis.
#'     With \code{reps > 1}, replace \eqn{CD_W} by \eqn{\overline{CD}_W} and add
#'     the screening term once.
#'   }
#'   \item{CD* (Pesaran-Xie correction; standardized-PCA variant)}{
#'     Each unit is demeaned and divided by its sample standard deviation before
#'     extracting \code{n_pc} principal components. Write \eqn{A} for this
#'     \eqn{T\times N} standardized matrix. The estimated loading matrix
#'     \eqn{\widehat\Gamma} is normalized so that
#'     \eqn{\widehat\Gamma'\widehat\Gamma/N=I}, and the factors are
#'     \eqn{\widehat F=A\widehat\Gamma/N}. For the factor-filtered residuals
#'     \eqn{r_{it}} of \eqn{A-\widehat F\widehat\Gamma'}, define
#'     \deqn{\widehat\sigma_i=\sqrt{T^{-1}\sum_t r_{it}^2},\quad
#'     \widehat\varphi=N^{-1}\sum_i\widehat\gamma_i/\widehat\sigma_i,\quad
#'     \widehat a_i=1-\widehat\sigma_i\widehat\varphi'\widehat\gamma_i.}
#'     The implemented bias correction is
#'     \deqn{\widehat\theta=1-N^{-1}\sum_i\widehat a_i^2,\qquad
#'     CD^*=\frac{CD(r)+\sqrt{T/2}\,\widehat\theta}{1-\widehat\theta}.}
#'     This uses the Pesaran-Xie correction algebra, but their PCA procedure
#'     operates on unstandardized observations. Standardizing before PCA can
#'     change the estimated factor space and is an implementation variant.
#'     Setting \code{n_pc = 0} returns classical CD without factor removal.
#'     For fitted models, the PCA input is
#'     \eqn{\widehat v_{it}=y_{it}-x_{it}'\widehat\beta_i}, where \eqn{x_{it}}
#'     includes all economic and deterministic design columns, including any
#'     intercept, trend, and constructed lags. The fitted CSA contribution is
#'     retained, following the input construction in Pesaran-Xie's regression
#'     procedure. CD, CDw, and CDw+ use the full regression residuals instead.
#'   }
#' }
#'
#' CD* requires a nondegenerate bias-correction denominator. Near-zero
#' denominators can produce severe size distortions, including proportional
#' loading/error-scale designs after standardization. Numerical rank checks
#' do not establish the validity of the asymptotic approximation.
#'
#' All p-values use the two-sided standard-normal approximation
#' \eqn{2\Phi(-|\mathrm{statistic}|)}. The implemented tests do not estimate a
#' serial-correlation correction. Stationarity, time-series dependence, factor
#' strength, and relative panel dimensions must satisfy the relevant theory;
#' accepting a numeric variable does not establish these conditions. For CD,
#' omitting pairs or using very short overlaps can also affect calibration.
#' Classical CD can be biased on CCE residuals because common time parameters
#' have been estimated. CDw addresses this problem under Juodis-Reese's
#' assumptions, including \eqn{\sqrt{T}/N\to0}.
#'
#' ## Missing data and balance
#'
#' Time periods containing no finite residual for any retained unit are outside
#' the effective residual sample and are always removed before balance is
#' assessed. Partially observed periods are handled according to \code{na.action}.
#'
#' \describe{
#'   \item{CD}{Uses pairwise-complete observations by default. Each pairwise
#'   correlation uses available overlaps.}
#'   \item{CDw, CDw+}{Require a balanced sample; explicitly select complete times if desired.}
#'   \item{CD*}{Requires a balanced panel. Explicitly setting \code{na.action = "drop.incomplete.times"}
#'   removes any time period with missing observations. With \code{na.action = "pairwise"},
#'   CD* returns \code{NA} and a warning when missing values are present.}
#' }
#'
#' @references
#' \insertRef{Pesaran2015}{csdm}
#'
#' \insertRef{Pesaran2021}{csdm}
#'
#' \insertRef{JuodisReese2021}{csdm}
#'
#' \insertRef{FanLiaoYao2015}{csdm}
#'
#' \insertRef{PesaranXie2021}{csdm}
#'
#' \insertRef{PesaranXie2026}{csdm}
#'
#' @examples
#' # Simulate independent and dependent panels
#' set.seed(1)
#' E_indep <- matrix(rnorm(100), nrow = 10)
#' E_dep <- matrix(rnorm(10), nrow = 10, ncol = 10, byrow = TRUE)
#'
#' # Compute all tests
#' cd_test(E_indep, type = "all")
#' cd_test(E_dep, type = "CD")
#'
#' # Specific test with parameters
#' cd_test(E_indep, type = "CDstar", n_pc = 2)
#'
#' # Test raw panel variables separately
#' panel <- expand.grid(id = 1:10, year = 1:10)
#' panel$x <- rnorm(nrow(panel))
#' panel$y <- rnorm(nrow(panel))
#' cd_test(panel, x, "y", id = "id", time = "year")
#'
#' # From a fitted csdm model
#' data(PWT_60_07, package = "csdm")
#' df <- PWT_60_07
#' ids <- unique(df$id)[1:10]
#' df_small <- df[df$id %in% ids & df$year >= 1970, ]
#' fit <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small,
#'   id = "id",
#'   time = "year",
#'   model = "cce",
#'   csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"))
#' )
#' cd_test(fit, type = "all")
#'
#' @export
cd_test <- function(object, ...) {
  UseMethod("cd_test")
}

#' @rdname cd_test
#' @param type Which test(s) to compute: one of \code{"CD"}, \code{"CDw"}, \code{"CDw+"},
#'   \code{"CDstar"}, or \code{"all"} (default: \code{"CD"}).
#' @param n_pc Number of principal components for CD* (default 4).
#' @param seed Integer seed for weight draws. Seeded calls restore the caller's RNG state; NULL uses the current RNG stream.
#' @param reps Positive integer number of independent weight draws for CDw and
#'   CDw+ (default 1). The draws are summed and divided by \code{sqrt(reps)}.
#'   CD and CD* are unaffected. Weighted test results record the number of draws
#'   in their \code{reps} field. For data inputs, \code{seed} is applied separately
#'   to each variable's test, making seeded results independent of selection order.
#' @param min_overlap Minimum number of finite observations for retaining a unit
#'   before time filtering (all tests), and minimum overlapping time periods
#'   for including a pair in classical CD (default 2).
#' @param na.action How to handle missing data: \code{"drop.incomplete.times"}
#'   removes time periods with any missing observations to create a balanced panel for CD*;
#'   \code{"pairwise"} (default) uses pairwise correlations for CD; unbalanced CDw/CDw+ requests error and CD* warns.
#' @export
#' @method cd_test default
cd_test.default <- function(object,
                            type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                            n_pc = 4L,
                            seed = NULL,
                            min_overlap = 2L,
                            na.action = c("pairwise", "drop.incomplete.times"),
                            reps = 1L,
                            ...) {
  result <- .cd_test_matrix(object, type = type, n_pc = n_pc, seed = seed,
                            min_overlap = min_overlap, na.action = na.action,
                            reps = reps, ...)
  result$call <- match.call()
  result
}

# The fitted-model method supplies a different panel to CD*, so it can request
# the other tests without first removing factors from the full CCE residuals.
.cd_test_matrix <- function(object,
                            type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                            n_pc = 4L, seed = NULL, min_overlap = 2L,
                            na.action = c("pairwise", "drop.incomplete.times"),
                            reps = 1L, .compute_star = TRUE, ...) {
  type <- match.arg(type)
  na.action <- match.arg(na.action)

  # Validate input
  if (!is.matrix(object) || !is.numeric(object)) {
    stop("cd_test.default: object must be a numeric matrix (N x T).")
  }
  if (nrow(object) < 2 || ncol(object) < 2) {
    stop("cd_test: At least 2 units and 2 time periods required.")
  }

  min_overlap <- .csdm_integer(min_overlap, "min_overlap")
  if (min_overlap < 2L) stop("'min_overlap' must be at least two.")
  reps <- .csdm_integer(reps, "reps")
  if (reps < 1L) stop("'reps' must be at least one.", call. = FALSE)
  if (any(is.infinite(object))) stop("Input values may contain NA, but not infinite values.")
  if (!is.null(seed)) {
    seed <- .csdm_integer(seed, "seed")
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
  }
  usable <- apply(object, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) >= min_overlap && stats::sd(x) > 0
  })
  if (sum(usable) < 2L) stop("At least two nonconstant units with sufficient observations are required.")
  excluded_units <- which(!usable)
  E <- object[usable, , drop = FALSE]

  # Periods with no finite values are not part of the diagnostic sample.
  time_labels <- if (is.null(colnames(E))) seq_len(ncol(E)) else colnames(E)
  empty_times <- colSums(is.finite(E)) == 0L
  excluded_times <- time_labels[empty_times]
  kept_time_indices <- which(!empty_times)
  if (any(empty_times)) E <- E[, !empty_times, drop = FALSE]
  if (ncol(E) < 2L) stop("At least two time periods with residual observations are required.")

  # Handle missing data
  if (na.action == "drop.incomplete.times" && anyNA(E)) {
    # Drop time periods with any missing observations
    complete_times <- colSums(is.na(E)) == 0
    n_dropped <- sum(!complete_times)
    if (sum(complete_times) < 2) {
      stop("cd_test: After dropping incomplete time periods, fewer than 2 periods remain.")
    }
    excluded <- empty_times
    excluded[kept_time_indices[!complete_times]] <- TRUE
    excluded_times <- time_labels[excluded]
    E <- E[, complete_times, drop = FALSE]
    if (n_dropped > 0) {
      message(sprintf("cd_test: Dropped %d incomplete time period%s (%.1f%%). Balanced panel: %d units x %d periods.",
                      n_dropped, if(n_dropped > 1) "s" else "",
                      100 * n_dropped / ncol(object), nrow(E), ncol(E)))
    }
  }

  N <- nrow(E)
  Tt <- ncol(E)

  # Convert to T x N for correlation computation (Stata convention)
  data_tn <- t(E)

  out <- list()

  # 1. CD (classic Pesaran)
  if (type %in% c("CD", "all")) {
    cd_res <- .cd_compute_classic(data_tn, N, Tt, min_overlap)
    out$CD <- list(
      statistic = cd_res$statistic,
      p.value = cd_res$p.value,
      N = N,
      T = Tt,
      pairs_used = cd_res$pairs_used
    )
  }

  if (type %in% c("CDw", "CDw+", "all")) {
    if (anyNA(E)) stop("CDw/CDw+ currently require a balanced sample; use explicit drop.incomplete.times or classical CD.")
    if (!is.null(seed)) set.seed(seed)
    centered <- sweep(data_tn, 2L, colMeans(data_tn))
    variance <- mean(centered^2)
    draws <- vapply(seq_len(reps), function(g) {
      w <- sample(c(-1, 1), N, replace = TRUE)
      weighted <- sweep(centered, 2L, w, "*")
      pair_sum <- sum(rowSums(weighted)^2 - rowSums(weighted^2)) / 2
      sqrt(2 / (Tt * N * (N - 1))) * pair_sum / variance
    }, numeric(1))
    statistic <- sum(draws) / sqrt(reps)
    out$CDw <- list(statistic = statistic,
      p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
      N = N, T = Tt, pairs_used = choose(N, 2), reps = reps)
    if (type %in% c("CDw+", "all")) {
      correlation <- stats::cor(centered)
      rho <- abs(correlation[upper.tri(correlation)])
      threshold <- 2 * sqrt(log(N) / Tt)
      enhancement <- sum(rho[rho > threshold])
      statistic <- statistic + enhancement
      out$CDw_plus <- list(statistic = statistic,
        p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
        N = N, T = Tt, pairs_used = choose(N, 2), reps = reps,
        threshold = threshold, enhancement = enhancement)
    }
  }

  # 4. CD* (bias-corrected with PCA factor removal)
  if (.compute_star && type %in% c("CDstar", "all")) {
    # CD* requires balanced panel
    if (anyNA(E)) {
      warning("cd_test: CD* requires a balanced panel; returning NA due to missing values.")
      out$CDstar <- list(
        statistic = NA_real_,
        p.value = NA_real_,
        N = N,
        T = Tt,
        n_pc = as.integer(n_pc)
      )
    } else {
      cdstar_res <- .cd_compute_star(data_tn, N, Tt, n_pc)
      out$CDstar <- list(
        statistic = cdstar_res$statistic,
        p.value = cdstar_res$p.value,
        N = N,
        T = Tt,
        n_pc = as.integer(n_pc)
      )
    }
  }

  res <- list(
    tests = out,
    type = if (type == "all") names(out) else type,
    N = N,
    T = Tt,
    na.action = na.action,
    excluded_units = excluded_units,
    excluded_times = excluded_times,
    kept_times = colnames(E),
    call = match.call()
  )
  class(res) <- "cd_test"
  res
}

#' @rdname cd_test
#' @export
#' @method cd_test csdm_fit
cd_test.csdm_fit <- function(object,
                              type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                              n_pc = 4L,
                              seed = NULL,
                              min_overlap = 2L,
                              na.action = c("pairwise", "drop.incomplete.times"),
                              reps = 1L,
                              ...) {
  type <- match.arg(type)
  na.action <- match.arg(na.action)
  E <- .csdm_get_residuals(object, type = "auto", strict = TRUE)
  if (type %in% c("CDstar", "all")) {
    V <- .cdstar_fit_input(object, E)
    star <- cd_test.default(V, type = "CDstar", n_pc = n_pc, seed = seed,
                            min_overlap = min_overlap, na.action = na.action,
                            reps = reps, ...)
    star$tests$CDstar$input <- "partial_residuals"
    # CD* can retain different units from tests on the full residuals.
    star$tests$CDstar$excluded_units <- star$excluded_units
    star$tests$CDstar$excluded_times <- star$excluded_times
    star$tests$CDstar$kept_times <- star$kept_times
    if (type == "CDstar") {
      star$call <- match.call()
      return(star)
    }
  }
  result <- .cd_test_matrix(E, type = type, n_pc = n_pc, seed = seed,
                            min_overlap = min_overlap, na.action = na.action,
                            reps = reps, .compute_star = FALSE, ...)
  if (type == "all") {
    result$tests$CDstar <- star$tests$CDstar
    result$type <- names(result$tests)
  }
  result$call <- match.call()
  result
}

# Pesaran-Xie start PCA from y minus the economic/deterministic component,
# retaining the common-factor component rather than using full CCE residuals.
.cdstar_fit_input <- function(object, E) {
  X <- object$model_matrix
  y <- if (!is.null(object$model_frame)) stats::model.response(object$model_frame)
  rows <- object$sample$row[object$sample$used]
  if (!is.matrix(X) || !is.numeric(y) || nrow(X) != length(rows) ||
      length(y) != length(rows) || !is.matrix(object$coef_i)) {
    stop("CD* requires stored economic design and unit coefficients; refit the model.")
  }
  keys <- object$data[rows, c(object$id, object$time), drop = FALSE]
  ids <- as.character(keys[[object$id]])
  partial <- numeric(length(rows))
  for (uid in unique(ids)) {
    selected <- which(ids == uid)
    beta <- object$coef_i[uid, colnames(X)]
    partial[selected] <- y[selected] - as.numeric(X[selected, , drop = FALSE] %*% beta)
  }
  V <- E
  V[] <- NA_real_
  times <- keys[[object$time]]
  if (isTRUE(object$meta$pdata) && is.factor(times)) {
    times <- suppressWarnings(as.numeric(as.character(times)))
  }
  cells <- cbind(match(ids, rownames(V)),
                 match(as.character(times), colnames(V)))
  if (anyNA(cells)) stop("Stored CD* sample indexes are inconsistent; refit the model.")
  V[cells] <- partial
  V
}

#' @rdname cd_test
#' @param x An object of class \code{cd_test} or \code{cd_test_list}.
#' @param digits Number of digits to print (default 3).
#' @export
#' @method print cd_test
print.cd_test <- function(x, digits = 3, ...) {
  if (is.null(x$tests) || length(x$tests) == 0) {
    cat("cd_test: no results\n")
    return(invisible(x))
  }

  cat("Cross-sectional dependence tests\n")
  different_samples <- any(vapply(x$tests, function(test) {
    !is.null(test$N) && !is.null(test$T) && (test$N != x$N || test$T != x$T)
  }, logical(1)))
  if (!different_samples && !is.null(x$N) && !is.null(x$T)) {
    cat(sprintf("N = %d, T = %d\n", x$N, x$T))
  }
  cat("\n")

  tests <- x$tests
  type_all <- is.character(x$type) && length(x$type) > 1
  if (!type_all && is.character(x$type) && length(x$type) == 1) {
    type_all <- identical(x$type, "all")
  }
  type_map <- c(
    "CD" = "CD",
    "CDw" = "CDw",
    "CDw+" = "CDw_plus",
    "CDstar" = "CDstar"
  )

  get_num <- function(test, key) {
    if (!is.null(test[[key]])) as.numeric(test[[key]]) else NA_real_
  }

  if (type_all) {
    test_names <- names(tests)
    display_names <- gsub("CDw_plus", "CDw+", test_names, fixed = TRUE)
    stat <- vapply(tests, get_num, key = "statistic", FUN.VALUE = numeric(1))
    pval <- vapply(tests, get_num, key = "p.value", FUN.VALUE = numeric(1))
    out <- data.frame(
      statistic = stat,
      p.value = pval,
      row.names = display_names,
      stringsAsFactors = FALSE
    )
    if (different_samples) {
      out$N <- vapply(tests, get_num, key = "N", FUN.VALUE = numeric(1))
      out$T <- vapply(tests, get_num, key = "T", FUN.VALUE = numeric(1))
    }
  } else {
    key <- type_map[[as.character(x$type[1])]]
    if (is.null(key) || is.null(tests[[key]])) {
      key <- names(tests)[1]
    }
    cd_test <- tests[[key]]
    test_name <- as.character(x$type[1])
    out <- data.frame(
      statistic = get_num(cd_test, "statistic"),
      p.value = get_num(cd_test, "p.value"),
      row.names = test_name,
      stringsAsFactors = FALSE
    )
  }

  fmt_num <- function(x) {
    ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
  }
  out$statistic <- fmt_num(out$statistic)
  out$p.value <- fmt_num(out$p.value)
  print(out, row.names = TRUE, right = TRUE)
  invisible(x)
}

# Internal helper: compute classic CD statistic with pairwise overlap
.cd_compute_classic <- function(data_tn, N, Tt, min_overlap = 2L) {
  # data_tn: T x N matrix
  # Returns: list(statistic, p.value, pairs_used)

  corr_mat <- stats::cor(data_tn, use = "pairwise.complete.obs")

  # Count valid pairs (those with sufficient overlap)
  pairs_used <- 0
  cd_sum <- 0

  for (i in seq_len(N - 1)) {
    for (j in (i + 1):N) {
      ok <- is.finite(data_tn[, i]) & is.finite(data_tn[, j])
      T_ij <- sum(ok)
      if (T_ij >= min_overlap) {
        rho_ij <- corr_mat[i, j]
        if (is.finite(rho_ij)) {
          cd_sum <- cd_sum + sqrt(T_ij) * rho_ij
          pairs_used <- pairs_used + 1
        }
      }
    }
  }

  if (pairs_used == 0) {
    warning("cd_test: No unit pairs have sufficient overlap for CD computation.")
    return(list(statistic = NA_real_, p.value = NA_real_, pairs_used = 0L))
  }

  # The Pesaran normalization uses the retained number of units.
  cd_stat <- sqrt(2 / (N * (N - 1))) * cd_sum
  cd_p <- 2 * stats::pnorm(abs(cd_stat), lower.tail = FALSE)

  list(statistic = cd_stat, p.value = cd_p, pairs_used = as.integer(pairs_used))
}

# Internal helper: compute CD* with PCA factor removal
.cd_compute_star <- function(data_tn, N, Tt, n_pc) {
  n_pc <- .csdm_integer(n_pc, "n_pc")
  if (n_pc >= min(N, Tt - 1L)) stop("'n_pc' must be below min(N, T - 1).")
  scales <- apply(data_tn, 2L, stats::sd)
  if (any(!is.finite(scales) | scales <= 0)) stop("CD* requires nonconstant complete unit series.")
  data_std <- sweep(sweep(data_tn, 2L, colMeans(data_tn)), 2L, scales, "/")
  if (n_pc == 0L) return(.cd_compute_classic(data_std, N, Tt))
  decomposition <- svd(data_std, nu = n_pc, nv = 0L)
  if (decomposition$d[n_pc] <= decomposition$d[1L] * 1e-7) stop("Requested factors exceed numerical rank.")
  factors <- cbind(1, decomposition$u[, seq_len(n_pc), drop = FALSE])
  beta <- qr.coef(qr(factors), data_std)
  residual <- data_std - factors %*% beta
  sigma <- sqrt(colMeans(residual^2))
  if (any(!is.finite(sigma) | sigma <= sqrt(.Machine$double.eps))) stop("CD* residual scales are degenerate.")
  loadings <- beta[-1L, , drop = FALSE]
  gamma <- sweep(loadings, 1L, sqrt(rowMeans(loadings^2)), "/")
  phi <- rowMeans(sweep(gamma, 2L, sigma, "/"))
  a <- as.numeric(1 - t(sweep(gamma, 2L, sigma, "*")) %*% phi)
  correction <- mean(a^2)
  if (!is.finite(correction) || correction <= sqrt(.Machine$double.eps)) stop("CD* bias correction is degenerate.")
  cd <- .cd_compute_classic(residual, N, Tt)$statistic
  statistic <- (cd + sqrt(Tt / 2) * (1 - correction)) / correction
  list(statistic = statistic, p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE))
}
