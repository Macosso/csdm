
test_that("R2_mg from residual-matrix matches manual computation in unbalanced panel", {
  set.seed(1234)
  df <- data.frame(
    id = rep(1:3, each = 10),
    time = rep(1:10, times = 3)
  )
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 2 + 0.5 * df$x1 - 0.2 * df$x2 + rnorm(nrow(df), sd = 0.7)
  # Drop some cells to create unbalancedness
  df$x1[df$id == 2 & df$time %in% c(2, 3, 4)] <- NA
  df$x2[df$id == 3 & df$time %in% c(7, 8)] <- NA
  df$y[df$id == 1 & df$time == 10] <- NA

  fit <- csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "mg")
  E <- fit$residuals_e
  yname <- "y"
  ids_levels <- rownames(E)
  time_levels <- colnames(E)
  # Strictly align Y to E
  Y <- matrix(NA_real_, nrow = length(ids_levels), ncol = length(time_levels),
              dimnames = list(ids_levels, time_levels))
  ii <- match(as.character(df$id), ids_levels)
  tt <- match(as.character(df$time), time_levels)
  keep <- is.finite(ii) & is.finite(tt)
  if (any(keep)) {
    Y[cbind(ii[keep], tt[keep])] <- as.numeric(df[[yname]][keep])
  }
  # Manual R2_i and R2_mg
  R2_i <- rep(NA_real_, nrow(E)); names(R2_i) <- ids_levels
  for (r in seq_len(nrow(E))) {
    er <- E[r, ]; yr <- Y[r, ]
    ok <- is.finite(er) & is.finite(yr)
    if (sum(ok) < 2L) next
    sse <- sum((er[ok])^2)
    yc <- yr[ok]
    sst <- sum((yc - mean(yc))^2)
    if (!is.finite(sst) || sst <= 0) next
    R2_i[[r]] <- 1 - sse / sst
  }
  R2_mg <- if (all(is.na(R2_i))) NA_real_ else mean(R2_i, na.rm = TRUE)
  # Compare to fit$stats
  expect_equal(fit$stats$R2_i, R2_i, tolerance = 1e-12)
  expect_equal(fit$stats$R2_mg, R2_mg, tolerance = 1e-12)
})
test_that("mg summary includes stats, table columns, and footer lists", {
  df <- data.frame(
    id = rep(1:4, each = 20),
    time = rep(1:20, times = 4)
  )
  set.seed(201)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 1 + 0.6 * df$x1 - 0.1 * df$x2 + rnorm(nrow(df), sd = 0.5)

  fit <- csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "mg")
  out <- utils::capture.output(summary(fit))

  expect_true(any(grepl("R-squared \\(mg\\)", out)))
  expect_true(any(grepl("CD Statistic", out)))
  expect_true(any(grepl("p-value", out)))
  expect_true(any(grepl("Mean Group Variables:", out)))
  expect_true(any(grepl("Cross Sectional Averaged Variables:", out)))

  sm <- summary(fit)
  expect_true(is.list(sm$stats))
  expect_true(is.numeric(sm$stats$cd_stat) || is.na(sm$stats$cd_stat))
  expect_true(is.numeric(sm$stats$cd_p_value) || is.na(sm$stats$cd_p_value))

  expect_true(is.list(sm$tables))
  expect_true("mean_group" %in% names(sm$tables))
  req_cols <- c("Coef.", "Std. Err.", "z", "P>|z|", "CI 2.5%", "CI 97.5%")
  expect_true(all(req_cols %in% colnames(sm$tables$mean_group)))

  expect_true(is.list(sm$lists))
  expect_true(all(c("mean_group_variables", "csa_vars", "csa_lags") %in% names(sm$lists)))
  expect_identical(sm$lists$csa_vars, "none")
})


test_that("R2_mg uses unit-level regression sample (mg)", {
  df <- data.frame(
    id = rep(1:4, each = 20),
    time = rep(1:20, times = 4)
  )
  set.seed(204)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 1 + 0.6 * df$x1 - 0.1 * df$x2 + rnorm(nrow(df), sd = 0.5)

  # Force different per-unit estimation samples
  df$x1[df$id == 2 & df$time %in% c(3, 4, 5, 6)] <- NA_real_
  df$x2[df$id == 3 & df$time %in% c(10, 11)] <- NA_real_

  fit <- csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "mg")
  expect_true(is.list(fit$stats))
  expect_true(!is.null(fit$stats$R2_i))

  # Recompute unit R2 from the exact na.omit regression sample
  expected_r2 <- setNames(rep(NA_real_, length(unique(df$id))), as.character(sort(unique(df$id))))
  for (uid in sort(unique(df$id))) {
    sub <- df[df$id == uid, , drop = FALSE]
    m <- stats::model.frame(y ~ x1 + x2, data = sub, na.action = stats::na.omit)
    if (nrow(m) < 2L) next
    f <- stats::lm(y ~ x1 + x2, data = sub, na.action = stats::na.omit)
    y_used <- stats::model.response(m)
    e_used <- stats::residuals(f)
    ok <- is.finite(y_used) & is.finite(e_used)
    if (sum(ok) < 2L) next
    sse <- sum((e_used[ok])^2)
    yc <- y_used[ok]
    sst <- sum((yc - mean(yc))^2)
    if (!is.finite(sst) || sst <= 0) next
    expected_r2[[as.character(uid)]] <- 1 - sse / sst
  }

  # Compare against stored unit-level R2 (allow for names ordering)
  common <- intersect(names(expected_r2), names(fit$stats$R2_i))
  expect_equal(fit$stats$R2_i[common], expected_r2[common], tolerance = 1e-12)
  expect_equal(fit$stats$R2_mg, mean(fit$stats$R2_i, na.rm = TRUE), tolerance = 1e-12)
})


test_that("cce summary includes CSA footer (lags) and stats", {
  df <- data.frame(
    id = rep(1:4, each = 20),
    time = rep(1:20, times = 4)
  )
  set.seed(202)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 1 + 0.6 * df$x1 - 0.1 * df$x2 + rnorm(nrow(df), sd = 0.5)

  fit <- csdm(
    y ~ x1 + x2,
    data = df,
    id = "id",
    time = "time",
    model = "cce",
    csa = csdm_csa(vars = "_all", lags = 0)
  )

  out <- utils::capture.output(summary(fit))
  expect_true(any(grepl("Cross Sectional Averaged Variables:", out)))

  sm <- summary(fit)
  expect_true(is.list(sm$lists))
  expect_true(sm$lists$csa_vars %in% c("_all", "none") || is.character(sm$lists$csa_vars))
})


test_that("dcce summary includes lag terms in footer when ylags>0", {
  df <- data.frame(
    id = rep(1:4, each = 20),
    time = rep(1:20, times = 4)
  )
  set.seed(203)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 1 + 0.6 * df$x1 - 0.1 * df$x2 + rnorm(nrow(df), sd = 0.5)

  fit <- csdm(
    y ~ x1 + x2,
    data = df,
    id = "id",
    time = "time",
    model = "dcce",
    csa = csdm_csa(vars = "_all", lags = 1),
    lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 0)
  )

  sm <- summary(fit)
  expect_true(is.list(sm$lists))
  expect_true(any(grepl("lag1_y", sm$lists$mean_group_variables)))
})
