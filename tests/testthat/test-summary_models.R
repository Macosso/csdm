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
