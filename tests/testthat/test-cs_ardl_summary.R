test_that("cs_ardl summary prints CS-ARDL blocks", {
  df <- data.frame(
    id = rep(1:4, each = 30),
    time = rep(1:30, times = 4)
  )
  set.seed(101)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))

  # Build a persistent y within each unit
  df$y <- NA_real_
  for (uid in unique(df$id)) {
    idx <- which(df$id == uid)
    idx <- idx[order(df$time[idx])]
    yv <- numeric(length(idx))
    yv[1] <- rnorm(1)
    for (t in 2:length(idx)) {
      yv[t] <- 0.7 * yv[t - 1] + 0.3 * df$x1[idx[t]] - 0.2 * df$x2[idx[t]] + rnorm(1, sd = 0.3)
    }
    df$y[idx] <- yv
  }

  fit <- csdm(
    y ~ x1 + x2,
    data = df,
    id = "id",
    time = "time",
    model = "cs_ardl",
    csa = csdm_csa(vars = "_all", lags = 1),
    lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 1)
  )

  out <- utils::capture.output(summary(fit))
  expect_true(any(grepl("Short Run Est\\.", out)))
  expect_true(any(grepl("Adjust\\. Term", out)))
  expect_true(any(grepl("Long Run Est\\.", out)))

  sm <- summary(fit)
  expect_true(is.list(sm$cs_ardl))
  expect_true(all(c("short_run", "adjust", "long_run") %in% names(sm$cs_ardl)))
})


test_that("cs_ardl coef() includes lr_ terms and uses xtdcce2 sign convention", {
  df <- data.frame(
    id = rep(1:4, each = 30),
    time = rep(1:30, times = 4)
  )
  set.seed(102)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))

  df$y <- NA_real_
  for (uid in unique(df$id)) {
    idx <- which(df$id == uid)
    idx <- idx[order(df$time[idx])]
    yv <- numeric(length(idx))
    yv[1] <- rnorm(1)
    for (t in 2:length(idx)) {
      yv[t] <- 0.8 * yv[t - 1] + 0.2 * df$x1[idx[t]] + rnorm(1, sd = 0.25)
    }
    df$y[idx] <- yv
  }

  fit <- csdm(
    y ~ x1 + x2,
    data = df,
    id = "id",
    time = "time",
    model = "cs_ardl",
    csa = csdm_csa(vars = "_all", lags = 1),
    lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 1)
  )

  cf <- coef(fit)
  expect_true("lr_y" %in% names(cf))
  expect_true(all(paste0("lr_", c("x1", "x2")) %in% names(cf)))

  # Adjustment term uses xtdcce2 sign convention: -(1 - sum(alpha_y_lags))
  expect_equal(
    unname(cf[["lr_y"]]),
    unname(fit$cs_ardl$mg$adj[["estimate"]])
  )

  # Pragmatic sign check: with persistence < 1, denom > 0 so lr_y should be negative
  expect_true(is.finite(cf[["lr_y"]]))
  expect_true(cf[["lr_y"]] < 0)
})
