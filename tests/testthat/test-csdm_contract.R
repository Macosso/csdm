test_that("csdm() returns a stable csdm_fit contract", {
  df <- data.frame(
    id = rep(1:4, each = 10),
    time = rep(1:10, times = 4)
  )
  set.seed(1)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 1 + 0.5 * df$x1 - 0.2 * df$x2 + rnorm(nrow(df), sd = 0.5)

  fit_mg <- csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "mg")
  expect_s3_class(fit_mg, "csdm_fit")
  expect_identical(fit_mg$model, "mg")
  expect_true(is.numeric(fit_mg$coef_mg))
  expect_true(is.numeric(fit_mg$se_mg))
  expect_true(is.matrix(fit_mg$vcov_mg))
  expect_true(is.matrix(fit_mg$coef_i))
  expect_true(is.matrix(fit_mg$residuals_e))
  expect_true(all(c("N", "T") %in% names(fit_mg$meta)))

  # S3 methods should not error
  expect_error(print(fit_mg), NA)
  sm <- summary(fit_mg)
  expect_s3_class(sm, "summary.csdm_fit")
  expect_true(all(c("call", "formula", "model", "N", "T", "coef_table") %in% names(sm)))
  expect_error(print(sm), NA)
  expect_true(is.numeric(coef(fit_mg)))
  expect_true(is.matrix(vcov(fit_mg)))
  expect_true(is.matrix(residuals(fit_mg, type = "e")))
  expect_true(is.matrix(predict(fit_mg, type = "residuals")))
  expect_true(is.matrix(predict(fit_mg, type = "xb")))

  # Stable dimensions and names
  expect_identical(rownames(fit_mg$coef_i), as.character(sort(unique(df$id))))
  expect_equal(ncol(fit_mg$residuals_e), length(sort(unique(df$time))))
  expect_identical(colnames(fit_mg$residuals_e), as.character(sort(unique(df$time))))

  fit_cce <- csdm(
    y ~ x1 + x2,
    data = df,
    id = "id",
    time = "time",
    model = "cce",
    csa = csdm_csa(vars = "_all", lags = 0)
  )
  expect_s3_class(fit_cce, "csdm_fit")
  expect_identical(fit_cce$model, "cce")

  fit_dcce <- csdm(
    y ~ x1 + x2,
    data = df,
    id = "id",
    time = "time",
    model = "dcce",
    csa = csdm_csa(vars = "_all", lags = 1)
  )
  expect_s3_class(fit_dcce, "csdm_fit")
  expect_identical(fit_dcce$model, "dcce")

  expect_equal(ncol(fit_dcce$residuals_e), length(sort(unique(df$time))))
  expect_identical(colnames(fit_dcce$residuals_e), as.character(sort(unique(df$time))))

  expect_error(
    csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "cs_ardl"),
    "Not implemented yet"
  )
})


test_that("csdm() residual mapping remains aligned with mid-panel NA", {
  df <- data.frame(
    id = rep(1:4, each = 10),
    time = rep(1:10, times = 4)
  )
  set.seed(2)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 1 + 0.5 * df$x1 - 0.2 * df$x2 + rnorm(nrow(df), sd = 0.5)

  # One missing regressor in the middle for unit 2 at time 5
  df$x1[df$id == 2 & df$time == 5] <- NA_real_

  fit <- csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "mg")
  E <- fit$residuals_e

  expect_true(is.na(E["2", "5"]))
  expect_true(is.finite(E["2", "6"]))
  expect_equal(sum(is.finite(E["2", ])), 9)
})


test_that("dcce supports within-unit y-lags via lr(type='ardl')", {
  df <- data.frame(
    id = rep(1:4, each = 10),
    time = rep(1:10, times = 4)
  )
  set.seed(3)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y  <- 1 + 0.5 * df$x1 - 0.2 * df$x2 + rnorm(nrow(df), sd = 0.5)

  fit <- csdm(
    y ~ x1 + x2,
    data = df,
    id = "id",
    time = "time",
    model = "dcce",
    csa = csdm_csa(vars = "_all", lags = 1),
    lr = csdm_lr(type = "ardl", ylags = 1)
  )

  expect_true("lag1_y" %in% names(coef(fit)))
  expect_true(is.na(fit$residuals_e["1", "1"]))
})
