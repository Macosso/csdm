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

  expect_error(
    csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "cs_ardl"),
    "Not implemented yet"
  )
})
