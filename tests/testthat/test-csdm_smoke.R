test_that("csdm returns a csdm_fit", {
  set.seed(1)
  df <- expand.grid(
    id = paste0("i", 1:4),
    time = 1:10,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y <- 1 + 0.3 * df$x1 - 0.1 * df$x2 + rnorm(nrow(df), sd = 0.1)

  fit <- csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "mg")
  expect_s3_class(fit, "csdm_fit")
  expect_true(is.matrix(residuals(fit)))
})
