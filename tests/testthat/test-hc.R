test_that("HC0 through HC3 match sandwich and are sign-invariant", {
  skip_if_not_installed("sandwich")
  d <- panel_fixture()
  m <- lm(y ~ x + z, d)
  X <- model.matrix(m); u <- residuals(m)
  for (type in c("HC0", "HC1", "HC2", "HC3")) {
    V <- sandwich_vcov(X, u, type)
    expect_equal(V, sandwich::vcovHC(m, type = type), tolerance = 1e-9)
    expect_equal(V, sandwich_vcov(X, -u, type))
  }
  expect_equal(unname(sandwich_vcov(matrix(1, 4, 1), c(-1, 1, -1, 1))), matrix(.25))
  expect_error(sandwich_vcov(cbind(X, X[, 2]), u), "full-rank")
})
