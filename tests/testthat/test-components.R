test_that("CS-ARDL coefficients and full covariance align", {
  d <- panel_fixture()
  a <- csdm(y ~ x + z, d, "id", "time", model = "cs_ardl",
    csa = csdm_csa("_none"), lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 1))
  for (block in names(a$components)) {
    b <- coef(a, component = block); V <- vcov(a, component = block)
    expect_identical(names(b), colnames(V))
    B <- a$components[[block]]$unit_coefficients
    expect_equal(b, colMeans(B))
    expect_equal(V, cov(B) / nrow(B))
  }
  expect_true(all(is.finite(confint(a))))
  expect_length(a$cs_ardl$ar_roots, 8)
})

test_that("component samples are explicit and consistent", {
  B <- cbind(a = 1:4, ratio = c(1, NA, 3, 4))
  all <- csdm:::.csdm_parameter_component(B)
  level <- csdm:::.csdm_parameter_component(B[, "a", drop = FALSE])
  expect_equal(all$n_used, 3)
  expect_equal(level$n_used, 4)
  expect_equal(all$coefficients, colMeans(B[c(1, 3, 4), ]))
  expect_equal(all$vcov, cov(B[c(1, 3, 4), ]) / 3)
})
