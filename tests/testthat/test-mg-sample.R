test_that("MG mean and covariance use exactly the same units", {
  d <- panel_fixture(); d$z[d$id == 1] <- 0
  a <- suppressWarnings(suppressMessages(csdm(y ~ x + z, d, "id", "time")))
  B <- a$coef_i[a$meta$included_units, , drop = FALSE]
  expect_equal(coef(a), colMeans(B))
  expect_equal(vcov(a), cov(B) / nrow(B))
  expect_equal(a$meta$N_observed, 8)
  expect_equal(a$meta$N_used, 7)
  expect_error(csdm(y ~ x + z, d[d$id == 2, ], "id", "time"), "At least two")
})
