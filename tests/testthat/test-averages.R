test_that("CSA missingness and leave-one-out honor their contracts", {
  d <- data.frame(id = 1:3, time = 1, x = c(NA, 2, 4))
  expect_equal(cross_sectional_avg(d, "id", "time", "x", leave_out = TRUE)$csa_x, c(3, 4, 2))
  expect_true(all(is.na(cross_sectional_avg(d, "id", "time", "x", na.rm = FALSE)$csa_x)))
  expect_equal(cross_sectional_avg(d, "id", "time", "x", leave_out = TRUE, na.rm = FALSE)$csa_x, c(3, NA, NA))
  expect_equal(cross_sectional_avg(d, "id", "time", "x", weights = c(1, 1, 3))$csa_x, rep(3.5, 3))
  expect_true(all(is.na(cross_sectional_avg(d, "id", "time", "x", weights = rep(0, 3))$csa_x)))
  expect_error(cross_sectional_avg(d, "id", "time", "x", weights = c(1, NA, 1)), "Weights")
  d$csa_x <- 0
  expect_error(cross_sectional_avg(d, "id", "time", "x"), "already exist")
})
