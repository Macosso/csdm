test_that("evaluated responses and dot expansion use user data", {
  d <- panel_fixture()
  a <- suppressMessages(csdm(I(y^2) ~ x + z, d, "id", "time"))
  d$yy <- d$y^2
  b <- suppressMessages(csdm(yy ~ x + z, d, "id", "time"))
  expect_equal(coef(a), coef(b))
  expect_equal(a$stats$R2_i, b$stats$R2_i)
  d$yy <- NULL
  a <- suppressMessages(csdm(y ~ . - id - time, d, "id", "time"))
  expect_false(any(grepl(".csdm_", names(coef(a)), fixed = TRUE)))
})

test_that("subset and NA handling affect the actual sample", {
  d <- panel_fixture()
  a <- suppressMessages(csdm(y ~ x, d, "id", "time", subset = time > 10))
  b <- suppressMessages(csdm(y ~ x, d[d$time > 10, ], "id", "time"))
  expect_equal(coef(a), coef(b))
  expect_true(all(a$sample$row %in% which(d$time > 10)))
  d$x[30] <- NA
  expect_error(csdm(y ~ x, d, "id", "time", na.action = na.fail), "Missing")
})

test_that("default CSA averages transformed design variables", {
  d <- panel_fixture()
  a <- suppressMessages(csdm(y ~ I(x^2), d, "id", "time", model = "cce"))
  d$xx <- d$x^2
  b <- suppressMessages(csdm(y ~ xx, d, "id", "time", model = "cce"))
  expect_equal(unname(coef(a)), unname(coef(b)))
  expect_identical(a$meta$csa$resolved_vars, c("y", "I(x^2)"))
  d$f <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  expect_no_error(suppressMessages(csdm(y ~ x * f, d, "id", "time", model = "cce")))
})
