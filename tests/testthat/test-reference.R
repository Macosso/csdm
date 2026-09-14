test_that("MG and CCE agree with independent unit regressions", {
  d <- panel_fixture()
  averages <- aggregate(cbind(y, x, z) ~ time, d, mean)
  names(averages)[-1L] <- c("cy", "cx", "cz")
  a <- merge(d, averages, by = "time")
  for (model in c("mg", "cce")) {
    formula <- if (model == "mg") y ~ x + z else y ~ x + z + cy + cx + cz
    refs <- lapply(split(a, a$id), function(s) lm(formula, s))
    b <- t(vapply(refs, function(m) coef(m)[1:3], numeric(3)))
    fit <- suppressMessages(csdm(y ~ x + z, d, "id", "time", model = model))
    expect_equal(coef(fit), colMeans(b), tolerance = 1e-10)
    expect_equal(vcov(fit), cov(b) / nrow(b), tolerance = 1e-10)
    for (i in seq_along(refs)) {
      expect_equal(unname(residuals(fit)[i, ]), unname(residuals(refs[[i]])), tolerance = 1e-10)
      expect_equal(unname(predict(fit)[i, ]), unname(fitted(refs[[i]])), tolerance = 1e-10)
    }
  }
})

test_that("CS-ARDL averages independently transformed unit ratios", {
  d <- panel_fixture()
  refs <- lapply(split(d, d$id), function(s) {
    s <- s[order(s$time), ]
    s$ly <- c(NA, head(s$y, -1))
    s$lx <- c(NA, head(s$x, -1))
    s$lz <- c(NA, head(s$z, -1))
    coef(lm(y ~ x + z + ly + lx + lz, s))
  })
  b <- do.call(rbind, refs)
  fit <- suppressMessages(csdm(y ~ x + z, d, "id", "time", model = "cs_ardl",
    csa = csdm_csa("_none"), lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 1)))
  expect_equal(unname(fit$coef_mg), unname(colMeans(b)), tolerance = 1e-10)
  expected <- c(lr_y = mean(b[, "ly"] - 1),
    lr_x = mean((b[, "x"] + b[, "lx"]) / (1 - b[, "ly"])),
    lr_z = mean((b[, "z"] + b[, "lz"]) / (1 - b[, "ly"])))
  expect_equal(tail(coef(fit), 3), expected, tolerance = 1e-10)
})
