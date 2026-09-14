test_that("equal and unequal fixed weights have correct mean covariance", {
  B <- cbind(a = 1:4, b = c(1, 5, 2, 7))
  expect_equal(pooled_vcov(B), cov(B) / 4)
  w <- c(.1, .2, .3, .4)
  centered <- sweep(B, 2, colSums(w * B))
  reference <- Reduce("+", lapply(1:4, function(i) w[i] * tcrossprod(centered[i, ])))
  reference <- reference * sum(w^2) / (1 - sum(w^2))
  dimnames(reference) <- list(colnames(B), colnames(B))
  expect_equal(pooled_vcov(B, w), reference)
  B[1, 2] <- NA
  expect_error(pooled_vcov(B), "Pairwise")
  expect_equal(pooled_vcov(B, pairwise = FALSE), cov(B[-1, ]) / 3)
  expect_error(pooled_vcov(B, rep(0, 4)), "Weights")
})
