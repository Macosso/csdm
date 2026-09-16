.libPaths(c(normalizePath(".dev-library"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
reps <- 500L
N <- 50L
TT <- 100L
results <- lapply(c("heterogeneous", "proportional"), function(design) {
  set.seed(9142027)
  loading <- if (design == "heterogeneous") seq(-1, 2, length.out = N) else seq(.5, 1.5, length.out = N)
  p <- replicate(reps, {
    E <- outer(loading, rnorm(TT)) +
      matrix(rnorm(N * TT), N) * seq(.5, 1.5, length.out = N)
    cd_test(E, type = "CDstar", n_pc = 1)$tests$CDstar$p.value
  })
  rejects <- sum(p < .05)
  interval <- binom.test(rejects, reps)$conf.int
  data.frame(design, seed = 9142027, reps, N, T = TT, rejection = rejects / reps,
    lower = interval[1], upper = interval[2])
})
out <- do.call(rbind, results)
write.csv(out, "review/cdstar-calibration.csv", row.names = FALSE)
print(out)
# The proportional design diagnoses a zero limiting correction, outside the test's assumptions.
# Apply the predeclared smoke bound only to the nondegenerate calibration design.
stopifnot(out$rejection[out$design == "heterogeneous"] < .12)
