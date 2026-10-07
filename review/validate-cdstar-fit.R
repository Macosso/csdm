# Selected null calibration for CD* on fitted partial residuals.
dir.create(".dev-check/config", recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CONFIG_DIR = normalizePath(".dev-check/config"))
.libPaths(c(normalizePath(".dev-library"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
set.seed(10072027)
simulations <- 500L
N <- 50L; TT <- 100L
d <- expand.grid(time = seq_len(TT), id = seq_len(N))
loading <- seq(-1, 2, length.out = N)
sigma <- seq(.5, 1.5, length.out = N)
p <- replicate(simulations, {
  factor <- rnorm(TT)
  d$x <- .4 * factor[d$time] + rnorm(nrow(d))
  d$y <- 1 + .6 * d$x + loading[d$id] * factor[d$time] +
    sigma[d$id] * rnorm(nrow(d))
  fit <- csdm(y ~ x, d, "id", "time", model = "cce")
  c(partial = cd_test(fit, type = "CDstar", n_pc = 1L)$tests$CDstar$p.value,
    full = cd_test(residuals(fit), type = "CDstar", n_pc = 1L)$tests$CDstar$p.value)
})
results <- lapply(rownames(p), function(input) {
  rejects <- sum(p[input, ] < .05)
  interval <- binom.test(rejects, simulations)$conf.int
  data.frame(input, seed = 10072027, simulations, N, T = TT, n_pc = 1L,
               rejection = rejects / simulations, lower = interval[1], upper = interval[2])
})
out <- do.call(rbind, results)
print(out, row.names = FALSE)
write.csv(out, "review/cdstar-fit-calibration.csv", row.names = FALSE)
# Broad predeclared smoke bound for the corrected input in this selected design.
# The full-residual variant is retained as a comparison, not an acceptance case.
stopifnot(all(is.finite(p)), out$rejection[out$input == "partial"] < .12)
