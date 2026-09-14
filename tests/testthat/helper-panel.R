panel_fixture <- function(n = 8L, periods = 60L) {
  set.seed(914)
  d <- expand.grid(time = seq_len(periods), id = seq_len(n))
  d$x <- rnorm(nrow(d))
  d$z <- rnorm(nrow(d))
  d$y <- 1 + 0.6 * d$x - 0.2 * d$z + rnorm(nrow(d))
  d
}
