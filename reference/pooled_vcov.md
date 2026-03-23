# Variance-covariance of Mean-Group (MG) averages

Computes the covariance matrix of the MG estimator \\\bar{\beta} =
N^{-1} \sum_i \hat\beta_i\\, using the cross-sectional covariance of
unit-specific slopes and dividing by the effective sample sizes per
coefficient (handles missing entries per unit).

## Usage

``` r
pooled_vcov(beta_i, weights = NULL, pairwise = TRUE)
```

## Arguments

- beta_i:

  Numeric matrix of unit-specific coefficients (\\N x K\\); rows =
  units, columns = coefficients. May contain `NA`s.

- weights:

  Optional numeric vector of length N with nonnegative weights summing
  to 1. If `NULL`, uses equal weights.

- pairwise:

  Logical; use pairwise-complete covariances across units (default
  `TRUE`).

## Value

A `K x K` covariance matrix for the MG mean, with `dimnames` inherited
from `colnames(beta_i)`.

## Details

For equal weights, diagonal entries coincide with
\\\widehat{\mathrm{Var}}(\hat\beta\_{ij}) / N\_{\text{eff},j}\\.
Off-diagonals are scaled analogously using pairwise effective N. If
`weights` are provided, computes \\\mathrm{Var}(\sum_i w_i
\hat\beta_i)\\ using a weighted covariance.
