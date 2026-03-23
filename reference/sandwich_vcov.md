# Heteroskedasticity-robust (HC) sandwich variance-covariance for OLS

Computes White/Huber HC0-HC3 sandwich vcov for an OLS design.

## Usage

``` r
sandwich_vcov(X, u, type = c("HC0", "HC1", "HC2", "HC3"))
```

## Arguments

- X:

  Numeric design matrix (n x k) used in OLS.

- u:

  Numeric residual vector (length n).

- type:

  Character; one of `"HC0"`, `"HC1"`, `"HC2"`, `"HC3"`.

## Value

A `k x k` variance-covariance matrix.
