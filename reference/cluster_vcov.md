# Cluster-robust variance-covariance for OLS

Computes one- or two-way cluster-robust vcov for an OLS design using the
Liang-Zeger "meat" and Cameron-Gelbach-Miller inclusion-exclusion for
two-way clustering.

## Usage

``` r
cluster_vcov(X, u, cluster, df_correction = TRUE, type = c("oneway", "twoway"))
```

## Arguments

- X:

  Numeric design matrix (n x k) used in OLS.

- u:

  Numeric residual vector (length n).

- cluster:

  One of:

  - a vector (length n) of cluster ids for one-way clustering; or

  - a data.frame/list with two vectors (each length n) for two-way
    clustering.

- df_correction:

  Logical; apply small-sample corrections. Default `TRUE`.

- type:

  Character, one of `"oneway"` or `"twoway"`.

## Value

A `k x k` variance-covariance matrix.
