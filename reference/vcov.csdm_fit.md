# Extract coefficient covariance matrix from a fitted csdm model

Extract coefficient covariance matrix from a fitted csdm model

## Usage

``` r
# S3 method for class 'csdm_fit'
vcov(object, ...)
```

## Arguments

- object:

  A fitted object of class `csdm_fit`.

- ...:

  Currently unused.

## Value

A numeric variance-covariance matrix aligned with `coef(object)` for
models where this is available.

## See also

[`coef.csdm_fit()`](coef.csdm_fit.md),
[`summary.csdm_fit()`](summary.csdm_fit.md)
