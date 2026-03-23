# Extract model coefficients from a fitted csdm model

Returns estimated mean-group coefficients from a `csdm_fit` object. For
`model = "cs_ardl"`, the returned vector includes short-run mean-group
coefficients, the adjustment coefficient (named `lr_<y>`), and long-run
coefficients when available.

## Usage

``` r
# S3 method for class 'csdm_fit'
coef(object, ...)
```

## Arguments

- object:

  A fitted object of class `csdm_fit`.

- ...:

  Currently unused.

## Value

A named numeric vector of estimated coefficients.

## See also

[`summary.csdm_fit()`](summary.csdm_fit.md),
[`vcov.csdm_fit()`](vcov.csdm_fit.md)
