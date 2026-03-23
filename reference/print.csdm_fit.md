# Compact print method for fitted csdm models

Prints a concise overview of a fitted `csdm_fit` object, including the
model type, formula, panel dimensions, and a coefficient table with
standard errors when available.

## Usage

``` r
# S3 method for class 'csdm_fit'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A fitted object of class `csdm_fit`.

- digits:

  Number of printed digits.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.

## See also

[`summary.csdm_fit()`](summary.csdm_fit.md),
[`coef.csdm_fit()`](coef.csdm_fit.md),
[`residuals.csdm_fit()`](residuals.csdm_fit.md)
