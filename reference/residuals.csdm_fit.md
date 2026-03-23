# Extract residual matrix from a fitted csdm model

Returns residuals as an \\N x T\\ matrix (rows are units, columns are
time). This method is designed for panel diagnostics and downstream
tools such as [`cd_test()`](cd_test.md).

## Usage

``` r
# S3 method for class 'csdm_fit'
residuals(object, type = c("e", "u"), ...)
```

## Arguments

- object:

  A fitted object of class `csdm_fit`.

- type:

  Residual type. Currently only `"e"` is implemented.

- ...:

  Currently unused.

## Value

A numeric matrix of residuals with dimensions \\N x T\\.

## See also

[`get_residuals()`](get_residuals.md), [`cd_test()`](cd_test.md),
[`predict.csdm_fit()`](predict.csdm_fit.md)
