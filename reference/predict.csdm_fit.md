# Predict method for csdm models

Produces fitted values (index `"xb"`) when available, or returns model
residuals. Prediction on new data is not yet implemented.

## Usage

``` r
# S3 method for class 'csdm_fit'
predict(object, newdata = NULL, type = c("xb", "residuals"), ...)
```

## Arguments

- object:

  A fitted object of class `csdm_fit`.

- newdata:

  Optional new data (not yet supported).

- type:

  One of `"xb"` for fitted values or `"residuals"`.

- ...:

  Currently unused.

## Value

A numeric matrix of fitted values or residuals, depending on `type`.

## See also

[`residuals.csdm_fit()`](residuals.csdm_fit.md),
[`summary.csdm_fit()`](summary.csdm_fit.md)
