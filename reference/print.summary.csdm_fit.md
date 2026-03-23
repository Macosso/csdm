# Print method for csdm summary objects

Formats and prints a `summary.csdm_fit` object. Output adapts to model
type and includes coefficient tables, selected goodness-of-fit
diagnostics, and compact model metadata.

## Usage

``` r
# S3 method for class 'summary.csdm_fit'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `summary.csdm_fit` object.

- digits:

  Number of digits to print.

- ...:

  Further arguments passed to methods.

## Value

Invisibly returns `x`.

## Details

The printout includes classic Pesaran CD diagnostics from the summary
object. For a full CD diagnostic panel (CD, CDw, CDw+, CD\*), use
[`cd_test()`](cd_test.md) on the fitted model.

## See also

[`summary.csdm_fit()`](summary.csdm_fit.md), [`cd_test()`](cd_test.md)
