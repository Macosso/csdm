# Extract residual matrices for panel diagnostics

Unified accessor that returns an \\N x T\\ residual matrix suitable for
cross-sectional dependence diagnostics and post-estimation analysis.

## Usage

``` r
get_residuals(object, type = c("auto", "cce", "pca", "pca_std"), strict = TRUE)
```

## Arguments

- object:

  A fitted model object supported by this package (e.g., class
  `csdm_fit`), or directly a numeric matrix of residuals shaped as \\N x
  T\\.

- type:

  Character string selecting which residuals to return when available:
  one of `"auto"`, `"cce"`, `"pca"`, or `"pca_std"`.

  - `"auto"`: prefer standardized PCA residuals if present, otherwise
    PCA residuals, otherwise CCE residuals, otherwise `object` if it is
    a matrix.

  - `"cce"`: residuals from the CCE-augmented per-unit regressions.

  - `"pca"`: residuals after removing estimated factors from \\\hat
    v\_{it}\\.

  - `"pca_std"`: `"pca"` residuals standardized by unit-specific scale
    (recommended for CD).

- strict:

  Logical; if `TRUE`, error on unsupported objects. If `FALSE`, return
  `NULL` when residuals cannot be found.

## Value

A numeric matrix of residuals with rows = units and columns = time,
preserving `rownames` and `colnames` when available; or `NULL` if
nothing suitable is found and `strict = FALSE`.

## Details

### Residual types

- cce:

  Residuals from the cross-sectionally augmented unit regressions.

- pca:

  Residuals after principal-component factor removal.

- pca_std:

  PCA residuals standardized by unit-specific scale.

- auto:

  Priority rule: `pca_std` -\> `pca` -\> `cce` -\> generic residual
  slots.

### Assumptions and usage

The returned matrix is intended for diagnostics that operate on
unit-time panels, including [`cd_test()`](cd_test.md). Missing values
are preserved unless downstream routines explicitly filter or balance
the panel.

## Examples

``` r
data(PWT_60_07, package = "csdm")
df <- PWT_60_07
ids <- unique(df$id)[1:10]
df_small <- df[df$id %in% ids & df$year >= 1970, ]

fit <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df_small,
  id = "id",
  time = "year",
  model = "cce",
  csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"))
)

E <- get_residuals(fit, type = "auto")
dim(E)
#> [1] 10 38
```
