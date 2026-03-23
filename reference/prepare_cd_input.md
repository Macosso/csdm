# Prepare residual matrices for CD and CD\* diagnostics

Cleans and transforms an \\N x T\\ residual matrix for cross-sectional
dependence testing. Operations include:

1.  Dropping time periods with fewer than `min_per_time` finite
    observations.

2.  Optional row-wise standardization to unit variance over available
    times.

3.  Optional demeaning across units at each time (recommended for CD).

## Usage

``` r
prepare_cd_input(
  E,
  standardize = c("row", "none"),
  demean_time = TRUE,
  min_per_time = 2L
)
```

## Arguments

- E:

  A numeric matrix of residuals (\\N x T\\); rows are units, columns are
  time; may be unbalanced (contain `NA`).

- standardize:

  One of `"row"`, `"none"`. If `"row"`, scale each row by its observed
  standard deviation.

- demean_time:

  Logical; if `TRUE`, subtract the cross-sectional mean at each time
  from available residuals in that column.

- min_per_time:

  Integer; drop time columns with fewer than this many finite
  observations.

## Value

A list with:

- Z:

  Processed residual matrix (\\N x T^\*\\) after
  filtering/standardizing/demeaning.

- kept_t:

  Integer indices of kept time columns (relative to the original `E`).

- m_t:

  Integer vector of cross-sectional counts per kept time (number of
  finite rows).

- row_sds:

  Numeric vector of row standard deviations used (invisibly `NA` if
  `standardize="none"`).

- col_means:

  Numeric vector of time means subtracted when `demean_time=TRUE`.

## Details

### Transformation steps

1.  Time periods with fewer than `min_per_time` finite observations are
    removed.

2.  If `standardize = "row"`, each unit is scaled by its observed
    standard deviation.

3.  If `demean_time = TRUE`, each time slice is demeaned across
    available units.

### Why this preprocessing matters

CD-type tests are sensitive to scale heterogeneity and sparse columns in
unbalanced panels. This helper creates a better-conditioned input matrix
while preserving as much usable information as possible.

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
E <- get_residuals(fit)
prep <- prepare_cd_input(E, standardize = "row", demean_time = TRUE, min_per_time = 3)
dim(prep$Z)
#> [1] 10 38
```
