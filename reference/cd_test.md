# Cross-sectional dependence (CD) tests for panel residuals

Computes Pesaran CD, CDw, CDw+, and CD\* tests for cross-sectional
dependence in panel residuals. The implementation supports residual
matrices or fitted `csdm_fit` objects and provides consistent handling
of unbalanced panels.

## Usage

``` r
cd_test(object, ...)

# Default S3 method
cd_test(
  object,
  type = c("CD", "CDw", "CDw+", "CDstar", "all"),
  n_pc = 4L,
  seed = NULL,
  min_overlap = 2L,
  na.action = c("drop.incomplete.times", "pairwise"),
  ...
)

# S3 method for class 'csdm_fit'
cd_test(
  object,
  type = c("CD", "CDw", "CDw+", "CDstar", "all"),
  n_pc = 4L,
  seed = NULL,
  min_overlap = 2L,
  na.action = c("drop.incomplete.times", "pairwise"),
  ...
)

# S3 method for class 'cd_test'
print(x, digits = 3, ...)
```

## Arguments

- object:

  A `csdm_fit` model object or a numeric matrix of residuals (N x T).

- ...:

  Additional arguments passed to methods.

- type:

  Which test(s) to compute: one of `"CD"`, `"CDw"`, `"CDw+"`,
  `"CDstar"`, or `"all"` (default: `"CD"`).

- n_pc:

  Number of principal components for CD\* (default 4).

- seed:

  Integer seed for weight draws in CDw/CDw+ (default NULL = no seed
  set).

- min_overlap:

  Minimum number of overlapping time periods required for a unit pair to
  be included in CD/CDw/CDw+ (default 2).

- na.action:

  How to handle missing data: `"drop.incomplete.times"` (default)
  removes time periods with any missing observations to create a
  balanced panel for CD\*; `"pairwise"` uses pairwise correlations for
  CD/CDw/CDw+ and warns for CD\*.

- x:

  An object of class `cd_test`.

- digits:

  Number of digits to print (default 3).

## Value

An object of class `cd_test` with fields `tests`, `type`, `N`, `T`,
`na.action`, and `call`. The `tests` list contains one or more test
results, each with `statistic` and `p.value`.

## Details

### Notation

Let \\E\\ be the residual matrix with \\N\\ cross-sectional units and
\\T\\ time periods. For each unit pair \\(i,j)\\, let \\T\_{ij}\\ be the
number of overlapping time periods and \\\rho\_{ij}\\ the pairwise
correlation.

### Test statistics

- CD (Pesaran, 2015):

  \$\$CD = \sqrt{\frac{2}{N(N-1)}} \sum\_{i\<j} \sqrt{T\_{ij}} \\
  \rho\_{ij}\$\$

- CDw (Juodis and Reese, 2021):

  Random sign flips \\w_i \in \\-1,1\\\\ are applied to residuals before
  computing correlations. The statistic is CD applied to the
  sign-flipped data.

- CDw+ (Fan, Liao, and Yao, 2015):

  Power enhancement adds a sparse thresholding term to CDw. The
  threshold is \$\$c_N = \sqrt{\frac{2 \log(N)}{T}}\$\$ and the power
  term sums \\\sqrt{T\_{ij}} \|\rho\_{ij}\|\\ for pairs exceeding the
  threshold.

- CD\* (Pesaran and Xie, 2021):

  CD is computed on residuals after removing `n_pc` principal components
  from \\E\\. This provides a bias-corrected test under multifactor
  errors.

### Missing data and balance

- CD, CDw, CDw+:

  Always use pairwise-complete observations. Each pairwise correlation
  uses available overlaps.

- CD\*:

  Requires a balanced panel. By default,
  `na.action = "drop.incomplete.times"` removes any time period with
  missing observations. With `na.action = "pairwise"`, CD\* returns `NA`
  and a warning when missing values are present.

## References

Pesaran MH (2015). “Testing weak cross-sectional dependence in large
panels.” *Econometric Reviews*, **34**(6-10), 1089–1117.

Pesaran MH (2021). “General diagnostic tests for cross-sectional
dependence in panels.” *Empirical Economics*, **60**(1), 13–50.

Juodis A, Reese S (2021). “The incidental parameters problem in testing
for remaining cross-sectional correlation.” *Journal of Business and
Economic Statistics*, **40**(3), 1191–1203.

Fan J, Liao Y, Yao J (2015). “Power Enhancement in High-Dimensional
Cross-Section Tests.” *Econometrica*, **83**(4), 1497–1541.

Pesaran MH, Xie Y (2021). “A bias-corrected CD test for error
cross-sectional dependence in panel models.” *Econometric Reviews*,
**41**(6), 649–677.

## Examples

``` r
# Simulate independent and dependent panels
set.seed(1)
E_indep <- matrix(rnorm(100), nrow = 10)
E_dep <- matrix(rnorm(10), nrow = 10, ncol = 10, byrow = TRUE)

# Compute all tests
cd_test(E_indep, type = "all")
#> Cross-sectional dependence tests
#> N = 10, T = 10
#> 
#>        statistic p.value
#> CD        -1.325   0.185
#> CDw       -0.722   0.470
#> CDw+      32.982   0.000
#> CDstar     0.528   0.598
cd_test(E_dep, type = "all")
#> Cross-sectional dependence tests
#> N = 10, T = 10
#> 
#>        statistic p.value
#> CD        21.213   0.000
#> CDw       -2.357   0.018
#> CDw+     139.945   0.000
#> CDstar     0.369   0.712

# Specific test with parameters
cd_test(E_indep, type = "CDstar", n_pc = 2)
#> Cross-sectional dependence tests
#> N = 10, T = 10
#> 
#>        statistic p.value
#> CDstar    -0.219   0.826

# From a fitted csdm model
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
cd_test(fit, type = "all")
#> Cross-sectional dependence tests
#> N = 10, T = 38
#> 
#>        statistic p.value
#> CD        -3.581   0.000
#> CDw        1.170   0.242
#> CDw+      68.347   0.000
#> CDstar    -3.124   0.002
```
