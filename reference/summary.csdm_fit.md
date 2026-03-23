# Summarize csdm model estimation results

Computes post-estimation summaries for `csdm_fit` objects, including
mean-group coefficient inference, model-level diagnostics, and
model-specific summary tables (for example, short-run and long-run
blocks for CS-ARDL).

## Usage

``` r
# S3 method for class 'csdm_fit'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  A fitted model object of class `csdm_fit`.

- digits:

  Number of digits to print.

- ...:

  Further arguments passed to methods.

## Value

An object of class `summary.csdm_fit` with core metadata
(call/formula/model/N/T), coefficient tables, fit statistics, and
model-specific components for printing and downstream inspection.

## Details

### Reported inference

For each coefficient \\\hat\beta_k\\, the summary reports standard
errors, \\z\\-statistics, and two-sided normal-approximation p-values:
\$\$z_k = \frac{\hat\beta_k}{\operatorname{se}(\hat\beta_k)}, \qquad p_k
= 2\\1-\Phi(\|z_k\|)\\.\$\$

### Diagnostics

The printed summary shows the classic Pesaran CD diagnostic by default.
Extended diagnostics (CDw, CDw+, CD\*) are available through
[`cd_test()`](cd_test.md).

## See also

[`print.summary.csdm_fit()`](print.summary.csdm_fit.md),
[`cd_test()`](cd_test.md), [`coef.csdm_fit()`](coef.csdm_fit.md),
[`vcov.csdm_fit()`](vcov.csdm_fit.md)

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
s <- summary(fit)
s
#> csdm summary: Static Common Correlated Error Model (CCE)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 10, T: 38
#> Number of obs: 380
#> R-squared (mg): 0.9643
#> CD = -3.5806, p = 3e-04
#> (For additional CD diagnostics, use cd_test())
#> 
#> Mean Group:
#>               Coef. Std. Err.       z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept)  0.5424    2.7204  0.1994 0.8420         -4.7895   5.8743
#> log_hc      -0.8807    1.1671 -0.7546 0.4505         -3.1682   1.4069
#> log_ck       0.1597    0.1263  1.2642 0.2061         -0.0879   0.4072
#> log_ngd      0.7779    0.5174  1.5034 0.1327         -0.2363   1.7920
#> 
#> Mean Group Variables: log_hc, log_ck, log_ngd
#> Cross Sectional Averaged Variables: log_rgdpo, log_hc, log_ck, log_ngd (lags=0)
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
