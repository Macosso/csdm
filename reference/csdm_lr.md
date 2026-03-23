# Specification: Long-run configuration

Specification: Long-run configuration

## Usage

``` r
csdm_lr(
  vars = NULL,
  type = c("none", "ecm", "ardl", "csdl"),
  ylags = 0,
  xdlags = 0,
  options = list()
)
```

## Arguments

- vars:

  Reserved for future use.

- type:

  One of c("none","ecm","ardl","csdl").

- ylags:

  Integer \>= 0. Within-unit lags of the dependent variable to include
  when supported by the chosen model/type.

- xdlags:

  Integer \>= 0. Scalar distributed lags to apply to each RHS regressor
  when supported by the chosen model/type.

- options:

  Reserved for future use.

## Value

A spec object (list) used by csdm().

## Examples

``` r
# Long-run / dynamic configuration (ARDL-style lags)
lr <- csdm_lr(type = "ardl", ylags = 1)
lr
#> $vars
#> NULL
#> 
#> $type
#> [1] "ardl"
#> 
#> $ylags
#> [1] 1
#> 
#> $xdlags
#> [1] 0
#> 
#> $options
#> list()
#> 
#> attr(,"class")
#> [1] "csdm_lr_spec"

# Minimal end-to-end DCCE example (kept small for speed)
data(PWT_60_07, package = "csdm")
df <- PWT_60_07
keep_ids <- unique(df$id)[1:10]
df_small <- df[df$id %in% keep_ids & df$year >= 1970, ]
fit <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df_small,
  id = "id",
  time = "year",
  model = "dcce",
  csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), lags = 3),
  lr = csdm_lr(type = "ardl", ylags = 1)
)
#> cd_test: Dropped 3 incomplete time periods (7.9%). Balanced panel: 10 units x 35 periods.
#> cd_test: Dropped 3 incomplete time periods (7.9%). Balanced panel: 10 units x 35 periods.
summary(fit)
#> csdm summary: Dynamic Common Correlated Error Model (DCCE)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 10, T: 38
#> Number of obs: 350
#> R-squared (mg): 0.9844
#> CD = -3.368, p = 8e-04
#> (For additional CD diagnostics, use cd_test())
#> 
#> Mean Group:
#>                 Coef. Std. Err.      z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept)    4.1315    5.3606 0.7707 0.4409         -6.3751  14.6381
#> log_hc         0.2549    1.0228 0.2492 0.8032         -1.7498   2.2596
#> log_ck         0.4697    0.2359 1.9914 0.0464       *  0.0074   0.9320
#> log_ngd        0.3952    1.4560 0.2714 0.7861         -2.4585   3.2490
#> lag1_log_rgdpo 0.0736    0.0697 1.0559 0.2910         -0.0630   0.2101
#> 
#> Mean Group Variables: log_hc, log_ck, log_ngd, lag1_log_rgdpo
#> Cross Sectional Averaged Variables: log_rgdpo, log_hc, log_ck, log_ngd (lags=3)
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
