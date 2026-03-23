# Panel Model Estimation with Cross-Sectional Dependence

Estimate heterogeneous panel data models with optional cross-sectional
augmentation and dynamic structure. The interface supports Mean Group
(MG), Common Correlated Effects (CCE), Dynamic CCE (DCCE), and
Cross-Sectionally Augmented ARDL (CS-ARDL) estimators with a consistent
specification workflow for cross-sectional averages, lag structure, and
variance-covariance estimation.

## Usage

``` r
csdm(
  formula,
  data,
  id,
  time,
  model = c("mg", "cce", "dcce", "cs_ardl", "cs_ecm", "cs_dl"),
  csa = csdm_csa(),
  lr = csdm_lr(),
  pooled = csdm_pooled(),
  trend = c("none", "unit", "pooled"),
  fullsample = FALSE,
  mgmissing = FALSE,
  vcov = csdm_vcov(),
  ...
)
```

## Arguments

- formula:

  Model formula of the form `y ~ x1 + x2`.

- data:

  A `data.frame` (or
  [`plm::pdata.frame`](https://rdrr.io/pkg/plm/man/pdata.frame.html))
  containing the variables in `formula`.

- id, time:

  Column names (strings) for the unit and time indexes. If `data` is a
  `pdata.frame`, these are taken from its index and the provided values
  are ignored.

- model:

  Estimator to fit. One of `"mg"`, `"cce"`, `"dcce"`, or `"cs_ardl"`.

- csa:

  Cross-sectional-average specification, created by
  [`csdm_csa()`](csdm_csa.md).

- lr:

  Long-run or dynamic specification, created by
  [`csdm_lr()`](csdm_lr.md).

- pooled:

  Pooled specification (reserved for future use), created by
  [`csdm_pooled()`](csdm_pooled.md).

- trend:

  One of `"none"` or `"unit"` (adds a linear unit trend). `"pooled"` is
  reserved and not implemented.

- fullsample:

  Logical; reserved for future extensions.

- mgmissing:

  Logical; reserved for future extensions.

- vcov:

  Variance-covariance specification, created by
  [`csdm_vcov()`](csdm_vcov.md).

- ...:

  Reserved for future extensions.

## Value

An object of class `csdm_fit` containing estimated coefficients,
residuals, variance-covariance estimates, model metadata, and
diagnostics. Use [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`cd_test()`](cd_test.md) to access standard outputs.

## Details

Let \\i = 1, \ldots, N\\ index cross-sectional units and \\t = 1,
\ldots, T\\ index time. A baseline heterogeneous panel model is

\$\$y\_{it} = \alpha_i + \beta_i^T x\_{it} + u\_{it}.\$\$

Here \\\alpha_i\\ is a unit-specific intercept, \\x\_{it}\\ is a vector
of regressors, \\\beta_i\\ is a vector of unit-specific slopes, and
\\u\_{it}\\ is an error term that may exhibit cross-sectional
dependence.

Cross-sectional averages are specified through
[`csdm_csa()`](csdm_csa.md) and dynamic or long-run structure is
specified through [`csdm_lr()`](csdm_lr.md). This keeps the model
interface consistent across estimators while allowing the degree of
cross-sectional augmentation and lag structure to vary by application.

**Implemented estimators**

**MG (Pesaran and Smith, 1995)**

The Mean Group estimator fits separate regressions for each unit and
averages the resulting coefficients:

\$\$\hat{\beta}\_{MG} = \frac{1}{N}\sum\_{i=1}^N \hat{\beta}\_i.\$\$

This estimator accommodates slope heterogeneity but does not explicitly
model cross-sectional dependence.

**CCE (Pesaran, 2006)**

Regressions are augmented with cross-sectional averages to proxy
unobserved common factors:

\$\$y\_{it} = \alpha_i + \beta_i^T x\_{it} + \gamma_i^T \bar{z}\_{t} +
v\_{it}.\$\$

A common choice is

\$\$\bar{z}\_t = (\bar{y}\_t, \bar{x}\_t),\$\$

with

\$\$\bar{x}\_t = \frac{1}{N}\sum\_{i=1}^N x\_{it}, \qquad \bar{y}\_t =
\frac{1}{N}\sum\_{i=1}^N y\_{it}.\$\$

More generally, \\\bar{z}\_t\\ collects the cross-sectional averages
specified in `csa`.

**DCCE (Chudik and Pesaran, 2015)**

Dynamic CCE extends CCE by allowing lagged dependent variables and
lagged cross-sectional averages:

\$\$y\_{it} = \alpha_i + \sum\_{p=1}^{P} \phi\_{ip} y\_{i,t-p} +
\beta_i^T x\_{it} + \sum\_{q=0}^{Q} \delta\_{iq}^T \bar{z}\_{t-q} +
e\_{it}.\$\$

In the package implementation, lagged dependent variables and
distributed lags of regressors are controlled through `lr`, while
contemporaneous and lagged cross-sectional averages are controlled
through `csa`.

**CS-ARDL (Chudik and Pesaran, 2015)**

In the package implementation, `model = "cs_ardl"` is obtained by first
estimating a cross-sectionally augmented ARDL-style regression in
levels, using the same dynamic specification as `model = "dcce"`, and
then transforming the unit-specific coefficients into adjustment and
long-run parameters.

The underlying unit-level regression is of the form

\$\$y\_{it} = \alpha_i + \sum\_{p=1}^{P} \phi\_{ip} y\_{i,t-p} +
\sum\_{q=0}^{Q} \beta\_{iq}^T x\_{i,t-q} + \sum\_{s=0}^{S}
\omega\_{is}^T \bar{z}\_{t-s} + e\_{it}.\$\$

From this dynamic specification, the package recovers the implied
error-correction form

\$\$\Delta y\_{it} = \alpha_i + \varphi_i \left(y\_{i,t-1} - \theta_i^T
x\_{i,t-1}\right) + \sum\_{j=1}^{P-1} \lambda\_{ij} \Delta y\_{i,t-j} +
\sum\_{j=0}^{Q-1} \psi\_{ij}^T \Delta x\_{i,t-j} + \sum\_{s=0}^{S}
\tilde{\omega}\_{is}^T \bar{z}\_{t-s} + e\_{it},\$\$

where \\\varphi_i\\ is the adjustment coefficient and \\\theta_i\\ is
the implied long-run relationship. In the current implementation, these
quantities are computed from the estimated lag polynomials rather than
from a direct ECM regression.

**Identification and assumptions**

MG requires sufficient time-series variation within each unit.

CCE relies on cross-sectional averages acting as proxies for latent
common factors, together with adequate cross-sectional and time
dimensions.

DCCE additionally requires enough time periods to support lagged
dependent variables, distributed lags, and lagged cross-sectional
averages.

CS-ARDL requires sufficient time length for the distributed-lag
structure and is intended for applications where both short-run dynamics
and long-run relationships are of interest in the presence of common
factors.

## References

Pesaran MH, Smith R (1995). “Estimating long-run relationships from
dynamic heterogeneous panels.” *Journal of Econometrics*, **68**(1),
79–113.

Pesaran MH (2006). “Estimation and inference in large heterogeneous
panels with multifactor error structure.” *Econometrica*, **74**(4),
967–1012.

Chudik A, Pesaran MH (2015). “Common correlated effects estimation of
heterogeneous dynamic panel data models with weakly exogenous
regressors.” *Journal of Econometrics*, **188**(2), 393–420.

## Examples

``` r
library(csdm)
data(PWT_60_07, package = "csdm")
df <- PWT_60_07

# Keep examples fast but fully runnable
keep_ids <- unique(df$id)[1:10]
df_small <- df[df$id %in% keep_ids & df$year >= 1970, ]

# Mean Group (MG)
mg <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df_small, id = "id", time = "year", model = "mg"
)
summary(mg)
#> csdm summary: Mean Group Model (MG)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 10, T: 38
#> Number of obs: 380
#> R-squared (mg): 0.9151
#> CD = 1.1106, p = 0.2667
#> (For additional CD diagnostics, use cd_test())
#> 
#> Mean Group:
#>              Coef. Std. Err.      z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept) 5.9905    1.1639 5.1470 0.0000     ***  3.7093   8.2717
#> log_hc      0.6417    1.1367 0.5645 0.5724         -1.5863   2.8697
#> log_ck      0.2707    0.1420 1.9066 0.0566       . -0.0076   0.5490
#> log_ngd     0.5797    0.4456 1.3008 0.1933         -0.2938   1.4531
#> 
#> Mean Group Variables: log_hc, log_ck, log_ngd
#> Cross Sectional Averaged Variables: none (lags=0)
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Common Correlated Effects (CCE)
cce <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df_small, id = "id", time = "year", model = "cce",
  csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"))
)
summary(cce)
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

# Dynamic CCE (DCCE)
dcce <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df_small, id = "id", time = "year", model = "dcce",
  csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), lags = 3),
  lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 0)
)
#> cd_test: Dropped 3 incomplete time periods (7.9%). Balanced panel: 10 units x 35 periods.
#> cd_test: Dropped 3 incomplete time periods (7.9%). Balanced panel: 10 units x 35 periods.
summary(dcce)
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

# CS-ARDL
cs_ardl <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df_small, id = "id", time = "year", model = "cs_ardl",
  csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), lags = 3),
  lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 0)
)
#> cd_test: Dropped 3 incomplete time periods (7.9%). Balanced panel: 10 units x 35 periods.
#> cd_test: Dropped 3 incomplete time periods (7.9%). Balanced panel: 10 units x 35 periods.
summary(cs_ardl)
#> csdm summary: Cross-Sectional ARDL (CS-ARDL)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 10, T: 38
#> Number of obs: 350
#> R-squared (mg): 0.9844
#> 
#> CD = -3.368, p = 8e-04
#> (For additional CD diagnostics, use cd_test())
#> 
#> Short Run Est.
#>                 Coef. Std. Err.      z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept)    4.1315    5.3606 0.7707 0.4409         -6.3751  14.6381
#> log_hc         0.2549    1.0228 0.2492 0.8032         -1.7498   2.2596
#> log_ck         0.4697    0.2359 1.9914 0.0464       *  0.0074   0.9320
#> log_ngd        0.3952    1.4560 0.2714 0.7861         -2.4585   3.2490
#> lag1_log_rgdpo 0.0736    0.0697 1.0559 0.2910         -0.0630   0.2101
#> 
#> Adjust. Term
#>                Coef. Std. Err.        z P>|z| Signif. CI 2.5% CI 97.5%
#> lr_log_rgdpo -0.9264    0.0697 -13.2994     0     ***  -1.063  -0.7899
#> 
#> Long Run Est.
#>             Coef. Std. Err.      z  P>|z| Signif. CI 2.5% CI 97.5% n_used
#> lr_log_hc  0.3172    0.9712 0.3266 0.7440         -1.5863   2.2207     10
#> lr_log_ck  0.4165    0.2469 1.6871 0.0916       . -0.0674   0.9004     10
#> lr_log_ngd 0.3572    1.6161 0.2210 0.8251         -2.8103   3.5247     10
#> 
#> Mean Group Variables: lag1_log_rgdpo, log_hc, log_ck, log_ngd
#> Cross Sectional Averaged Variables: log_rgdpo, log_hc, log_ck, log_ngd (lags=3)
#> Long Run Variables: log_hc, log_ck, log_ngd
#> Cointegration variable(s): log_rgdpo
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
