# Introduction to csdm: Panel Data Models with Cross-Sectional Dependence

![](../reference/figures/logo.svg)

------------------------------------------------------------------------

## Overview

The `csdm` package implements econometric methods for panel data with
cross-sectional dependence (CSD). In many applications, observations
across units (e.g., countries, firms, regions) are not
independent—macroeconomic shocks, trade relationships, or spillovers
create correlation across cross-sectional units. The `csdm` package
provides robust estimators that account for this dependence structure,
plus diagnostic tests to detect and characterize it.

This vignette demonstrates four core estimation methods and related
inference tools on real panel data from the Penn World Table (PWT).

## Methodology: Four Estimators

### Model Specification

The [`csdm()`](../reference/csdm.md) interface estimates heterogeneous
panel data models with optional cross-sectional augmentation and dynamic
structure. Let $i = 1,\ldots,N$ index cross-sectional units and
$t = 1,\ldots,T$ index time. A baseline heterogeneous panel model is

$$y_{it} = \alpha_{i} + \beta_{i}\prime x_{it} + u_{it},\qquad i = 1,\ldots,N,\; t = 1,\ldots,T$$

where:

- $y_{it}$ is the outcome variable for unit $i$ at time $t$
- $\alpha_{i}$ is a unit-specific intercept
- $\beta_{i}$ is a $(k \times 1)$ vector of unit-specific slopes
- $x_{it}$ is a $(k \times 1)$ vector of explanatory variables
- $u_{it}$ is the error term, which may exhibit cross-sectional
  dependence

The inner product $\beta_{i}\prime x_{it}$ is scalar-valued.
Heterogeneous slopes allow each unit to respond differently to the
regressors. In many applications, cross-sectional dependence arises
because the error term contains unobserved common factors. The
estimators implemented in [`csdm()`](../reference/csdm.md) differ in how
they handle this dependence and whether they allow for dynamic
adjustment.

------------------------------------------------------------------------

### 1. Mean Group (MG) Estimator

The Mean Group estimator fits separate regressions for each unit and
averages the resulting coefficients:

$${\widehat{\beta}}_{MG} = \frac{1}{N}\sum\limits_{i = 1}^{N}{\widehat{\beta}}_{i}$$

**Key idea**: Estimation is performed unit by unit, with no pooling of
slope coefficients across cross-sectional units.

**Interpretation**:

- ${\widehat{\beta}}_{MG}$ is the cross-sectional average of the
  unit-specific estimates
- all slope coefficients are allowed to differ across units

**Properties**:

- accommodates slope heterogeneity
- requires sufficient time-series information within each unit
- does not explicitly model cross-sectional dependence

**Use case**: A natural benchmark when the main concern is heterogeneous
slopes and no explicit factor structure is imposed.

------------------------------------------------------------------------

### 2. Common Correlated Effects (CCE) Estimator

The CCE estimator augments each unit regression with cross-sectional
averages to proxy unobserved common factors:

$$y_{it} = \alpha_{i} + \beta_{i}\prime x_{it} + \gamma_{i}\prime{\bar{z}}_{t} + v_{it}$$

where ${\bar{z}}_{t}$ collects the cross-sectional averages specified
through [`csdm_csa()`](../reference/csdm_csa.md), for example

$${\bar{z}}_{t} = \left( {\bar{y}}_{t},{\bar{x}}_{t} \right),\qquad{\bar{x}}_{t} = \frac{1}{N}\sum\limits_{i = 1}^{N}x_{it},\qquad{\bar{y}}_{t} = \frac{1}{N}\sum\limits_{i = 1}^{N}y_{it}.$$

**Key idea**: Cross-sectional averages serve as proxies for latent
common factors that induce dependence across units.

**Interpretation**:

- $\beta_{i}$ measures the unit-specific effect conditional on the
  included cross-sectional averages
- $\gamma_{i}$ captures unit-specific exposure to the common components
  with ${\bar{z}}_{t}$ as a proxy.

**Properties**:

- allows heterogeneous slopes
- augments the regression with cross-sectional averages supplied through
  `csa`
- suitable when cross-sectional dependence is driven by latent common
  shocks

**Use case**: When dependence across units is believed to reflect common
unobserved factors.

------------------------------------------------------------------------

### 3. Dynamic CCE (DCCE) Estimator

The DCCE estimator extends CCE to dynamic settings by including lagged
dependent variables, optional distributed lags of regressors, and lagged
cross-sectional averages:

$$y_{it} = \alpha_{i} + \sum\limits_{p = 1}^{P}\phi_{ip}y_{i,t - p} + \sum\limits_{q = 0}^{Q}\beta_{iq}\prime x_{i,t - q} + \sum\limits_{s = 0}^{S}\delta_{is}\prime{\bar{z}}_{t - s} + e_{it}$$

where the dynamic structure is controlled through
[`csdm_lr()`](../reference/csdm_lr.md) and the cross-sectional averages
and their lags are controlled through
[`csdm_csa()`](../reference/csdm_csa.md).

**Key idea**: Dynamics are introduced directly in the unit equation,
while lagged cross-sectional averages help absorb common factor
dependence over time.

**Interpretation**:

- $\phi_{ip}$ captures unit-specific persistence
- $\beta_{iq}$ captures contemporaneous and lagged effects of regressors
- $\delta_{is}$ captures the effect of contemporaneous and lagged common
  components

**Properties**:

- allows heterogeneous dynamic adjustment across units
- combines lagged dependent variables, optional distributed lags, and
  cross-sectional augmentation
- requires enough time periods to support the chosen lag structure

**Use case**: When the outcome is persistent over time and
cross-sectional dependence remains important.

------------------------------------------------------------------------

### 4. Cross-Sectionally Augmented ARDL (CS-ARDL)

In the current [`csdm()`](../reference/csdm.md) implementation,
`model = "cs_ardl"` is obtained by first estimating a cross-sectionally
augmented ARDL-style regression in levels, using the same dynamic
specification as `model = "dcce"`, and then transforming the estimated
unit-specific coefficients into adjustment and long-run parameters.

The underlying unit-level regression is

$$y_{it} = \alpha_{i} + \sum\limits_{p = 1}^{P}\phi_{ip}y_{i,t - p} + \sum\limits_{q = 0}^{Q}\beta_{iq}\prime x_{i,t - q} + \sum\limits_{s = 0}^{S}\omega_{is}\prime{\bar{z}}_{t - s} + e_{it}$$

From this dynamic specification, the implied error-correction form is

$$\Delta y_{it} = \alpha_{i} + \varphi_{i}\left( y_{i,t - 1} - \theta_{i}\prime x_{i,t - 1} \right) + \sum\limits_{j = 1}^{P - 1}\lambda_{ij}\Delta y_{i,t - j} + \sum\limits_{j = 0}^{Q - 1}\psi_{ij}\prime\Delta x_{i,t - j} + \sum\limits_{s = 0}^{S}{\widetilde{\omega}}_{is}\prime{\bar{z}}_{t - s} + e_{it}$$

where the dynamic structure is controlled through
[`csdm_lr()`](../reference/csdm_lr.md) and the cross-sectional averages
are supplied through [`csdm_csa()`](../reference/csdm_csa.md).

**Key idea**: `cs_ardl` reports the implied short-run and long-run
quantities from a cross-sectionally augmented ARDL fit.

**Interpretation**:

- $\theta_{i}$ is the unit-specific long-run relationship
- $\varphi_{i}$ is the implied speed of adjustment back toward
  equilibrium
- $\psi_{ij}$ captures short-run effects of changes in regressors
- ${\widetilde{\omega}}_{is}$ captures the role of common
  cross-sectional components

**Properties**:

- supports heterogeneous short-run and long-run dynamics
- combines ARDL-style dynamics with cross-sectional augmentation
- recovers adjustment and long-run coefficients from estimated lag
  polynomials rather than fitting a separate ECM directly

**Use case**: When the objective is to study long-run relationships
together with heterogeneous short-run adjustment in panels affected by
common factors.

------------------------------------------------------------------------

### Cross-Sectional Averages and Dynamic Structure

Two helper specifications control the main extensions in
[`csdm()`](../reference/csdm.md):

- [`csdm_csa()`](../reference/csdm_csa.md) defines which variables enter
  as cross-sectional averages and how many lags of those averages are
  included
- [`csdm_lr()`](../reference/csdm_lr.md) defines the dynamic or long-run
  structure, such as lagged dependent variables and distributed lags

This design keeps the estimation interface consistent across the four
estimators while allowing the model specification to vary by
application.

------------------------------------------------------------------------

### Summary

| Estimator | Heterogeneous Slopes | Cross-Sectional Averages | Dynamics | Long-Run Structure |
|-----------|----------------------|--------------------------|----------|--------------------|
| MG        | Yes                  | No                       | No       | No                 |
| CCE       | Yes                  | Yes                      | No       | No                 |
| DCCE      | Yes                  | Yes                      | Yes      | No                 |
| CS-ARDL   | Yes                  | Yes                      | Yes      | Yes                |

## Data: Penn World Table Subset

The `PWT_60_07` dataset contains macroeconomic indicators for 93
countries covering 1960–2007 (48 years). Key variables include:

- `id`: Country identifier
- `year`: Calendar year (1960–2007)
- `log_rgdpo`: Log real GDP per capita
- `log_hc`: Log human capital index
- `log_ck`: Log capital stock
- `log_ngd`: Log government debt (control variable)

``` r
data(PWT_60_07, package = "csdm")
head(PWT_60_07, 10)
#>    id year log_rgdpo    log_hc   log_ck   log_ngd
#> 1   1 1960  7.780284 0.7042058 11.33559        NA
#> 2   1 1961  7.792448 0.7096307 11.39625 -2.714639
#> 3   1 1962  7.800655 0.7150558 11.45449 -2.719555
#> 4   1 1963  7.751311 0.7204807 11.44691 -2.723755
#> 5   1 1964  7.786854 0.7259058 11.48774 -2.727269
#> 6   1 1965  7.879184 0.7313308 11.54133 -2.730137
#> 7   1 1966  7.885357 0.7386493 11.57771 -2.737341
#> 8   1 1967  7.903803 0.7459679 11.61129 -2.744694
#> 9   1 1968  7.923327 0.7532863 11.64887 -2.745329
#> 10  1 1969  8.004197 0.7606049 11.71465 -2.739602
str(PWT_60_07)
#> Classes 'tbl_df', 'tbl' and 'data.frame':    4464 obs. of  6 variables:
#>  $ id       : num  1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "format.stata")= chr "%12.0g"
#>  $ year     : num  1960 1961 1962 1963 1964 ...
#>   ..- attr(*, "label")= chr "(firstnm) year"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ log_rgdpo: num  7.78 7.79 7.8 7.75 7.79 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ log_hc   : num  0.704 0.71 0.715 0.72 0.726 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ log_ck   : num  11.3 11.4 11.5 11.4 11.5 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ log_ngd  : num  NA -2.71 -2.72 -2.72 -2.73 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"

# For computational speed in this vignette, use a subset: 
# first 15 countries, 1970-2007
first_15_ids <- unique(PWT_60_07$id)[1:15]
df <- subset(PWT_60_07, id %in% first_15_ids & year >= 1970 & year <= 2007)
```

The panel is relatively balanced. We will use growth regressions:
modeling log-GDP (`log_rgdpo`) as a function of human capital
(`log_hc`), capital stock (`log_ck`), and government debt (`log_ngd`),
and test cross-sectional dependence in residuals.

## Package installation

To install the `csdm` package from CRAN, run:

``` r
install.packages("csdm")
```

To install the latest development version from GitHub, run:

``` r
install.packages("remotes")
remotes::install_github("Macosso/csdm")
```

## Model Estimation: Four Examples

All models are fitted with [`csdm()`](../reference/csdm.md), which
automatically detects the input structure and applies the appropriate
methodology. The key arguments are `id` and `time` to specify the
cross-sectional and time-period identifiers, and `model` to choose the
estimator. For CCE and DCCE, additional arguments (`csa` and `lr`)
specify treatment of cross-sectional averages and dynamics.

### Example 1: Mean Group (MG) Estimation

``` r
# MG: Separate regression per country, then average coefficients
fit_mg <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df,
  id = "id", 
  time = "year",
  model = "mg"
)

print(fit_mg)
#> csdm fit (mg)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#>             Estimate Std.Error
#> (Intercept)   6.5609    0.8648
#> log_hc        0.1725    0.8530
#> log_ck        0.3056    0.1359
#> log_ngd       0.7800    0.3777
summary(fit_mg)
#> csdm summary: Mean Group Model (MG)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#> Number of obs: 570
#> R-squared (mg): 0.9449
#> CD = 3.2379, p = 0.0012
#> (For additional CD diagnostics, use cd_test())
#> 
#> Mean Group:
#>              Coef. Std. Err.      z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept) 6.5609    0.8648 7.5867 0.0000     ***  4.8659   8.2558
#> log_hc      0.1725    0.8530 0.2022 0.8398         -1.4993   1.8443
#> log_ck      0.3056    0.1359 2.2485 0.0245       *  0.0392   0.5720
#> log_ngd     0.7800    0.3777 2.0652 0.0389       *  0.0398   1.5203
#> 
#> Mean Group Variables: log_hc, log_ck, log_ngd
#> Cross Sectional Averaged Variables: none (lags=0)
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

**Interpretation**: The MG estimate suggests that on average across
countries, increases in human capital, capital stock, and changes in
debt are associated with changes in real GDP. The standard errors
reflect cross-country heterogeneity in these relationships.

### Example 2: Common Correlated Effects (CCE)

``` r
# CCE: Add cross-sectional means to control for common shocks
fit_cce <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df,
  id = "id", 
  time = "year",
  model = "cce",
  csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"))
)

print(fit_cce)
#> csdm fit (cce)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#>             Estimate Std.Error
#> (Intercept)   1.9003    2.1195
#> log_hc       -1.4921    1.0152
#> log_ck        0.1367    0.0956
#> log_ngd       0.8075    0.2972
summary(fit_cce)
#> csdm summary: Static Common Correlated Error Model (CCE)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#> Number of obs: 570
#> R-squared (mg): 0.9777
#> CD = -2.6758, p = 0.0075
#> (For additional CD diagnostics, use cd_test())
#> 
#> Mean Group:
#>               Coef. Std. Err.       z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept)  1.9003    2.1195  0.8965 0.3700         -2.2540   6.0545
#> log_hc      -1.4921    1.0152 -1.4697 0.1416         -3.4819   0.4977
#> log_ck       0.1367    0.0956  1.4298 0.1528         -0.0507   0.3240
#> log_ngd      0.8075    0.2972  2.7175 0.0066      **  0.2251   1.3899
#> 
#> Mean Group Variables: log_hc, log_ck, log_ngd
#> Cross Sectional Averaged Variables: log_rgdpo, log_hc, log_ck, log_ngd (lags=0)
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

**Interpretation**: After accounting for global shocks (via
cross-sectional averages of all variables), the CCE coefficients and
standard errors may differ from MG. This indicates whether common
factors (e.g., technology, energy prices) are a major source of
cross-sectional dependence.

### Example 3: Dynamic CCE (DCCE)

``` r
# DCCE: Include dynamics and cross-sectional means
# Use lagged dependent variable to capture dynamic adjustment
fit_dcce <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df,
  id = "id", 
  time = "year",
  model = "dcce",
  csa = csdm_csa(
    vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), 
    lags = 3
  ),
  lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 0)
)

print(fit_dcce)
#> csdm fit (dcce)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#>                Estimate Std.Error
#> (Intercept)      9.2888    8.0386
#> log_hc          -1.9558    1.5659
#> log_ck           0.6666    0.2424
#> log_ngd         -0.3178    1.1271
#> lag1_log_rgdpo  -0.0745    0.0555
summary(fit_dcce)
#> csdm summary: Dynamic Common Correlated Error Model (DCCE)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#> Number of obs: 525
#> R-squared (mg): 0.9861
#> CD = -2.392, p = 0.0168
#> (For additional CD diagnostics, use cd_test())
#> 
#> Mean Group:
#>                  Coef. Std. Err.       z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept)     9.2888    8.0386  1.1555 0.2479         -6.4666  25.0441
#> log_hc         -1.9558    1.5659 -1.2490 0.2117         -5.0249   1.1132
#> log_ck          0.6666    0.2424  2.7499 0.0060      **  0.1915   1.1417
#> log_ngd        -0.3178    1.1271 -0.2819 0.7780         -2.5268   1.8913
#> lag1_log_rgdpo -0.0745    0.0555 -1.3424 0.1795         -0.1834   0.0343
#> 
#> Mean Group Variables: log_hc, log_ck, log_ngd, lag1_log_rgdpo
#> Cross Sectional Averaged Variables: log_rgdpo, log_hc, log_ck, log_ngd (lags=3)
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

**Interpretation**: The DCCE model includes lagged GDP to capture
dynamic adjustment. The lagged coefficient typically lies between
0.8–0.95, indicating strong income persistence. The coefficients on
other variables represent short-run elasticities; to compute long-run
effects, divide by $\left( 1 - \text{lag coefficient} \right)$.

### Example 4: Cross-Sectionally Augmented ARDL (CS-ARDL)

``` r
# CS-ARDL: Separate short-run and long-run dynamics
# Includes lagged dependent and lagged regressors
fit_csardl <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df,
  id = "id", 
  time = "year",
  model = "cs_ardl",
  csa = csdm_csa(
    vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"), 
    lags = 3
  ),
  lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 0)
)

print(fit_csardl)
#> csdm fit (cs_ardl)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#>                Estimate Std.Error
#> (Intercept)      9.2888    8.0386
#> log_hc          -1.9558    1.5659
#> log_ck           0.6666    0.2424
#> log_ngd         -0.3178    1.1271
#> lag1_log_rgdpo  -0.0745    0.0555
summary(fit_csardl)
#> csdm summary: Cross-Sectional ARDL (CS-ARDL)
#> Formula: log_rgdpo ~ log_hc + log_ck + log_ngd
#> N: 15, T: 38
#> Number of obs: 525
#> R-squared (mg): 0.9861
#> 
#> CD = -2.392, p = 0.0168
#> (For additional CD diagnostics, use cd_test())
#> 
#> Short Run Est.
#>                  Coef. Std. Err.       z  P>|z| Signif. CI 2.5% CI 97.5%
#> (Intercept)     9.2888    8.0386  1.1555 0.2479         -6.4666  25.0441
#> log_hc         -1.9558    1.5659 -1.2490 0.2117         -5.0249   1.1132
#> log_ck          0.6666    0.2424  2.7499 0.0060      **  0.1915   1.1417
#> log_ngd        -0.3178    1.1271 -0.2819 0.7780         -2.5268   1.8913
#> lag1_log_rgdpo -0.0745    0.0555 -1.3424 0.1795         -0.1834   0.0343
#> 
#> Adjust. Term
#>                Coef. Std. Err.        z P>|z| Signif. CI 2.5% CI 97.5%
#> lr_log_rgdpo -1.0745    0.0555 -19.3513     0     *** -1.1834  -0.9657
#> 
#> Long Run Est.
#>              Coef. Std. Err.       z  P>|z| Signif. CI 2.5% CI 97.5% n_used
#> lr_log_hc  -1.8378    1.4743 -1.2465 0.2126         -4.7274   1.0518     15
#> lr_log_ck   0.6106    0.2076  2.9409 0.0033      **  0.2037   1.0175     15
#> lr_log_ngd -0.4098    1.2364 -0.3314 0.7403         -2.8330   2.0134     15
#> 
#> Mean Group Variables: lag1_log_rgdpo, log_hc, log_ck, log_ngd
#> Cross Sectional Averaged Variables: log_rgdpo, log_hc, log_ck, log_ngd (lags=3)
#> Long Run Variables: log_hc, log_ck, log_ngd
#> Cointegration variable(s): log_rgdpo
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

**Interpretation**: The CS-ARDL model returns both short-run
coefficients (immediate response to shocks) and long-run coefficients
(equilibrium effect after full adjustment). The long-run elasticities
are often larger than short-run responses, consistent with gradual
accumulation effects in capital and human capital.

## Cross-Sectional Dependence Testing

After fitting a model, we can test whether residuals exhibit
cross-sectional dependence using the Pesaran CD test and related
variants. CSD tests detect whether residuals $u_{it}$ are correlated
across units—a key assumption violation that can bias standard errors.

### Four CD Test Types

All CD tests have null hypothesis: **residuals are cross-sectionally
independent**.

#### 1. Pesaran CD Test

The Pesaran CD statistic is:

$$CD = \sqrt{\frac{2}{N(N - 1)}}\sum\limits_{i = 1}^{N - 1}\sum\limits_{j = i + 1}^{N}{\widehat{\rho}}_{ij}\sqrt{T}$$

where ${\widehat{\rho}}_{ij}$ is the cross-sectional correlation between
residuals of units $i$ and $j$. The test statistic is approximately
standard normal under the null.

**Interpretation**: Large $|CD|$ rejects independence; both positive and
negative correlations are flagged. This is the most general CD test and
works even when $N$ is fixed and $\left. T\rightarrow\infty \right.$.

#### 2. Pesaran CD Weighted (CDw)

The CDw statistic applies unit-level random sign weights to the
cross-sectional correlations:

$$CD_{w} = \sqrt{\frac{2}{N(N - 1)}}\sum\limits_{i = 1}^{N - 1}\sum\limits_{j = i + 1}^{N}w_{i}w_{j}{\widehat{\rho}}_{ij}\sqrt{T}$$

where weights $w_{i} \in \{ - 1,1\}$ are independent random sign flips
assigned at the unit level and held fixed within a replication. This
random-weighting scheme improves the behavior of the test in the
presence of heteroskedasticity.

#### 3. Pesaran CD Weighted Plus (CDw+)

CDw+ uses the same unit-level random sign weights but applies them to a
bias-adjusted version of the CD statistic:

$$CD_{w}^{+} = \sqrt{\frac{2}{N(N - 1)}}\sum\limits_{i = 1}^{N - 1}\sum\limits_{j = i + 1}^{N}w_{i}w_{j}{\widehat{\rho}}_{ij}^{+}\sqrt{T}$$

where ${\widehat{\rho}}_{ij}^{+}$ denotes the adjusted cross-sectional
correlation. CDw+ is designed to improve robustness in large panels with
heteroskedasticity.

#### 4. Pesaran CD\*, Fan-Liao-Yao (FLY)

The CD\* statistic is a semiparametric refinement for large $N$ and $T$:

$$CD^{*} = \frac{1}{\sqrt{N(N - 1)}}\sum\limits_{i = 1}^{N - 1}\sum\limits_{j = i + 1}^{N}\left( {\widehat{\rho}}_{ij}^{2} - \tau_{T} \right)$$

where $\tau_{T}$ is a variance adjustment. FLY-type tests are designed
for large panel dimensions and provide robustness against certain forms
of weak cross-sectional dependence.

### Running CD Tests with Seed Selection

The [`cd_test()`](../reference/cd_test.md) function accepts the fitted
model and computes all test variants. Tests use a **random seed** to
initialize pseudo-random computations (for `cdw` and `cdw+`); setting a
`seed` ensures reproducibility of numerical results across runs.

``` r
# Test MG residuals for CSD
cd_mg <- cd_test(fit_mg, type = "CD")
print(cd_mg)
#> Cross-sectional dependence tests
#> N = 15, T = 38
#> 
#>    statistic p.value
#> CD     3.238   0.001

# Test CCE residuals for CSD
set.seed(1234)
cd_cce <- cd_test(fit_cce, type = "all")
print(cd_cce)
#> Cross-sectional dependence tests
#> N = 15, T = 38
#> 
#>        statistic p.value
#> CD        -2.676   0.007
#> CDw       -1.532   0.126
#> CDw+     150.431   0.000
#> CDstar    -2.835   0.005

# Test DCCE residuals for CSD
set.seed(1234)
cd_dcce <- cd_test(fit_dcce, type = "CDw")
print(cd_dcce)
#> Cross-sectional dependence tests
#> N = 15, T = 35
#> 
#>     statistic p.value
#> CDw     1.092   0.275

# Test CS-ARDL residuals for CSD
set.seed(1234)
cd_csardl <- cd_test(fit_csardl, type = "all")
print(cd_csardl)
#> Cross-sectional dependence tests
#> N = 15, T = 35
#> 
#>        statistic p.value
#> CD        -2.392   0.017
#> CDw        1.092   0.275
#> CDw+     130.969   0.000
#> CDstar    -1.769   0.077
```

**Interpreting Results**:

- **CD statistic p-value \< 0.05**: Reject null of CSD independence;
  residuals are correlated across units.
- **CDw, CDw+, CD\* variants**: Provide robustness checks; if all reject
  the null, CSD is strongly evidenced.
- **Magnitude**: Large $|CD|$ statistics (e.g., $|CD| > 3$) indicate
  substantial and economically meaningful dependence.

In practice, models that do not account for cross-sectional dependence
(like MG without augmentation) typically show significant CD test
rejections, justifying the use of CSD-robust methods like CCE and DCCE.

## References

Chudik, A., & Pesaran, M. H. (2013). Large panel data models with
cross-sectional dependence: A survey \[Globalization Institute Working
Papers\]. *Federal Reserve Bank of Dallas*, (153).

Chudik, A., & Pesaran, M. H. (2015). Common correlated effects
estimation of heterogeneous dynamic panel data models with weakly
exogenous regressors. *Journal of Econometrics*, *188*(2), 393–420.

Ditzen, J. (2018). Estimating dynamic common-correlated effects in
STATA. *The STATA Journal*, *18*(3), 585–617.
<https://doi.org/10.1177/1536867X1801800306>

Fan, J., Liao, Y., & Yao, J. (2015). Power enhancement in
high-dimensional cross-section tests. *Econometrica*, *83*(4),
1497–1541.

Juodis, A., & Reese, S. (2021). The incidental parameters problem in
testing for remaining cross-sectional correlation. *Journal of Business
and Economic Statistics*, *40*(3), 1191–1203.

Pesaran, M. H. (2006). Estimation and inference in large heterogeneous
panels with multifactor error structure. *Econometrica*, *74*(4),
967–1012.

Pesaran, M. H. (2007). A simple unit root test in the presence of
cross-section dependence. *Journal of Applied Econometrics*, *22*(2),
265–312.

Pesaran, M. H. (2015). Testing weak cross-sectional dependence in large
panels. *Econometric Reviews*, *34*(6-10), 1089–1117.

Pesaran, M. H. (2021). General diagnostic tests for cross-sectional
dependence in panels. *Empirical Economics*, *60*(1), 13–50.

Pesaran, M. H., & Smith, R. (1995). Estimating long-run relationships
from dynamic heterogeneous panels. *Journal of Econometrics*, *68*(1),
79–113.

Pesaran, M. H., & Xie, Y. (2021). A bias-corrected CD test for error
cross-sectional dependence in panel models. *Econometric Reviews*,
*41*(6), 649–677.

------------------------------------------------------------------------

*For further details on the theoretical foundations and implementation
of CSD-robust methods, see the documentation for
[`?csdm`](../reference/csdm.md), [`?cd_test`](../reference/cd_test.md),
and [`?summary.csdm_fit`](../reference/summary.csdm_fit.md).*
