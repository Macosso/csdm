<!-- badges: start -->
[![R-CMD-check](https://github.com/Macosso/csdm/workflows/R-CMD-check/badge.svg)](https://github.com/Macosso/csdm/actions)
<!-- badges: end -->

# csdm

`csdm` provides Mean Group (MG), Common Correlated Effects (CCE), and Dynamic CCE (DCCE) estimators via `csdm()`.

This package does **not** execute Stata code and does **not** claim numerical equivalence to `xtdcce2`. The goal of the mapping below is to provide a conceptual translation and runnable R usage.

## Replicating xtdcce2 syntax (Stata) in csdm (R)

### Quick mapping

| Stata concept / xtdcce2 option (descriptive) | `csdm` syntax |
|---|---|
| Dynamic dependent-variable lag (e.g., `L.y`) | `lr = csdm_lr(type="ardl", ylags=1)` |
| Cross-sectional averages controls | `csa = csdm_csa(vars=c(...), lags=...)` |
| CSA lags (e.g., 3) | `csa = csdm_csa(..., lags=3)` |
| Scalar distributed lags for all RHS regressors | `lr = csdm_lr(type="ardl", xdlags=K)` |
| Pesaran CD test | `cd_test(fit)` |

### Worked end-to-end example (DCCE + CSA + y-lag)

```r
library(csdm)
data(PWT_60_07, package="csdm")
df <- PWT_60_07

fit <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df,
  id = "id",
  time = "year",
  model = "dcce",
  csa = csdm_csa(
    vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"),
    lags = 3
  ),
  lr = csdm_lr(type = "ardl", ylags = 1)
)

summary(fit)
cd_test(fit)
```

### Distributed lags for RHS regressors (`xdlags`)

```r
fit_dx <- csdm(
  log_rgdpo ~ log_hc + log_ck + log_ngd,
  data = df,
  id = "id",
  time = "year",
  model = "dcce",
  csa = csdm_csa(
    vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"),
    lags = 3
  ),
  lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 1)
)

names(coef(fit_dx))
```

### Notes

- `xdlags` currently supports only **simple RHS variable names** (no transformations or interactions like `log(x)` or `x1:x2`).
- Generated lag columns follow the naming convention `lag1_<varname>`, `lag2_<varname>`, ...
