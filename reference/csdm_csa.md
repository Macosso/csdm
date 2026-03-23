# Specification: Cross-sectional averages (CSA)

Specification: Cross-sectional averages (CSA)

## Usage

``` r
csdm_csa(
  vars = "_all",
  lags = 0,
  scope = c("estimation", "global", "cluster"),
  cluster = NULL
)
```

## Arguments

- vars:

  Character. One of "\_all", "\_none", or a character vector of variable
  names.

- lags:

  Integer. Either a scalar integer \>= 0 applied to all CSA variables,
  or a named integer vector giving per-variable maximum lags.

- scope:

  Character vector. One or more of c("estimation","global","cluster").

- cluster:

  Reserved for future use.

## Value

A spec object (list) used by csdm().

## Examples

``` r
# Cross-sectional averages (CSA) configuration for DCCE
csa <- csdm_csa(
  vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"),
  lags = 3
)
csa
#> $vars
#> [1] "log_rgdpo" "log_hc"    "log_ck"    "log_ngd"  
#> 
#> $lags
#> [1] 3
#> 
#> $scope
#> [1] "estimation" "global"     "cluster"   
#> 
#> $cluster
#> NULL
#> 
#> attr(,"class")
#> [1] "csdm_csa_spec"
```
