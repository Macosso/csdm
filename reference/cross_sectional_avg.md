# Cross-sectional averages by time (with optional leave-one-out)

Computes cross-sectional averages (CSAs) of specified variables for each
time period, optionally in a leave-one-out (LOO) fashion per
observation. Supports unbalanced panels and observation weights.

## Usage

``` r
cross_sectional_avg(
  data,
  id = NULL,
  time = NULL,
  vars,
  leave_out = FALSE,
  weights = NULL,
  suffix = "csa",
  return_mode = c("attach", "time"),
  na.rm = TRUE
)
```

## Arguments

- data:

  A `data.frame` or
  [`plm::pdata.frame`](https://rdrr.io/pkg/plm/man/pdata.frame.html).

- id, time:

  Character scalar names of unit and time columns when `data` is a plain
  `data.frame`. If `data` is a `pdata.frame`, these are inferred from
  its index and can be omitted.

- vars:

  Character vector of column names to average cross-sectionally.

- leave_out:

  Logical; if `TRUE`, computes LOO means for each row: \\\bar{x}\_{-i,t}
  = (\sum\_{j \neq i} w\_{jt} x\_{jt}) / (\sum\_{j \neq i} w\_{jt})\\.
  If `FALSE`, computes standard time means: \\\bar{x}\_{t} = (\sum_j
  w\_{jt} x\_{jt}) / (\sum_j w\_{jt})\\.

- weights:

  Optional. Either:

  - a numeric vector of length `nrow(data)`, or

  - the name of a column in `data` with nonnegative weights.

  If `NULL`, uses equal weights (1 for observed values).

- suffix:

  Character suffix to append to CSA columns (default `"csa"`).

- return_mode:

  One of `"attach"` or `"time"`.

  - `"attach"` returns the original `data` with CSA columns added.

  - `"time"` returns a unique-time table `[time, csa_*]`.

- na.rm:

  Logical; if `TRUE`, excludes `NA`s from sums and denominators. If
  `FALSE`, any `NA` in a time slice yields `NA` for that time's CSA for
  that variable.

## Value

A `data.frame`:

- If `return_mode="attach"`: original data + CSA columns named
  `paste0(suffix, "_", vars)`.

- If `return_mode="time"`: unique time rows with CSA columns.

## Details

Efficiently computes, for each `v in vars` and time `t`, \$\$\bar v_t =
\frac{\sum_i w\_{it}\\ 1\_{\\v\_{it}\text{ finite}\\}\\ v\_{it}} {\sum_i
w\_{it}\\ 1\_{\\v\_{it}\text{ finite}\\}}\$\$ For `leave_out=TRUE`, each
row's CSA excludes its own contribution; if the denominator becomes
\\\le 0\\ (e.g., only one finite observation at that time), the LOO mean
is set to `NA` for that row/variable.
