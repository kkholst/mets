# Bivariate Survival Data on Rectangular Regions

Restricts bivariate survival data to a rectangular time region defined
by left-truncation and right-censoring boundaries, for use with
piecewise twostage models.

## Usage

``` r
surv_boxarea(
  left.trunc,
  right.cens,
  data,
  timevar = "time",
  status = "status",
  id = "id",
  covars = NULL,
  covars.pairs = NULL,
  num = NULL,
  silent = 1,
  boxtimevar = "boxtime"
)
```

## Arguments

- left.trunc:

  vector of length 2 giving left truncation times.

- right.cens:

  vector of length 2 giving right censoring times.

- data:

  a data.frame with the survival data.

- timevar:

  name of the time variable.

- status:

  name of the status variable.

- id:

  name of the cluster identifier.

- covars:

  optional covariate names.

- covars.pairs:

  optional pair-level covariate names.

- num:

  within-cluster numbering variable name.

- silent:

  verbosity level (1=silent).

- boxtimevar:

  name for the created box-time variable.

## Value

A data.frame in long format restricted to the specified region.
