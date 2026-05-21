# Force Same Censoring Within Clusters

Enforces the same censoring time within clusters (pairs) by censoring
both subjects at the minimum censoring time when censoring occurs first.

## Usage

``` r
force.same.cens(
  data,
  id = "id",
  time = "time",
  cause = "cause",
  entrytime = NULL,
  cens.code = 0
)

force_same_cens(
  data,
  id = "id",
  time = "time",
  cause = "cause",
  entrytime = NULL,
  cens.code = 0
)
```

## Arguments

- data:

  a data.frame in long format.

- id:

  name of the cluster/pair identifier column.

- time:

  name of the time variable.

- cause:

  name of the cause/status variable.

- entrytime:

  optional name of left-truncation time variable.

- cens.code:

  value indicating censoring in the cause variable.

## Value

A data.frame with enforced same censoring.
