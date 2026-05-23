# Extend Cumulative Hazard Functions to Common Time Range

Extends a collection of cumulative hazard functions so that they all
cover the same time range. Cumulative hazards that end before the
maximum observed time are extrapolated linearly using either a
user-supplied rate or the average hazard rate estimated from the
existing cumulative hazard.

## Usage

``` r
extendCums(cumA, cumB, extend = NULL)
```

## Arguments

- cumA:

  a cumulative hazard matrix (two columns: time, cumulative hazard) or a
  list of such matrices.

- cumB:

  an optional second cumulative hazard matrix. If provided it is
  combined with `cumA` into a list.

- extend:

  numeric vector of hazard rates used for extrapolation, one per
  cumulative hazard. If `NULL` (default), rates are estimated as the
  ratio of the last cumulative hazard value to the last time point.

## Value

A list of cumulative hazard matrices extended to the common maximum
time.

## See also

[`rchaz`](http://kkholst.github.io/mets/reference/rchaz.md),
[`rcrisk`](http://kkholst.github.io/mets/reference/rcrisk.md),
[`mets-simulation`](http://kkholst.github.io/mets/reference/mets-simulation.md)
