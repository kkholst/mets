# Influence Functions for phreg objects

Computes the influence functions (IID decomposition) for the regression
coefficients and/or the baseline cumulative hazard at a specific time
point.

## Usage

``` r
# S3 method for class 'phreg'
IC(x, type = "robust", all = FALSE, time = NULL, baseline = NULL, ...)
```

## Arguments

- x:

  Object of class `"phreg"`.

- type:

  Type of influence function: `"robust"` (default) or `"martingale"`.

- all:

  Logical; if `TRUE`, returns both beta and baseline influence
  functions.

- time:

  Time point for baseline influence function (required if baseline is
  requested).

- baseline:

  Arguments for baseline estimation.

- ...:

  Additional arguments.

## Value

A matrix of influence functions. If `all=TRUE`, columns correspond to
regression coefficients and baseline cumulative hazard. Attributes
include `coef` and `time`.

## See also

[`phreg`](http://kkholst.github.io/mets/reference/phreg.md),
[`iidBaseline`](http://kkholst.github.io/mets/reference/iidBaseline.md)

## Author

Thomas Scheike
