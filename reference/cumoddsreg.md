# Cumulative Odds Regression for Discrete Time Data

A wrapper function for `interval_logitsurv_discrete` that simplifies the
interface for discrete time-to-event data where the event time is
observed exactly or as a factor level. It converts a factor response
into an interval-censored format internally.

## Usage

``` r
cumoddsreg(formula, data, ...)
```

## Arguments

- formula:

  Formula with a factor response on the left-hand side (representing the
  event time) and covariates on the right.

- data:

  Data frame.

- ...:

  Arguments passed to `interval_logitsurv_discrete`.

## Value

An object of class `"cumoddsreg"` with the same structure as
`interval_logitsurv_discrete`.

## See also

[`interval_logitsurv_discrete`](http://kkholst.github.io/mets/reference/interval_logitsurv_discrete.md)

## Author

Thomas Scheike
