# Robust Baseline Hazard Standard Errors

Computes robust (sandwich) standard errors for the cumulative baseline
hazard from a `phreg` object.

Summarizes cumulative baseline hazard estimates from a phreg object,
optionally with robust standard errors.

Computes confidence intervals using log or plain transformations, with
optional restrictions to positive values or probability scale.

## Usage

``` r
robust.basehaz.phreg(x, type = "robust", fixbeta = NULL, ...)

summarybase.phreg(object, robust = FALSE, ...)

conftype(
  x,
  std.err,
  conf.type = c("log", "plain"),
  restrict = c("positive", "prob", "none"),
  conf.int = 0.95
)
```

## Arguments

- x:

  point estimate(s).

- type:

  type of standard error (default `"robust"`).

- fixbeta:

  if non-NULL, fixes beta at given value.

- ...:

  additional arguments.

- object:

  a `phreg` object.

- robust:

  logical; if TRUE, uses robust standard errors.

- std.err:

  standard error(s).

- conf.type:

  type of transformation: `"log"` or `"plain"`.

- restrict:

  restriction: `"positive"`, `"prob"`, or `"none"`.

- conf.int:

  confidence level (default 0.95).

## Value

A list with `cumhaz`, `se.cumhaz`, and `strata`.

An object of class `"summary.recurrent"`.

A list with `upper`, `lower`, `conf.type`, and `conf.int`.
