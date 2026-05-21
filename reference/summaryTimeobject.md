# Summarize a Time-Varying Estimate with Confidence Bands

Extracts estimates at specified time points with confidence intervals.

## Usage

``` r
summaryTimeobject(
  mutimes,
  mu,
  se.mu = NULL,
  times = NULL,
  type = "log",
  level = 0.95,
  ...
)
```

## Arguments

- mutimes:

  vector of estimation time points.

- mu:

  vector or matrix of estimates.

- se.mu:

  vector or matrix of standard errors (optional).

- times:

  time points at which to evaluate (default: `mutimes`).

- type:

  confidence interval type: `"log"` or `"plain"`.

- level:

  confidence level (default 0.95).

- ...:

  additional arguments.

## Value

A data.frame with columns: times, mean, se-mean, CI-2.5%, CI-97.5%.
