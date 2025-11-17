# Simulate from the Aalen Frailty model

Simulate observations from Aalen Frailty model with Gamma distributed
frailty and constant intensity.

## Usage

``` r
simAalenFrailty(
  n = 5000,
  theta = 0.3,
  K = 2,
  beta0 = 1.5,
  beta = 1,
  cens = 1.5,
  cuts = 0,
  ...
)
```

## Arguments

- n:

  Number of observations in each cluster

- theta:

  Dependence paramter (variance of frailty)

- K:

  Number of clusters

- beta0:

  Baseline (intercept)

- beta:

  Effect (log hazard ratio) of covariate

- cens:

  Censoring rate

- cuts:

  time cuts

- ...:

  Additional arguments

## Author

Klaus K. Holst
