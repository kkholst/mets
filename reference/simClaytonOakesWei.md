# Simulate from the Clayton-Oakes frailty model

Simulate observations from the Clayton-Oakes copula model with Weibull
type baseline and Cox marginals.

## Usage

``` r
simClaytonOakesWei(
  K,
  n,
  eta,
  beta,
  stoptime,
  weiscale = 1,
  weishape = 2,
  left = 0,
  pairleft = 0
)
```

## Arguments

- K:

  Number of clusters

- n:

  Number of observations in each cluster

- eta:

  1/variance

- beta:

  Effect (log hazard ratio) of covariate

- stoptime:

  Stopping time

- weiscale:

  weibull scale parameter

- weishape:

  weibull shape parameter

- left:

  Left truncation

- pairleft:

  pairwise (1) left truncation or individual (0)

## Author

Klaus K. Holst
