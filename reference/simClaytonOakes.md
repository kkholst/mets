# Simulate from the Clayton-Oakes frailty model

Simulate observations from the Clayton-Oakes copula model with piecewise
constant marginals.

## Usage

``` r
simClaytonOakes(
  K,
  n,
  eta,
  beta,
  stoptime,
  lam = 1,
  left = 0,
  pairleft = 0,
  trunc.prob = 0.5,
  same = 0
)
```

## Arguments

- K:

  Number of clusters

- n:

  Number of observations in each cluster

- eta:

  variance

- beta:

  Effect (log hazard ratio) of covariate

- stoptime:

  Stopping time

- lam:

  constant hazard

- left:

  Left truncation

- pairleft:

  pairwise (1) left truncation or individual (0)

- trunc.prob:

  Truncation probability

- same:

  if 1 then left-truncation is same also for univariate truncation

## Author

Thomas Scheike and Klaus K. Holst
