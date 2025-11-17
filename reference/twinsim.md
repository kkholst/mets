# Simulate twin data

Simulate twin data from a linear normal ACE/ADE/AE model.

## Usage

``` r
twinsim(
  nMZ = 100,
  nDZ = nMZ,
  b1 = c(),
  b2 = c(),
  mu = 0,
  acde = c(1, 1, 0, 1),
  randomslope = NULL,
  threshold = 0,
  cens = FALSE,
  wide = FALSE,
  ...
)
```

## Arguments

- nMZ:

  Number of monozygotic twin pairs

- nDZ:

  Number of dizygotic twin pairs

- b1:

  Effect of covariates (labelled x1,x2,...) of type 1. One distinct
  covariate value for each twin/individual.

- b2:

  Effect of covariates (labelled g1,g2,...) of type 2. One covariate
  value for each twin pair.

- mu:

  Intercept parameter.

- acde:

  Variance of random effects (in the order A,C,D,E)

- randomslope:

  Logical indicating wether to include random slopes of the variance
  components w.r.t. x1,x2,...

- threshold:

  Treshold used to define binary outcome y0

- cens:

  Logical variable indicating whether to censor outcome

- wide:

  Logical indicating if wide data format should be returned

- ...:

  Additional arguments parsed on to lower-level functions

## See also

[`twinlm`](http://kkholst.github.io/mets/reference/twinlm.md)

## Author

Klaus K. Holst
