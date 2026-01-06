# Fittting of Two-stage recurrent events random effects model based on Cox/Cox or Cox/Ghosh-Lin marginals

Fittting of Two-stage recurrent events random effects model based on
Cox/Cox or Cox/Ghosh-Lin marginals. Random effects model fore recurrent
events with terminal event. Marginal models fitted first and given to
twostageREC function.

## Usage

``` r
twostageREC(
  margsurv,
  recurrent,
  data = parent.frame(),
  theta = NULL,
  model = c("full", "shared", "non-shared"),
  ghosh.lin = NULL,
  theta.des = NULL,
  var.link = 0,
  method = "NR",
  no.opt = FALSE,
  weights = NULL,
  se.cluster = NULL,
  fnu = NULL,
  nufix = 0,
  nu = NULL,
  numderiv = 1,
  derivmethod = c("simple", "Richardson"),
  ...
)
```

## Arguments

- margsurv:

  marginal model for terminal event

- recurrent:

  marginal model for recurrent events

- data:

  used for fitting

- theta:

  starting value for total variance of gamma frailty

- model:

  can fully shared "full", partly shared "shared", or non-shared where
  the random effect acts only on recurrent events

- ghosh.lin:

  to force use ghosh.lin marginals based on recurrent (taking baseline
  and coefficients)

- theta.des:

  regression design for variance

- var.link:

  possible link function 1 is exponential link

- method:

  NR

- no.opt:

  to not optimize

- weights:

  possible weights

- se.cluster:

  to combine influence functions for naive variance based on these
  clusters GEE style

- fnu:

  a function to make transformation for nu (amount shared)

- nufix:

  to fix the amount shared

- nu:

  starting value for amount shared

- numderiv:

  uses numerical derivatives for some derivatives

- derivmethod:

  method for numerical derivative

- ...:

  arguments for

## References

Scheike (2026), Two-stage recurrent events random effects models, LIDA,
to appear
