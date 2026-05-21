# Fitting of Two-Stage Recurrent Events Random Effects Model

Fits a two-stage random effects model for recurrent events with a
terminal event. Marginal models (Cox or Ghosh-Lin) are fitted first and
passed to this function.

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

  Marginal model for the terminal event (object of class `"phreg"`).

- recurrent:

  Marginal model for recurrent events (object of class `"phreg"` or
  `"recreg"`).

- data:

  Data frame used for fitting.

- theta:

  Starting value for total variance of gamma frailty.

- model:

  Model type: `"full"` (fully shared), `"shared"` (partly shared), or
  `"non-shared"`.

- ghosh.lin:

  Logical; if `TRUE`, forces use of Ghosh-Lin marginals based on the
  recurrent model.

- theta.des:

  Regression design for variance parameters.

- var.link:

  Link function for variance (1 for exponential).

- method:

  Optimization method (default "NR").

- no.opt:

  Logical; if `TRUE`, skips optimization.

- weights:

  Weights.

- se.cluster:

  Clusters for SE calculation (GEE style).

- fnu:

  Function to transform \\\nu\\ (amount shared).

- nufix:

  Logical; if `TRUE`, fixes the amount shared.

- nu:

  Starting value for the amount shared.

- numderiv:

  Logical; if `TRUE`, uses numerical derivatives.

- derivmethod:

  Method for numerical derivative.

- ...:

  Arguments for the optimizer.

## Value

An object of class `"twostageREC"` containing:

- coef:

  Estimated coefficients.

- var:

  Variance-covariance matrix.

- theta:

  Dependence parameters.

- model:

  Model type.

## Details

Supports:

- Cox/Cox marginals.

- Cox/Ghosh-Lin marginals.

- Fully shared, partly shared, or non-shared random effects.

## References

Scheike (2026), Two-stage recurrent events random effects models, LIDA,
to appear.

## Author

Thomas Scheike
