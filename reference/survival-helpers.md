# Survival Twostage Helpers

Helper functions for the twostage survival dependence models.

## Usage

``` r
survival.twostage(x, ...)

matplot.mets.twostage(object, ...)

alpha2spear(theta, link = 1)

alpha2kendall(theta, link = 0)

piecewise_twostage(
  cut1,
  cut2,
  data = parent.frame(),
  timevar = "time",
  status = "status",
  id = "id",
  covars = NULL,
  covars.pairs = NULL,
  num = NULL,
  method = "optimize",
  Nit = 100,
  detail = 0,
  silent = 1,
  weights = NULL,
  control = list(),
  theta = NULL,
  theta.des = NULL,
  var.link = 1,
  step = 0.5,
  model = "plackett",
  data.return = 0
)

piecewise_data(
  cut1,
  cut2,
  data = parent.frame(),
  timevar = "time",
  status = "status",
  id = "id",
  covars = NULL,
  covars.pairs = NULL,
  num = NULL,
  silent = 1
)
```

## Arguments

- x:

  a marginal model object (for `survival.twostage`).

- ...:

  additional arguments.

- object:

  a twostage model object (for matplot method).

- theta:

  initial dependence parameter values.

- link:

  if 1, parameters are on log scale (for alpha2kendall/alpha2spear).

- cut1:

  vector of cut points for the first time axis.

- cut2:

  vector of cut points for the second time axis.

- data:

  a data.frame with the survival data.

- timevar:

  character name of the time variable.

- status:

  character name of the status variable.

- id:

  name of the cluster identifier column.

- covars:

  optional character vector of covariate names.

- covars.pairs:

  optional covariates at the pair level.

- num:

  character name of the within-cluster number variable.

- method:

  optimization method.

- Nit:

  maximum number of iterations.

- detail:

  level of detail in output.

- silent:

  level of verbosity (1=silent).

- weights:

  optional weights.

- control:

  optimization control list.

- theta.des:

  theta design matrix.

- var.link:

  if 1, log-link for variance parameters.

- step:

  step size for optimization.

- model:

  dependence model: `"plackett"` or `"clayton.oakes"`.

- data.return:

  if 1, return data with model fits.

## Details

`survival.twostage` is an alias for `survival_twostage`.

`alpha2kendall` converts the Clayton-Oakes alpha parameter to Kendall's
tau.

`alpha2spear` converts the Clayton-Oakes alpha parameter to Spearman's
rho.

`piecewise_twostage` fits twostage models on piecewise time intervals.

`piecewise_data` prepares data for piecewise twostage analysis.

`matplot.mets.twostage` produces matplot of twostage baseline estimates.

## Author

Klaus K. Holst, Thomas Scheike
