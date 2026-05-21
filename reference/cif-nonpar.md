# Non-parametric Cumulative Incidence Functions

Functions for computing and visualizing non-parametric cumulative
incidence estimates, as well as dependence measures (odds ratio,
relative risk) for bivariate competing risks data.

## Usage

``` r
rr_cif(
  cif,
  data,
  cause = NULL,
  cif2 = NULL,
  times = NULL,
  cause1 = 1,
  cause2 = 1,
  cens.code = NULL,
  cens.model = "KM",
  Nit = 40,
  detail = 0,
  clusters = NULL,
  theta = NULL,
  theta.des = NULL,
  step = 1,
  sym = 0,
  weights = NULL,
  same.cens = FALSE,
  censoring.weights = NULL,
  silent = 1,
  par.func = NULL,
  dpar.func = NULL,
  dimpar = NULL,
  score.method = "nlminb",
  entry = NULL,
  estimator = 1,
  trunkp = 1,
  admin.cens = NULL,
  ...
)

or_cif(
  cif,
  data,
  cause = NULL,
  cif2 = NULL,
  times = NULL,
  cause1 = 1,
  cause2 = 1,
  cens.code = NULL,
  cens.model = "KM",
  Nit = 40,
  detail = 0,
  clusters = NULL,
  theta = NULL,
  theta.des = NULL,
  step = 1,
  sym = 0,
  weights = NULL,
  same.cens = FALSE,
  censoring.weights = NULL,
  silent = 1,
  par.func = NULL,
  dpar.func = NULL,
  dimpar = NULL,
  score.method = "nlminb",
  entry = NULL,
  estimator = 1,
  trunkp = 1,
  admin.cens = NULL,
  ...
)

random.cif(cif, ...)

Grandom.cif(cif, ...)

predictPairPlack(cif1, cif2, status1, status2, theta)

npc(T, cause, same.cens = TRUE, sep = FALSE)

nonparcuminc(t, status, cens = 0)

plotcr(
  x,
  col,
  lty,
  legend = TRUE,
  which = 1:2,
  cause = 1:2,
  ask = prod(par("mfcol")) < length(which) && dev.interactive(),
  ...
)
```

## Arguments

- cif:

  a cumulative incidence model object (from timereg).

- data:

  a data.frame with the variables.

- cause:

  causes to plot.

- cif2:

  optional second CIF model if different from first.

- times:

  time points for evaluation.

- cause1:

  cause for first coordinate.

- cause2:

  cause for second coordinate.

- cens.code:

  censoring code value.

- cens.model:

  censoring model type (default `"KM"`).

- Nit:

  maximum number of iterations.

- detail:

  level of output detail.

- clusters:

  cluster variable name or vector.

- theta:

  dependence parameter(s).

- theta.des:

  design matrix for theta.

- step:

  step size for optimization.

- sym:

  if 1, symmetric dependence structure.

- weights:

  optional weights.

- same.cens:

  logical; if TRUE, uses joint censoring weights.

- censoring.weights:

  optional pre-computed censoring weights.

- silent:

  verbosity level.

- par.func:

  optional parameter function.

- dpar.func:

  optional derivative of parameter function.

- dimpar:

  dimension of parameter vector.

- score.method:

  optimization method (default `"nlminb"`).

- entry:

  optional entry time variable.

- estimator:

  estimator type.

- trunkp:

  truncation probability.

- admin.cens:

  administrative censoring time.

- ...:

  additional arguments.

- cif1:

  CIF values for subject 1 (for `predictPairPlack`).

- status1:

  status for subject 1.

- status2:

  status for subject 2.

- T:

  matrix with columns: time1, time2, status1, status2 (for `npc`).

- sep:

  logical; if TRUE, uses separate censoring models for each subject.

- t:

  vector of event/censoring times (for `nonparcuminc`).

- status:

  vector of status codes (for `nonparcuminc`).

- cens:

  censoring code (default 0).

- x:

  data matrix or competing risks object.

- col:

  colors for curves.

- lty:

  line types for curves.

- legend:

  logical; if TRUE, add legend.

- which:

  which plots to show.

- ask:

  logical; if TRUE, prompt before new page.

## Value

For `npc`: matrix with columns (time, cumulative incidence). For
`nonparcuminc`: matrix with time and cause-specific cumulative
incidences.

## Details

`npc` computes bivariate non-parametric cumulative incidence using
inverse-probability-of-censoring weights.

`nonparcuminc` computes univariate non-parametric cumulative incidence
for multiple causes.

`plotcr` plots cumulative incidence curves for competing risks using the
prodlim package.

`or_cif` fits an odds-ratio model for bivariate cumulative incidence.

`rr_cif` fits a relative-risk model for bivariate cumulative incidence.

`random.cif` and `Grandom.cif` are aliases for `random_cif` and
`Grandom_cif` (random effects CIF models).

`predictPairPlack` predicts pairwise joint probabilities under a
Plackett (odds-ratio) dependence model.

`matplot.mets.twostage` produces matrix-plots of concordance over time
from a twostage object.

## Author

Klaus K. Holst, Thomas Scheike
