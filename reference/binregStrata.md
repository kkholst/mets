# Stratified Binomial Regression with cross-validation option

Fits a stratified binomial regression model for the cumulative incidence
function (CIF), restricted mean survival time (RMST), or restricted mean
time lost (RMTL), using inverse probability of censoring weighting
(IPCW). Each stratum may either be fit on its own observations or, when
`cv.fold = TRUE`, on the complement of its observations (i.e. a
leave-one-stratum-out / cross-validation style fit). Standard errors are
computed via an iid decomposition that accounts for the estimation of
the censoring distribution.

## Usage

``` r
binregStrata(
  formula,
  data,
  cause = 1,
  time = NULL,
  beta = NULL,
  type = c("II", "I"),
  offset = NULL,
  weights = NULL,
  cens.weights = NULL,
  cens.model = ~+1,
  se = TRUE,
  kaplan.meier = TRUE,
  cens.code = 0,
  no.opt = FALSE,
  method = "nr",
  augmentation = NULL,
  outcome = c("cif", "rmst", "rmtl"),
  model = c("default", "logit", "exp", "lin"),
  Ydirect = NULL,
  strata = NULL,
  cv.fold = FALSE,
  low.memory = FALSE,
  ...
)
```

## Arguments

- formula:

  A formula with a [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)
  or `Event` object on the left-hand side, and covariates, optional
  [`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html),
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html),
  [`offset()`](https://rdrr.io/r/stats/offset.html), and
  [`weights()`](https://rdrr.io/r/stats/weights.html) terms on the
  right-hand side.

- data:

  A `data.frame` containing the variables in `formula`.

- cause:

  Value(s) of the status variable that identify the event of interest
  (cause 1 by default).

- time:

  The fixed time-point at which the CIF, RMST, or RMTL is evaluated.
  Required.

- beta:

  Optional starting values for the regression coefficients. Either a
  vector (recycled across strata) or a matrix with one row per stratum.

- type:

  Character; either `"II"` (default) or `"I"`. When `"II"`, an
  additional augmentation term and an extra component of the
  censoring-related iid decomposition are computed.

- offset:

  Optional offset term; can also be specified via the formula.

- weights:

  Optional case weights; can also be specified via the formula.

- cens.weights:

  Optional pre-computed inverse probability of censoring weights. If
  supplied, the internal censoring model fit is skipped and `se` is
  forced to `FALSE`.

- cens.model:

  A formula for the censoring model (right-hand side only, fitted via
  [`phreg`](http://kkholst.github.io/mets/reference/phreg.md)). Defaults
  to the Kaplan-Meier (`~+1`).

- se:

  Logical. If `TRUE` (default), standard errors are computed, including
  the contribution from estimating the censoring distribution.

- kaplan.meier:

  Logical. If `TRUE` (default), uses the Kaplan-Meier estimator for
  censoring weights when the censoring model has no covariates;
  otherwise a Cox-type estimator is used internally when covariates are
  present.

- cens.code:

  The code in the status variable denoting censoring (default 0).

- no.opt:

  Logical. If `TRUE`, no optimization is performed and the object is
  evaluated at the supplied `beta`.

- method:

  Optimization method. Currently only `"nr"` (Newton-Raphson via
  [`lava::NR`](https://kkholst.github.io/lava/reference/NR.html)) is
  supported.

- augmentation:

  Optional augmentation term(s) added to the score equations, either a
  vector (recycled across strata) or a matrix with one row per stratum.

- outcome:

  Character; the outcome scale to model: `"cif"` (default), `"rmst"`, or
  `"rmtl"`.

- model:

  Character; the link function: `"default"` (chooses `"logit"` for CIF
  and `"exp"` for RMST/RMTL), `"logit"`, `"exp"`, or `"lin"`.

- Ydirect:

  Optional, directly supplied response vector, overriding the response
  otherwise constructed from `outcome`.

- strata:

  Optional numeric vector defining strata, used when strata are not
  specified via
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) in the
  formula.

- cv.fold:

  Logical. If `TRUE`, for each stratum `s` the regression coefficients
  are estimated using all observations *not* in stratum `s` (a
  cross-validation-style fit), rather than the observations within
  stratum `s` (the default, `FALSE`).

- low.memory:

  Logical. If `TRUE` only saves key quanties, coef, var, iid and some
  design structures.

- ...:

  Additional arguments, currently passed as optimizer `control` settings
  (e.g. `tol`, `stepsize`) to the Newton-Raphson routine.

## Value

An object of class `c("binregStrata","binreg")`, a list with components
including:

- coef:

  Estimated regression coefficients (stacked across strata).

- iid:

  Subject-level iid decomposition of the coefficients, adjusted for
  censoring estimation when `se = TRUE`.

- var, robvar:

  Robust (sandwich) variance-covariance matrix.

- se.coef, se.robust:

  Standard errors derived from `var`.

- hessian, ihessian:

  Per-stratum Hessians and their (pseudo-)inverses.

- coef_list:

  Per-stratum coefficient vectors.

- time, cause, cens.code, model, outcome:

  Echoed arguments.

- Y, Yipcw, IPCW:

  The constructed response, IPCW-weighted response, and IPCW weights.

- cens.weights, formC, cens.strata, cens.nstrata:

  Censoring model weights and related information (when internally
  estimated).

- nstrata, strata_orig, strata_call, cv.fold:

  Stratification information and whether the cross-validation fold
  scheme was used.

- n, nevent, ncluster, nid, id, name.id:

  Sample size and identifier information.

- call, design:

  The matched call and the design object.

## Details

When `cv.fold = FALSE` (the standard case), the coefficients for stratum
`s` are estimated using only observations with `strata == s`. When
`cv.fold = TRUE`, the coefficients for stratum `s` are instead estimated
using observations with `strata != s`, i.e. the complement; this also
affects the construction of the censoring-related iid terms (`MGCiid`)
and the per-stratum subject-level influence functions used for the
robust variance.

Censoring weights are computed via a Cox model
([`phreg`](http://kkholst.github.io/mets/reference/phreg.md)) on
`cens.model`, evaluated at `pmin(exit, time)`, unless `cens.weights` is
supplied directly. When `se = TRUE`, the additional variability coming
from the estimation of the censoring distribution is incorporated into
the standard errors via a martingale-type iid decomposition (`MGCiid`).

## See also

[`phreg`](http://kkholst.github.io/mets/reference/phreg.md),
[`binreg`](http://kkholst.github.io/mets/reference/binreg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(mets)
data(bmt)
out <- binregStrata(Event(time, cause) ~ tcell + strata(platelet),
                     data = bmt, cause = 1, time = 50)
summary(out)
} # }
```
