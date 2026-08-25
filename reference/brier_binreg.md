# Cross-validated Brier score for competing risks, RMST and RMTL regression using stratified leave-fold-out fits (`binreg` `binregStrata`)

Fits candidate models via
[`binregStrata`](http://kkholst.github.io/mets/reference/binregStrata.md)
with `cv.fold = TRUE` and evaluates each using apparent (in-sample) and
*k*-fold cross-validated IPCW Brier scores relative to a null model.
Supports all outcome types available in `binreg`/`binregStrata`:
cause-specific cumulative incidence (competing risks), restricted mean
survival time (RMST), and restricted mean time lost (RMTL).

## Usage

``` r
brier_binreg(
  formula,
  data,
  time = NULL,
  fold = 5,
  rhs = NULL,
  rhs0 = ~+1,
  cens.model = ~+1,
  brier.cens.model = ~+1,
  brier.args = list(),
  cv.split.index = NULL,
  ...
)
```

## Arguments

- formula:

  Two-sided formula of the form
  `Event(time, status) ~ covariates + cluster(id)` or
  `Event(time, status, cause) ~ covariates + cluster(id)`. The
  `cluster(id)` term is optional; if absent an internal identifier
  `id__` is added automatically assuming independent observations.

- data:

  A `data.frame` containing all variables in `formula` and `rhs`.

- time:

  Scalar or numeric vector of evaluation time points. When a vector is
  supplied the function loops over time points reusing the same fold
  structure and returns a named list of `binregCV` objects with class
  `"binregCV_list"`.

- fold:

  Integer. Number of cross-validation folds / strata (default 5).

- rhs:

  Named list of one-sided formulas specifying candidate models, e.g.
  `list(small = ~Z1, full = ~Z1 + Z2)`. If `NULL` (default), the
  right-hand side of `formula` (minus any
  [`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html) term) is
  used as a single candidate named `"model1"`.

- rhs0:

  One-sided formula for the null (intercept-only) model (default `~+1`).

- cens.model:

  Censoring model passed to `binregStrata` for the CV fold fits (default
  `~+1`). Since folds are handled internally by `binregStrata`'s
  `strata`/`cv.fold` mechanism, this does *not* need to include
  `strata(fold)`.

- brier.cens.model:

  Censoring model passed to `resmeanIPCW` when computing Brier scores.
  Default `~+1` gives Kaplan-Meier weights. Supply a stratified model,
  e.g. `~strata(Z1, Z2)`, for covariate-adjusted IPCW weights.

- brier.args:

  Named list of additional arguments forwarded to every `resmeanIPCW`
  call (beyond `formula`, `data`, `time`, `Ydirect`, `cens.model`, and
  `model` which are set internally). Typically empty.

- cv.split.index:

  List from `folds(n, fold)` or a vector of fold indices (1:fold) for
  each observation, to fix the cross-validation folds. `NULL` (default)
  generates folds internally via `folds(n, fold)`.

- ...:

  Additional arguments passed to every `binreg` / `binregStrata` call,
  e.g. `cause = 2`, `outcome = "rmst"`, or `offset`.

## Value

An object of class `"binregCV"` (single time point) or `"binregCV_list"`
(multiple time points), with the following structure:

- `null`:

  List with `coef`, `se`, and `iid` for the null model. `coef` and `se`
  are named vectors with entries `bbrier` (apparent, log scale) and
  `bbrierCV` (cross-validated, log scale). `iid` contains the
  influence-function matrix for `bbrierCV` with delta-method correction
  for the stratified CV prediction step and the log transform.

- `models`:

  Named list, one element per entry in `rhs`. Each element contains
  `coef` (`brier`, `brierCV`, log scale), `se`, `dcoef` (differences
  from null: log Brier ratios), `sed`, and `iid` (influence-function
  matrix for `brierCV`, log scale, with CV delta-method correction).

- `iid_mat`:

  Matrix of dimension \\n \times (1 + \code{length(rhs)})\\ with columns
  `null` and one per entry of `rhs`, each the log-scale CV Brier
  influence function for the corresponding model. Used by
  [`estimate.binregCV`](http://kkholst.github.io/mets/reference/estimate.binregCV.md)
  for joint inference.

- `cv.split.index`:

  List of fold index vectors used.

- `rhs`, `rhs0`:

  Processed candidate and null formulas (cluster terms removed).

- `model`:

  Link function read from the `binreg` fit: `"logit"`, `"exp"`, or
  `"lin"`.

- `call`, `time`, `fold`, `cens.model`, `brier.cens.model`:

  Call and settings for printing and downstream use.

## Details

Cross-validation is delegated directly to
[`binregStrata`](http://kkholst.github.io/mets/reference/binregStrata.md)
via its `strata` and `cv.fold = TRUE` arguments: each stratum (fold)
gets its own coefficient vector \\\beta_s\\, estimated on the complement
of that fold. Out-of-fold predictions are obtained by routing each
observation to the \\\beta_s\\ corresponding to its own fold, which was
trained on the other folds.

The fold membership is stored in an internal `fold` column (1-based) in
`data` and passed to `binregStrata` as a 0-based `strata` vector.
Censoring weights for the underlying competing-risks model are estimated
according to `cens.model` as usual; since strata already define the
folds, `cens.model` need not include `strata(fold)` explicitly.

Brier scores are computed on the **linear scale** via
`resmeanIPCW(..., model = "lin")` and then log-transformed using
`lava::estimate(..., f = log)` with the delta method. All stored
coefficients, standard errors, and influence functions are therefore on
the **log scale**. Use `summary(fit, transform = exp)` or
`estimate(fit, f = exp)` to back-transform to the original Brier scale.

The delta-method correction for the cross-validated influence function
accounts for the block structure of `binregStrata`'s influence function
(\\n\_{id} \times (\code{nstrata} \times p)\\, one \\p\\-column block
per stratum). For each stratum \\s\\, only observations belonging to
fold \\s\\ have a non-zero derivative of their out-of-fold prediction
with respect to \\\beta_s\\, since \\\beta_s\\ is exactly the
coefficient used to predict fold \\s\\. The full gradient is therefore
block-sparse, with one non-zero \\p\\-length block per stratum, and the
correction reduces to a single matrix product against `iid(outff_fit)`.

Left-truncation (delayed entry) is not supported and will raise an
error. Specifically, `Event(entry, time, status)` is rejected; use
`Event(time, status)` or `Event(time, status, cause)`.

## Outcome types

The function supports all outcome types provided by `binreg`/
`binregStrata`:

- Competing risks CIF:

  Default; use `Event(time, status, cause)` and pass `cause` via `...`.
  The Brier score measures squared error for the cause-specific CIF.

- RMST:

  Pass `outcome = "rmst"` via `...`. The Brier score measures squared
  error for the restricted mean survival time up to `time`.

- RMTL:

  Pass `outcome = "rmtl"` via `...`. The Brier score measures squared
  error for the restricted mean time lost.

## Scale and inference

All output is on the log scale. To obtain Brier scores on the original
scale with confidence intervals, use:


      summary(fit, transform = exp)
      estimate(fit, f = exp)

Differences on the log scale correspond to log Brier ratios;
exponentiate to obtain Brier ratios. For absolute differences on the
original scale, compute contrasts after `estimate(fit, f = exp)`.

## See also

[`estimate.binregCV`](http://kkholst.github.io/mets/reference/estimate.binregCV.md)
for joint inference across models and time points;
[`summary.binregCV`](http://kkholst.github.io/mets/reference/summary.binregCV.md)
and
[`print.summary_binregCV_multi`](http://kkholst.github.io/mets/reference/print.summary_binregCV_multi.md)
for formatted output;
[`binregStrata`](http://kkholst.github.io/mets/reference/binregStrata.md)
for the underlying stratified regression model;
[`resmeanIPCW`](http://kkholst.github.io/mets/reference/resmeanIPCW.md)
for IPCW mean estimation;
[`estimate`](https://kkholst.github.io/lava/reference/estimate.default.html)
for delta-method inference.

## Examples

``` r
if (FALSE) { # \dontrun{
library(mets)
data(bmt)
bmt$id <- seq_len(nrow(bmt))

## --- competing risks CIF/Survival (default) ---
fit <- brier_binreg(
  Event(time, cause) ~ tcell + platelet + age + cluster(id),
  data = bmt, time = 50,
  rhs = list(small = ~age, full = ~tcell + platelet + age)
)
summary(fit)                        ## log scale
summary(fit, transform = exp)       ## Brier scale

## multiple time points
fitm <- brier_binreg(
  Event(time, cause) ~ tcell + platelet + age + cluster(id),
  data = bmt, time = c(50, 60),
  rhs = list(small = ~age, full = ~tcell + platelet + age)
)
summary(fitm, transform = exp)

## joint inference across models via lava::estimate
e <- estimate(fitm)
e
estimate(e, f = exp)               ## back-transform

## stratified IPCW weights for Brier score
fit_strat <- brier_binreg(
  Event(time, cause) ~ tcell + platelet + age + cluster(id),
  data = bmt, time = 50,
  rhs = list(small = ~age, full = ~tcell + platelet + age),
  brier.cens.model = ~strata(tcell)
)

## --- RMST --- survival case only
fit_rmst <- brier_binreg(
  Event(time, cause) ~ tcell + platelet + age + cluster(id),
  data = bmt, time = 50, cause = 1:2,
  rhs    = list(small = ~age, full = ~tcell + platelet + age),
  outcome = "rmst"
)
summary(fit_rmst, transform = exp)

## --- RMTL --- competing risks or survival (RMTL)
fit_rmtl <- brier_binreg(
  Event(time, cause) ~ tcell + platelet + age + cluster(id),
  data = bmt, time = 50,
  rhs    = list(small = ~age, full = ~tcell + platelet + age),
  outcome = "rmtl"
)
summary(fit_rmtl, transform = exp)

## --- no IPCW adjustment, regression ---------
fit_reg <- brier_binreg(
  time ~ tcell + platelet + age + cluster(id),
  data = bmt, 
  rhs    = list(small = ~age, full = ~tcell + platelet + age)
)
summary(fit_reg, transform = exp)

## --- no IPCW adjustment  binomial-regression  ---------------
bmt$binY  <- as.numeric(bmt$time>30)
fit_bin <- brier_binreg(
  binY ~ tcell + platelet + age + cluster(id),
  data = bmt, model="logit",
  rhs    = list(small = ~age, full = ~tcell + platelet + age)
)
summary(fit_bin, transform = exp)


} # }
```
