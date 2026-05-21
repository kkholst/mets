# Simulate recurrent events from a two-stage Cox or Ghosh-Lin model

Simulates recurrent event data from a fitted two-stage model, where the
recurrent event process and the terminal event are each described by a
separate fitted model. The recurrent event model may be either a Cox
proportional hazards model (`phreg`) or a Ghosh-Lin marginal rate model
(`recreg`); the terminal event model must be a Cox model (`phreg`).

## Usage

``` r
sim_recurrent_ts(
  cox1,
  coxd = NULL,
  n = 1,
  data = NULL,
  type = c("default", "cox-cox", "gl-cox"),
  id = "id",
  varz = 1,
  share = 1,
  cens = 0.001,
  scale1 = 1,
  scaled = 1,
  dependence = NULL,
  r1 = NULL,
  rd = NULL,
  rc = NULL,
  strata1 = NULL,
  stratad = NULL,
  death.code = 3,
  ...
)
```

## Arguments

- cox1:

  A fitted `phreg` object for the recurrent event rate, or a fitted
  `recreg` (Ghosh-Lin) object. The model type is detected automatically
  from the class.

- coxd:

  A fitted `phreg` object for the terminal event. May be `NULL` if no
  terminal event is modelled.

- n:

  Number of subjects to simulate. Default is `1`.

- data:

  The data frame on which `cox1` and `coxd` were fitted, used to draw
  covariate values for the simulated subjects. If `NULL`, covariates
  must be supplied via `r1`, `rd`, `strata1`, and `stratad`.

- type:

  Simulation type: `"default"` (auto-detected from class of `cox1`),
  `"cox-cox"`, or `"gl-cox"`.

- id:

  Name of the subject identifier variable in `data`. Default is `"id"`.

- varz:

  Variance of the frailty distribution in the two-stage model. Default
  is `1`.

- share:

  Proportion of the shared frailty assigned to the recurrent event
  process in the partial-sharing model. Default is `1`.

- cens:

  Rate of exponential censoring. Default is `0.001`.

- scale1:

  Scalar multiplier for the baseline cumulative hazard of the recurrent
  event process. Default is `1`.

- scaled:

  Scalar multiplier for the baseline cumulative hazard of the terminal
  event. Default is `1`.

- dependence:

  If non-`NULL`, falls back to
  [`sim_recurrent_list`](http://kkholst.github.io/mets/reference/sim_recurrentII.md)
  with this frailty structure (see
  [`sim_recurrentII`](http://kkholst.github.io/mets/reference/sim_recurrentII.md)
  for valid values). Default is `NULL`.

- r1:

  Optional numeric vector of length `n` of subject-specific relative
  risks for the recurrent event, used when `data = NULL`.

- rd:

  Optional numeric vector of length `n` of subject-specific relative
  risks for the terminal event, used when `data = NULL`.

- rc:

  Optional numeric vector of length `n` of subject-specific censoring
  rate multipliers.

- strata1:

  Optional integer vector of length `n` specifying the stratum index
  (0-based) for the recurrent event model, used when `data = NULL`.

- stratad:

  Optional integer vector of length `n` specifying the stratum index
  (0-based) for the terminal event model, used when `data = NULL`.

- death.code:

  Integer status code used for the terminal event in the output `status`
  column. Default is `3`.

- ...:

  Further arguments passed to
  [`sim_GLcox`](http://kkholst.github.io/mets/reference/sim_GLcox.md),
  including `nmin` and `nmax` for the linear approximation grid.

## Value

A data frame in counting-process format with one row per event interval
per subject. Column names match those in the original model formula
(entry, exit, and status variables). Additional columns include `id` and
the covariates drawn from `data` (if supplied). The terminal event is
coded as `death.code` in the status variable; recurrent events are coded
as `1`.

## Details

Covariates are drawn by bootstrap from `data` (if supplied), and
subject-specific relative risks and strata are derived from the fitted
model objects. Stratified baselines are fully supported. When
`dependence` is `NULL` (default), the simulation uses the two-stage
structure from
[`sim_GLcox`](http://kkholst.github.io/mets/reference/sim_GLcox.md);
setting `dependence` to an integer falls back to
[`sim_recurrent_list`](http://kkholst.github.io/mets/reference/sim_recurrentII.md)
with the corresponding frailty model.

## References

Scheike, T. H. (2026). Two-stage recurrent events random effects models.
*Lifetime Data Analysis*.

## See also

[`recurrent_marginal`](http://kkholst.github.io/mets/reference/recurrent_marginal.md),
[`sim_recurrent_list`](http://kkholst.github.io/mets/reference/sim_recurrentII.md),
[`sim_GLcox`](http://kkholst.github.io/mets/reference/sim_GLcox.md)

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
hf    <- hfactioncpx12
hf$x  <- as.numeric(hf$treatment)
n     <- 100

## Cox-Cox two-stage model
xr <- phreg(Surv(entry, time, status == 1) ~ x + cluster(id), data = hf)
dr <- phreg(Surv(entry, time, status == 2) ~ x + cluster(id), data = hf)
simcoxcox <- sim_recurrent_ts(xr, dr, n = n, data = hf, death.code = 2)

## Ghosh-Lin/Cox two-stage model
recGL  <- recreg(Event(entry, time, status) ~ x + cluster(id), hf, death.code = 2)
simglcox <- sim_recurrent_ts(recGL, dr, n = n, data = hf, death.code = 2)
```
