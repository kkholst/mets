# Simulation of two-stage recurrent events data based on Cox/Cox or Cox/Ghosh-Lin structure

Simulation of two-stage recurrent events data based on Cox/Cox or
Cox/Ghosh-Lin structure. type=3 will generate Cox/Cox twostage mode,
type=2 will generate Ghosh-Lin/Cox model. If the variance is var.z=0,
then generates data without any dependence or frailty. If
model="twostage" then default is to generate data from Ghosh-Lin/Cox
model, and if type=3 then will generate data with marginal Cox models
(Cox/Cox). Simulation based on linear aproximation of hazard for
two-stage models based on grid on time-scale. Must be sufficientyly
fine.

## Usage

``` r
simGLcox(
  n,
  base1,
  drcumhaz,
  var.z = 0,
  r1 = NULL,
  rd = NULL,
  rc = NULL,
  fz = NULL,
  fdz = NULL,
  model = c("twostage", "frailty", "shared"),
  type = NULL,
  share = 1,
  cens = NULL,
  nmin = 100,
  nmax = 1000
)
```

## Arguments

- n:

  number of id's

- base1:

  baseline for cox/ghosh-lin models

- drcumhaz:

  baseline for terminal event

- var.z:

  variance of gamma frailty

- r1:

  relative risk term for baseline

- rd:

  relative risk term for terminal event

- rc:

  relative risk term for censorings

- fz:

  possible transformation (function) of frailty term

- fdz:

  possible transformation (function) of frailty term for death

- model:

  twostage, frailty, shared (partly shared two-stage model)

- type:

  type of simulation, default is decided based on model

- share:

  to fit patly shared random effects model

- cens:

  censoring rate for exponential censoring

- nmin:

  default 100, at least nmin or number of rows of the two-baselines
  max(nmin,nrow(base1),nrow(drcumhaz)) points in time-grid from 0 to
  maximum time for base1

- nmax:

  default 1000, at most nmax points in time-grid

## Details

Must specify baselines of recurrent events and terminal event and
possible covariate effects.

## References

Scheike (2025), Two-stage recurrent events random effects models, LIDA,
to appear
