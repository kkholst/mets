# Simulation of Two-Stage Recurrent Events Data

Simulates data based on Cox/Cox or Cox/Ghosh-Lin structures for
recurrent events with a terminal event.

## Usage

``` r
sim_GLcox(
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

  Number of IDs.

- base1:

  Baseline cumulative hazard for recurrent events.

- drcumhaz:

  Baseline cumulative hazard for the terminal event.

- var.z:

  Variance of the gamma frailty.

- r1:

  Relative risk term for the recurrent event baseline.

- rd:

  Relative risk term for the terminal event.

- rc:

  Relative risk term for censoring.

- fz:

  Function for transformation of the frailty term.

- fdz:

  Function for transformation of the frailty term for death.

- model:

  Model type: `"twostage"`, `"frailty"`, or `"shared"`.

- type:

  Type of simulation (default depends on `model`).

- share:

  Proportion of shared random effects (for partially shared models).

- cens:

  Censoring rate for exponential censoring.

- nmin:

  Minimum number of points in the time grid (default 100).

- nmax:

  Maximum number of points in the time grid (default 1000).

## Value

A data frame with simulated recurrent events and terminal events,
including frailty terms.

## Details

- `type=3`: Generates data from a Cox/Cox two-stage model.

- `type=2`: Generates data from a Ghosh-Lin/Cox model.

If `var.z=0`, data is generated without dependence or frailty. If
`model="twostage"`, the default is to generate data from a Ghosh-Lin/Cox
model. If `type=3`, data is generated with marginal Cox models
(Cox/Cox).

Simulation is based on a linear approximation of the hazard for
two-stage models on a time grid. The grid must be sufficiently fine.

## References

Scheike (2025), Two-stage recurrent events random effects models, LIDA,
to appear.

## Author

Thomas Scheike
