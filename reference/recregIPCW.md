# IPCW Estimator for Recurrent Events

Computes the Inverse Probability of Censoring Weighted (IPCW) estimator
for the mean number of recurrent events. Supports various estimators
including the Ghosh-Lin and Lawless-Cook estimators.

## Usage

``` r
recregIPCW(
  formula,
  data = data,
  cause = 1,
  cens.code = 0,
  death.code = 2,
  cens.model = ~1,
  km = TRUE,
  times = NULL,
  beta = NULL,
  offset = NULL,
  type = c("II", "I"),
  marks = NULL,
  weights = NULL,
  model = c("exp", "lin"),
  no.opt = FALSE,
  augmentation = NULL,
  method = "nr",
  se = TRUE,
  ...
)
```

## Arguments

- formula:

  Formula with an 'Event' outcome.

- data:

  Data frame.

- cause:

  Cause of interest (default is 1).

- cens.code:

  Censoring code (default is 0).

- death.code:

  Death code (default is 2).

- cens.model:

  Formula for the censoring model (default is `~1`).

- km:

  Logical; if `TRUE`, uses Kaplan-Meier for censoring weights; otherwise
  uses Cox model.

- times:

  Time points for estimation (required).

- beta:

  Initial values for coefficients (optional).

- offset:

  Offsets.

- type:

  Type of estimator: `"II"` (default) or `"I"`.

- marks:

  Mark values.

- weights:

  Weights.

- model:

  Model type for the mean: `"exp"` (default) or `"lin"`.

- no.opt:

  Logical; if `TRUE`, skips optimization.

- augmentation:

  Augmentation terms.

- method:

  Optimization method (default is "nr").

- se:

  Logical; if `TRUE`, computes standard errors.

- ...:

  Additional arguments.

## Value

An object of class `"binreg"` (extending `"resmean"`) containing:

- coef:

  Estimated coefficients.

- var:

  Variance-covariance matrix.

- iid:

  Influence functions.

- times:

  Time points.

- Y:

  Observed counts.

## See also

[`recreg`](http://kkholst.github.io/mets/reference/recreg.md)

## Author

Thomas Scheike
