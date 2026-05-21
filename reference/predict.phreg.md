# Predictions from Proportional Hazards Model

Computes predictions for survival probability, cumulative hazard, or
risk at specified time points for new data or existing data. Includes
standard errors and confidence intervals.

## Usage

``` r
# S3 method for class 'phreg'
predict(
  object,
  newdata,
  times = NULL,
  individual.time = FALSE,
  tminus = FALSE,
  se = TRUE,
  robust = FALSE,
  conf.type = "log",
  conf.int = 0.95,
  km = FALSE,
  ...
)
```

## Arguments

- object:

  Object of class `"phreg"`.

- newdata:

  Data frame for prediction. If `NULL`, predictions are made for the
  original data.

- times:

  Time points for prediction. Defaults to all unique event times in the
  model.

- individual.time:

  Logical; if `TRUE`, uses one time per subject (requires `newdata` and
  `times` to be same length).

- tminus:

  Logical; if `TRUE`, predicts at \\t-\\ (strictly before \\t\\).

- se:

  Logical; if `TRUE`, computes standard errors and confidence intervals.

- robust:

  Logical; if `TRUE`, uses robust standard errors (default for most
  functions).

- conf.type:

  Transformation for survival estimates: `"log"` (default) or `"plain"`.

- conf.int:

  Confidence level (default 0.95).

- km:

  Logical; if `TRUE`, uses Kaplan-Meier product-limit for baseline;
  otherwise uses exponential of cumulative baseline.

- ...:

  Additional arguments for plotting functions.

## Value

An object of class `"predictphreg"` containing:

- surr:

  Matrix of survival probabilities.

- cumhaz:

  Matrix of cumulative hazards.

- cif:

  Matrix of cumulative incidence functions (if applicable).

- times:

  Vector of time points.

- surv.upper, surv.lower:

  Confidence bounds for survival.

- RR:

  Relative risks.

## Author

Thomas Scheike
