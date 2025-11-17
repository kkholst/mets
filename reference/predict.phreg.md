# Predictions from proportional hazards model

Predictions from proportional hazards model

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

  phreg object

- newdata:

  data.frame

- times:

  Time where to predict variable, default is all time-points from the
  object sorted

- individual.time:

  to use one (individual) time per subject, and then newdata and times
  have same length and makes only predictions for these individual
  times.

- tminus:

  to make predictions in T- that is strictly before given times, useful
  for IPCW techniques

- se:

  with standard errors and upper and lower confidence intervals.

- robust:

  to get robust se's also default for most functions (uses robse.cumhaz
  otherwise se.cumhaz).

- conf.type:

  transformation for suvival estimates, default is log

- conf.int:

  significance level

- km:

  to use Kaplan-Meier product-limit for baseline \$\$S\_{s0}(t)= (1 -
  dA\_{s0}(t))\$\$, otherwise take exp of cumulative baseline.

- ...:

  Additional arguments to plot functions
