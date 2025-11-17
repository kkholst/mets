# While-Alive estimands for recurrent events

Considers the ratio of means \$\$E(N(min(D,t)))/E(min(D,t))\$\$ and the
the mean of the events per time unit \$\$E(N(min(D,t))/min(D,t))\$\$
both based on IPCW etimation. RMST estimator equivalent to Kaplan-Meier
based estimator.

## Usage

``` r
WA_recurrent(
  formula,
  data,
  time = NULL,
  cens.code = 0,
  cause = 1,
  death.code = 2,
  trans = NULL,
  cens.formula = NULL,
  augmentR = NULL,
  augmentC = NULL,
  type = NULL,
  marks = NULL,
  ...
)
```

## Arguments

- formula:

  Event formula first covariate on rhs must be a factor giving the
  treatment

- data:

  data frame

- time:

  for estimation

- cens.code:

  of censorings

- cause:

  of events

- death.code:

  of terminal events

- trans:

  possible power for mean of events per time-unit

- cens.formula:

  censoring model, default is to use strata(treatment)

- augmentR:

  covariates for model of mean ratio

- augmentC:

  covariates for censoring augmentation

- type:

  augmentation for call of binreg, when augmentC is given default is "I"
  and otherwise "II"

- marks:

  possible marks for composite outcome situation for model for counts
  with marks

- ...:

  arguments for binregATE

## Author

Thomas Scheike
