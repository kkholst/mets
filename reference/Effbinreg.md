# Efficient IPCW for binary data

Simple version of comp.risk function of timereg for just one time-point
thus fitting the model \$\$E(T \leq t \| X ) = expit( X^T beta) \$\$

## Usage

``` r
Effbinreg(
  formula,
  data,
  cause = 1,
  time = NULL,
  beta = NULL,
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
  h = NULL,
  MCaugment = NULL,
  ...
)
```

## Arguments

- formula:

  formula with outcome (see `coxph`)

- data:

  data frame

- cause:

  cause of interest

- time:

  time of interest

- beta:

  starting values

- offset:

  offsets for partial likelihood

- weights:

  for score equations

- cens.weights:

  censoring weights

- cens.model:

  only stratified cox model without covariates

- se:

  to compute se's based on IPCW

- kaplan.meier:

  uses Kaplan-Meier for IPCW in contrast to exp(-Baseline)

- cens.code:

  gives censoring code

- no.opt:

  to not optimize

- method:

  for optimization

- augmentation:

  to augment binomial regression

- h:

  h for estimating equation

- MCaugment:

  iid of h and censoring model

- ...:

  Additional arguments to lower level funtions

- model:

  exp or linear

## Details

Based on binomial regresion IPCW response estimating equation: \$\$ X (
\Delta (T \leq t)/G_c(T_i-) - expit( X^T beta)) = 0 \$\$ for IPCW
adjusted responses.

Based on binomial regresion IPCW response estimating equation: \$\$ h(X)
X ( \Delta (T \leq t)/G_c(T_i-) - expit( X^T beta)) = 0 \$\$ for IPCW
adjusted responses where \$h\$ is given as an argument together with iid
of censoring with h. By using appropriately the h argument we can also
do the efficient IPCW estimator estimator this works the prepsurv and
prepcif for survival or competing risks data. In this case also the
censoring martingale should be given for variance calculation and this
also comes out of the prepsurv or prepcif functions. (Experimental
version at this stage).

Variance is based on \$\$ \sum w_i^2 \$\$ also with IPCW adjustment, and
naive.var is variance under known censoring model.

Censoring model may depend on strata.

## Author

Thomas Scheike
