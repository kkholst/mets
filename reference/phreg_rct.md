# Lu-Tsiatis More Efficient Log-Rank for Randomized studies with baseline covariates

Efficient implementation of the Lu-Tsiatis improvement using baseline
covariates, extended to competing risks and recurrent events. Results
almost equivalent with the speffSurv function of the speff2trial
function in the survival case. A dynamic censoring augmentation
regression is also computed to gain even more from the censoring
augmentation. Furhter, we also deal with twostage randomizations. The
function was implemented to deal with recurrent events (start,stop) +
cluster, and more examples in vignette.

## Usage

``` r
phreg_rct(
  formula,
  data,
  cause = 1,
  cens.code = 0,
  typesR = c("R0", "R1", "R01"),
  typesC = c("C", "dynC"),
  weights = NULL,
  augmentR0 = NULL,
  augmentR1 = NULL,
  augmentC = NULL,
  treat.model = ~+1,
  RCT = TRUE,
  treat.var = NULL,
  km = TRUE,
  level = 0.95,
  cens.model = NULL,
  estpr = 1,
  pi0 = 0.5,
  base.augment = FALSE,
  return.augmentR0 = FALSE,
  mlogit = FALSE,
  ...
)
```

## Arguments

- formula:

  formula with 'Surv' or 'Event' outcome (see `coxph`) and treatment
  (randomization 0/1)

- data:

  data frame

- cause:

  to use for competing risks, recurrent events data

- cens.code:

  to use for competing risks, recurrent events data

- typesR:

  augmentations used for randomization

- typesC:

  augmentations used for censoring

- weights:

  weights for score equation

- augmentR0:

  formula for the randomization augmentation (~age+sex)

- augmentR1:

  formula for the randomization augmentation (~age+sex)

- augmentC:

  formula for the censoring augmentation (~age+sex)

- treat.model:

  propensity score model, default is ~+1, assuming an RCT study

- RCT:

  if false will use propensity score adjustment for marginal model

- treat.var:

  in case of twostage randomization, this variable is 1 for the
  treatment times, if start,stop then default assumes that only one
  treatment at first record

- km:

  use Kaplan-Meier for the censoring weights (stratified on treatment)

- level:

  of confidence intervals

- cens.model:

  default is censoring model ~strata(treatment) but any model can be
  used to make censoring martingales

- estpr:

  estimates propensity scores

- pi0:

  possible fixed propensity scores for randomizations

- base.augment:

  TRUE to covariate augment baselines (only for R0 augmentation)

- return.augmentR0:

  to return augmentation data

- mlogit:

  if TRUE then forces use of this function for propensity scores,
  default for binary treatment is glm

- ...:

  Additional arguments to phreg function

## References

Lu, Tsiatis (2008), Improving the efficiency of the log-rank test using
auxiliary covariates, Biometrika, 679–694

Scheike, Nerstroem and Martinussen (2025), Randomized clinical trials
and the proportional hazards model for recurrent events.

## Author

Thomas Scheike

## Examples

``` r
## Lu, Tsiatis simulation
data <- mets:::simLT(0.7,100)
dfactor(data) <- Z.f~Z

out <- phreg_rct(Surv(time,status)~Z.f,data=data,augmentR0=~X,augmentC=~factor(Z):X)
summary(out)
#>                 Estimate   Std.Err       2.5%     97.5%   P-value
#> Marginal-Z.f1  0.1781155 0.2639114 -0.3391412 0.6953723 0.4997350
#> R0_C:Z.f1     -0.0175242 0.2071406 -0.4235122 0.3884638 0.9325790
#> R0_dynC:Z.f1   0.1109733 0.2022974 -0.2855223 0.5074689 0.5833039
#> attr(,"class")
#> [1] "summary.phreg_rct"
```
