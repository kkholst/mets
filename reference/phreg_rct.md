# Lu-Tsiatis More Efficient Log-Rank for Randomized Studies with Baseline Covariates

Efficient implementation of the Lu-Tsiatis improvement using baseline
covariates, extended to competing risks and recurrent events. The
results are almost equivalent to the `speffSurv` function of the
`speff2trial` package in the survival case. A dynamic censoring
augmentation regression is also computed to gain additional efficiency
from the censoring augmentation. The function handles two-stage
randomizations and recurrent events (start,stop) with cluster structure.

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

  Formula with `Surv` or `Event` outcome (see `coxph`) and treatment
  variable (randomization 0/1). The treatment variable must be the first
  covariate on the right-hand side.

- data:

  Data frame containing all variables referenced in the formula.

- cause:

  Numeric code for the event of interest in competing risks or recurrent
  events.

- cens.code:

  Numeric code for censoring in competing risks or recurrent events.

- typesR:

  Character vector specifying augmentation types for randomization
  (options: "R0" for baseline, "R1" for post-baseline, "R01" for both).

- typesC:

  Character vector specifying augmentation types for censoring (options:
  "C" for static, "dynC" for dynamic).

- weights:

  Weights to be used for `phreg`.

- augmentR0:

  Formula for the first randomization augmentation (e.g., `~age+sex`).

- augmentR1:

  Formula for the second randomization augmentation (e.g., `~age+sex`).

- augmentC:

  Formula for the censoring augmentation (e.g., `~age+sex`).

- treat.model:

  Propensity score model formula (default is `~+1`, assuming RCT).

- RCT:

  Logical; if FALSE, uses propensity score adjustment for marginal
  model.

- treat.var:

  Variable indicating treatment times in two-stage randomization.

- km:

  Logical; use Kaplan-Meier for censoring weights (stratified on
  treatment).

- level:

  Confidence level for intervals (default 0.95).

- cens.model:

  Censoring model formula (default is `~strata(treatment)`).

- estpr:

  Numeric code (1/0); estimate propensity scores or not (default TRUE).

- pi0:

  Fixed propensity scores for randomizations (if not estimating).

- base.augment:

  Logical; covariate augment baselines (only for R0 augmentation).

- return.augmentR0:

  Logical; return augmentation data.

- mlogit:

  Logical; use multinomial logistic regression for propensity scores.

- ...:

  Additional arguments passed to `phreg` function.

## Value

An object of class `"phreg_rct"` containing:

- coefs:

  Coefficient estimates for the treatment effect.

- iid:

  Influence function (IID) decomposition for variance estimation.

- AugR0, AugR1, AugCdyn, AugClt:

  Augmentation terms for different strategies.

- cumhaz:

  Cumulative hazards.

- var:

  Variance-covariance matrix.

- se:

  Standard errors.

- call:

  Original function call.

- formula:

  Formula used.

- data:

  The data used (if requested).

The object includes results for different augmentation combinations (R0,
R1, R01, C, dynC).

## Details

The method improves the efficiency of the log-rank test by utilizing
auxiliary baseline covariates to reduce variance, particularly useful in
randomized clinical trials (RCTs) where covariate adjustment can
increase power.

## References

Lu, T. and Tsiatis, A. A. (2008), Improving the efficiency of the
log-rank test using auxiliary covariates, Biometrika, 95, 679–694.

Scheike, T. H., Nerstroem, C. and Martinussen, T. (2026), Randomized
clinical trials and the proportional hazards model for recurrent events,
TEST.

## Author

Thomas Scheike

## Examples

``` r
## Lu, Tsiatis simulation
data <- mets:::simLT(0.7,100)
dfactor(data) <- Z.f~Z

out <- phreg_rct(Surv(time,status)~Z.f,data=data,augmentR0=~X,augmentC=~factor(Z):X)
summary(out)
#>                Estimate   Std.Err        2.5%     97.5%    P-value
#> Marginal-Z.f1 0.6721162 0.2763128  0.13055314 1.2136793 0.01499718
#> R0_C:Z.f1     0.4990308 0.2230568  0.06184748 0.9362141 0.02527090
#> R0_dynC:Z.f1  0.3950996 0.2172938 -0.03078852 0.8209876 0.06902238
#> attr(,"class")
#> [1] "summary.phreg_rct"
```
