# Two-Stage Randomization for Survival or Competing Risks Data

Estimates the average treatment effect \\E(Y(i,j))\\ of treatment regime
\\(i,j)\\ under two-stage randomization. The estimator can be augmented
using information from both randomizations and dynamic censoring
augmentation to improve efficiency.

## Usage

``` r
binregTSR(
  formula,
  data,
  cause = 1,
  time = NULL,
  cens.code = 0,
  response.code = NULL,
  augmentR0 = NULL,
  treat.model0 = ~+1,
  augmentR1 = NULL,
  treat.model1 = ~+1,
  augmentC = NULL,
  cens.model = ~+1,
  estpr = c(1, 1),
  response.name = NULL,
  offset = NULL,
  weights = NULL,
  cens.weights = NULL,
  beta = NULL,
  kaplan.meier = TRUE,
  no.opt = FALSE,
  method = "nr",
  augmentation = NULL,
  outcome = c("cif", "rmst", "rmst-cause"),
  model = "exp",
  Ydirect = NULL,
  return.dataw = 0,
  pi0 = 0.5,
  pi1 = 0.5,
  cens.time.fixed = 1,
  outcome.iid = 1,
  meanCs = 0,
  ...
)
```

## Arguments

- formula:

  Formula with outcome (see `coxph`), typically
  `Event(entry,time,status)~+1+cluster(id)`.

- data:

  Data frame containing all variables.

- cause:

  Cause of interest for competing risks (default 1).

- time:

  Time point for estimation.

- cens.code:

  Censoring code (default 0).

- response.code:

  Code of status indicating response at which 2nd randomization occurs.

- augmentR0:

  Covariates for augmentation model of the first randomization.

- treat.model0:

  Logistic treatment model for the first randomization.

- augmentR1:

  Covariates for augmentation model of the second randomization.

- treat.model1:

  Logistic treatment model for the second randomization.

- augmentC:

  Covariates for censoring augmentation model.

- cens.model:

  Stratification for censoring model based on observed covariates.

- estpr:

  Logical; estimate randomization probabilities using model (default
  TRUE).

- response.name:

  Name of response variable (reads from `treat.model1` if NULL).

- offset:

  Not implemented.

- weights:

  Not implemented.

- cens.weights:

  Can be provided externally.

- beta:

  Starting values for optimization.

- kaplan.meier:

  Logical; use Kaplan-Meier for censoring weights rather than exp
  cumulative hazard.

- no.opt:

  Not implemented.

- method:

  Not implemented.

- augmentation:

  Not implemented.

- outcome:

  Outcome type: `"cif"` (cumulative incidence), `"rmst"` (restricted
  mean survival time), or `"rmst-cause"` (restricted mean time lost for
  cause).

- model:

  Not implemented, uses linear regression for augmentation.

- Ydirect:

  Use this Y instead of outcome constructed inside the program.

- return.dataw:

  Logical; return weighted data for all treatment regimes.

- pi0:

  Known randomization probabilities for first randomization.

- pi1:

  Known randomization probabilities for second randomization.

- cens.time.fixed:

  Logical; use time-dependent weights for censoring estimation.

- outcome.iid:

  Logical; get iid contribution from outcome model.

- meanCs:

  Logical; indicates censoring augmentation is centered by
  `CensAugment.times/n`.

- ...:

  Additional arguments to lower-level functions.

## Value

An object of class `"binregTSR"` containing:

- riskG:

  Simple estimator results (coefficient, SE).

- riskG0:

  First randomization augmentation results.

- riskG1:

  Second randomization augmentation results.

- riskG01:

  Both randomizations augmentation results.

- riskG.iid:

  Influence functions for all estimators.

- varG:

  Variance-covariance matrices.

- MGc:

  Censoring martingale contributions.

- CensAugment.times:

  Censoring augmentation terms.

- dynCens.coef:

  Dynamic censoring coefficients.

- dataW:

  Weighted data (if `return.dataw=TRUE`).

## Details

The method solves the estimating equation: \$\$ \frac{I(\min(T_i,t) \<
G_i)}{G_c(\min(T_i,t))} I(T \leq t, \epsilon=1) - AUG_0 - AUG_1 +
AUG_C - p(i,j) = 0 \$\$ where:

- \\AUG_0 = \frac{A_0(i) - \pi_0(i)}{\pi_0(i)} X_0 \gamma_0\\ uses
  covariates from `augmentR0`

- \\AUG_1 = \frac{A_0(i)}{\pi_0(i)} \frac{A_1(j) - \pi_1(j)}{\pi_1(j)}
  X_1 \gamma_1\\ uses covariates from `augmentR1`

- \\AUG_C = \int_0^t \gamma_c(s)^T (e(s) - \bar e(s)) \frac{1}{G_c(s)}
  dM_c(s)\\ is the censoring augmentation

Standard errors are estimated using the influence function of all
estimators, enabling tests of differences to be computed subsequently.
The method handles both survival data and competing risks data, and
supports multiple treatment levels.

## References

Scheike, T. H. (2024). Two-stage randomization analysis for survival
data. mets package documentation.

## See also

[`binreg`](http://kkholst.github.io/mets/reference/binreg.md),
[`phreg_rct`](http://kkholst.github.io/mets/reference/phreg_rct.md)

## Author

Thomas Scheike

## Examples

``` r
ddf <- mets:::gsim(200,covs=1,null=0,cens=1,ce=2)

bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),ddf$datat,time=2,cause=c(1),
        cens.code=0,treat.model0=A0.f~+1,treat.model1=A1.f~A0.f,
        augmentR1=~X11+X12+TR,augmentR0=~X01+X02,
        augmentC=~A01+A02+X01+X02+A11t+A12t+X11+X12+TR,
        response.code=2)
summary(bb) 
#> Simple estimator :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.7965990 0.06876063
#> A0.f=1, response*A1.f=2 0.7290403 0.05904496
#> A0.f=2, response*A1.f=1 0.2652459 0.08526775
#> A0.f=2, response*A1.f=2 0.3463623 0.07507682
#> 
#> First Randomization Augmentation :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.8020553 0.07037422
#> A0.f=1, response*A1.f=2 0.7380517 0.05897235
#> A0.f=2, response*A1.f=1 0.2643500 0.08356227
#> A0.f=2, response*A1.f=2 0.3392007 0.07533492
#> 
#> Second Randomization Augmentation :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.8270658 0.05679816
#> A0.f=1, response*A1.f=2 0.7414130 0.06070947
#> A0.f=2, response*A1.f=1 0.2497165 0.09169221
#> A0.f=2, response*A1.f=2 0.3645886 0.07180896
#> 
#> 1st and 2nd Randomization Augmentation :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.8352750 0.05697325
#> A0.f=1, response*A1.f=2 0.7531912 0.05912032
#> A0.f=2, response*A1.f=1 0.2495341 0.08954463
#> A0.f=2, response*A1.f=2 0.3630083 0.07086160
#> 
```
