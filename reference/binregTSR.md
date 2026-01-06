# 2 Stage Randomization for Survival Data or competing Risks Data

Under two-stage randomization we can estimate the average treatment
effect E(Y(i,j)) of treatment regime (i,j). The estimator can be
agumented in different ways: using the two randomizations and the
dynamic censoring augmetatation. The treatment's must be given as
factors.

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

  formula with outcome (see `coxph`)

- data:

  data frame

- cause:

  cause of interest

- time:

  time of interest

- cens.code:

  gives censoring code

- response.code:

  code of status of survival data that indicates a response at which 2nd
  randomization is performed

- augmentR0:

  augmentation model for 1st randomization

- treat.model0:

  logistic treatment model for 1st randomization

- augmentR1:

  augmentation model for 2nd randomization

- treat.model1:

  logistic treatment model for 2ndrandomization

- augmentC:

  augmentation model for censoring

- cens.model:

  stratification for censoring model based on observed covariates

- estpr:

  estimate randomization probabilities using model

- response.name:

  can give name of response variable, otherwise reads this as first
  variable of treat.model1

- offset:

  not implemented

- weights:

  not implemented

- cens.weights:

  can be given

- beta:

  starting values

- kaplan.meier:

  for censoring weights, rather than exp cumulative hazard

- no.opt:

  not implemented

- method:

  not implemented

- augmentation:

  not implemented

- outcome:

  can be c("cif","rmst","rmst-cause")

- model:

  not implemented, uses linear regression for augmentation

- Ydirect:

  use this Y instead of outcome constructed inside the program (e.g.
  I(T\< t, epsilon=1)), see binreg for more on this

- return.dataw:

  to return weighted data for all treatment regimes

- pi0:

  set up known randomization probabilities

- pi1:

  set up known randomization probabilities

- cens.time.fixed:

  to use time-dependent weights for censoring estimation using weights

- outcome.iid:

  to get iid contribution from outcome model (here linear regression
  working models).

- meanCs:

  \(0\) indicates that censoring augmentation is centered by
  CensAugment.times/n

- ...:

  Additional arguments to lower level funtions

## Details

The solved estimating eqution is \$\$ ( I(min(T_i,t) \< G_i)/G_c(min(T_i
,t)) I(T \leq t, \epsilon=1 ) - AUG_0 - AUG_1 + AUG_C - p(i,j)) = 0 \$\$
where using the covariates from augmentR0 \$\$ AUG_0 = \frac{A_0(i) -
\pi_0(i)}{ \pi_0(i)} X_0 \gamma_0\$\$ and using the covariates from
augmentR1 \$\$ AUG_1 = \frac{A_0(i)}{\pi_0(i)} \frac{A_1(j) - \pi_1(j)}{
\pi_1(j)} X_1 \gamma_1\$\$ and the censoring augmentation is \$\$ AUG_C
= \int_0^t \gamma_c(s)^T (e(s) - \bar e(s)) \frac{1}{G_c(s) } dM_c(s)
\$\$ where \$\$ \gamma_c(s)\$\$ is chosen to minimize the variance given
the dynamic covariates specified by augmentC.

In the observational case, we can use propensity score modelling and
outcome modelling (using linear regression).

Standard errors are estimated using the influence function of all
estimators and tests of differences can therefore be computed
subsequently.

## Author

Thomas Scheike

## Examples

``` r
library(mets)
ddf <- mets:::gsim(200,covs=1,null=0,cens=1,ce=2)

bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),ddf$datat,time=2,cause=c(1),
        cens.code=0,treat.model0=A0.f~+1,treat.model1=A1.f~A0.f,
        augmentR1=~X11+X12+TR,augmentR0=~X01+X02,
        augmentC=~A01+A02+X01+X02+A11t+A12t+X11+X12+TR,
        response.code=2)
summary(bb) 
#> Simple estimator :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.6513604 0.05894890
#> A0.f=1, response*A1.f=2 0.8454786 0.06940480
#> A0.f=2, response*A1.f=1 0.3381049 0.07717948
#> A0.f=2, response*A1.f=2 0.1536713 0.05861389
#> 
#> First Randomization Augmentation :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.6509532 0.05841946
#> A0.f=1, response*A1.f=2 0.8466905 0.06961735
#> A0.f=2, response*A1.f=1 0.3372760 0.07745806
#> A0.f=2, response*A1.f=2 0.1539456 0.05837829
#> 
#> Second Randomization Augmentation :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.6575295 0.05738919
#> A0.f=1, response*A1.f=2 0.8256192 0.06904964
#> A0.f=2, response*A1.f=1 0.3298011 0.07878161
#> A0.f=2, response*A1.f=2 0.1492153 0.05848652
#> 
#> 1st and 2nd Randomization Augmentation :
#>                              coef           
#> A0.f=1, response*A1.f=1 0.6614013 0.05563921
#> A0.f=1, response*A1.f=2 0.8295814 0.06803427
#> A0.f=2, response*A1.f=1 0.3304605 0.07864288
#> A0.f=2, response*A1.f=2 0.1502968 0.05795214
#> 
```
