# Average Treatment effect for censored competing risks data using Binomial Regression

Under the standard causal assumptions we can estimate the average
treatment effect E(Y(1) - Y(0)). We need Consistency, ignorability (
Y(1), Y(0) indep A given X), and positivity.

## Usage

``` r
binregATE(
  formula,
  data,
  cause = 1,
  time = NULL,
  beta = NULL,
  treat.model = ~+1,
  cens.model = ~+1,
  offset = NULL,
  weights = NULL,
  cens.weights = NULL,
  se = TRUE,
  type = c("II", "I"),
  kaplan.meier = TRUE,
  cens.code = 0,
  no.opt = FALSE,
  method = "nr",
  augmentation = NULL,
  outcome = c("cif", "rmst", "rmtl"),
  model = c("default", "logit", "exp", "lin"),
  Ydirect = NULL,
  typeATE = "II",
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

- treat.model:

  logistic treatment model given covariates

- cens.model:

  only stratified cox model without covariates

- offset:

  offsets for partial likelihood

- weights:

  for score equations

- cens.weights:

  censoring weights

- se:

  to compute se's with IPCW adjustment, otherwise assumes that IPCW
  weights are known

- type:

  "II" adds augmentation term, and "I" classic binomial regression

- kaplan.meier:

  uses Kaplan-Meier for IPCW in contrast to exp(-Baseline)

- cens.code:

  gives censoring code

- no.opt:

  to not optimize

- method:

  for optimization

- augmentation:

  for augment binomial regression

- outcome:

  can do CIF regression "cif"=F(t\|X), "rmst"=E( min(T, t) \| X) , or E(
  I(epsilon==cause) ( t - mint(T,t)) ) \| X) depending on the number of
  the number of causes.

- model:

  exp or linear model for E( min(T, t) \| X)=exp(X^t beta), or E(
  I(epsilon==cause) ( t - mint(T,t)) ) \| X)=exp(X^t beta)

- Ydirect:

  use this outcome Y with IPCW vesion

- typeATE:

  "II" to censor augment the estimating equation

- ...:

  Additional arguments to lower level funtions (binreg that fits outcome
  model)

## Details

The first covariate in the specification of the competing risks
regression model must be the treatment variable that should be coded as
a factor. If the factor has more than two levels then it uses the mlogit
for propensity score modelling. If there are no censorings this is the
same as ordinary logistic regression modelling.

Estimates the ATE using the the standard binary double robust estimating
equations that are IPCW censoring adjusted. Rather than binomial
regression we also consider a IPCW weighted version of standard logistic
regression logitIPCWATE.

typeATE="II" will augment the estimating equation with \$\$ (A/\pi(X))
\int E( O(t) \| T \geq t, S(X))/ G_c(t,S(X)) d \hat M_c(s) \$\$ when
estimating the mean outcome for the treated.

## References

Blanche PF, Holt A, Scheike T (2022). “On logistic regression with right
censored data, with or without competing risks, and its use for
estimating treatment effects.” Lifetime data analysis, 29, 441–482.

## Author

Thomas Scheike

## Examples

``` r
library(mets); data(bmt)
dfactor(bmt)  <-  ~.

brs <- binregATE(Event(time,cause)~tcell.f+platelet+age,bmt,time=50,cause=1,
  treat.model=tcell.f~platelet+age)
summary(brs)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.199112  0.130982 -0.455831  0.057607  0.1285
#> tcell.f1    -0.637221  0.356617 -1.336177  0.061735  0.0740
#> platelet    -0.344504  0.245974 -0.826604  0.137596  0.1613
#> age          0.437222  0.107263  0.226991  0.647454  0.0000
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.81946 0.63392 1.0593
#> tcell.f1     0.52876 0.26285 1.0637
#> platelet     0.70857 0.43753 1.1475
#> age          1.54840 1.25482 1.9107
#> 
#> Average Treatment effects (G-formula) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.428748  0.027511  0.374828  0.482668  0.0000
#> treat1     0.289931  0.065905  0.160761  0.419102  0.0000
#> treat:1-0 -0.138817  0.071772 -0.279487  0.001854  0.0531
#> 
#> Average Treatment effects (double robust) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.428159  0.027612  0.374040  0.482278  0.0000
#> treat1     0.250414  0.064791  0.123426  0.377402  0.0001
#> treat:1-0 -0.177745  0.070145 -0.315226 -0.040263  0.0113
#> 
#> 
head(brs$riskDR.iid)
#>          iidriska      iidriska
#> [1,] -0.001158914 -3.544289e-05
#> [2,] -0.001200978  7.591687e-05
#> [3,] -0.001326400  3.360021e-04
#> [4,] -0.001320260  3.247937e-04
#> [5,] -0.001140662 -9.113840e-05
#> [6,] -0.001398172  4.595460e-04
head(brs$riskG.iid)
#>        riskGa.iid    riskGa.iid
#> [1,] -0.001190622 -0.0001527653
#> [2,] -0.001242318  0.0001090009
#> [3,] -0.001355147  0.0006917410
#> [4,] -0.001350560  0.0006678242
#> [5,] -0.001164390 -0.0002837983
#> [6,] -0.001403991  0.0009473264

brsi <- binregATE(Event(time,cause)~tcell.f+tcell.f*platelet+tcell.f*age,bmt,time=50,cause=1,
  treat.model=tcell.f~platelet+age)
summary(brsi)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>                    Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)       -0.161478  0.133240 -0.422623  0.099667  0.2255
#> tcell.f1          -1.029704  0.513385 -2.035920 -0.023489  0.0449
#> platelet          -0.490119  0.270795 -1.020867  0.040630  0.0703
#> age                0.445946  0.112205  0.226028  0.665864  0.0001
#> tcell.f1:platelet  0.956586  0.694975 -0.405540  2.318713  0.1687
#> tcell.f1:age      -0.154714  0.427213 -0.992036  0.682609  0.7172
#> 
#> exp(coeffients):
#>                   Estimate    2.5%   97.5%
#> (Intercept)        0.85089 0.65533  1.1048
#> tcell.f1           0.35711 0.13056  0.9768
#> platelet           0.61255 0.36028  1.0415
#> age                1.56197 1.25361  1.9462
#> tcell.f1:platelet  2.60280 0.66662 10.1626
#> tcell.f1:age       0.85666 0.37082  1.9790
#> 
#> Average Treatment effects (G-formula) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.427718  0.027557  0.373706  0.481729  0.0000
#> treat1     0.265587  0.069738  0.128903  0.402272  0.0001
#> treat:1-0 -0.162130  0.074795 -0.308726 -0.015535  0.0302
#> 
#> Average Treatment effects (double robust) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.428133  0.027614  0.374012  0.482255  0.0000
#> treat1     0.253093  0.066785  0.122196  0.383990  0.0002
#> treat:1-0 -0.175040  0.072033 -0.316223 -0.033857  0.0151
#> 
#> 
head(brs$riskDR.iid)
#>          iidriska      iidriska
#> [1,] -0.001158914 -3.544289e-05
#> [2,] -0.001200978  7.591687e-05
#> [3,] -0.001326400  3.360021e-04
#> [4,] -0.001320260  3.247937e-04
#> [5,] -0.001140662 -9.113840e-05
#> [6,] -0.001398172  4.595460e-04
head(brs$riskG.iid)
#>        riskGa.iid    riskGa.iid
#> [1,] -0.001190622 -0.0001527653
#> [2,] -0.001242318  0.0001090009
#> [3,] -0.001355147  0.0006917410
#> [4,] -0.001350560  0.0006678242
#> [5,] -0.001164390 -0.0002837983
#> [6,] -0.001403991  0.0009473264

```
