# Percentage of years lost due to cause regression

Estimates the percentage of the years lost that is due to a cause and
how covariates affects this percentage by doing ICPW regression.

## Usage

``` r
binregRatio(
  formula,
  data,
  cause = 1,
  time = NULL,
  beta = NULL,
  type = c("II", "I"),
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
  outcome = c("cif", "rmtl"),
  model = c("logit", "exp", "lin"),
  Ydirect = NULL,
  ...
)
```

## Arguments

- formula:

  formula with outcome (see `coxph`)

- data:

  data frame

- cause:

  cause of interest (numeric variable)

- time:

  time of interest

- beta:

  starting values

- type:

  "II" adds augmentation term, and "I" classical outcome IPCW regression

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

- outcome:

  can do CIF regression "cif"=F(t\|X), "rmtl"=E( t- min(T, t) \| X)"

- model:

  logit, exp or lin(ear)

- Ydirect:

  use this Y instead of outcome constructed inside the program, should
  be a matrix with two column for numerator and denominator.

- ...:

  Additional arguments to lower level funtions

## Details

Let the years lost be \$\$Y1= t- min(T ,) \$\$ and the years lost due to
cause 1 \$\$Y2= I(epsilon==1) ( t- min(T ,t) \$\$ , then we model the
ratio \$\$logit( E(Y2 \| X)/E(Y1 \| X)) = X^T \beta \$\$. Estimation is
based on on binomial regresion IPCW response estimating equation: \$\$ X
( \Delta^{ipcw}(t) Y2 expit(X^T \beta) - Y1 ) = 0 \$\$ where
\$\$\Delta^{ipcw}(t) = I((min(t,T)\< C)/G_c(min(t,T)-)\$\$ is IPCW
adjustment of the response \$\$Y(t)= I(T \leq t, \epsilon=1 )\$\$.

(type="I") sovlves this estimating equation using a stratified
Kaplan-Meier for the censoring distribution. For (type="II") the default
an additional censoring augmentation term \$\$X \int E(Y(t)\|
T\>s)/G_c(s) d \hat M_c\$\$ is added.

The variance is based on the squared influence functions that are also
returned as the iid component. naive.var is variance under known
censoring model.

Censoring model may depend on strata (cens.model=~strata(gX)).

## References

Scheike & Tanaka (2025), Restricted mean time lost ratio regression:
Percentage of restricted mean time lost due to specific cause, WIP

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt); bmt$time <- bmt$time+runif(408)*0.001

rmst30 <- rmstIPCW(Event(time,cause!=0)~platelet+tcell+age,bmt,time=30,cause=1)
rmst301 <- rmstIPCW(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=1)
rmst302 <- rmstIPCW(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=2)

estimate(rmst30)
#>             Estimate Std.Err     2.5%    97.5%   P-value
#> (Intercept)   2.6104 0.05744  2.49778  2.72294 0.000e+00
#> platelet      0.2450 0.08435  0.07970  0.41035 3.674e-03
#> tcell         0.1456 0.11877 -0.08719  0.37838 2.202e-01
#> age          -0.1728 0.03794 -0.24718 -0.09844 5.255e-06
estimate(rmst301)
#>             Estimate Std.Err    2.5%      97.5%    P-value
#> (Intercept)   2.4449 0.07291  2.3020  2.5878403 1.655e-246
#> platelet     -0.4519 0.16738 -0.7799 -0.1238183  6.940e-03
#> tcell        -0.4935 0.25140 -0.9862 -0.0007153  4.967e-02
#> age           0.2746 0.06538  0.1465  0.4027874  2.658e-05
estimate(rmst302)
#>             Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept)  1.45850  0.1362  1.1915 1.7255 9.617e-27
#> platelet    -0.01886  0.2224 -0.4548 0.4171 9.324e-01
#> tcell        0.40810  0.2648 -0.1108 0.9270 1.232e-01
#> age          0.04566  0.1191 -0.1878 0.2791 7.014e-01

## percentage of total cumulative incidence due to cause 1
rmtlratioI <- rmtlRatio(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=1)
summary(rmtlratioI)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  0.977101  0.178189  0.627856  1.326345  0.0000
#> platelet    -0.412142  0.323288 -1.045775  0.221490  0.2024
#> tcell       -0.855816  0.419188 -1.677408 -0.034223  0.0412
#> age          0.217651  0.168416 -0.112439  0.547741  0.1962
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  2.65674 1.87359 3.7672
#> platelet     0.66223 0.35142 1.2479
#> tcell        0.42494 0.18686 0.9664
#> age          1.24315 0.89365 1.7293
#> 
#> 

pp <- predict(rmtlratioI,bmt[1:5,])
pp
#>        pred         se     lower     upper
#> 1 0.7349115 0.03513428 0.6660483 0.8037747
#> 2 0.7529172 0.03812709 0.6781881 0.8276463
#> 3 0.7839582 0.05014699 0.6856701 0.8822463
#> 4 0.7828518 0.04965019 0.6855374 0.8801661
#> 5 0.7243286 0.03566672 0.6544218 0.7942353

newdata <- data.frame(platelet=1,tcell=1,age=1)
## percentage of total cumulative incidence due to cause 1
cifratio <- binregRatio(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=1)
summary(cifratio)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)  0.91183  0.17554  0.56778  1.25589  0.0000
#> platelet    -0.46563  0.31922 -1.09128  0.16003  0.1447
#> tcell       -1.06111  0.41568 -1.87583 -0.24639  0.0107
#> age          0.21200  0.16115 -0.10384  0.52784  0.1883
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  2.48888 1.76434 3.5109
#> platelet     0.62774 0.33579 1.1735
#> tcell        0.34607 0.15323 0.7816
#> age          1.23615 0.90137 1.6953
#> 
#> 
pp <- predict(cifratio,newdata)
pp
#>        pred       se     lower    upper
#> 1 0.4006137 0.106969 0.1909544 0.610273

rmtlratioI <- binregRatio(Event(time,cause)~platelet+tcell+age,bmt,
                               time=30,cause=1,outcome="rmtl")
summary(rmtlratioI)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  0.977101  0.178189  0.627856  1.326345  0.0000
#> platelet    -0.412142  0.323288 -1.045775  0.221490  0.2024
#> tcell       -0.855816  0.419188 -1.677408 -0.034223  0.0412
#> age          0.217651  0.168416 -0.112439  0.547741  0.1962
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  2.65674 1.87359 3.7672
#> platelet     0.66223 0.35142 1.2479
#> tcell        0.42494 0.18686 0.9664
#> age          1.24315 0.89365 1.7293
#> 
#> 

pp <- predict(rmtlratioI,newdata)
pp
#>        pred        se    lower     upper
#> 1 0.4817067 0.1092932 0.267492 0.6959214
```
