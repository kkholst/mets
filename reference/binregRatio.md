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
#>             Estimate Std.Err     2.5%   97.5%   P-value
#> (Intercept)   2.6104 0.05744  2.49777  2.7229 0.000e+00
#> platelet      0.2450 0.08434  0.07971  0.4103 3.672e-03
#> tcell         0.1459 0.11873 -0.08680  0.3786 2.191e-01
#> age          -0.1729 0.03794 -0.24723 -0.0985 5.217e-06
estimate(rmst301)
#>             Estimate Std.Err    2.5%      97.5%    P-value
#> (Intercept)   2.4449 0.07291  2.3020  2.5878272 1.617e-246
#> platelet     -0.4519 0.16739 -0.7800 -0.1238273  6.940e-03
#> tcell        -0.4937 0.25146 -0.9866 -0.0008631  4.960e-02
#> age           0.2747 0.06538  0.1466  0.4028551  2.648e-05
estimate(rmst302)
#>             Estimate Std.Err    2.5%  97.5%   P-value
#> (Intercept)  1.45849  0.1362  1.1915 1.7255 9.623e-27
#> platelet    -0.01890  0.2224 -0.4548 0.4170 9.323e-01
#> tcell        0.40759  0.2648 -0.1113 0.9265 1.237e-01
#> age          0.04577  0.1191 -0.1877 0.2792 7.008e-01

## percentage of total cumulative incidence due to cause 1
rmtlratioI <- rmtlRatio(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=1)
summary(rmtlratioI)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  0.977110  0.178192  0.627860  1.326361  0.0000
#> platelet    -0.412122  0.323292 -1.045763  0.221520  0.2024
#> tcell       -0.855519  0.419219 -1.677174 -0.033865  0.0413
#> age          0.217628  0.168419 -0.112468  0.547724  0.1963
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  2.65677 1.87360 3.7673
#> platelet     0.66224 0.35142 1.2480
#> tcell        0.42506 0.18690 0.9667
#> age          1.24312 0.89363 1.7293
#> 
#> 

pp <- predict(rmtlratioI,bmt)
ppb <- cbind(pp,bmt)

## percentage of total cumulative incidence due to cause 1
cifratio <- binregRatio(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=1)
summary(cifratio)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)  0.91179  0.17554  0.56774  1.25585  0.0000
#> platelet    -0.46503  0.31924 -1.09072  0.16066  0.1452
#> tcell       -1.06112  0.41571 -1.87590 -0.24634  0.0107
#> age          0.21199  0.16115 -0.10386  0.52784  0.1884
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  2.48878 1.76427 3.5108
#> platelet     0.62812 0.33597 1.1743
#> tcell        0.34607 0.15322 0.7817
#> age          1.23614 0.90135 1.6953
#> 
#> 
pp <- predict(cifratio,bmt)

rmtlratioI <- binregRatio(Event(time,cause)~platelet+tcell+age,bmt,
                               time=30,cause=1,outcome="rmtl")
summary(rmtlratioI)
#>    n events
#>  408    154
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  0.977110  0.178192  0.627860  1.326361  0.0000
#> platelet    -0.412122  0.323292 -1.045763  0.221520  0.2024
#> tcell       -0.855519  0.419219 -1.677174 -0.033865  0.0413
#> age          0.217628  0.168419 -0.112468  0.547724  0.1963
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  2.65677 1.87360 3.7673
#> platelet     0.66224 0.35142 1.2480
#> tcell        0.42506 0.18690 0.9667
#> age          1.24312 0.89363 1.7293
#> 
#> 

pp <- predict(rmtlratioI,bmt)
ppb <- cbind(pp,bmt)
```
