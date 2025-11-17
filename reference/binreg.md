# Binomial Regression for censored competing risks data

Simple version of comp.risk function of timereg for just one time-point
thus fitting the model \$\$P(T \leq t, \epsilon=1 \| X ) = expit( X^T
beta) \$\$

## Usage

``` r
binreg(
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
  outcome = c("cif", "rmst", "rmtl"),
  model = c("default", "logit", "exp", "lin"),
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

  "II" adds augmentation term, and "I" classic binomial regression

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

  can do CIF regression "cif"=F(t\|X), "rmst"=E( min(T, t) \| X) , or
  years-lost "rmtl"=E( I(epsilon==cause) ( t - mint(T,t)) ) \| X)

- model:

  link functions used, with defaults logit for cif, exp for rmst or
  rmtl, but can be logit, exp or lin (for identity link)

- Ydirect:

  use this Y instead of outcome constructed inside the program (e.g.
  I(T\< t, epsilon=1)), then uses IPCW vesion of the Y, set outcome to
  "rmst" to fit using the model specified by model

- ...:

  Additional arguments to lower level funtions

## Details

Based on binomial regresion IPCW response estimating equation: \$\$ X (
\Delta^{ipcw}(t) I(T \leq t, \epsilon=1 ) - expit( X^T beta)) = 0 \$\$
where \$\$\Delta^{ipcw}(t) = I((min(t,T)\< C)/G_c(min(t,T)-)\$\$ is IPCW
adjustment of the response \$\$Y(t)= I(T \leq t, \epsilon=1 )\$\$.

(type="I") sovlves this estimating equation using a stratified
Kaplan-Meier for the censoring distribution. For (type="II") the default
an additional censoring augmentation term \$\$X \int E(Y(t)\|
T\>s)/G_c(s) d \hat M_c\$\$ is added.

logitIPCW instead considers \$\$ X I(min(T_i,t) \< G_i)/G_c(min(T_i ,t))
( I(T \leq t, \epsilon=1 ) - expit( X^T beta)) = 0 \$\$ a standard
logistic regression with weights that adjust for IPCW.

The variance is based on the squared influence functions that are also
returned as the iid component. naive.var is variance under known
censoring model.

Censoring model may depend on strata (cens.model=~strata(gX)).

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt); bmt$time <- bmt$time+runif(408)*0.001
# logistic regresion with IPCW binomial regression 
out <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50)
summary(out)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.180339  0.126754 -0.428773  0.068095  0.1548
#> tcell       -0.418248  0.345428 -1.095275  0.258779  0.2260
#> platelet    -0.437644  0.240971 -0.909938  0.034650  0.0693
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.83499 0.65131 1.0705
#> tcell        0.65820 0.33445 1.2953
#> platelet     0.64556 0.40255 1.0353
#> 
#> 
head(iid(out))
#>              [,1]        [,2]        [,3]
#> [1,] -0.006946389 0.004004234 0.006177038
#> [2,] -0.006946389 0.004004234 0.006177038
#> [3,] -0.006946389 0.004004234 0.006177038
#> [4,] -0.006946389 0.004004234 0.006177038
#> [5,] -0.006946389 0.004004234 0.006177038
#> [6,] -0.006946389 0.004004234 0.006177038

predict(out,data.frame(tcell=c(0,1),platelet=c(1,1)),se=TRUE)
#>        pred         se     lower     upper
#> 1 0.3502403 0.04847344 0.2552323 0.4452483
#> 2 0.2618778 0.06969023 0.1252850 0.3984707

outs <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50,cens.model=~strata(tcell,platelet))
summary(outs)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.180703  0.127413 -0.430429  0.069022  0.1561
#> tcell       -0.365924  0.350632 -1.053150  0.321302  0.2967
#> platelet    -0.433487  0.240270 -0.904408  0.037433  0.0712
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.83468 0.65023 1.0715
#> tcell        0.69356 0.34884 1.3789
#> platelet     0.64824 0.40478 1.0381
#> 
#> 

## glm with IPCW weights 
outl <- logitIPCW(Event(time,cause)~tcell+platelet,bmt,time=50)
summary(outl)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.241513  0.131457 -0.499164  0.016138  0.0662
#> tcell       -0.344686  0.368367 -1.066671  0.377299  0.3494
#> platelet    -0.292897  0.262665 -0.807711  0.221918  0.2648
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.78544 0.60704 1.0163
#> tcell        0.70844 0.34415 1.4583
#> platelet     0.74610 0.44588 1.2485
#> 
#> 

##########################################
### risk-ratio of different causes #######
##########################################
data(bmt)
bmt$id <- 1:nrow(bmt)
bmt$status <- bmt$cause
bmt$strata <- 1
bmtdob <- bmt
bmtdob$strata <-2
bmtdob <- dtransform(bmtdob,status=1,cause==2)
bmtdob <- dtransform(bmtdob,status=2,cause==1)
###
bmtdob <- rbind(bmt,bmtdob)
dtable(bmtdob,cause+status~strata)
#> strata: 1
#> 
#>       status   0   1   2
#> cause                   
#> 0            160   0   0
#> 1              0 161   0
#> 2              0   0  87
#> ------------------------------------------------------------ 
#> strata: 2
#> 
#>       status   0   1   2
#> cause                   
#> 0            160   0   0
#> 1              0   0 161
#> 2              0  87   0

cif1 <- cif(Event(time,cause)~+1,bmt,cause=1)
cif2 <- cif(Event(time,cause)~+1,bmt,cause=2)
plot(cif1)
plot(cif2,add=TRUE,col=2)


cifs1 <- binreg(Event(time,cause)~tcell+platelet+age,bmt,cause=2,time=50)
cifs2 <- binreg(Event(time,cause)~tcell+platelet+age,bmt,cause=2,time=50)
summary(cifs1)
#>    n events
#>  408     85
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.321995  0.157777 -1.631233 -1.012758  0.0000
#> tcell        0.747294  0.352319  0.056761  1.437826  0.0339
#> platelet    -0.019451  0.277005 -0.562372  0.523469  0.9440
#> age         -0.072282  0.141683 -0.349976  0.205412  0.6099
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.26660 0.19569 0.3632
#> tcell        2.11128 1.05840 4.2115
#> platelet     0.98074 0.56986 1.6879
#> age          0.93027 0.70471 1.2280
#> 
#> 
summary(cifs2)
#>    n events
#>  408     85
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.321995  0.157777 -1.631233 -1.012758  0.0000
#> tcell        0.747294  0.352319  0.056761  1.437826  0.0339
#> platelet    -0.019451  0.277005 -0.562372  0.523469  0.9440
#> age         -0.072282  0.141683 -0.349976  0.205412  0.6099
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.26660 0.19569 0.3632
#> tcell        2.11128 1.05840 4.2115
#> platelet     0.98074 0.56986 1.6879
#> age          0.93027 0.70471 1.2280
#> 
#> 

cifdob <- binreg(Event(time,status)~-1+factor(strata)+
   tcell*factor(strata)+platelet*factor(strata)+age*factor(strata)
   +cluster(id),bmtdob,cause=1,time=50,cens.model=~strata(strata))
summary(cifdob)
#>    n events
#>  816    245
#> 
#>  408 clusters
#> coeffients:
#>                           Estimate   Std.Err      2.5%     97.5% P-value
#> factor(strata)1          -0.198956  0.130987 -0.455685  0.057773  0.1288
#> factor(strata)2          -1.321995  0.157777 -1.631233 -1.012758  0.0000
#> tcell                    -0.637010  0.356610 -1.335954  0.061934  0.0741
#> platelet                 -0.344885  0.246012 -0.827058  0.137289  0.1609
#> age                       0.437259  0.107266  0.227021  0.647497  0.0000
#> factor(strata)2:tcell     1.384304  0.600983  0.206399  2.562209  0.0213
#> factor(strata)2:platelet  0.325433  0.432030 -0.521331  1.172197  0.4513
#> factor(strata)2:age      -0.509541  0.208100 -0.917410 -0.101672  0.0143
#> 
#> exp(coeffients):
#>                          Estimate    2.5%   97.5%
#> factor(strata)1           0.81959 0.63401  1.0595
#> factor(strata)2           0.26660 0.19569  0.3632
#> tcell                     0.52887 0.26291  1.0639
#> platelet                  0.70830 0.43733  1.1472
#> age                       1.54846 1.25486  1.9108
#> factor(strata)2:tcell     3.99205 1.22924 12.9644
#> factor(strata)2:platelet  1.38463 0.59373  3.2291
#> factor(strata)2:age       0.60077 0.39955  0.9033
#> 
#> 

riskratio <- function(p) {
  Z <- rbind(c(1,0,1,1,0,0,0,0), c(0,1,1,1,0,1,1,0))
  lp <- c(Z %*% p)
  p <- lava::expit(lp)
  return(p[1]/p[2])
}

lava::estimate(cifdob,f=riskratio)
#>    Estimate Std.Err   2.5% 97.5% P-value
#> p1   0.6604  0.2738 0.1238 1.197 0.01585
```
