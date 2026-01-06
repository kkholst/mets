# Cross-odds-ratio, OR or RR risk regression for competing risks

Fits a parametric model for the log-cross-odds-ratio for the predictive
effect of for the cumulative incidence curves for \\T_1\\ experiencing
cause i given that \\T_2\\ has experienced a cause k : \$\$
\log(COR(i\|k)) = h(\theta,z_1,i,z_2,k,t)=\_{default} \theta^T z = \$\$
with the log cross odds ratio being \$\$ COR(i\|k) = \frac{O(T_1 \leq
t,cause_1=i \| T_2 \leq t,cause_2=k)}{ O(T_1 \leq t,cause_1=i)} \$\$ the
conditional odds divided by the unconditional odds, with the odds being,
respectively \$\$ O(T_1 \leq t,cause_1=i \| T_2 \leq t,cause_1=k) =
\frac{ P_x(T_1 \leq t,cause_1=i \| T_2 \leq t,cause_2=k)}{ P_x((T_1 \leq
t,cause_1=i)^c \| T_2 \leq t,cause_2=k)} \$\$ and \$\$ O(T_1 \leq
t,cause_1=i) = \frac{P_x(T_1 \leq t,cause_1=i )}{P_x((T_1 \leq
t,cause_1=i)^c )}. \$\$ Here \\B^c\\ is the complement event of \\B\\,
\\P_x\\ is the distribution given covariates (\\x\\ are subject specific
and \\z\\ are cluster specific covariates), and \\h()\\ is a function
that is the simple identity \\\theta^T z\\ by default.

## Usage

``` r
cor.cif(
  cif,
  data,
  cause = NULL,
  times = NULL,
  cause1 = 1,
  cause2 = 1,
  cens.code = NULL,
  cens.model = "KM",
  Nit = 40,
  detail = 0,
  clusters = NULL,
  theta = NULL,
  theta.des = NULL,
  step = 1,
  sym = 0,
  weights = NULL,
  par.func = NULL,
  dpar.func = NULL,
  dimpar = NULL,
  score.method = "nlminb",
  same.cens = FALSE,
  censoring.weights = NULL,
  silent = 1,
  ...
)
```

## Arguments

- cif:

  a model object from the timereg::comp.risk function with the marginal
  cumulative incidence of cause1, i.e., the event of interest, and whose
  odds the comparision is compared to the conditional odds given cause2

- data:

  a data.frame with the variables.

- cause:

  specifies the causes related to the death times, the value cens.code
  is the censoring value. When missing it comes from marginal cif.

- times:

  time-vector that specifies the times used for the estimating euqations
  for the cross-odds-ratio estimation.

- cause1:

  specificies the cause considered.

- cause2:

  specificies the cause that is conditioned on.

- cens.code:

  specificies the code for the censoring if NULL then uses the one from
  the marginal cif model.

- cens.model:

  specified which model to use for the ICPW, KM is Kaplan-Meier
  alternatively it may be "cox"

- Nit:

  number of iterations for Newton-Raphson algorithm.

- detail:

  if 0 no details are printed during iterations, if 1 details are given.

- clusters:

  specifies the cluster structure.

- theta:

  specifies starting values for the cross-odds-ratio parameters of the
  model.

- theta.des:

  specifies a regression design for the cross-odds-ratio parameters.

- step:

  specifies the step size for the Newton-Raphson algorithm.

- sym:

  specifies if symmetry is used in the model.

- weights:

  weights for estimating equations.

- par.func:

  parfunc

- dpar.func:

  dparfunc

- dimpar:

  dimpar

- score.method:

  "nlminb", can also use "nr".

- same.cens:

  if true then censoring within clusters are assumed to be the same
  variable, default is independent censoring.

- censoring.weights:

  these probabilities are used for the bivariate censoring dist.

- silent:

  1 to suppress output about convergence related issues.

- ...:

  Not used.

## Value

returns an object of type 'cor'. With the following arguments:

- theta:

  estimate of proportional odds parameters of model.

- var.theta:

  variance for gamma.

- hess:

  the derivative of the used score.

- score:

  scores at final stage.

- score:

  scores at final stage.

- theta.iid:

  matrix of iid decomposition of parametric effects.

## Details

The OR dependence measure is given by \$\$ OR(i,k) = \log ( \frac{O(T_1
\leq t,cause_1=i \| T_2 \leq t,cause_2=k)}{ O(T_1 \leq t,cause_1=i) \|
T_2 \leq t,cause_2=k)} \$\$ This measure is numerically more stabile
than the COR measure, and is symetric in i,k.

The RR dependence measure is given by \$\$ RR(i,k) = \log ( \frac{P(T_1
\leq t,cause_1=i , T_2 \leq t,cause_2=k)}{ P(T_1 \leq t,cause_1=i) P(T_2
\leq t,cause_2=k)} \$\$ This measure is numerically more stabile than
the COR measure, and is symetric in i,k.

The model is fitted under symmetry (sym=1), i.e., such that it is
assumed that \\T_1\\ and \\T_2\\ can be interchanged and leads to the
same cross-odd-ratio (i.e. \\COR(i\|k) = COR(k\|i))\\, as would be
expected for twins or without symmetry as might be the case with mothers
and daughters (sym=0).

\\h()\\ may be specified as an R-function of the parameters, see example
below, but the default is that it is simply \\\theta^T z\\.

## References

Cross odds ratio Modelling of dependence for Multivariate Competing
Risks Data, Scheike and Sun (2012), Biostatistics.

A Semiparametric Random Effects Model for Multivariate Competing Risks
Data, Scheike, Zhang, Sun, Jensen (2010), Biometrika.

## Author

Thomas Scheike

## Examples

``` r
 ## Reduce Ex.Timings
library("timereg")
data(multcif);
multcif$cause[multcif$cause==0] <- 2
zyg <- rep(rbinom(200,1,0.5),each=2)
theta.des <- model.matrix(~-1+factor(zyg))

times=seq(0.05,1,by=0.05) # to speed up computations use only these time-points
add <- timereg::comp.risk(Event(time,cause)~+1+cluster(id),data=multcif,cause=1,
               n.sim=0,times=times,model="fg",max.clust=NULL)
add2 <- timereg::comp.risk(Event(time,cause)~+1+cluster(id),data=multcif,cause=2,
               n.sim=0,times=times,model="fg",max.clust=NULL)

out1 <- cor.cif(add,data=multcif,cause1=1,cause2=1)
summary(out1)
#> Cross odds ratio dependence for competing risks
#> 
#> Effect of cause1=1 on cause2=1 under symmetry=0
#> 
#>           log-Coef.   SE     z P-val Cross odds ratio   SE
#> intercept      1.04 9.33 0.111 0.911             2.82 26.4
#> 

out2 <- cor.cif(add,data=multcif,cause1=1,cause2=1,theta.des=theta.des)
summary(out2)
#> Cross odds ratio dependence for competing risks
#> 
#> Effect of cause1=1 on cause2=1 under symmetry=0
#> 
#>              log-Coef.    SE      z P-val Cross odds ratio   SE
#> factor(zyg)0     0.882  8.96 0.0984 0.922             2.42 21.6
#> factor(zyg)1     1.230 21.80 0.0562 0.955             3.41 74.4
#> 

##out3 <- cor.cif(add,data=multcif,cause1=1,cause2=2,cif2=add2)
##summary(out3)
###########################################################
# investigating further models using parfunc and dparfunc
###########################################################
set.seed(100)
prt<-simnordic.random(2000,cordz=2,cormz=5)
prt$status <-prt$cause
table(prt$status)
#> 
#>    0    1    2 
#> 2102  597 5301 

times <- seq(40,100,by=10)
cifmod <- timereg::comp.risk(Event(time,cause)~+1+cluster(id),data=prt,
                    cause=1,n.sim=0,
                    times=times,conservative=1,max.clust=NULL,model="fg")
theta.des <- model.matrix(~-1+factor(zyg),data=prt)

parfunc <- function(par,t,pardes)
{
par <- pardes %*% c(par[1],par[2]) +
       pardes %*% c( par[3]*(t-60)/12,par[4]*(t-60)/12)
par
}
head(parfunc(c(0.1,1,0.1,1),50,theta.des))
#>        [,1]
#> 1 0.1666667
#> 2 0.1666667
#> 3 0.1666667
#> 4 0.1666667
#> 5 0.1666667
#> 6 0.1666667

dparfunc <- function(par,t,pardes)
{
dpar <- cbind(pardes, t(t(pardes) * c( (t-60)/12,(t-60)/12)) )
dpar
}
head(dparfunc(c(0.1,1,0.1,1),50,theta.des))
#>   factor(zyg)MZ factor(zyg)DZ factor(zyg)MZ factor(zyg)DZ
#> 1             0             1             0    -0.8333333
#> 2             0             1             0    -0.8333333
#> 3             0             1             0    -0.8333333
#> 4             0             1             0    -0.8333333
#> 5             0             1             0    -0.8333333
#> 6             0             1             0    -0.8333333

names(prt)
#>  [1] "time"      "cause"     "x"         "country"   "id"        "cens"     
#>  [7] "stime"     "entry"     "truncated" "zyg"       "status"   
or1 <- or.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,
              same.cens=TRUE,theta=c(0.6,1.1,0.1,0.1),
              par.func=parfunc,dpar.func=dparfunc,dimpar=4,
              score.method="nr",detail=1)
#> [1] "Fisher-Scoring ===================: it= 1"
#> theta:[1] 0.6 1.1 0.1 0.1
#> score:[1]  4.4038139 -0.5697823 11.6772976 -1.6949449
#> hess:          [,1]      [,2]       [,3]       [,4]
#> [1,] -2.380951  0.000000  -6.184041   0.000000
#> [2,]  0.000000 -3.262648   0.000000  -8.285153
#> [3,] -6.184041  0.000000 -17.409502   0.000000
#> [4,]  0.000000 -8.285153   0.000000 -23.052048
#> [1] "Fisher-Scoring ===================: it= 2"
#> theta:[1]  1.98846481  1.23830300  0.27754528 -0.02323445
#> score:[1] -1.25603529 -0.01146765 -2.79344759 -0.04626781
#> hess:          [,1]      [,2]       [,3]       [,4]
#> [1,] -3.955147  0.000000  -9.393342   0.000000
#> [2,]  0.000000 -3.009163   0.000000  -7.530575
#> [3,] -9.393342  0.000000 -25.237970   0.000000
#> [4,]  0.000000 -7.530575   0.000000 -20.789841
#> [1] "Fisher-Scoring ===================: it= 3"
#> theta:[1]  1.51716906  1.25710704  0.34227293 -0.03227122
#> score:[1]  0.0084015577  0.0001720138  0.0266070641 -0.0006890471
#> hess:          [,1]      [,2]       [,3]       [,4]
#> [1,] -3.976728  0.000000  -9.711857   0.000000
#> [2,]  0.000000 -3.004696   0.000000  -7.507249
#> [3,] -9.711857  0.000000 -26.396637   0.000000
#> [4,]  0.000000 -7.507249   0.000000 -20.707731
#> [1] "Fisher-Scoring ===================: it= 4"
#> theta:[1]  1.51373019  1.25859721  0.34454613 -0.03284474
#> score:[1]  1.000615e-04  1.638158e-05  2.676517e-04 -3.827292e-05
#> hess:          [,1]      [,2]      [,3]       [,4]
#> [1,] -3.974979  0.000000  -9.70733   0.000000
#> [2,]  0.000000 -3.004947   0.00000  -7.507026
#> [3,] -9.707330  0.000000 -26.38239   0.000000
#> [4,]  0.000000 -7.507026   0.00000 -20.705898
#> [1] "Fisher-Scoring ===================: it= 5"
#> theta:[1]  1.51373411  1.25870404  0.34455483 -0.03288531
#> score:[1]  7.061280e-08  1.183464e-06  1.291720e-06 -2.693329e-06
#> hess:          [,1]      [,2]       [,3]       [,4]
#> [1,] -3.974972  0.000000  -9.707299   0.000000
#> [2,]  0.000000 -3.004967   0.000000  -7.507016
#> [3,] -9.707299  0.000000 -26.382278   0.000000
#> [4,]  0.000000 -7.507016   0.000000 -20.705785
#> [1] "Fisher-Scoring ===================: it= 6"
#> theta:[1]  1.51373311  1.25871167  0.34455525 -0.03288821
#> score:[1]  3.249126e-08  8.427829e-08  5.920297e-08 -1.934684e-07
#> hess:          [,1]      [,2]       [,3]       [,4]
#> [1,] -3.974972  0.000000  -9.707298   0.000000
#> [2,]  0.000000 -3.004969   0.000000  -7.507015
#> [3,] -9.707298  0.000000 -26.382277   0.000000
#> [4,]  0.000000 -7.507015   0.000000 -20.705777
summary(or1)
#> OR for dependence for competing risks
#> 
#> OR of cumulative incidence for cause1= 1  and cause2= 1
#>        log-ratio Coef.    SE      z    P-val Ratio    SE
#> R-func          1.5100 0.348  4.350 1.39e-05 4.540 1.580
#> R-func          1.2600 0.408  3.090 2.03e-03 3.520 1.440
#> R-func          0.3450 0.129  2.680 7.43e-03 1.410 0.182
#> R-func         -0.0329 0.111 -0.296 7.67e-01 0.968 0.107
#> 

 cor1 <- cor.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,
                 same.cens=TRUE,theta=c(0.5,1.0,0.1,0.1),
                 par.func=parfunc,dpar.func=dparfunc,dimpar=4,
                 control=list(trace=TRUE),detail=1)
#>   0:     764.39429: 0.500000  1.00000 0.100000 0.100000
#>   1:     760.41547: 0.638110 0.960514 0.467772 -0.0132205
#>   2:     760.40387: 0.653394 0.974160 0.451035 0.0181984
#>   3:     760.38389: 0.690757 0.972239 0.448591 0.00142152
#>   4:     760.37606: 0.718635 0.978903 0.420006 0.00837761
#>   5:     760.25923:  1.06302  1.01086 0.300325 -0.0208226
#>   6:     760.25074:  1.17117  1.03359 0.291034 0.00480685
#>   7:     760.21950:  1.27004  1.04554 0.238147 -0.0108041
#>   8:     760.20998:  1.37486  1.06362 0.202831 -0.0306896
#>   9:     760.20957:  1.45490  1.09158 0.167676 -0.0343349
#>  10:     760.20850:  1.41961  1.09023 0.182575 -0.0388101
#>  11:     760.20845:  1.41809  1.09720 0.183351 -0.0412561
#>  12:     760.20842:  1.41724  1.11066 0.183663 -0.0460526
#>  13:     760.20842:  1.41733  1.11052 0.183624 -0.0459859
#> $par
#> [1]  1.4173313  1.1105171  0.1836238 -0.0459859
#> 
#> $objective
#> [1] 760.2084
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 13
#> 
#> $evaluations
#> function gradient 
#>       19       71 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
#> iid decomposition
summary(cor1)
#> Cross odds ratio dependence for competing risks
#> 
#> Effect of cause1=1 on cause2=1 under symmetry=0
#> 
#>        log-Coef.   SE        z P-val Cross odds ratio  SE
#> R-func     1.420 66.5  0.02130 0.983            4.130 274
#> R-func     1.110 71.3  0.01560 0.988            3.040 217
#> R-func     0.184 23.3  0.00789 0.994            1.200  28
#> R-func    -0.046 22.0 -0.00209 0.998            0.955  21
#> 

### piecewise contant OR model
gparfunc <- function(par,t,pardes)
{
  cuts <- c(0,80,90,120)
  grop <- diff(t<cuts)
paru  <- (pardes[,1]==1) * sum(grop*par[1:3]) +
    (pardes[,2]==1) * sum(grop*par[4:6])
paru
}

dgparfunc <- function(par,t,pardes)
{
  cuts <- c(0,80,90,120)
  grop <- diff(t<cuts)
par1 <- matrix(c(grop),nrow(pardes),length(grop),byrow=TRUE)
parmz <- par1* (pardes[,1]==1)
pardz <- (pardes[,2]==1) * par1
dpar <- cbind( parmz,pardz)
dpar
}
head(dgparfunc(rep(0.1,6),50,theta.des))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    0    0    1    0    0
#> [2,]    0    0    0    1    0    0
#> [3,]    0    0    0    1    0    0
#> [4,]    0    0    0    1    0    0
#> [5,]    0    0    0    1    0    0
#> [6,]    0    0    0    1    0    0
head(gparfunc(rep(0.1,6),50,theta.des))
#>   1   2   3   4   5   6 
#> 0.1 0.1 0.1 0.1 0.1 0.1 

or1g <- or.cif(cifmod,data=prt,cause1=1,cause2=1,
               theta.des=theta.des, same.cens=TRUE,
               par.func=gparfunc,dpar.func=dgparfunc,
               dimpar=6,score.method="nr",detail=1)
#> [1] "Fisher-Scoring ===================: it= 1"
#> theta:[1] 0.1 0.1 0.1 0.1 0.1 0.1
#> score:[1] 0.1272695 0.6721568 3.2995329 0.1193734 0.3414362 1.1416387
#> hess:            [,1]       [,2]       [,3]        [,4]       [,5]       [,6]
#> [1,] -0.05264815  0.0000000  0.0000000  0.00000000  0.0000000  0.0000000
#> [2,]  0.00000000 -0.2202109  0.0000000  0.00000000  0.0000000  0.0000000
#> [3,]  0.00000000  0.0000000 -0.8603209  0.00000000  0.0000000  0.0000000
#> [4,]  0.00000000  0.0000000  0.0000000 -0.05264815  0.0000000  0.0000000
#> [5,]  0.00000000  0.0000000  0.0000000  0.00000000 -0.2202109  0.0000000
#> [6,]  0.00000000  0.0000000  0.0000000  0.00000000  0.0000000 -0.8603209
#> [1] "Fisher-Scoring ===================: it= 2"
#> theta:[1] 2.517359 3.152333 3.935235 2.367380 1.650497 1.426992
#> score:[1] -0.4338385 -1.0147952 -2.5460332 -0.3722739 -0.3313216 -0.7014266
#> hess:           [,1]       [,2]      [,3]       [,4]       [,5]      [,6]
#> [1,] -0.4976024  0.0000000  0.000000  0.0000000  0.0000000  0.000000
#> [2,]  0.0000000 -0.8549109  0.000000  0.0000000  0.0000000  0.000000
#> [3,]  0.0000000  0.0000000 -1.600519  0.0000000  0.0000000  0.000000
#> [4,]  0.0000000  0.0000000  0.000000 -0.4751796  0.0000000  0.000000
#> [5,]  0.0000000  0.0000000  0.000000  0.0000000 -0.8248922  0.000000
#> [6,]  0.0000000  0.0000000  0.000000  0.0000000  0.0000000 -2.408282
#> [1] "Fisher-Scoring ===================: it= 3"
#> theta:[1] 1.645501 1.965314 2.344480 1.583942 1.248842 1.135736
#> score:[1] -0.02982302  0.04411753  0.63062461 -0.02738202 -0.01381585 -0.01990217
#> hess:           [,1]       [,2]      [,3]      [,4]       [,5]      [,6]
#> [1,] -0.3202401  0.0000000  0.000000  0.000000  0.0000000  0.000000
#> [2,]  0.0000000 -0.9097494  0.000000  0.000000  0.0000000  0.000000
#> [3,]  0.0000000  0.0000000 -2.766973  0.000000  0.0000000  0.000000
#> [4,]  0.0000000  0.0000000  0.000000 -0.305418  0.0000000  0.000000
#> [5,]  0.0000000  0.0000000  0.000000  0.000000 -0.6707065  0.000000
#> [6,]  0.0000000  0.0000000  0.000000  0.000000  0.0000000 -2.098342
#> [1] "Fisher-Scoring ===================: it= 4"
#> theta:[1] 1.552374 2.013808 2.572392 1.494288 1.228243 1.126251
#> score:[1] -4.297262e-04 -1.168380e-04  3.731020e-03 -5.158185e-04 -4.497370e-05
#> [6]  2.252054e-06
#> hess:           [,1]       [,2]      [,3]       [,4]       [,5]      [,6]
#> [1,] -0.2978324  0.0000000  0.000000  0.0000000  0.0000000  0.000000
#> [2,]  0.0000000 -0.9191812  0.000000  0.0000000  0.0000000  0.000000
#> [3,]  0.0000000  0.0000000 -2.695765  0.0000000  0.0000000  0.000000
#> [4,]  0.0000000  0.0000000  0.000000 -0.2839372  0.0000000  0.000000
#> [5,]  0.0000000  0.0000000  0.000000  0.0000000 -0.6619257  0.000000
#> [6,]  0.0000000  0.0000000  0.000000  0.0000000  0.0000000 -2.087228
#> [1] "Fisher-Scoring ===================: it= 5"
#> theta:[1] 1.550931 2.013681 2.573776 1.492471 1.228175 1.126252
#> score:[1]  1.115720e-06 -7.362229e-10  1.590683e-06 -1.109555e-06 -4.887707e-10
#> [6] -2.862365e-09
#> hess:           [,1]       [,2]      [,3]       [,4]       [,5]     [,6]
#> [1,] -0.2974862  0.0000000  0.000000  0.0000000  0.0000000  0.00000
#> [2,]  0.0000000 -0.9191579  0.000000  0.0000000  0.0000000  0.00000
#> [3,]  0.0000000  0.0000000 -2.695165  0.0000000  0.0000000  0.00000
#> [4,]  0.0000000  0.0000000  0.000000 -0.2835044  0.0000000  0.00000
#> [5,]  0.0000000  0.0000000  0.000000  0.0000000 -0.6618966  0.00000
#> [6,]  0.0000000  0.0000000  0.000000  0.0000000  0.0000000 -2.08723
#> [1] "Fisher-Scoring ===================: it= 6"
#> theta:[1] 1.550935 2.013681 2.573776 1.492467 1.228175 1.126252
#> score:[1] -3.212719e-09 -6.107094e-15  8.052016e-10 -1.957429e-09 -2.185535e-15
#> [6]  7.569806e-11
#> hess:           [,1]       [,2]      [,3]       [,4]       [,5]     [,6]
#> [1,] -0.2974871  0.0000000  0.000000  0.0000000  0.0000000  0.00000
#> [2,]  0.0000000 -0.9191578  0.000000  0.0000000  0.0000000  0.00000
#> [3,]  0.0000000  0.0000000 -2.695165  0.0000000  0.0000000  0.00000
#> [4,]  0.0000000  0.0000000  0.000000 -0.2835035  0.0000000  0.00000
#> [5,]  0.0000000  0.0000000  0.000000  0.0000000 -0.6618966  0.00000
#> [6,]  0.0000000  0.0000000  0.000000  0.0000000  0.0000000 -2.08723
summary(or1g)
#> OR for dependence for competing risks
#> 
#> OR of cumulative incidence for cause1= 1  and cause2= 1
#>        log-ratio Coef.    SE    z    P-val Ratio    SE
#> R-func            1.55 0.362 4.28 1.87e-05  4.72 1.710
#> R-func            2.01 0.325 6.20 5.64e-10  7.49 2.430
#> R-func            2.57 0.352 7.31 2.71e-13 13.10 4.620
#> R-func            1.49 0.369 4.05 5.12e-05  4.45 1.640
#> R-func            1.23 0.307 3.99 6.49e-05  3.41 1.050
#> R-func            1.13 0.274 4.11 4.01e-05  3.08 0.846
#> 
names(or1g)
#>  [1] "theta"      "score"      "hess"       "hessi"      "var.theta" 
#>  [6] "theta.iid"  "score1"     "thetanames" "brierscore" "p11"       
#> [11] "call"      
head(or1g$theta.iid)
#>      [,1] [,2] [,3]         [,4]         [,5]         [,6]
#> [1,]    0    0    0 0.0007341431 0.0008025713 0.0008395676
#> [2,]    0    0    0 0.0007341431 0.0008025713 0.0008395676
#> [3,]    0    0    0 0.0007341431 0.0008025713 0.0008395676
#> [4,]    0    0    0 0.0007341431 0.0008025713 0.0008395676
#> [5,]    0    0    0 0.0007341431 0.0008025713 0.0008395676
#> [6,]    0    0    0 0.0007341431 0.0008025713 0.0008395676
```
