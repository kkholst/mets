# IPCW Cumulative Cost

## Overview

We describe how to perform regression modelling for cumulative cost
\begin{align\*} {\cal U}(t) & = \int_0^t Z(s) dN(s) \end{align\*} where
N(s) is a counting process that registers the times at which costs are
realised and accumulated, and Z(t) is the cost (or mark) at the event
times. The counting process can be a mix of random and fixed times, and
the data are represented in counting process format with the marks/costs
attached to the event times. There are many additional uses of such
cumulative processes; for example, when considering time lost in a
recurrent events setting, which we return to below.

We can estimate the marginal mean of the cumulative process
\begin{align\*} \nu(t) & = E ( {\cal U}(t) ) \end{align\*} possibly for
strata with standard errors based on the derived influence function.

We provide semi-parametric regression modelling using the proportional
model \begin{align\*} E ( {\cal U}(t) \| X) & = \Lambda_0(t) \exp( X^T
\beta). \end{align\*}

In addition for a fixed time-point t \in \[0,\tau\] we can estimate the
mean given covariates \begin{align\*} E ( {\cal U}(t) \| X) & = \exp(
X^T \beta) \end{align\*} where \tau is some maximum follow-up time.

- These quantities are estimated under independent right-censoring given
  X, using IPCW-adjusted estimating equations,
  - similarly to the Ghosh-Lin model for recurrent events.
- A terminal event can be specified.

We also estimate the probability of exceeding thresholds over time
\begin{align\*} P ( {\cal U}(t) \> k ) & = \mu_k(t), \end{align\*} and
in the setting with a terminal event this is based on a derived
competing risks data structure that keeps track of the competing
terminal event.

Regression modelling of this quantity is also possible using competing
risks regression models, for example via the `cifreg` function in mets.

## HF-action data

Using the HF-action data, we simulate a severity score for each event.

``` r

library(mets)
data(hfactioncpx12)
hf <- hfactioncpx12
set.seed(1)
hf$severity <- abs((5+rnorm(741)*2))[hf$id]

## marginal mean using formula  
outNZ <- recurrent_marginal(Event(entry,time,status)~strata(treatment)+cluster(id)
             +marks(severity),hf,cause=1,death.code=2)
plot(outNZ,se=TRUE)
summary(outNZ,times=3) 
#> [[1]]
#>     new.time    mean        se CI-2.5% CI-97.5% strata
#> 682        3 10.3249 0.6044746  9.2056  11.5803      0
#> 
#> [[2]]
#>     new.time    mean        se  CI-2.5% CI-97.5% strata
#> 601        3 9.33613 0.6429995 8.157232 10.68541      1

outN <- recurrent_marginal(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,
               cause=1,death.code=2)
plot(outN,se=TRUE,add=TRUE)
```

![Marginal mean cumulative cost by
treatment.](figure/cumulative-cost-unnamed-chunk-2-1.png)

``` r

summary(outN,times=3) 
#> [[1]]
#>     new.time     mean        se  CI-2.5% CI-97.5% strata
#> 682        3 2.118496 0.1138572 1.906692 2.353829      0
#> 
#> [[2]]
#>     new.time     mean        se  CI-2.5% CI-97.5% strata
#> 601        3 1.924062 0.1216577 1.699801 2.177912      1
```

For comparison we also compute the IPCW estimates with and without marks
at time 3 using the linear model, and note that they are identical.
Standard errors are based on different formulae that are asymptotically
equivalent, and we note that they are very similar.

``` r

outNZ3 <- recregIPCW(Event(entry,time,status)~-1+treatment+cluster(id)+marks(severity),data=hf,
          cause=1,death.code=2,time=3,cens.model=~strata(treatment),model="lin")
summary(outNZ3)
#>    n events
#>  741   1281
#> 
#>  741 clusters
#> coeffients:
#>            Estimate  Std.Err     2.5%    97.5% P-value
#> treatment0 10.32490  0.60444  9.14022 11.50959       0
#> treatment1  9.33613  0.64295  8.07597 10.59629       0
head(iid(outNZ3))
#>            [,1]      [,2]
#> 1 -0.0013037837 0.0000000
#> 2  0.0089249278 0.0000000
#> 3  0.0000000000 0.0211858
#> 4 -0.0056613816 0.0000000
#> 5 -0.0005288663 0.0000000
#> 6 -0.0342554114 0.0000000

outN3 <- recregIPCW(Event(entry,time,status)~-1+treatment+cluster(id),data=hf,cause=1,death.code=2,time=3,
         cens.model=~strata(treatment),model="lin")
summary(outN3)
#>    n events
#>  741   1281
#> 
#>  741 clusters
#> coeffients:
#>            Estimate Std.Err    2.5%   97.5% P-value
#> treatment0  2.11850 0.11385 1.89535 2.34164       0
#> treatment1  1.92406 0.12165 1.68564 2.16248       0
head(iid(outN3))
#>            [,1]        [,2]
#> 1  0.0004542472 0.000000000
#> 2  0.0009756994 0.000000000
#> 3  0.0000000000 0.009301496
#> 4 -0.0029668336 0.000000000
#> 5 -0.0001120764 0.000000000
#> 6 -0.0070693971 0.000000000
```

We also apply the semiparametric proportional cost model with IPCW
adjustment:

``` r

propNZ <- recreg(Event(entry,time,status)~treatment+marks(severity)+cluster(id),data=hf,cause=1,death.code=2)
summary(propNZ) 
#> 
#>     n events
#>  2132   1391
#> 
#>  741 clusters
#> coefficients:
#>             Estimate      S.E.   dU^-1/2 P-value
#> treatment1 -0.101014  0.089005  0.024304  0.2564
#> 
#> exp(coefficients):
#>            Estimate    2.5%  97.5%
#> treatment1  0.90392 0.75922 1.0762
plot(propNZ,main="Baselines")
     
GL <- recreg(Event(entry,time,status)~treatment+cluster(id),hf,cause=1,death.code=2)
summary(GL)
#> 
#>     n events
#>  2132   1391
#> 
#>  741 clusters
#> coefficients:
#>             Estimate      S.E.   dU^-1/2 P-value
#> treatment1 -0.110404  0.078656  0.053776  0.1604
#> 
#> exp(coefficients):
#>            Estimate    2.5%  97.5%
#> treatment1  0.89547 0.76754 1.0447
plot(GL,add=TRUE,col=2)
```

![Proportional cost model baseline
estimates.](figure/cumulative-cost-unnamed-chunk-4-1.png) Those treated
have 14% lower cumulative severity and 11% lower expected number of
events.

## Exceed threshold

Finally, we estimate the probability of exceeding cumulative severity
thresholds of 1, 5, and 10:

``` r

ooNZ <- prob_exceed_recurrent(Event(entry,time,status)~strata(treatment)+cluster(id)+marks(severity),data=hf,
                  cause=1,death.code=2,exceed=c(1,5,10,20))
plot(ooNZ,strata=1)
plot(ooNZ,strata=2,add=TRUE)
```

![Probability of exceeding cumulative cost
thresholds.](figure/cumulative-cost-unnamed-chunk-5-1.png)

``` r

summary(ooNZ,times=3)
#> $`0`
#> $`0`$prob
#>            times           
#>                3 2.99865085
#> N<1            3 0.06474819
#> exceed>=1      3 0.93525181
#> exceed>=5      3 0.90118229
#> exceed>=10     3 0.74537955
#> exceed>=20     3 0.49282394
#> 
#> $`0`$se
#>            times           
#>                3 2.99865085
#> N<1            3 0.02657120
#> exceed>=1      3 0.02657120
#> exceed>=5      3 0.03013468
#> exceed>=10     3 0.04466101
#> exceed>=20     3 0.05149605
#> 
#> $`0`$lower
#>      times          
#> [1,]     3 2.9986509
#> [2,]     3 0.1154033
#> [3,]     3 0.8845967
#> [4,]     3 0.8440133
#> [5,]     3 0.6627899
#> [6,]     3 0.4015580
#> 
#> $`0`$upper
#>      times           
#> [1,]     3 2.99865085
#> [2,]     3 0.01119234
#> [3,]     3 0.98880766
#> [4,]     3 0.96222364
#> [5,]     3 0.83826056
#> [6,]     3 0.60483284
#> 
#> 
#> $`1`
#> $`1`$prob
#>            times           
#>                3 2.99865085
#> N<1            3 0.01770732
#> exceed>=1      3 0.98229268
#> exceed>=5      3 0.94934557
#> exceed>=10     3 0.79560979
#> exceed>=20     3 0.53319063
#> 
#> $`1`$se
#>            times           
#>                3 2.99865085
#> N<1            3 0.01097663
#> exceed>=1      3 0.01097663
#> exceed>=5      3 0.01645312
#> exceed>=10     3 0.03727391
#> exceed>=20     3 0.05545584
#> 
#> $`1`$lower
#>      times           
#> [1,]     3 2.99865085
#> [2,]     3 0.03898725
#> [3,]     3 0.96101275
#> [4,]     3 0.91763959
#> [5,]     3 0.72580801
#> [6,]     3 0.43486168
#> 
#> $`1`$upper
#>      times          
#> [1,]     3 2.9986509
#> [2,]     3 0.0000000
#> [3,]     3 1.0000000
#> [4,]     3 0.9821470
#> [5,]     3 0.8721245
#> [6,]     3 0.6537533
```

## Cumulative time lost for recurrent events

The cumulative time lost for recurrent events is defined as
\begin{align\*} {\cal M}(t) = E\left\[ \int_0^\tau (\tau-s) dN(s)
\right\] = \int_0^\tau \mu(s) ds \end{align\*} where \mu(t) = E( N(t) )
is the marginal mean of the recurrent events at time t.

``` r

hf$lost5 <- 5-hf$time

RecLost <- recregIPCW(Event(entry,time,status)~-1+treatment+cluster(id)+marks(lost5),data=hf,
           cause=1,death.code=2,time=5,cens.model=~strata(treatment),model="lin")
summary(RecLost)
#>    n events
#>  741   1391
#> 
#>  741 clusters
#> coeffients:
#>            Estimate Std.Err    2.5%   97.5% P-value
#> treatment0  8.58300 0.42951 7.74118 9.42482       0
#> treatment1  7.66234 0.46400 6.75292 8.57177       0
head(iid(RecLost))
#>            [,1]       [,2]
#> 1  0.0016920221 0.00000000
#> 2  0.0073388996 0.00000000
#> 3  0.0000000000 0.02120478
#> 4 -0.0095548150 0.00000000
#> 5 -0.0005696809 0.00000000
#> 6 -0.0201750011 0.00000000
```

## SessionInfo

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /home/kkzh/.asdf/installs/r/4.6.0/lib/R/lib/libRblas.so 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Copenhagen
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] timereg_2.0.7  survival_3.8-6 mets_1.3.10   
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6              knitr_1.51             rlang_1.2.0           
#>  [4] xfun_0.57              otel_0.2.0             future.apply_1.20.2   
#>  [7] listenv_0.10.1         lava_1.9.1             stats4_4.6.0          
#> [10] grid_4.6.0             evaluate_1.0.5         yaml_2.3.12           
#> [13] mvtnorm_1.3-7          numDeriv_2016.8-1.1    compiler_4.6.0        
#> [16] codetools_0.2-20       Rcpp_1.1.1-1.1         ucminf_1.2.3          
#> [19] future_1.70.0          lattice_0.22-9         digest_0.6.39         
#> [22] parallelly_1.47.0      parallel_4.6.0         splines_4.6.0         
#> [25] Matrix_1.7-5           tools_4.6.0            RcppArmadillo_15.2.6-1
#> [28] globals_0.19.1
```
