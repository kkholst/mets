# GOF for Cox covariates in PH regression

Cumulative residuals after model matrix for Cox PH regression p-values
based on Lin, Wei, Ying resampling.

## Usage

``` r
gofM.phreg(
  formula,
  data,
  offset = NULL,
  weights = NULL,
  modelmatrix = NULL,
  n.sim = 1000,
  silent = 1,
  ...
)
```

## Arguments

- formula:

  formula for cox regression

- data:

  data for model

- offset:

  offset

- weights:

  weights

- modelmatrix:

  matrix for cumulating residuals

- n.sim:

  number of simulations for score processes

- silent:

  to keep it absolutely silent, otherwise timing estimate will be
  prduced for longer jobs.

- ...:

  Additional arguments to lower level funtions

## Details

That is, computes \$\$ U(t) = \int_0^t M^t d \hat M \$\$ and resamples
its asymptotic distribution.

This will show if the residuals are consistent with the model.
Typically, M will be a design matrix for the continous covariates that
gives for example the quartiles, and then the plot will show if for the
different quartiles of the covariate the risk prediction is consistent
over time (time x covariate interaction).

## Author

Thomas Scheike and Klaus K. Holst

## Examples

``` r
library(mets)
data(TRACE)
set.seed(1)
TRACEsam <- blocksample(TRACE,idvar="id",replace=FALSE,100)

dcut(TRACEsam)  <- ~. 
mm <- model.matrix(~-1+factor(wmicat.4),data=TRACEsam)
m1 <- gofM.phreg(Surv(time,status==9)~vf+chf+wmi,data=TRACEsam,modelmatrix=mm)
summary(m1)
#> Cumulative residuals versus modelmatrix :
#>                            Sup_t |U(t)|  pval
#> factor(wmicat.4)[0.4,1.1]      5.788752 0.021
#> factor(wmicat.4)(1.1,1.4]      2.633143 0.591
#> factor(wmicat.4)(1.4,1.72]     5.657370 0.038
#> factor(wmicat.4)(1.72,2]       1.227485 0.948
#> 
#> Cumulative score process versus covariates (discrete z via model.matrix):
#>         Sup_z |U(tau,z)|  pval
#> matrixZ         1.801626 0.812
if (interactive()) {
par(mfrow=c(2,2))
plot(m1)
}

m1 <- gofM.phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACEsam,modelmatrix=mm) 
summary(m1)
#> Cumulative residuals versus modelmatrix :
#>                            Sup_t |U(t)|  pval
#> factor(wmicat.4)[0.4,1.1]      6.049754 0.022
#> factor(wmicat.4)(1.1,1.4]      2.211031 0.757
#> factor(wmicat.4)(1.4,1.72]     5.040704 0.085
#> factor(wmicat.4)(1.72,2]       1.482109 0.890
#> 
#> Cumulative score process versus covariates (discrete z via model.matrix):
#>         Sup_z |U(tau,z)|  pval
#> matrixZ         1.228032 0.934

## cumulative sums in covariates, via design matrix mm 
mm <- cumContr(TRACEsam$wmi,breaks=10,equi=TRUE)
m1 <- gofM.phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACEsam,
      modelmatrix=mm,silent=0)
#> Cumulative score process test for modelmatrix:
#>        Sup_t |U(t)| pval
#> <=0.56         0.84 0.35
#> <=0.72         1.93 0.39
#> <=0.88         5.06 0.04
#> <=1.04         2.63 0.58
#> <=1.2          4.27 0.17
#> <=1.36         4.04 0.17
#> <=1.52         3.20 0.29
#> <=1.68         1.83 0.81
#> <=1.84         1.25 0.82
#> <=2            0.00 1.00
summary(m1)
#> Cumulative residuals versus modelmatrix :
#>        Sup_t |U(t)|  pval
#> <=0.56    0.8423856 0.352
#> <=0.72    1.9258405 0.394
#> <=0.88    5.0552596 0.037
#> <=1.04    2.6334412 0.579
#> <=1.2     4.2747875 0.174
#> <=1.36    4.0393945 0.169
#> <=1.52    3.2004455 0.292
#> <=1.68    1.8330704 0.810
#> <=1.84    1.2509429 0.823
#> <=2       0.0000000 1.000
#> 
#> Cumulative score process versus covariates (discrete z via model.matrix):
#>         Sup_z |U(tau,z)|  pval
#> matrixZ          3.73439 0.283
```
