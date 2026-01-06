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
#> factor(wmicat.4)[0.4,1.1]      5.679902 0.035
#> factor(wmicat.4)(1.1,1.4]      2.557399 0.593
#> factor(wmicat.4)(1.4,1.72]     5.411115 0.058
#> factor(wmicat.4)(1.72,2]       1.384392 0.931
#> 
#> Cumulative score process versus covariates (discrete z via model.matrix):
#>         Sup_z |U(tau,z)|  pval
#> matrixZ         1.771578 0.812

## cumulative sums in covariates, via design matrix mm 
mm <- cumContr(TRACEsam$wmi,breaks=10,equi=TRUE)
m1 <- gofM.phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACEsam,
      modelmatrix=mm,silent=0)
#> Cumulative score process test for modelmatrix:
#>        Sup_t |U(t)| pval
#> <=0.56         0.88 0.33
#> <=0.72         2.25 0.28
#> <=0.88         5.28 0.03
#> <=1.04         2.68 0.55
#> <=1.2          4.12 0.20
#> <=1.36         4.09 0.16
#> <=1.52         3.34 0.27
#> <=1.68         1.78 0.84
#> <=1.84         1.17 0.85
#> <=2            0.00 1.00
summary(m1)
#> Cumulative residuals versus modelmatrix :
#>        Sup_t |U(t)|  pval
#> <=0.56    0.8821214 0.330
#> <=0.72    2.2521756 0.276
#> <=0.88    5.2766972 0.029
#> <=1.04    2.6807567 0.552
#> <=1.2     4.1217296 0.201
#> <=1.36    4.0926272 0.161
#> <=1.52    3.3447151 0.271
#> <=1.68    1.7766131 0.844
#> <=1.84    1.1652990 0.853
#> <=2       0.0000000 1.000
#> 
#> Cumulative score process versus covariates (discrete z via model.matrix):
#>         Sup_z |U(tau,z)|  pval
#> matrixZ         3.931466 0.225
```
