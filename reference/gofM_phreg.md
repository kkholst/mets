# Goodness-of-Fit for Cox Covariates (Model Matrix)

Tests the functional form of covariates in a Cox PH model by computing
cumulative residuals against a user-specified model matrix. This helps
detect non-linear effects or time-varying coefficients (interaction with
time).

## Usage

``` r
gofM_phreg(
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

  Formula for the Cox regression model.

- data:

  Data frame.

- offset:

  Offset vector.

- weights:

  Weights vector.

- modelmatrix:

  Matrix for cumulating residuals. Typically constructed using
  [`cumContr()`](http://kkholst.github.io/mets/reference/gofZ_phreg.md)
  or manually (e.g., quartiles of a continuous covariate).

- n.sim:

  Number of simulations (default 1000).

- silent:

  Logical; suppresses timing estimates if TRUE.

- ...:

  Additional arguments passed to `phreg`.

## Value

An object of class `"gof.phreg"` containing:

- jumptimes:

  Event times.

- supUsim:

  Simulated supremum values.

- res:

  Matrix with observed supremum and p-values for each column of
  `modelmatrix`.

- score:

  Cumulative score process.

- simUt:

  Simulated processes.

- Utlast, pval.last:

  Supremum and p-value for the final time point (covariate direction).

- type:

  Type of test ("modelmatrix").

## Details

The test statistic is: \$\$ U(t) = \int_0^t M^T d \hat M \$\$ where
\\M\\ is the model matrix (e.g., a set of basis functions for a
continuous covariate) and \\\hat M\\ are the martingale residuals.

P-values are based on the Lin, Wei, and Ying (1993) resampling method.
The plot shows whether the residuals are consistent with the model
across the range of the covariate.

## References

Lin, D. Y., Wei, L. J., & Ying, Z. (1993). Checking the Cox model with
cumulative sums of martingale-based residuals. Biometrika, 80(3),
557-572.

## See also

[`gof.phreg`](http://kkholst.github.io/mets/reference/gof.phreg.md),
[`gofZ_phreg`](http://kkholst.github.io/mets/reference/gofZ_phreg.md),
[`cumContr`](http://kkholst.github.io/mets/reference/gofZ_phreg.md)

## Author

Thomas Scheike and Klaus K. Holst

## Examples

``` r
data(TRACE)
set.seed(1)
TRACEsam <- blocksample(TRACE, idvar="id", replace=FALSE, 100)
dcut(TRACEsam) <- ~. 
mm <- model.matrix(~-1+factor(wmicat.4), data=TRACEsam)
m1 <- gofM_phreg(Surv(time,status==9)~vf+chf+wmi, data=TRACEsam, modelmatrix=mm)
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

## Cumulative sums in covariates via design matrix
mm <- mets:::cumContr(TRACEsam$wmi, breaks=10, equi=TRUE)
m1 <- gofM_phreg(Surv(time,status==9)~strata(vf)+chf+wmi, data=TRACEsam,
         modelmatrix=mm, silent=0)
#> Cumulative score process test for modelmatrix:
#>        Sup_t |U(t)| pval
#> <=0.56         0.88 0.33
#> <=0.72         2.25 0.26
#> <=0.88         5.28 0.03
#> <=1.04         2.68 0.54
#> <=1.2          4.12 0.17
#> <=1.36         4.09 0.17
#> <=1.52         3.34 0.24
#> <=1.68         1.78 0.83
#> <=1.84         1.17 0.84
#> <=2            0.00 1.00
summary(m1)
#> Cumulative residuals versus modelmatrix :
#>        Sup_t |U(t)|  pval
#> <=0.56    0.8821214 0.327
#> <=0.72    2.2521756 0.263
#> <=0.88    5.2766972 0.031
#> <=1.04    2.6807567 0.542
#> <=1.2     4.1217296 0.174
#> <=1.36    4.0926272 0.173
#> <=1.52    3.3447151 0.242
#> <=1.68    1.7766131 0.831
#> <=1.84    1.1652990 0.839
#> <=2       0.0000000 1.000
#> 
#> Cumulative score process versus covariates (discrete z via model.matrix):
#>         Sup_z |U(tau,z)| pval
#> matrixZ         3.931466 0.23
```
