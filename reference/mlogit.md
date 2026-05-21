# Multinomial Regression Based on phreg

Fits a multinomial regression model for a categorical outcome with \\K\\
levels: \$\$ P_i = \frac{ \exp( X^\top \beta_i ) }{ \sum\_{j=1}^K \exp(
X^\top \beta_j ) } \$\$ for \\i=1, \dots, K\\, where \\\beta_1 = 0\\
(baseline category).

## Usage

``` r
mlogit(formula, data, offset = NULL, weights = NULL, fix.X = FALSE, ...)
```

## Arguments

- formula:

  Formula with the outcome (similar to `coxph`). The outcome must be a
  factor.

- data:

  Data frame containing the variables.

- offset:

  Offsets for the partial likelihood.

- weights:

  Weights for the score equations.

- fix.X:

  Logical; if `TRUE`, forces the same coefficients for all categories
  (except intercepts).

- ...:

  Additional arguments passed to lower-level functions.

## Value

An object of class `"mlogit"` (extending `"phreg"`) containing:

- coef:

  Matrix of estimated coefficients (rows correspond to categories,
  columns to covariates).

- var:

  Robust variance-covariance matrix.

- iid:

  Influence functions for the coefficients.

- nlev:

  Number of levels in the outcome.

- px:

  Number of covariates.

## Details

This ensures that \\\sum_j P_j = 1\\. The model is fitted using the
`phreg` function by expanding the data into a long format with strata
for each category.

The coefficients represent the log-Relative-Risk (log-RR) relative to
the baseline group (the first level of the factor, which can be reset
using `relevel`).

Standard errors are computed based on the sandwich form: \$\$ D U^{-1}
\left( \sum U_i^2 \right) D U^{-1} \$\$ where \\U\\ is the score vector
and \\D\\ is the derivative matrix.

Influence functions (possibly robust) can be obtained via the
[`iid()`](https://kkholst.github.io/lava/reference/iid.html) function.
The response variable should be a factor.

Can also fit a cumulative odds model as a special case of
`interval_logitsurv_discrete`.

## Author

Thomas Scheike

## Examples

``` r
data(bmt)
bmt$id <- sample(200,408,replace=TRUE)
dfactor(bmt) <- cause1f~cause
drelevel(bmt,ref=3) <- cause3f~cause
dlevels(bmt)
#> cause1f #levels=:3 
#> [1] "0" "1" "2"
#> -----------------------------------------
#> cause3f #levels=:3 
#> [1] "2" "0" "1"
#> -----------------------------------------

mreg <- mlogit(cause1f~+1+cluster(id),bmt)
summary(mreg)
#> 
#>     n events
#>  1224    408
#> 
#>  1224 clusters
#> coefficients:
#>              Estimate       S.E.    dU^-1/2 P-value
#> Intercept2  0.0062305  0.1048289  0.1116297  0.9526
#> Intercept3 -0.6092657  0.1439164  0.1332076  0.0000
#> 
#> exp(coefficients):
#>            Estimate    2.5%  97.5%
#> Intercept2  1.00625 0.81936 1.2358
#> Intercept3  0.54375 0.41011 0.7209
#> 
head(iid(mreg))
#>      Intercept2    Intercept3
#> 1  2.484472e-02 -5.367040e-17
#> 2  6.211180e-03 -1.353524e-17
#> 3  6.172360e-03 -6.250000e-03
#> 4  1.242236e-02  1.149425e-02
#> 5  6.211180e-03  1.149425e-02
#> 6 -3.881988e-05 -6.250000e-03
dim(iid(mreg))
#> [1] 173   2

mreg <- mlogit(cause1f~tcell+platelet,bmt)
summary(mreg)
#> 
#>     n events
#>  1224    408
#> 
#>  1224 clusters
#> coefficients:
#>            Estimate     S.E.  dU^-1/2 P-value
#> Intercept2  0.25002  0.13906  0.13858  0.0722
#> tcell2     -0.28389  0.36431  0.36285  0.4358
#> platelet2  -0.68611  0.24797  0.24956  0.0057
#> Intercept3 -0.56565  0.16921  0.17078  0.0008
#> tcell3      0.50505  0.36481  0.36226  0.1662
#> platelet3  -0.35890  0.29130  0.28727  0.2179
#> 
#> exp(coefficients):
#>            Estimate    2.5%  97.5%
#> Intercept2  1.28406 0.97773 1.6864
#> tcell2      0.75285 0.36864 1.5375
#> platelet2   0.50353 0.30971 0.8187
#> Intercept3  0.56799 0.40767 0.7914
#> tcell3      1.65707 0.81062 3.3874
#> platelet3   0.69844 0.39462 1.2362
#> 
head(iid(mreg))
#>      Intercept2       tcell2    platelet2 Intercept3      tcell3   platelet3
#> 1 -0.0001355312 0.0009504519 8.857909e-05  0.0185731 -0.01020386 -0.01569271
#> 2 -0.0001355312 0.0009504519 8.857909e-05  0.0185731 -0.01020386 -0.01569271
#> 3 -0.0001355312 0.0009504519 8.857909e-05  0.0185731 -0.01020386 -0.01569271
#> 4 -0.0001355312 0.0009504519 8.857909e-05  0.0185731 -0.01020386 -0.01569271
#> 5 -0.0001355312 0.0009504519 8.857909e-05  0.0185731 -0.01020386 -0.01569271
#> 6 -0.0001355312 0.0009504519 8.857909e-05  0.0185731 -0.01020386 -0.01569271
dim(iid(mreg))
#> [1] 408   6

mreg3 <- mlogit(cause3f~tcell+platelet,bmt)
summary(mreg3)
#> 
#>     n events
#>  1224    408
#> 
#>  1224 clusters
#> coefficients:
#>            Estimate     S.E.  dU^-1/2 P-value
#> Intercept2  0.56565  0.16921  0.17078  0.0008
#> tcell2     -0.50505  0.36481  0.36226  0.1662
#> platelet2   0.35890  0.29130  0.28727  0.2179
#> Intercept3  0.81567  0.16346  0.16467  0.0000
#> tcell3     -0.78894  0.39244  0.38890  0.0444
#> platelet3  -0.32721  0.30423  0.30139  0.2821
#> 
#> exp(coefficients):
#>            Estimate    2.5%  97.5%
#> Intercept2  1.76059 1.26364 2.4530
#> tcell2      0.60347 0.29521 1.2336
#> platelet2   1.43175 0.80893 2.5341
#> Intercept3  2.26070 1.64100 3.1144
#> tcell3      0.45433 0.21053 0.9804
#> platelet3   0.72093 0.39713 1.3087
#> 

## inverse information standard errors 
lava::estimate(coef=mreg3$coef,vcov=mreg3$II)
#>            Estimate Std.Err    2.5%    97.5%   P-value
#> Intercept2   0.5657  0.1708  0.2309  0.90038 9.259e-04
#> tcell2      -0.5051  0.3623 -1.2151  0.20497 1.633e-01
#> platelet2    0.3589  0.2873 -0.2041  0.92194 2.115e-01
#> Intercept3   0.8157  0.1647  0.4929  1.13842 7.294e-07
#> tcell3      -0.7889  0.3889 -1.5512 -0.02672 4.249e-02
#> platelet3   -0.3272  0.3014 -0.9179  0.26350 2.776e-01

## predictions based on seen response or not 
## all probabilities
head(predict(mreg,response=FALSE))
#>           0         1         2
#> 1 0.3506254 0.4502227 0.1991519
#> 2 0.3506254 0.4502227 0.1991519
#> 3 0.3506254 0.4502227 0.1991519
#> 4 0.3506254 0.4502227 0.1991519
#> 5 0.3506254 0.4502227 0.1991519
#> 6 0.3506254 0.4502227 0.1991519
head(predict(mreg))
#>        pred         se      lower    upper
#> 1 0.1991519 0.06681508 0.06819674 0.330107
#> 2 0.1991519 0.06681508 0.06819674 0.330107
#> 3 0.1991519 0.06681508 0.06819674 0.330107
#> 4 0.1991519 0.06681508 0.06819674 0.330107
#> 5 0.1991519 0.06681508 0.06819674 0.330107
#> 6 0.1991519 0.06681508 0.06819674 0.330107
## using newdata 
newdata <- data.frame(tcell=c(1,1,1),platelet=c(0,1,1),cause1f=c("2","2","0"))
## only probability of seen response 
predict(mreg,newdata)
#>        pred         se      lower     upper
#> 1 0.3236700 0.13603308 0.05705011 0.5902900
#> 2 0.3065921 0.10966904 0.09164473 0.5215395
#> 3 0.4663875 0.07376361 0.32181348 0.6109615
## without response
predict(mreg,newdata,response=FALSE)
#>           0         1         2
#> 1 0.3438904 0.3324396 0.3236700
#> 2 0.4663875 0.2270204 0.3065921
#> 3 0.4663875 0.2270204 0.3065921
## given indexx of P(Y=j)
predict(mreg,newdata,Y=c(1,2,3))
#>        pred         se       lower     upper
#> 1 0.3438904 0.06916812  0.20832337 0.4794574
#> 2 0.2270204 0.12225406 -0.01259315 0.4666340
#> 3 0.3065921 0.10966904  0.09164473 0.5215395
##  reponse not given 
newdata <- data.frame(tcell=c(1,1,1),platelet=c(0,1,1))
predict(mreg,newdata)
#>           0         1         2
#> 1 0.3438904 0.3324396 0.3236700
#> 2 0.4663875 0.2270204 0.3065921
#> 3 0.4663875 0.2270204 0.3065921
```
