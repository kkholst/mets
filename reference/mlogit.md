# Multinomial regression based on phreg regression

Fits multinomial regression model \$\$ P_i = \frac{ \exp( X^\beta_i ) }{
\sum\_{j=1}^K \exp( X^\beta_j ) }\$\$ for \$\$i=1,..,K\$\$ where
\$\$\beta_1 = 0\$\$, such that \$\$\sum_j P_j = 1\$\$ using phreg
function. Thefore the ratio \$\$\frac{P_i}{P_1} = \exp( X^\beta_i )\$\$

## Usage

``` r
mlogit(formula, data, offset = NULL, weights = NULL, fix.X = FALSE, ...)
```

## Arguments

- formula:

  formula with outcome (see `coxph`)

- data:

  data frame

- offset:

  offsets for partial likelihood

- weights:

  for score equations

- fix.X:

  to have same coefficients for all categories

- ...:

  Additional arguments to lower level funtions

## Details

Coefficients give log-Relative-Risk relative to baseline group (first
level of factor, so that it can reset by relevel command). Standard
errors computed based on sandwhich form \$\$ DU^-1 \sum U_i^2 DU^-1\$\$.

Can also get influence functions (possibly robust) via iid() function,
response should be a factor.

Can fit cumulative odds model as a special case of
interval.logitsurv.discrete

## Author

Thomas Scheike

## Examples

``` r
library(mets)
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
#> coeffients:
#>                   Estimate       S.E.    dU^-1/2 P-value
#> XX(Intercept).2  0.0062305  0.1187375  0.1116297  0.9582
#> XX(Intercept).3 -0.6092657  0.1414469  0.1332076  0.0000
#> 
#> exp(coeffients):
#>                 Estimate    2.5%  97.5%
#> XX(Intercept).2  1.00625 0.79733 1.2699
#> XX(Intercept).3  0.54375 0.41210 0.7175
#> 
head(iid(mreg))
#>   XX(Intercept).2 XX(Intercept).3
#> 1   -6.288820e-03   -1.005747e-03
#> 2    1.692758e-17    1.149425e-02
#> 3    6.172360e-03   -6.250000e-03
#> 4    6.211180e-03    2.298851e-02
#> 5    1.242236e-02   -2.801584e-17
#> 6    3.246739e-17    2.298851e-02
dim(iid(mreg))
#> [1] 174   2

mreg <- mlogit(cause1f~tcell+platelet,bmt)
summary(mreg)
#> 
#>     n events
#>  1224    408
#> 
#>  1224 clusters
#> coeffients:
#>                 Estimate     S.E.  dU^-1/2 P-value
#> XX(Intercept).2  0.25002  0.13906  0.13858  0.0722
#> XXtcell.2       -0.28389  0.36431  0.36285  0.4358
#> XXplatelet.2    -0.68611  0.24797  0.24956  0.0057
#> XX(Intercept).3 -0.56565  0.16921  0.17078  0.0008
#> XXtcell.3        0.50505  0.36481  0.36226  0.1662
#> XXplatelet.3    -0.35890  0.29130  0.28727  0.2179
#> 
#> exp(coeffients):
#>                 Estimate    2.5%  97.5%
#> XX(Intercept).2  1.28406 0.97773 1.6864
#> XXtcell.2        0.75285 0.36864 1.5375
#> XXplatelet.2     0.50353 0.30971 0.8187
#> XX(Intercept).3  0.56799 0.40767 0.7914
#> XXtcell.3        1.65707 0.81062 3.3874
#> XXplatelet.3     0.69844 0.39462 1.2362
#> 
head(iid(mreg))
#>   XX(Intercept).2    XXtcell.2 XXplatelet.2 XX(Intercept).3   XXtcell.3
#> 1   -0.0001355312 0.0009504519 8.857909e-05       0.0185731 -0.01020386
#> 2   -0.0001355312 0.0009504519 8.857909e-05       0.0185731 -0.01020386
#> 3   -0.0001355312 0.0009504519 8.857909e-05       0.0185731 -0.01020386
#> 4   -0.0001355312 0.0009504519 8.857909e-05       0.0185731 -0.01020386
#> 5   -0.0001355312 0.0009504519 8.857909e-05       0.0185731 -0.01020386
#> 6   -0.0001355312 0.0009504519 8.857909e-05       0.0185731 -0.01020386
#>   XXplatelet.3
#> 1  -0.01569271
#> 2  -0.01569271
#> 3  -0.01569271
#> 4  -0.01569271
#> 5  -0.01569271
#> 6  -0.01569271
dim(iid(mreg))
#> [1] 408   6

mreg3 <- mlogit(cause3f~tcell+platelet,bmt)
summary(mreg3)
#> 
#>     n events
#>  1224    408
#> 
#>  1224 clusters
#> coeffients:
#>                 Estimate     S.E.  dU^-1/2 P-value
#> XX(Intercept).2  0.56565  0.16921  0.17078  0.0008
#> XXtcell.2       -0.50505  0.36481  0.36226  0.1662
#> XXplatelet.2     0.35890  0.29130  0.28727  0.2179
#> XX(Intercept).3  0.81567  0.16346  0.16467  0.0000
#> XXtcell.3       -0.78894  0.39244  0.38890  0.0444
#> XXplatelet.3    -0.32721  0.30423  0.30139  0.2821
#> 
#> exp(coeffients):
#>                 Estimate    2.5%  97.5%
#> XX(Intercept).2  1.76059 1.26364 2.4530
#> XXtcell.2        0.60347 0.29521 1.2336
#> XXplatelet.2     1.43175 0.80893 2.5341
#> XX(Intercept).3  2.26070 1.64100 3.1144
#> XXtcell.3        0.45433 0.21053 0.9804
#> XXplatelet.3     0.72093 0.39713 1.3087
#> 

## inverse information standard errors 
lava::estimate(coef=mreg3$coef,vcov=mreg3$II)
#>                 Estimate Std.Err    2.5%    97.5%   P-value
#> XX(Intercept).2   0.5657  0.1708  0.2309  0.90038 9.259e-04
#> XXtcell.2        -0.5051  0.3623 -1.2151  0.20497 1.633e-01
#> XXplatelet.2      0.3589  0.2873 -0.2041  0.92194 2.115e-01
#> XX(Intercept).3   0.8157  0.1647  0.4929  1.13842 7.294e-07
#> XXtcell.3        -0.7889  0.3889 -1.5512 -0.02672 4.249e-02
#> XXplatelet.3     -0.3272  0.3014 -0.9179  0.26350 2.776e-01

## predictions based on seen response or not 
newdata <- data.frame(tcell=c(1,1,1),platelet=c(0,1,1),cause1f=c("2","1","0"))
## all probabilities
predict(mreg,newdata,response=FALSE)
#>           0         1         2
#> 1 0.3438904 0.3324396 0.3236700
#> 2 0.4663875 0.2270204 0.3065921
#> 3 0.4663875 0.2270204 0.3065921
## only probability of seen response 
predict(mreg,newdata)
#>        pred         se       lower     upper
#> 1 0.3236700 0.13603308  0.05705011 0.5902900
#> 2 0.2270204 0.12225406 -0.01259315 0.4666340
#> 3 0.4663875 0.07376361  0.32181348 0.6109615
```
