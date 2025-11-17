# Fast Cox PH regression

Fast Cox PH regression Robust variance is default variance with the
summary.

## Usage

``` r
phreg(formula, data, offset = NULL, weights = NULL, ...)
```

## Arguments

- formula:

  formula with 'Surv' outcome (see `coxph`)

- data:

  data frame

- offset:

  offsets for Cox model

- weights:

  weights for Cox score equations

- ...:

  Additional arguments to lower level funtions

## Details

influence functions (iid) will follow numerical order of given cluster
variable so ordering after \$id will give iid in order of data-set.

## Author

Klaus K. Holst, Thomas Scheike

## Examples

``` r
library(mets)
data(TRACE)
dcut(TRACE) <- ~.
out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4)+cluster(id),data=TRACE)
summary(out1)
#> 
#>     n events
#>  1878    958
#> coeffients:
#>     Estimate     S.E.  dU^-1/2 P-value
#> vf  0.452306 0.136473 0.111038   9e-04
#> chf 0.931822 0.074226 0.074650   0e+00
#> 
#> exp(coeffients):
#>     Estimate   2.5%  97.5%
#> vf    1.5719 1.2030 2.0540
#> chf   2.5391 2.1954 2.9367
#> 

par(mfrow=c(1,2))
plot(out1)

## computing robust variance for baseline
rob1 <- robust.phreg(out1)
plot(rob1,se=TRUE,robust=TRUE)


## iid decomposition, with scaled influence functions
## for regression parameters
head(iid(out1))
#>               vf           chf
#> 3  -0.0004533127 -0.0027020566
#> 7   0.0089952016  0.0001774394
#> 13  0.0024358693 -0.0016241766
#> 15 -0.0006893327  0.0016569786
#> 17 -0.0018064125  0.0006870047
#> 22  0.0001540691  0.0008916881
## making iid decomposition of baseline at a specific time-point
Aiiid <- iid(out1,time=30)
head(Aiiid)
#>          strata0       strata1       strata2       strata3
#> 3   0.0015058392  0.0008448281  0.0018373121  4.104713e-04
#> 7  -0.0001202832 -0.0004714924 -0.0002518205 -1.588681e-04
#> 13 -0.0022290969  0.0003821440  0.0003022758  2.068804e-04
#> 15 -0.0008476467 -0.0004731874  0.0006257034 -2.374801e-04
#> 17  0.0004864903 -0.0001256304 -0.0001095240 -7.608883e-05
#> 22 -0.0004972827 -0.0002790038 -0.0006202510 -1.355228e-04
## both iid decompositions
dd <- iidBaseline(out1,time=30)
head(dd$beta.iid)
#>               [,1]          [,2]
#> [1,] -0.0004533127 -0.0027020566
#> [2,]  0.0089952016  0.0001774394
#> [3,]  0.0024358693 -0.0016241766
#> [4,] -0.0006893327  0.0016569786
#> [5,] -0.0018064125  0.0006870047
#> [6,]  0.0001540691  0.0008916881
head(dd$base.iid)
#>          strata0       strata1       strata2       strata3
#> 3   0.0015058392  0.0008448281  0.0018373121  4.104713e-04
#> 7  -0.0001202832 -0.0004714924 -0.0002518205 -1.588681e-04
#> 13 -0.0022290969  0.0003821440  0.0003022758  2.068804e-04
#> 15 -0.0008476467 -0.0004731874  0.0006257034 -2.374801e-04
#> 17  0.0004864903 -0.0001256304 -0.0001095240 -7.608883e-05
#> 22 -0.0004972827 -0.0002790038 -0.0006202510 -1.355228e-04
```
