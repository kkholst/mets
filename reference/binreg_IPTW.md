# IPTW logistic regression, Inverse Probabibilty of Treatment Weighted binreg

Fits logistic regression model with treatment weights \$\$ w(A)= \sum_a
I(A=a)/P(A=a\|X) \$\$, computes standard errors via influence functions
that are returned as the IID argument. Propensity scores are fitted
using either logistic regression (glm) or the multinomial model (mlogit)
when more than two categories for treatment. The treatment needs to be a
factor and is identified on the rhs of the "treat.model". Can handle
right censored binreg type estimating equations with IPTW weights.

## Usage

``` r
binreg_IPTW(
  formula,
  data,
  treat.model = NULL,
  weights = NULL,
  estpr = 1,
  pi0 = 0.5,
  ...
)
```

## Arguments

- formula:

  for binreg

- data:

  data-frame for estimation

- treat.model:

  propensity score model (binary or multinomial)

- weights:

  may be given, and then uses weights\*w(A) as the weights

- estpr:

  to estimate propensity scores and get infuence function contribution
  to uncertainty

- pi0:

  fixed simple weights

- ...:

  arguments for binreg call

## Details

Also works with cluster argument.

## Author

Thomas Scheike

## Examples

``` r

data(bmt)
dfactor(bmt) <- platelet.f~platelet
## logistic modelling of cumulative incidence 
gg <- binreg_IPTW(Event(time,cause)~platelet.f+age,bmt,
         treat.model=platelet.f~age,time=30)
summary(gg)
#>    n   events
#>  408 157.4166
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.270287  0.126209 -0.517652 -0.022922  0.0322
#> platelet.f1 -0.715767  0.241921 -1.189924 -0.241611  0.0031
#> age          0.399177  0.110405  0.182787  0.615568  0.0003
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.76316 0.59592 0.9773
#> platelet.f1  0.48882 0.30424 0.7854
#> age          1.49060 1.20056 1.8507
#> 
#> 
head(iid(gg))
#>           [,1]        [,2]          [,3]
#> 1 -0.006674904 0.006826456 -0.0004443729
#> 2 -0.006969323 0.007595318 -0.0025144319
#> 3 -0.007438561 0.009032129 -0.0069378873
#> 4 -0.007423086 0.008978572 -0.0067597973
#> 5 -0.006498370 0.006397056  0.0006225813
#> 6 -0.007587962 0.009582294 -0.0088280331

## logistic modelling  
gg <- binreg_IPTW(tcell~platelet.f+age,bmt,
         treat.model=platelet.f~age)
summary(gg)
#>    n events
#>  408     54
#> 
#>  408 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept) -2.43527  0.22211 -2.87058 -1.99995   0e+00
#> platelet.f1  1.13878  0.30683  0.53740  1.74016   2e-04
#> age          0.59404  0.16633  0.26803  0.92004   4e-04
#> 
#> exp(coeffients):
#>             Estimate     2.5%  97.5%
#> (Intercept) 0.087575 0.056666 0.1353
#> platelet.f1 3.122963 1.711558 5.6983
#> age         1.811283 1.307384 2.5094
#> 
#> 
```
