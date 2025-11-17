# IPTW Cox, Inverse Probaibilty of Treatment Weighted Cox regression

Fits Cox model with treatment weights \$\$ w(A)= \sum_a
I(A=a)/\pi(a\|X)\$\$, where \$\$\pi(a\|X)=P(A=a\|X)\$\$. Computes
standard errors via influence functions that are returned as the IID
argument. Propensity scores are fitted using either logistic regression
(glm) or the multinomial model (mlogit) when there are than treatment
categories. The treatment needs to be a factor and is identified on the
rhs of the "treat.model". Recurrent events can be considered with
start,stop structure and then cluster(id) must be specified. Robust
standard errors are computed in all cases.

## Usage

``` r
phreg_IPTW(
  formula,
  data,
  treat.model = NULL,
  treat.var = NULL,
  weights = NULL,
  estpr = 1,
  pi0 = 0.5,
  se.cluster = NULL,
  ...
)
```

## Arguments

- formula:

  for phreg

- data:

  data frame for risk averaging

- treat.model:

  propensity score model (binary or multinomial)

- treat.var:

  a 1/0 variable that indicates when treatment is given and the
  propensity score is computed

- weights:

  may be given, and then uses weights\*w(A) as the weights

- estpr:

  (=1, default) to estimate propensity scores and get infuence function
  contribution to uncertainty

- pi0:

  fixed simple weights

- se.cluster:

  to compute GEE type standard errors when additional cluster structure
  is present

- ...:

  arguments for phreg call

## Details

Time-dependent propensity score weights can also be computed when
treat.var is used, it must be 1 at the time of first (A_0) and 2nd
treatment (A_1), then uses weights \$\$w_0(A_0) \* w_1(A_1)^{t\>T_r}\$\$
where \$\$T_r\$\$ is time of 2nd randomization.

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data <- mets:::simLT(0.7,100,beta=0.3,betac=0,ce=1,betao=0.3)
dfactor(data) <- Z.f~Z
out <- phreg_IPTW(Surv(time,status)~Z.f,data=data,treat.model=Z.f~X)
summary(out)
#> 
#>    n events
#>  100     45
#> 
#>  100 clusters
#> coeffients:
#>      Estimate     S.E.  dU^-1/2 P-value
#> Z.f1 -0.47160  0.25676  0.20863  0.0662
#> 
#> exp(coeffients):
#>      Estimate    2.5%  97.5%
#> Z.f1  0.62400 0.37725 1.0321
#> 
```
