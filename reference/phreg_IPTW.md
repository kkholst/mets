# IPTW Cox Regression (Inverse Probability of Treatment Weighted)

Fits a Cox model with treatment weights \$\$w(A) = \sum_a
I(A=a)/\pi(a\|X)\$\$, where \$\$\pi(a\|X) = P(A=a\|X)\$\$.

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

  Formula for `phreg`.

- data:

  Data frame for risk averaging.

- treat.model:

  Propensity score model (binary or multinomial).

- treat.var:

  A 1/0 variable indicating when treatment is given (for time-dependent
  weights).

- weights:

  Optional weights to multiply with the IPTW weights.

- estpr:

  (=1, default) to estimate propensity scores and include their
  uncertainty in the influence function.

- pi0:

  Fixed simple weights (if `estpr=0`).

- se.cluster:

  To compute GEE-type standard errors when additional cluster structure
  is present.

- ...:

  Arguments for `phreg` call.

## Value

An object of class `"phreg"` with additional IPTW components:

- IID:

  Influence functions including propensity score uncertainty.

- iptw:

  IPTW weights used.

- naive.var:

  Naive variance ignoring propensity score uncertainty.

## Details

Standard errors are computed via influence functions that are returned
as the IID argument. Propensity scores are fitted using either logistic
regression (`glm`) or the multinomial model (`mlogit`) when there are
more than two treatment categories.

The treatment variable must be a factor and is identified on the RHS of
the `treat.model`. Recurrent events can be considered with a start-stop
structure, requiring `cluster(id)`. Robust standard errors are computed
in all cases.

Time-dependent propensity score weights can be computed when `treat.var`
is used. This weight be 1 at the time of first (A_0) and 2nd treatment
(A_1), then uses weights \$\$w_0(A_0) \* w_1(A_1)^{t\>T_r}\$\$ where
\$\$T_r\$\$ is time of 2nd randomization. The weights are constructed
using a `glm` or `mlogit` model based on the data where `treat.var=1`.
The propensity score can be constructed for any number of treatments in
a similar manner.

## Author

Thomas Scheike

## Examples

``` r
data <- mets:::simLT(0.7, 100, beta = 0.3, betac = 0, ce = 1, betao = 0.3)
dfactor(data) <- Z.f ~ Z
out <- phreg_IPTW(Surv(time, status) ~ Z.f, data = data, treat.model = Z.f ~ X)
summary(out)
#> 
#>    n events
#>  100     38
#> 
#>  100 clusters
#> coefficients:
#>      Estimate     S.E.  dU^-1/2 P-value
#> Z.f1 -0.84702  0.27941  0.23943  0.0024
#> 
#> exp(coefficients):
#>      Estimate    2.5%  97.5%
#> Z.f1  0.42869 0.24792 0.7413
#> 
```
