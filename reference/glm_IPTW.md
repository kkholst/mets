# IPTW GLM, Inverse Probaibilty of Treatment Weighted GLM

Fits GLM model with treatment weights \$\$ w(A)= \sum_a I(A=a)/P(A=a\|X)
\$\$, computes standard errors via influence functions that are returned
as the IID argument. Propensity scores are fitted using either logistic
regression (glm) or the multinomial model (mlogit) when more than two
categories for treatment. The treatment needs to be a factor and is
identified on the rhs of the "treat.model".

## Usage

``` r
glm_IPTW(
  formula,
  data,
  treat.model = NULL,
  family = binomial(),
  id = NULL,
  weights = NULL,
  estpr = 1,
  pi0 = 0.5,
  ...
)
```

## Arguments

- formula:

  for glm

- data:

  data frame for risk averaging

- treat.model:

  propensity score model (binary or multinomial)

- family:

  of glm (logistic regression)

- id:

  cluster id for standard errors

- weights:

  may be given, and then uses weights\*w(A) as the weights

- estpr:

  to estimate propensity scores and get infuence function contribution
  to uncertainty

- pi0:

  fixed simple weights

- ...:

  arguments for glm call

## Details

Also works with cluster argument.

## Author

Thomas Scheike
