# Classic twin model for quantitative traits

Fits a classical twin model for quantitative traits.

## Usage

``` r
twinlm(
  formula,
  data,
  id,
  zyg,
  DZ,
  group = NULL,
  group.equal = FALSE,
  strata = NULL,
  weights = NULL,
  type = c("ace"),
  twinnum = "twinnum",
  binary = FALSE,
  ordinal = 0,
  keep = weights,
  estimator = NULL,
  constrain = TRUE,
  control = list(),
  messages = 1,
  ...
)
```

## Arguments

- formula:

  Formula specifying effects of covariates on the response

- data:

  `data.frame` with one observation pr row. In addition a column with
  the zygosity (DZ or MZ given as a factor) of each individual much be
  specified as well as a twin id variable giving a unique pair of
  numbers/factors to each twin pair

- id:

  The name of the column in the dataset containing the twin-id variable.

- zyg:

  The name of the column in the dataset containing the zygosity variable

- DZ:

  Character defining the level in the zyg variable corresponding to the
  dyzogitic twins. If this argument is missing, the reference level
  (i.e. the first level) will be interpreted as the dyzogitic twins

- group:

  Optional. Variable name defining group for interaction analysis (e.g.,
  gender)

- group.equal:

  If TRUE marginals of groups are asummed to be the same

- strata:

  Strata variable name

- weights:

  Weights matrix if needed by the chosen estimator. For use with Inverse
  Probability Weights

- type:

  Character defining the type of analysis to be performed. Can be a
  subset of "aced" (additive genetic factors, common environmental
  factors, unique environmental factors, dominant genetic factors).
  Other choices are:

  - "0" (or "sat"): Saturated model where twin 1 and twin 2 within each
    twin pair may have a different marginal distribution.

  - "1" (or "flex","zyg"): Within twin pairs the marginal distribution
    is the same, but the marginal distribution may differ between MZ and
    DZ twins. A free correlation structure within MZ and DZ twins.

  - "2" (or "u", "eqmarg"): All individuals have the same marginals but
    a free correlation structure within MZ and DZ twins.

  The default value is an additive polygenic model `type="ace"`.

- twinnum:

  The name of the column in the dataset numbering the twins (1,2). If it
  does not exist in `data` it will automatically be created.

- binary:

  If `TRUE` a liability model is fitted. Note that if the
  right-hand-side of the formula is a factor, character vector, og
  logical variable, then the liability model is automatically chosen
  (wrapper of the `bptwin` function).

- ordinal:

  If non-zero (number of bins) a liability model is fitted.

- keep:

  Vector of variables from `data` that are not specified in `formula`,
  to be added to data.frame of the SEM

- estimator:

  Choice of estimator/model

- constrain:

  Development argument

- control:

  Control argument parsed on to the optimization routine

- messages:

  Control amount of messages shown

- ...:

  Additional arguments parsed on to lower-level functions

## Value

Returns an object of class `twinlm`.

## See also

[`bptwin`](http://kkholst.github.io/mets/reference/bptwin.md),
[`twinlm.time`](http://kkholst.github.io/mets/reference/bptwin.md),
`twinlm.strata`,
[`twinsim`](http://kkholst.github.io/mets/reference/twinsim.md)

## Author

Klaus K. Holst

## Examples

``` r
## Simulate data
set.seed(1)
d <- twinsim(1000,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
## E(y|z1,z2) = z1 - z2. var(A) = var(C) = var(E) = 1

## E.g to fit the data to an ACE-model without any confounders we simply write
ace <- twinlm(y ~ 1, data=d, DZ="DZ", zyg="zyg", id="id")
ace
#>        Estimate Std. Error Z value  Pr(>|z|)
#> y     -0.019439   0.041817 -0.4649     0.642
#> sd(A)  0.902004   0.203739  4.4273 9.544e-06
#> sd(C)  1.137025   0.132852  8.5586 < 2.2e-16
#> sd(E)  1.728992   0.037408 46.2194 < 2.2e-16
#> 
#> MZ-pairs DZ-pairs 
#>     1000     1000 
#> 
#> Variance decomposition:
#>   Estimate 2.5%    97.5%  
#> A 0.15966  0.01867 0.30065
#> C 0.25370  0.13920 0.36820
#> E 0.58664  0.53677 0.63650
#> 
#> 
#>                          Estimate 2.5%    97.5%  
#> Broad-sense heritability 0.15966  0.01867 0.30065
#> 
#>                        Estimate 2.5%    97.5%  
#> Correlation within MZ: 0.41336  0.36229 0.46196
#> Correlation within DZ: 0.33353  0.27933 0.38561
#> 
#> 'log Lik.' -8779.953 (df=4)
#> AIC: 17567.91 
#> BIC: 17590.31 
## An AE-model could be fitted as
ae <- twinlm(y ~ 1, data=d, DZ="DZ", zyg="zyg", id="id", type="ae")
## LRT:
lava::compare(ae,ace)
#> 
#>  - Likelihood ratio test -
#> 
#> data:  
#> chisq = 17.207, df = 1, p-value = 3.353e-05
#> sample estimates:
#> log likelihood (model 1) log likelihood (model 2) 
#>                -8788.556                -8779.953 
#> 
## AIC
AIC(ae)-AIC(ace)
#> [1] 15.20656
## To adjust for the covariates we simply alter the formula statement
ace2 <- twinlm(y ~ x1+x2, data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
## Summary/GOF
summary(ace2)
#>        Estimate Std. Error  Z value Pr(>|z|)
#> y     -0.026049   0.034844  -0.7476   0.4547
#> sd(A)  1.066060   0.072890  14.6256   <2e-16
#> sd(C)  0.980740   0.073569  13.3309   <2e-16
#> sd(E)  0.979980   0.021887  44.7736   <2e-16
#> y~x1   1.006963   0.021900  45.9807   <2e-16
#> y~x2  -0.993802   0.021962 -45.2512   <2e-16
#> 
#> MZ-pairs DZ-pairs 
#>     1000     1000 
#> 
#> Variance decomposition:
#>   Estimate 2.5%    97.5%  
#> A 0.37156  0.27300 0.47012
#> C 0.31446  0.22643 0.40250
#> E 0.31398  0.28381 0.34414
#> 
#> 
#>                          Estimate 2.5%    97.5%  
#> Broad-sense heritability 0.37156  0.27300 0.47012
#> 
#>                        Estimate 2.5%    97.5%  
#> Correlation within MZ: 0.68602  0.65467 0.71502
#> Correlation within DZ: 0.50024  0.45538 0.54257
#> 
#> 'log Lik.' -7449.697 (df=6)
#> AIC: 14911.39 
#> BIC: 14945 
 ## Reduce Ex.Timings
## An interaction could be analyzed as:
ace3 <- twinlm(y ~ x1+x2 + x1:I(x2<0), data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
ace3
#>                     Estimate Std. Error  Z value Pr(>|z|)
#> y                  -0.026089   0.034847  -0.7487   0.4541
#> sd(A)               1.065529   0.072975  14.6012   <2e-16
#> sd(C)               0.981354   0.073598  13.3340   <2e-16
#> sd(E)               0.980018   0.021889  44.7711   <2e-16
#> y~x1                1.010637   0.030279  33.3778   <2e-16
#> y~x2               -0.993859   0.021964 -45.2498   <2e-16
#> y~x1:I(x2 < 0)TRUE -0.007626   0.043403  -0.1757   0.8605
#> 
#> MZ-pairs DZ-pairs 
#>     1000     1000 
#> 
#> Variance decomposition:
#>   Estimate 2.5%    97.5%  
#> A 0.37117  0.27253 0.46981
#> C 0.31484  0.22673 0.40296
#> E 0.31399  0.28382 0.34415
#> 
#> 
#>                          Estimate 2.5%    97.5%  
#> Broad-sense heritability 0.37117  0.27253 0.46981
#> 
#>                        Estimate 2.5%    97.5%  
#> Correlation within MZ: 0.68601  0.65466 0.71501
#> Correlation within DZ: 0.50043  0.45553 0.54279
#> 
#> 'log Lik.' -7449.682 (df=7)
#> AIC: 14913.36 
#> BIC: 14952.57 
## Categorical variables are also supported 
d2 <- transform(d,x2cat=cut(x2,3,labels=c("Low","Med","High")))
ace4 <- twinlm(y ~ x1+x2cat, data=d2, DZ="DZ", zyg="zyg", id="id", type="ace")
```
