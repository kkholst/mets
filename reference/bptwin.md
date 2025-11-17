# Liability model for twin data

Liability-threshold model for twin data

## Usage

``` r
bptwin(
  x,
  data,
  id,
  zyg,
  DZ,
  group = NULL,
  num = NULL,
  weights = NULL,
  weights.fun = function(x) ifelse(any(x <= 0), 0, max(x)),
  strata = NULL,
  messages = 1,
  control = list(trace = 0),
  type = "ace",
  eqmean = TRUE,
  pairs.only = FALSE,
  samecens = TRUE,
  allmarg = samecens & !is.null(weights),
  stderr = TRUE,
  robustvar = TRUE,
  p,
  indiv = FALSE,
  constrain,
  varlink,
  ...
)
```

## Arguments

- x:

  Formula specifying effects of covariates on the response.

- data:

  `data.frame` with one observation pr row. In addition a column with
  the zygosity (DZ or MZ given as a factor) of each individual much be
  specified as well as a twin id variable giving a unique pair of
  numbers/factors to each twin pair.

- id:

  The name of the column in the dataset containing the twin-id variable.

- zyg:

  The name of the column in the dataset containing the zygosity
  variable.

- DZ:

  Character defining the level in the zyg variable corresponding to the
  dyzogitic twins.

- group:

  Optional. Variable name defining group for interaction analysis (e.g.,
  gender)

- num:

  Optional twin number variable

- weights:

  Weight matrix if needed by the chosen estimator (IPCW)

- weights.fun:

  Function defining a single weight each individual/cluster

- strata:

  Strata

- messages:

  Control amount of messages shown

- control:

  Control argument parsed on to the optimization routine. Starting
  values may be parsed as '`start`'.

- type:

  Character defining the type of analysis to be performed. Should be a
  subset of "acde" (additive genetic factors, common environmental
  factors, dominant genetic factors, unique environmental factors).

- eqmean:

  Equal means (with type="cor")?

- pairs.only:

  Include complete pairs only?

- samecens:

  Same censoring

- allmarg:

  Should all marginal terms be included

- stderr:

  Should standard errors be calculated?

- robustvar:

  If TRUE robust (sandwich) variance estimates of the variance are used

- p:

  Parameter vector p in which to evaluate log-Likelihood and score
  function

- indiv:

  If TRUE the score and log-Likelihood contribution of each twin-pair

- constrain:

  Development argument

- varlink:

  Link function for variance parameters

- ...:

  Additional arguments to lower level functions

## See also

[`twinlm`](http://kkholst.github.io/mets/reference/twinlm.md),
`twinlm.time`,
[`twinlm.strata`](http://kkholst.github.io/mets/reference/twinlm.md),
[`twinsim`](http://kkholst.github.io/mets/reference/twinsim.md)

## Author

Klaus K. Holst

## Examples

``` r
data(twinstut)
b0 <- bptwin(stutter~sex,
             data=droplevels(subset(twinstut,zyg%in%c("mz","dz"))),
             id="tvparnr",zyg="zyg",DZ="dz",type="ae")
#> Warning: setting environment(<primitive function>) is not possible and trying it is deprecated
#> Warning: setting environment(<primitive function>) is not possible and trying it is deprecated
#> Warning: setting environment(<primitive function>) is not possible and trying it is deprecated
summary(b0)
#> 
#>             Estimate  Std.Err        Z   p-value    
#> (Intercept) -3.70371  0.24449 -15.1485 < 2.2e-16 ***
#> sexmale      0.83310  0.08255  10.0920 < 2.2e-16 ***
#> log(var(A))  1.18278  0.17179   6.8851 5.774e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>  Total MZ/DZ Complete pairs MZ/DZ
#>  8777/12511  3255/4058           
#> 
#>                    Estimate 2.5%    97.5%  
#> A                  0.76545  0.70500 0.82590
#> E                  0.23455  0.17410 0.29500
#> MZ Tetrachoric Cor 0.76545  0.69793 0.81948
#> DZ Tetrachoric Cor 0.38272  0.35210 0.41253
#> 
#> MZ:
#>                      Estimate 2.5%     97.5%   
#> Concordance           0.01560  0.01273  0.01912
#> Casewise Concordance  0.42830  0.36248  0.49677
#> Marginal              0.03643  0.03294  0.04027
#> Rel.Recur.Risk       11.75741  9.77237 13.74246
#> log(OR)               3.52382  3.13466  3.91298
#> DZ:
#>                      Estimate 2.5%    97.5%  
#> Concordance          0.00558  0.00465 0.00670
#> Casewise Concordance 0.15327  0.13749 0.17050
#> Marginal             0.03643  0.03294 0.04027
#> Rel.Recur.Risk       4.20744  3.78588 4.62900
#> log(OR)              1.69996  1.57262 1.82730
#> 
#>                          Estimate 2.5%    97.5%  
#> Broad-sense heritability 0.76545  0.70500 0.82590
#> 
```
