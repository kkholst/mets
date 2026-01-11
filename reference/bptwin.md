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
        data=droplevels(
          subset(twinstut, zyg%in%c("mz","dz") & tvparnr<5e3)
        ),
        id="tvparnr",zyg="zyg",DZ="dz",type="ae")
summary(b0)
#> 
#>             Estimate  Std.Err       Z   p-value    
#> (Intercept) -3.84320  0.67034 -5.7333 9.852e-09 ***
#> sexmale      0.80018  0.19370  4.1310 3.612e-05 ***
#> log(var(A))  1.12281  0.46483  2.4155   0.01571 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>  Total MZ/DZ Complete pairs MZ/DZ
#>  1483/2934   542/897             
#> 
#>                    Estimate 2.5%    97.5%  
#> A                  0.75451  0.58576 0.92326
#> E                  0.24549  0.07674 0.41424
#> MZ Tetrachoric Cor 0.75451  0.53102 0.87986
#> DZ Tetrachoric Cor 0.37726  0.28992 0.45836
#> 
#> MZ:
#>                      Estimate 2.5%     97.5%   
#> Concordance           0.01126  0.00629  0.02008
#> Casewise Concordance  0.39595  0.23173  0.58753
#> Marginal              0.02844  0.02205  0.03661
#> Rel.Recur.Risk       13.92073  7.23136 20.61010
#> log(OR)               3.59486  2.53322  4.65650
#> DZ:
#>                      Estimate 2.5%    97.5%  
#> Concordance          0.00378  0.00230 0.00621
#> Casewise Concordance 0.13291  0.09607 0.18105
#> Marginal             0.02844  0.02205 0.03661
#> Rel.Recur.Risk       4.67300  3.32186 6.02414
#> log(OR)              1.77247  1.40907 2.13587
#> 
#>                          Estimate 2.5%    97.5%  
#> Broad-sense heritability 0.75451  0.58576 0.92326
#> 
```
