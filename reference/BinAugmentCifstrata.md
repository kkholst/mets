# Augmentation for Binomial regression based on stratified NPMLE Cif (Aalen-Johansen)

Computes the augmentation term for each individual as well as the sum
\$\$ A = \int_0^t H(u,X) \frac{1}{S^\*(u,s)} \frac{1}{G_c(u)} dM_c(u)
\$\$ with \$\$ H(u,X) = F_1^\*(t,s) - F_1^\*(u,s) \$\$ using a KM for
\$\$G_c(t)\$\$ and a working model for cumulative baseline related to
\$\$F_1^\*(t,s)\$\$ and \$\$s\$\$ is strata, \$\$S^\*(t,s) = 1 -
F_1^\*(t,s) - F_2^\*(t,s)\$\$.

## Usage

``` r
BinAugmentCifstrata(
  formula,
  data = data,
  cause = 1,
  cens.code = 0,
  km = TRUE,
  time = NULL,
  weights = NULL,
  offset = NULL,
  ...
)
```

## Arguments

- formula:

  formula with 'Event', strata model for CIF given by strata, and
  strataC specifies censoring strata

- data:

  data frame

- cause:

  of interest

- cens.code:

  code of censoring

- km:

  to use Kaplan-Meier

- time:

  of interest

- weights:

  weights for estimating equations

- offset:

  offsets for logistic regression

- ...:

  Additional arguments to binreg function.

## Details

Standard errors computed under assumption of correct \$\$G_c(s)\$\$
model.

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt)
dcut(bmt,breaks=2) <- ~age 
out1<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
      strata(platelet,agecat.2),data=bmt,cause=1,time=40)
summary(out1)
#>    n events
#>  408    157
#> 
#>  408 clusters
#> coeffients:
#>                      Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)          -0.51295  0.17090 -0.84791 -0.17799  0.0027
#> platelet             -0.63011  0.23585 -1.09237 -0.16785  0.0075
#> agecat.2(0.203,1.94]  0.55926  0.21211  0.14353  0.97500  0.0084
#> 
#> exp(coeffients):
#>                      Estimate    2.5%  97.5%
#> (Intercept)           0.59873 0.42831 0.8370
#> platelet              0.53253 0.33542 0.8455
#> agecat.2(0.203,1.94]  1.74938 1.15434 2.6512
#> 
#> 

out2<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
    strata(platelet,agecat.2)+strataC(platelet),data=bmt,cause=1,time=40)
summary(out2)
#>    n events
#>  408    157
#> 
#>  408 clusters
#> coeffients:
#>                      Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)          -0.51346  0.17109 -0.84879 -0.17814  0.0027
#> platelet             -0.63636  0.23653 -1.09996 -0.17276  0.0071
#> agecat.2(0.203,1.94]  0.56280  0.21229  0.14672  0.97889  0.0080
#> 
#> exp(coeffients):
#>                      Estimate    2.5%  97.5%
#> (Intercept)           0.59842 0.42793 0.8368
#> platelet              0.52922 0.33288 0.8413
#> agecat.2(0.203,1.94]  1.75559 1.15803 2.6615
#> 
#> 
```
