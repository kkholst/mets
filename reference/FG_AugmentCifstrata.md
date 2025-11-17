# Augmentation for Fine-Gray model based on stratified NPMLE Cif (Aalen-Johansen)

Computes the augmentation term for each individual as well as the sum
\$\$ A(\beta) = \int H(t,X,\beta) \frac{F_2^\*(t,s)}{S^\*(t,s)}
\frac{1}{G_c(t)} dM_c \$\$ with \$\$ H(t,X,\beta) = \int_t^\infty (X -
E(\beta,t) ) G_c(t) d\Lambda_1^\*i(t,s) \$\$ using a KM for
\$\$G_c(t)\$\$ and a working model for cumulative baseline related to
\$\$F_1^\*(t,s)\$\$ and \$\$s\$\$ is strata, \$\$S^\*(t,s) = 1 -
F_1^\*(t,s) - F_2^\*(t,s)\$\$, and \$\$E(\beta^p,t)\$\$ is given.
Assumes that no strata for baseline of ine-Gay model that is augmented.

## Usage

``` r
FG_AugmentCifstrata(
  formula,
  data = data,
  E = NULL,
  cause = NULL,
  cens.code = 0,
  km = TRUE,
  case.weights = NULL,
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

- E:

  from FG-model

- cause:

  of interest

- cens.code:

  code of censoring

- km:

  to use Kaplan-Meier

- case.weights:

  weights for FG score equations (that follow dN_1)

- weights:

  weights for FG score equations

- offset:

  offsets for FG model

- ...:

  Additional arguments to lower level funtions

## Details

After a couple of iterations we end up with a solution of \$\$ \int (X -
E(\beta) ) Y_1(t) w(t) dM_1 + A(\beta) \$\$ the augmented FG-score.

Standard errors computed under assumption of correct \$\$G_c\$\$ model.

## Author

Thomas Scheike

## Examples

``` r
set.seed(100)
rho1 <- 0.2; rho2 <- 10
n <- 400
beta=c(0.0,-0.1,-0.5,0.3)
dats <- simul.cifs(n,rho1,rho2,beta,rc=0.2)
dtable(dats,~status)
#> 
#> status
#>   0   1   2 
#>  14  54 332 
#> 
dsort(dats) <- ~time
fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
summary(fg)
#> 
#>    n events
#>  400     54
#> 
#>  400 clusters
#> coeffients:
#>     Estimate      S.E.   dU^-1/2 P-value
#> Z1  0.028262  0.135312  0.136188  0.8346
#> Z2 -0.149224  0.271222  0.272379  0.5822
#> 
#> exp(coeffients):
#>    Estimate    2.5%  97.5%
#> Z1  1.02866 0.78903 1.3411
#> Z2  0.86138 0.50621 1.4657
#> 

fgaugS <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fg$E)
summary(fgaugS)
#> 
#>    n events
#>  400     54
#> 
#>  400 clusters
#> coeffients:
#>     Estimate      S.E.   dU^-1/2 P-value
#> Z1  0.028262  0.130145  0.136188  0.8281
#> Z2 -0.149224  0.260989  0.272379  0.5675
#> 
#> exp(coeffients):
#>    Estimate    2.5%  97.5%
#> Z1  1.02866 0.79707 1.3276
#> Z2  0.86138 0.51646 1.4366
#> 
fgaugS2 <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fgaugS$E)
summary(fgaugS2)
#> 
#>    n events
#>  400     54
#> 
#>  400 clusters
#> coeffients:
#>     Estimate      S.E.   dU^-1/2 P-value
#> Z1  0.028262  0.130145  0.136188  0.8281
#> Z2 -0.149224  0.260989  0.272379  0.5675
#> 
#> exp(coeffients):
#>    Estimate    2.5%  97.5%
#> Z1  1.02866 0.79707 1.3276
#> Z2  0.86138 0.51646 1.4366
#> 
```
