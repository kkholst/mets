# Twostage survival model fitted by pseudo MLE

Fits Clayton-Oakes clustered survival data using marginals that are on
Cox form in the likelihood for the dependence parameter as in Glidden
(2000). The dependence can be modelled via a

1.  Regression design on dependence parameter.

We allow a regression structure for the indenpendent gamma distributed
random effects and their variances that may depend on cluster
covariates. So \$\$ \theta = h( z_j^T \alpha) \$\$ where \\z\\ is
specified by theta.des . The link function can be the exp when
var.link=1

## Usage

``` r
twostageMLE(
  margsurv,
  data = parent.frame(),
  theta = NULL,
  theta.des = NULL,
  var.link = 0,
  method = "NR",
  no.opt = FALSE,
  weights = NULL,
  se.cluster = NULL,
  ...
)
```

## Arguments

- margsurv:

  Marginal model from phreg

- data:

  data frame

- theta:

  Starting values for variance components

- theta.des:

  design for dependence parameters, when pairs are given this is could
  be a (pairs) x (numer of parameters) x (max number random effects)
  matrix

- var.link:

  Link function for variance if 1 then uses exp link

- method:

  type of opitmizer, default is Newton-Raphson "NR"

- no.opt:

  to not optimize, for example to get score and iid for specific theta

- weights:

  cluster specific weights, but given with length equivalent to
  data-set, weights for score equations

- se.cluster:

  specifies how the influence functions are summed before squared when
  computing the variance. Note that the id from the marginal model is
  used to construct MLE, and then these scores can be summed with the
  se.cluster argument.

- ...:

  arguments to be passed to optimizer

## References

Measuring early or late dependence for bivariate twin data Scheike,
Holst, Hjelmborg (2015), LIDA

Twostage modelling of additive gamma frailty models for survival data.
Scheike and Holst, working paper

Shih and Louis (1995) Inference on the association parameter in copula
models for bivariate survival data, Biometrics, (1995).

Glidden (2000), A Two-Stage estimator of the dependence parameter for
the Clayton Oakes model, LIDA, (2000).

## Author

Thomas Scheike

## Examples

``` r
data(diabetes)
dd <- phreg(Surv(time,status==1)~treat+cluster(id),diabetes)
oo <- twostageMLE(dd,data=diabetes)
summary(oo)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> $estimates
#>                 Coef.        SE       z       P-val Kendall tau         SE
#> dependence1 0.9526614 0.3543033 2.68883 0.007170289    0.322645 0.08127892
#> 
#> $type
#> NULL
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

theta.des <- model.matrix(~-1+factor(adult),diabetes)

oo <-twostageMLE(dd,data=diabetes,theta.des=theta.des)
summary(oo)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> $estimates
#>                    Coef.        SE        z      P-val Kendall tau         SE
#> factor(adult)1 0.9117633 0.4000030 2.279391 0.02264381   0.3131310 0.09435851
#> factor(adult)2 1.0570600 0.7014182 1.507032 0.13180233   0.3457767 0.15010636
#> 
#> $type
#> NULL
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```
