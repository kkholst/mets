# Random effects model for competing risks data

Fits a random effects model describing the dependence in the cumulative
incidence curves for subjects within a cluster. Given the gamma
distributed random effects it is assumed that the cumulative incidence
curves are indpendent, and that the marginal cumulative incidence curves
are on the form \$\$ P(T \leq t, cause=1 \| x,z) = P_1(t,x,z) = 1- exp(
-x^T A(t) exp(z^T \beta)) \$\$ We allow a regression structure for the
random effects variances that may depend on cluster covariates.

## Usage

``` r
random.cif(
  cif,
  data,
  cause = NULL,
  cif2 = NULL,
  cause1 = 1,
  cause2 = 1,
  cens.code = NULL,
  cens.model = "KM",
  Nit = 40,
  detail = 0,
  clusters = NULL,
  theta = NULL,
  theta.des = NULL,
  sym = 1,
  step = 1,
  same.cens = FALSE,
  var.link = 0,
  score.method = "nr",
  entry = NULL,
  trunkp = 1,
  ...
)
```

## Arguments

- cif:

  a model object from the comp.risk function with the marginal
  cumulative incidence of cause2, i.e., the event that is conditioned
  on, and whose odds the comparision is made with respect to

- data:

  a data.frame with the variables.

- cause:

  specifies the causes related to the death times, the value cens.code
  is the censoring value.

- cif2:

  specificies model for cause2 if different from cause1.

- cause1:

  cause of first coordinate.

- cause2:

  cause of second coordinate.

- cens.code:

  specificies the code for the censoring if NULL then uses the one from
  the marginal cif model.

- cens.model:

  specified which model to use for the ICPW, KM is Kaplan-Meier
  alternatively it may be "cox"

- Nit:

  number of iterations for Newton-Raphson algorithm.

- detail:

  if 0 no details are printed during iterations, if 1 details are given.

- clusters:

  specifies the cluster structure.

- theta:

  specifies starting values for the cross-odds-ratio parameters of the
  model.

- theta.des:

  specifies a regression design for the cross-odds-ratio parameters.

- sym:

  1 for symmetry 0 otherwise

- step:

  specifies the step size for the Newton-Raphson algorith.m

- same.cens:

  if true then censoring within clusters are assumed to be the same
  variable, default is independent censoring.

- var.link:

  if var.link=1 then var is on log-scale.

- score.method:

  default uses "nlminb" optimzer, alternatively, use the "nr" algorithm.

- entry:

  entry-age in case of delayed entry. Then two causes must be given.

- trunkp:

  gives probability of survival for delayed entry, and related to
  entry-ages given above.

- ...:

  extra arguments.

## Value

returns an object of type 'cor'. With the following arguments:

- theta:

  estimate of proportional odds parameters of model.

- var.theta:

  variance for gamma.

- hess:

  the derivative of the used score.

- score:

  scores at final stage.

- score:

  scores at final stage.

- theta.iid:

  matrix of iid decomposition of parametric effects.

## References

A Semiparametric Random Effects Model for Multivariate Competing Risks
Data, Scheike, Zhang, Sun, Jensen (2010), Biometrika.

Cross odds ratio Modelling of dependence for Multivariate Competing
Risks Data, Scheike and Sun (2012), work in progress.

## Author

Thomas Scheike

## Examples

``` r
 ## Reduce Ex.Timings
 d <- simnordic.random(5000,delayed=TRUE,cordz=0.5,cormz=2,lam0=0.3,country=TRUE)
 times <- seq(50,90,by=10)
 add1 <- timereg::comp.risk(Event(time,cause)~-1+factor(country)+cluster(id),data=d,
 times=times,cause=1,max.clust=NULL)

 ### making group indidcator 
 mm <- model.matrix(~-1+factor(zyg),d)

 out1<-random.cif(add1,data=d,cause1=1,cause2=1,theta=1,same.cens=TRUE)
 summary(out1)
#> Random effect variance for variation due to clusters
#> 
#> Cause 1 and cause 1
#> 
#> 
#>              Coef.        SE        z        P-val Cross odds ratio        SE
#> intercept 1.104589 0.1416994 7.795297 6.439294e-15         2.104589 0.1416994

 out2<-random.cif(add1,data=d,cause1=1,cause2=1,theta=1,
       theta.des=mm,same.cens=TRUE)
 summary(out2)
#> Random effect variance for variation due to clusters
#> 
#> Cause 1 and cause 1
#> 
#> 
#>                   Coef.        SE        z        P-val Cross odds ratio
#> factor(zyg)MZ 2.1534461 0.3028887 7.109694 1.163070e-12         3.153446
#> factor(zyg)DZ 0.3138631 0.1303997 2.406931 1.608721e-02         1.313863
#>                      SE
#> factor(zyg)MZ 0.3028887
#> factor(zyg)DZ 0.1303997

#########################################
##### 2 different causes
#########################################

 add2 <- timereg::comp.risk(Event(time,cause)~-1+factor(country)+cluster(id),data=d,
                  times=times,cause=2,max.clust=NULL)
 out3 <- random.cif(add1,data=d,cause1=1,cause2=2,cif2=add2,sym=1,same.cens=TRUE)
 summary(out3) ## negative dependence
#> Random effect variance for variation due to clusters
#> 
#> Cause 1 and cause 2
#> 
#> 
#>                Coef.         SE         z P-val Cross odds ratio         SE
#> intercept -0.3341633 0.02629939 -12.70612     0        0.6658367 0.02629939

 out4 <- random.cif(add1,data=d,cause1=1,cause2=2,cif2=add2,theta.des=mm,sym=1,same.cens=TRUE)
 summary(out4) ## negative dependence
#> Random effect variance for variation due to clusters
#> 
#> Cause 1 and cause 2
#> 
#> 
#>                    Coef.         SE          z        P-val Cross odds ratio
#> factor(zyg)MZ -0.4502417 0.03083045 -14.603799 0.000000e+00        0.5497583
#> factor(zyg)DZ -0.2007758 0.04502517  -4.459191 8.226949e-06        0.7992242
#>                       SE
#> factor(zyg)MZ 0.03083045
#> factor(zyg)DZ 0.04502517
```
