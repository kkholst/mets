# Additive Random effects model for competing risks data for polygenetic modelling

Fits a random effects model describing the dependence in the cumulative
incidence curves for subjects within a cluster. Given the gamma
distributed random effects it is assumed that the cumulative incidence
curves are indpendent, and that the marginal cumulative incidence curves
are on additive form \$\$ P(T \leq t, cause=1 \| x,z) = P_1(t,x,z) = 1-
exp( -x^T A(t) - t z^T \beta) \$\$

## Usage

``` r
Grandom.cif(
  cif,
  data,
  cause = NULL,
  cif2 = NULL,
  times = NULL,
  cause1 = 1,
  cause2 = 1,
  cens.code = NULL,
  cens.model = "KM",
  Nit = 40,
  detail = 0,
  clusters = NULL,
  theta = NULL,
  theta.des = NULL,
  weights = NULL,
  step = 1,
  sym = 0,
  same.cens = FALSE,
  censoring.weights = NULL,
  silent = 1,
  var.link = 0,
  score.method = "nr",
  entry = NULL,
  estimator = 1,
  trunkp = 1,
  admin.cens = NULL,
  random.design = NULL,
  ...
)
```

## Arguments

- cif:

  a model object from the timereg::comp.risk function with the marginal
  cumulative incidence of cause2, i.e., the event that is conditioned
  on, and whose odds the comparision is made with respect to

- data:

  a data.frame with the variables.

- cause:

  specifies the causes related to the death times, the value cens.code
  is the censoring value.

- cif2:

  specificies model for cause2 if different from cause1.

- times:

  time points

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

- weights:

  weights for score equations.

- step:

  specifies the step size for the Newton-Raphson algorith.m

- sym:

  1 for symmetri and 0 otherwise

- same.cens:

  if true then censoring within clusters are assumed to be the same
  variable, default is independent censoring.

- censoring.weights:

  Censoring probabilities

- silent:

  debug information

- var.link:

  if var.link=1 then var is on log-scale.

- score.method:

  default uses "nlminb" optimzer, alternatively, use the "nr" algorithm.

- entry:

  entry-age in case of delayed entry. Then two causes must be given.

- estimator:

  estimator

- trunkp:

  gives probability of survival for delayed entry, and related to
  entry-ages given above.

- admin.cens:

  Administrative censoring

- random.design:

  specifies a regression design of 0/1's for the random effects.

- ...:

  extra arguments.

## Value

returns an object of type 'random.cif'. With the following arguments:

- theta:

  estimate of parameters of model.

- var.theta:

  variance for gamma.

- hess:

  the derivative of the used score.

- score:

  scores at final stage.

- theta.iid:

  matrix of iid decomposition of parametric effects.

## Details

We allow a regression structure for the indenpendent gamma distributed
random effects and their variances that may depend on cluster
covariates.

random.design specificies the random effects for each subject within a
cluster. This is a matrix of 1's and 0's with dimension n x d. With d
random effects. For a cluster with two subjects, we let the
random.design rows be \\v_1\\ and \\v_2\\. Such that the random effects
for subject 1 is \$\$v_1^T (Z_1,...,Z_d)\$\$, for d random effects. Each
random effect has an associated parameter \\(\lambda_1,...,\lambda_d)\\.
By construction subjects 1's random effect are Gamma distributed with
mean \\\lambda_1/v_1^T \lambda\\ and variance \\\lambda_1/(v_1^T
\lambda)^2\\. Note that the random effect \\v_1^T (Z_1,...,Z_d)\\ has
mean 1 and variance \\1/(v_1^T \lambda)\\.

The parameters \\(\lambda_1,...,\lambda_d)\\ are related to the
parameters of the model by a regression construction \\pard\\ (d x k),
that links the \\d\\ \\\lambda\\ parameters with the (k) underlying
\\\theta\\ parameters \$\$ \lambda = pard \theta \$\$

## References

A Semiparametric Random Effects Model for Multivariate Competing Risks
Data, Scheike, Zhang, Sun, Jensen (2010), Biometrika.

Cross odds ratio Modelling of dependence for Multivariate Competing
Risks Data, Scheike and Sun (2013), Biostatitistics.

Scheike, Holst, Hjelmborg (2014), LIDA, Estimating heritability for
cause specific hazards based on twin data

## Author

Thomas Scheike

## Examples

``` r
 ## Reduce Ex.Timings
 d <- simnordic.random(5000,delayed=TRUE,
       cordz=1.0,cormz=2,lam0=0.3,country=TRUE)
 times <- seq(50,90,by=10)
 addm <- timereg::comp.risk(Event(time,cause)~-1+factor(country)+cluster(id),data=d,
 times=times,cause=1,max.clust=NULL)

 ### making group indidcator 
 mm <- model.matrix(~-1+factor(zyg),d)

 out1m<-random.cif(addm,data=d,cause1=1,cause2=1,theta=1,
       theta.des=mm,same.cens=TRUE)
 summary(out1m)
#> Random effect variance for variation due to clusters
#> 
#> Cause 1 and cause 1
#> 
#> 
#>                   Coef.        SE        z        P-val Cross odds ratio
#> factor(zyg)MZ 2.2995032 0.3309552 6.948079 3.703038e-12         3.299503
#> factor(zyg)DZ 0.8440647 0.1826187 4.622005 3.800487e-06         1.844065
#>                      SE
#> factor(zyg)MZ 0.3309552
#> factor(zyg)DZ 0.1826187
 
 ## this model can also be formulated as a random effects model 
 ## but with different parameters
 out2m<-Grandom.cif(addm,data=d,cause1=1,cause2=1,
        theta=c(0.5,1),step=1.0,
        random.design=mm,same.cens=TRUE)
 summary(out2m)
#> Random effect parameters for additive gamma random effects 
#> 
#> Cause 1 and cause 1
#> 
#> 
#>      Coef.     SE    z   P-val
#> [1,] 0.435 0.0626 6.95 3.7e-12
#> [2,] 1.180 0.2560 4.62 3.8e-06
#> $estimate
#>      Coef.     SE    z   P-val
#> [1,] 0.435 0.0626 6.95 3.7e-12
#> [2,] 1.180 0.2560 4.62 3.8e-06
#> 
#> $h
#>    Estimate Std.Err   2.5% 97.5%   P-value
#> p1   0.3671 0.09538 0.1801 0.554 0.0001189
#> p2   1.0000 0.00000 1.0000 1.000 0.0000000
#> 
 1/out2m$theta
#>           [,1]
#> [1,] 2.2995032
#> [2,] 0.8440647
 out1m$theta
#>           [,1]
#> [1,] 2.2995032
#> [2,] 0.8440647
 
 ####################################################################
 ################### ACE modelling of twin data #####################
 ####################################################################
 ### assume that zygbin gives the zygosity of mono and dizygotic twins
 ### 0 for mono and 1 for dizygotic twins. We now formulate and AC model
 zygbin <- d$zyg=="DZ"

 n <- nrow(d)
 ### random effects for each cluster
 des.rv <- cbind(mm,(zygbin==1)*rep(c(1,0)),(zygbin==1)*rep(c(0,1)),1)
 ### design making parameters half the variance for dizygotic components
 pardes <- rbind(c(1,0), c(0.5,0),c(0.5,0), c(0.5,0), c(0,1))

 outacem <-Grandom.cif(addm,data=d,cause1=1,cause2=1,
    same.cens=TRUE,theta=c(0.35,0.15),
            step=1.0,theta.des=pardes,random.design=des.rv)
 summary(outacem)
#> Random effect parameters for additive gamma random effects 
#> 
#> Cause 1 and cause 1
#> 
#> 
#>        Coef.     SE      z    P-val
#> [1,]  0.4660 0.0787  5.920 3.14e-09
#> [2,] -0.0313 0.0815 -0.384 7.01e-01
#> $estimate
#>        Coef.     SE      z    P-val
#> [1,]  0.4660 0.0787  5.920 3.14e-09
#> [2,] -0.0313 0.0815 -0.384 7.01e-01
#> 
#> $h
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> p1  1.07202  0.1921  0.6956 1.4485 2.383e-08
#> p2 -0.07202  0.1921 -0.4485 0.3044 7.077e-01
#> 

```
