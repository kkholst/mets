# Twostage Survival Model for Multivariate Survival Data

Fits Clayton-Oakes or bivariate Plackett models for bivariate survival
data using marginals that are on Cox form. The dependence can be
modelled via:

1.  Regression design on dependence parameter.

2.  Random effects, additive gamma model.

## Usage

``` r
survival_twostage(
  margsurv,
  data = NULL,
  method = "nr",
  detail = 0,
  clusters = NULL,
  silent = 1,
  weights = NULL,
  theta = NULL,
  theta.des = NULL,
  var.link = 1,
  baseline.iid = 1,
  model = "clayton.oakes",
  marginal.trunc = NULL,
  marginal.survival = NULL,
  strata = NULL,
  se.clusters = NULL,
  numDeriv = 1,
  random.design = NULL,
  pairs = NULL,
  dim.theta = NULL,
  numDeriv.method = "simple",
  additive.gamma.sum = NULL,
  var.par = 1,
  no.opt = FALSE,
  ...
)
```

## Arguments

- margsurv:

  Marginal model object.

- data:

  Data frame (must be given).

- method:

  Scoring method: `"nr"` for lava NR optimizer.

- detail:

  Detail level for output.

- clusters:

  Cluster variable.

- silent:

  Debug information level.

- weights:

  Weights for score equations.

- theta:

  Starting values for variance components.

- theta.des:

  Design for dependence parameters; when pairs are given, the indices of
  the theta-design for this pair are given in pairs as column 5.

- var.link:

  Link function for variance: 1 for exponential link.

- baseline.iid:

  To adjust for baseline estimation, using `phreg` function on same
  data.

- model:

  Model type: `"clayton.oakes"` or `"plackett"`.

- marginal.trunc:

  Marginal left truncation probabilities.

- marginal.survival:

  Optional vector of marginal survival probabilities.

- strata:

  Strata for fitting (see examples).

- se.clusters:

  Clusters for SE calculation with IID.

- numDeriv:

  To get numDeriv version of second derivative; otherwise uses sum of
  squared scores for each pair.

- random.design:

  Random effect design for additive gamma model; when pairs are given,
  the indices of the pairs' random.design rows are given as columns 3:4.

- pairs:

  Matrix with rows of indices (two-columns) for the pairs considered in
  the pairwise composite score; useful for case-control sampling when
  marginal is known.

- dim.theta:

  Dimension of the theta parameter for pairs situation.

- numDeriv.method:

  Method for numerical derivative (e.g., `"simple"` to speed up things).

- additive.gamma.sum:

  For `two.stage=0`, this is specification of the `lamtot` in the models
  via a matrix that is multiplied onto the parameters theta (dimensions
  = number of random effects \\\times\\ number of theta parameters);
  when `NULL`, sums all parameters.

- var.par:

  Is 1 for the default parametrization with the variances of the random
  effects; `var.par=0` specifies that the \\\lambda_j\\'s are used as
  parameters.

- no.opt:

  For not optimizing.

- ...:

  Additional arguments to maximizer NR of lava.

## Value

An object of class `"mets.twostage"` containing:

- theta:

  Estimated dependence parameters.

- coef:

  Coefficients.

- score:

  Score vector.

- hess:

  Hessian matrix.

- hessi:

  Inverse Hessian matrix.

- var.theta:

  Variance of theta parameters.

- robvar.theta:

  Robust variance of theta parameters.

- loglike:

  Log-likelihood value.

- theta.iid:

  Influence functions for theta.

- model:

  Model type used.

## Details

If clusters contain more than two subjects, a composite likelihood based
on pairwise bivariate models is used (for full MLE see `twostageMLE`).

The two-stage model is constructed such that, given gamma distributed
random effects, the survival functions are assumed independent, and the
marginal survival functions are on Cox form: \$\$ P(T \> t\| x) =
S(t\|x)= \exp( -\exp(x^T \beta) A_0(t) ) \$\$

One possibility is to model the variance within clusters via a
regression design, specifying a regression structure for the independent
gamma distributed random effect for each cluster, such that the variance
is given by: \$\$ \theta = h( z_j^T \alpha) \$\$ where \\z\\ is
specified by `theta.des`, and a possible link function `var.link=1` will
use the exponential link \\h(x)=\exp(x)\\, and `var.link=0` the identity
link \\h(x)=x\\.

The reported standard errors are based on the estimated information from
the likelihood assuming that the marginals are known (unlike
`twostageMLE` and for the additive gamma model below).

Can also fit a structured additive gamma random effects model, such as
the ACE, ADE model for survival data. In this case, the `random.design`
specifies the random effects for each subject within a cluster. This is
a matrix of 1's and 0's with dimension \\n \times d\\ (for \\d\\ random
effects).

For a cluster with two subjects, the `random.design` rows are \\v_1\\
and \\v_2\\. Such that the random effect for subject 1 is \$\$v_1^T
(Z_1,...,Z_d)\$\$, for \\d\\ random effects. Each random effect has an
associated parameter \\(\lambda_1,...,\lambda_d)\\. By construction,
subject 1's random effect is Gamma distributed with mean
\\\lambda_j/v_1^T \lambda\\ and variance \\\lambda_j/(v_1^T
\lambda)^2\\. Note that the random effect \\v_1^T (Z_1,...,Z_d)\\ has
mean 1 and variance \\1/(v_1^T \lambda)\\. It is assumed that
\\lamtot=v_1^T \lambda\\ is fixed within clusters as it would be for the
ACE model.

Based on these parameters, the relative contribution (the heritability,
\\h\\) is equivalent to the expected values of the random effects:
\\\lambda_j/v_1^T \lambda\\.

The DEFAULT parametrization (`var.par=1`) uses the variances of the
random effects: \$\$ \theta_j = \lambda_j/(v_1^T \lambda)^2 \$\$ For
alternative parametrizations, specify how the parameters relate to
\\\lambda_j\\ with the argument `var.par=0`.

For both types of models, the basic model assumptions are that given the
random effects of the clusters, the survival distributions within a
cluster are independent and on the form: \$\$ P(T \> t\| x,z) = \exp( -Z
\cdot \text{Laplace}^{-1}(lamtot,lamtot,S(t\|x)) ) \$\$ with the inverse
Laplace of the gamma distribution with mean 1 and variance \\1/lamtot\\.

The parameters \\(\lambda_1,...,\lambda_d)\\ are related to the
parameters of the model by a regression construction `pard` (\\d \times
k\\), that links the \\d\\ \\\lambda\\ parameters with the \\k\\
underlying \\\theta\\ parameters: \$\$ \lambda = theta.des \times \theta
\$\$ here using `theta.des` to specify these low-dimension associations.
Default is a diagonal matrix. This can be used to make structural
assumptions about the variances of the random-effects as is needed for
the ACE model for example.

The `case.control` option can be used with the pair specification of the
pairwise parts of the estimating equations. Here it is assumed that the
second subject of each pair is the proband.

## References

- Twostage estimation of additive gamma frailty models for survival
  data. Scheike (2019), work in progress.

- Shih and Louis (1995) Inference on the association parameter in copula
  models for bivariate survival data, *Biometrics*.

- Glidden (2000), A Two-Stage estimator of the dependence parameter for
  the Clayton Oakes model, LIDA.

- Measuring early or late dependence for bivariate twin data. Scheike,
  Holst, Hjelmborg (2015), LIDA.

- Estimating heritability for cause specific mortality based on twins
  studies. Scheike, Holst, Hjelmborg (2014), LIDA.

- Additive Gamma frailty models for competing risks data. Eriksson and
  Scheike (2015), *Biometrics*.

## Author

Thomas Scheike

## Examples

``` r
data(diabetes)

# Marginal Cox model  with treat as covariate
margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
### Clayton-Oakes, MLE
fitco1<-twostageMLE(margph,data=diabetes,theta=1.0)
summary(fitco1)
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

 ## Reduce Ex.Timings
### Plackett model
mph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
fitp <- survival_twostage(mph,data=diabetes,theta=3.0,Nit=40,
               clusters=diabetes$id,var.link=1,model="plackett")
summary(fitp)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>             log-Coef.        SE        z        P-val Spearman Corr.         SE
#> dependence1   1.14188 0.3026537 3.772891 0.0001613666      0.3648216 0.08869229
#> 
#> $or
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> dependence1    3.133  0.9481 1.274 4.991 0.0009528
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

### Clayton-Oakes
fitco2 <- survival_twostage(mph,data=diabetes,theta=0.0,detail=0,
                 clusters=diabetes$id,var.link=1,model="clayton.oakes")
summary(fitco2)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>              log-Coef.        SE          z     P-val Kendall tau         SE
#> dependence1 -0.0484957 0.3718487 -0.1304178 0.8962359    0.322645 0.08126576
#> 
#> $vargam
#>             Estimate Std.Err   2.5% 97.5%  P-value
#> dependence1   0.9527  0.3542 0.2584 1.647 0.007161
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
fitco3 <- survival_twostage(margph,data=diabetes,theta=1.0,detail=0,
             clusters=diabetes$id,var.link=0,model="clayton.oakes")
summary(fitco3)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> $estimates
#>                 Coef.     SE        z       P-val Kendall tau         SE
#> dependence1 0.9526614 0.3543 2.688855 0.007169754    0.322645 0.08127816
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

### without covariates but with stratafied
marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
fitpa <- survival_twostage(marg,data=diabetes,theta=1.0,
                clusters=diabetes$id,model="clayton.oakes")
summary(fitpa)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>               log-Coef.        SE          z     P-val Kendall tau         SE
#> dependence1 -0.05684062 0.3721207 -0.1527478 0.8785971    0.320824 0.08108359
#> 
#> $vargam
#>             Estimate Std.Err   2.5% 97.5%  P-value
#> dependence1   0.9447  0.3516 0.2557 1.634 0.007203
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"


### Piecewise constant cross hazards ratio modelling
########################################################

d <- subset(sim_ClaytonOakes(1000,2,0.5,0,stoptime=2,left=0),!truncated)
udp <- piecewise_twostage(c(0,0.5,2),data=d,method="optimize",
                          id="cluster",timevar="time",
                          status="status",model="clayton.oakes",silent=0)
#> Data-set  1 out of  4 
#>   Number of joint events: 262 of  1000
#> Data-set  2 out of  4 
#>   Number of joint events: 136 of  608
#> Data-set  3 out of  4 
#>   Number of joint events: 128 of  601
#> Data-set  4 out of  4 
#>   Number of joint events: 306 of  471
summary(udp)
#> [1] 1
#> Dependence parameter for Clayton-Oakes model 
#> log-coefficient for dependence parameter (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  0.593 (0.100)  0.819 (0.125)
#> 0.5 - 2  0.578 (0.133)  0.741 (0.079)
#> 
#> Kendall's tau (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  0.475 (0.025)  0.531 (0.031)
#> 0.5 - 2  0.471 (0.033)  0.512 (0.020)
#> 


 ## Reduce Ex.Timings
### Same model using the strata option, a bit slower
########################################################
## makes the survival pieces for different areas in the plane
##ud1=surv_boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
##ud2=surv_boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
##ud3=surv_boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
##ud4=surv_boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")

## everything done in one call
ud <- piecewise_data(c(0,0.5,2),data=d,timevar="time",status="status",id="cluster")
ud$strata <- factor(ud$strata);
ud$intstrata <- factor(ud$intstrata)

## makes strata specific id variable to identify pairs within strata
## se's computed based on the id variable across strata "cluster"
ud$idstrata <- ud$id+(as.numeric(ud$strata)-1)*2000

marg2 <- timereg::aalen(Surv(boxtime,status)~-1+factor(num):factor(intstrata),
               data=ud,n.sim=0,robust=0)
tdes <- model.matrix(~-1+factor(strata),data=ud)
fitp2 <- survival_twostage(marg2,data=ud,se.clusters=ud$cluster,clusters=ud$idstrata,
                model="clayton.oakes",theta.des=tdes,step=0.5)
summary(fitp2)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>                           log-Coef.         SE        z        P-val
#> factor(strata)0-0.5,0-0.5 0.7677352 0.10124764 7.582747 3.375078e-14
#> factor(strata)0-0.5,0.5-2 0.5719969 0.14361250 3.982919 6.807407e-05
#> factor(strata)0.5-2,0-0.5 0.3717046 0.15778635 2.355746 1.848555e-02
#> factor(strata)0.5-2,0.5-2 0.6915363 0.07753039 8.919552 0.000000e+00
#>                           Kendall tau         SE
#> factor(strata)0-0.5,0-0.5   0.5186384 0.02527674
#> factor(strata)0-0.5,0.5-2   0.4697494 0.03577171
#> factor(strata)0.5-2,0-0.5   0.4203242 0.03844492
#> factor(strata)0.5-2,0.5-2   0.4995973 0.01938258
#> 
#> $vargam
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> dependence1    2.155  0.2182 1.727 2.582 5.250e-23
#> dependence2    1.772  0.2545 1.273 2.271 3.327e-12
#> dependence3    1.450  0.2288 1.002 1.899 2.332e-10
#> dependence4    1.997  0.1548 1.693 2.300 4.609e-38
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
ud$stratas <- ud$strata;
ud$stratas[ud$strata=="0.5-2,0-0.5"] <- "0-0.5,0.5-2"
tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
fitp3 <- survival_twostage(marg2,data=ud,clusters=ud$idstrata,se.cluster=ud$cluster,
                model="clayton.oakes",theta.des=tdes2,step=0.5)
summary(fitp3)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>                            log-Coef.         SE        z        P-val
#> factor(stratas)0-0.5,0-0.5 0.7677352 0.10124764 7.582747 3.375078e-14
#> factor(stratas)0-0.5,0.5-2 0.4740701 0.12104406 3.916509 8.984052e-05
#> factor(stratas)0.5-2,0.5-2 0.6915363 0.07753039 8.919552 0.000000e+00
#>                            Kendall tau         SE
#> factor(stratas)0-0.5,0-0.5   0.5186384 0.02527674
#> factor(stratas)0-0.5,0.5-2   0.4454487 0.02990081
#> factor(stratas)0.5-2,0.5-2   0.4995973 0.01938258
#> 
#> $vargam
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> dependence1    2.155  0.2182 1.727 2.582 5.250e-23
#> dependence2    1.607  0.1945 1.225 1.988 1.439e-16
#> dependence3    1.997  0.1548 1.693 2.300 4.609e-38
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

### same model using strata option, a bit slower
fitp4 <- survival_twostage(marg2,data=ud,clusters=ud$cluster,se.cluster=ud$cluster,
                model="clayton.oakes",theta.des=tdes2,step=0.5,strata=ud$strata)
summary(fitp4)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>                            log-Coef.         SE        z        P-val
#> factor(stratas)0-0.5,0-0.5 0.7677352 0.10124764 7.582747 3.375078e-14
#> factor(stratas)0-0.5,0.5-2 0.4740701 0.12104406 3.916509 8.984052e-05
#> factor(stratas)0.5-2,0.5-2 0.6915363 0.07753039 8.919552 0.000000e+00
#>                            Kendall tau         SE
#> factor(stratas)0-0.5,0-0.5   0.5186384 0.02527674
#> factor(stratas)0-0.5,0.5-2   0.4454487 0.02990081
#> factor(stratas)0.5-2,0.5-2   0.4995973 0.01938258
#> 
#> $vargam
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> dependence1    2.155  0.2182 1.727 2.582 5.250e-23
#> dependence2    1.607  0.1945 1.225 1.988 1.439e-16
#> dependence3    1.997  0.1548 1.693 2.300 4.609e-38
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"


 ## Reduce Ex.Timings
### structured random effects model additive gamma ACE
### simulate structured two-stage additive gamma ACE model
data <- sim_ClaytonOakes_twin_ace(2000,2,1,0,3)
out <- twin_polygen_design(data,id="cluster")
pardes <- out$pardes
pardes
#>      [,1] [,2]
#> [1,]  1.0    0
#> [2,]  0.5    0
#> [3,]  0.5    0
#> [4,]  0.5    0
#> [5,]  0.0    1
des.rv <- out$des.rv
head(des.rv)
#>   MZ DZ DZns1 DZns2 env
#> 1  1  0     0     0   1
#> 2  1  0     0     0   1
#> 3  1  0     0     0   1
#> 4  1  0     0     0   1
#> 5  1  0     0     0   1
#> 6  1  0     0     0   1
aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data,robust=0)
ts <- survival_twostage(aa,data=data,clusters=data$cluster,detail=0,
         theta=c(2,1),var.link=0,step=0.5,
         random.design=des.rv,theta.des=pardes)
summary(ts)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> $estimates
#>                 Coef.        SE         z        P-val Kendall tau         SE
#> dependence1 2.1349984 0.2124813 10.047936 0.000000e+00   0.5163239 0.02485421
#> dependence2 0.7883763 0.1698435  4.641782 3.454178e-06   0.2827367 0.04368940
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err   2.5%  97.5%   P-value
#> [1,]   0.7303 0.05871 0.6153 0.8454 1.584e-35
#> [2,]   0.2697 0.05871 0.1546 0.3847 4.355e-06
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1    2.923  0.1387 2.652 3.195 1.276e-98
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```
