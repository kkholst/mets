# Twostage survival model for multivariate survival data

Fits Clayton-Oakes or bivariate Plackett models for bivariate survival
data using marginals that are on Cox form. The dependence can be
modelled via

1.  Regression design on dependence parameter.

2.  Random effects, additive gamma model.

If clusters contain more than two subjects, we use a composite
likelihood based on the pairwise bivariate models, for full MLE see
twostageMLE.

The two-stage model is constructed such that given the gamma distributed
random effects it is assumed that the survival functions are indpendent,
and that the marginal survival functions are on Cox form (or additive
form) \$\$ P(T \> t\| x) = S(t\|x)= exp( -exp(x^T \beta) A_0(t) ) \$\$

One possibility is to model the variance within clusters via a
regression design, and then one can specify a regression structure for
the independent gamma distributed random effect for each cluster, such
that the variance is given by \$\$ \theta = h( z_j^T \alpha) \$\$ where
\\z\\ is specified by theta.des, and a possible link function var.link=1
will will use the exponential link \\h(x)=exp(x)\\, and var.link=0 the
identity link \\h(x)=x\\. The reported standard errors are based on the
estimated information from the likelihood assuming that the marginals
are known (unlike the twostageMLE and for the additive gamma model
below).

Can also fit a structured additive gamma random effects model, such as
the ACE, ADE model for survival data. In this case the random.design
specificies the random effects for each subject within a cluster. This
is a matrix of 1's and 0's with dimension n x d. With d random effects.
For a cluster with two subjects, we let the random.design rows be
\\v_1\\ and \\v_2\\. Such that the random effects for subject 1 is
\$\$v_1^T (Z_1,...,Z_d)\$\$, for d random effects. Each random effect
has an associated parameter \\(\lambda_1,...,\lambda_d)\\. By
construction subjects 1's random effect are Gamma distributed with mean
\\\lambda_j/v_1^T \lambda\\ and variance \\\lambda_j/(v_1^T
\lambda)^2\\. Note that the random effect \\v_1^T (Z_1,...,Z_d)\\ has
mean 1 and variance \\1/(v_1^T \lambda)\\. It is here asssumed that
\\lamtot=v_1^T \lambda\\ is fixed within clusters as it would be for the
ACE model below.

Based on these parameters the relative contribution (the heritability,
h) is equivalent to the expected values of the random effects:
\\\lambda_j/v_1^T \lambda\\

The DEFAULT parametrization (var.par=1) uses the variances of the random
effecs \$\$ \theta_j = \lambda_j/(v_1^T \lambda)^2 \$\$ For alternative
parametrizations one can specify how the parameters relate to
\\\lambda_j\\ with the argument var.par=0.

For both types of models the basic model assumptions are that given the
random effects of the clusters the survival distributions within a
cluster are independent and ' on the form \$\$ P(T \> t\| x,z) = exp( -Z
\cdot Laplace^{-1}(lamtot,lamtot,S(t\|x)) ) \$\$ with the inverse
laplace of the gamma distribution with mean 1 and variance 1/lamtot.

The parameters \\(\lambda_1,...,\lambda_d)\\ are related to the
parameters of the model by a regression construction \\pard\\ (d x k),
that links the \\d\\ \\\lambda\\ parameters with the (k) underlying
\\\theta\\ parameters \$\$ \lambda = theta.des \theta \$\$ here using
theta.des to specify these low-dimension association. Default is a
diagonal matrix. This can be used to make structural assumptions about
the variances of the random-effects as is needed for the ACE model for
example.

The case.control option that can be used with the pair specification of
the pairwise parts of the estimating equations. Here it is assumed that
the second subject of each pair is the proband.

## Usage

``` r
survival.twostage(
  margsurv,
  data = parent.frame(),
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

  Marginal model

- data:

  data frame

- method:

  Scoring method "nr", for lava NR optimizer

- detail:

  Detail

- clusters:

  Cluster variable

- silent:

  Debug information

- weights:

  Weights

- theta:

  Starting values for variance components

- theta.des:

  design for dependence parameters, when pairs are given the indeces of
  the theta-design for this pair, is given in pairs as column 5

- var.link:

  Link function for variance: exp-link.

- baseline.iid:

  to adjust for baseline estimation, using phreg function on same data.

- model:

  model

- marginal.trunc:

  marginal left truncation probabilities

- marginal.survival:

  optional vector of marginal survival probabilities

- strata:

  strata for fitting, see example

- se.clusters:

  for clusters for se calculation with iid

- numDeriv:

  to get numDeriv version of second derivative, otherwise uses sum of
  squared scores for each pair

- random.design:

  random effect design for additive gamma model, when pairs are given
  the indeces of the pairs random.design rows are given as columns 3:4

- pairs:

  matrix with rows of indeces (two-columns) for the pairs considered in
  the pairwise composite score, useful for case-control sampling when
  marginal is known.

- dim.theta:

  dimension of the theta parameter for pairs situation.

- numDeriv.method:

  uses simple to speed up things and second derivative not so important.

- additive.gamma.sum:

  for two.stage=0, this is specification of the lamtot in the models via
  a matrix that is multiplied onto the parameters theta
  (dimensions=(number random effects x number of theta parameters), when
  null then sums all parameters.

- var.par:

  is 1 for the default parametrization with the variances of the random
  effects, var.par=0 specifies that the \\\lambda_j\\'s are used as
  parameters.

- no.opt:

  for not optimizng

- ...:

  Additional arguments to maximizer NR of lava. and ascertained sampling

## References

Twostage estimation of additive gamma frailty models for survival data.
Scheike (2019), work in progress

Shih and Louis (1995) Inference on the association parameter in copula
models for bivariate survival data, Biometrics, (1995).

Glidden (2000), A Two-Stage estimator of the dependence parameter for
the Clayton Oakes model, LIDA, (2000).

Measuring early or late dependence for bivariate twin data Scheike,
Holst, Hjelmborg (2015), LIDA

Estimating heritability for cause specific mortality based on twins
studies Scheike, Holst, Hjelmborg (2014), LIDA

Additive Gamma frailty models for competing risks data, Biometrics
(2015) Eriksson and Scheike (2015),

## Author

Thomas Scheike

## Examples

``` r
library(mets)
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

### Plackett model
mph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
fitp <- survival.twostage(mph,data=diabetes,theta=3.0,Nit=40,
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
fitco2 <- survival.twostage(mph,data=diabetes,theta=0.0,detail=0,
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
fitco3 <- survival.twostage(margph,data=diabetes,theta=1.0,detail=0,
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
fitpa <- survival.twostage(marg,data=diabetes,theta=1.0,
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

fitcoa <- survival.twostage(marg,data=diabetes,theta=1.0,clusters=diabetes$id,
                 model="clayton.oakes")
summary(fitcoa)
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

d <- subset(simClaytonOakes(2000,2,0.5,0,stoptime=2,left=0),!truncated)
udp <- piecewise.twostage(c(0,0.5,2),data=d,method="optimize",
                          id="cluster",timevar="time",
                          status="status",model="clayton.oakes",silent=0)
#> Data-set  1 out of  4 
#>   Number of joint events: 530 of  2000
#> Data-set  2 out of  4 
#>   Number of joint events: 265 of  1218
#> Data-set  3 out of  4 
#>   Number of joint events: 246 of  1197
#> Data-set  4 out of  4 
#>   Number of joint events: 613 of  945
summary(udp)
#> [1] 1
#> Dependence parameter for Clayton-Oakes model 
#> log-coefficient for dependence parameter (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  0.655 (0.067)  0.621 (0.101)
#> 0.5 - 2  0.629 (0.097)  0.747 (0.058)
#> 
#> Kendall's tau (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  0.49  (0.017)  0.482 (0.025)
#> 0.5 - 2  0.484 (0.024)  0.513 (0.015)
#> 

 ## Reduce Ex.Timings
### Same model using the strata option, a bit slower
########################################################
## makes the survival pieces for different areas in the plane
##ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
##ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
##ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
##ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")

## everything done in one call
ud <- piecewise.data(c(0,0.5,2),data=d,timevar="time",status="status",id="cluster")
ud$strata <- factor(ud$strata);
ud$intstrata <- factor(ud$intstrata)

## makes strata specific id variable to identify pairs within strata
## se's computed based on the id variable across strata "cluster"
ud$idstrata <- ud$id+(as.numeric(ud$strata)-1)*2000

marg2 <- timereg::aalen(Surv(boxtime,status)~-1+factor(num):factor(intstrata),
               data=ud,n.sim=0,robust=0)
tdes <- model.matrix(~-1+factor(strata),data=ud)
fitp2 <- survival.twostage(marg2,data=ud,se.clusters=ud$cluster,clusters=ud$idstrata,
                model="clayton.oakes",theta.des=tdes,step=0.5)
summary(fitp2)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>                           log-Coef.         SE        z        P-val
#> factor(strata)0-0.5,0-0.5 0.8276088 0.06812882 12.14770 0.0000000000
#> factor(strata)0-0.5,0.5-2 0.4049929 0.11343320  3.57032 0.0003565449
#> factor(strata)0.5-2,0-0.5 0.4176728 0.11408832  3.66096 0.0002512718
#> factor(strata)0.5-2,0.5-2 0.6903006 0.05740192 12.02574 0.0000000000
#>                           Kendall tau         SE
#> factor(strata)0-0.5,0-0.5   0.5335648 0.01695545
#> factor(strata)0-0.5,0.5-2   0.4284558 0.02777768
#> factor(strata)0.5-2,0-0.5   0.4315636 0.02798774
#> factor(strata)0.5-2,0.5-2   0.4992884 0.01435045
#> 
#> $vargam
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> dependence1    2.288  0.1559 1.982 2.593 8.908e-49
#> dependence2    1.499  0.1701 1.166 1.833 1.189e-18
#> dependence3    1.518  0.1732 1.179 1.858 1.865e-18
#> dependence4    1.994  0.1145 1.770 2.219 5.715e-68
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
fitp3 <- survival.twostage(marg2,data=ud,clusters=ud$idstrata,se.cluster=ud$cluster,
                model="clayton.oakes",theta.des=tdes2,step=0.5)
summary(fitp3)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>                            log-Coef.         SE         z        P-val
#> factor(stratas)0-0.5,0-0.5 0.8276088 0.06812882 12.147705 0.000000e+00
#> factor(stratas)0-0.5,0.5-2 0.4110824 0.08945478  4.595422 4.318746e-06
#> factor(stratas)0.5-2,0.5-2 0.6903006 0.05740192 12.025741 0.000000e+00
#>                            Kendall tau         SE
#> factor(stratas)0-0.5,0-0.5   0.5335648 0.01695545
#> factor(stratas)0-0.5,0.5-2   0.4299477 0.02192471
#> factor(stratas)0.5-2,0.5-2   0.4992884 0.01435045
#> 
#> $vargam
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> dependence1    2.288  0.1559 1.982 2.593 8.908e-49
#> dependence2    1.508  0.1349 1.244 1.773 5.177e-29
#> dependence3    1.994  0.1145 1.770 2.219 5.715e-68
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

### same model using strata option, a bit slower
fitp4 <- survival.twostage(marg2,data=ud,clusters=ud$cluster,se.cluster=ud$cluster,
                model="clayton.oakes",theta.des=tdes2,step=0.5,strata=ud$strata)
summary(fitp4)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link 
#> $estimates
#>                            log-Coef.         SE         z        P-val
#> factor(stratas)0-0.5,0-0.5 0.8276088 0.06812882 12.147705 0.000000e+00
#> factor(stratas)0-0.5,0.5-2 0.4110824 0.08945478  4.595422 4.318746e-06
#> factor(stratas)0.5-2,0.5-2 0.6903006 0.05740192 12.025741 0.000000e+00
#>                            Kendall tau         SE
#> factor(stratas)0-0.5,0-0.5   0.5335648 0.01695545
#> factor(stratas)0-0.5,0.5-2   0.4299477 0.02192471
#> factor(stratas)0.5-2,0.5-2   0.4992884 0.01435045
#> 
#> $vargam
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> dependence1    2.288  0.1559 1.982 2.593 8.908e-49
#> dependence2    1.508  0.1349 1.244 1.773 5.177e-29
#> dependence3    1.994  0.1145 1.770 2.219 5.715e-68
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"


 ## Reduce Ex.Timings
### structured random effects model additive gamma ACE
### simulate structured two-stage additive gamma ACE model
data <- simClaytonOakes.twin.ace(4000,2,1,0,3)
out <- twin.polygen.design(data,id="cluster")
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
ts <- survival.twostage(aa,data=data,clusters=data$cluster,detail=0,
         theta=c(2,1),var.link=0,step=0.5,
         random.design=des.rv,theta.des=pardes)
summary(ts)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> $estimates
#>                Coef.        SE         z P-val Kendall tau         SE
#> dependence1 1.886703 0.1355750 13.916305     0    0.485425 0.01794927
#> dependence2 1.018832 0.1090384  9.343785     0    0.337492 0.02392940
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err   2.5%  97.5%   P-value
#> [1,]   0.6493 0.03797 0.5749 0.7238 1.416e-65
#> [2,]   0.3507 0.03797 0.2762 0.4251 2.569e-20
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%    P-value
#> p1    2.906 0.09363 2.722 3.089 1.902e-211
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

```
