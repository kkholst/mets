# Fits Clayton-Oakes or bivariate Plackett (OR) models for binary data using marginals that are on logistic form. If clusters contain more than two times, the algoritm uses a compososite likelihood based on all pairwise bivariate models.

The pairwise pairwise odds ratio model provides an alternative to the
alternating logistic regression (ALR).

## Usage

``` r
binomial_twostage(
  margbin,
  data = parent.frame(),
  method = "nr",
  detail = 0,
  clusters = NULL,
  silent = 1,
  weights = NULL,
  theta = NULL,
  theta.des = NULL,
  var.link = 0,
  var.par = 1,
  var.func = NULL,
  iid = 1,
  notaylor = 1,
  model = "plackett",
  marginal.p = NULL,
  beta.iid = NULL,
  Dbeta.iid = NULL,
  strata = NULL,
  max.clust = NULL,
  se.clusters = NULL,
  numDeriv = 0,
  random.design = NULL,
  pairs = NULL,
  dim.theta = NULL,
  additive.gamma.sum = NULL,
  pair.ascertained = 0,
  case.control = 0,
  no.opt = FALSE,
  twostage = 1,
  beta = NULL,
  ...
)
```

## Arguments

- margbin:

  Marginal binomial model

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

  Weights for log-likelihood, can be used for each type of outcome in
  2x2 tables.

- theta:

  Starting values for variance components

- theta.des:

  design for dependence parameters, when pairs are given the indeces of
  the theta-design for this pair, is given in pairs as column 5

- var.link:

  Link function for variance

- var.par:

  parametrization

- var.func:

  when alternative parametrizations are used this function can specify
  how the paramters are related to the \\\lambda_j\\'s.

- iid:

  Calculate i.i.d. decomposition when iid\>=1, when iid=2 then avoids
  adding the uncertainty for marginal paramters for additive gamma model
  (default).

- notaylor:

  Taylor expansion

- model:

  model

- marginal.p:

  vector of marginal probabilities

- beta.iid:

  iid decomposition of marginal probability estimates for each subject,
  if based on GLM model this is computed.

- Dbeta.iid:

  derivatives of marginal model wrt marginal parameters, if based on GLM
  model this is computed.

- strata:

  strata for fitting: considers only pairs where both are from same
  strata

- max.clust:

  max clusters

- se.clusters:

  clusters for iid decomposition for roubst standard errors

- numDeriv:

  uses Fisher scoring aprox of second derivative if 0, otherwise
  numerical derivatives

- random.design:

  random effect design for additive gamma model, when pairs are given
  the indeces of the pairs random.design rows are given as columns 3:4

- pairs:

  matrix with rows of indeces (two-columns) for the pairs considered in
  the pairwise composite score, useful for case-control sampling when
  marginal is known.

- dim.theta:

  dimension of theta when pairs and pairs specific design is given. That
  is when pairs has 6 columns.

- additive.gamma.sum:

  this is specification of the lamtot in the models via a matrix that is
  multiplied onto the parameters theta (dimensions=(number random
  effects x number of theta parameters), when null then sums all
  parameters. Default is a matrix of 1's

- pair.ascertained:

  if pairs are sampled only when there are events in the pair i.e.
  Y1+Y2\>=1.

- case.control:

  if data is case control data for pair call, and here 2nd column of
  pairs are probands (cases or controls)

- no.opt:

  for not optimizing

- twostage:

  default twostage=1, to fit MLE use twostage=0

- beta:

  is starting value for beta for MLE version

- ...:

  for NR of lava

## Details

The reported standard errors are based on a cluster corrected score
equations from the pairwise likelihoods assuming that the marginals are
known. This gives correct standard errors in the case of the Odds-Ratio
model (Plackett distribution) for dependence, but incorrect standard
errors for the Clayton-Oakes types model (that is also called
"gamma"-frailty). For the additive gamma version of the standard errors
are adjusted for the uncertainty in the marginal models via an iid
deomposition using the iid() function of lava. For the clayton oakes
model that is not speicifed via the random effects these can be fixed
subsequently using the iid influence functions for the marginal model,
but typically this does not change much.

For the Clayton-Oakes version of the model, given the gamma distributed
random effects it is assumed that the probabilities are indpendent, and
that the marginal survival functions are on logistic form \$\$
logit(P(Y=1\|X)) = \alpha + x^T \beta \$\$ therefore conditional on the
random effect the probability of the event is \$\$ logit(P(Y=1\|X,Z)) =
exp( -Z \cdot Laplace^{-1}(lamtot,lamtot,P(Y=1\|x)) ) \$\$

Can also fit a structured additive gamma random effects model, such the
ACE, ADE model for survival data:

Now random.design specificies the random effects for each subject within
a cluster. This is a matrix of 1's and 0's with dimension n x d. With d
random effects. For a cluster with two subjects, we let the
random.design rows be \\v_1\\ and \\v_2\\. Such that the random effects
for subject 1 is \$\$v_1^T (Z_1,...,Z_d)\$\$, for d random effects. Each
random effect has an associated parameter \\(\lambda_1,...,\lambda_d)\\.
By construction subjects 1's random effect are Gamma distributed with
mean \\\lambda_j/v_1^T \lambda\\ and variance \\\lambda_j/(v_1^T
\lambda)^2\\. Note that the random effect \\v_1^T (Z_1,...,Z_d)\\ has
mean 1 and variance \\1/(v_1^T \lambda)\\. It is here asssumed that
\\lamtot=v_1^T \lambda\\ is fixed over all clusters as it would be for
the ACE model below.

The DEFAULT parametrization uses the variances of the random effecs
(var.par=1) \$\$ \theta_j = \lambda_j/(v_1^T \lambda)^2 \$\$

For alternative parametrizations (var.par=0) one can specify how the
parameters relate to \\\lambda_j\\ with the function

Based on these parameters the relative contribution (the heritability,
h) is equivalent to the expected values of the random effects
\\\lambda_j/v_1^T \lambda\\

Given the random effects the probabilities are independent and on the
form \$\$ logit(P(Y=1\|X)) = exp( -
Laplace^{-1}(lamtot,lamtot,P(Y=1\|x)) ) \$\$ with the inverse laplace of
the gamma distribution with mean 1 and variance lamtot.

The parameters \\(\lambda_1,...,\lambda_d)\\ are related to the
parameters of the model by a regression construction \\pard\\ (d x k),
that links the \\d\\ \\\lambda\\ parameters with the (k) underlying
\\\theta\\ parameters \$\$ \lambda = theta.des \theta \$\$ here using
theta.des to specify these low-dimension association. Default is a
diagonal matrix.

## References

Two-stage binomial modelling

## Author

Thomas Scheike

## Examples

``` r
data(twinstut)
twinstut0 <- subset(twinstut, tvparnr<4000)
twinstut <- twinstut0
twinstut$binstut <- (twinstut$stutter=="yes")*1
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
margbin <- glm(binstut~factor(sex)+age,data=twinstut,family=binomial())
bin <- binomial_twostage(margbin,data=twinstut,var.link=1,
         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0)
summary(bin)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                    theta        se
#> factor(zyg)dz -0.2853738 0.9894082
#> factor(zyg)mz  3.3391390 0.5590195
#> factor(zyg)os  0.4920396 0.7634939
#> 
#> $or
#>               Estimate Std.Err   2.5%  97.5% P-value
#> factor(zyg)dz   0.7517  0.7438 -0.706  2.209 0.31216
#> factor(zyg)mz  28.1948 15.7615 -2.697 59.087 0.07364
#> factor(zyg)os   1.6356  1.2488 -0.812  4.083 0.19027
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

twinstut$cage <- scale(twinstut$age)
theta.des <- model.matrix( ~-1+factor(zyg)+cage,data=twinstut)
bina <- binomial_twostage(margbin,data=twinstut,var.link=1,
             clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bina)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                    theta        se
#> factor(zyg)dz -0.2684851 0.9930894
#> factor(zyg)mz  3.4239727 0.5773886
#> factor(zyg)os  0.4778091 0.7628390
#> cage           0.2519096 0.4821619
#> 
#> $or
#>               Estimate Std.Err     2.5%  97.5% P-value
#> factor(zyg)dz   0.7645  0.7593 -0.72357  2.253 0.31395
#> factor(zyg)mz  30.6911 17.7207 -4.04082 65.423 0.08328
#> factor(zyg)os   1.6125  1.2301 -0.79843  4.024 0.18989
#> cage            1.2865  0.6203  0.07073  2.502 0.03808
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

theta.des <- model.matrix( ~-1+factor(zyg)+factor(zyg)*cage,data=twinstut)
bina <- binomial_twostage(margbin,data=twinstut,var.link=1,
             clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bina)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link 
#> $estimates
#>                          theta        se
#> factor(zyg)dz      -0.27701246 1.0713974
#> factor(zyg)mz       3.52484148 0.5951743
#> factor(zyg)os       0.48859941 0.7644050
#> cage                0.06420907 3.6225641
#> factor(zyg)mz:cage  0.49441325 3.6865646
#> factor(zyg)os:cage -0.12312666 3.6921748
#> 
#> $or
#>                    Estimate Std.Err     2.5%  97.5% P-value
#> factor(zyg)dz        0.7580  0.8122  -0.8338  2.350 0.35063
#> factor(zyg)mz       33.9484 20.2052  -5.6531 73.550 0.09292
#> factor(zyg)os        1.6300  1.2460  -0.8121  4.072 0.19080
#> cage                 1.0663  3.8628  -6.5046  8.637 0.78251
#> factor(zyg)mz:cage   1.6395  6.0443 -10.2070 13.486 0.78619
#> factor(zyg)os:cage   0.8842  3.2644  -5.5140  7.282 0.78651
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

### use of clayton oakes binomial additive gamma model
###########################################################
 ## Reduce Ex.Timings
data <- sim_binClaytonOakes_family_ace(1000,2,1,beta=NULL,alpha=NULL)
margbin <- glm(ybin~x,data=data,family=binomial())
margbin
#> 
#> Call:  glm(formula = ybin ~ x, family = binomial(), data = data)
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.5001       0.2206  
#> 
#> Degrees of Freedom: 3999 Total (i.e. Null);  3998 Residual
#> Null Deviance:       5193 
#> Residual Deviance: 5182  AIC: 5186

head(data)
#>   ybin x   type cluster
#> 1    1 0 mother       1
#> 2    0 0 father       1
#> 3    1 0  child       1
#> 4    1 1  child       1
#> 5    0 1 mother       2
#> 6    0 1 father       2
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)

### make ace random effects design
out <- ace_family_design(data,member="type",id="cluster")
out$pardes
#>       [,1] [,2]
#>  [1,] 0.25    0
#>  [2,] 0.25    0
#>  [3,] 0.25    0
#>  [4,] 0.25    0
#>  [5,] 0.25    0
#>  [6,] 0.25    0
#>  [7,] 0.25    0
#>  [8,] 0.25    0
#>  [9,] 0.00    1
head(out$des.rv)
#>      m1 m2 m3 m4 f1 f2 f3 f4 env
#> [1,]  1  1  1  1  0  0  0  0   1
#> [2,]  0  0  0  0  1  1  1  1   1
#> [3,]  1  1  0  0  1  1  0  0   1
#> [4,]  1  0  1  0  1  0  1  0   1
#> [5,]  1  1  1  1  0  0  0  0   1
#> [6,]  0  0  0  0  1  1  1  1   1

bints <- binomial_twostage(margbin,data=data,
     clusters=data$cluster,detail=0,var.par=1,
     theta=c(2,1),var.link=0,
     random.design=out$des.rv,theta.des=out$pardes)
summary(bints)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> $estimates
#>                theta        se
#> dependence1 1.477113 0.4185208
#> dependence2 1.104406 0.1564451
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err   2.5%  97.5%   P-value
#> dependence1   0.5722 0.08722 0.4012 0.7431 5.364e-11
#> dependence2   0.4278 0.08722 0.2569 0.5988 9.337e-07
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1    2.582  0.3953 1.807 3.356 6.533e-11
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

data <- sim_binClaytonOakes_twin_ace(1000,2,1,beta=NULL,alpha=NULL)
out  <- twin_polygen_design(data,id="cluster",zygname="zygosity")
out$pardes
#>      [,1] [,2]
#> [1,]  1.0    0
#> [2,]  0.5    0
#> [3,]  0.5    0
#> [4,]  0.5    0
#> [5,]  0.0    1
head(out$des.rv)
#>   MZ DZ DZns1 DZns2 env
#> 1  1  0     0     0   1
#> 2  1  0     0     0   1
#> 3  1  0     0     0   1
#> 4  1  0     0     0   1
#> 5  1  0     0     0   1
#> 6  1  0     0     0   1
margbin <- glm(ybin~x,data=data,family=binomial())

bintwin <- binomial_twostage(margbin,data=data,
     clusters=data$cluster,var.par=1,
     theta=c(2,1),random.design=out$des.rv,theta.des=out$pardes)
summary(bintwin)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> $estimates
#>                 theta        se
#> dependence1 0.8423038 0.6254036
#> dependence2 1.5659818 0.4357446
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>             Estimate Std.Err     2.5%  97.5%  P-value
#> dependence1   0.3498  0.2223 -0.08598 0.7855 0.115669
#> dependence2   0.6502  0.2223  0.21451 1.0860 0.003446
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err 2.5% 97.5%   P-value
#> p1    2.408  0.3871 1.65 3.167 4.944e-10
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
concordanceTwinACE(bintwin)
#> $MZ
#>                      Estimate Std.Err   2.5%  97.5%    P-value
#> concordance            0.5131 0.01834 0.4771 0.5490 3.838e-172
#> casewise concordance   0.8090 0.01586 0.7779 0.8401  0.000e+00
#> marginal               0.6342 0.01691 0.6011 0.6674 6.370e-308
#> 
#> $DZ
#>                      Estimate Std.Err   2.5%  97.5%    P-value
#> concordance            0.4917 0.01935 0.4538 0.5296 1.884e-142
#> casewise concordance   0.7753 0.01742 0.7412 0.8094  0.000e+00
#> marginal               0.6342 0.01691 0.6011 0.6674 6.370e-308
#> 

```
