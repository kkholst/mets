# Analysis of multivariate survival data

## Overview

When looking at multivariate survival data with the aim of learning
about the dependence that is present, possibly after correcting for some
covariates different approaches are available in the mets package

- Binary models and adjust for censoring with inverse probabilty of
  censoring weighting
  - biprobit model
- Bivariate surival models of Clayton-Oakes type
  - With regression structure on dependence parameter
  - With additive gamma distributed random effects
  - Special functionality for polygenic random effects modelling such as
    ACE, ADE ,AE and so forth.
- Plackett OR model model
  - With regression structure on OR dependence parameter
- Cluster stratified Cox

Typically it can be hard or impossible to specify random effects models
with special structure among the parameters of the random effects. This
is possible for our specification of the random effects models.

To be concrete about the model structure assume that we have paired
survival data
$\left( T_{1},\delta_{1},T_{2},\delta_{2},X_{1},X_{2} \right)$ where the
censored survival responses are
$\left( T_{1},\delta_{1},T_{2},\delta_{2} \right)$ and the covariates
are $\left( X_{1},X_{2} \right)$.

The basic models assumes that each subject has a marginal on Cox-form
$$\lambda_{s{(k,i)}}(t)\exp\left( X_{ki}^{T}\beta \right)$$ where
$s(k,i)$ is a strata variable.

The constructed likelihood is a composite likehood based on

- all pairs within each cluster (when the pairs argument is not used)

- the specified pairs when the pairs argument is used.

In addition to the clusters specified for the construction of the
composite likelihood these can be added further summed using the
se.clusters argument that sums the influence functions over the
se.clusters. When se.clusters are not specified it is the same as the
cluster argument.

When the clusters of the dependence parametaters are the same as those
of the marginal model and the phreg function is used then the standard
errors are corrected for the uncertainty from the marginal models,
otherwise the returned standard errors are computed as if the marginals
are known.

## Gamma distributed frailties

The focus of this vignette is describe how to work on bivariate survival
data using the addtive gamma-random effects models. We present two
different ways of specifying different dependence structures.

- Univariate models with a single random effect for each cluster and
  with a regression design on the variance.

- Multivariate models with multiple random effects for each cluster.

The univariate models are then given a given cluster random effects
$Z_{k}$ with parameter $\theta$ the joint survival function is given by
the Clayton copula and on the form
$$\psi(\theta,\psi^{- 1}\left( \theta,S_{1}\left( t,X_{k1} \right) \right) + \psi^{- 1}\left( \theta,S_{1}\left( t,X_{k1} \right) \right)$$
where $\psi$ is the Laplace transform of a gamma distributed random
variable with mean 1 and variance $\theta$.

We then model the variance within clusters by a cluster specific
regression design such that $$\theta = h\left( z_{j}^{T}\alpha \right)$$
where $z$ is the regression design (specified by theta.des in the
software), and $h$ is link function, that is either $exp$ or the
identity.

This model can be fitted using a pairwise likelihood or the
pseudo-likelihood using either

- twostage

- twostageMLE

To make the twostage approach possible we need a model with specific
structure for the marginals. Therefore given the random effect of the
clusters the survival distributions within a cluster are independent and
on the form
$$P\left( T_{j} > t|X_{j},Z \right) = exp\left( - Z \cdot \Psi^{- 1}\left( \nu^{- 1},S\left( t|X_{j} \right) \right) \right)$$
with $\Psi$ the laplace of the gamma distribution with mean 1 and
variance $1/\nu$.

## Additive Gamma frailties

For the multivariate models we are given a multivarite random effect
each cluster $Z = \left( Z_{1},...,Z_{d} \right)$ with d random effects.
The total random effect for each subject $j$ in a cluster is then
specified using a regression design on these random effects, with a
regression vector $V_{j}$ such that the total random effect is
$V_{j}^{T}\left( Z_{1},...,Z_{d} \right)$. The elements of $V_{J}$ are
1/0. The random effects $\left( Z_{1},...,Z_{d} \right)$ has associated
parameters $\left( \lambda_{1},...,\lambda_{d} \right)$ and $Z_{j}$ is
Gamma distributed with

- mean $\lambda_{j}/V_{1}^{T}\lambda$

- variance $\lambda_{j}/\left( V_{1}^{T}\lambda \right)^{2}$

The key assumption to make the two-stage fitting possible is that
$$\begin{array}{r}
{\nu = V_{j}^{T}\lambda}
\end{array}$$ is constant within clusters. The consequence of this is
that the total random effect for each subject within a cluster,
$V_{j}^{T}\left( Z_{1},...,Z_{d} \right)$, is gamma distributed with
variance $1/\nu$.

The DEFAULT parametrization (var.par=1) uses the variances of the random
effecs $$\begin{array}{r}
{\theta_{j} = \lambda_{j}/\nu^{2}}
\end{array}$$ For alternative parametrizations one can specify that the
parameters are $\theta_{j} = \lambda_{j}$ with the argument var.par=0.

Finally the parameters $\left( \theta_{1},...,\theta_{d} \right)$ are
related to the parameters of the model by a regression construction $M$
(d x k), that links the $d$$\theta$ parameters with the $k$ underlying
$\alpha$ parameters $$\begin{aligned}
\theta & {= M\alpha.}
\end{aligned}$$ The default is a diagonal matrix for $M$. This can be
used to make structural assumptions about the variances of the
random-effects as is needed for the ACE model for example. In the
software \$ M \$ is called theta.des

Assume that the marginal survival distribution for subject $i$ within
cluster $k$ is given by $S_{X_{k,i}}(t)$ given covariates $X_{k,i}$.

Now given the random effects of the cluster $Z_{k}$ and the
covariates$X_{k,i}$$i = 1,\ldots,n_{k}$ we assume that subjects within
the cluster are independent with survival distributions
$$\begin{array}{r}
{\exp\left( - \left( V_{k,i}Z_{k} \right)\Psi^{- 1}\left( \nu,S_{X_{k,i}}(t) \right) \right).}
\end{array}$$

A consequence of this is that the hazards given the covariates $X_{k,i}$
and the random effects $Z_{k}$ are given by $$\begin{array}{r}
{\lambda_{k,i}\left( t;X_{k,i},Z_{k,i} \right) = \left( V_{k,i}V_{k} \right)D_{3}\Psi^{- 1}\left( \nu,S_{X_{k,i}}(t) \right)D_{t}S_{X_{k,i}}(t)}
\end{array}$$ where $D_{t}$ and $D_{3}$ denotes the partial derivatives
with respect to $t$ and the third argument, respectively.

Further, we can express the multivariate survival distribution as
$$\begin{aligned}
{S\left( t_{1},\ldots,t_{m} \right)} & {= \exp\left( - \sum\limits_{i = 1}^{m}\left( V_{i}Z \right)\Psi^{- 1}\left( \eta_{l},\nu_{l},S_{X_{k,i}}\left( t_{i} \right) \right) \right)} \\
 & {= \prod\limits_{l = 1}^{p}\Psi\left( \eta_{l},\eta,\sum\limits_{i = 1}^{m}Q_{k,i}\Psi^{- 1}\left( \eta,\eta,S_{X_{k,i}}\left( t_{i} \right) \right) \right).}
\end{aligned}$$ In the case of considering just pairs, we write this
function as $C\left( S_{k,i}(t),S_{k,j}(t) \right)$.

In addition to survival times from this model, we assume that we
independent right censoring present $U_{k,i}$ such that the given
$V_{k}$ and the
covariates$X_{k,i}$$i = 1,\ldots,n_{k}$$\left( U_{k,1},\ldots,U_{k,n_{k}} \right)$
of $\left( T_{k,1},\ldots,T_{k,n_{k}} \right)$, and the conditional
censoring distribution do not depend on $V_{k}$.

One consequence of the model strucure is that the Kendall’s can be
computed for two-subjects $(i,j)$ across two clusters `1'' and`2’’ as
$$\begin{array}{r}
{E\left( \frac{\left( V_{1i}Z_{1} - V_{1j}Z_{2} \right)\left( V_{2i}Z_{1} - V_{2j}Z_{2} \right)}{\left( V_{1i}Z_{1} + V_{2i}Z_{2} \right)\left( V_{1j}Z_{1} + V_{2j}Z_{2} \right)} \right)}
\end{array}$$ under the assumption that that we compare pairs with
equivalent marginals, $S_{X_{1,i}}(t) = S_{X_{2,i}}(t)$ and
$S_{X_{1,j}}(t) = S_{X_{2,j}}(t)$, and that
$S_{X_{1,i}}(\infty) = S_{X_{1,j}}(\infty) = 0$. Here we also use that
$\eta$ is the same across clusters. The Kendall’s tau would be the same
for due to the same additive structure for the frailty terms, and the
random effects thus have the same interpretation in terms of Kendall’s
tau.

## Univariate gamma (clayton-oakes) model twostage models

We start by fitting simple Clayton-Oakes models for the data, that is
with an overall random effect that is Gamma distrubuted with variance
$\theta$. We can fit the model by a pseudo-MLE (twostageMLE) and a
pairwise composite likelihood approach (twostage).

The pseudo-liklihood and the composite pairwise likelhood should give
the same for this model since we have paired data. In addition the
log-parametrization is illustrated with the var.link=1 option. In
addition it is specified that we want a “clayton.oakes” model. We note
that the standard errors differs because the twostage does not include
the variance due to the baseline parameters for this type of modelling,
so here it is better to use the twostageMLE.

``` r
 library(mets)
 data(diabetes)
 set.seed(100)
 
 # Marginal Cox model  with treat as covariate
 margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
 # Clayton-Oakes, MLE 
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
 
 # Clayton-Oakes
 fitco2 <- survival.twostage(margph,data=diabetes,theta=0.0,
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
 fitco3 <- survival.twostage(margph,data=diabetes,theta=1.0,
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
```

Note, the standard errors are slightly different when comparing fitco1
with fitco3 since the survival.twostage uses numerical derivatives for
the hessian and the derivative in the direction of the marginal model.

The marginal models can be either structured Cox model or as here with a
baseline for each strata. This gives quite similar results to those
before.

``` r
  # without covariates but marginal model stratified 
  marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
 fitco<-twostageMLE(marg,data=diabetes,theta=1.0)
 summary(fitco)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z       P-val Kendall tau         SE
#> dependence1 0.9447446 0.3516204 2.686831 0.007213349    0.320824 0.08109776
#> 
#> $type
#> NULL
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

  fitcoa <- survival.twostage(marg,data=diabetes,theta=1.0,clusters=diabetes$id,
           model="clayton.oakes",var.link=0)
  summary(fitcoa)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z       P-val Kendall tau         SE
#> dependence1 0.9447446 0.3516129 2.686889 0.007212097    0.320824 0.08109601
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

## Piecewise constant Clayton-Oakes model

Let the cross-hazard ratio (CHR) be defined as $$\begin{array}{r}
{\eta\left( t_{1},t_{2} \right) = \frac{\lambda_{1}\left( t_{1}|T_{2} = t_{2} \right)}{\lambda_{1}\left( t_{1}|T_{2} \geq t_{2} \right)} = \frac{\lambda_{2}\left( t_{2}|T_{1} = t_{1} \right)}{\lambda_{2}\left( t_{2}|T_{1} \geq t_{1} \right)}}
\end{array}$$ where $\lambda_{1}$ and $\lambda_{2}$ are the conditional
hazard functions of $T_{1}$ and $T_{2}$ given covariates. For the
Clayton-Oakes model this ratio is
$\eta\left( t_{1},t_{2} \right) = 1 + \theta$, and as a consequence we
see that if the co-twin is dead at any time we would increase our risk
assessment on the hazard scale with the constant
$\eta\left( t_{1},t_{2} \right)$. The Clayton-Oakes model also has the
nice property that Kendall’s tau is linked directly to the dependence
parameter $\theta$ and is $1/(1 + 2/\theta)$.

A very useful extension of the model the constant cross-hazard ratio
(CHR) model is the piecewise constant cross-hazard ratio (CHR) for
bivariate survival data , and this model was extended to competing risks
in .

In the survival setting we let the CHR $$\begin{aligned}
{\eta\left( t_{1},t_{2} \right)} & {= \sum\eta_{i,j}I\left( t_{1} \in I_{i},t_{2} \in I_{j} \right)}
\end{aligned}$$

The model lets the CHR by constant in different part of the plane. This
can be thought of also as having a separate Clayton-Oakes model for each
of the regions specified in the plane here by the cut-points
$c(0,0.5,2)$ thus defining 9 regions.

This provides a constructive goodness of fit test for the whether the
Clayton-Oakes model is valid. Indeed if valid the parameter should be
the same in all regions.

First we generate some data from the Clayton-Oakes model with variance
$0.5$ and 2000 pairs. And fit the related model.

``` r
 d <- simClaytonOakes(200,2,0.5,0,3)
  margph <- phreg(Surv(time,status)~x+cluster(cluster),data=d)
 # Clayton-Oakes, MLE 
 fitco1<-twostageMLE(margph,data=d)
 summary(fitco1)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                Coef.        SE       z        P-val Kendall tau         SE
#> dependence1 2.283394 0.3164797 7.21498 5.393463e-13   0.5330806 0.03449846
#> 
#> $type
#> NULL
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

Now we cut the region at the cut-points $c(0,0.5,2)$ thus defining 9
regions and fit a separate model for each region.  
We see that the parameter is indeed rather constant over the 9 regions.
A formal test can be constructed.

``` r
 udp <- piecewise.twostage(c(0,0.5,2),data=d,id="cluster",timevar="time",status="status",model="clayton.oakes",silent=0)
#> Data-set  1 out of  4
#>   Number of joint events: 55 of  200
#> Data-set  2 out of  4
#>   Number of joint events: 25 of  116
#> Data-set  3 out of  4
#>   Number of joint events: 28 of  120
#> Data-set  4 out of  4
#>   Number of joint events: 50 of  91
 summary(udp)
#> [1] 1
#> Dependence parameter for Clayton-Oakes model 
#> log-coefficient for dependence parameter (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  0.566 (0.222)  1.278 (0.223)
#> 0.5 - 2  0.519 (0.268)  1.092 (0.205)
#> 
#> Kendall's tau (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  0.468 (0.055)  0.642 (0.051)
#> 0.5 - 2  0.457 (0.066)  0.598 (0.049)
```

## Multivariate gamma twostage models

To illustrate how the multivariate models can be used, we first set up
some twin data with ACE structure. That is two shared random effects,
one being the genes $\sigma_{g}^{2}$ and one the environmental effect
$\sigma_{e}^{2}$. Monozygotic twins share all genes whereas the
dizygotic twins only share half the genes. This can be expressed via 5
random effect for each twin pair (for example). We start by setting this
up.

The pardes matrix tells how the the parameters of the 5 random effects
are related, and the matrix her first has one random effect with
parameter $\theta_{1}$ (here the $\sigma_{g}^{2}$ ), then the next 3
random effects have parameters $0.5\theta_{1}$ (here $0.5\sigma_{g}^{2}$
), and the last random effect that is given by its own parameter
$\theta_{2}$ (here $\sigma_{e}^{2}$ ).

``` r
 data <- simClaytonOakes.twin.ace(200,2,1,0,3)

 out <- twin.polygen.design(data,id="cluster",zyg="DZ",zygname="zyg",type="ace")
 pardes <- out$pardes
 pardes 
#>      [,1] [,2]
#> [1,]  1.0    0
#> [2,]  0.5    0
#> [3,]  0.5    0
#> [4,]  0.5    0
#> [5,]  0.0    1
```

The last part of the model structure is to decide how the random effects
are shared for the different pairs (MZ and DZ), this is specfied by the
random effects design ($V_{1}$ and $V_{2}$) for each pair. This is here
specified by an overall designmatrix for each subject (since they enter
all pairs with the same random effects design).

For an MZ pair the two share the full gene random effect and the full
environmental random effect. In contrast the DZ pairs share the 2nd
random effect with half the gene-variance and have both a non-shared
gene-random effect with half the variance, and finally a fully shared
environmental random effect.

``` r
 des.rv <- out$des.rv
 # MZ
 head(des.rv,2)
#>   MZ DZ DZns1 DZns2 env
#> 1  1  0     0     0   1
#> 2  1  0     0     0   1
 # DZ 
 tail(des.rv,2)
#>     MZ DZ DZns1 DZns2 env
#> 399  0  1     1     0   1
#> 400  0  1     0     1   1
```

Now we call the twostage function. We see that we essentially recover
the true values, and note that the output also compares the sizes of the
genetic and environmental random effect. This number is sometimes called
the heritability. In addition the total variance for each subject is
also computed and is here around $3$, as we indeed constructed.

``` r
### data <- simClaytonOakes.twin.ace(2000,2,1,0,3)
### out <- twin.polygen.design(data,id="cluster",zyg="DZ",zygname="zyg",type="ace")
 aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data)
 ts <- twostage(aa,data=data,clusters=data$cluster,
      theta=c(2,1),var.link=0,random.design=out$des.rv,theta.des=out$pardes)
 summary(ts)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                Coef.        SE        z       P-val Kendall tau         SE
#> dependence1 2.051732 0.6981882 2.938652 0.003296431   0.5063840 0.08505916
#> dependence2 1.025707 0.5674045 1.807717 0.070650519   0.3389974 0.12395643
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err     2.5%  97.5%   P-value
#> [1,]   0.6667    0.19  0.29431 1.0391 0.0004498
#> [2,]   0.3333    0.19 -0.03909 0.7057 0.0793902
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1    3.077  0.3982 2.297 3.858 1.092e-14
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

- A nice feature of the procdure is that it scales linearly in the
  number of observations
  - 1 mill pairs had a running time of around 100 seconds.

``` r
run <- 0
if (run==1) {
 data <- simClaytonOakes.twin.ace(1000000,2,1,0,3)

 out <- twin.polygen.design(data,id="cluster",zyg="DZ",zygname="zyg",type="ace")
 pardes <- out$pardes
 aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data)
 system.time(
 ts <- twostage(aa,data=data,clusters=data$cluster,
      theta=c(2,1),var.link=0,random.design=out$des.rv,theta.des=out$pardes)
 )
 summary(ts)
}
```

The estimates can be transformed into Kendall’s tau estimates for MZ and
DZ twins. The Kendall’s tau in the above output reflects how a gamma
distributed random effect in the normal Clayton-Oakes model is related
to the Kendall’s tau. In this setting the Kendall’s of MZ and DZ,
however, should reflect both random effects.

We do this based on simulations. The Kendall’s tau of the MZ is around
0.60, and for DZ around 0.33. Both are quite high and this is due to a
large shared environmental effect and large genetic effect.

``` r
kendall.ClaytonOakes.twin.ace(ts$theta[1],ts$theta[2],K=10000) 
#> $mz.kendall
#> [1] 0.6186546
#> 
#> $dz.kendall
#> [1] 0.3350827
```

## Family data

For family data, things are quite similar since we use only the pairwise
structure. We show how the designs are specified.

First we simulate data from an ACE model. 2000 families with two-parents
that share only the environment, and two-children that share genes with
their parents.

``` r
library(mets)
set.seed(1000)
data <- simClaytonOakes.family.ace(200,2,1,0,3)
head(data)
#>        time status x cluster   type   mintime lefttime truncated
#> 1 0.4139376      1 0       1 mother 0.4139376        0         0
#> 2 1.6435053      1 1       1 father 0.4139376        0         0
#> 3 1.2485275      1 0       1  child 0.4139376        0         0
#> 4 1.1079118      1 1       1  child 0.4139376        0         0
#> 5 3.0000000      0 0       2 mother 1.6718726        0         0
#> 6 2.0085566      1 1       2 father 1.6718726        0         0
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)
```

To set up the random effects some functions can be used. We here set up
the ACE model that has 9 random effects with one shared environmental
effect (the last random effect) and 4 genetic random effects for each
parent, with variance $\sigma_{g}^{2}/4$.

The random effect is again set-up with an overall designmatrix because
it is again the same for each subject for all comparisons across family
members. We below demonstrate how the model can be specified in various
other ways.

Each child share 2 genetic random effects with each parent, and also
share 2 genetic random effects with his/her sibling.

``` r
out <- ace.family.design(data,member="type",id="cluster")
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
head(out$des.rv,4)
#>      m1 m2 m3 m4 f1 f2 f3 f4 env
#> [1,]  1  1  1  1  0  0  0  0   1
#> [2,]  0  0  0  0  1  1  1  1   1
#> [3,]  1  1  0  0  1  1  0  0   1
#> [4,]  1  0  1  0  1  0  1  0   1
```

Then we fit the model

``` r
pa <- phreg(Surv(time,status)~+1+cluster(cluster),data=data)

# make ace random effects design
ts <- twostage(pa,data=data,clusters=data$cluster,var.par=1,var.link=0,theta=c(2,1),
        random.design=out$des.rv,theta.des=out$pardes)
summary(ts)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z        P-val Kendall tau         SE
#> dependence1 2.0094622 0.5369170 3.742594 0.0001821303   0.5011800 0.06679822
#> dependence2 0.6503464 0.2136795 3.043561 0.0023379628   0.2453817 0.06083976
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err    2.5% 97.5%   P-value
#> [1,]   0.7555 0.09823 0.56297 0.948 1.458e-14
#> [2,]   0.2445 0.09823 0.05199 0.437 1.280e-02
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1     2.66  0.4459 1.786 3.534 2.437e-09
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

The model can also be fitted by specifying the pairs that one wants for
the pairwise likelhood. This is done by specifying the pairs argument.
We start by considering all pairs as we also did before.

All pairs can be written up by calling the familycluster.index function.

There are xx pairs to consider, and the first 6 pairs for the first
family is written out here.

``` r
# now specify fitting via specific pairs 
# first all pairs 
mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=12)
#>  [1] 1 2 1 3 1 4 2 3 2 4 3 4
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
head(pairs,n=6)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    3
#> [3,]    1    4
#> [4,]    2    3
#> [5,]    2    4
#> [6,]    3    4
```

Then fitting the model using only specified pairs

``` r
ts <- twostage(pa,data=data,clusters=data$cluster, theta=c(2,1),var.link=0,step=1.0,
        random.design=out$des.rv, theta.des=out$pardes,pairs=pairs)
summary(ts)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z        P-val Kendall tau         SE
#> dependence1 2.0094622 0.5369170 3.742594 0.0001821303   0.5011800 0.06679822
#> dependence2 0.6503464 0.2136795 3.043561 0.0023379628   0.2453817 0.06083976
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err    2.5% 97.5%   P-value
#> [1,]   0.7555 0.09823 0.56297 0.948 1.458e-14
#> [2,]   0.2445 0.09823 0.05199 0.437 1.280e-02
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1     2.66  0.4459 1.786 3.534 2.437e-09
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

Now we only use a random sample of the pairs by sampling these. The
pairs picked still refers to the data given in the data argument, and
clusters (families) are also specified as before.

``` r
ssid <- sort(sample(1:nrow(pairs),200))
tsd <- twostage(pa,data=data,clusters=data$cluster,
    theta=c(2,1)/10,var.link=0,random.design=out$des.rv,
   theta.des=out$pardes,pairs=pairs[ssid,])
summary(tsd)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z     P-val Kendall tau         SE
#> dependence1 0.9237395 0.6070493 1.521688 0.1280873   0.3159445 0.14202885
#> dependence2 0.4087673 0.2660432 1.536469 0.1244233   0.1696998 0.09170489
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err    2.5%  97.5% P-value
#> [1,]   0.6932  0.2523  0.1987 1.1878 0.00601
#> [2,]   0.3068  0.2523 -0.1878 0.8013 0.22410
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err   2.5% 97.5% P-value
#> p1    1.333  0.4803 0.3912 2.274 0.00553
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

Sometimes one only has the data from the pairs in addition to for
example a cohort estimate of the marginal surival models. We now
demonstrate how this is dealt with. Everything is essentially as before
but need to organize the design differently compared to before we
specified the design  
for everybody in the cohort. In addition we do not here bring in the
uncertainty from the baseline in the estimates, even though this is
formally possible, but when the data of the marginal model and twostage
data are not the same, we have to specify that we do not want the
decomposition for the uncertainty due to the baseline (baseline.iid=0).

``` r
ids <- sort(unique(c(pairs[ssid,])))

pairsids <- c(pairs[ssid,])
pair.new <- matrix(fast.approx(ids,c(pairs[ssid,])),ncol=2)
head(pair.new)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4
#> [3,]    3    6
#> [4,]    5    6
#> [5,]    7    8
#> [6,]    9   11

# this requires that pair.new refers to id's in dataid (survival, status and so forth)
# random.design and theta.des are constructed to be the array 3 dims via individual specfication from ace.family.design
dataid <- dsort(data[ids,],"cluster")
outid <- ace.family.design(dataid,member="type",id="cluster")
outid$pardes
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
head(outid$des.rv)
#>      m1 m2 m3 m4 f1 f2 f3 f4 env
#> [1,]  0  0  0  0  1  1  1  1   1
#> [2,]  1  0  1  0  1  0  1  0   1
#> [3,]  1  1  1  1  0  0  0  0   1
#> [4,]  0  0  0  0  1  1  1  1   1
#> [5,]  1  1  0  0  1  1  0  0   1
#> [6,]  1  0  1  0  1  0  1  0   1
```

Now fitting the model using only the pair data.

``` r
tsdid <- twostage(pa,data=dataid,clusters=dataid$cluster,theta=c(2,1)/10,var.link=0,baseline.iid=0,
          random.design=outid$des.rv,theta.des=outid$pardes,pairs=pair.new)
summary(tsdid)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z     P-val Kendall tau         SE
#> dependence1 0.9237395 0.6040859 1.529153 0.1262266   0.3159445 0.14133553
#> dependence2 0.4087673 0.2590853 1.577732 0.1146271   0.1696998 0.08930652
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err    2.5%  97.5%  P-value
#> [1,]   0.6932  0.2487  0.2058 1.1807 0.005311
#> [2,]   0.3068  0.2487 -0.1807 0.7942 0.217376
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err   2.5% 97.5%  P-value
#> p1    1.333  0.4784 0.3949  2.27 0.005345
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"

paid <- phreg(Surv(time,status)~+1+cluster(cluster),data=dataid)
tsdidb <- twostage(paid,data=dataid,clusters=dataid$cluster,theta=c(2,1)/10,
  var.link=0,random.design=outid$des.rv,theta.des=outid$pardes,pairs=pair.new)
summary(tsdidb)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z     P-val Kendall tau         SE
#> dependence1 0.9405413 0.6151098 1.529062 0.1262490   0.3198531 0.14227484
#> dependence2 0.3381704 0.2557053 1.322501 0.1860015   0.1446303 0.09354431
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err    2.5%  97.5%  P-value
#> [1,]   0.7355  0.2499  0.2458 1.2253 0.003247
#> [2,]   0.2645  0.2499 -0.2253 0.7542 0.289925
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err   2.5% 97.5%  P-value
#> p1    1.279  0.4865 0.3251 2.232 0.008583
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
coef(tsdid)
#>                 Coef.        SE        z     P-val Kendall tau         SE
#> dependence1 0.9237395 0.6040859 1.529153 0.1262266   0.3159445 0.14133553
#> dependence2 0.4087673 0.2590853 1.577732 0.1146271   0.1696998 0.08930652
```

Estimates changed because we used either the marginal from the
full-data, in which case the standard errors did not reflect the
uncertainty from the baseline, or the marginal estimated from only the
sub-sample in which case the marginals were slightly different.

### Pairwise specification of random effects and variances

Now we illustrate how one can also directly specify the random.design
and theta.design for each pair, rather than taking an overall
specification that can be used for the whole family via the rows of the
des.rv for the relevant pairs. This can be much simpler in some
situations.

``` r
pair.types <-  matrix(dataid[c(t(pair.new)),"type"],byrow=T,ncol=2)
head(pair.new)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4
#> [3,]    3    6
#> [4,]    5    6
#> [5,]    7    8
#> [6,]    9   11
head(pair.types)
#>      [,1]     [,2]    
#> [1,] "father" "child" 
#> [2,] "mother" "father"
#> [3,] "mother" "child" 
#> [4,] "child"  "child" 
#> [5,] "father" "child" 
#> [6,] "mother" "child"

theta.des  <- rbind( c(rbind(c(1,0),c(1,0),c(0,1),c(0,0))),
        c(rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))))
random.des <- rbind( 
        c(1,0,1,0),c(0,1,1,0),
        c(1,1,0,1),c(1,0,1,1))
mf <- 1*(pair.types[,1]=="mother" & pair.types[,2]=="father")
##          pair, rv related to pairs,  theta.des related to pair 
pairs.new <- cbind(pair.new,(mf==1)*1+(mf==0)*3,(mf==1)*2+(mf==0)*4,(mf==1)*1+(mf==0)*2,(mf==1)*3+(mf==0)*4)
```

pairs.new is matix with

- columns 1:2 giving the indeces of the data points

- columns 3:4 giving the indeces of the random.design for the different
  pairs

- columns 5 giving the indeces of the theta.des written as rows

- columns 6 giving the number of random variables for this pair

Looking at the first three rows. We see that the composite likehood is
based on data-points (1,2), (3,4) and (5,6), these are (mother, father),
(mother, child), and (father, child), respectively.

``` r
head(pairs.new[1:3,])
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    2    3    4    2    4
#> [2,]    3    4    1    2    1    3
#> [3,]    3    6    3    4    2    4
head(dataid)
#>       time status x cluster   type   mintime lefttime truncated number child
#> 2 1.643505      1 1       1 father 0.4139376        0         0      2     0
#> 4 1.107912      1 1       1  child 0.4139376        0         0      4     0
#> 5 3.000000      0 0       2 mother 1.6718726        0         0      1     0
#> 6 2.008557      1 1       2 father 1.6718726        0         0      2     0
#> 7 1.671873      1 1       2  child 1.6718726        0         0      3     1
#> 8 1.763022      1 1       2  child 1.6718726        0         0      4     0
```

The random effects for these are specified from random effects with
design read from the random.design, using the rows (1,2), (3,4) and
(3,4), respecively, and with random effects that have variances given by
theta.des rows, 1,2, and 2 respectively in the three cases. For the
first pair (1,2), the random vectors and their variances are given by,
(mother, father) pair,

``` r
random.des[1,]
#> [1] 1 0 1 0
random.des[2,]
#> [1] 0 1 1 0
matrix(theta.des[1,],4,2)
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    0
#> [3,]    0    1
#> [4,]    0    0
```

thus sharing only the third random effect with variance $\sigma_{e}^{2}$
and having two non-shared random effects with variances
$\sigma_{g}^{2}$, and finally a last 4th random effect with variance $0$
that thus could have been omitted.

The length of all rows of theta.des are the maximum number of random
effects $\times$ the number of parameters. These two numbers are given
in the call. In this case 4 $\times$ 2. So theta.des has rows of length
$8$, possibly including some 0’s for rows not relevant due to fewer
random effects, as is the case here for pairs that do not share genetic
effects.

Now considering the parent and their child, they are thus sharing the
first random effect with variance $0.5\sigma_{g}^{2}$ then there are two
non-shared random effects with variances $0.5\sigma_{g}^{2}$, and
finally a shared environment with variance $\sigma_{e}^{2}$.

``` r
head(dataid)
#>       time status x cluster   type   mintime lefttime truncated number child
#> 2 1.643505      1 1       1 father 0.4139376        0         0      2     0
#> 4 1.107912      1 1       1  child 0.4139376        0         0      4     0
#> 5 3.000000      0 0       2 mother 1.6718726        0         0      1     0
#> 6 2.008557      1 1       2 father 1.6718726        0         0      2     0
#> 7 1.671873      1 1       2  child 1.6718726        0         0      3     1
#> 8 1.763022      1 1       2  child 1.6718726        0         0      4     0
matrix(theta.des[2,],4,2)
#>      [,1] [,2]
#> [1,]  0.5    0
#> [2,]  0.5    0
#> [3,]  0.5    0
#> [4,]  0.0    1
random.des[3,]
#> [1] 1 1 0 1
random.des[4,]
#> [1] 1 0 1 1
```

And fitting again the same model as before

``` r
tsdid2 <- twostage(pa,data=dataid,clusters=dataid$cluster,
       theta=c(2,1)/10,var.link=0,step=1.0,random.design=random.des,
       baseline.iid=0, theta.des=theta.des,pairs=pairs.new,dim.theta=2)
summary(tsdid2)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z     P-val Kendall tau         SE
#> dependence1 0.9237395 0.6040859 1.529153 0.1262266   0.3159445 0.14133553
#> dependence2 0.4087673 0.2590853 1.577732 0.1146271   0.1696998 0.08930652
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err    2.5%  97.5%  P-value
#> [1,]   0.6932  0.2487  0.2058 1.1807 0.005311
#> [2,]   0.3068  0.2487 -0.1807 0.7942 0.217376
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err   2.5% 97.5%  P-value
#> p1    1.333  0.4784 0.3949  2.27 0.005345
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
tsd$theta
#>           [,1]
#> [1,] 0.9237395
#> [2,] 0.4087673
tsdid2$theta
#>           [,1]
#> [1,] 0.9237395
#> [2,] 0.4087673
tsdid$theta
#>           [,1]
#> [1,] 0.9237395
#> [2,] 0.4087673
```

Finally the same model structure can be setup based on a Kinship
coefficient.

``` r
kinship  <- rep(0.5,nrow(pair.types))
kinship[pair.types[,1]=="mother" & pair.types[,2]=="father"] <- 0
head(kinship,n=10)
#>  [1] 0.5 0.0 0.5 0.5 0.5 0.5 0.5 0.0 0.5 0.0

out <- make.pairwise.design(pair.new,kinship,type="ace") 
```

Same output as before

``` r
tsdid3 <- twostage(pa,data=dataid,clusters=dataid$cluster,
   theta=c(2,1)/10,var.link=0,step=1.0,random.design=out$random.design,
   baseline.iid=0,theta.des=out$theta.des,pairs=out$new.pairs,dim.theta=2)
summary(tsdid3)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.        SE        z     P-val Kendall tau         SE
#> dependence1 0.9237395 0.6040859 1.529153 0.1262266   0.3159445 0.14133553
#> dependence2 0.4087673 0.2590853 1.577732 0.1146271   0.1696998 0.08930652
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err    2.5%  97.5%  P-value
#> [1,]   0.6932  0.2487  0.2058 1.1807 0.005311
#> [2,]   0.3068  0.2487 -0.1807 0.7942 0.217376
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err   2.5% 97.5%  P-value
#> p1    1.333  0.4784 0.3949  2.27 0.005345
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
tsdid2$theta
#>           [,1]
#> [1,] 0.9237395
#> [2,] 0.4087673
tsdid$theta
#>           [,1]
#> [1,] 0.9237395
#> [2,] 0.4087673
```

Now fitting the AE model without the “C” component for shared
environment:

``` r
outae <- make.pairwise.design(pair.new,kinship,type="ae") 
tsdid4 <- twostage(pa,data=dataid,clusters=dataid$cluster,
   theta=c(2,1)/10,var.link=0,random.design=outae$random.design,
   baseline.iid=0,theta.des=outae$theta.des,pairs=outae$new.pairs,dim.theta=1)
summary(tsdid4)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                Coef.        SE        z        P-val Kendall tau        SE
#> dependence1 1.932485 0.5702673 3.388735 0.0007021583   0.4914157 0.0737521
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> $h
#>      Estimate Std.Err 2.5% 97.5% P-value
#> [1,]        1       0    1     1       0
#> 
#> $vare
#> NULL
#> 
#> $vartot
#>    Estimate Std.Err   2.5% 97.5%   P-value
#> p1    1.932  0.5703 0.8148  3.05 0.0007022
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

### Pairwise dependence modelling

We now illustate how to estimate all pairwise associations between
different family-members using the twostage function. They key is
specify the pairs for the composite likelihood directly and then the
associated design-matrix directly. This needs to be done looking at both
subjects in the pairs, and the design therefore follows the pairs and is
a matrix where its row is given as the third column in the pairs
argument.

First we get the pairs to be considered

``` r
library(mets)
set.seed(1000)
data <- simClaytonOakes.family.ace(200,2,1,0,3)
head(data)
#>        time status x cluster   type   mintime lefttime truncated
#> 1 0.4139376      1 0       1 mother 0.4139376        0         0
#> 2 1.6435053      1 1       1 father 0.4139376        0         0
#> 3 1.2485275      1 0       1  child 0.4139376        0         0
#> 4 1.1079118      1 1       1  child 0.4139376        0         0
#> 5 3.0000000      0 0       2 mother 1.6718726        0         0
#> 6 2.0085566      1 1       2 father 1.6718726        0         0
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)

mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=20)
#>  [1] 1 2 1 3 1 4 2 3 2 4 3 4 5 6 5 7 5 8 6 7
pairs <- mm$pairs
dim(pairs)
#> [1] 1200    2
head(pairs,12)
#>       [,1] [,2]
#>  [1,]    1    2
#>  [2,]    1    3
#>  [3,]    1    4
#>  [4,]    2    3
#>  [5,]    2    4
#>  [6,]    3    4
#>  [7,]    5    6
#>  [8,]    5    7
#>  [9,]    5    8
#> [10,]    6    7
#> [11,]    6    8
#> [12,]    7    8
```

Now we construct the design matrix related to the pairs

``` r
 dtypes <- interaction( data[pairs[,1],"type"], data[pairs[,2],"type"])
 dtypes <- droplevels(dtypes)
 table(dtypes)
#> dtypes
#>   child.child  father.child  mother.child mother.father 
#>           200           400           400           200
 dm <- model.matrix(~-1+factor(dtypes))
 head(dm)
#>   factor(dtypes)child.child factor(dtypes)father.child
#> 1                         0                          0
#> 2                         0                          0
#> 3                         0                          0
#> 4                         0                          1
#> 5                         0                          1
#> 6                         1                          0
#>   factor(dtypes)mother.child factor(dtypes)mother.father
#> 1                          0                           1
#> 2                          1                           0
#> 3                          1                           0
#> 4                          0                           0
#> 5                          0                           0
#> 6                          0                           0
```

Then we fit the model:

``` r
pa <- phreg(Surv(time,status)~cluster(cluster),data)

tsp <- twostage(pa,data=data,theta.des=dm,pairs=cbind(pairs,1:nrow(dm)),se.clusters=data$clust)
summary(tsp)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects 
#> With log-link
#> $estimates
#>                               log-Coef.        SE           z       P-val
#> factor(dtypes)child.child    0.01865073 0.2087690  0.08933666 0.928814359
#> factor(dtypes)father.child  -0.29471237 0.2167688 -1.35956992 0.173966062
#> factor(dtypes)mother.child  -0.24879997 0.1949185 -1.27643105 0.201803188
#> factor(dtypes)mother.father -1.42550068 0.5159729 -2.76274312 0.005731786
#>                             Kendall tau         SE
#> factor(dtypes)child.child     0.3374907 0.04667882
#> factor(dtypes)father.child    0.2713351 0.04285787
#> factor(dtypes)mother.child    0.2805072 0.03933901
#> factor(dtypes)mother.father   0.1072975 0.04942234
#> 
#> $vargam
#>             Estimate Std.Err      2.5%  97.5%   P-value
#> dependence1   1.0188  0.2127  0.601943 1.4357 1.668e-06
#> dependence2   0.7447  0.1614  0.428334 1.0612 3.965e-06
#> dependence3   0.7797  0.1520  0.481851 1.0776 2.892e-07
#> dependence4   0.2404  0.1240 -0.002714 0.4835 5.261e-02
#> 
#> $type
#> [1] "clayton.oakes"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

We note that mother-father have the smallest correlation since they only
share the environmental random effect, whereas other pairs have a
similar correlation (in fact the same due to the way the data was
simulated) consisting of half the genetic effect and the environmental
effect that was also shared.

## Univariate plackett model twostage models

The copula known as the Plackett distribution, see , is on the form
$$\begin{array}{r}
{C(u,v;\theta) = \begin{cases}
\frac{S - \left( S^{2} - 4uv\theta(\theta - a) \right)}{2(\theta - 1)} & {{\mspace{6mu}\text{if}\mspace{6mu}}\theta \neq 1} \\
{uv} & {{\mspace{6mu}\text{if}\mspace{6mu}}\theta = 1}
\end{cases}}
\end{array}$$ with $S = 1 + (\theta - 1)(u + v)$. With marginals $S_{i}$
we now define the bivariate survival function as
$C\left( u_{1},u_{2} \right) = H\left( S_{1}\left( t_{1} \right),S_{2}\left( t_{2} \right) \right)$
with $u_{i} = S_{i}\left( t_{i} \right)$.

The dependence parameter $\theta$ has the nice interpretation that the
it is equivalent to the odds-ratio of all $2 \times 2$ tables for
surviving past any cut of the plane $\left( t_{1},t_{2} \right)$, that
is
$$\theta = \frac{P\left( T_{1} > t_{1}|T_{2} > t_{2} \right)P\left( T_{1} \leq t_{1}|T_{2} > t_{2} \right)}{P\left( T_{1} > t_{1}|T_{2} \leq t_{2} \right)P\left( T_{1} \leq t_{1}|T_{2} \leq t_{2} \right)}.$$

One additional nice feature of the odds-ratio measure it that it is
directly linked to the Spearman correlation, $\rho$, that can be
computed as $$\begin{array}{r}
{\frac{\theta + 1}{\theta - 1} - \frac{2\theta}{(\theta - 1)^{2}}\log(\theta)}
\end{array}$$ when $\theta \neq 1$, if $\theta = 1$ then $\rho = 0$.

This model has a more free parameter than the Clayton-Oakes model.

``` r
 library(mets)
 data(diabetes)
 
 # Marginal Cox model  with treat as covariate
 margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
 # Clayton-Oakes, MLE 
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
 
 # Plackett model
 mph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
 fitp <- survival.twostage(mph,data=diabetes,theta=3.0,
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
 
 # without covariates but with stratafied 
 marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
 fitpa <- survival.twostage(marg,data=diabetes,theta=1.0,
                 clusters=diabetes$id)
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
```

With a regression design

``` r
 mm <- model.matrix(~-1+factor(adult),diabetes)
 fitp <- survival.twostage(mph,data=diabetes,theta=3.0,Nit=40,
                clusters=diabetes$id,var.link=1,model="plackett",
        theta.des=mm)
 summary(fitp)
#> Dependence parameter for Odds-Ratio (Plackett) model 
#> With log-link
#> $estimates
#>                log-Coef.        SE        z       P-val Spearman Corr.
#> factor(adult)1  1.098332 0.3630524 3.025272 0.002484097      0.3519987
#> factor(adult)2  1.231962 0.5233261 2.354100 0.018567607      0.3909504
#>                       SE
#> factor(adult)1 0.1074107
#> factor(adult)2 0.1501987
#> 
#> $or
#>             Estimate Std.Err    2.5% 97.5% P-value
#> dependence1    2.999   1.089  0.8650 5.133 0.00588
#> dependence2    3.428   1.794 -0.0881 6.944 0.05602
#> 
#> $type
#> [1] "plackett"
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

``` r
 # Piecewise constant cross hazards ratio modelling

 d <- subset(simClaytonOakes(1000,2,0.5,0,stoptime=2,left=0),!truncated)
 udp <- piecewise.twostage(c(0,0.5,2),data=d,id="cluster",timevar="time",
                           status="status",model="plackett",silent=0)
#> Data-set  1 out of  4
#>   Number of joint events: 254 of  1000
#> Data-set  2 out of  4
#>   Number of joint events: 106 of  609
#> Data-set  3 out of  4
#>   Number of joint events: 136 of  640
#> Data-set  4 out of  4
#>   Number of joint events: 328 of  503
 summary(udp)
#> [1] 1
#> Dependence parameter for Plackett model 
#> log-coefficient for dependence parameter (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  1.794 (0.106)  1.888 (0.139)
#> 0.5 - 2  1.582 (0.141)  2.18  (0.109)
#> 
#> Spearman Correlation (SE) 
#>          0 - 0.5        0.5 - 2      
#> 0 - 0.5  0.54  (0.026)  0.563 (0.033)
#> 0.5 - 2  0.487 (0.037)  0.628 (0.023)
```

## SessionInfo

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] mets_1.3.9
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.5              knitr_1.51             rlang_1.1.7           
#>  [4] xfun_0.55              textshaping_1.0.4      jsonlite_2.0.0        
#>  [7] listenv_0.10.0         future.apply_1.20.1    lava_1.8.2            
#> [10] htmltools_0.5.9        ragg_1.5.0             sass_0.4.10           
#> [13] rmarkdown_2.30         grid_4.5.2             evaluate_1.0.5        
#> [16] jquerylib_0.1.4        fastmap_1.2.0          numDeriv_2016.8-1.1   
#> [19] yaml_2.3.12            mvtnorm_1.3-3          lifecycle_1.0.5       
#> [22] timereg_2.0.7          compiler_4.5.2         codetools_0.2-20      
#> [25] fs_1.6.6               htmlwidgets_1.6.4      Rcpp_1.1.1            
#> [28] future_1.68.0          lattice_0.22-7         systemfonts_1.3.1     
#> [31] digest_0.6.39          R6_2.6.1               parallelly_1.46.1     
#> [34] parallel_4.5.2         splines_4.5.2          Matrix_1.7-4          
#> [37] bslib_0.9.0            tools_4.5.2            RcppArmadillo_15.2.3-1
#> [40] globals_0.18.0         survival_3.8-3         pkgdown_2.2.0         
#> [43] cachem_1.1.0           desc_1.4.3
```
