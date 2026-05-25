# Marginal modelling of clustered survival data

## Overview

A basic component of our modelling of multivariate survival data is that
many models are built around marginals on Cox form. The marginal Cox
model can be fitted efficiently in the mets package, in particular the
handling of strata and robust standard errors is optimized.

The basic models assumes that each subject has a marginal on Cox-form
\lambda\_{g(k,i)}(t) \exp( X\_{ki}^T \beta). where g(k,i) gives the
strata for the subject.

We here discuss and show how to get robust standard errors of

- the regression parameters

- the baseline

and how to do goodness of fit test using

- cumulative residuals score test

First we generate some data from the Clayton-Oakes model, with 5 members
in each cluster and a variance parameter at 2

``` r

 library(mets)
 options(warn=-1)
 set.seed(1000) # to control output in simulations for p-values below.
 n <- 1000
 k <- 5
 theta <- 2
 data <- sim_ClaytonOakes(n,k,theta,0.3,3)
```

The data has one subject per row.

- time : time of event
- status : 1 for event and 0 for censoring
- x : x is a binary covariate
- cluster : cluster

Now we fit the model and produce robust standard errors for both
regression parameters and baseline.

First, recall that the baseline for strata g is asymptotically
equivalent to \begin{align} \hat A_g(t) - A_g(t) & = \sum\_{k \in g}
\left( \sum\_{i: ki \in g} \int_0^t \frac{1}{S\_{0,g}} dM\_{ki}^g -
P^g(t) \beta_k \right) \end{align} with P^g(t) = \int_0^t E_g(s) d \hat
\Lambda_g(s) the derivative of \int_0^t 1/S\_{0,g}(s) dN\_{\cdot g} wrt
to \beta, and \begin{align} \hat \beta - \beta & = I(\tau)^{-1} \sum_k (
\sum_i \int_0^\tau (Z\_{ki} - E\_{g}) dM\_{ki}^g ) = \sum_k \beta\_{k}
\end{align} with \begin{align} M\_{ki}^g(t) & = N\_{ki}(t) - \int_0^t
Y\_{ki}(s) \exp( Z\_{ki} \beta) d \Lambda\_{g(k,i)}(t), \\ \beta\_{k} &
= I(\tau)^{-1} \sum_i \int_0^\tau (Z\_{ki} - E\_{g}) dM\_{ki}^g
\end{align} the basic 0-mean processes, that are martingales in the iid
setting, and I(t) is the derivative of the total score, \hat
U(t,\beta)), with respect to \beta evaluated at time t.

The variance of the baseline of strata g is estimated by \begin{align}
\sum\_{k \in g} ( \sum\_{i: ki \in g} \int_0^t \frac{1}{S\_{0,g(k,i)}}
d\hat M\_{ki}^g - P^g(t) \beta_k )^2 \end{align} that can be computed
using the particular structure of \begin{align} d \hat M\_{ik}^g(t) & =
dN\_{ik}(t) - \frac{1}{S\_{0,g(i,k)}} \exp(Z\_{ik} \beta) dN\_{g.}(t)
\end{align}

This robust variance of the baseline and the iid decomposition for \beta
is computed in mets as:

``` r

   out <- phreg(Surv(time,status)~x+cluster(cluster),data=data)
   summary(out)
#> 
#>     n events
#>  5000   4854
#> coefficients:
#>   Estimate     S.E.  dU^-1/2 P-value
#> x 0.287859 0.028177 0.028897       0
#> 
#> exp(coefficients):
#>   Estimate   2.5%  97.5%
#> x   1.3336 1.2619 1.4093
   # robust standard errors attached to output
   rob <- robust_phreg(out)
```

We can get the iid decomposition of the \hat \beta - \beta by

``` r

   # making iid decomposition of regression parameters
   betaiid <- IC(out)
   head(betaiid)
#>             x
#> 1 -0.34616008
#> 2 -1.44918926
#> 3 -0.03898156
#> 4  0.42156050
#> 5  0.34253904
#> 6 -0.07706668
   # robust standard errors
   crossprod(betaiid/NROW(betaiid))^.5
#>            x
#> x 0.02817714
   # same as 
```

We now look at the plot with robust standard errors

``` r

  plot(rob,se=TRUE,robust=TRUE,col=3)
```

![Cumulative baseline hazard with robust standard
errors.](figure/marginal-cox-unnamed-chunk-5-1.png)

We can also make survival prediction with robust standard errors using
the phreg.

``` r

  pp <-  predict(out,data[1:20,],se=TRUE,robust=TRUE)
  plot(pp,se=TRUE,whichx=1:10)
```

![Survival predictions with robust standard
errors.](figure/marginal-cox-unnamed-chunk-6-1.png)

Finally, just to check that we can recover the model we also estimate
the dependence parameter

``` r

tt <- twostageMLE(out,data=data)
summary(tt)
#> Dependence parameter for Clayton-Oakes model
#> Variance of Gamma distributed random effects
#> $estimates
#>                 Coef.         SE        z P-val Kendall tau        SE
#> dependence1 0.5316753 0.03497789 15.20032     0   0.2100093 0.0109146
#> 
#> $type
#> NULL
#> 
#> attr(,"class")
#> [1] "summary.mets.twostage"
```

## Goodness of fit

The observed score process is given by \begin{align} U(t,\hat \beta) & =
\sum_k \sum_i \int_0^t (Z\_{ki} - \hat E_g ) d \hat M\_{ki}^g
\end{align} where g is strata g(k,i). The observed score has the iid
decomposition \begin{align} \hat U(t) = \sum_k \sum_i \int_0^t
(Z\_{ki} - E_g) dM\_{ki}^g - I(t) \sum_k \beta_k \end{align} where
\beta_k is the iid decomposition of the score process for the true \beta
\begin{align} \beta_k & = I(\tau)^{-1} \sum_i \int_0^\tau (Z\_{ki} - E_g
) d M\_{ki}^g \end{align} and I(t) is the derivative of the total score,
\hat U(t,\beta)), with respect to \beta evaluated at time t.

This observed score can be resampled given it is on iid form in terms of
clusters.

Now using the cumulative score process for checking proportional hazards

``` r

gout <- gof(out)
gout
#> Cumulative score process test for Proportionality:
#>   Sup|U(t)|  pval
#> x  30.24353 0.401
```

The p-value reflects whether the observed score process is consistent
with the model.

``` r

  plot(gout)
```

![Cumulative score process for goodness of
fit.](figure/marginal-cox-unnamed-chunk-9-1.png)

### Computational aspects

The score processes can be resampled as in Lin, Wei, Ying (1993) using
the martingale structure, such that the observed score process is
resampled by \begin{align} \sum_k \sum_i \int_0^t g\_{ki} (Z\_{ki} -
E_g) dN\_{ki} - I(t) I^{-1}(\tau) g\_{ki} \int_0^{\tau} (Z\_{ki} - E_g)
dN\_{ki} . \end{align} where g\_{ki} are i.i.d. standard normals.

Based on the zero mean processes we more generally with clusters can
resample the score process. For resampling of score process we need
\begin{align} U(t,\beta) & = \sum_k \sum_i g_k \int_0^t (Z\_{ki} - E_g )
dM\_{ki}^g \end{align} where g is strata. We write g_k as g\_{ki} and
thus repeating g_k within each cluster.

Computations are done using that \begin{align\*} \int_0^t (Z\_{ki} -
E\_{g}) dM\_{ki}^g & = \int_0^t (Z\_{ki} - E\_{g}) dN\_{ki}^g -
\int_0^{t} (Z\_{ki} - E\_{g}) Y\_{ki}(u) d\Lambda^g(u) \end{align\*}
therefore and summing the compensator part with the g\_{ki} multipliers
then gives for each strata g \begin{align\*} & \int_0^t
\frac{S\_{1g}^w(u)}{S\_{0g}(u)} dN\_{g.}(v) - \int_0^t E\_{g}(u)
\frac{S\_{0g}^w(u)}{S\_{0g}(u)} dN\_{g.}(v) \end{align\*} with  
\begin{align\*} S\_{jg}^w(t) & = \sum\_{ki \in g} \exp(Z\_{ki} \beta)
Z\_{ki}^j Y\_{ki}(t) g\_{ki} \\ S\_{jg}(t) & = \sum\_{ki \in g}
\exp(Z\_{ki} \beta) Z\_{ki}^j Y\_{ki}(t). \end{align\*}

## Cluster stratified Cox models

For clustered data it is possible to estimate the regression coefficient
within clusters by using Cox’s partial likelihood stratified on
clusters.

Note, here that the data is generated with a different subject specific
structure, so we will not recover \beta = 0.3 and the model will not be
a proportional Cox model; we would therefore also expect to reject
“proportionality” with the goodness-of-fit test.

The model can be thought of as \lambda\_{g(k,i)} (t) \exp( X\_{ki}^T
\beta) where \lambda_g(t) is some cluster specific baseline.

The regression coefficient \beta can be estimated by using the partial
likelihood for clusters.

``` r

 out <- phreg(Surv(time,status)~x+strata(cluster),data=data)
 summary(out)
#> 
#>     n events
#>  5000   4854
#> coefficients:
#>   Estimate     S.E.  dU^-1/2 P-value
#> x 0.406307 0.032925 0.039226       0
#> 
#> exp(coefficients):
#>   Estimate   2.5%  97.5%
#> x   1.5013 1.4074 1.6013
```

## SessionInfo

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /home/kkzh/.asdf/installs/r/4.6.0/lib/R/lib/libRblas.so 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Copenhagen
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] timereg_2.0.7  survival_3.8-6 mets_1.3.10   
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6              knitr_1.51             rlang_1.2.0           
#>  [4] xfun_0.57              otel_0.2.0             future.apply_1.20.2   
#>  [7] listenv_0.10.1         lava_1.9.1             stats4_4.6.0          
#> [10] grid_4.6.0             evaluate_1.0.5         yaml_2.3.12           
#> [13] mvtnorm_1.3-7          numDeriv_2016.8-1.1    compiler_4.6.0        
#> [16] codetools_0.2-20       Rcpp_1.1.1-1.1         ucminf_1.2.3          
#> [19] future_1.70.0          lattice_0.22-9         digest_0.6.39         
#> [22] parallelly_1.47.0      parallel_4.6.0         splines_4.6.0         
#> [25] Matrix_1.7-5           tools_4.6.0            RcppArmadillo_15.2.6-1
#> [28] globals_0.19.1
```
