# Marginal modelling of clustered survival data

## Overview

A basic component for our modelling of multivariate survival data is
that many models are build around marginals that on Cox form. The
marginal Cox model can be fitted efficiently in the mets package, in
particular the handling of strata and robust standard errors is
optimized.

The basic models assumes that each subject has a marginal on Cox-form
$$\lambda_{g{(k,i)}}(t)\exp\left( X_{ki}^{T}\beta \right).$$ where
$g(k,i)$ gives the strata for the subject.

We here discuss and show how to get robust standard errors of

- the regression parameters

- the baseline

and how to do goodness of fit test using

- cumulative residuals score test

First we generate some data from the Clayton-Oakes model, with $5$
members in each cluster and a variance parameter at $2$

``` r
 library(mets)
 options(warn=-1)
 set.seed(1000) # to control output in simulatins for p-values below.
 n <- 1000
 k <- 5
 theta <- 2
 data <- simClaytonOakes(n,k,theta,0.3,3)
```

The data is on has one subject per row.

- time : time of event
- status : 1 for event and 0 for censoring
- x : x is a binary covariate
- cluster : cluster

Now we fit the model and produce robust standard errors for both
regression parameters and baseline.

First, recall that the baseline for strata $g$ is asymptotically
equivalent to $$\begin{aligned}
{{\widehat{A}}_{g}(t) - A_{g}(t)} & {= \sum\limits_{k \in g}\left( \sum\limits_{i:ki \in g}\int_{0}^{t}\frac{1}{S_{0,g}}dM_{ki}^{g} - P^{g}(t)\beta_{k} \right)}
\end{aligned}$$ with
$P^{g}(t) = \int_{0}^{t}E_{g}(s)d{\widehat{\Lambda}}_{g}(s)$ the
derivative of $\int_{0}^{t}1/S_{0,g}(s)dN_{\cdot g}$ wrt to $\beta$, and
$$\begin{aligned}
{\widehat{\beta} - \beta} & {= I(\tau)^{- 1}\sum\limits_{k}\left( \sum\limits_{i}\int_{0}^{\tau}\left( Z_{ki} - E_{g} \right)dM_{ki}^{g} \right) = \sum\limits_{k}\beta_{k}}
\end{aligned}$$ with $$\begin{aligned}
{M_{ki}^{g}(t)} & {= N_{ki}(t) - \int_{0}^{t}Y_{ki}(s)\exp\left( Z_{ki}\beta \right)d\Lambda_{g{(k,i)}}(t),} \\
\beta_{k} & {= I(\tau)^{- 1}\sum\limits_{i}\int_{0}^{\tau}\left( Z_{ki} - E_{g} \right)dM_{ki}^{g}}
\end{aligned}$$ the basic 0-mean processes, that are martingales in the
iid setting, and $I(t)$ is the derivative of the total score,
$\widehat{U}(t,\beta))$, with respect to $\beta$ evaluated at time $t$.

The variance of the baseline of strata g is estimated by
$$\begin{array}{r}
{\sum\limits_{k \in g}\left( \sum\limits_{i:ki \in g}\int_{0}^{t}\frac{1}{S_{0,g{(k,i)}}}d{\widehat{M}}_{ki}^{g} - P^{g}(t)\beta_{k} \right)^{2}}
\end{array}$$ that can be computed using the particular structure of
$$\begin{aligned}
{d{\widehat{M}}_{ik}^{g}(t)} & {= dN_{ik}(t) - \frac{1}{S_{0,g{(i,k)}}}\exp\left( Z_{ik}\beta \right)dN_{g.}(t)}
\end{aligned}$$

This robust variance of the baseline and the iid decomposition for
$\beta$ is computed in mets as:

``` r
   out <- phreg(Surv(time,status)~x+cluster(cluster),data=data)
   summary(out)
#> 
#>     n events
#>  5000   4854
#> coeffients:
#>   Estimate     S.E.  dU^-1/2 P-value
#> x 0.287859 0.028177 0.028897       0
#> 
#> exp(coeffients):
#>   Estimate   2.5%  97.5%
#> x   1.3336 1.2619 1.4093
   # robust standard errors attached to output
   rob <- robust.phreg(out)
```

We can get the iid decomposition of the $\widehat{\beta} - \beta$ by

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

![](marginal-cox_files/figure-html/unnamed-chunk-5-1.png)

We can also make survival prediction with robust standard errors using
the phreg.

``` r
  pp <-  predict(out,data[1:20,],se=TRUE,robust=TRUE)
  plot(pp,se=TRUE,whichx=1:10)
```

![](marginal-cox_files/figure-html/unnamed-chunk-6-1.png)

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

The observed score process is given by $$\begin{aligned}
{U\left( t,\widehat{\beta} \right)} & {= \sum\limits_{k}\sum\limits_{i}\int_{0}^{t}\left( Z_{ki} - {\widehat{E}}_{g} \right)d{\widehat{M}}_{ki}^{g}}
\end{aligned}$$ where $g$ is strata $g(k,i)$. The observed score has the
iid decomposition $$\begin{array}{r}
{\widehat{U}(t) = \sum\limits_{k}\sum\limits_{i}\int_{0}^{t}\left( Z_{ki} - E_{g} \right)dM_{ki}^{g} - I(t)\sum\limits_{k}\beta_{k}}
\end{array}$$ where $\beta_{k}$ is the iid decomposition of the score
process for the true $\beta$$$\begin{aligned}
\beta_{k} & {= I(\tau)^{- 1}\sum\limits_{i}\int_{0}^{\tau}\left( Z_{ki} - E_{g} \right)dM_{ki}^{g}}
\end{aligned}$$ and $I(t)$ is the derivative of the total score,
$\widehat{U}(t,\beta))$, with respect to $\beta$ evaluated at time $t$.

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

The p-value reflects wheter the observed score process is consistent
with the model.

``` r
  plot(gout)
```

![](marginal-cox_files/figure-html/unnamed-chunk-9-1.png)

### Computational aspects

The score processes can be resampled as in Lin, Wei, Ying (1993) using
the martingale structure, such that the observed score process is
resampled by $$\begin{array}{r}
{\sum\limits_{k}\sum\limits_{i}\int_{0}^{t}g_{ki}\left( Z_{ki} - E_{g} \right)dN_{ki} - I(t)I^{- 1}(\tau)g_{ki}\int_{0}^{\tau}\left( Z_{ki} - E_{g} \right)dN_{ki}.}
\end{array}$$ where $g_{ki}$ are i.i.d. standard normals.

Based on the zero mean processes we more generally with clusters can
resample the score process. For resampling of score process we need
$$\begin{aligned}
{U(t,\beta)} & {= \sum\limits_{k}\sum\limits_{i}g_{k}\int_{0}^{t}\left( Z_{ki} - E_{g} \right)dM_{ki}^{g}}
\end{aligned}$$ where $g$ is strata. We write $g_{k}$ as $g_{ki}$ and
thus repeating $g_{k}$ within each cluster.

Computations are done using that $$\begin{aligned}
{\int_{0}^{t}\left( Z_{ki} - E_{g} \right)dM_{ki}^{g}} & {= \int_{0}^{t}\left( Z_{ki} - E_{g} \right)dN_{ki}^{g} - \int_{0}^{t}\left( Z_{ki} - E_{g} \right)Y_{ki}(u)d\Lambda^{g}(u)}
\end{aligned}$$ therefore and summing the compensator part with the
$g_{ki}$ multipliers then gives for each strata $g$$$\begin{aligned}
 & {\int_{0}^{t}\frac{S_{1g}^{w}(u)}{S_{0g}(u)}dN_{g.}(v) - \int_{0}^{t}E_{g}(u)\frac{S_{0g}^{w}(u)}{S_{0g}(u)}dN_{g.}(v)}
\end{aligned}$$ with  
$$\begin{aligned}
{S_{jg}^{w}(t)} & {= \sum\limits_{ki \in g}\exp\left( Z_{ki}\beta \right)Z_{ki}^{j}Y_{ki}(t)g_{ki}} \\
{S_{jg}(t)} & {= \sum\limits_{ki \in g}\exp\left( Z_{ki}\beta \right)Z_{ki}^{j}Y_{ki}(t).}
\end{aligned}$$

## Cluster stratified Cox models

For clustered data it is possible to estimate the regression coefficient
within clusters by using Cox’s partial likelihood stratified on
clusters.

Note, here that the data is generated with a different subject specific
structure, so we will not recover the $\beta$ at 0.3 and the model will
not be a proportional Cox model, we we would also expect to reject
“proportionality” with the gof-test.

The model can be thought of as
$$\lambda_{g{(k,i)}}(t)\exp\left( X_{ki}^{T}\beta \right)$$ where
$\lambda_{g}(t)$ is some cluster specific baseline.

The regression coefficient $\beta$ can be estimated by using the partial
likelihood for clusters.

``` r
 out <- phreg(Surv(time,status)~x+strata(cluster),data=data)
 summary(out)
#> 
#>     n events
#>  5000   4854
#> coeffients:
#>   Estimate     S.E.  dU^-1/2 P-value
#> x 0.406307 0.032925 0.039226       0
#> 
#> exp(coeffients):
#>   Estimate   2.5%  97.5%
#> x   1.5013 1.4074 1.6013
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
