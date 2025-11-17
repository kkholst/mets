# Binomial Regression for Survival and Competing Risks Data

## Binomial Regression for censored data

The binreg function can fit a logistic link model with IPCW adjustment
for a specific time-point, and can thus be used for describing survival
or competing risks data. The function can be used for large data and is
completely scalable, that is, linear in the data. A nice feature is that
influcence functions are computed and are available, and can thus be
used for all other settings based on these parameters.

In addition and to summarize

- the censoring weights can be strata dependent
- predictions can be computed with standard errors
- computation time linear in data
  - including standard errors
- clusters can be given and then cluster corrected standard errors are
  computed

## Details

The binreg function does direct binomial regression for one time-point,
$t$, fitting the model $$\begin{aligned}
{P\left( T \leq t,\epsilon = 1|X \right)} & {= \text{expit}\left( X^{T}\beta \right) = F_{1}(t,X,\beta)}
\end{aligned}$$ to an IPCW adjusted estmating equation (EE) with
response $Y(t) = I(T \leq t,\epsilon = 1)$$$\begin{aligned}
{U\left( \beta,{\widehat{G}}_{c} \right) =} & {X\left( Y(t)\frac{\Delta(t)}{{\widehat{G}}_{c}\left( T_{i} \land t \right)} - \text{expit}\left( X^{T}\beta \right) \right) = 0,}
\end{aligned}$$ with $G_{c}(t) = P(C > t)$, the censoring survival
distribution, and with
$\Delta(t) = I\left( C_{i} > T_{i} \land t \right)$ the indicator of
being uncensored at time $t$ (type=“I”). The default type=“II” is to
augment with a censoring term, that is solve $$\begin{aligned}
 & {U\left( \beta,{\widehat{G}}_{c} \right) + \int_{0}^{t}X\frac{\widehat{E}\left( Y(t)|T > u \right)}{{\widehat{G}}_{c}(u)}d{\widehat{M}}_{c}(u) = 0}
\end{aligned}$$ where $M_{c}(u)$ is the censoring martingale, this
typically improves the performance. This is equivlent to the
pseudo-value approach (see Overgaard (2025)).

The influence function for the type=“II” estimator is $$\begin{array}{r}
{U\left( \beta,G_{c} \right) + \int_{0}^{t}X\frac{E\left( Y|T > u \right)}{G_{c}(u)}dM_{c}(u) - \int_{0}^{t}\frac{E\left( X|T > u \right)E\left( Y|T > u \right)}{G_{c}(u)}dM_{c}(u) - \int_{0}^{t}\frac{E\left( XY|T > u \right)}{G_{c}(u)}dM_{c}(u)}
\end{array}$$ and for type=“I” $$\begin{aligned}
 & {U(\beta) + \int_{0}^{t}\frac{E\left( XY|T > u \right)}{G_{c}(u)}dM_{c}(u).}
\end{aligned}$$ The means $E\left( XY(t)|T > u \right)$ and
$E\left( Y(t)|T > u \right)$ are estimated by IPCW estimators among
survivors to get estimates of the influence functions.

The function logitIPCW instead considers $$\begin{aligned}
{U^{glm}\left( \beta,{\widehat{G}}_{c} \right) =} & {\frac{\Delta(t)}{{\widehat{G}}_{c}\left( T_{i} \land t \right)}X\left( Y(t) - \text{expit}\left( X^{T}\beta \right) \right) = 0.}
\end{aligned}$$ This score equation is quite similar to those of the
binreg, and exactly the same when the censoring model is
fully-nonparametric.

The logitIPCW has influence function $$\begin{aligned}
 & {U^{glm}\left( \beta,G_{c} \right) + \int_{0}^{t}\frac{E\left( X\left( Y - F_{1}(t,\beta) \right)|T > u \right)}{G_{c}(u)}dM_{c}(u)}
\end{aligned}$$

Which estimator performs the best depends on the censoring distribution
and it seems that the binreg with type=“II” performs overall quite
nicely (see Blanche et al (2023) and Overgaard (2024)). For the full
estimated censoring model all estimators have the same influence
function (see Blanche et al (2023)).

Additional functions logitATE, and binregATE computes the average
treatment effect. We demonstrate this in another vignette.

The functions logitATE/binregATE can be used there is no censoring and
we thus have simple binary outcome.

The variance is based on sandwich formula with IPCW adjustment (using
the influence functions), and naive.var is the variance under known
censoring model. The influence functions are stored in the output.
Clusters can be specified to get cluster corrected standard errors.

## Examples

``` r
 library(mets)
 options(warn=-1)
 set.seed(1000) # to control output in random noise just below.
 data(bmt)
 bmt$time <- bmt$time+runif(nrow(bmt))*0.01

 # logistic regresion with IPCW binomial regression 
 out <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50)
 summary(out)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.180338  0.126748 -0.428760  0.068084  0.1548
#> tcell       -0.418545  0.345480 -1.095675  0.258584  0.2257
#> platelet    -0.437644  0.240978 -0.909952  0.034665  0.0694
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.83499 0.65132 1.0705
#> tcell        0.65800 0.33431 1.2951
#> platelet     0.64556 0.40254 1.0353
```

We can also compute predictions using the estimates

``` r
 predict(out,data.frame(tcell=c(0,1),platelet=c(1,1)),se=TRUE)
#>        pred         se     lower     upper
#> 1 0.3502406 0.04847582 0.2552280 0.4452533
#> 2 0.2618207 0.06969334 0.1252217 0.3984196
```

Further the censoring model can depend on strata

``` r
 outs <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50,cens.model=~strata(tcell,platelet))
 summary(outs)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.180697  0.127414 -0.430424  0.069030  0.1561
#> tcell       -0.365928  0.350632 -1.053154  0.321299  0.2967
#> platelet    -0.433494  0.240270 -0.904415  0.037428  0.0712
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.83469 0.65023 1.0715
#> tcell        0.69355 0.34884 1.3789
#> platelet     0.64824 0.40478 1.0381
```

## Absolute risk differences and ratio

Now for illustrations I wish to consider the absolute risk difference
depending on tcell

``` r
 outs <- binreg(Event(time,cause)~tcell,bmt,time=50,cens.model=~strata(tcell))
 summary(outs)
#>    n events
#>  408    160
#> 
#>  408 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept) -0.30054  0.11153 -0.51914 -0.08194  0.0070
#> tcell       -0.51741  0.33981 -1.18342  0.14860  0.1278
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.74042 0.59503 0.9213
#> tcell        0.59606 0.30623 1.1602
```

the risk difference is

``` r
ps <-  predict(outs,data.frame(tcell=c(0,1)),se=TRUE)
ps
#>        pred         se     lower     upper
#> 1 0.4254253 0.02726306 0.3719897 0.4788609
#> 2 0.3061988 0.06819019 0.1725461 0.4398516
sum( c(1,-1) * ps[,1])
#> [1] 0.1192264
```

Getting the standard errors are easy enough since the two-groups are
independent. In the case where we in addition had adjusted for other
covariates, however, we would need the to apply the delta-theorem thus
using the relevant covariances along the lines of

``` r
dd <- data.frame(tcell=c(0,1))
p <- predict(outs,dd)

riskdifratio <- function(p,contrast=c(1,-1)) {
   outs$coef <- p
   p <- predict(outs,dd)[,1]
   pd <- sum(contrast*p)
   r1 <- p[1]/p[2]
   r2 <- p[2]/p[1]
   return(c(pd,r1,r2))
}
     
estimate(outs,f=riskdifratio,dd,null=c(0,1,1))
#>      Estimate Std.Err     2.5%  97.5% P-value
#> [p1]   0.1192 0.07344 -0.02471 0.2632 0.10448
#> [p2]   1.3894 0.32197  0.75833 2.0204 0.22652
#> [p3]   0.7197 0.16679  0.39284 1.0467 0.09291
#> 
#>  Null Hypothesis: 
#>   [p1] = 0
#>   [p2] = 1
#>   [p3] = 1 
#>  
#> chisq = 12.0249, df = 3, p-value = 0.007298
```

same as

``` r
run <- 0
if (run==1) {
library(prodlim)
pl <- prodlim(Hist(time,cause)~tcell,bmt)
spl <- summary(pl,times=50,asMatrix=TRUE)
spl
}
```

## Augmenting the Binomial Regression

Rather than using a larger censoring model we can also compute an
augmentation term and then fit the binomial regression model based on
this augmentation term. Here we compute the augmentation based on
stratified non-parametric estimates of $F_{1}\left( t,S(X) \right)$,
where $S(X)$ gives strata based on $X$ as a working model.

Computes the augmentation term for each individual as well as the sum
$$\begin{aligned}
A & {= \int_{0}^{t}H(u,X)\frac{1}{S^{*}(u,s)}\frac{1}{G_{c}(u)}dM_{c}(u)}
\end{aligned}$$ with $$\begin{aligned}
{H(u,X)} & {= F_{1}^{*}\left( t,S(X) \right) - F_{1}^{*}\left( u,S(X) \right)}
\end{aligned}$$ using a KM for $G_{c}(t)$ and a working model for
cumulative baseline related to $F_{1}^{*}(t,s)$ and $s$ is strata,
$S^{*}(t,s) = 1 - F_{1}^{*}(t,s) - F_{2}^{*}(t,s)$.

Standard errors computed under assumption of correct but estimated
$G_{c}(s)$ model.

``` r
 data(bmt)
 dcut(bmt,breaks=2) <- ~age 
 out1<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
              strata(platelet,agecat.2),data=bmt,cause=1,time=40)
 summary(out1)
#>    n events
#>  408    157
#> 
#>  408 clusters
#> coeffients:
#>                      Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)          -0.51295  0.17090 -0.84791 -0.17799  0.0027
#> platelet             -0.63011  0.23585 -1.09237 -0.16785  0.0075
#> agecat.2(0.203,1.94]  0.55926  0.21211  0.14353  0.97500  0.0084
#> 
#> exp(coeffients):
#>                      Estimate    2.5%  97.5%
#> (Intercept)           0.59873 0.42831 0.8370
#> platelet              0.53253 0.33542 0.8455
#> agecat.2(0.203,1.94]  1.74938 1.15434 2.6512

 out2<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
     strata(platelet,agecat.2)+strataC(platelet),data=bmt,cause=1,time=40)
 summary(out2)
#>    n events
#>  408    157
#> 
#>  408 clusters
#> coeffients:
#>                      Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)          -0.51346  0.17109 -0.84879 -0.17814  0.0027
#> platelet             -0.63636  0.23653 -1.09996 -0.17276  0.0071
#> agecat.2(0.203,1.94]  0.56280  0.21229  0.14672  0.97889  0.0080
#> 
#> exp(coeffients):
#>                      Estimate    2.5%  97.5%
#> (Intercept)           0.59842 0.42793 0.8368
#> platelet              0.52922 0.33288 0.8413
#> agecat.2(0.203,1.94]  1.75559 1.15803 2.6615
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
#>  [1] cli_3.6.5           knitr_1.50          rlang_1.1.6        
#>  [4] xfun_0.54           textshaping_1.0.4   jsonlite_2.0.0     
#>  [7] listenv_0.10.0      future.apply_1.20.0 lava_1.8.2         
#> [10] htmltools_0.5.8.1   ragg_1.5.0          sass_0.4.10        
#> [13] rmarkdown_2.30      grid_4.5.2          evaluate_1.0.5     
#> [16] jquerylib_0.1.4     fastmap_1.2.0       numDeriv_2016.8-1.1
#> [19] yaml_2.3.10         mvtnorm_1.3-3       lifecycle_1.0.4    
#> [22] timereg_2.0.7       compiler_4.5.2      codetools_0.2-20   
#> [25] fs_1.6.6            htmlwidgets_1.6.4   Rcpp_1.1.0         
#> [28] future_1.68.0       lattice_0.22-7      systemfonts_1.3.1  
#> [31] digest_0.6.38       R6_2.6.1            parallelly_1.45.1  
#> [34] parallel_4.5.2      splines_4.5.2       Matrix_1.7-4       
#> [37] bslib_0.9.0         tools_4.5.2         globals_0.18.0     
#> [40] survival_3.8-3      pkgdown_2.2.0       cachem_1.1.0       
#> [43] desc_1.4.3
```
