# Discrete Interval Censored Survival Models

## Discrete Inteval Censored survival times

We consider the cumulative odds model for the probability of dying
before time t: $$\begin{aligned}
{\text{logit}\left( P\left( T \leq t|x \right) \right)} & {= \log\left( G(t) \right) + x^{T}\beta} \\
{P\left( T \leq t|x \right)} & {= \frac{G(t)exp\left( x^{T}\beta \right)}{1 + G(t)exp\left( x^{T}\beta \right)}} \\
{P\left( T > t|x \right)} & {= \frac{1}{1 + G(t)exp\left( x^{T}\beta \right)}}
\end{aligned}$$

Input are intervals given by $\rbrack t_{l},t_{r}\rbrack$ where t_r can
be infinity for right-censored intervals. When the data is discrete, in
contrast to grouping of continuous data, $\rbrack 0,1\rbrack$ then the
intervals $\rbrack j,j + 1\rbrack$ will be equvilant to an observation
at j+1 (see below example).

Likelihood is maximized: $$\begin{array}{r}
{\prod\limits_{i}P\left( T_{i} > t_{il}|x \right) - P\left( T_{i} > t_{ir}|x \right).}
\end{array}$$

This model is also called the cumulative odds model $$\begin{aligned}
{P\left( T \leq t|x \right)} & {= \frac{G(t)exp\left( x^{T}\beta \right)}{1 + G(t)exp\left( x^{T}\beta \right)}.}
\end{aligned}$$ and $\beta$ says something about the OR of probability
of being before $t$.

The baseline is parametrized as $$\begin{aligned}
{G(t)} & {= \sum\limits_{j \leq t}\exp\left( \alpha_{j} \right)}
\end{aligned}$$

An important consequence of the model is that for all cut-points $t$ we
have the same OR parameters for the OR of being early or later than $t$.

## Discrete TTP

First we look at some time to pregnancy data (simulated discrete
survival data) that is right-censored, and set it up to fit the
cumulative odds model by constructing the intervals appropriately:

``` r
library(mets)

data(ttpd) 
dtable(ttpd,~entry+time2)
#> 
#>       time2   1   2   3   4   5   6 Inf
#> entry                                  
#> 0           316   0   0   0   0   0   0
#> 1             0 133   0   0   0   0   0
#> 2             0   0 150   0   0   0   0
#> 3             0   0   0  23   0   0   0
#> 4             0   0   0   0  90   0   0
#> 5             0   0   0   0   0  68   0
#> 6             0   0   0   0   0   0 220
out <- interval.logitsurv.discrete(Interval(entry,time2)~X1+X2+X3+X4,ttpd)
summary(out)
#> $baseline
#>       Estimate Std.Err   2.5%   97.5%   P-value
#> time1  -2.0064  0.1523 -2.305 -1.7079 1.273e-39
#> time2  -2.1749  0.1599 -2.488 -1.8614 4.118e-42
#> time3  -1.4581  0.1544 -1.761 -1.1554 3.636e-21
#> time4  -2.9260  0.2453 -3.407 -2.4453 8.379e-33
#> time5  -1.2051  0.1706 -1.539 -0.8706 1.633e-12
#> time6  -0.9102  0.1860 -1.275 -0.5457 9.843e-07
#> 
#> $logor
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> X1   0.9913  0.1179 0.76024 1.2223 4.100e-17
#> X2   0.6962  0.1162 0.46847 0.9238 2.064e-09
#> X3   0.3466  0.1159 0.11941 0.5738 2.788e-03
#> X4   0.3223  0.1151 0.09668 0.5478 5.111e-03
#> 
#> $or
#>    Estimate     2.5%    97.5%
#> X1 2.694610 2.138791 3.394874
#> X2 2.006032 1.597554 2.518953
#> X3 1.414239 1.126834 1.774950
#> X4 1.380231 1.101503 1.729490

dfactor(ttpd) <- entry.f~entry
out <- cumoddsreg(entry.f~X1+X2+X3+X4,ttpd)
summary(out)
#> $baseline
#>       Estimate Std.Err   2.5%   97.5%   P-value
#> time1  -2.0064  0.1523 -2.305 -1.7079 1.273e-39
#> time2  -2.1749  0.1599 -2.488 -1.8614 4.118e-42
#> time3  -1.4581  0.1544 -1.761 -1.1554 3.636e-21
#> time4  -2.9260  0.2453 -3.407 -2.4453 8.379e-33
#> time5  -1.2051  0.1706 -1.539 -0.8706 1.633e-12
#> time6  -0.9102  0.1860 -1.275 -0.5457 9.843e-07
#> 
#> $logor
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> X1   0.9913  0.1179 0.76024 1.2223 4.100e-17
#> X2   0.6962  0.1162 0.46847 0.9238 2.064e-09
#> X3   0.3466  0.1159 0.11941 0.5738 2.788e-03
#> X4   0.3223  0.1151 0.09668 0.5478 5.111e-03
#> 
#> $or
#>    Estimate     2.5%    97.5%
#> X1 2.694610 2.138791 3.394874
#> X2 2.006032 1.597554 2.518953
#> X3 1.414239 1.126834 1.774950
#> X4 1.380231 1.101503 1.729490
```

We note that the probability of dying is increased considerably for all
covariates.

Now using this discrete survival model we simulate some data from this
model

``` r
set.seed(1000) # to control output in simulatins for p-values below.
n <- 200
Z <- matrix(rbinom(n*4,1,0.5),n,4)
outsim <- simlogitSurvd(out$coef,Z)
outsim <- transform(outsim,left=time,right=time+1)
outsim <- dtransform(outsim,right=Inf,status==0)
outss <- interval.logitsurv.discrete(Interval(left,right)~+X1+X2+X3+X4,outsim)
summary(outss)
#> $baseline
#>       Estimate Std.Err    2.5%   97.5%   P-value
#> time1  -2.0154  0.3698 -2.7402 -1.2906 5.036e-08
#> time2  -1.5474  0.3473 -2.2281 -0.8666 8.385e-06
#> time3  -0.8119  0.3411 -1.4804 -0.1434 1.729e-02
#> time4  -2.0085  0.5102 -3.0084 -1.0086 8.248e-05
#> time5  -0.2185  0.3858 -0.9746  0.5376 5.711e-01
#> time6   0.2637  0.4618 -0.6415  1.1689 5.681e-01
#> 
#> $logor
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> X1  1.27893  0.2804  0.7293 1.8286 5.106e-06
#> X2  0.39293  0.2635 -0.1235 0.9094 1.359e-01
#> X3 -0.09008  0.2524 -0.5847 0.4045 7.211e-01
#> X4  0.20766  0.2627 -0.3072 0.7225 4.292e-01
#> 
#> $or
#>    Estimate      2.5%    97.5%
#> X1 3.592796 2.0735647 6.225116
#> X2 1.481310 0.8838237 2.482711
#> X3 0.913858 0.5572845 1.498582
#> X4 1.230798 0.7355301 2.059553

pred <- predictlogitSurvd(out,se=TRUE)
plotSurvd(pred,se=TRUE)
```

![](interval-discrete-survival_files/figure-html/unnamed-chunk-3-1.png)

Finally, we look at some data and compare with the icenReg package that
can also fit the proportional odds model for continous or discrete data.
We make the data fully interval censored/discrete by letting also exact
obsevations be only observed to be in an interval.

We consider the interval censored survival times for time from onset of
diabetes to to diabetic nephronpathy, then modify it to observe only
that the event times are in certain intervals.

``` r
test <- 0 
if (test==1) {

require(icenReg)
data(IR_diabetes)
IRdia <- IR_diabetes
## removing fully observed data in continuous version, here making it a discrete observation 
IRdia <- dtransform(IRdia,left=left-1,left==right)
dtable(IRdia,~left+right,level=1)

ints <- with(IRdia,dInterval(left,right,cuts=c(0,5,10,20,30,40,Inf),show=TRUE) )
}
```

We note that the gender effect is equivalent for the two approaches.

``` r
if (test==1) {
ints$Ileft <- ints$left
ints$Iright <- ints$right
IRdia <- cbind(IRdia,data.frame(Ileft=ints$Ileft,Iright=ints$Iright))
dtable(IRdia,~Ileft+Iright)
# 
#       Iright   1   2   3   4   5 Inf
# Ileft                               
# 0             10   1  34  25   4   0
# 1              0  55  19  17   1   1
# 2              0   0 393  16   4   0
# 3              0   0   0 127   1   0
# 4              0   0   0   0  21   0
# 5              0   0   0   0   0   2

outss <- interval.logitsurv.discrete(Interval(Ileft,Iright)~+gender,IRdia)
#            Estimate Std.Err    2.5%    97.5%   P-value
# time1        -3.934  0.3316 -4.5842 -3.28418 1.846e-32
# time2        -2.042  0.1693 -2.3742 -1.71038 1.710e-33
# time3         1.443  0.1481  1.1530  1.73340 1.911e-22
# time4         3.545  0.2629  3.0295  4.06008 1.976e-41
# time5         6.067  0.7757  4.5470  7.58784 5.217e-15
# gendermale   -0.385  0.1691 -0.7165 -0.05351 2.283e-02
summary(outss)
outss$ploglik
# [1] -646.1946

fit <- ic_sp(cbind(Ileft, Iright) ~ gender, data = IRdia, model = "po")
# 
# Model:  Proportional Odds
# Dependency structure assumed: Independence
# Baseline:  semi-parametric 
# Call: ic_sp(formula = cbind(Ileft, Iright) ~ gender, data = IRdia, 
#     model = "po")
# 
#            Estimate Exp(Est)
# gendermale    0.385     1.47
# 
# final llk =  -646.1946 
# Iterations =  6 
# Bootstrap Samples =  0 
# WARNING: only  0  bootstrap samples used for standard errors. 
# Suggest using more bootstrap samples for inference
summary(fit)

## sometimes NR-algorithm needs modifications of stepsize to run 
## outss <- interval.logitsurv.discrete(Interval(Ileft,Iright)~+gender,IRdia,control=list(trace=TRUE,stepsize=1.0))
}
```

Also agrees with the cumulative link regression of the ordinal package,
although the baseline is parametrized differently. In additon the clm is
describing the probability of surviving rather than the probabibility of
dying.

``` r

data(ttpd) 
dtable(ttpd,~entry+time2)
#> 
#>       time2   1   2   3   4   5   6 Inf
#> entry                                  
#> 0           316   0   0   0   0   0   0
#> 1             0 133   0   0   0   0   0
#> 2             0   0 150   0   0   0   0
#> 3             0   0   0  23   0   0   0
#> 4             0   0   0   0  90   0   0
#> 5             0   0   0   0   0  68   0
#> 6             0   0   0   0   0   0 220
ttpd <- dfactor(ttpd,fentry~entry)
out <- cumoddsreg(fentry~X1+X2+X3+X4,ttpd)
summary(out)
#> $baseline
#>       Estimate Std.Err   2.5%   97.5%   P-value
#> time1  -2.0064  0.1523 -2.305 -1.7079 1.273e-39
#> time2  -2.1749  0.1599 -2.488 -1.8614 4.118e-42
#> time3  -1.4581  0.1544 -1.761 -1.1554 3.636e-21
#> time4  -2.9260  0.2453 -3.407 -2.4453 8.379e-33
#> time5  -1.2051  0.1706 -1.539 -0.8706 1.633e-12
#> time6  -0.9102  0.1860 -1.275 -0.5457 9.843e-07
#> 
#> $logor
#>    Estimate Std.Err    2.5%  97.5%   P-value
#> X1   0.9913  0.1179 0.76024 1.2223 4.100e-17
#> X2   0.6962  0.1162 0.46847 0.9238 2.064e-09
#> X3   0.3466  0.1159 0.11941 0.5738 2.788e-03
#> X4   0.3223  0.1151 0.09668 0.5478 5.111e-03
#> 
#> $or
#>    Estimate     2.5%    97.5%
#> X1 2.694610 2.138791 3.394874
#> X2 2.006032 1.597554 2.518953
#> X3 1.414239 1.126834 1.774950
#> X4 1.380231 1.101503 1.729490

out$ploglik
#> [1] -1676.456

if (test==1) {
### library(ordinal)
### out1 <- clm(fentry~X1+X2+X3+X4,data=ttpd)
### summary(out1)

# formula: fentry ~ X1 + X2 + X3 + X4
# data:    ttpd
# 
#  link  threshold nobs logLik   AIC     niter max.grad cond.H 
#  logit flexible  1000 -1676.46 3372.91 6(2)  1.17e-12 5.3e+02
# 
# Coefficients:
#    Estimate Std. Error z value Pr(>|z|)    
# X1  -0.9913     0.1171  -8.465  < 2e-16 ***
# X2  -0.6962     0.1156  -6.021 1.74e-09 ***
# X3  -0.3466     0.1150  -3.013  0.00259 ** 
# X4  -0.3223     0.1147  -2.810  0.00495 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Threshold coefficients:
#     Estimate Std. Error z value
# 0|1  -2.0064     0.1461 -13.733
# 1|2  -1.3940     0.1396  -9.984
# 2|3  -0.7324     0.1347  -5.435
# 3|4  -0.6266     0.1343  -4.667
# 4|5  -0.1814     0.1333  -1.361
# 5|6   0.2123     0.1342   1.582
}
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
#>  [1] cli_3.6.5           knitr_1.51          rlang_1.1.6        
#>  [4] xfun_0.55           textshaping_1.0.4   jsonlite_2.0.0     
#>  [7] listenv_0.10.0      future.apply_1.20.1 lava_1.8.2         
#> [10] htmltools_0.5.9     ragg_1.5.0          sass_0.4.10        
#> [13] rmarkdown_2.30      grid_4.5.2          evaluate_1.0.5     
#> [16] jquerylib_0.1.4     fastmap_1.2.0       numDeriv_2016.8-1.1
#> [19] yaml_2.3.12         mvtnorm_1.3-3       lifecycle_1.0.4    
#> [22] timereg_2.0.7       compiler_4.5.2      codetools_0.2-20   
#> [25] fs_1.6.6            htmlwidgets_1.6.4   Rcpp_1.1.0         
#> [28] future_1.68.0       lattice_0.22-7      systemfonts_1.3.1  
#> [31] digest_0.6.39       R6_2.6.1            parallelly_1.46.0  
#> [34] parallel_4.5.2      splines_4.5.2       Matrix_1.7-4       
#> [37] bslib_0.9.0         tools_4.5.2         globals_0.18.0     
#> [40] survival_3.8-3      pkgdown_2.2.0       cachem_1.1.0       
#> [43] desc_1.4.3
```
