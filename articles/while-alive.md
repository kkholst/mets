# While Alive estimands for Recurrent Events

## While Alive estimands for Recurrent Events

We consider two while-alive estimands for recurrent events data
\begin{align\*} \frac{E(N(D \wedge t))}{E(D \wedge t)} \end{align\*} and
the mean of the subject specific events per time-unit \begin{align\*} E(
\frac{N(D \wedge t)}{D \wedge t} ) \end{align\*} for two
treatment-groups in the case of an RCT. For the mean of events per
time-unit it has been seen that when the sample size is small one can
improve the finite sample properties by employing a transformation such
as square or cube-root, and thus consider \begin{align\*} E( (\frac{N(D
\wedge t)}{D \wedge t})^.33 ) \end{align\*}

``` r

data(hfactioncpx12)

dtable(hfactioncpx12,~status)
#> 
#> status
#>    0    1    2 
#>  617 1391  124
dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,death.code=2)
summary(dd)
#> While-Alive summaries:  
#> 
#> RMST,  E(min(D,t)) 
#>            Estimate Std.Err  2.5% 97.5% P-value
#> treatment0    1.859 0.02108 1.817 1.900       0
#> treatment1    1.924 0.01502 1.894 1.953       0
#>  
#>                           Estimate Std.Err    2.5%    97.5% P-value
#> [treatment0] - [treat.... -0.06517 0.02588 -0.1159 -0.01444  0.0118
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#>  
#> chisq = 6.3405, df = 1, p-value = 0.0118
#> mean events, E(N(min(D,t))): 
#>            Estimate Std.Err  2.5% 97.5%   P-value
#> treatment0    1.572 0.09573 1.384 1.759 1.375e-60
#> treatment1    1.453 0.10315 1.251 1.656 4.376e-45
#>  
#>                           Estimate Std.Err    2.5%  97.5% P-value
#> [treatment0] - [treat....   0.1185  0.1407 -0.1574 0.3943     0.4
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#>  
#> chisq = 0.7085, df = 1, p-value = 0.4
#> _______________________________________________________ 
#> Ratio of means E(N(min(D,t)))/E(min(D,t)) 
#>    Estimate Std.Err   2.5%  97.5%   P-value
#> p1   0.8457 0.05264 0.7425 0.9488 4.411e-58
#> p2   0.7555 0.05433 0.6490 0.8619 5.963e-44
#>  
#>             Estimate Std.Err     2.5%  97.5% P-value
#> [p1] - [p2]  0.09022 0.07565 -0.05805 0.2385   0.233
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p1] - [p2] = 0 
#>  
#> chisq = 1.4222, df = 1, p-value = 0.233
#> _______________________________________________________ 
#> Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
#>        Estimate Std.Err   2.5%  97.5%   P-value
#> treat0   1.0725  0.1222 0.8331 1.3119 1.645e-18
#> treat1   0.7552  0.0643 0.6291 0.8812 7.508e-32
#>  
#>                     Estimate Std.Err    2.5%  97.5% P-value
#> [treat0] - [treat1]   0.3173  0.1381 0.04675 0.5879 0.02153
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treat0] - [treat1] = 0 
#>  
#> chisq = 5.2837, df = 1, p-value = 0.02153

dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
           death.code=2,trans=.333)
summary(dd,type="log")
#> While-Alive summaries, log-scale:  
#> 
#> RMST,  E(min(D,t)) 
#>            Estimate  Std.Err   2.5%  97.5% P-value
#> treatment0   0.6199 0.011340 0.5977 0.6421       0
#> treatment1   0.6543 0.007807 0.6390 0.6696       0
#>  
#>                           Estimate Std.Err     2.5%     97.5% P-value
#> [treatment0] - [treat.... -0.03446 0.01377 -0.06145 -0.007478 0.01231
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#>  
#> chisq = 6.2656, df = 1, p-value = 0.01231
#> mean events, E(N(min(D,t))): 
#>            Estimate Std.Err   2.5%  97.5%   P-value
#> treatment0   0.4523 0.06090 0.3329 0.5716 1.119e-13
#> treatment1   0.3739 0.07097 0.2348 0.5130 1.376e-07
#>  
#>                           Estimate Std.Err    2.5%  97.5% P-value
#> [treatment0] - [treat....  0.07835 0.09352 -0.1049 0.2616  0.4022
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#>  
#> chisq = 0.7018, df = 1, p-value = 0.4022
#> _______________________________________________________ 
#> Ratio of means E(N(min(D,t)))/E(min(D,t)) 
#>    Estimate Std.Err    2.5%    97.5%   P-value
#> p1  -0.1676 0.06224 -0.2896 -0.04563 7.081e-03
#> p2  -0.2804 0.07192 -0.4214 -0.13947 9.651e-05
#>  
#>             Estimate Std.Err     2.5%  97.5% P-value
#> [p1] - [p2]   0.1128 0.09511 -0.07361 0.2992  0.2356
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p1] - [p2] = 0 
#>  
#> chisq = 1.4067, df = 1, p-value = 0.2356
#> _______________________________________________________ 
#> Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
#>        Estimate Std.Err    2.5%   97.5%   P-value
#> treat0  -0.3833 0.04939 -0.4801 -0.2865 8.487e-15
#> treat1  -0.5380 0.05666 -0.6491 -0.4270 2.191e-21
#>  
#>                     Estimate Std.Err     2.5%  97.5% P-value
#> [treat0] - [treat1]   0.1548 0.07517 0.007459 0.3021 0.03948
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [treat0] - [treat1] = 0 
#>  
#> chisq = 4.2403, df = 1, p-value = 0.03948
```

We see that the ratio of means is not very different between groups, but
that the subject-specific mean of events per time-unit shows that those
on the active treatment have fewer events per time-unit on average.

We can also fit a regression model for the mean of the subject-specific
events per time-unit, here using the exponential link.

``` r

hfactioncpx12$age <- rnorm(741)[hfactioncpx12$id]
hfactioncpx12$sex <- rbinom(741,1,0.5)[hfactioncpx12$id]

dd <- WA_reg(Event(entry,time,status)~treatment+age+sex+cluster(id),hfactioncpx12,time=2,death.code=2)
summary(dd)
#>    n events
#>  741     86
#> 
#>  741 clusters
#> coeffients:
#>               Estimate    Std.Err       2.5%      97.5% P-value
#> (Intercept)  0.0878839  0.0937519 -0.0958665  0.2716342  0.3485
#> treatment1  -0.3522855  0.1368005 -0.6204094 -0.0841615  0.0100
#> age         -0.0036586  0.0499082 -0.1014768  0.0941597  0.9416
#> sex         -0.0316276  0.1479585 -0.3216210  0.2583658  0.8307
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  1.09186 0.90859 1.3121
#> treatment1   0.70308 0.53772 0.9193
#> age          0.99635 0.90350 1.0987
#> sex          0.96887 0.72497 1.2948
```

## Composite outcomes involving death and marks

The event count can be generalised in various ways by using outcomes
other than N(D \wedge t), for example, \begin{align\*} \tilde N(D \wedge
t) = \int_0^t I(D \geq s) M(s) dN(s) + \sum_j M_j I(D \leq t,\epsilon=j)
) \end{align\*} where M(s) are the marks associated with N(s) and M_j
are marks associated with the different causes of the terminal event.
This provides an extension of the weighted composite outcomes measure of
Mao & Lin (2022).

The marks (or weights) can be stochastic, for example when counting
hospital expenses, and are stored as a vector in the data frame. The
marks for the event times (defined through the causes) are then used.

Here we weight death by 2 and otherwise count recurrent events as before
(with weight 1):

``` r

hfactioncpx12$marks <- runif(nrow(hfactioncpx12))

##ddmg <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
##cause=1:2,death.code=2,marks=hfactioncpx12$marks)
##summary(ddmg)

ddm <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
cause=1:2,death.code=2,marks=hfactioncpx12$status)
```

## SessionInfo

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#> [1] mets_1.3.10
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6              knitr_1.51             rlang_1.2.0           
#>  [4] xfun_0.57              textshaping_1.0.5      jsonlite_2.0.0        
#>  [7] listenv_0.10.1         future.apply_1.20.2    lava_1.9.1            
#> [10] htmltools_0.5.9        ragg_1.5.2             sass_0.4.10           
#> [13] rmarkdown_2.31         grid_4.6.0             evaluate_1.0.5        
#> [16] jquerylib_0.1.4        fastmap_1.2.0          numDeriv_2016.8-1.1   
#> [19] yaml_2.3.12            mvtnorm_1.3-7          lifecycle_1.0.5       
#> [22] timereg_2.0.7          compiler_4.6.0         codetools_0.2-20      
#> [25] fs_2.1.0               htmlwidgets_1.6.4      Rcpp_1.1.1-1.1        
#> [28] future_1.70.0          lattice_0.22-9         systemfonts_1.3.2     
#> [31] digest_0.6.39          R6_2.6.1               parallelly_1.47.0     
#> [34] parallel_4.6.0         splines_4.6.0          Matrix_1.7-5          
#> [37] bslib_0.11.0           tools_4.6.0            RcppArmadillo_15.2.6-1
#> [40] globals_0.19.1         survival_3.8-6         pkgdown_2.2.0         
#> [43] cachem_1.1.0           desc_1.4.3
```
