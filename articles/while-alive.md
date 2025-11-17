# While Alive estimands for Recurrent Events

## While Alive estimands for Recurrent Events

We consider two while-alive estimands for recurrent events data
$$\begin{array}{r}
\frac{E\left( N(D \land t) \right)}{E(D \land t)}
\end{array}$$ and the mean of the subject specific events per time-unit
$$\begin{array}{r}
{E\left( \frac{N(D \land t)}{D \land t} \right)}
\end{array}$$ for two treatment-groups in the case of an RCT. For the
mean of events per time-unit it has been seen that when the sample size
is small one can improve the finite sample properties by employing a
transformation such as square or cube-root, and thus consider
$$\begin{array}{r}
{E\left( \left( \frac{N(D \land t)}{D \land t} \right)^{.33} \right)}
\end{array}$$

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
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> mean events, E(N(min(D,t))): 
#>            Estimate Std.Err  2.5% 97.5%   P-value
#> treatment0    1.572 0.09573 1.384 1.759 1.375e-60
#> treatment1    1.453 0.10315 1.251 1.656 4.376e-45
#>  
#>                           Estimate Std.Err    2.5%  97.5% P-value
#> [treatment0] - [treat....   0.1185  0.1407 -0.1574 0.3943     0.4
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> _______________________________________________________ 
#> Ratio of means E(N(min(D,t)))/E(min(D,t)) 
#>            Estimate Std.Err   2.5%  97.5%   P-value
#> treatment0   0.8457 0.05264 0.7425 0.9488 4.411e-58
#> treatment1   0.7555 0.05433 0.6490 0.8619 5.963e-44
#>  
#>                           Estimate Std.Err     2.5%  97.5% P-value
#> [treatment0] - [treat....  0.09022 0.07565 -0.05805 0.2385   0.233
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> _______________________________________________________ 
#> Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
#>        Estimate Std.Err   2.5%  97.5%   P-value
#> treat0   1.0725  0.1222 0.8331 1.3119 1.645e-18
#> treat1   0.7552  0.0643 0.6291 0.8812 7.508e-32
#>  
#>                     Estimate Std.Err    2.5%  97.5% P-value
#> [treat0] - [treat1]   0.3173  0.1381 0.04675 0.5879 0.02153
#> 
#>  Null Hypothesis: 
#>   [treat0] - [treat1] = 0

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
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> mean events, E(N(min(D,t))): 
#>            Estimate Std.Err   2.5%  97.5%   P-value
#> treatment0   0.4523 0.06090 0.3329 0.5716 1.119e-13
#> treatment1   0.3739 0.07097 0.2348 0.5130 1.376e-07
#>  
#>                           Estimate Std.Err    2.5%  97.5% P-value
#> [treatment0] - [treat....  0.07835 0.09352 -0.1049 0.2616  0.4022
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> _______________________________________________________ 
#> Ratio of means E(N(min(D,t)))/E(min(D,t)) 
#>            Estimate Std.Err    2.5%    97.5%   P-value
#> treatment0  -0.1676 0.06224 -0.2896 -0.04563 7.081e-03
#> treatment1  -0.2804 0.07192 -0.4214 -0.13947 9.651e-05
#>  
#>                           Estimate Std.Err     2.5%  97.5% P-value
#> [treatment0] - [treat....   0.1128 0.09511 -0.07361 0.2992  0.2356
#> 
#>  Null Hypothesis: 
#>   [treatment0] - [treatment1] = 0 
#> _______________________________________________________ 
#> Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
#>        Estimate Std.Err    2.5%   97.5%   P-value
#> treat0  -0.3833 0.04939 -0.4801 -0.2865 8.487e-15
#> treat1  -0.5380 0.05666 -0.6491 -0.4270 2.191e-21
#>  
#>                     Estimate Std.Err     2.5%  97.5% P-value
#> [treat0] - [treat1]   0.1548 0.07517 0.007459 0.3021 0.03948
#> 
#>  Null Hypothesis: 
#>   [treat0] - [treat1] = 0
```

We see that the ratio of means are not very different, but that the
subject specific mean of events per time-unit shows that those on the
active treatment has fewer events per time-unit on average.

## Composite outcomes involving death and marks

The number of events can be generalized in various ways by using other
outcomes than $N(D \land t)$, for example,  
$$\begin{array}{r}
{\widetilde{N}(D \land t) = \int_{0}^{t}I(D \geq s)M(s)dN(s) + \sum\limits_{j}M_{j}I(D \leq t,\epsilon = j))}
\end{array}$$ where $M(s)$ are the marks related to $N(s)$ and are
$M_{j}$ marks associated with the different causes of the terminal
event. This provides an extension of the weighted composite outcomes
measure of Mao & Lin (2022).

The marks (or here weights) can be stochastic if we are couting hosptial
expenses, for example, and is vector on the data-frame. The marks for
the event times (defined through the causes) will then be used.

Here weighting death with weight 2 and otherwise couting the recurrent
of events as before (with weight 1)

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
