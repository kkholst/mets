# Haplotype Discrete Survival Models

## Haplotype Analysis for discrete TTP

Cycle-specific logistic regression of haplo-type effects with known
haplo-type probabilities. Given observed genotype G and unobserved
haplotypes H we here mix out over the possible haplotypes using that
$P\left( H|G \right)$ is given as input.

$$\begin{aligned}
{S\left( t|x,G \right)} & {= E(S\left( t\left| x,H) \right|G \right) = \sum\limits_{h \in G}P\left( h|G \right)S\left( t|z,h \right)}
\end{aligned}$$ so survival can be computed by mixing out over possible
h given g.

Survival is based on logistic regression for the discrete hazard
function of the form $$\begin{aligned}
{\text{logit}\left( P\left( T = t|T > = t,x,h \right) \right)} & {= \alpha_{t} + x(h)beta}
\end{aligned}$$ where x(h) is a regression design of x and haplotypes
$h = \left( h_{1},h_{2} \right)$.

Simple binomial data can be fitted using this function.

For standard errors we assume that haplotype probabilities are known.

We are particularly interested in the types haplotypes:

``` r
types <- c("DCGCGCTCACG","DTCCGCTGACG","ITCAGTTGACG","ITCCGCTGAGG")

## some haplotypes frequencies for simulations 
data(haplo)
hapfreqs <- haplo$hapfreqs 
print(hapfreqs)
#>             index   haplotype     freq
#> DCGAGCTCACG     1 DCGAGCTCACG 0.010681
#> DCGCGCTCACG     2 DCGCGCTCACG 0.138387
#> DTGAGCTCACG     3 DTGAGCTCACG 0.000310
#> DTGAGCTCACA     4 DTGAGCTCACA 0.006800
#> DTGAGCTCGCG     5 DTGAGCTCGCG 0.034517
#> DTGACCTCACG     6 DTGACCTCACG 0.001336
#> DTGCGCTCACG     7 DTGCGCTCACG 0.009969
#> DTGCGCTCACA     8 DTGCGCTCACA 0.011833
#> DTGCGCTCGCG     9 DTGCGCTCGCG 0.302389
#> DTGCGCCCGCG    10 DTGCGCCCGCG 0.001604
#> DTGCCCTCACG    11 DTGCCCTCACG 0.003912
#> DTCAGCTGACG    12 DTCAGCTGACG 0.001855
#> DTCCGCTGACG    13 DTCCGCTGACG 0.103394
#> DTCCCCTGACG    14 DTCCCCTGACG 0.000310
#> ITCAGTTGACG    15 ITCAGTTGACG 0.048124
#> ITCCGCTGAGG    16 ITCCGCTGAGG 0.291273
#> ITCCGTTGACG    17 ITCCGTTGACG 0.031089
#> ITCCGTCGACG    18 ITCCGTCGACG 0.001502
#> ITCCCCTGAGG    19 ITCCCCTGAGG 0.000653
```

Among the types of interest we look up the frequencies and choose a
baseline

``` r
www <-which(hapfreqs$haplotype %in% types)
hapfreqs$freq[www]
#> [1] 0.138387 0.103394 0.048124 0.291273

baseline=hapfreqs$haplotype[9]
baseline
#> [1] "DTGCGCTCGCG"
```

We have cycle specific data with $id$ and outcome $y$

``` r
haploX  <- haplo$haploX
dlist(haploX,.~id|id %in% c(1,4,7))
#> id: 1
#>   y X1 X2 X3 X4 times lbnr__id Count1
#> 1 0 0  0  0  0  1     1        0     
#> 2 0 0  0  0  0  2     2        0     
#> 3 0 0  0  0  0  3     3        0     
#> 4 0 0  0  0  0  4     4        0     
#> 5 0 0  0  0  0  5     5        0     
#> 6 0 0  0  0  0  6     6        0     
#> ------------------------------------------------------------ 
#> id: 4
#>    y X1 X2 X3 X4 times lbnr__id Count1
#> 19 1 0  0  0  0  1     1        0     
#> ------------------------------------------------------------ 
#> id: 7
#>    y X1 X2 X3 X4 times lbnr__id Count1
#> 37 0 1  0  0  0  1     1        0     
#> 38 0 1  0  0  0  2     2        0     
#> 39 1 1  0  0  0  3     3        0
```

and a list of possible haplo-types for each id and how likely they are
$p$ (the sum of within each id is 1):

``` r
ghaplos <- haplo$ghaplos
head(ghaplos)
#>    id      haplo1      haplo2          p
#> 1   1 DTGCGCTCGCG DTGAGCTCGCG 1.00000000
#> 19  2 ITCCGTTGACG DTGAGCTCGCG 0.06867716
#> 21  2 ITCAGTTGACG DTGCGCTCGCG 0.93132284
#> 51  3 ITCCGTTGACG DTGAGCTCGCG 0.06867716
#> 53  3 ITCAGTTGACG DTGCGCTCGCG 0.93132284
#> 66  4 DTGCGCTCGCG DTGCGCTCGCG 1.00000000
```

The first id=1 has the haplotype fully observed, but id=2 has two
possible haplotypes consistent with the observed genotype for this id,
the probabiblities are 7% and 93%, respectively.

With the baseline given above we can specify a regression design that
gives an effect if a “type” is present (sm=0), or an additive effect of
haplotypes (sm=1):

``` r
designftypes <- function(x,sm=0) {
hap1=x[1]
hap2=x[2]
if (sm==0) y <- 1*( (hap1==types) | (hap2==types))
if (sm==1) y <- 1*(hap1==types) + 1*(hap2==types)
return(y)
}
```

To fit the model we start by constructing a time-design (named X) and
takes the haplotype distributions for each id

``` r
haploX$time <- haploX$times
Xdes <- model.matrix(~factor(time),haploX)
colnames(Xdes) <- paste("X",1:ncol(Xdes),sep="")
X <- dkeep(haploX,~id+y+time)
X <- cbind(X,Xdes)
Haplos <- dkeep(ghaplos,~id+"haplo*"+p)
desnames=paste("X",1:6,sep="")   # six X's related to 6 cycles 
head(X)
#>   id y time X1 X2 X3 X4 X5 X6
#> 1  1 0    1  1  0  0  0  0  0
#> 2  1 0    2  1  1  0  0  0  0
#> 3  1 0    3  1  0  1  0  0  0
#> 4  1 0    4  1  0  0  1  0  0
#> 5  1 0    5  1  0  0  0  1  0
#> 6  1 0    6  1  0  0  0  0  1
```

Now we can fit the model with the design given by the designfunction

``` r
out <- haplo.surv.discrete(X=X,y="y",time.name="time",
      Haplos=Haplos,desnames=desnames,designfunc=designftypes) 
names(out$coef) <- c(desnames,types)
out$coef
#>          X1          X2          X3          X4          X5          X6 
#> -1.82153345 -0.61608261 -0.17143057 -1.27152045 -0.28635976 -0.19349091 
#> DCGCGCTCACG DTCCGCTGACG ITCAGTTGACG ITCCGCTGAGG 
#>  0.79753613  0.65747412  0.06119231  0.31666905
summary(out)
#>             Estimate Std.Err     2.5%   97.5%   P-value
#> X1          -1.82153  0.1619 -2.13892 -1.5041 2.355e-29
#> X2          -0.61608  0.1895 -0.98748 -0.2447 1.149e-03
#> X3          -0.17143  0.1799 -0.52398  0.1811 3.406e-01
#> X4          -1.27152  0.2631 -1.78719 -0.7559 1.346e-06
#> X5          -0.28636  0.2030 -0.68425  0.1115 1.584e-01
#> X6          -0.19349  0.2134 -0.61184  0.2249 3.647e-01
#> DCGCGCTCACG  0.79754  0.1494  0.50465  1.0904 9.445e-08
#> DTCCGCTGACG  0.65747  0.1621  0.33971  0.9752 5.007e-05
#> ITCAGTTGACG  0.06119  0.2145 -0.35931  0.4817 7.755e-01
#> ITCCGCTGAGG  0.31667  0.1361  0.04989  0.5834 1.999e-02
```

Haplotypes “DCGCGCTCACG” “DTCCGCTGACG” gives increased hazard of
pregnancy

The data was generated with these true coefficients

``` r
tcoef=c(-1.93110204,-0.47531630,-0.04118204,-1.57872602,-0.22176426,-0.13836416,
0.88830288,0.60756224,0.39802821,0.32706859)

cbind(out$coef,tcoef)
#>                               tcoef
#> X1          -1.82153345 -1.93110204
#> X2          -0.61608261 -0.47531630
#> X3          -0.17143057 -0.04118204
#> X4          -1.27152045 -1.57872602
#> X5          -0.28635976 -0.22176426
#> X6          -0.19349091 -0.13836416
#> DCGCGCTCACG  0.79753613  0.88830288
#> DTCCGCTGACG  0.65747412  0.60756224
#> ITCAGTTGACG  0.06119231  0.39802821
#> ITCCGCTGAGG  0.31666905  0.32706859
```

The design fitted can be found in the output

``` r
head(out$X,10)
#>    X1 X2 X3 X4 X5 X6 haplo1 haplo2 haplo3 haplo4
#> 1   1  0  0  0  0  0      0      0      0      0
#> 2   1  1  0  0  0  0      0      0      0      0
#> 3   1  0  1  0  0  0      0      0      0      0
#> 4   1  0  0  1  0  0      0      0      0      0
#> 5   1  0  0  0  1  0      0      0      0      0
#> 6   1  0  0  0  0  1      0      0      0      0
#> 8   1  0  0  0  0  0      0      0      1      0
#> 10  1  1  0  0  0  0      0      0      1      0
#> 12  1  0  1  0  0  0      0      0      1      0
#> 14  1  0  0  1  0  0      0      0      1      0
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
