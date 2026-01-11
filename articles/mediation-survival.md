# Mediation Analysis for survival data

## Overview

Fit

- binomial-regression IPCW, binreg
- additive Lin-Ying model, aalenMets
- cox model phreg
- standard logistic regression via binreg

in the context of mediation analysis using mediation weights as in the
medFlex package. We thus fit natural effects models, that for example on
the binary scale might state that $$\begin{array}{r}
{\text{logit}(P\left( Y\left( x,M\left( x^{*} \right) \right) = 1|Z \right) = \beta_{0} + \beta_{1}x + \beta_{2}x^{*} + \beta_{3}^{T}Z,}
\end{array}$$ in this case the the Natural Direct Effect (NDE) for fixed
covariates $Z$ is $$\begin{array}{r}
{\text{OR}_{1,0|Z}^{\text{NDE}} = \frac{\text{odds}\left( Y\left( 1,M(x) \right)|Z \right)}{\text{odds}\left( Y\left( 0,M(x) \right)|Z \right)} = \exp\left( \beta_{1} \right),}
\end{array}$$ and the Natural Inderect Effect (NIE) for fixed covariates
$Z$ is $$\begin{array}{r}
{\text{OR}_{1,0|Z}^{\text{NIE}} = \frac{\text{odds}\left( Y\left( x,M(1) \right)|Z \right)}{\text{odds}\left( Y\left( x,M(0) \right)|Z \right)} = \exp\left( \beta_{2} \right).}
\end{array}$$ See the medFlex package for additional discussion of the
parametrization.

The mediator can be

- binomial using glm-binomial.
- multnomial via the mlogit function of mets

Both mediator and exposure must be coded as factors.

In the below example these are

- mediator: gp.f
- exposure : dnr.f

and the outcome model is concerned with the risk/hazard of cause=2.

The key is that the standard errors are computed using the i.i.d
influence functions and a Taylor expansion to deal with the uncertainty
from the mediation weights.

## Simulated Data

First we simulate some data that mimics that of Kumar et al 2012. This
is data from multiple myeloma patients treated with allogeneic stem cell
transplantation from the Center for International Blood and Marrow
Transplant Research (CIBMTR) Kumar et al (2012), “Trends in allogeneic
stem cell transplantation for multiple myeloma: a CIBMTR analysis”. The
data used in this paper consist of patients transplanted from 1995 to
2005, and we compared the outcomes between transplant periods: 2001-2005
(N=488) versus 1995-2000 (N=375). The two competing events were relapse
(cause 2) and treatment-related mortality (TRM, cause 1)) defined as
death without relapse. considered the following risk covariates:
transplant time period (gp (main interest of the study): 1 for
transplanted in 2001-2005 versus 0 for transplanted in 1995-2000), donor
type (dnr: 1 for Unrelated or other related donor (N=280) versus 0 for
HLA-identical sibling (N=584)), prior autologous transplant (preauto: 1
for Auto+Allo transplant (N=399) versus 0 for allogeneic transplant
alone (N=465)) and time to transplant (ttt24: 1 for more than 24 months
(N=289) versus 0 for less than or equal to 24 months (N=575))).

The interest is then on the effect of the period (gp) and the possible
mediation via the amount of unrealted or related donors (dnr). A
somewhat artificial example ! All adjusted for other important
counfounders.

``` r
 library(mets)
 runb <- 0
 options(warn=-1)
 set.seed(1000) # to control output in simulatins for p-values below.

n <- 200; k.boot <- 10; 

dat <- kumarsimRCT(n,rho1=0.5,rho2=0.5,rct=2,censpar=c(0,0,0,0),
          beta = c(-0.67, 0.59, 0.55, 0.25, 0.98, 0.18, 0.45, 0.31),
    treatmodel = c(-0.18, 0.56, 0.56, 0.54),restrict=1)
dfactor(dat) <- dnr.f~dnr
dfactor(dat) <- gp.f~gp
drename(dat) <- ttt24~"ttt24*"
dat$id <- 1:n
dat$ftime <- 1
```

## Mediation Weights

Then compute the mediation weights based on a mediation model

``` r
weightmodel <- fit <- glm(gp.f~dnr.f+preauto+ttt24,data=dat,family=binomial)
wdata <- medweight(fit,data=dat)
```

## Binomial Regression

A simple multvariate regression of the probaibility of relapse at 50
months with both exposure and mediator (given the other covariates)

``` r
aaMss2 <- binreg(Event(time,status)~gp+dnr+preauto+ttt24+cluster(id),data=dat,time=50,cause=2)
summary(aaMss2)
#>    n events
#>  200     97
#> 
#>  200 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept) -1.01508  0.31869 -1.63971 -0.39046  0.0014
#> gp           1.08533  0.34216  0.41471  1.75594  0.0015
#> dnr          0.51969  0.35757 -0.18113  1.22051  0.1461
#> preauto      0.39417  0.35936 -0.31017  1.09851  0.2727
#> ttt24        0.50469  0.38681 -0.25344  1.26283  0.1920
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.36237 0.19404 0.6767
#> gp           2.96041 1.51394 5.7889
#> dnr          1.68151 0.83433 3.3889
#> preauto      1.48316 0.73332 2.9997
#> ttt24        1.65648 0.77612 3.5354
```

## Binomial regression IPCW Mediation Analysis

We first look at the probability of relapse at 50 months

``` r
### binomial regression ###########################################################
aaMss <- binreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
        time=50,weights=wdata$weights,cause=2)
summary(aaMss)
#>    n events
#>  400    194
#> 
#>  200 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.535534  0.256218 -1.037712 -0.033356  0.0366
#> dnr.f01      0.375817  0.348618 -0.307462  1.059095  0.2810
#> dnr.f11      0.275383  0.071199  0.135836  0.414931  0.0001
#> preauto      0.588221  0.350437 -0.098624  1.275066  0.0932
#> ttt24        0.266179  0.363603 -0.446469  0.978827  0.4641
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.58536 0.35426 0.9672
#> dnr.f01      1.45618 0.73531 2.8838
#> dnr.f11      1.31704 1.14549 1.5143
#> preauto      1.80078 0.90608 3.5789
#> ttt24        1.30497 0.63988 2.6613

ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  400    194
#> 
#>  200 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.535534  0.254832 -1.034995 -0.036073  0.0356
#> dnr.f01      0.375817  0.317732 -0.246927  0.998560  0.2369
#> dnr.f11      0.275383  0.117175  0.045726  0.505041  0.0188
#> preauto      0.588221  0.346523 -0.090951  1.267394  0.0896
#> ttt24        0.266179  0.366361 -0.451875  0.984233  0.4675
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.58536 0.35523 0.9646
#> dnr.f01      1.45618 0.78120 2.7144
#> dnr.f11      1.31704 1.04679 1.6571
#> preauto      1.80078 0.91306 3.5516
#> ttt24        1.30497 0.63643 2.6758
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}
```

So the NDE is $1.40(0.72,2.76)$ and the NIE is $1.32(1.05,1.66)$.

## Mediation Analysis

We here also illustrate how to use the other models mentioned above.

``` r
### lin-ying model ################################################################
aaMss <- aalenMets(Surv(time/100,status==2)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
           weights=wdata$weights)
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  400    196
#> coeffients:
#>          Estimate   Std.Err      2.5%     97.5% P-value
#> dnr.f01  1.169592  0.739323 -0.279454  2.618637  0.1137
#> dnr.f11  0.206757  0.131289 -0.050565  0.464078  0.1153
#> preauto  0.617537  0.504302 -0.370877  1.605950  0.2207
#> ttt24    0.457736  0.517822 -0.557175  1.472648  0.3767
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### cox model ###############################################################################
aaMss <- phreg(Surv(time,status==2)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
           weights=wdata$weights)
summary(aaMss)
#> 
#>    n events
#>  400    196
#> coeffients:
#>         Estimate     S.E.  dU^-1/2 P-value
#> dnr.f01 0.414565 0.213724 0.157231  0.0524
#> dnr.f11 0.100656 0.039308 0.144971  0.0104
#> preauto 0.284460 0.232166 0.162375  0.2205
#> ttt24   0.185561 0.226044 0.160886  0.4117
#> 
#> exp(coeffients):
#>         Estimate    2.5%  97.5%
#> dnr.f01  1.51371 0.99568 2.3013
#> dnr.f11  1.10590 1.02389 1.1945
#> preauto  1.32904 0.84318 2.0949
#> ttt24    1.20389 0.77300 1.8750
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  400    196
#> coeffients:
#>            Estimate     Std.Err        2.5%       97.5% P-value
#> dnr.f01  0.41456472  0.20869639  0.00552731  0.82360212  0.0470
#> dnr.f11  0.10065575  0.05121458  0.00027702  0.20103448  0.0494
#> preauto  0.28445952  0.23037280 -0.16706288  0.73598192  0.2169
#> ttt24    0.18556110  0.22549763 -0.25640614  0.62752835  0.4106
#> 
#> exp(coeffients):
#>         Estimate    2.5%  97.5%
#> dnr.f01  1.51371 1.00554 2.2787
#> dnr.f11  1.10590 1.00028 1.2227
#> preauto  1.32904 0.84615 2.0875
#> ttt24    1.20389 0.77383 1.8730
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### Fine-Gray #############################################################3
aaMss <- cifreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
        weights=wdata$weights,propodds=NULL,cause=2)
summary(aaMss)
#> 
#>    n events
#>  400    196
#> 
#>  200 clusters
#> coeffients:
#>         Estimate    S.E. dU^-1/2 P-value
#> dnr.f01  0.18943 0.21986 0.15855  0.3889
#> dnr.f11  0.18730 0.04083 0.14503  0.0000
#> preauto  0.41452 0.22783 0.16098  0.0688
#> ttt24    0.17304 0.22892 0.16308  0.4497
#> 
#> exp(coeffients):
#>         Estimate    2.5%  97.5%
#> dnr.f01  1.20856 0.78545 1.8596
#> dnr.f11  1.20599 1.11324 1.3065
#> preauto  1.51364 0.96849 2.3656
#> ttt24    1.18892 0.75910 1.8621
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  400    196
#> 
#>  200 clusters
#> coeffients:
#>          Estimate   Std.Err      2.5%     97.5% P-value
#> dnr.f01  0.189426  0.233939 -0.269087  0.647939  0.4181
#> dnr.f11  0.187298  0.047733  0.093744  0.280853  0.0001
#> preauto  0.414517  0.230676 -0.037600  0.866634  0.0723
#> ttt24    0.173042  0.230810 -0.279338  0.625422  0.4534
#> 
#> exp(coeffients):
#>         Estimate    2.5%  97.5%
#> dnr.f01  1.20856 0.76408 1.9116
#> dnr.f11  1.20599 1.09828 1.3243
#> preauto  1.51364 0.96310 2.3789
#> ttt24    1.18892 0.75628 1.8690
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### logit model  #############################################################3
aaMss <- cifreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
        weights=wdata$weights,cause=2)
summary(aaMss)
#> 
#>    n events
#>  400    196
#> 
#>  200 clusters
#> coeffients:
#>         Estimate     S.E.  dU^-1/2 P-value
#> dnr.f01 0.357168 0.339848 0.158937  0.2933
#> dnr.f11 0.272392 0.064166 0.145076  0.0000
#> preauto 0.657010 0.326082 0.160361  0.0439
#> ttt24   0.191333 0.353606 0.167443  0.5884
#> 
#> exp(coeffients):
#>         Estimate    2.5%  97.5%
#> dnr.f01  1.42928 0.73424 2.7822
#> dnr.f11  1.31310 1.15792 1.4891
#> preauto  1.92902 1.01806 3.6551
#> ttt24    1.21086 0.60549 2.4215
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  400    196
#> 
#>  200 clusters
#> coeffients:
#>          Estimate   Std.Err      2.5%     97.5% P-value
#> dnr.f01  0.357168  0.351089 -0.330953  1.045289  0.3090
#> dnr.f11  0.272392  0.068131  0.138857  0.405927  0.0001
#> preauto  0.657010  0.328207  0.013736  1.300284  0.0453
#> ttt24    0.191333  0.356086 -0.506583  0.889250  0.5910
#> 
#> exp(coeffients):
#>         Estimate    2.5%  97.5%
#> dnr.f01  1.42928 0.71824 2.8442
#> dnr.f11  1.31310 1.14896 1.5007
#> preauto  1.92902 1.01383 3.6703
#> ttt24    1.21086 0.60255 2.4333
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### binomial outcome  ############################
aaMss <- binreg(Event(ftime,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
        time=50,weights=wdata$weights,cens.weights=1,cause=2)
summary(aaMss)
#>    n events
#>  400    196
#> 
#>  200 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.674433  0.235285 -1.135583 -0.213284  0.0042
#> dnr.f01      0.221834  0.318264 -0.401952  0.845620  0.4858
#> dnr.f11      0.262722  0.060281  0.144572  0.380871  0.0000
#> preauto      0.578077  0.319091 -0.047331  1.203484  0.0700
#> ttt24        0.214442  0.328183 -0.428784  0.857669  0.5135
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.50944 0.32123 0.8079
#> dnr.f01      1.24836 0.66901 2.3294
#> dnr.f11      1.30046 1.15555 1.4636
#> preauto      1.78261 0.95377 3.3317
#> ttt24        1.23917 0.65130 2.3577
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  400    196
#> 
#>  200 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.674433  0.235022 -1.135069 -0.213798  0.0041
#> dnr.f01      0.221834  0.286717 -0.340122  0.783789  0.4391
#> dnr.f11      0.262722  0.107508  0.052011  0.473432  0.0145
#> preauto      0.578077  0.315260 -0.039822  1.195975  0.0667
#> ttt24        0.214442  0.329107 -0.430596  0.859480  0.5147
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.50944 0.32140 0.8075
#> dnr.f01      1.24836 0.71168 2.1898
#> dnr.f11      1.30046 1.05339 1.6055
#> preauto      1.78261 0.96096 3.3068
#> ttt24        1.23917 0.65012 2.3619
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}
```

## Multinomial regression

Also works with mediator with more than two levels

- meditor: wmi in 4 categories
- exposure: age in 4 categories

``` r
library(mets)
data(tTRACE)
dcut(tTRACE) <- ~. 

weightmodel <- fit <- mlogit(wmicat.4 ~agecat.4+vf+chf,data=tTRACE,family=binomial)
wdata <- medweight(fit,data=tTRACE)

aaMss <- binreg(Event(time,status)~agecat.40+ agecat.41+ vf+chf+cluster(id),data=wdata,
        time=7,weights=wdata$weights,cause=9)
summary(aaMss)
MultMed <- mediatorSurv(aaMss,fit,data=tTRACE,wdata=wdata)
summary(MultMed)
```

``` r
summary(MultMed)
#>     n events
#>  4000   2016
#> 
#>  1000 clusters
#> coeffients:
#>                       Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)          -1.839306  0.174541 -2.181401 -1.497211  0.0000
#> agecat.40(60.1,68.6]  0.916646  0.223488  0.478618  1.354674  0.0000
#> agecat.40(68.6,75.6]  1.363830  0.222418  0.927898  1.799762  0.0000
#> agecat.40(75.6,96.3]  2.277415  0.249815  1.787786  2.767044  0.0000
#> agecat.41(60.1,68.6]  0.121100  0.053334  0.016567  0.225633  0.0232
#> agecat.41(68.6,75.6]  0.119374  0.053193  0.015118  0.223631  0.0248
#> agecat.41(75.6,96.3]  0.095356  0.053874 -0.010234  0.200947  0.0767
#> vf                    0.712461  0.293627  0.136962  1.287960  0.0152
#> chf                   1.166578  0.154721  0.863331  1.469825  0.0000
#> 
#> exp(coeffients):
#>                      Estimate    2.5%   97.5%
#> (Intercept)           0.15893 0.11288  0.2238
#> agecat.40(60.1,68.6]  2.50089 1.61384  3.8755
#> agecat.40(68.6,75.6]  3.91114 2.52919  6.0482
#> agecat.40(75.6,96.3]  9.75144 5.97621 15.9115
#> agecat.41(60.1,68.6]  1.12874 1.01671  1.2531
#> agecat.41(68.6,75.6]  1.12679 1.01523  1.2506
#> agecat.41(75.6,96.3]  1.10005 0.98982  1.2226
#> vf                    2.03900 1.14678  3.6254
#> chf                   3.21099 2.37105  4.3485
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
