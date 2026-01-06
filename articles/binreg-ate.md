# Average treatment effect (ATE) for Competing risks and binary outcomes

The binregATE function can fit a logistic link model with IPCW
adjustment for a specific time-point, and can thus be used for
describing survival or competing risks data. The function can be used
for large data and is completely scalable, that is, linear in the data.
A nice feature is that influcence functions are computed and are
available.

In addition and to summarize

- the censoring weights can be strata dependent
- predictions can be computed with standard errors
- computation time is linear in data
  - including standard errors
- clusters can be given and then cluster corrected standard errors are
  computed

## Average treatment effect

First we simulate some data that mimics that of Kumar et al 2012. This
is data from multiple myeloma patients treated with allogeneic stem cell
transplantation from the Center for International Blood and Marrow
Transplant Research (CIBMTR) Kumar et al (2012), “Trends in allogeneic
stem cell transplantation for multiple myeloma: a CIBMTR analysis”. The
data used in this paper consist of patients transplanted from 1995 to
2005, and we compared the outcomes between transplant periods: 2001-2005
(N=488) versus 1995-2000 (N=375). The two competing events were relapse
(cause 2) and treatment-related mortality (TRM, cause 1) defined as
death without relapse. considered the following risk covariates:
transplant time period (gp (main interest of the study): 1 for
transplanted in 2001-2005 versus 0 for transplanted in 1995-2000), donor
type (dnr: 1 for Unrelated or other related donor (N=280) versus 0 for
HLA-identical sibling (N=584)), prior autologous transplant (preauto: 1
for Auto+Allo transplant (N=399) versus 0 for allogeneic transplant
alone (N=465)) and time to transplant (ttt24: 1 for more than 24 months
(N=289) versus 0 for less than or equal to 24 months (N=575))).

We here generate similar data by assuming that the two cumlative
incidence curves are logistic and we have censoring that depends on the
covariates via a Cox model. All this is wrapped in the kumarsim
function. The simulation does not deal with possible violations of the
bound that $F_{1} + F_{2} < 1$. But as we increase the sample size we
still see that we recover the parameters of cause 2.

``` r
library(mets) 
set.seed(100)
###
n <- 400
kumar <- kumarsim(n,depcens=1)
kumar$cause <- kumar$status
kumar$ttt24 <- kumar[,6]
dtable(kumar,~cause)
#> 
#> cause
#>   0   1   2 
#> 182  72 146
dfactor(kumar) <- gp.f~gp
kumar$id <- 1:n
kumar$idc <- sample(100,n,TRUE)
kumar$ids <- sample(n,n)
kumar$id2 <- sample(n,n)
kumar2 <- kumar[order(kumar$id2),]
kumar$int <- interaction(kumar$gp,kumar$dnr)
kumar2$int <- interaction(kumar2$gp,kumar2$dnr)
clust <- 0

b2 <- binregATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,
    treat.model=gp.f~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr))
summary(b2)
#>    n events
#>  400    137
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.310227  0.221121 -1.743616 -0.876838  0.0000
#> gp.f1        0.898030  0.271691  0.365526  1.430535  0.0009
#> dnr          0.323059  0.271545 -0.209160  0.855278  0.2342
#> preauto      0.269177  0.278181 -0.276049  0.814402  0.3332
#> ttt24        0.496202  0.276045 -0.044836  1.037240  0.0722
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.26976 0.17489 0.4161
#> gp.f1        2.45476 1.44127 4.1809
#> dnr          1.38135 0.81127 2.3520
#> preauto      1.30889 0.75878 2.2578
#> ttt24        1.64247 0.95615 2.8214
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    0.293021 0.039148 0.216292 0.369750   0e+00
#> treat1    0.496481 0.041744 0.414665 0.578297   0e+00
#> treat:1-0 0.203460 0.060294 0.085285 0.321635   7e-04
#> 
#> Average Treatment effects (double robust) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    0.311984 0.042914 0.227875 0.396094   0e+00
#> treat1    0.512310 0.042454 0.429103 0.595518   0e+00
#> treat:1-0 0.200326 0.060585 0.081582 0.319070   9e-04

b5 <- binregATE(Event(time,cause)~int+preauto+ttt24,kumar,cause=2,
        treat.model=int~preauto+ttt24,cens.code=0,time=60)
summary(b5)
#>    n events
#>  400    142
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.115768  0.285147 -1.674646 -0.556890  0.0001
#> int1.0       0.674582  0.333967  0.020018  1.329146  0.0434
#> int0.1       0.034751  0.494229 -0.933920  1.003421  0.9439
#> int1.1       1.340278  0.431989  0.493594  2.186961  0.0019
#> preauto      0.298967  0.272324 -0.234779  0.832713  0.2723
#> ttt24        0.476755  0.286955 -0.085667  1.039177  0.0966
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.32766 0.18737 0.5730
#> int1.0       1.96321 1.02022 3.7778
#> int0.1       1.03536 0.39301 2.7276
#> int1.1       3.82010 1.63819 8.9081
#> preauto      1.34846 0.79075 2.2995
#> ttt24        1.61084 0.91790 2.8269
#> 
#> Average Treatment effects (G-formula) :
#>                 Estimate    Std.Err       2.5%      97.5% P-value
#> treat0.0       0.3124362  0.0586931  0.1973999  0.4274724  0.0000
#> treat1.0       0.4676985  0.0422018  0.3849846  0.5504125  0.0000
#> treat0.1       0.3197899  0.0896402  0.1440983  0.4954815  0.0004
#> treat1.1       0.6274521  0.0726051  0.4851488  0.7697554  0.0000
#> treat:1.0-0.0  0.1552624  0.0738341  0.0105502  0.2999745  0.0355
#> treat:0.1-0.0  0.0073537  0.1048541 -0.1981565  0.2128640  0.9441
#> treat:1.1-0.0  0.3150160  0.0978046  0.1233225  0.5067094  0.0013
#> treat:0.1-1.0 -0.1479086  0.1003040 -0.3445009  0.0486837  0.1403
#> treat:1.1-1.0  0.1597536  0.0845518 -0.0059649  0.3254721  0.0588
#> treat:1.1-0.1  0.3076622  0.1178093  0.0767602  0.5385643  0.0090
#> 
#> Average Treatment effects (double robust) :
#>                 Estimate    Std.Err       2.5%      97.5% P-value
#> treat0.0       0.3580588  0.0674640  0.2258318  0.4902858  0.0000
#> treat1.0       0.4871368  0.0501492  0.3888461  0.5854274  0.0000
#> treat0.1       0.3102783  0.1210473  0.0730299  0.5475266  0.0104
#> treat1.1       0.7645588  0.1887893  0.3945387  1.1345790  0.0001
#> treat:1.0-0.0  0.1290779  0.0831680 -0.0339284  0.2920843  0.1207
#> treat:0.1-0.0 -0.0477806  0.1400295 -0.3222333  0.2266722  0.7329
#> treat:1.1-0.0  0.4065000  0.2040550  0.0065596  0.8064404  0.0464
#> treat:0.1-1.0 -0.1768585  0.1315359 -0.4346641  0.0809471  0.1788
#> treat:1.1-1.0  0.2774221  0.1937486 -0.1023182  0.6571624  0.1522
#> treat:1.1-0.1  0.4542806  0.2102654  0.0421680  0.8663931  0.0307
```

We note that the estimates found using the large censoring model are
very different from those using the simple Kaplan-Meier weights that are
severely biased for these data. This is due to a strong censoring
dependence.

The average treatment is around $0.17 = E\left( Y(1) - Y(0) \right)$ at
time 60 for the transplant period, under the standard causal
assumptions. The 1/0 treatment variable used for the causal computation
is found as the right hand side (rhs) of the treat.model or as the first
argument on the rhs of the response model.

The binregATE default uses binreg with its default to fit the working
model and is recommended, but the logitIPCW and logitIPCWATE can also be
used and are GLM-type IPCW weighted models (see binreg help
page/vignette).

``` r
ib2 <- logitIPCWATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,
    treat.model=gp.f~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr))
summary(ib2)
#>    n events
#>  400    137
#> 
#>  400 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept) -1.27557  0.21917 -1.70514 -0.84601  0.0000
#> gp.f1        0.83337  0.25904  0.32566  1.34108  0.0013
#> dnr          0.36815  0.27656 -0.17390  0.91020  0.1831
#> preauto      0.48521  0.30772 -0.11791  1.08832  0.1148
#> ttt24        0.18958  0.30947 -0.41698  0.79614  0.5401
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.27927 0.18175 0.4291
#> gp.f1        2.30105 1.38494 3.8232
#> dnr          1.44506 0.84038 2.4848
#> preauto      1.62452 0.88878 2.9693
#> ttt24        1.20874 0.65904 2.2170
#> 
#> Average Treatment effects (G-formula) :
#>         Estimate  Std.Err     2.5%    97.5% P-value
#> treat-1 0.493059 0.040018 0.414625 0.571494  0.0000
#> treat-0 0.302695 0.040274 0.223759 0.381632  0.0000
#> p1      0.190364 0.058248 0.076200 0.304527  0.0011
#> 
#> Average Treatment effects (double robust) :
#>         Estimate  Std.Err     2.5%    97.5% P-value
#> treat-1 0.524459 0.043164 0.439860 0.609058   0e+00
#> treat-0 0.309211 0.043507 0.223938 0.394483   0e+00
#> p1      0.215248 0.061268 0.095165 0.335331   4e-04

ib5 <- logitIPCW(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,cens.code=0,
     time=60,cens.model=~strata(gp,dnr))
summary(ib5)
#>    n events
#>  400    142
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.142446  0.225267 -1.583960 -0.700931  0.0000
#> gp.f1        0.716372  0.267369  0.192339  1.240405  0.0074
#> dnr          0.551548  0.325499 -0.086418  1.189515  0.0902
#> preauto      0.748078  0.343802  0.074239  1.421917  0.0296
#> ttt24       -0.184844  0.373708 -0.917298  0.547610  0.6209
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.31904 0.20516 0.4961
#> gp.f1        2.04699 1.21208 3.4570
#> dnr          1.73594 0.91721 3.2855
#> preauto      2.11294 1.07706 4.1451
#> ttt24        0.83123 0.39960 1.7291

ibs <- logitIPCW(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,cens.code=0,time=60)
summary(ibs)
#>    n events
#>  400    142
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.294250  0.220215 -1.725864 -0.862636  0.0000
#> gp.f1        1.633772  0.294651  1.056267  2.211277  0.0000
#> dnr          0.020551  0.309235 -0.585538  0.626639  0.9470
#> preauto      0.547275  0.305442 -0.051381  1.145931  0.0732
#> ttt24        0.288357  0.309124 -0.317515  0.894228  0.3509
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.27410 0.17802 0.4220
#> gp.f1        5.12316 2.87562 9.1274
#> dnr          1.02076 0.55681 1.8713
#> preauto      1.72854 0.94992 3.1454
#> ttt24        1.33423 0.72796 2.4454

check <- 0
if (check==1) {
require(riskRegression)
e.wglm <- wglm( regressor.event=~gp.f+dnr+preauto+ttt24, formula.censor = Surv(time,cause==0)~+1, times = 60, data = kumar, product.limit=TRUE,cause=2)
summary(e.wglm)$coef
estimate(ibs)

es.wglm <- wglm( regressor.event=~gp.f+dnr+preauto+ttt24, 
formula.censor = Surv(time,cause==0)~strata(gp,dnr), times = 60, 
data = kumar, product.limit=TRUE,cause=2)
summary(es.wglm)$coef
estimate(ib5)
}
```

The cluster argument should not be used for the logitIPCWATE, but works
for binregATE.

## Average treatment for Competing risks data

The binreg function does direct binomial regression for one time-point,
$t$, fitting the model $$\begin{aligned}
{P\left( T \leq t,\epsilon = 1|X \right)} & {= \text{expit}\left( X^{T}\beta \right),}
\end{aligned}$$ for possible right censored data. The estimation
procedure is based on IPCW adjusted estimating equation (EE)
$$\begin{aligned}
{U(\beta) =} & {X\left( \Delta(t)I(T \leq t,\epsilon = 1)/G_{c}(T \land t) - \text{expit}\left( X^{T}beta \right) \right) = 0}
\end{aligned}$$ where $G_{c}(t) = P(C > t)$, the censoring survival
distribution, and with $\Delta(t) = I(C > T \land t)$ the indicator of
being uncensored at time $t$.

The function logitIPCW instead considers the EE the EE $$\begin{aligned}
{U(\beta) =} & {X\frac{\Delta(t)}{G_{c}(T \land t)}\left( I(T \leq t,\epsilon = 1) - \text{expit}\left( X^{T}beta \right) \right) = 0.}
\end{aligned}$$ The two score equations are quite similar, and exactly
the same when the censoring model is fully-nonparametric given $X$.

- It seems that the binreg estimating equations most often is preferable
  to use, and the estimating equation used is also augmented in the
  default implementation (see the binreg vignette).

Additional functions logitATE, and binregATE computes the average
treatment effect. We demonstrate their use below.

The functions binregATE (recommended) and logitATE also works when there
is no censoring and we thus have simple binary outcome.

Variance is based on sandwich formula with IPCW adjustment, and
naive.var is variance under a known censoring model. The influence
functions are stored in the output. Further, the standard errors can be
cluster corrected by specifying the relevant cluster for the working
outcome model.

- We estimate the average treatment effect of our binary response
  $I(T \leq t,\epsilon = 1)$
  - Using a working logistic model for the resonse (possibly with a
    cluster specification)
  - Using a working logistic model for treatment given covariates
    - The binregATE can also handle a factor with more than two levels
      and then uses the mlogit multinomial regression function (of
      mets).
  - Using a working model for censoring given covariates, this must be a
    stratified Kaplan-Meier.

If there are no censoring then the censoring weights are simply set to
1.

The average treatment effect is $$\begin{array}{r}
{E\left( Y(1) - Y(0) \right)}
\end{array}$$ using counterfactual outcomes.

We compute the simple G-estimator $$\begin{array}{r}
{\sum m_{a}\left( X_{i} \right)}
\end{array}$$ to estimate the risk $E\left( Y(a) \right)$.

The DR-estimator instead uses the estimating equations that are double
robust wrt

- A working logistic model for the resonse
- A working logistic model for treatment given covariates

This is estimated using the estimator $$\begin{array}{r}
{\sum\left\lbrack \frac{A_{i}Y_{i}}{\pi_{A}\left( X_{i} \right)} - \frac{A_{i} - \pi_{A}\left( X_{i} \right)}{\pi_{A}\left( X_{i} \right)}m_{1}\left( X_{i} \right) \right\rbrack - \left\lbrack \frac{\left( 1 - A_{i} \right)Y_{i}}{1 - \pi_{A}\left( X_{i} \right)} + \frac{A_{i} - \pi_{A}\left( X_{i} \right)}{1 - \pi_{A}\left( X_{i} \right)}m_{0}\left( X_{i} \right) \right\rbrack}
\end{array}$$ where

- $A_{i}$ is treatment indicator
- $\pi_{A}\left( X_{i} \right) = P\left( A_{i} = 1|X_{i} \right)$ is
  treatment model
- $Y_{i}$ outcome, that in case of censoring is censoring adjusted
  ${\widetilde{Y}}_{i}\Delta(t)/G_{c}\left( T_{i} - \land t \right)$
- ${\widetilde{Y}}_{i} = I\left( T_{i} \leq t,\epsilon_{i} = 1 \right)$
  oucome before censoring.
- $m_{j}\left( X_{i} \right) = P\left( Y_{i} = 1|A_{i} = j,X_{i} \right)$
  is outcome model, using binomial regression.

The standard errors are then based on an iid decomposition using
taylor-expansions for the parameters of the treatment-model and the
outcome-model, and the censoring probability.

We need that the censoring model is correct, so it can be important to
use a sufficiently large censorng model as we also illustrate below.

- The censoring model can be specified by strata (used for phreg

We also compute standard marginalization for average treatment effect
(called differenceG) $$\begin{array}{r}
{\sum\left\lbrack m_{1}\left( X_{i} \right) - m_{0}\left( X_{i} \right) \right\rbrack}
\end{array}$$ and again standard errors are based on the related
influcence functions and are also returned.

For large data where there are more than 2 treatment groups the
computations can be memory extensive when there are many covariates due
to the multinomial-regression model used for the propensity scores.
Otherwise the function (binregATE) will run for large data.

The ATE functions need that the treatment that is given as the first
variable on the right hand side of the outcome model is a factor. The
variable is also indentified from the left hand side of the treatment
model (treat.model), that per default assumes that treatment does not
depend on any covariates.

## Average treatment effect for binary or continuous responses

In the binary case a binary outcome is specified instead of the survival
outcome, and as a consequence no-censoring adjustment is done

- the binary/numeric outcome must be a variable in the data-frame

Running the code (can also use binregATE koding cause without censorings
values, so setting cens.code=2, and time large)

``` r
kumar$cause2 <- 1*(kumar$cause==2)

b3 <- logitATE(cause2~gp.f+dnr+preauto+ttt24,kumar,treat.model=gp.f~dnr+preauto+ttt24)
summary(b3)
#>    n events
#>  400    400
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.219960  0.200830 -1.613580 -0.826341  0.0000
#> gp.f1        0.387595  0.243738 -0.090123  0.865314  0.1118
#> dnr          0.633992  0.241410  0.160837  1.107147  0.0086
#> preauto      0.139356  0.248680 -0.348049  0.626761  0.5752
#> ttt24        0.449527  0.243475 -0.027675  0.926730  0.0648
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.29524 0.19917 0.4376
#> gp.f1        1.47343 0.91382 2.3758
#> dnr          1.88512 1.17449 3.0257
#> preauto      1.14953 0.70606 1.8715
#> ttt24        1.56757 0.97270 2.5262
#> 
#> Average Treatment effects (G-formula) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.316084  0.037954  0.241695  0.390472  0.0000
#> treat1     0.400809  0.033395  0.335356  0.466262  0.0000
#> treat:1-0  0.084726  0.052651 -0.018468  0.187919  0.1076
#> 
#> Average Treatment effects (double robust) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.343451  0.042284  0.260577  0.426325  0.0000
#> treat1     0.421615  0.033616  0.355729  0.487500  0.0000
#> treat:1-0  0.078164  0.053961 -0.027599  0.183926  0.1475

###library(targeted)
###b3a <- ate(cause2~gp.f|dnr+preauto+ttt24| dnr+preauto+ttt24,kumar,family=binomial)
###summary(b3a)

## calculate also relative risk
estimate(coef=b3$riskDR,vcov=b3$var.riskDR,f=function(p) p[1]/p[2])
#>        Estimate Std.Err   2.5% 97.5%   P-value
#> treat0   0.8146  0.1194 0.5807 1.049 8.831e-12
```

Or with continuous response using normal estimating equations

``` r
b3 <- normalATE(time~gp.f+dnr+preauto+ttt24,kumar,treat.model=gp.f~dnr+preauto+ttt24)
summary(b3)
#>    n events
#>  400    400
#> 
#>  400 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)  43.4758   4.0796  35.4800  51.4716  0.0000
#> gp.f1       -27.2211   4.1379 -35.3313 -19.1109  0.0000
#> dnr           2.9485   4.4768  -5.8259  11.7228  0.5101
#> preauto      -2.8172   3.6191  -9.9105   4.2760  0.4363
#> ttt24        -1.9787   4.1285 -10.0704   6.1131  0.6318
#> 
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0     42.3170   3.8362  34.7982  49.8357       0
#> treat1     15.0958   1.3284  12.4922  17.6994       0
#> treat:1-0 -27.2211   4.1379 -35.3313 -19.1109       0
#> 
#> Average Treatment effects (double robust) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0     42.3200   4.0417  34.3985  50.2415       0
#> treat1     15.1545   1.2580  12.6888  17.6201       0
#> treat:1-0 -27.1655   4.2405 -35.4768 -18.8543       0
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
