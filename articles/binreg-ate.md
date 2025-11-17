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
#> 122 126 152
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
#>  400    144
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.074326  0.214347 -1.494439 -0.654214  0.0000
#> gp.f1        0.669400  0.259837  0.160129  1.178670  0.0100
#> dnr          0.238437  0.263384 -0.277787  0.754661  0.3653
#> preauto      0.587886  0.275985  0.046965  1.128806  0.0332
#> ttt24        0.085280  0.275382 -0.454458  0.625018  0.7568
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.34153 0.22437 0.5199
#> gp.f1        1.95306 1.17366 3.2500
#> dnr          1.26926 0.75746 2.1269
#> preauto      1.80018 1.04809 3.0920
#> ttt24        1.08902 0.63479 1.8683
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    0.339665 0.040722 0.259851 0.419479  0.0000
#> treat1    0.496244 0.041593 0.414724 0.577764  0.0000
#> treat:1-0 0.156579 0.060890 0.037238 0.275921  0.0101
#> 
#> Average Treatment effects (double robust) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    0.311513 0.041198 0.230767 0.392260  0.0000
#> treat1    0.483995 0.043077 0.399566 0.568423  0.0000
#> treat:1-0 0.172481 0.059766 0.055343 0.289620  0.0039

b5 <- binregATE(Event(time,cause)~int+preauto+ttt24,kumar,cause=2,
        treat.model=int~preauto+ttt24,cens.code=0,time=60)
summary(b5)
#>    n events
#>  400    150
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.067627  0.285570 -1.627333 -0.507921  0.0002
#> int1.0       0.846566  0.320929  0.217557  1.475575  0.0083
#> int0.1       0.619052  0.494731 -0.350604  1.588707  0.2108
#> int1.1       0.896555  0.410239  0.092501  1.700609  0.0289
#> preauto      0.373155  0.273876 -0.163631  0.909942  0.1730
#> ttt24        0.468170  0.300260 -0.120328  1.056668  0.1189
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.34382 0.19645 0.6017
#> int1.0       2.33163 1.24304 4.3735
#> int0.1       1.85717 0.70426 4.8974
#> int1.1       2.45115 1.09691 5.4773
#> preauto      1.45231 0.84906 2.4842
#> ttt24        1.59707 0.88663 2.8768
#> 
#> Average Treatment effects (G-formula) :
#>                Estimate   Std.Err      2.5%     97.5% P-value
#> treat0.0       0.330393  0.058669  0.215404  0.445382  0.0000
#> treat1.0       0.529404  0.041929  0.447225  0.611584  0.0000
#> treat0.1       0.474113  0.101946  0.274303  0.673923  0.0000
#> treat1.1       0.541508  0.071126  0.402103  0.680913  0.0000
#> treat:1.0-0.0  0.199011  0.072251  0.057402  0.340621  0.0059
#> treat:0.1-0.0  0.143720  0.116578 -0.084770  0.372209  0.2176
#> treat:1.1-0.0  0.211115  0.096285  0.022399  0.399831  0.0283
#> treat:0.1-1.0 -0.055292  0.112083 -0.274970  0.164387  0.6218
#> treat:1.1-1.0  0.012104  0.082916 -0.150409  0.174617  0.8839
#> treat:1.1-0.1  0.067395  0.125308 -0.178205  0.312995  0.5907
#> 
#> Average Treatment effects (double robust) :
#>                Estimate   Std.Err      2.5%     97.5% P-value
#> treat0.0       0.269296  0.059958  0.151781  0.386811  0.0000
#> treat1.0       0.548493  0.054747  0.441191  0.655795  0.0000
#> treat0.1       0.505759  0.276892 -0.036940  1.048458  0.0678
#> treat1.1       0.533535  0.178811  0.183071  0.883998  0.0028
#> treat:1.0-0.0  0.279198  0.079076  0.124211  0.434184  0.0004
#> treat:0.1-0.0  0.236463  0.286973 -0.325992  0.798919  0.4099
#> treat:1.1-0.0  0.264239  0.188996 -0.106186  0.634664  0.1621
#> treat:0.1-1.0 -0.042734  0.301903 -0.634453  0.548985  0.8874
#> treat:1.1-1.0 -0.014959  0.192427 -0.392108  0.362190  0.9380
#> treat:1.1-0.1  0.027775  0.291101 -0.542772  0.598323  0.9240
```

We note that the estimates found using the large censoring model are
very different from those using the simple Kaplan-Meier weights that are
severely biased for these data. This is due to a stong censoring
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
#>  400    144
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -1.086262  0.218545 -1.514603 -0.657921  0.0000
#> gp.f1        0.627397  0.256308  0.125042  1.129753  0.0144
#> dnr          0.216831  0.267804 -0.308055  0.741717  0.4181
#> preauto      0.833926  0.298578  0.248723  1.419128  0.0052
#> ttt24       -0.050684  0.302852 -0.644263  0.542895  0.8671
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.33748 0.21990 0.5179
#> gp.f1        1.87273 1.13320 3.0949
#> dnr          1.24213 0.73488 2.0995
#> preauto      2.30234 1.28239 4.1335
#> ttt24        0.95058 0.52505 1.7210
#> 
#> Average Treatment effects (G-formula) :
#>         Estimate  Std.Err     2.5%    97.5% P-value
#> treat-1 0.499553 0.039520 0.422095 0.577011  0.0000
#> treat-0 0.354175 0.042699 0.270486 0.437864  0.0000
#> p1      0.145378 0.059607 0.028551 0.262205  0.0147
#> 
#> Average Treatment effects (double robust) :
#>         Estimate  Std.Err     2.5%    97.5% P-value
#> treat-1 0.496378 0.043208 0.411693 0.581063  0.0000
#> treat-0 0.312278 0.041300 0.231331 0.393225  0.0000
#> p1      0.184100 0.059946 0.066608 0.301593  0.0021

ib5 <- logitIPCW(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,cens.code=0,
     time=60,cens.model=~strata(gp,dnr))
summary(ib5)
#>    n events
#>  400    150
#> 
#>  400 clusters
#> coeffients:
#>               Estimate    Std.Err       2.5%      97.5% P-value
#> (Intercept) -1.0134105  0.2136059 -1.4320705 -0.5947505  0.0000
#> gp.f1        0.6820524  0.2685519  0.1557004  1.2084044  0.0111
#> dnr          0.1294593  0.2763652 -0.4122064  0.6711251  0.6395
#> preauto      0.4963073  0.3529323 -0.1954273  1.1880419  0.1597
#> ttt24        0.5577291  0.2879451 -0.0066329  1.1220912  0.0528
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.36298 0.23881 0.5517
#> gp.f1        1.97793 1.16848 3.3481
#> dnr          1.13821 0.66219 1.9564
#> preauto      1.64264 0.82248 3.2807
#> ttt24        1.74670 0.99339 3.0713

ibs <- logitIPCW(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,cens.code=0,time=60)
summary(ibs)
#>    n events
#>  400    150
#> 
#>  400 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept) -0.95461  0.21119 -1.36854 -0.54069  0.0000
#> gp.f1        0.92146  0.26818  0.39584  1.44707  0.0006
#> dnr          0.15315  0.28159 -0.39876  0.70507  0.5865
#> preauto      0.42757  0.29321 -0.14711  1.00224  0.1448
#> ttt24        0.39881  0.29279 -0.17505  0.97266  0.1732
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.38496 0.25448 0.5823
#> gp.f1        2.51296 1.48564 4.2507
#> dnr          1.16551 0.67115 2.0240
#> preauto      1.53352 0.86320 2.7244
#> ttt24        1.49005 0.83941 2.6450

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
#>                Estimate     Std.Err        2.5%       97.5% P-value
#> (Intercept) -1.07282963  0.19388190 -1.45283117 -0.69282808  0.0000
#> gp.f1        0.26229718  0.23212023 -0.19265012  0.71724448  0.2585
#> dnr          0.47162640  0.24030425  0.00063871  0.94261408  0.0497
#> preauto      0.26934293  0.24156220 -0.20411028  0.74279614  0.2648
#> ttt24        0.40066829  0.24437556 -0.07829902  0.87963559  0.1011
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.34204 0.23391 0.5002
#> gp.f1        1.29991 0.82477 2.0488
#> dnr          1.60260 1.00064 2.5667
#> preauto      1.30910 0.81537 2.1018
#> ttt24        1.49282 0.92469 2.4100
#> 
#> Average Treatment effects (G-formula) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.345704  0.037959  0.271306  0.420102  0.0000
#> treat1     0.404802  0.033367  0.339404  0.470200  0.0000
#> treat:1-0  0.059098  0.052111 -0.043037  0.161232  0.2568
#> 
#> Average Treatment effects (double robust) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     0.331726  0.039466  0.254375  0.409078  0.0000
#> treat1     0.398074  0.034407  0.330637  0.465511  0.0000
#> treat:1-0  0.066348  0.052260 -0.036080  0.168775  0.2042

###library(targeted)
###b3a <- ate(cause2~gp.f|dnr+preauto+ttt24| dnr+preauto+ttt24,kumar,family=binomial)
###summary(b3a)

## calculate also relative risk
estimate(coef=b3$riskDR,vcov=b3$var.riskDR,f=function(p) p[1]/p[2])
#>        Estimate Std.Err   2.5% 97.5%   P-value
#> treat0   0.8333  0.1223 0.5936 1.073 9.589e-12
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
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  27.80015   3.30304  21.32631  34.27399  0.0000
#> gp.f1       -16.06100   3.25515 -22.44098  -9.68102  0.0000
#> dnr          -6.51649   3.16550 -12.72075  -0.31223  0.0395
#> preauto       1.66745   2.93033  -4.07590   7.41080  0.5693
#> ttt24        -1.12151   3.27807  -7.54640   5.30338  0.7323
#> 
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0     26.2456   2.9998  20.3660  32.1251       0
#> treat1     10.1846   1.0464   8.1336  12.2355       0
#> treat:1-0 -16.0610   3.2552 -22.4410  -9.6810       0
#> 
#> Average Treatment effects (double robust) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0     28.0289   3.2124  21.7328  34.3251       0
#> treat1     10.7094   1.0007   8.7481  12.6707       0
#> treat:1-0 -17.3195   3.3652 -23.9153 -10.7238       0
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
