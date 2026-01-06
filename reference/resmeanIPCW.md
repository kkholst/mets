# Restricted IPCW mean for censored survival data

Simple and fast version for IPCW regression for just one time-point thus
fitting the model \$\$E( min(T, t) \| X ) = exp( X^T beta) \$\$ or in
the case of competing risks data \$\$E( I(epsilon=1) (t - min(T ,t)) \|
X ) = exp( X^T beta) \$\$ thus given years lost to cause, see `binreg`
for the arguments.

## Usage

``` r
resmeanIPCW(formula, data, outcome = c("rmst", "rmtl"), ...)
```

## Arguments

- formula:

  formula with outcome on Event form

- data:

  data frame

- outcome:

  can do either rmst regression ('rmst') or years-lost regression
  ('rmtl')

- ...:

  Additional arguments to binreg

## Details

When the status is binary assumes it is a survival setting and default
is to consider outcome Y=min(T,t), if status has more than two levels,
then computes years lost due to the specified cause, thus using the
response \$\$ Y = (t-min(T,t)) I(status=cause) \$\$

Based on binomial regresion IPCW response estimating equation: \$\$ X (
\Delta(min(T,t)) Y /G_c(min(T,t)) - exp( X^T beta)) = 0 \$\$ for IPCW
adjusted responses. Here \$\$ \Delta(min(T,t)) = I ( min(T ,t) \leq C )
\$\$ is indicator of being uncensored. Concretely, the uncensored
observations at time t will count those with an event (of any type)
before t and those with a censoring time at t or further out. One should
therefore be a bit careful when data has been constructed such that some
of the event times T are equivalent to t.

Can also solve the binomial regresion IPCW response estimating equation:
\$\$ h(X) X ( \Delta(min(T,t)) Y /G_c(min(T,t)) - exp( X^T beta)) = 0
\$\$ for IPCW adjusted responses where \$h\$ is given as an argument
together with iid of censoring with h.

By using appropriately the h argument we can also do the efficient IPCW
estimator estimator.

Variance is based on \$\$ \sum w_i^2 \$\$ also with IPCW adjustment, and
naive.var is variance under known censoring model.

When Ydirect is given it solves : \$\$ X ( \Delta(min(T,t)) Ydirect
/G_c(min(T,t)) - exp( X^T beta)) = 0 \$\$ for IPCW adjusted responses.

The actual influence (type="II") function is based on augmenting with
\$\$ X \int_0^t E(Y \| T\>s) /G_c(s) dM_c(s) \$\$ and alternatively just
solved directly (type="I") without any additional terms.

Censoring model may depend on strata.

## References

Scheike, T. and Holst, K. K. Restricted mean time lost for survival and
competing risks data using mets in R, WIP

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001
# E( min(T;t) | X ) = exp( a+b X) with IPCW estimation 
out <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
                time=50,cens.model=~strata(platelet),model="exp")
summary(out)
#>    n events
#>  408    245
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  3.014393  0.065110  2.886780  3.142007  0.0000
#> tcell        0.106872  0.138216 -0.164027  0.377771  0.4394
#> platelet     0.246943  0.097333  0.056173  0.437712  0.0112
#> age         -0.185971  0.043566 -0.271359 -0.100583  0.0000
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 20.37672 17.93546 23.1503
#> tcell        1.11279  0.84872  1.4590
#> platelet     1.28011  1.05778  1.5492
#> age          0.83030  0.76234  0.9043
#> 
#> 

## weighted GLM version   RMST
out2 <- logitIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
            time=50,cens.model=~strata(platelet),model="exp",outcome="rmst")
summary(out2)
#>    n   events
#>  408 7467.653
#> 
#>  408 clusters
#> coeffients:
#>               Estimate    Std.Err       2.5%      97.5% P-value
#> (Intercept)  3.0201981  0.0748709  2.8734539  3.1669423  0.0000
#> tcell        0.0249363  0.1868557 -0.3412942  0.3911667  0.8938
#> platelet     0.2530618  0.1262735  0.0055703  0.5005532  0.0451
#> age         -0.1803245  0.0526438 -0.2835045 -0.0771445  0.0006
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 20.49535 17.69804 23.7348
#> tcell        1.02525  0.71085  1.4787
#> platelet     1.28796  1.00559  1.6496
#> age          0.83500  0.75314  0.9258
#> 
#> 

### time-lost
outtl <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
                time=50,cens.model=~strata(platelet),model="exp",outcome="rmtl")
summary(outtl)
#>    n events
#>  408    245
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  3.361764  0.047705  3.268264  3.455265  0.0000
#> tcell       -0.085417  0.127458 -0.335230  0.164396  0.5028
#> platelet    -0.239822  0.100764 -0.437317 -0.042327  0.0173
#> age          0.175688  0.044858  0.087769  0.263607  0.0001
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 28.84003 26.26571 31.6667
#> tcell        0.91813  0.71517  1.1787
#> platelet     0.78677  0.64577  0.9586
#> age          1.19207  1.09174  1.3016
#> 
#> 

### same as Kaplan-Meier for full censoring model 
bmt$int <- with(bmt,strata(tcell,platelet))
out <- resmeanIPCW(Event(time,cause!=0)~-1+int,bmt,time=30,
                             cens.model=~strata(platelet,tcell),model="lin")
estimate(out)
#>                        Estimate Std.Err  2.5% 97.5%   P-value
#> inttcell=0, platelet=0    13.60  0.8316 11.97 15.23 3.824e-60
#> inttcell=0, platelet=1    18.90  1.2696 16.41 21.39 4.002e-50
#> inttcell=1, platelet=0    16.19  2.4061 11.48 20.91 1.705e-11
#> inttcell=1, platelet=1    17.77  2.4536 12.96 22.58 4.464e-13
out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)
rm1 <- resmean.phreg(out1,times=30)
summary(rm1)
#>                     strata times    rmean  se.rmean    lower    upper
#> tcell=0, platelet=0      0    30 13.60295 0.8315412 12.06701 15.33440
#> tcell=0, platelet=1      1    30 18.90126 1.2693314 16.57019 21.56026
#> tcell=1, platelet=0      2    30 16.19122 2.4006180 12.10806 21.65132
#> tcell=1, platelet=1      3    30 17.76607 2.4422046 13.57005 23.25956
#>                     years.lost
#> tcell=0, platelet=0   16.39705
#> tcell=0, platelet=1   11.09874
#> tcell=1, platelet=0   13.80878
#> tcell=1, platelet=1   12.23393

### years lost regression
outl <- resmeanIPCW(Event(time,cause!=0)~-1+int,bmt,time=30,outcome="years-lost",
                             cens.model=~strata(platelet,tcell),model="lin")
estimate(outl)
#>                        Estimate Std.Err   2.5% 97.5%   P-value
#> inttcell=0, platelet=0    16.40  0.8315 14.767 18.03 1.485e-86
#> inttcell=0, platelet=1    11.10  1.2693  8.611 13.59 2.255e-18
#> inttcell=1, platelet=0    13.81  2.4006  9.104 18.51 8.811e-09
#> inttcell=1, platelet=1    12.23  2.4423  7.447 17.02 5.467e-07

## competing risks years-lost for cause 1  
out <- resmeanIPCW(Event(time,cause)~-1+int,bmt,time=30,cause=1,
                            cens.model=~strata(platelet,tcell),model="lin")
estimate(out)
#>                        Estimate Std.Err   2.5%  97.5%   P-value
#> inttcell=0, platelet=0   12.105  0.8508 10.438 13.773 6.162e-46
#> inttcell=0, platelet=1    6.884  1.1741  4.583  9.185 4.536e-09
#> inttcell=1, platelet=0    7.261  2.3533  2.648 11.873 2.033e-03
#> inttcell=1, platelet=1    5.780  2.0925  1.679  9.882 5.737e-03
## same as integrated cumulative incidence 
rmc1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=30)
summary(rmc1)
#> $estimate
#>                     strata times    intF_1   intF_2 se.intF_1 se.intF_2
#> tcell=0, platelet=0      0    30 12.105108 4.291939 0.8508088 0.6161454
#> tcell=0, platelet=1      1    30  6.884191 4.214551 1.1741028 0.9057053
#> tcell=1, platelet=0      2    30  7.260755 6.548026 2.3532855 1.9703338
#> tcell=1, platelet=1      3    30  5.780372 6.453554 2.0924973 2.0815288
#>                     total.years.lost lower_intF_1 upper_intF_1 lower_intF_2
#> tcell=0, platelet=0         16.39705    10.547314    13.892982     3.239337
#> tcell=0, platelet=1         11.09874     4.928105     9.616694     2.765849
#> tcell=1, platelet=0         13.80878     3.846791    13.704556     3.630610
#> tcell=1, platelet=1         12.23393     2.843285    11.751441     3.429671
#>                     upper_intF_2
#> tcell=0, platelet=0     5.686578
#> tcell=0, platelet=1     6.422058
#> tcell=1, platelet=0    11.809764
#> tcell=1, platelet=1    12.143543
#> 
```
