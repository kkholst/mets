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
#> (Intercept)  3.014403  0.065110  2.886790  3.142016  0.0000
#> tcell        0.106717  0.138215 -0.164179  0.377613  0.4401
#> platelet     0.246958  0.097329  0.056196  0.437720  0.0112
#> age         -0.185997  0.043562 -0.271378 -0.100617  0.0000
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 20.37692 17.93564 23.1505
#> tcell        1.11262  0.84859  1.4588
#> platelet     1.28013  1.05780  1.5492
#> age          0.83028  0.76233  0.9043
#> 
#> 

## weighted GLM version   RMST
out2 <- logitIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
            time=50,cens.model=~strata(platelet),model="exp",outcome="rmst")
summary(out2)
#>    n   events
#>  408 7467.648
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  3.020196  0.067437  2.888023  3.152369  0.0000
#> tcell        0.024937  0.186148 -0.339907  0.389781  0.8934
#> platelet     0.253065  0.098422  0.060162  0.445967  0.0101
#> age         -0.180325  0.052152 -0.282541 -0.078108  0.0005
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 20.49531 17.95777 23.3914
#> tcell        1.02525  0.71184  1.4767
#> platelet     1.28797  1.06201  1.5620
#> age          0.83500  0.75387  0.9249
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
#> (Intercept)  3.361766  0.047705  3.268265  3.455266  0.0000
#> tcell       -0.085417  0.127458 -0.335230  0.164396  0.5028
#> platelet    -0.239824  0.100765 -0.437319 -0.042330  0.0173
#> age          0.175688  0.044858  0.087769  0.263607  0.0001
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 28.84007 26.26574 31.6667
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
#> inttcell=0, platelet=0    13.60  0.8316 11.97 15.23 3.830e-60
#> inttcell=0, platelet=1    18.90  1.2696 16.41 21.39 4.001e-50
#> inttcell=1, platelet=0    16.19  2.4061 11.48 20.91 1.705e-11
#> inttcell=1, platelet=1    17.77  2.4536 12.96 22.58 4.462e-13
out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)
rm1 <- resmean.phreg(out1,times=30)
summary(rm1)
#>                     strata times    rmean  se.rmean    lower    upper
#> tcell=0, platelet=0      0    30 13.60291 0.8315431 12.06697 15.33436
#> tcell=0, platelet=1      1    30 18.90126 1.2693300 16.57019 21.56025
#> tcell=1, platelet=0      2    30 16.19118 2.4006274 12.10801 21.65130
#> tcell=1, platelet=1      3    30 17.76614 2.4421920 13.57013 23.25959
#>                     years.lost
#> tcell=0, platelet=0   16.39709
#> tcell=0, platelet=1   11.09874
#> tcell=1, platelet=0   13.80882
#> tcell=1, platelet=1   12.23386

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
#> tcell=0, platelet=0      0    30 12.105152 4.291936 0.8508119 0.6161445
#> tcell=0, platelet=1      1    30  6.884185 4.214559 1.1741025 0.9057040
#> tcell=1, platelet=0      2    30  7.260785 6.548039 2.3532953 1.9703393
#> tcell=1, platelet=1      3    30  5.780325 6.453536 2.0924809 2.0815229
#>                     total.years.lost lower_intF_1 upper_intF_1 lower_intF_2
#> tcell=0, platelet=0         16.39709    10.547353    13.893032     3.239335
#> tcell=0, platelet=1         11.09874     4.928100     9.616688     2.765858
#> tcell=1, platelet=0         13.80882     3.846807    13.704613     3.630616
#> tcell=1, platelet=1         12.23386     2.843261    11.751348     3.429661
#>                     upper_intF_2
#> tcell=0, platelet=0     5.686572
#> tcell=0, platelet=1     6.422060
#> tcell=1, platelet=0    11.809793
#> tcell=1, platelet=1    12.143509
#> 
```
