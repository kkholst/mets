# Restricted IPCW Mean for Censored Survival Data

Provides a fast implementation of Inverse Probability of Censoring
Weighting (IPCW) regression for a single time point. It fits the model:
\$\$ E( \min(T, t) \| X ) = \exp( X^T \beta) \$\$ or, in the case of
competing risks data: \$\$ E( I(\epsilon=1) (t - \min(T, t)) \| X ) =
\exp( X^T \beta) \$\$ which represents the "Years Lost Due to Cause"
(RMTL).

## Usage

``` r
resmeanIPCW(formula, data, outcome = c("rmst", "rmtl"), ...)
```

## Arguments

- formula:

  Formula with an `Event` outcome (e.g., `Event(time, cause)`).

- data:

  Data frame containing the variables.

- outcome:

  Outcome type: `"rmst"` (Restricted Mean Survival Time) or `"rmtl"`
  (Restricted Mean Time Lost).

- ...:

  Additional arguments passed to `binreg`, such as `time`, `cause`,
  `cens.model`, `model`, `type`, etc.

## Value

An object of class `"binreg"` containing:

- coef:

  Coefficient estimates.

- se:

  Standard errors.

- var:

  Variance-covariance matrix.

- iid:

  Influence function decomposition.

- naive.var:

  Variance under known censoring model (if applicable).

- time:

  Time point used.

- outcome:

  Type of outcome analyzed.

## Details

The method solves the binomial regression IPCW response estimating
equation: \$\$ X \left( \frac{\Delta(\min(T,t)) Y}{G_c(\min(T,t))} -
\exp( X^T \beta) \right) = 0 \$\$ where \\\Delta(\min(T,t)) =
I(\min(T,t) \leq C)\\ is the indicator of being uncensored at the time
of interest.

When the status variable is binary, the outcome is assumed to be \\Y =
\min(T,t)\\ (RMST). If the status has more than two levels (competing
risks), the outcome is \\Y = (t - \min(T,t))
I(\text{status}=\text{cause})\\ (RMTL for a specific cause).

The function supports:

- **IPCW Adjustment**: Weights by the inverse of the censoring survival
  probability.

- **Augmentation**: Can include an augmentation term (type="II" or
  "III") to improve efficiency and robustness (Double Robust
  estimation).

- **Variance Estimation**: Based on the influence function, including
  adjustments for the estimation of the censoring model.

## References

Scheike, T. and Holst, K. K. (2024). Restricted mean time lost for
survival and competing risks data using mets in R. WIP.

## See also

[`binreg`](http://kkholst.github.io/mets/reference/binreg.md),
[`resmeanATE`](http://kkholst.github.io/mets/reference/resmeanATE.md),
`rmstIPCW`

## Author

Thomas Scheike

## Examples

``` r
data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001

# E( min(T;t) | X ) = exp( a+b X) with IPCW estimation 
out <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age, bmt,
                time=50, cens.model=~strata(platelet), model="exp")
summary(out)
#>    n events
#>  408    245
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  3.014368  0.065112  2.886750  3.141986  0.0000
#> tcell        0.106635  0.138266 -0.164362  0.377631  0.4406
#> platelet     0.247037  0.097339  0.056256  0.437819  0.0112
#> age         -0.185916  0.043567 -0.271305 -0.100527  0.0000
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 20.37621 17.93493 23.1498
#> tcell        1.11253  0.84844  1.4588
#> platelet     1.28023  1.05787  1.5493
#> age          0.83034  0.76238  0.9044
#> 
#> 

## Weighted GLM version (RMST)
out2 <- logitIPCW(Event(time,cause!=0)~tcell+platelet+age, bmt,
            time=50, cens.model=~strata(platelet), model="exp", outcome="rmst")
summary(out2)
#>    n   events
#>  408 7467.644
#> 
#>  408 clusters
#> coeffients:
#>               Estimate    Std.Err       2.5%      97.5% P-value
#> (Intercept)  3.0201548  0.0748731  2.8734062  3.1669033  0.0000
#> tcell        0.0248183  0.1868530 -0.3414069  0.3910434  0.8943
#> platelet     0.2531265  0.1262746  0.0056329  0.5006201  0.0450
#> age         -0.1803174  0.0526445 -0.2834988 -0.0771360  0.0006
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 20.49446 17.69720 23.7339
#> tcell        1.02513  0.71077  1.4785
#> platelet     1.28805  1.00565  1.6497
#> age          0.83501  0.75314  0.9258
#> 
#> 

### Time-lost (RMTL)
outtl <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age, bmt,
                time=50, cens.model=~strata(platelet), model="exp", outcome="rmtl")
summary(outtl)
#>    n events
#>  408    245
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  3.361796  0.047702  3.268301  3.455291  0.0000
#> tcell       -0.085197  0.127470 -0.335034  0.164641  0.5039
#> platelet    -0.239899  0.100767 -0.437399 -0.042399  0.0173
#> age          0.175625  0.044854  0.087713  0.263537  0.0001
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 28.84093 26.26666 31.6675
#> tcell        0.91833  0.71531  1.1790
#> platelet     0.78671  0.64571  0.9585
#> age          1.19199  1.09167  1.3015
#> 
#> 

### Same as Kaplan-Meier for full censoring model 
bmt$int <- with(bmt, strata(tcell, platelet))
out <- resmeanIPCW(Event(time,cause!=0)~-1+int, bmt, time=30,
                             cens.model=~strata(platelet, tcell), model="lin")
estimate(out)
#>                        Estimate Std.Err  2.5% 97.5%   P-value
#> inttcell=0, platelet=0    13.60  0.8316 11.97 15.23 3.830e-60
#> inttcell=0, platelet=1    18.90  1.2696 16.41 21.39 4.003e-50
#> inttcell=1, platelet=0    16.19  2.4061 11.48 20.91 1.704e-11
#> inttcell=1, platelet=1    17.77  2.4536 12.96 22.58 4.464e-13
out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet), data=bmt)
rm1 <- resmean_phreg(out1, times=30)
summary(rm1)
#>                     strata times    rmean  se.rmean    lower    upper
#> tcell=0, platelet=0      0    30 13.60291 0.8315434 12.06697 15.33436
#> tcell=0, platelet=1      1    30 18.90123 1.2693318 16.57016 21.56023
#> tcell=1, platelet=0      2    30 16.19125 2.4006087 12.10811 21.65132
#> tcell=1, platelet=1      3    30 17.76607 2.4422072 13.57004 23.25956
#>                     years.lost
#> tcell=0, platelet=0   16.39709
#> tcell=0, platelet=1   11.09877
#> tcell=1, platelet=0   13.80875
#> tcell=1, platelet=1   12.23393

### Years lost regression
outl <- resmeanIPCW(Event(time,cause!=0)~-1+int, bmt, time=30, outcome="years-lost",
                             cens.model=~strata(platelet, tcell), model="lin")
estimate(outl)
#>                        Estimate Std.Err   2.5% 97.5%   P-value
#> inttcell=0, platelet=0    16.40  0.8315 14.767 18.03 1.486e-86
#> inttcell=0, platelet=1    11.10  1.2693  8.611 13.59 2.254e-18
#> inttcell=1, platelet=0    13.81  2.4006  9.104 18.51 8.810e-09
#> inttcell=1, platelet=1    12.23  2.4423  7.447 17.02 5.467e-07

## Competing risks years-lost for cause 1  
out <- resmeanIPCW(Event(time,cause)~-1+int, bmt, time=30, cause=1,
                            cens.model=~strata(platelet, tcell), model="lin")
estimate(out)
#>                        Estimate Std.Err   2.5%  97.5%   P-value
#> inttcell=0, platelet=0   12.105  0.8508 10.438 13.773 6.162e-46
#> inttcell=0, platelet=1    6.884  1.1741  4.583  9.185 4.536e-09
#> inttcell=1, platelet=0    7.261  2.3533  2.648 11.873 2.033e-03
#> inttcell=1, platelet=1    5.780  2.0925  1.679  9.882 5.737e-03
## Same as integrated cumulative incidence 
rmc1 <- cif_yearslost(Event(time,cause)~strata(tcell,platelet), data=bmt, times=30)
summary(rmc1)
#> $estimate
#> $estimate$intF_1
#>                     strata times    intF_1 se.intF_1 lower_intF_1 upper_intF_1
#> tcell=0, platelet=0      0    30 12.105150 0.8508117    10.547351    13.893030
#> tcell=0, platelet=1      1    30  6.884207 1.1741044     4.928119     9.616714
#> tcell=1, platelet=0      2    30  7.260710 2.3532715     3.846767    13.704473
#> tcell=1, platelet=1      3    30  5.780407 2.0925106     2.843302    11.751515
#> 
#> $estimate$intF_2
#>                     strata times   intF_2 se.intF_2 lower_intF_2 upper_intF_2
#> tcell=0, platelet=0      0    30 4.291939 0.6161454     3.239337     5.686578
#> tcell=0, platelet=1      1    30 4.214566 0.9057057     2.765862     6.422072
#> tcell=1, platelet=0      2    30 6.548043 1.9703358     3.630623    11.809783
#> tcell=1, platelet=1      3    30 6.453525 2.0815194     3.429656    12.143489
#> 
#> 
#> $total.years.lost
#> [1] 16.39709 11.09877 13.80875 12.23393
#> 
```
