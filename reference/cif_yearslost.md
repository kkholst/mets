# Restricted Mean Time Lost for Competing Risks

Computes the Restricted Mean Time Lost (RMTL) for competing risks based
on the integrated Aalen-Johansen estimator.

## Usage

``` r
cif_yearslost(formula, data = data, cens.code = 0, times = NULL, ...)
```

## Arguments

- formula:

  Formula for `phreg` object with `strata` to indicate strata, or `+1`
  if no strata.

- data:

  Data frame for calculations.

- cens.code:

  Censoring code (needed to separate event codes from censorings).

- times:

  Possible times for which to report restricted mean. Summary displays
  estimates for these times.

- ...:

  Additional arguments passed to `phreg`.

## Value

An object of class `"resmean_phreg"` containing:

- cumhaz:

  Matrix of cumulative hazards (years lost).

- se.cumhaz:

  Standard errors.

- intF1times:

  Years lost at specified times.

- causes:

  Vector of cause codes.

## Details

A set of time points can be given to be returned in the summary, but the
function computes years-lost for all event times and can be
plotted/viewed. The RMTL for a specific time-point can also be computed
using the `rmstIPCW` function.

## Author

Thomas Scheike

## Examples

``` r
data(bmt)
bmt$time <- bmt$time + runif(408) * 0.001

## Years lost decomposed into causes
drm1 <- cif_yearslost(Event(time, cause) ~ strata(tcell, platelet), data = bmt, times = c(40, 50))
par(mfrow = c(1, 2))
plot(drm1, cause = 1, se = 1)
plot(drm1, cause = 2, se = 1)

summary(drm1)
#> $estimate
#> $estimate$intF_1
#>                       strata times    intF_1 se.intF_1 lower_intF_1
#> tcell.0..platelet.0        0    40 16.718647  1.162628    14.588407
#> tcell.0..platelet.1        1    40  9.728016  1.609499     7.033849
#> tcell.1..platelet.0        2    40  9.953058  3.221203     5.278056
#> tcell.1..platelet.1        3    40  8.302397  2.871793     4.214767
#> tcell.0..platelet.0.1      0    50 21.367831  1.476647    18.661101
#> tcell.0..platelet.1.1      1    50 12.979253  2.047517     9.527304
#> tcell.1..platelet.0.1      2    50 12.645366  4.089961     6.708456
#> tcell.1..platelet.1.1      3    50 11.809339  3.673686     6.418465
#>                       upper_intF_1
#> tcell.0..platelet.0       19.15995
#> tcell.0..platelet.1       13.45413
#> tcell.1..platelet.0       18.76891
#> tcell.1..platelet.1       16.35436
#> tcell.0..platelet.0.1     24.46716
#> tcell.0..platelet.1.1     17.68192
#> tcell.1..platelet.0.1     23.83638
#> tcell.1..platelet.1.1     21.72801
#> 
#> $estimate$intF_2
#>                       strata times    intF_2 se.intF_2 lower_intF_2
#> tcell.0..platelet.0        0    40  6.121405 0.8509979     4.661408
#> tcell.0..platelet.1        1    40  6.388328 1.2998315     4.287395
#> tcell.1..platelet.0        2    40 10.497731 2.8144210     6.207118
#> tcell.1..platelet.1        3    40  9.264319 2.9840973     4.927606
#> tcell.0..platelet.0.1      0    50  8.149712 1.0945194     6.263607
#> tcell.0..platelet.1.1      1    50  8.690047 1.7124397     5.905904
#> tcell.1..platelet.0.1      2    50 14.608620 3.7302702     8.856401
#> tcell.1..platelet.1.1      3    50 12.075080 3.8902302     6.421823
#>                       upper_intF_2
#> tcell.0..platelet.0       8.038689
#> tcell.0..platelet.1       9.518773
#> tcell.1..platelet.0      17.754191
#> tcell.1..platelet.1      17.417710
#> tcell.0..platelet.0.1    10.603763
#> tcell.0..platelet.1.1    12.786681
#> tcell.1..platelet.0.1    24.096897
#> tcell.1..platelet.1.1    22.705012
#> 
#> 
#> $total.years.lost
#> [1] 22.84005 16.11634 20.45079 17.56672 29.51754 21.66930 27.25399 23.88442
#> 
estimate(drm1, cause = 1)
#> [[1]]
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0   16.719   1.163 14.440 19.00 6.905e-47
#> tcell=0, platelet=1    9.728   1.609  6.573 12.88 1.502e-09
#> tcell=1, platelet=0    9.953   3.221  3.640 16.27 2.003e-03
#> tcell=1, platelet=1    8.302   2.872  2.674 13.93 3.840e-03
#> 
#> [[2]]
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    21.37   1.477 18.474 24.26 1.861e-47
#> tcell=0, platelet=1    12.98   2.048  8.966 16.99 2.312e-10
#> tcell=1, platelet=0    12.65   4.090  4.629 20.66 1.989e-03
#> tcell=1, platelet=1    11.81   3.674  4.609 19.01 1.306e-03
#> 
estimate(drm1, cause = 2)
#> [[1]]
#>                     Estimate Std.Err  2.5%  97.5%   P-value
#> tcell=0, platelet=0    6.121   0.851 4.453  7.789 6.329e-13
#> tcell=0, platelet=1    6.388   1.300 3.841  8.936 8.890e-07
#> tcell=1, platelet=0   10.498   2.814 4.982 16.014 1.915e-04
#> tcell=1, platelet=1    9.264   2.984 3.416 15.113 1.906e-03
#> 
#> [[2]]
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0     8.15   1.095 6.004 10.29 9.627e-14
#> tcell=0, platelet=1     8.69   1.712 5.334 12.05 3.882e-07
#> tcell=1, platelet=0    14.61   3.730 7.297 21.92 8.994e-05
#> tcell=1, platelet=1    12.08   3.890 4.450 19.70 1.910e-03
#> 

## Comparing populations
drm1 <- cif_yearslost(Event(time, cause) ~ strata(tcell, platelet), data = bmt, times = 40)
summary(drm1, contrast = list(1:4))
#> $testintF_1
#>             Estimate Std.Err   2.5% 97.5%   P-value
#> [p1] - [p2]    6.991   1.985 3.0991 10.88 0.0004302
#> [p1] - [p3]    6.766   3.425 0.0535 13.48 0.0482015
#> [p1] - [p4]    8.416   3.098 2.3439 14.49 0.0065978
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p1] - [p2] = 0
#>   [p1] - [p3] = 0
#>   [p1] - [p4] = 0 
#>  
#> chisq = 17.643, df = 3, p-value = 0.0005211
#> 
#> $testintF_2
#>             Estimate Std.Err    2.5% 97.5% P-value
#> [p1] - [p2]  -0.2669   1.554  -3.312 2.778  0.8636
#> [p1] - [p3]  -4.3763   2.940 -10.139 1.386  0.1366
#> [p1] - [p4]  -3.1429   3.103  -9.225 2.939  0.3111
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [p1] - [p2] = 0
#>   [p1] - [p3] = 0
#>   [p1] - [p4] = 0 
#>  
#> chisq = 3.0579, df = 3, p-value = 0.3828
#> 
#> $estimate
#> $estimate$intF_1
#>                     strata times    intF_1 se.intF_1 lower_intF_1 upper_intF_1
#> tcell=0, platelet=0      0    40 16.718647  1.162628    14.588407     19.15995
#> tcell=0, platelet=1      1    40  9.728016  1.609499     7.033849     13.45413
#> tcell=1, platelet=0      2    40  9.953058  3.221203     5.278056     18.76891
#> tcell=1, platelet=1      3    40  8.302397  2.871793     4.214767     16.35436
#> 
#> $estimate$intF_2
#>                     strata times    intF_2 se.intF_2 lower_intF_2 upper_intF_2
#> tcell=0, platelet=0      0    40  6.121405 0.8509979     4.661408     8.038689
#> tcell=0, platelet=1      1    40  6.388328 1.2998315     4.287395     9.518773
#> tcell=1, platelet=0      2    40 10.497731 2.8144210     6.207118    17.754191
#> tcell=1, platelet=1      3    40  9.264319 2.9840973     4.927606    17.417710
#> 
#> 
#> $total.years.lost
#> [1] 22.84005 16.11634 20.45079 17.56672
#> 
e1 <- estimate(drm1)
estimate(e1, rbind(c(1, -1, 0, 0)))
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    6.991   1.985 3.099 10.88 0.0004302
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [tcell=0, platelet=0] - [tcell=0, platelet=1] = 0 
#>  
#> chisq = 12.3964, df = 1, p-value = 0.0004302
```
