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
#> tcell.0..platelet.0        0    40 16.718643  1.162627    14.588405
#> tcell.0..platelet.1        1    40  9.731559  1.609585     7.037121
#> tcell.1..platelet.0        2    40  9.953077  3.221209     5.278067
#> tcell.1..platelet.1        3    40  8.302397  2.871799     4.214760
#> tcell.0..platelet.0.1      0    50 21.367845  1.476646    18.661116
#> tcell.0..platelet.1.1      1    50 12.985913  2.047935     9.533103
#> tcell.1..platelet.0.1      2    50 12.645385  4.089967     6.708466
#> tcell.1..platelet.1.1      3    50 11.809303  3.673690     6.418429
#>                       upper_intF_1
#> tcell.0..platelet.0       19.15994
#> tcell.0..platelet.1       13.45767
#> tcell.1..platelet.0       18.76894
#> tcell.1..platelet.1       16.35438
#> tcell.0..platelet.0.1     24.46718
#> tcell.0..platelet.1.1     17.68930
#> tcell.1..platelet.0.1     23.83641
#> tcell.1..platelet.1.1     21.72800
#> 
#> $estimate$intF_2
#>                       strata times    intF_2 se.intF_2 lower_intF_2
#> tcell.0..platelet.0        0    40  6.120930 0.8509983     4.660947
#> tcell.0..platelet.1        1    40  6.388230 1.2998312     4.287303
#> tcell.1..platelet.0        2    40 10.497757 2.8144225     6.207140
#> tcell.1..platelet.1        3    40  9.264228 2.9840677     4.927558
#> tcell.0..platelet.0.1      0    50  8.148258 1.0944538     6.262295
#> tcell.0..platelet.1.1      1    50  8.689754 1.7124193     5.905655
#> tcell.1..platelet.0.1      2    50 14.608646 3.7302732     8.856421
#> tcell.1..platelet.1.1      3    50 12.074989 3.8902005     6.421775
#>                       upper_intF_2
#> tcell.0..platelet.0       8.038235
#> tcell.0..platelet.1       9.518684
#> tcell.1..platelet.0      17.754217
#> tcell.1..platelet.1      17.417538
#> tcell.0..platelet.0.1    10.602202
#> tcell.0..platelet.1.1    12.786358
#> tcell.1..platelet.0.1    24.096928
#> tcell.1..platelet.1.1    22.704840
#> 
#> 
#> $total.years.lost
#> [1] 22.83957 16.11979 20.45083 17.56663 29.51610 21.67567 27.25403 23.88429
#> 
estimate(drm1, cause = 1)
#> [[1]]
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0   16.719   1.163 14.440 19.00 6.904e-47
#> tcell=0, platelet=1    9.732   1.610  6.577 12.89 1.485e-09
#> tcell=1, platelet=0    9.953   3.221  3.640 16.27 2.003e-03
#> tcell=1, platelet=1    8.302   2.872  2.674 13.93 3.840e-03
#> 
#> [[2]]
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    21.37   1.477 18.474 24.26 1.860e-47
#> tcell=0, platelet=1    12.99   2.048  8.972 17.00 2.283e-10
#> tcell=1, platelet=0    12.65   4.090  4.629 20.66 1.989e-03
#> tcell=1, platelet=1    11.81   3.674  4.609 19.01 1.306e-03
#> 
estimate(drm1, cause = 2)
#> [[1]]
#>                     Estimate Std.Err  2.5%  97.5%   P-value
#> tcell=0, platelet=0    6.121   0.851 4.453  7.789 6.355e-13
#> tcell=0, platelet=1    6.388   1.300 3.841  8.936 8.894e-07
#> tcell=1, platelet=0   10.498   2.814 4.982 16.014 1.915e-04
#> tcell=1, platelet=1    9.264   2.984 3.416 15.113 1.906e-03
#> 
#> [[2]]
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    8.148   1.094 6.003 10.29 9.691e-14
#> tcell=0, platelet=1    8.690   1.712 5.333 12.05 3.884e-07
#> tcell=1, platelet=0   14.609   3.730 7.297 21.92 8.994e-05
#> tcell=1, platelet=1   12.075   3.890 4.450 19.70 1.910e-03
#> 

## Comparing populations
drm1 <- cif_yearslost(Event(time, cause) ~ strata(tcell, platelet), data = bmt, times = 40)
summary(drm1, contrast = list(1:4))
#> $testintF_1
#>             Estimate Std.Err    2.5% 97.5%   P-value
#> [p1] - [p2]    6.987   1.986 3.09545 10.88 0.0004333
#> [p1] - [p3]    6.766   3.425 0.05347 13.48 0.0482026
#> [p1] - [p4]    8.416   3.098 2.34386 14.49 0.0065980
#> 
#> $testintF_2
#>             Estimate Std.Err    2.5% 97.5% P-value
#> [p1] - [p2]  -0.2673   1.554  -3.312 2.778  0.8634
#> [p1] - [p3]  -4.3768   2.940 -10.140 1.386  0.1366
#> [p1] - [p4]  -3.1433   3.103  -9.225 2.939  0.3111
#> 
#> $estimate
#> $estimate$intF_1
#>                     strata times    intF_1 se.intF_1 lower_intF_1 upper_intF_1
#> tcell=0, platelet=0      0    40 16.718643  1.162627    14.588405     19.15994
#> tcell=0, platelet=1      1    40  9.731559  1.609585     7.037121     13.45767
#> tcell=1, platelet=0      2    40  9.953077  3.221209     5.278067     18.76894
#> tcell=1, platelet=1      3    40  8.302397  2.871799     4.214760     16.35438
#> 
#> $estimate$intF_2
#>                     strata times    intF_2 se.intF_2 lower_intF_2 upper_intF_2
#> tcell=0, platelet=0      0    40  6.120930 0.8509983     4.660947     8.038235
#> tcell=0, platelet=1      1    40  6.388230 1.2998312     4.287303     9.518684
#> tcell=1, platelet=0      2    40 10.497757 2.8144225     6.207140    17.754217
#> tcell=1, platelet=1      3    40  9.264228 2.9840677     4.927558    17.417538
#> 
#> 
#> $total.years.lost
#> [1] 22.83957 16.11979 20.45083 17.56663
#> 
e1 <- estimate(drm1)
estimate(e1, rbind(c(1, -1, 0, 0)))
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    6.987   1.986 3.095 10.88 0.0004333
```
