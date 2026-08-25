# Restricted Mean for Stratified Kaplan-Meier or Cox Model

Computes the Restricted Mean Survival Time (RMST) for stratified
Kaplan-Meier or stratified Cox models with martingale standard errors.

## Usage

``` r
resmean_phreg(x, times = NULL, covs = NULL, ...)
```

## Arguments

- x:

  Object of class `"phreg"`.

- times:

  Possible times for which to report restricted mean. If `NULL`, reports
  for all event times.

- covs:

  Possible covariates for Cox model adjustment.

- ...:

  Additional arguments passed to lower-level functions.

## Value

An object of class `"resmean_phreg"` containing:

- rmst:

  Matrix of restricted mean survival times (all event times, used by
  `plot`).

- se.rmst:

  Standard errors for RMST (all event times).

- rmst_times:

  Restricted mean (and years lost) at the specified `times`, one row per
  strata/time (`NULL` if `times=NULL`).

- estimate:

  If `times` is given: a `lava` `"estimate"` object (or, if
  `length(times)>1`, a list of such objects, one per time) with the
  across-strata estimate and covariance at each time. `NULL` if
  `times=NULL`.

## Details

The standard error is computed using linear interpolation between
standard errors at jump-times. This allows plotting the restricted mean
as a function of time.

Years lost can be computed based on this and decomposed into years lost
for different causes using the `cif_yearslost` function.

## Author

Thomas Scheike

## Examples

``` r
data(bmt)
bmt$time <- bmt$time + runif(408) * 0.001
out1 <- phreg(Surv(time, cause != 0) ~ strata(tcell, platelet), data = bmt)

## No times given: the full rmst curve at all event times (what plot() shows)
rm1 <- resmean_phreg(out1)
rm1                     ## short description; use summary() for the full curve
#> 'resmean_phreg' object
#> Strata (strata(tcell, platelet)): tcell=0, platelet=0, tcell=0, platelet=1, tcell=1, platelet=0, tcell=1, platelet=1
#> 248 distinct event times, range 0.03 - 70.625
#> ('times' was not given at call time -- use summary(x) for the full curve at all event times (also shown by plot(x)), or refit with times=... for time-specific estimates)
head(summary(rm1))
#>                strata       time       rmst      se.rmst      lower      upper
#> 1 tcell=0, platelet=0 0.03016045 0.03016045 0.000000e+00 0.03016045 0.03016045
#> 2 tcell=0, platelet=0 0.03020276 0.03020260 1.662472e-07 0.03020227 0.03020292
#> 3 tcell=0, platelet=0 0.03031414 0.03031309 7.443638e-07 0.03031163 0.03031455
#> 4 tcell=0, platelet=0 0.03054388 0.03054012 2.201127e-06 0.03053581 0.03054444
#> 5 tcell=0, platelet=0 0.03059844 0.03059383 2.572038e-06 0.03058879 0.03059887
#> 6 tcell=0, platelet=0 0.03065250 0.03064682 2.960358e-06 0.03064102 0.03065263
#>     years.lost
#> 1 0.000000e+00
#> 2 1.665755e-07
#> 3 1.043543e-06
#> 4 3.757058e-06
#> 5 4.616327e-06
#> 6 5.680482e-06

## Several times: one lava 'estimate' object per time, stored on the object
rm1 <- resmean_phreg(out1, times = 10 * (1:6))
summary(rm1)
#> time10 
#>                     Estimate Std.Err  2.5% 97.5%    P-value
#> tcell=0, platelet=0    5.863  0.2566 5.360 6.366 1.457e-115
#> tcell=0, platelet=1    7.632  0.3424 6.961 8.303 4.584e-110
#> tcell=1, platelet=0    7.278  0.7093 5.887 8.668  1.060e-24
#> tcell=1, platelet=1    7.670  0.5625 6.568 8.773  2.422e-42
#> time20 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    9.889  0.5394  8.832 10.95 4.460e-75
#> tcell=0, platelet=1   13.506  0.8000 11.938 15.07 6.051e-64
#> tcell=1, platelet=0   12.103  1.5546  9.056 15.15 6.941e-15
#> tcell=1, platelet=1   12.788  1.4676  9.911 15.66 2.949e-18
#> time30 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    13.60  0.8315 11.97 15.23 3.774e-60
#> tcell=0, platelet=1    18.90  1.2693 16.41 21.39 3.785e-50
#> tcell=1, platelet=0    16.19  2.4006 11.49 20.90 1.534e-11
#> tcell=1, platelet=1    17.77  2.4422 12.98 22.55 3.475e-13
#> time40 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    17.16   1.124 14.96 19.36 1.169e-52
#> tcell=0, platelet=1    23.88   1.737 20.48 27.29 5.267e-43
#> tcell=1, platelet=0    19.55   3.203 13.27 25.83 1.039e-09
#> tcell=1, platelet=1    22.43   3.384 15.80 29.07 3.368e-11
#> time50 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    20.48   1.411 17.72 23.25 9.551e-48
#> tcell=0, platelet=1    28.33   2.196 24.03 32.64 4.499e-38
#> tcell=1, platelet=0    22.75   4.054 14.80 30.69 2.009e-08
#> tcell=1, platelet=1    26.12   4.231 17.82 34.41 6.705e-10
#> time60 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    23.74   1.704 20.40 27.08 3.851e-44
#> tcell=0, platelet=1    32.77   2.687 27.51 38.04 3.180e-34
#> tcell=1, platelet=0    25.94   4.948 16.25 35.64 1.575e-07
#> tcell=1, platelet=1    29.67   5.160 19.56 39.78 8.904e-09
e1 <- estimate(rm1)     ## == rm1$estimate: a "resmean_estimate" list, one per time
e1
#> time10 
#>                     Estimate Std.Err  2.5% 97.5%    P-value
#> tcell=0, platelet=0    5.863  0.2566 5.360 6.366 1.457e-115
#> tcell=0, platelet=1    7.632  0.3424 6.961 8.303 4.584e-110
#> tcell=1, platelet=0    7.278  0.7093 5.887 8.668  1.060e-24
#> tcell=1, platelet=1    7.670  0.5625 6.568 8.773  2.422e-42
#> time20 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    9.889  0.5394  8.832 10.95 4.460e-75
#> tcell=0, platelet=1   13.506  0.8000 11.938 15.07 6.051e-64
#> tcell=1, platelet=0   12.103  1.5546  9.056 15.15 6.941e-15
#> tcell=1, platelet=1   12.788  1.4676  9.911 15.66 2.949e-18
#> time30 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    13.60  0.8315 11.97 15.23 3.774e-60
#> tcell=0, platelet=1    18.90  1.2693 16.41 21.39 3.785e-50
#> tcell=1, platelet=0    16.19  2.4006 11.49 20.90 1.534e-11
#> tcell=1, platelet=1    17.77  2.4422 12.98 22.55 3.475e-13
#> time40 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    17.16   1.124 14.96 19.36 1.169e-52
#> tcell=0, platelet=1    23.88   1.737 20.48 27.29 5.267e-43
#> tcell=1, platelet=0    19.55   3.203 13.27 25.83 1.039e-09
#> tcell=1, platelet=1    22.43   3.384 15.80 29.07 3.368e-11
#> time50 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    20.48   1.411 17.72 23.25 9.551e-48
#> tcell=0, platelet=1    28.33   2.196 24.03 32.64 4.499e-38
#> tcell=1, platelet=0    22.75   4.054 14.80 30.69 2.009e-08
#> tcell=1, platelet=1    26.12   4.231 17.82 34.41 6.705e-10
#> time60 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    23.74   1.704 20.40 27.08 3.851e-44
#> tcell=0, platelet=1    32.77   2.687 27.51 38.04 3.180e-34
#> tcell=1, platelet=0    25.94   4.948 16.25 35.64 1.575e-07
#> tcell=1, platelet=1    29.67   5.160 19.56 39.78 8.904e-09

## Apply a contrast to every time at once (comparing the 4 strata)
summary(e1, rbind(c(1, -1, 0, 0)))
#> time10 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err   2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....   -1.769  0.4279 -2.607 -0.93 3.572e-05
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 17.0864, df = 1, p-value = 3.572e-05
#> time20 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err   2.5%  97.5%   P-value
#> [tcell=0, platelet=0]....   -3.618  0.9649 -5.509 -1.726 0.0001774
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 14.0565, df = 1, p-value = 0.0001774
#> time30 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err   2.5%  97.5%   P-value
#> [tcell=0, platelet=0]....   -5.298   1.517 -8.272 -2.324 0.0004801
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 12.1913, df = 1, p-value = 0.0004801
#> time40 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err   2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....   -6.723   2.069 -10.78 -2.668 0.001156
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 10.5593, df = 1, p-value = 0.001156
#> time50 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err   2.5% 97.5%  P-value
#> [tcell=0, platelet=0]....   -7.847    2.61 -12.96 -2.73 0.002648
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 9.0357, df = 1, p-value = 0.002648
#> time60 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err   2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....   -9.027   3.181 -15.26 -2.792 0.004545
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 8.0522, df = 1, p-value = 0.004545

## Restrict to a single time first ...
summary(rm1, time = 50)
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    20.48   1.411 17.72 23.25 9.551e-48
#> tcell=0, platelet=1    28.33   2.196 24.03 32.64 4.499e-38
#> tcell=1, platelet=0    22.75   4.054 14.80 30.69 2.009e-08
#> tcell=1, platelet=1    26.12   4.231 17.82 34.41 6.705e-10
## ... optionally with a contrast for that time only
summary(rm1, time = 50, rbind(c(1, -1, 0, 0)))
#> Call: estimate.default(x = est, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err   2.5% 97.5%  P-value
#> [tcell=0, platelet=0]....   -7.847    2.61 -12.96 -2.73 0.002648
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 9.0357, df = 1, p-value = 0.002648
estimate(rm1, time = 50, rbind(c(1, -1, 0, 0)))
#>                           Estimate Std.Err   2.5% 97.5%  P-value
#> [tcell=0, platelet=0]....   -7.847    2.61 -12.96 -2.73 0.002648

## All pairwise differences between the 4 strata, for one time or for all times
de1 <- estimate(e1, lava:::pairwise_diff(4))
de1
#> $time10
#>                             Estimate Std.Err   2.5%    97.5%   P-value
#> [tcell=0, platelet=0]....   -1.76863  0.4279 -2.607 -0.93002 3.572e-05
#> [tcell=0, platelet=0].....1 -1.41432  0.7543 -2.893  0.06401 6.078e-02
#> [tcell=0, platelet=0].....2 -1.80681  0.6182 -3.019 -0.59511 3.472e-03
#> [tcell=0, platelet=1]....    0.35431  0.7876 -1.189  1.89796 6.528e-01
#> [tcell=0, platelet=1].....1 -0.03818  0.6585 -1.329  1.25241 9.538e-01
#> [tcell=1, platelet=0]....   -0.39249  0.9052 -2.167  1.38172 6.646e-01
#> 
#> $time20
#>                             Estimate Std.Err   2.5%   97.5%   P-value
#> [tcell=0, platelet=0]....    -3.6175  0.9649 -5.509 -1.7264 0.0001774
#> [tcell=0, platelet=0].....1  -2.2141  1.6455 -5.439  1.0109 0.1784373
#> [tcell=0, platelet=0].....2  -2.8988  1.5636 -5.963  0.1658 0.0637469
#> [tcell=0, platelet=1]....     1.4034  1.7483 -2.023  4.8301 0.4221471
#> [tcell=0, platelet=1].....1   0.7187  1.6715 -2.557  3.9948 0.6672017
#> [tcell=1, platelet=0]....    -0.6847  2.1379 -4.875  3.5055 0.7487733
#> 
#> $time30
#>                             Estimate Std.Err   2.5%   97.5%   P-value
#> [tcell=0, platelet=0]....     -5.298   1.517 -8.272 -2.3242 0.0004801
#> [tcell=0, platelet=0].....1   -2.588   2.541 -7.568  2.3911 0.3082988
#> [tcell=0, platelet=0].....2   -4.163   2.580 -9.220  0.8934 0.1065949
#> [tcell=0, platelet=1]....      2.710   2.716 -2.612  8.0324 0.3182912
#> [tcell=0, platelet=1].....1    1.135   2.752 -4.259  6.5298 0.6800154
#> [tcell=1, platelet=0]....     -1.575   3.425 -8.287  5.1371 0.6456080
#> 
#> $time40
#>                             Estimate Std.Err    2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....     -6.723   2.069 -10.778 -2.668 0.001156
#> [tcell=0, platelet=0].....1   -2.389   3.394  -9.042  4.264 0.481587
#> [tcell=0, platelet=0].....2   -5.273   3.566 -12.261  1.715 0.139181
#> [tcell=0, platelet=1]....      4.334   3.644  -2.808 11.476 0.234246
#> [tcell=0, platelet=1].....1    1.450   3.804  -6.005  8.906 0.702980
#> [tcell=1, platelet=0]....     -2.884   4.659 -12.016  6.248 0.535941
#> 
#> $time50
#>                             Estimate Std.Err    2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....     -7.847   2.610 -12.963 -2.730 0.002648
#> [tcell=0, platelet=0].....1   -2.262   4.292 -10.675  6.151 0.598168
#> [tcell=0, platelet=0].....2   -5.632   4.460 -14.373  3.109 0.206672
#> [tcell=0, platelet=1]....      5.585   4.610  -3.452 14.621 0.225775
#> [tcell=0, platelet=1].....1    2.215   4.767  -7.128 11.558 0.642149
#> [tcell=1, platelet=0]....     -3.370   5.859 -14.854  8.114 0.565240
#> 
#> $time60
#>                             Estimate Std.Err    2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....     -9.027   3.181 -15.263 -2.792 0.004545
#> [tcell=0, platelet=0].....1   -2.199   5.233 -12.455  8.057 0.674306
#> [tcell=0, platelet=0].....2   -5.928   5.434 -16.578  4.723 0.275329
#> [tcell=0, platelet=1]....      6.828   5.630  -4.206 17.863 0.225183
#> [tcell=0, platelet=1].....1    3.100   5.817  -8.302 14.502 0.594158
#> [tcell=1, platelet=0]....     -3.729   7.149 -17.740 10.282 0.601953
#> 
#> attr(,"class")
#> [1] "estimate.list" "list"         
summary(rm1, time = 50, lava:::pairwise_diff(4))
#> Call: estimate.default(x = est, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                             Estimate Std.Err    2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....     -7.847   2.610 -12.963 -2.730 0.002648
#> [tcell=0, platelet=0].....1   -2.262   4.292 -10.675  6.151 0.598168
#> [tcell=0, platelet=0].....2   -5.632   4.460 -14.373  3.109 0.206672
#> [tcell=0, platelet=1]....      5.585   4.610  -3.452 14.621 0.225775
#> [tcell=0, platelet=1].....1    2.215   4.767  -7.128 11.558 0.642149
#> [tcell=1, platelet=0]....     -3.370   5.859 -14.854  8.114 0.565240
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=1]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=1]] = 0
#>   [[tcell=1, platelet=0] - [tcell=1, platelet=1]] = 0 
#>  
#> chisq = 9.6173, df = 3, p-value = 0.02212

par(mfrow = c(1, 2))
plot(rm1, se = 1)
plot(rm1, years.lost = TRUE, se = 1)


## Comparing populations (single time -> a plain estimate object)
rm1 <- resmean_phreg(out1, times = 40)
e1 <- estimate(rm1)
estimate(e1, rbind(c(1, -1, 0, 0)))
#>                           Estimate Std.Err   2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....   -6.723   2.069 -10.78 -2.668 0.001156
de1 <- estimate(e1, lava:::pairwise_diff(4))
de1
#>                             Estimate Std.Err    2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....     -6.723   2.069 -10.778 -2.668 0.001156
#> [tcell=0, platelet=0].....1   -2.389   3.394  -9.042  4.264 0.481587
#> [tcell=0, platelet=0].....2   -5.273   3.566 -12.261  1.715 0.139181
#> [tcell=0, platelet=1]....      4.334   3.644  -2.808 11.476 0.234246
#> [tcell=0, platelet=1].....1    1.450   3.804  -6.005  8.906 0.702980
#> [tcell=1, platelet=0]....     -2.884   4.659 -12.016  6.248 0.535941
summary(rm1, lava:::pairwise_diff(4))
#> Call: estimate.default(x = est, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                             Estimate Std.Err    2.5%  97.5%  P-value
#> [tcell=0, platelet=0]....     -6.723   2.069 -10.778 -2.668 0.001156
#> [tcell=0, platelet=0].....1   -2.389   3.394  -9.042  4.264 0.481587
#> [tcell=0, platelet=0].....2   -5.273   3.566 -12.261  1.715 0.139181
#> [tcell=0, platelet=1]....      4.334   3.644  -2.808 11.476 0.234246
#> [tcell=0, platelet=1].....1    1.450   3.804  -6.005  8.906 0.702980
#> [tcell=1, platelet=0]....     -2.884   4.659 -12.016  6.248 0.535941
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=1]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=1]] = 0
#>   [[tcell=1, platelet=0] - [tcell=1, platelet=1]] = 0 
#>  
#> chisq = 11.439, df = 3, p-value = 0.009574
```
