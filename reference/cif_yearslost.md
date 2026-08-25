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

  Matrix of years lost (all event times, used by `plot`), one column per
  cause.

- se.cumhaz:

  Standard errors (all event times).

- yearslost_times:

  Years lost per cause at the specified `times`, one row per strata/time
  (`NULL` if `times=NULL`).

- estimate:

  If `times` is given: a named list, one element per cause, each a
  `lava` `"estimate"` object (or, if `length(times)>1`, a list of such
  objects, one per time). Causes are estimated separately as their
  covariance is not available. `NULL` if `times=NULL`.

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

## No times given: full years-lost curve at all event times, one table per cause
drm1 <- cif_yearslost(Event(time, cause) ~ strata(tcell, platelet), data = bmt)
drm1                         ## short description; use summary() for the full curves
#> 'resmean_phreg' object
#> Competing risks, causes: 1, 2
#> Strata (strata(tcell, platelet)): tcell=0, platelet=0, tcell=0, platelet=1, tcell=1, platelet=0, tcell=1, platelet=1
#> 248 distinct event times, range 0.03 - 70.626
#> ('times' was not given at call time -- use summary(x) for the full curve at all event times (also shown by plot(x)), or refit with times=... for time-specific estimates)
tail(summary(drm1, cause=1))  ## just cause 1
#>                  strata     time years.lost se.years.lost     lower    upper
#> 243 tcell=0, platelet=0 40.13293  16.780445      1.166778 14.642584 19.23044
#> 244 tcell=0, platelet=1 40.26369   9.811029      1.621059  7.096955 13.56304
#> 245 tcell=1, platelet=1 41.77681   8.821595      3.013795  4.515897 17.23258
#> 246 tcell=0, platelet=0 52.27050  22.423447      1.548327 19.585171 25.67305
#> 247 tcell=0, platelet=0 68.28972  29.871106      2.056342 26.100821 34.18601
#> 248 tcell=0, platelet=0 70.62570  30.957149      2.130654 27.050559 35.42792
plot(drm1, se=1)


## Years lost decomposed into causes, at specific times
drm1 <- cif_yearslost(Event(time, cause) ~ strata(tcell, platelet), data = bmt, times = c(40, 50))
par(mfrow = c(1, 2))
plot(drm1, cause = 1, se = 1)
plot(drm1, cause = 2, se = 1)

summary(drm1)          ## both causes
#> cause1 
#> time40 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0   16.719   1.163 14.440 19.00 6.904e-47
#> tcell=0, platelet=1    9.732   1.610  6.577 12.89 1.485e-09
#> tcell=1, platelet=0    9.953   3.221  3.640 16.27 2.003e-03
#> tcell=1, platelet=1    8.302   2.872  2.674 13.93 3.840e-03
#> time50 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    21.37   1.477 18.474 24.26 1.860e-47
#> tcell=0, platelet=1    12.99   2.048  8.972 17.00 2.283e-10
#> tcell=1, platelet=0    12.65   4.090  4.629 20.66 1.989e-03
#> tcell=1, platelet=1    11.81   3.674  4.609 19.01 1.306e-03
#> cause2 
#> time40 
#>                     Estimate Std.Err  2.5%  97.5%   P-value
#> tcell=0, platelet=0    6.121   0.851 4.453  7.789 6.355e-13
#> tcell=0, platelet=1    6.388   1.300 3.841  8.936 8.894e-07
#> tcell=1, platelet=0   10.498   2.814 4.982 16.014 1.915e-04
#> tcell=1, platelet=1    9.264   2.984 3.416 15.113 1.906e-03
#> time50 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    8.148   1.094 6.003 10.29 9.691e-14
#> tcell=0, platelet=1    8.690   1.712 5.333 12.05 3.884e-07
#> tcell=1, platelet=0   14.609   3.730 7.297 21.92 8.994e-05
#> tcell=1, platelet=1   12.075   3.890 4.450 19.70 1.910e-03
summary(drm1, cause=1) ## just cause 1
#> time40 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0   16.719   1.163 14.440 19.00 6.904e-47
#> tcell=0, platelet=1    9.732   1.610  6.577 12.89 1.485e-09
#> tcell=1, platelet=0    9.953   3.221  3.640 16.27 2.003e-03
#> tcell=1, platelet=1    8.302   2.872  2.674 13.93 3.840e-03
#> time50 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    21.37   1.477 18.474 24.26 1.860e-47
#> tcell=0, platelet=1    12.99   2.048  8.972 17.00 2.283e-10
#> tcell=1, platelet=0    12.65   4.090  4.629 20.66 1.989e-03
#> tcell=1, platelet=1    11.81   3.674  4.609 19.01 1.306e-03

## Causes are stored (and estimated) separately: one "resmean_estimate" list per cause
e1 <- estimate(drm1, cause = 1)
e1
#> time40 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0   16.719   1.163 14.440 19.00 6.904e-47
#> tcell=0, platelet=1    9.732   1.610  6.577 12.89 1.485e-09
#> tcell=1, platelet=0    9.953   3.221  3.640 16.27 2.003e-03
#> tcell=1, platelet=1    8.302   2.872  2.674 13.93 3.840e-03
#> time50 
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    21.37   1.477 18.474 24.26 1.860e-47
#> tcell=0, platelet=1    12.99   2.048  8.972 17.00 2.283e-10
#> tcell=1, platelet=0    12.65   4.090  4.629 20.66 1.989e-03
#> tcell=1, platelet=1    11.81   3.674  4.609 19.01 1.306e-03
e2 <- estimate(drm1, cause = 2)
e2
#> time40 
#>                     Estimate Std.Err  2.5%  97.5%   P-value
#> tcell=0, platelet=0    6.121   0.851 4.453  7.789 6.355e-13
#> tcell=0, platelet=1    6.388   1.300 3.841  8.936 8.894e-07
#> tcell=1, platelet=0   10.498   2.814 4.982 16.014 1.915e-04
#> tcell=1, platelet=1    9.264   2.984 3.416 15.113 1.906e-03
#> time50 
#>                     Estimate Std.Err  2.5% 97.5%   P-value
#> tcell=0, platelet=0    8.148   1.094 6.003 10.29 9.691e-14
#> tcell=0, platelet=1    8.690   1.712 5.333 12.05 3.884e-07
#> tcell=1, platelet=0   14.609   3.730 7.297 21.92 8.994e-05
#> tcell=1, platelet=1   12.075   3.890 4.450 19.70 1.910e-03

## Apply a contrast to every time at once, for one cause
summary(e1, rbind(c(1, -1, 0, 0)))
#> time40 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    6.987   1.986 3.095 10.88 0.0004333
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 12.3829, df = 1, p-value = 0.0004333
#> time50 
#> Call: estimate.default(x = o, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    8.382   2.525 3.433 13.33 0.0009006
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 11.0215, df = 1, p-value = 0.0009006

## Restrict to a single time first ...
summary(drm1, cause = 1, time = 50)
#>                     Estimate Std.Err   2.5% 97.5%   P-value
#> tcell=0, platelet=0    21.37   1.477 18.474 24.26 1.860e-47
#> tcell=0, platelet=1    12.99   2.048  8.972 17.00 2.283e-10
#> tcell=1, platelet=0    12.65   4.090  4.629 20.66 1.989e-03
#> tcell=1, platelet=1    11.81   3.674  4.609 19.01 1.306e-03
## ... optionally with a contrast for that time only
summary(drm1, cause = 1, time = 50, rbind(c(1, -1, 0, 0)))
#> Call: estimate.default(x = est, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    8.382   2.525 3.433 13.33 0.0009006
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 11.0215, df = 1, p-value = 0.0009006
estimate(drm1, cause = 1, time = 50, rbind(c(1, -1, 0, 0)))
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    8.382   2.525 3.433 13.33 0.0009006

## All pairwise differences between the 4 strata, for one time or for all times
de1 <- estimate(e1, lava:::pairwise_diff(4))
de1
#> $time40
#>                             Estimate Std.Err     2.5%  97.5%   P-value
#> [tcell=0, platelet=0]....     6.9871   1.986  3.09545 10.879 0.0004333
#> [tcell=0, platelet=0].....1   6.7656   3.425  0.05347 13.478 0.0482026
#> [tcell=0, platelet=0].....2   8.4162   3.098  2.34386 14.489 0.0065980
#> [tcell=0, platelet=1]....    -0.2215   3.601 -7.27928  6.836 0.9509481
#> [tcell=0, platelet=1].....1   1.4292   3.292 -5.02326  7.882 0.6642032
#> [tcell=1, platelet=0]....     1.6507   4.315 -6.80752 10.109 0.7020894
#> 
#> $time50
#>                             Estimate Std.Err    2.5%  97.5%   P-value
#> [tcell=0, platelet=0]....     8.3819   2.525  3.4335 13.330 0.0009006
#> [tcell=0, platelet=0].....1   8.7225   4.348  0.1998 17.245 0.0448653
#> [tcell=0, platelet=0].....2   9.5585   3.959  1.7984 17.319 0.0157712
#> [tcell=0, platelet=1]....     0.3405   4.574 -8.6244  9.305 0.9406540
#> [tcell=0, platelet=1].....1   1.1766   4.206 -7.0669  9.420 0.7796702
#> [tcell=1, platelet=0]....     0.8361   5.498 -9.9391 11.611 0.8791232
#> 
#> attr(,"class")
#> [1] "estimate.list" "list"         
summary(drm1, cause = 1, time = 50, lava:::pairwise_diff(4))
#> Call: estimate.default(x = est, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                             Estimate Std.Err    2.5%  97.5%   P-value
#> [tcell=0, platelet=0]....     8.3819   2.525  3.4335 13.330 0.0009006
#> [tcell=0, platelet=0].....1   8.7225   4.348  0.1998 17.245 0.0448653
#> [tcell=0, platelet=0].....2   9.5585   3.959  1.7984 17.319 0.0157712
#> [tcell=0, platelet=1]....     0.3405   4.574 -8.6244  9.305 0.9406540
#> [tcell=0, platelet=1].....1   1.1766   4.206 -7.0669  9.420 0.7796702
#> [tcell=1, platelet=0]....     0.8361   5.498 -9.9391 11.611 0.8791232
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=1]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=1]] = 0
#>   [[tcell=1, platelet=0] - [tcell=1, platelet=1]] = 0 
#>  
#> chisq = 15.5277, df = 3, p-value = 0.001417

## Comparing populations (single time -> a plain estimate object per cause)
drm1 <- cif_yearslost(Event(time, cause) ~ strata(tcell, platelet), data = bmt, times = 40)
summary(drm1, cause = 1, contrast = rbind(c(1, -1, 0, 0)))
#> Call: estimate.default(x = est, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    6.987   1.986 3.095 10.88 0.0004333
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0 
#>  
#> chisq = 12.3829, df = 1, p-value = 0.0004333
e1 <- estimate(drm1, cause = 1)
estimate(e1, rbind(c(1, -1, 0, 0)))
#>                           Estimate Std.Err  2.5% 97.5%   P-value
#> [tcell=0, platelet=0]....    6.987   1.986 3.095 10.88 0.0004333
de1 <- estimate(e1, lava:::pairwise_diff(4))
de1
#>                             Estimate Std.Err     2.5%  97.5%   P-value
#> [tcell=0, platelet=0]....     6.9871   1.986  3.09545 10.879 0.0004333
#> [tcell=0, platelet=0].....1   6.7656   3.425  0.05347 13.478 0.0482026
#> [tcell=0, platelet=0].....2   8.4162   3.098  2.34386 14.489 0.0065980
#> [tcell=0, platelet=1]....    -0.2215   3.601 -7.27928  6.836 0.9509481
#> [tcell=0, platelet=1].....1   1.4292   3.292 -5.02326  7.882 0.6642032
#> [tcell=1, platelet=0]....     1.6507   4.315 -6.80752 10.109 0.7020894
summary(drm1, cause = 1, contrast = lava:::pairwise_diff(4))
#> Call: estimate.default(x = est, f = contrast)
#> ────────────────────────────────────────────────────────────
#>                             Estimate Std.Err     2.5%  97.5%   P-value
#> [tcell=0, platelet=0]....     6.9871   1.986  3.09545 10.879 0.0004333
#> [tcell=0, platelet=0].....1   6.7656   3.425  0.05347 13.478 0.0482026
#> [tcell=0, platelet=0].....2   8.4162   3.098  2.34386 14.489 0.0065980
#> [tcell=0, platelet=1]....    -0.2215   3.601 -7.27928  6.836 0.9509481
#> [tcell=0, platelet=1].....1   1.4292   3.292 -5.02326  7.882 0.6642032
#> [tcell=1, platelet=0]....     1.6507   4.315 -6.80752 10.109 0.7020894
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [[tcell=0, platelet=0] - [tcell=0, platelet=1]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=0] - [tcell=1, platelet=1]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=0]] = 0
#>   [[tcell=0, platelet=1] - [tcell=1, platelet=1]] = 0
#>   [[tcell=1, platelet=0] - [tcell=1, platelet=1]] = 0 
#>  
#> chisq = 17.6322, df = 3, p-value = 0.0005238
```
