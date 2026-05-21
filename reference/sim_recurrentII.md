# Simulate recurrent events with two event types and a terminal event

Simulates recurrent event data with up to two distinct event types and
an optional terminal event (death), based on user-supplied cumulative
hazard functions. Dependence between processes can be introduced via
shared or correlated gamma-distributed frailties.

## Usage

``` r
sim_recurrentII(
  n,
  cumhaz,
  cumhaz2,
  death.cumhaz = NULL,
  r1 = NULL,
  r2 = NULL,
  rd = NULL,
  rc = NULL,
  dependence = 0,
  var.z = 1,
  cor.mat = NULL,
  cens = NULL,
  gap.time = FALSE,
  max.recurrent = 100,
  ...
)
```

## Arguments

- n:

  Number of subjects to simulate.

- cumhaz:

  Two-column matrix `(time, cumhaz)` giving the cumulative hazard of the
  first type of recurrent event.

- cumhaz2:

  Two-column matrix `(time, cumhaz)` giving the cumulative hazard of the
  second type of recurrent event.

- death.cumhaz:

  Two-column matrix `(time, cumhaz)` giving the cumulative hazard of the
  terminal event. If `NULL`, no terminal event is simulated and
  follow-up ends at the end of `cumhaz`.

- r1:

  Optional numeric vector of length `n` with subject-specific relative
  risk multipliers for the first event type.

- r2:

  Optional numeric vector of length `n` with subject-specific relative
  risk multipliers for the second event type.

- rd:

  Optional numeric vector of length `n` with subject-specific relative
  risk multipliers for the terminal event.

- rc:

  Optional numeric vector of length `n` with subject-specific
  multipliers for the exponential censoring rate.

- dependence:

  Integer specifying the frailty structure. One of `0` (independence),
  `1` (shared gamma frailty), or `4` (shared frailty for recurrent
  events only). See Details.

- var.z:

  Variance of the gamma-distributed frailty. Default is `1`.

- cor.mat:

  Correlation matrix for the random effects. Used when `dependence = 2`
  (in `sim_recurrent_list`).

- cens:

  Rate of exponential censoring. If `NULL` (default), no additional
  censoring is applied.

- gap.time:

  Logical. If `TRUE`, event times are drawn as gap times (time since the
  last event) rather than calendar times. Default is `FALSE`.

- max.recurrent:

  Maximum number of recurrent events allowed per subject. Default is
  `100`.

- ...:

  Further arguments passed to `sim_recurrent_list`.

## Value

A data frame in counting-process format (one row per event interval per
subject) with columns:

- id:

  Subject identifier.

- start, entry:

  Interval start time.

- stop, time:

  Interval end time (event or censoring time).

- status:

  Event type at `stop`: `1` or `2` for a recurrent event of the
  corresponding type, `0` for censoring.

- death:

  Indicator for a terminal event (`1`) or censoring/survival (`0`).

Attributes `"cumhaz"`, `"death.cumhaz"`, `"rr"`, and `"rd"` store the
inputs used for simulation.

## Details

The simulation proceeds by sequentially drawing the next event time from
the specified cumulative hazards, taking the minimum of the two
recurrent event times, and stopping each subject at death or
administrative censoring.

Dependence between processes is controlled by `dependence`:

- `0`:

  Independence: all subjects have frailty fixed at 1.

- `1`:

  Shared frailty: all processes share a single gamma-distributed random
  effect with mean 1 and variance `var.z`.

- `4`:

  Recurrent-event frailty only: the two recurrent event processes share
  a gamma frailty but the terminal event is independent.

For more complex correlation structures across two event types and
death, use
[`sim_recurrentTS`](http://kkholst.github.io/mets/reference/sim_recurrentTS.md).

## See also

[`sim_recurrent`](http://kkholst.github.io/mets/reference/sim_recurrent.md),
`sim_recurrent_list`,
[`sim_recurrentTS`](http://kkholst.github.io/mets/reference/sim_recurrentTS.md)

## Author

Thomas Scheike

## Examples

``` r
data(CPH_HPN_CRBSI)
dr    <- CPH_HPN_CRBSI$terminal
base1 <- CPH_HPN_CRBSI$crbsi
base4 <- CPH_HPN_CRBSI$mechanical

## Single recurrent event type, with and without terminal event
rr <- sim_recurrent(5, base1)
dlist(rr, . ~ id, n = 0)
#> id: 1
#>        entry      time status dtime fdeath death     start      stop
#> 1     0.0000  169.3049      1  5110      0     0    0.0000  169.3049
#> 6   169.3049  481.1343      1  5110      0     0  169.3049  481.1343
#> 11  481.1343  952.4639      1  5110      0     0  481.1343  952.4639
#> 16  952.4639 1681.5204      1  5110      0     0  952.4639 1681.5204
#> 21 1681.5204 2291.2539      1  5110      0     0 1681.5204 2291.2539
#> 26 2291.2539 3236.2176      1  5110      0     0 2291.2539 3236.2176
#> 31 3236.2176 3243.0040      1  5110      0     0 3236.2176 3243.0040
#> 35 3243.0040 3712.4692      1  5110      0     0 3243.0040 3712.4692
#> 37 3712.4692 3802.3129      1  5110      0     0 3712.4692 3802.3129
#> 39 3802.3129 3877.4558      1  5110      0     0 3802.3129 3877.4558
#> 41 3877.4558 4008.8353      1  5110      0     0 3877.4558 4008.8353
#> 43 4008.8353 5110.0000      0  5110      0     0 4008.8353 5110.0000
#> ------------------------------------------------------------ 
#> id: 2
#>        entry      time status dtime fdeath death     start      stop
#> 2     0.0000  836.2402      1  5110      0     0    0.0000  836.2402
#> 7   836.2402 1582.7756      1  5110      0     0  836.2402 1582.7756
#> 12 1582.7756 1909.3792      1  5110      0     0 1582.7756 1909.3792
#> 17 1909.3792 2788.0148      1  5110      0     0 1909.3792 2788.0148
#> 22 2788.0148 4131.8937      1  5110      0     0 2788.0148 4131.8937
#> 27 4131.8937 4542.2007      1  5110      0     0 4131.8937 4542.2007
#> 32 4542.2007 5110.0000      0  5110      0     0 4542.2007 5110.0000
#> ------------------------------------------------------------ 
#> id: 3
#>        entry      time status dtime fdeath death     start      stop
#> 3     0.0000  477.0110      1  5110      0     0    0.0000  477.0110
#> 8   477.0110  684.0066      1  5110      0     0  477.0110  684.0066
#> 13  684.0066  752.7554      1  5110      0     0  684.0066  752.7554
#> 18  752.7554  786.7790      1  5110      0     0  752.7554  786.7790
#> 23  786.7790 1523.6981      1  5110      0     0  786.7790 1523.6981
#> 28 1523.6981 1697.5033      1  5110      0     0 1523.6981 1697.5033
#> 33 1697.5033 2121.4121      1  5110      0     0 1697.5033 2121.4121
#> 36 2121.4121 2211.0749      1  5110      0     0 2121.4121 2211.0749
#> 38 2211.0749 2248.7271      1  5110      0     0 2211.0749 2248.7271
#> 40 2248.7271 3559.0995      1  5110      0     0 2248.7271 3559.0995
#> 42 3559.0995 3823.9527      1  5110      0     0 3559.0995 3823.9527
#> 44 3823.9527 4697.2929      1  5110      0     0 3823.9527 4697.2929
#> 45 4697.2929 4697.6391      1  5110      0     0 4697.2929 4697.6391
#> 46 4697.6391 5031.7334      1  5110      0     0 4697.6391 5031.7334
#> 47 5031.7334 5110.0000      0  5110      0     0 5031.7334 5110.0000
#> ------------------------------------------------------------ 
#> id: 4
#>       entry     time status dtime fdeath death    start     stop
#> 4     0.000 1570.693      1  5110      0     0    0.000 1570.693
#> 9  1570.693 3057.897      1  5110      0     0 1570.693 3057.897
#> 14 3057.897 3389.234      1  5110      0     0 3057.897 3389.234
#> 19 3389.234 3913.541      1  5110      0     0 3389.234 3913.541
#> 24 3913.541 3976.522      1  5110      0     0 3913.541 3976.522
#> 29 3976.522 5110.000      0  5110      0     0 3976.522 5110.000
#> ------------------------------------------------------------ 
#> id: 5
#>       entry     time status dtime fdeath death    start     stop
#> 5     0.000 1253.563      1  5110      0     0    0.000 1253.563
#> 10 1253.563 1947.437      1  5110      0     0 1253.563 1947.437
#> 15 1947.437 2075.897      1  5110      0     0 1947.437 2075.897
#> 20 2075.897 2211.924      1  5110      0     0 2075.897 2211.924
#> 25 2211.924 3991.204      1  5110      0     0 2211.924 3991.204
#> 30 3991.204 5084.684      1  5110      0     0 3991.204 5084.684
#> 34 5084.684 5110.000      0  5110      0     0 5084.684 5110.000

rr <- sim_recurrent(5, base1, death.cumhaz = dr)
dlist(rr, . ~ id, n = 0)
#> id: 1
#>        entry      time status    dtime fdeath death     start      stop
#> 1     0.0000  350.6289      1 2996.717      1     0    0.0000  350.6289
#> 6   350.6289  362.1943      1 2996.717      1     0  350.6289  362.1943
#> 10  362.1943  884.5432      1 2996.717      1     0  362.1943  884.5432
#> 12  884.5432 1546.1629      1 2996.717      1     0  884.5432 1546.1629
#> 14 1546.1629 1701.7343      1 2996.717      1     0 1546.1629 1701.7343
#> 16 1701.7343 2996.7166      0 2996.717      1     1 1701.7343 2996.7166
#> ------------------------------------------------------------ 
#> id: 2
#>   entry     time status    dtime fdeath death start     stop
#> 2     0 209.9771      0 209.9771      1     1     0 209.9771
#> ------------------------------------------------------------ 
#> id: 3
#>      entry      time status    dtime fdeath death    start      stop
#> 3   0.0000  509.7645      1 4193.438      1     0   0.0000  509.7645
#> 7 509.7645 4193.4377      0 4193.438      1     1 509.7645 4193.4377
#> ------------------------------------------------------------ 
#> id: 4
#>        entry      time status    dtime fdeath death     start      stop
#> 4     0.0000  298.0493      1 2335.534      1     0    0.0000  298.0493
#> 8   298.0493  482.9285      1 2335.534      1     0  298.0493  482.9285
#> 11  482.9285  490.5588      1 2335.534      1     0  482.9285  490.5588
#> 13  490.5588  630.2717      1 2335.534      1     0  490.5588  630.2717
#> 15  630.2717  750.7360      1 2335.534      1     0  630.2717  750.7360
#> 17  750.7360  768.4209      1 2335.534      1     0  750.7360  768.4209
#> 18  768.4209 1510.5186      1 2335.534      1     0  768.4209 1510.5186
#> 19 1510.5186 1936.9710      1 2335.534      1     0 1510.5186 1936.9710
#> 20 1936.9710 2263.3131      1 2335.534      1     0 1936.9710 2263.3131
#> 21 2263.3131 2335.5340      0 2335.534      1     1 2263.3131 2335.5340
#> ------------------------------------------------------------ 
#> id: 5
#>      entry      time status    dtime fdeath death    start      stop
#> 5   0.0000  358.8214      1 4797.555      1     0   0.0000  358.8214
#> 9 358.8214 4797.5555      0 4797.555      1     1 358.8214 4797.5555

## Verify that estimated rates recover the true baselines (increase n for precision)
rr <- sim_recurrent(100, base1, death.cumhaz = dr)
par(mfrow = c(1, 3))
mets:::showfitsim(causes = 1, rr, dr, base1, base1)

## Shared frailty across all processes
rr <- sim_recurrent(100, base1, death.cumhaz = dr, dependence = 1, var.z = 0.4)
dtable(rr, ~death + status)
#> 
#>       status   0   1
#> death               
#> 0             26 221
#> 1             74   0

## Two event types; second type uses the mechanical complication rate
set.seed(100)
rr <- sim_recurrentII(100, base1, base4, death.cumhaz = dr)
dtable(rr, ~death + status)
#> 
#>       status   0   1   2
#> death                   
#> 0             10 295  39
#> 1             90   0   0
par(mfrow = c(2, 2))

mets:::showfitsim(causes = 2, rr, dr, base1, base4)

## Three event types and two causes of death via sim_recurrent_list
set.seed(100)
cumhaz <- list(base1, base1, base4)
drl    <- list(dr, base4)
rr     <- sim_recurrent_list(100, cumhaz, death.cumhaz = drl, dependence = 0)
dtable(rr, ~death + status)
#> 
#>       status   0   1   2   3
#> death                       
#> 0              4 232 268  33
#> 1             70   0   0   0
#> 2             26   0   0   0
mets:::showfitsimList(rr, cumhaz, drl)


```
