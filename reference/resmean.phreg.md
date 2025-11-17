# Restricted mean for stratified Kaplan-Meier or Cox model with martingale standard errors

Restricted mean for stratified Kaplan-Meier or stratified Cox with
martingale standard error. Standard error is computed using linear
interpolation between standard errors at jump-times. Plots gives
restricted mean at all times. Years lost can be computed based on this
and decomposed into years lost for different causes using the
cif.yearslost function that is based on integrating the cumulative
incidence functions. One particular feature of these functions are that
the restricted mean and years-lost are computed for all event times as
functions and can be plotted/viewed. When times are given and beyond the
last event time withn a strata the curves are extrapolated using the
estimates of cumulative incidence.

## Usage

``` r
resmean.phreg(x, times = NULL, covs = NULL, ...)
```

## Arguments

- x:

  phreg object

- times:

  possible times for which to report restricted mean

- covs:

  possible covariate for Cox model

- ...:

  Additional arguments to lower level funtions

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt); bmt$time <- bmt$time+runif(408)*0.001
out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)

rm1 <- resmean.phreg(out1,times=10*(1:6))
summary(rm1)
#>                       strata times     rmean  se.rmean     lower     upper
#> tcell.0..platelet.0        0    10  5.863348 0.2565969  5.381392  6.388468
#> tcell.0..platelet.1        1    10  7.631932 0.3423878  6.989521  8.333387
#> tcell.1..platelet.0        2    10  7.277571 0.7092918  6.012098  8.809410
#> tcell.1..platelet.1        3    10  7.670168 0.5624467  6.643348  8.855697
#> tcell.0..platelet.0.1      0    20  9.888956 0.5393851  8.886328 11.004709
#> tcell.0..platelet.1.1      1    20 13.506404 0.8000257 12.025979 15.169072
#> tcell.1..platelet.0.1      2    20 12.102998 1.5545668  9.409371 15.567732
#> tcell.1..platelet.1.1      3    20 12.787751 1.4675791 10.211899 16.013337
#> tcell.0..platelet.0.2      0    30 13.602943 0.8315414 12.067001 15.334387
#> tcell.0..platelet.1.2      1    30 18.901243 1.2693309 16.570176 21.560241
#> tcell.1..platelet.0.2      2    30 16.191191 2.4006237 12.108029 21.651307
#> tcell.1..platelet.1.2      3    30 17.766106 2.4421997 13.570089 23.259577
#> tcell.0..platelet.0.3      0    40 17.159949 1.1235958 15.093195 19.509709
#> tcell.0..platelet.1.3      1    40 23.883651 1.7373041 20.710189 27.543389
#> tcell.1..platelet.0.3      2    40 19.549224 3.2030928 14.179542 26.952361
#> tcell.1..platelet.1.3      3    40 22.433334 3.3838371 16.691634 30.150103
#> tcell.0..platelet.0.4      0    50 20.482459 1.4110533 17.895429 23.443479
#> tcell.0..platelet.1.4      1    50 28.330696 2.1961763 24.337316 32.979329
#> tcell.1..platelet.0.4      2    50 22.746027 4.0537082 16.040124 32.255470
#> tcell.1..platelet.1.4      3    50 26.115671 4.2306850 19.011146 35.875180
#> tcell.0..platelet.0.5      0    60 23.741473 1.7038101 20.626288 27.327145
#> tcell.0..platelet.1.5      1    60 32.771223 2.6865772 27.906884 38.483446
#> tcell.1..platelet.0.5      2    60 25.942830 4.9476022 17.851850 37.700879
#> tcell.1..platelet.1.5      3    60 29.671639 5.1599146 21.101684 41.722081
#>                       years.lost
#> tcell.0..platelet.0     4.136652
#> tcell.0..platelet.1     2.368068
#> tcell.1..platelet.0     2.722429
#> tcell.1..platelet.1     2.329832
#> tcell.0..platelet.0.1  10.111044
#> tcell.0..platelet.1.1   6.493596
#> tcell.1..platelet.0.1   7.897002
#> tcell.1..platelet.1.1   7.212249
#> tcell.0..platelet.0.2  16.397057
#> tcell.0..platelet.1.2  11.098757
#> tcell.1..platelet.0.2  13.808809
#> tcell.1..platelet.1.2  12.233894
#> tcell.0..platelet.0.3  22.840051
#> tcell.0..platelet.1.3  16.116349
#> tcell.1..platelet.0.3  20.450776
#> tcell.1..platelet.1.3  17.566666
#> tcell.0..platelet.0.4  29.517541
#> tcell.0..platelet.1.4  21.669304
#> tcell.1..platelet.0.4  27.253973
#> tcell.1..platelet.1.4  23.884329
#> tcell.0..platelet.0.5  36.258527
#> tcell.0..platelet.1.5  27.228777
#> tcell.1..platelet.0.5  34.057170
#> tcell.1..platelet.1.5  30.328361
par(mfrow=c(1,2))
plot(rm1,se=1)
plot(rm1,years.lost=TRUE,se=1)


## comparing populations, can also be done using rmstIPCW via influence functions
rm1 <- resmean.phreg(out1,times=40)
e1 <- estimate(rm1)
e1
#>    Estimate Std.Err  2.5% 97.5%   P-value
#> p1    17.16   1.124 14.96 19.36 1.169e-52
#> p2    23.88   1.737 20.48 27.29 5.269e-43
#> p3    19.55   3.203 13.27 25.83 1.039e-09
#> p4    22.43   3.384 15.80 29.07 3.367e-11
estimate(e1,rbind(c(1,-1,0,0)))
#>             Estimate Std.Err   2.5%  97.5%  P-value
#> [p1] - [p2]   -6.724   2.069 -10.78 -2.669 0.001155
#> 
#>  Null Hypothesis: 
#>   [p1] - [p2] = 0 

## years.lost decomposed into causes
drm1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=10*(1:6))
par(mfrow=c(1,2)); plot(drm1,cause=1,se=1); plot(drm1,cause=2,se=1);

summary(drm1)
#> $estimate
#>                       strata times    intF_1     intF_2 se.intF_1 se.intF_2
#> tcell.0..platelet.0        0    10  3.117716  1.0189361 0.2487175 0.1703642
#> tcell.0..platelet.1        1    10  1.710926  0.6571419 0.3238403 0.1870597
#> tcell.1..platelet.0        2    10  1.876136  0.8462935 0.6339068 0.4726441
#> tcell.1..platelet.1        3    10  1.358597  0.9712345 0.5303058 0.3617336
#> tcell.0..platelet.0.1      0    20  7.517560  2.5934839 0.5441129 0.3861266
#> tcell.0..platelet.1.1      1    20  4.230958  2.2626386 0.7414091 0.5327621
#> tcell.1..platelet.0.1      2    20  4.568444  3.3285585 1.4876774 1.1718411
#> tcell.1..platelet.1.1      3    20  3.569482  3.6427672 1.3002981 1.1906629
#> tcell.0..platelet.0.2      0    30 12.105122  4.2919342 0.8508098 0.6161442
#> tcell.0..platelet.1.2      1    30  6.884187  4.2145694 1.1741026 0.9057061
#> tcell.1..platelet.0.2      2    30  7.260751  6.5480579 2.3532866 1.9703434
#> tcell.1..platelet.1.2      3    30  5.780366  6.4535279 2.0924951 2.0815219
#> tcell.0..platelet.0.3      0    40 16.718632  6.1214190 1.1626270 0.8509985
#> tcell.0..platelet.1.3      1    40  9.728009  6.3883401 1.6094993 1.2998358
#> tcell.1..platelet.0.3      2    40  9.953059 10.4977173 3.2212046 2.8144104
#> tcell.1..platelet.1.3      3    40  8.302378  9.2642886 2.8717876 2.9840876
#> tcell.0..platelet.0.4      0    50 21.367817  8.1497247 1.4766454 1.0945204
#> tcell.0..platelet.1.4      1    50 12.979245  8.6900589 2.0475166 1.7124434
#> tcell.1..platelet.0.4      2    50 12.645367 14.6086064 4.0899618 3.7302619
#> tcell.1..platelet.1.4      3    50 11.809280 12.0750492 3.6736800 3.8902204
#> tcell.0..platelet.0.5      0    60 26.017001 10.2415261 1.7930399 1.3471514
#> tcell.0..platelet.1.5      1    60 16.237000 10.9917777 2.5097415 2.1389922
#> tcell.1..platelet.0.5      2    60 15.337674 18.7194955 4.9591172 4.6873852
#> tcell.1..platelet.1.5      3    60 15.442551 14.8858099 4.5899546 4.7978997
#>                       total.years.lost lower_intF_1 upper_intF_1 lower_intF_2
#> tcell.0..platelet.0           4.136652    2.6664378     3.645370    0.7342225
#> tcell.0..platelet.1           2.368068    1.1806408     2.479390    0.3761484
#> tcell.1..platelet.0           2.722429    0.9675227     3.638039    0.2832278
#> tcell.1..platelet.1           2.329832    0.6321767     2.919732    0.4680546
#> tcell.0..platelet.0.1        10.111044    6.5233073     8.663352    1.9371079
#> tcell.0..platelet.1.1         6.493596    3.0010893     5.964835    1.4262341
#> tcell.1..platelet.0.1         7.897002    2.4131328     8.648789    1.6694947
#> tcell.1..platelet.1.1         7.212249    1.7479449     7.289246    1.9195989
#> tcell.0..platelet.0.2        16.397057   10.5473267    13.892998    3.2393340
#> tcell.0..platelet.1.2        11.098757    4.9281021     9.616691    2.7658651
#> tcell.1..platelet.0.2        13.808809    3.8467866    13.704558    3.6306276
#> tcell.1..platelet.1.2        12.233894    2.8432825    11.751429    3.4296557
#> tcell.0..platelet.0.3        22.840051   14.5883949    19.159933    4.6614199
#> tcell.0..platelet.1.3        16.116349    7.0338426    13.454119    4.2874001
#> tcell.1..platelet.0.3        20.450776    5.2780561    18.768914    6.2071177
#> tcell.1..platelet.1.3        17.566666    4.2147552    16.354324    4.9275895
#> tcell.0..platelet.0.4        29.517541   18.6610883    24.467147    6.2636181
#> tcell.0..platelet.1.4        21.669304    9.5272969    17.681908    5.9059106
#> tcell.1..platelet.0.4        27.253973    6.7084554    23.836380    8.8563983
#> tcell.1..platelet.1.4        23.884329    6.4184193    21.727949    6.4218060
#> tcell.0..platelet.0.5        36.258527   22.7297215    29.779702    7.9140574
#> tcell.0..platelet.1.5        27.228777   11.9932300    21.982415    7.5062792
#> tcell.1..platelet.0.5        34.057170    8.1384414    28.905320   11.4591447
#> tcell.1..platelet.1.5        30.328361    8.6241905    27.651568    7.9144109
#>                       upper_intF_2
#> tcell.0..platelet.0       1.414055
#> tcell.0..platelet.1       1.148046
#> tcell.1..platelet.0       2.528752
#> tcell.1..platelet.1       2.015356
#> tcell.0..platelet.0.1     3.472269
#> tcell.0..platelet.1.1     3.589546
#> tcell.1..platelet.0.1     6.636320
#> tcell.1..platelet.1.1     6.912774
#> tcell.0..platelet.0.2     5.686570
#> tcell.0..platelet.1.2     6.422076
#> tcell.1..platelet.0.2    11.809821
#> tcell.1..platelet.1.2    12.143500
#> tcell.0..platelet.0.3     8.038703
#> tcell.0..platelet.1.3     9.518797
#> tcell.1..platelet.0.3    17.754145
#> tcell.1..platelet.1.3    17.417653
#> tcell.0..platelet.0.4    10.603777
#> tcell.0..platelet.1.4    12.786703
#> tcell.1..platelet.0.4    24.096859
#> tcell.1..platelet.1.4    22.704955
#> tcell.0..platelet.0.5    13.253487
#> tcell.0..platelet.1.5    16.095748
#> tcell.1..platelet.0.5    30.579901
#> tcell.1..platelet.1.5    27.997957
#> 

## comparing populations, can also be done using rmstIPCW via influence functions
drm1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=40)
summary(drm1)
#> $estimate
#>                     strata times    intF_1    intF_2 se.intF_1 se.intF_2
#> tcell=0, platelet=0      0    40 16.718632  6.121419  1.162627 0.8509985
#> tcell=0, platelet=1      1    40  9.728009  6.388340  1.609499 1.2998358
#> tcell=1, platelet=0      2    40  9.953059 10.497717  3.221205 2.8144104
#> tcell=1, platelet=1      3    40  8.302378  9.264289  2.871788 2.9840876
#>                     total.years.lost lower_intF_1 upper_intF_1 lower_intF_2
#> tcell=0, platelet=0         22.84005    14.588395     19.15993     4.661420
#> tcell=0, platelet=1         16.11635     7.033843     13.45412     4.287400
#> tcell=1, platelet=0         20.45078     5.278056     18.76891     6.207118
#> tcell=1, platelet=1         17.56667     4.214755     16.35432     4.927590
#>                     upper_intF_2
#> tcell=0, platelet=0     8.038703
#> tcell=0, platelet=1     9.518797
#> tcell=1, platelet=0    17.754145
#> tcell=1, platelet=1    17.417653
#> 
## first cause 
e1 <- estimate(drm1)
estimate(e1,rbind(c(1,-1,0,0)))
#>             Estimate Std.Err  2.5% 97.5%   P-value
#> [p1] - [p2]    6.991   1.985 3.099 10.88 0.0004302
#> 
#>  Null Hypothesis: 
#>   [p1] - [p2] = 0 
```
