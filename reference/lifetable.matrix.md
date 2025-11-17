# Life table

Create simple life table

## Usage

``` r
# S3 method for class 'matrix'
lifetable(x, strata = list(), breaks = c(),
   weights=NULL, confint = FALSE, ...)

 # S3 method for class 'formula'
lifetable(x, data=parent.frame(), breaks = c(),
   weights=NULL, confint = FALSE, ...)
```

## Arguments

- x:

  time formula (Surv) or matrix/data.frame with columns time,status or
  entry,exit,status

- strata:

  strata

- breaks:

  time intervals

- weights:

  weights variable

- confint:

  if TRUE 95% confidence limits are calculated

- ...:

  additional arguments to lower level functions

- data:

  data.frame

## Author

Klaus K. Holst

## Examples

``` r
library(timereg)
data(TRACE)

d <- with(TRACE,lifetable(Surv(time,status==9)~sex+vf,breaks=c(0,0.2,0.5,8.5)))
lava::estimate(glm(events ~ offset(log(atrisk))+factor(int.end)*vf + sex*vf,
            data=d,poisson))
#>                        Estimate  Std.Err    2.5%   97.5%    P-value
#> (Intercept)           -0.444337 0.010009 -0.4640 -0.4247  0.000e+00
#> factor(int.end)0.5    -1.197746 0.024205 -1.2452 -1.1503  0.000e+00
#> factor(int.end)8.5    -1.871838 0.008473 -1.8884 -1.8552  0.000e+00
#> vf                     1.830440 0.064658  1.7037  1.9572 2.631e-176
#> sex                   -0.239036 0.006105 -0.2510 -0.2271  0.000e+00
#> factor(int.end)0.5:vf -1.746744 0.218350 -2.1747 -1.3188  1.247e-15
#> factor(int.end)8.5:vf -1.927748 0.097242 -2.1183 -1.7372  1.838e-87
#> vf:sex                 0.009668 0.081323 -0.1497  0.1691  9.054e-01
```
