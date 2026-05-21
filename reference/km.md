# Kaplan-Meier with Robust Standard Errors

Computes Kaplan-Meier estimates with robust standard errors. Robust
variance is the default and is obtained from the `predict` call.

## Usage

``` r
km(formula, data = data, km = TRUE, ...)
```

## Arguments

- formula:

  Formula with 'Surv' or 'Event' outcome.

- data:

  Data frame.

- km:

  Logical; if `TRUE`, returns Kaplan-Meier estimates; otherwise returns
  Nelson-Aalen based estimates.

- ...:

  Additional arguments passed to `phreg`.

## Value

An object of class `"km"` (extends `"predictphreg"`) containing:

- surr:

  Survival probabilities.

- se.surv:

  Standard errors.

- lower, upper:

  Confidence intervals.

## Author

Thomas Scheike

## Examples

``` r
data(sTRACE)
sTRACE$cluster <- sample(1:100, 500, replace = TRUE)
out1 <- km(Surv(time, status == 9) ~ strata(vf, chf), data = sTRACE)
out2 <- km(Surv(time, status == 9) ~ strata(vf, chf) + cluster(cluster), data = sTRACE)

summary(out1, times = 1:3)
#> Predictions of type 'surv'
#>   Showing subjects: 1, 2, 3, 4
#>   Showing times:    1, 2, 3
#> 
#> -- Subject 1 --
#>  time   surv     se  lower  upper
#>     1 0.9215 0.0177 0.8875 0.9568
#>     2 0.8604 0.0228 0.8169 0.9062
#>     3 0.8167 0.0254 0.7683 0.8681
#> 
#> -- Subject 2 --
#>  time   surv     se  lower  upper
#>     1 0.7175 0.0291 0.6627 0.7768
#>     2 0.6035 0.0316 0.5447 0.6687
#>     3 0.5022 0.0323 0.4428 0.5697
#> 
#> -- Subject 3 --
#>  time   surv     se  lower upper
#>     1 0.6667 0.1566 0.4207     1
#>     2 0.6667 0.1566 0.4207     1
#>     3 0.6667 0.1566 0.4207     1
#> 
#> -- Subject 4 --
#>  time   surv     se  lower  upper
#>     1 0.5217 0.0978 0.3614 0.7533
#>     2 0.3913 0.0942 0.2441 0.6274
#>     3 0.3913 0.0942 0.2441 0.6274
#> 
summary(out2, times = 1:3)
#> Predictions of type 'surv'
#>   Showing subjects: 1, 2, 3, 4
#>   Showing times:    1, 2, 3
#> 
#> -- Subject 1 --
#>  time   surv     se  lower  upper
#>     1 0.9215 0.0157 0.8913 0.9528
#>     2 0.8604 0.0220 0.8184 0.9045
#>     3 0.8167 0.0239 0.7712 0.8649
#> 
#> -- Subject 2 --
#>  time   surv     se  lower  upper
#>     1 0.7175 0.0284 0.6639 0.7753
#>     2 0.6035 0.0305 0.5467 0.6663
#>     3 0.5022 0.0331 0.4414 0.5715
#> 
#> -- Subject 3 --
#>  time   surv     se  lower upper
#>     1 0.6667 0.1566 0.4207     1
#>     2 0.6667 0.1566 0.4207     1
#>     3 0.6667 0.1566 0.4207     1
#> 
#> -- Subject 4 --
#>  time   surv     se  lower  upper
#>     1 0.5217 0.0934 0.3674 0.7409
#>     2 0.3913 0.0899 0.2494 0.6140
#>     3 0.3913 0.0899 0.2494 0.6140
#> 

par(mfrow = c(1, 2))
plot(out1, se = TRUE)
plot(out2, se = TRUE)
```
