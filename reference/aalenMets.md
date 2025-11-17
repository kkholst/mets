# Fast additive hazards model with robust standard errors

Fast Lin-Ying additive hazards model with a possibly stratified
baseline. Robust variance is default variance with the summary.

## Usage

``` r
aalenMets(formula, data = data, no.baseline = FALSE, ...)
```

## Arguments

- formula:

  formula with 'Surv' outcome (see `coxph`)

- data:

  data frame

- no.baseline:

  to fit model without baseline hazard

- ...:

  Additional arguments to phreg

## Details

influence functions (iid) will follow numerical order of given cluster
variable so ordering after \$id will give iid in order of data-set.

## Author

Thomas Scheike

## Examples

``` r
data(bmt); bmt$time <- bmt$time+runif(408)*0.001
out <- aalenMets(Surv(time,cause==1)~tcell+platelet+age,data=bmt)
summary(out)
#> 
#>    n events
#>  408    161
#> 
#>  408 clusters
#> coeffients:
#>            Estimate       S.E.    dU^-1/2 P-value
#> tcell    -0.0129630  0.0041295  0.2303820  0.0017
#> platelet -0.0087405  0.0028057  0.1664364  0.0018
#> age       0.0066215  0.0013881  0.0789239  0.0000
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> tcell     0.98712 0.97916 0.9951
#> platelet  0.99130 0.98586 0.9968
#> age       1.00664 1.00391 1.0094
#> 

## out2 <- timereg::aalen(Surv(time,cause==1)~const(tcell)+const(platelet)+const(age),data=bmt)
## summary(out2)
```
