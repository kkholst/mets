# Fast Additive Hazards Model with Robust Standard Errors

Fits a fast Lin-Ying additive hazards model with a possibly stratified
baseline. Robust variance is the default variance estimate in the
summary.

## Usage

``` r
aalenMets(formula, data = data, no.baseline = FALSE, ...)
```

## Arguments

- formula:

  Formula with a 'Surv' outcome (similar to `coxph`).

- data:

  Data frame.

- no.baseline:

  Logical; if `TRUE`, fits the model without baseline hazard estimation.

- ...:

  Additional arguments passed to `phreg`.

## Value

An object of class `"aalenMets"` (extends `"phreg"`) containing:

- coef:

  Estimated coefficients.

- var:

  Robust variance-covariance matrix.

- iid:

  Influence functions.

- intZHZ:

  Integrated ZHZ matrix.

- gamma:

  Coefficient estimates.

## Details

Influence functions (IID) follow the numerical order of the given
cluster variable. Ordering by `$id` aligns the IID terms with the
dataset order.

## Author

Thomas Scheike

## Examples

``` r
data(bmt)
bmt$time <- bmt$time + runif(408) * 0.001
out <- aalenMets(Surv(time, cause == 1) ~ tcell + platelet + age, data = bmt)
summary(out)
#> 
#>    n events
#>  408    161
#> 
#>  408 clusters
#> coefficients:
#>            Estimate       S.E.    dU^-1/2 P-value
#> tcell    -0.0129507  0.0041288  0.2304083  0.0017
#> platelet -0.0087471  0.0028053  0.1664259  0.0018
#> age       0.0066174  0.0013879  0.0789210  0.0000
#> 
#> exp(coefficients):
#>          Estimate    2.5%  97.5%
#> tcell     0.98713 0.97918 0.9952
#> platelet  0.99129 0.98586 0.9968
#> age       1.00664 1.00390 1.0094
#> 

## Comparison with timereg::aalen
## out2 <- timereg::aalen(
##   Surv(time, cause == 1) ~ const(tcell) + const(platelet) + const(age),
##   data = bmt
## )
## summary(out2)
```
