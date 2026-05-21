# While-Alive Regression for Recurrent Events

Performs regression analysis for the "While-Alive" mean of events per
time unit, defined as \\Z(t) = N(\min(D,t)) / \min(D,t)\\. This function
models how covariates affect the rate of recurrent events per unit of
time alive.

## Usage

``` r
WA_reg(
  formula,
  data,
  time = NULL,
  cens.code = 0,
  cause = 1,
  death.code = 2,
  marks = NULL,
  ...,
  trans = 1
)
```

## Arguments

- formula:

  Formula with regression design. The first covariate on the RHS must be
  the treatment factor. Can include other covariates and `cluster(id)`.

- data:

  Data frame.

- time:

  Time point \\t\\ for estimation.

- cens.code:

  Censoring code.

- cause:

  Event cause code.

- death.code:

  Death code.

- marks:

  Marks for composite outcomes.

- ...:

  Additional arguments passed to `binreg`.

- trans:

  Power transformation for the outcome (default 1).

## Value

An object of class `"binreg"` containing coefficient estimates, standard
errors, confidence intervals, and influence functions for the regression
of the event rate per time alive.

## Details

The estimation is based on IPCW (Inverse Probability of Censoring
Weighting) and calls `binreg` after constructing the outcome variable.
It supports double robust estimation if covariate augmentation is
specified.

## References

Ragni, A., Martinussen, T., & Scheike, T. H. (2023). Nonparametric
estimation of the Patient Weighted While-Alive Estimand. arXiv preprint.

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
hfactioncpx12$age <- rnorm(741)[hfactioncpx12$id] 
dtable(hfactioncpx12,~status)
#> 
#> status
#>    0    1    2 
#>  617 1391  124 
#> 
## exp-link regression 
dd <- WA_reg(Event(entry,time,status)~treatment+age+cluster(id),data=hfactioncpx12,
                    time=2,death.code=2)
summary(dd)
#>    n events
#>  741     86
#> 
#>  741 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  0.070899  0.112756 -0.150098  0.291896  0.5295
#> treatment1  -0.354038  0.143843 -0.635964 -0.072112  0.0138
#> age          0.053433  0.070949 -0.085624  0.192491  0.4514
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  1.07347 0.86062 1.3390
#> treatment1   0.70185 0.52942 0.9304
#> age          1.05489 0.91794 1.2123
#> 
#> 
```
