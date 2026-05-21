# Average Treatment Effect for Restricted Mean Time

Estimates the Average Treatment Effect (ATE) for Restricted Mean
Survival Time (RMST) or Restricted Mean Time Lost (RMTL) in censored
competing risks data using IPCW.

## Usage

``` r
resmeanATE(formula, data, model = "exp", outcome = c("rmst", "rmtl"), ...)
```

## Arguments

- formula:

  Formula with an `Event` outcome. The first covariate must be the
  treatment factor.

- data:

  Data frame.

- model:

  Link function: `"exp"` (exponential) or `"lin"` (identity).

- outcome:

  Outcome type: `"rmst"` or `"rmtl"`.

- ...:

  Additional arguments passed to `binregATE`, such as `time`,
  `treat.model`, `augmentR0`, `augmentC`, etc.

## Value

An object of class `"binregATE"` containing:

- riskG:

  Simple IPCW estimator results.

- riskDR:

  Double Robust estimator results.

- riskG.iid, riskDR.iid:

  Influence functions.

- coef:

  Treatment effect estimates.

- se:

  Standard errors.

## Details

Under standard causal assumptions (Consistency, Ignorability,
Positivity), the ATE is estimated as \\E(Y(1) - Y(0))\\, where \\Y(a)\\
is the potential outcome under treatment \\a\\. The method uses double
robust estimating equations that are IPCW-adjusted for censoring.

The first covariate in the formula must be the treatment effect (a
factor). If the factor has more than two levels, multinomial logistic
regression (mlogit) is used for propensity score modeling.

## References

Scheike, T. and Holst, K. K. (2024). Restricted mean time lost for
survival and competing risks data using mets in R. WIP.

Scheike, T. and Tanaka, S. (2025). Restricted mean time lost ratio
regression: Percentage of restricted mean time lost due to specific
cause. WIP.

## See also

[`binregATE`](http://kkholst.github.io/mets/reference/binregATE.md),
[`ratioATE`](http://kkholst.github.io/mets/reference/ratioATE.md)

## Author

Thomas Scheike

## Examples

``` r
data(bmt); bmt$event <- bmt$cause!=0; dfactor(bmt) <- tcell~tcell

out <- resmeanATE(Event(time,event)~tcell+platelet, data=bmt, time=40, 
                  treat.model=tcell~platelet, outcome="rmtl")
summary(out)
#>    n events
#>  408    241
#> 
#>  408 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept)  3.121409  0.048611  3.026134  3.216684  0.0000
#> tcell1      -0.022691  0.131519 -0.280465  0.235082  0.8630
#> platelet    -0.318830  0.106383 -0.527336 -0.110324  0.0027
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 22.67831 20.61737 24.9453
#> tcell1       0.97756  0.75543  1.2650
#> platelet     0.72700  0.59018  0.8955
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    20.73597  0.95653 18.86121 22.61073  0.0000
#> treat1    20.27074  2.47786 15.41422 25.12726  0.0000
#> treat:1-0 -0.46523  2.67400 -5.70617  4.77572  0.8619
#> 
#> Average Treatment effects (double robust) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    20.71846  0.95824 18.84035 22.59657  0.0000
#> treat1    19.65975  2.50069 14.75848 24.56102  0.0000
#> treat:1-0 -1.05871  2.67373 -6.29913  4.18171  0.6921
#> 
#> 

out1 <- resmeanATE(Event(time,cause)~tcell+platelet, data=bmt, cause=1, time=40,
                   treat.model=tcell~platelet, outcome="rmtl")
summary(out1)
#>    n events
#>  408    157
#> 
#>  408 clusters
#> coeffients:
#>             Estimate  Std.Err     2.5%    97.5% P-value
#> (Intercept)  2.80626  0.06962  2.66981  2.94271  0.0000
#> tcell1      -0.37413  0.24769 -0.85960  0.11133  0.1309
#> platelet    -0.49164  0.16493 -0.81490 -0.16837  0.0029
#> 
#> exp(coeffients):
#>             Estimate     2.5%   97.5%
#> (Intercept) 16.54790 14.43717 18.9672
#> tcell1       0.68788  0.42333  1.1178
#> platelet     0.61162  0.44268  0.8450
#> 
#> Average Treatment effects (G-formula) :
#>           Estimate  Std.Err     2.5%    97.5% P-value
#> treat0    14.53165  0.95705 12.65587 16.40742  0.0000
#> treat1     9.99609  2.37815  5.33499 14.65718  0.0000
#> treat:1-0 -4.53556  2.57515 -9.58276  0.51164  0.0782
#> 
#> Average Treatment effects (double robust) :
#>            Estimate   Std.Err      2.5%     97.5% P-value
#> treat0     14.51355   0.95800  12.63590  16.39120  0.0000
#> treat1      9.36465   2.41708   4.62727  14.10203  0.0001
#> treat:1-0  -5.14890   2.59798 -10.24084  -0.05696  0.0475
#> 
#> 
```
