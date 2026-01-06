# G-estimator for binomial regression model (Standardized estimates)

Computes G-estimator \$\$ \hat F(t,A=a) = n^{-1} \sum_i \hat
F(t,A=a,Z_i) \$\$. Assumes that the first covariate is \$A\$. Gives
influence functions of these risk estimates and SE's are based on these.
If first covariate is a factor then all contrast are computed, and if
continuous then considered covariate values are given by Avalues.

## Usage

``` r
binregG(x, data, Avalues = NULL, varname = NULL)
```

## Arguments

- x:

  binreg object

- data:

  data frame for risk averaging

- Avalues:

  values to compare for first covariate A, assumes that first variable
  is factor and take all levels

- varname:

  if given then averages for this variable, default is first variable

## References

Blanche PF, Holt A, Scheike T (2022). “On logistic regression with right
censored data, with or without competing risks, and its use for
estimating treatment effects.” Lifetime data analysis, 29, 441–482.

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt); bmt$time <- bmt$time+runif(408)*0.001
bmt$event <- (bmt$cause!=0)*1

b1 <- binreg(Event(time,cause)~age+tcell+platelet,bmt,cause=1,time=50)
sb1 <- binregG(b1,bmt,Avalues=c(0,1,2))
summary(sb1)
#> G-estimator :
#>       Estimate Std.Err   2.5%  97.5%   P-value
#> risk0   0.4058 0.02588 0.3551 0.4565 1.982e-55
#> risk1   0.5119 0.03706 0.4393 0.5846 2.057e-43
#> risk2   0.6168 0.05516 0.5087 0.7250 4.993e-29
#> 
#> Average Treatment effect: difference (G-estimator) :
#>      Estimate Std.Err    2.5%  97.5%   P-value
#> pa     0.1061 0.02623 0.05471 0.1575 5.222e-05
#> pa.1   0.2110 0.04960 0.11381 0.3082 2.096e-05
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>       Estimate    Std.Err      2.5%     97.5%      P-value
#> [pa] 0.2323087 0.05277448 0.1288726 0.3357448 1.073002e-05
#> [pa] 0.4187166 0.08402886 0.2540231 0.5834101 6.260295e-07
#> ratio: 
#>      Estimate     2.5%    97.5%
#> [pa] 1.261509 1.137545 1.398982
#> [pa] 1.520010 1.289202 1.792139
#> 
```
