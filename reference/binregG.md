# G-estimator for binomial regression model (Standardized estimates)

Computes G-estimator \$\$ \hat F(t,A=a) = n^{-1} \sum_i \hat
F(t,A=a,Z_i) \$\$. Assumes that the first covariate is \$A\$. Gives
influence functions of these risk estimates and SE's are based on these.
If first covariate is a factor then all contrast are computed, and if
continuous then considered covariate values are given by Avalues.

## Usage

``` r
binregG(x, data, Avalues = c(0, 1), varname = NULL)
```

## Arguments

- x:

  binreg object

- data:

  data frame for risk averaging

- Avalues:

  values to compare for first covariate A

- varname:

  if given then averages for this variable, default is first variable

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
#> risk0   0.4058 0.02588 0.3551 0.4566 1.981e-55
#> risk1   0.5120 0.03706 0.4393 0.5846 2.092e-43
#> risk2   0.6169 0.05517 0.5087 0.7250 5.030e-29
#> 
#> Average Treatment effect: difference (G-estimator) :
#>    Estimate Std.Err   2.5%  97.5%   P-value
#> p1   0.1061 0.02623 0.0547 0.1575 5.227e-05
#> p2   0.2110 0.04960 0.1138 0.3082 2.098e-05
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>       Estimate    Std.Err      2.5%     97.5%      P-value
#> [p1] 0.2322843 0.05276981 0.1288573 0.3357112 1.073368e-05
#> [p2] 0.4186733 0.08402102 0.2539951 0.5833514 6.261926e-07
#> ratio: 
#>      Estimate     2.5%    97.5%
#> [p1] 1.261478 1.137528 1.398935
#> [p2] 1.519944 1.289165 1.792034
#> 
```
