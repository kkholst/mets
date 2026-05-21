# G-Estimator for Binomial Regression Model (Standardized Estimates)

Computes the G-estimator (G-formula) for standardized risk estimates
based on a fitted `binreg` object. The G-estimator standardizes
predictions over the covariate distribution in the data: \$\$ \hat F(t,
A=a) = n^{-1} \sum\_{i=1}^n \hat F(t, A=a, Z_i) \$\$

## Usage

``` r
binregG(x, data, Avalues = NULL, varname = NULL)
```

## Arguments

- x:

  An object of class `"binreg"` obtained from
  [`binreg()`](http://kkholst.github.io/mets/reference/binreg.md) or
  [`logitIPCW()`](http://kkholst.github.io/mets/reference/binreg.md).

- data:

  A data frame containing the covariates to be used for averaging the
  risk estimates. This should ideally be the same data used to fit the
  model, or a representative sample.

- Avalues:

  Numeric or factor vector specifying the values of the first covariate
  (\\A\\) for which to compute standardized risks.

  - If the first covariate is a factor and `Avalues` is `NULL`, all
    levels of the factor are used.

  - If the first covariate is continuous, `Avalues` must be provided.

- varname:

  Optional character string specifying the name of the variable to be
  treated as the treatment/exposure variable. If `NULL`, the first
  variable in the model formula is used.

## Value

An object of class `"survivalG"` containing:

- risk:

  A table of standardized risk estimates for each value of `Avalues`.

- risk.iid:

  Influence functions for the standardized risk estimates.

- difference:

  Pairwise differences in risks between levels of `A`.

- ratio:

  Risk ratios between levels of `A`.

- vcov:

  Variance-covariance matrix of the risk estimates.

- model:

  The link function used in the original model.

## Details

This function assumes that the first covariate in the original model
formula represents the treatment or exposure variable (\\A\\). It
calculates the marginal risk for specified values of \\A\\ by averaging
the conditional predictions over the observed covariate distribution
\\Z\\.

The function returns influence functions for these risk estimates,
allowing for the computation of standard errors and confidence
intervals.

If the first covariate is a factor, contrasts between all levels are
computed automatically. If it is continuous, specific values must be
provided via `Avalues`.

## References

- Blanche PF, Holt A, Scheike T (2022). "On logistic regression with
  right censored data, with or without competing risks, and its use for
  estimating treatment effects." *Lifetime Data Analysis*, 29, 441–482.

## See also

[`binreg`](http://kkholst.github.io/mets/reference/binreg.md),
[`binregATE`](http://kkholst.github.io/mets/reference/binregATE.md)

## Author

Thomas Scheike

## Examples

``` r
data(bmt); bmt$time <- bmt$time+runif(408)*0.001
bmt$event <- (bmt$cause!=0)*1

b1 <- binreg(Event(time,cause)~age+tcell+platelet,bmt,cause=1,time=50)
sb1 <- binregG(b1,bmt,Avalues=c(0,1,2))
summary(sb1)
#> G-estimator :
#>       Estimate Std.Err   2.5%  97.5%   P-value
#> risk0   0.4058 0.02588 0.3551 0.4566 1.980e-55
#> risk1   0.5120 0.03706 0.4393 0.5846 2.087e-43
#> risk2   0.6169 0.05517 0.5087 0.7250 5.018e-29
#> 
#> Average Treatment effect: difference (G-estimator) :
#>      Estimate Std.Err    2.5%  97.5%   P-value
#> pa     0.1061 0.02623 0.05471 0.1575 5.224e-05
#> pa.1   0.2110 0.04960 0.11380 0.3082 2.096e-05
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>     Estimate    Std.Err      2.5%     97.5%      P-value
#> pa 0.2322908 0.05276944 0.1288646 0.3357170 1.072603e-05
#> pa 0.4186840 0.08402001 0.2540078 0.5833602 6.255859e-07
#> ratio: 
#>    Estimate     2.5%    97.5%
#> pa 1.261486 1.137536 1.398943
#> pa 1.519960 1.289182 1.792050
#> 
```
