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
#> risk0   0.4058 0.02588 0.3551 0.4565 1.972e-55
#> risk1   0.5119 0.03705 0.4393 0.5845 2.062e-43
#> risk2   0.6168 0.05516 0.5087 0.7249 5.047e-29
#> 
#> Average Treatment effect: difference (G-estimator) :
#> Call: estimate.default(coef = Gest$Gest, vcov = vv)
#> ────────────────────────────────────────────────────────────
#>      Estimate Std.Err   2.5%  97.5%   P-value
#> pa     0.4058 0.02588 0.3551 0.4565 1.972e-55
#> pa.1   0.5119 0.03705 0.4393 0.5845 2.062e-43
#> pa.2   0.6168 0.05516 0.5087 0.7249 5.047e-29
#> ────────────────────────────────────────────────────────────
#> Null Hypothesis: 
#>   [pa] = 0
#>   [pa] = 0
#>   [pa] = 0 
#>  
#> chisq = 259.0715, df = 3, p-value < 2.2e-16
#> 
#> Average Treatment effect: ratio (G-estimator) :
#> log-ratio: 
#>     Estimate    Std.Err      2.5%     97.5%      P-value
#> pa 0.4058263 0.02587658 0.3551092 0.4565435 1.972151e-55
#> pa 0.5119228 0.03705489 0.4392966 0.5845491 2.062121e-43
#> pa 0.6168006 0.05516466 0.5086799 0.7249214 5.047007e-29
#> ratio: 
#>    Estimate     2.5%    97.5%
#> pa 1.500542 1.426336 1.578608
#> pa 1.668496 1.551615 1.794182
#> pa 1.852990 1.663094 2.064569
#> 
```
