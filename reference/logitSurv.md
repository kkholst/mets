# Proportional Odds Survival Model

Fits a semiparametric proportional odds model where: \$\$
\mbox{logit}(S(t\|x)) = \log(\Lambda(t)) + x \beta \$\$ Thus, covariate
effects represent the odds ratio (OR) of survival.

## Usage

``` r
logitSurv(formula, data, offset = NULL, weights = NULL, ...)
```

## Arguments

- formula:

  Formula with 'Surv' outcome (similar to `coxph`).

- data:

  Data frame.

- offset:

  Offsets for \\\exp(x \beta)\\ terms.

- weights:

  Weights for score equations.

- ...:

  Additional arguments passed to lower-level functions.

## Value

An object of class `"phreg"` with `propodds=1`.

## Details

This is equivalent to using a hazards model: \$\$ Z \lambda(t) \exp(x
\beta) \$\$ where \\Z\\ is gamma distributed with mean and variance 1.

## References

Eriksson, Frank, Li, Jianing, Scheike, Thomas, and Zhang, Mei-Jie
(2015). "The proportional odds cumulative incidence model for competing
risks." *Biometrics*, 71(3), 687–695.

## Author

Thomas Scheike

## Examples

``` r
data(TRACE)
dcut(TRACE) <- ~.
out1 <- logitSurv(Surv(time, status == 9) ~ vf + chf + strata(wmicat.4), data = TRACE)
summary(out1)
#> 
#>     n events
#>  1878    958
#> coefficients:
#>     Estimate    S.E. dU^-1/2 P-value
#> vf   0.30049 0.22633 0.11154  0.1843
#> chf  1.26008 0.10095 0.07316  0.0000
#> 
#> exp(coefficients):
#>     Estimate    2.5%  97.5%
#> vf   1.35052 0.86667 2.1045
#> chf  3.52570 2.89277 4.2971
#> 
gof(out1)
#> Cumulative score process test for Proportionality:
#>     Sup|U(t)|  pval
#> vf   42.28659 0.000
#> chf  20.75308 0.485
plot(out1)
```
