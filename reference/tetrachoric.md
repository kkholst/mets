# Estimate parameters from odds-ratio

Calculate tetrachoric correlation of probabilities from odds-ratio

## Usage

``` r
tetrachoric(P, OR, approx = 0, ...)
```

## Arguments

- P:

  Joint probabilities or marginals (if OR is given)

- OR:

  Odds-ratio

- approx:

  If TRUE an approximation of the tetrachoric correlation is used

- ...:

  Additional arguments

## Examples

``` r
tetrachoric(0.3,1.25) # Marginal p1=p2=0.3, OR=2
#> [1] 0.08173353
P <- matrix(c(0.1,0.2,0.2,0.5),2)
prod(diag(P))/prod(lava::revdiag(P))
#> [1] 1.25
##mets:::assoc(P)
tetrachoric(P)
#> [1] 0.08173353
or2prob(2,0.1)
#>            [,1]       [,2]
#> [1,] 0.01690481 0.08309519
#> [2,] 0.08309519 0.81690481
#> attr(,"marg")
#> [1] 0.1 0.1
or2prob(2,c(0.1,0.2))
#>            [,1]      [,2]
#> [1,] 0.03153416 0.1684658
#> [2,] 0.06846584 0.7315342
#> attr(,"marg")
#> [1] 0.1 0.2
```
