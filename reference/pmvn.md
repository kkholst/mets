# Multivariate normal distribution function

Multivariate normal distribution function

## Usage

``` r
pmvn(lower, upper, mu, sigma, cor = FALSE)
```

## Arguments

- lower:

  lower limits

- upper:

  upper limits

- mu:

  mean vector

- sigma:

  variance matrix or vector of correlation coefficients

- cor:

  if TRUE sigma is treated as standardized (correlation matrix)

## Examples

``` r
lower <- rbind(c(0,-Inf),c(-Inf,0))
upper <- rbind(c(Inf,0),c(0,Inf))
mu <- rbind(c(1,1),c(-1,1))
sigma <- diag(2)+1
pmvn(lower=lower,upper=upper,mu=mu,sigma=sigma)
#> [1] 0.1265479 0.5361516
```
