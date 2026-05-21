# Casewise Concordance from Concordant/Discordant Counts

Computes casewise concordance probability and confidence interval from
counts of concordant and discordant pairs using a binomial GLM.

## Usage

``` r
casewise_bin(nc, nd)
```

## Arguments

- nc:

  number of concordant pairs.

- nd:

  number of discordant pairs.

## Value

A list with `p.casewise` (estimated probability) and `ci.casewise`
(confidence interval).
