# Cumulative incidence with robust standard errors

Cumulative incidence with robust standard errors

## Usage

``` r
cif(formula, data = data, cause = 1, cens.code = 0, death.code = NULL, ...)
```

## Arguments

- formula:

  formula with 'Event' outcome and strata (only!)

- data:

  data frame

- cause:

  NULL looks at all, otherwise specify which cause to consider

- cens.code:

  censoring code "0" is default, and death is cens.code!=0

- death.code:

  alternative to cens.code give codes of death

- ...:

  Additional arguments to lower level funtions

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt)
bmt$cluster <- sample(1:100,408,replace=TRUE)
out1 <- cif(Event(time,cause)~+1,data=bmt,cause=1)
out2 <- cif(Event(time,cause)~+1+cluster(cluster),data=bmt,cause=1)

par(mfrow=c(1,2))
plot(out1,se=TRUE)
plot(out2,se=TRUE)
```
