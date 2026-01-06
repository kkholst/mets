# Estimates the casewise concordance based on Concordance and marginal estimate using binreg

Estimates the casewise concordance based on Concordance and marginal
estimate using binreg

## Usage

``` r
binregCasewise(concbreg, margbreg, zygs = c("DZ", "MZ"), newdata = NULL, ...)
```

## Arguments

- concbreg:

  Concordance

- margbreg:

  Marginal estimate

- zygs:

  order of zygosity for estimation of concordance and casewise.

- newdata:

  to give instead of zygs.

- ...:

  to pass to estimate function

## Details

Uses cluster iid for the two binomial-regression estimates standard
errors better than those of casewise that are often conservative.

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(prt)
prt <- force.same.cens(prt,cause="status")

dd <- bicompriskData(Event(time, status)~strata(zyg)+id(id), data=prt, cause=c(2, 2))
newdata <- data.frame(zyg=c("DZ","MZ"),id=1)

## concordance 
bcif1 <- binreg(Event(time,status)~-1+factor(zyg)+cluster(id), data=dd,
                time=80, cause=1, cens.model=~strata(zyg))
pconc <- predict(bcif1,newdata)

## marginal estimates 
mbcif1 <- binreg(Event(time,status)~cluster(id), data=prt, time=80, cause=2)
mc <- predict(mbcif1,newdata)

cse <- binregCasewise(bcif1,mbcif1)
cse
#> $coef
#>     Estimate      2.5%     97.5%
#> p1 0.1586277 0.1445496 0.1740770
#> p2 0.4041311 0.3682646 0.4434908
#> 
#> $logcoef
#>    Estimate Std.Err   2.5%   97.5%   P-value
#> p1   -1.841 0.04742 -1.934 -1.7483 0.000e+00
#> p2   -0.906 0.04742 -0.999 -0.8131 2.208e-81
#> 
```
