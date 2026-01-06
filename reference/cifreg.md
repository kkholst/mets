# CIF regression

CIF logistic-link for propodds=1 default and CIF Fine-Gray (cloglog)
regression for propodds=NULL. The FG model can also be called using the
cifregFG function that has propodds=NULL.

## Usage

``` r
cifreg(
  formula,
  data,
  propodds = 1,
  cause = 1,
  cens.code = 0,
  no.codes = NULL,
  death.code = NULL,
  ...
)
```

## Arguments

- formula:

  formula with 'Event' outcome

- data:

  data frame

- propodds:

  to fit logit link model, and propodds=NULL to fit Fine-Gray model

- cause:

  of interest

- cens.code:

  code of censoring

- no.codes:

  certain event codes to be ignored when finding competing causes, can
  be used with administrative censoring.

- death.code:

  can also specify death.code (in addition to cause) to overrule default
  which takes all remaining codes (minus cause,cens.code,no.codes)

- ...:

  Additional arguments to recreg

## Details

For FG model: \$\$ \int (X - E ) Y_1(t) w(t) dM_1 \$\$ is computed and
summed over clusters and returned multiplied with inverse of second
derivative as iid.naive. Here \$\$w(t) = G(t) (I(T_i \wedge t \<
C_i)/G_c(T_i \wedge t))\$\$ and \$\$E(t) = S_1(t)/S_0(t)\$\$ and
\$\$S_j(t) = \sum X_i^j Y\_{i1}(t) w_i(t) \exp(X_i^T \beta)\$\$.

The iid decomposition of the beta's, however, also have a censoring term
that is also is computed and added (still scaled with inverse second
derivative) \$\$ \int (X - E ) Y_1(t) w(t) dM_1 + \int q(s)/p(s) dM_c
\$\$ and returned as the iid

For logistic link standard errors are slightly to small since
uncertainty from recursive baseline is not considered, so for smaller
data-sets it is recommended to use the prop.odds.subdist of timereg that
is also more efficient due to use of different weights for the
estimating equations. Alternatively, one can also bootstrap the standard
errors.

## Author

Thomas Scheike

## Examples

``` r
## data with no ties
library(mets)
data(bmt,package="mets")
bmt$time <- bmt$time+runif(nrow(bmt))*0.01
bmt$id <- 1:nrow(bmt)

## logistic link  OR interpretation
or=cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
summary(or)
#> 
#>    n events
#>  408    161
#> 
#>  408 clusters
#> coeffients:
#>           Estimate      S.E.   dU^-1/2 P-value
#> tcell    -0.709569  0.331982  0.274926  0.0326
#> platelet -0.455401  0.236010  0.187919  0.0537
#> age       0.391196  0.098040  0.083671  0.0001
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> tcell     0.49186 0.25660 0.9428
#> platelet  0.63419 0.39933 1.0072
#> age       1.47875 1.22023 1.7920
#> 
par(mfrow=c(1,2))
plot(or)
nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
por <- predict(or,nd)
plot(por)


## approximate standard errors 
por <-mets:::predict.phreg(or,nd)
plot(por,se=1)

## Fine-Gray model
fg=cifregFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
summary(fg)
#> 
#>    n events
#>  408    161
#> 
#>  408 clusters
#> coeffients:
#>           Estimate      S.E.   dU^-1/2 P-value
#> tcell    -0.596588  0.270514  0.275786  0.0274
#> platelet -0.425558  0.180762  0.187724  0.0186
#> age       0.343748  0.080265  0.086283  0.0000
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> tcell     0.55069 0.32407 0.9358
#> platelet  0.65340 0.45848 0.9312
#> age       1.41022 1.20494 1.6505
#> 
##fg=recreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,death.code=2)
##summary(fg)
plot(fg)

nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
pfg <- predict(fg,nd,se=1)
plot(pfg,se=1)

## bt <- iidBaseline(fg,time=30)
## bt <- IIDrecreg(fg$cox.prep,fg,time=30)

## not run to avoid timing issues
## gofFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)

sfg <- cifregFG(Event(time,cause)~strata(tcell)+platelet+age,data=bmt,cause=1)
summary(sfg)
#> 
#>    n events
#>  408    161
#> 
#>  408 clusters
#> coeffients:
#>           Estimate      S.E.   dU^-1/2 P-value
#> platelet -0.424172  0.180849  0.187824   0.019
#> age       0.341941  0.079859  0.086283   0.000
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> platelet  0.65431 0.45903 0.9327
#> age       1.40768 1.20372 1.6462
#> 
plot(sfg)


### predictions with CI based on iid decomposition of baseline and beta
### these are used in the predict function above
fg <- cifregFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
Biid <- iidBaseline(fg,time=20)
pfg1 <- FGprediid(Biid,nd)
pfg1
#>           pred    se-log     lower     upper
#> [1,] 0.2693322 0.2276224 0.1723994 0.4207661
#> [2,] 0.4343766 0.0747748 0.3751613 0.5029385
```
