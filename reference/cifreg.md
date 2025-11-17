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

  certain event codes to be ignored when finding competing causes

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
#> tcell    -0.709525  0.331991  0.274927  0.0326
#> platelet -0.455274  0.236012  0.187919  0.0537
#> age       0.391162  0.098039  0.083670  0.0001
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> tcell     0.49188 0.25661 0.9429
#> platelet  0.63427 0.39938 1.0073
#> age       1.47870 1.22019 1.7920
#> 
par(mfrow=c(1,2))
plot(or)
nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
por <- predict(or,nd)
plot(por)


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
#> tcell    -0.596733  0.270480  0.275784  0.0274
#> platelet -0.425608  0.180733  0.187722  0.0185
#> age       0.343742  0.080269  0.086284  0.0000
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> tcell     0.55061 0.32405 0.9356
#> platelet  0.65337 0.45848 0.9311
#> age       1.41021 1.20493 1.6505
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
#> platelet -0.424327  0.180814  0.187824  0.0189
#> age       0.341961  0.079861  0.086283  0.0000
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> platelet  0.65421 0.45900 0.9325
#> age       1.40771 1.20374 1.6462
#> 
plot(sfg)

### predictions with CI based on iid decomposition of baseline and beta
### these are used in the predict function above
fg <- cifregFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
Biid <- iidBaseline(fg,time=20)
pfg1 <- FGprediid(Biid,nd)
pfg1
#>           pred     se-log     lower     upper
#> [1,] 0.2693191 0.22758875 0.1724024 0.4207179
#> [2,] 0.4344048 0.07477547 0.3751851 0.5029717
```
