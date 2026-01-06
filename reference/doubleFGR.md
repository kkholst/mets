# Double CIF Fine-Gray model with two causes

Estimation based on derived hazards and recursive estimating equations.
fits two parametrizations 1) \$\$ F_1(t,X) = 1 - \exp( - \exp( X^T \beta
) \Lambda_1(t)) \$\$ and \$\$ F_2(t,X_2) = 1 - \exp( - \exp( X_2^T
\beta_2 ) \Lambda_2(t)) \$\$ or restricted version 2) \$\$ F_1(t,X) =
1 - \exp( - \exp( X^T \beta ) \Lambda_1(t)) \$\$ and \$\$ F_2(t,X_2,X) =
( 1 - \exp( - \exp( X_2^T \beta_2 ) \Lambda_2(t)) ) (1 - F_1(\infty,X))
\$\$

## Usage

``` r
doubleFGR(formula, data, offset = NULL, weights = NULL, X2 = NULL, ...)
```

## Arguments

- formula:

  formula with 'Event'

- data:

  data frame

- offset:

  offsets for cox model

- weights:

  weights for Cox score equations

- X2:

  specifies the regression design for second CIF model

- ...:

  Additional arguments to lower level funtions

## Author

Thomas Scheike

## Examples

``` r
library(mets)
res <- 0
data(bmt)
bmt$age2 <- bmt$age
newdata <- bmt[1:19,]
## if (interactive()) par(mfrow=c(5,3))

## same X1 and X2
pr2 <- doubleFGR(Event(time,cause)~age+platelet,data=bmt,restrict=res)
##if (interactive()) {
##  bplotdFG(pr2,cause=1)
##  bplotdFG(pr2,cause=2,add=TRUE)
##}
##pp21 <- predictdFG(pr2,newdata=newdata)
##pp22 <- predictdFG(pr2,newdata=newdata,cause=2)
##if (interactive()) {
##  plot(pp21)
##  plot(pp22,add=TRUE,col=2)
##}
##pp21 <- predictdFG(pr2)
##pp22 <- predictdFG(pr2,cause=2)
##if (interactive()) {
 ## plot(pp21)
 ## plot(pp22,add=TRUE,col=2)
##}

pr2 <- doubleFGR(Event(time,cause)~strata(platelet),data=bmt,restrict=res)

## different X1 and X2
pr2 <- doubleFGR(Event(time,cause)~age+platelet+age2,data=bmt,X2=3,restrict=res)

### uden X1
pr2 <- doubleFGR(Event(time,cause)~age+platelet,data=bmt,X2=1:2,restrict=res)

### without X2
pr2 <- doubleFGR(Event(time,cause)~age+platelet,data=bmt,X2=0,restrict=res)
```
