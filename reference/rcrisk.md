# Simulation of Piecewise constant hazard models with two causes (Cox).

Simulates data from piecwise constant baseline hazard that can also be
of Cox type. Censor data at highest value of the break points for either
of the cumulatives.

## Usage

``` r
rcrisk(
  cumhaz1,
  cumhaz2,
  rr1,
  rr2,
  n = NULL,
  cens = NULL,
  rrc = NULL,
  extend = TRUE,
  ...
)
```

## Arguments

- cumhaz1:

  cumulative hazard of cause 1

- cumhaz2:

  cumulative hazard of cause 1

- rr1:

  number of simulations or vector of relative risk for simuations.

- rr2:

  number of simulations or vector of relative risk for simuations.

- n:

  number of simulation if rr not given

- cens:

  to censor further , rate or cumumlative hazard

- rrc:

  retlativ risk for censoring.

- extend:

  to extend the cumulative hazards to largest end-point

- ...:

  arguments for rchaz

## Author

Thomas Scheike

## Examples

``` r
data(bmt); 
cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
cox2 <- phreg(Surv(time,cause==2)~tcell+platelet,data=bmt)

X1 <- bmt[,c("tcell","platelet")]
n <- 100
xid <- sample(1:nrow(X1),n,replace=TRUE)
Z1 <- X1[xid,]
Z2 <- X1[xid,]
rr1 <- exp(as.matrix(Z1) %*% cox1$coef)
rr2 <- exp(as.matrix(Z2) %*% cox2$coef)

d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2)
dd <- cbind(d,Z1)

scox1 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
scox2 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
par(mfrow=c(1,2))
plot(cox1); plot(scox1,add=TRUE)
plot(cox2); plot(scox2,add=TRUE)

cbind(cox1$coef,scox1$coef,cox2$coef,scox2$coef)
#>                [,1]        [,2]       [,3]        [,4]
#> tcell    -0.4232606 -0.02131069  0.3991068 -0.82695114
#> platelet -0.5654438 -0.39726694 -0.2461474 -0.09356943
```
