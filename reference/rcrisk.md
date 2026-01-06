# Simulation of Piecewise constant hazard models with two causes (Cox).

Simulates data from piecwise constant baseline hazard that can also be
of Cox type. Censor data at highest value of the break points for either
of the cumulatives, see also sim.phregs

## Usage

``` r
rcrisk(
  cumA,
  cumB,
  rr1 = NULL,
  rr2 = NULL,
  n = NULL,
  cens = NULL,
  rrc = NULL,
  extend = TRUE,
  causes = NULL,
  ...
)
```

## Arguments

- cumA:

  cumulative hazard of cause 1, or list of multiple cumulative hazards

- cumB:

  cumulative hazard of cause 2 or NULL when cumA is a list

- rr1:

  number of simulations or vector of relative risk for simuations, or
  matrix with columns equal to number of hazards in list

- rr2:

  number of simulations or vector of relative risk for simuations.

- n:

  number of simulation if rr not given, must be given when rr is not
  given

- cens:

  to censor further , rate or cumumlative hazard

- rrc:

  retlativ risk for censoring.

- extend:

  to extend the cumulative hazards to largest end-point

- causes:

  to assign status values for each of the causes, vector of integers

- ...:

  arguments for rchaz

## Author

Thomas Scheike

## Examples

``` r
library(mets)
data(bmt); 
n <- 100
cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
cox2 <- phreg(Surv(time,cause==2)~tcell+platelet,data=bmt)

X1 <- bmt[,c("tcell","platelet")]
xid <- sample(1:nrow(X1),n,replace=TRUE)
Z1 <- X1[xid,]
Z2 <- X1[xid,]
rr1 <- exp(as.matrix(Z1) %*% cox1$coef)
rr2 <- exp(as.matrix(Z2) %*% cox2$coef)

d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2,cens=2/70)
dd <- cbind(d,Z1)

d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2,cens=cbind(c(1,30,68),c(.01,1,3)))
dd <- cbind(d,Z1)

par(mfrow=c(1,3))
scox0 <- phreg(Surv(time,status==0)~tcell+platelet,data=dd)
plot(scox0); lines(cbind(c(1,30,68),c(.01,1,3)),col=2)
##
scox1 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
scox2 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
plot(cox1); plot(scox1,add=TRUE,col=2)
plot(cox2); plot(scox2,add=TRUE,col=2)

cbind(cox1$coef,scox1$coef,cox2$coef,scox2$coef)
#>                [,1]        [,2]       [,3]      [,4]
#> tcell    -0.4232606 -12.1354452  0.3991068 1.2674840
#> platelet -0.5654438  -0.3064046 -0.2461474 0.5558998

# 3 causes and censoring 
d3 <-  rcrisk(list(cox1$cum,cox2$cum,cox1$cum),NULL,n=100,cens=cbind(c(1,30,68),c(.01,1,3)))
dtable(d3,~status)
#> 
#> status
#>  0  1  2  3 
#> 19 37 13 31 
#> 
```
