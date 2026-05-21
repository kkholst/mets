# Simulation of Illness-Death Model

Simulates data from a full illness-death model with reversible
transitions and multiple causes of death. Supports various dependence
structures via shared frailties.

## Usage

``` r
sim_multistate(
  n,
  cumhaz,
  cumhaz2,
  death.cumhaz,
  death.cumhaz2,
  rr12 = NULL,
  rr21 = NULL,
  rd13 = NULL,
  rd23 = NULL,
  rrc = NULL,
  gap.time = FALSE,
  max.recurrent = 100,
  dependence = 0,
  var.z = 0.5,
  cor.mat = NULL,
  cens = NULL,
  extend = TRUE,
  ...
)
```

## Arguments

- n:

  Number of IDs.

- cumhaz:

  Cumulative hazard from state 1 to 2.

- cumhaz2:

  Cumulative hazard from state 2 to 1.

- death.cumhaz:

  Cumulative hazard of death from state 1.

- death.cumhaz2:

  Cumulative hazard of death from state 2.

- rr12:

  Relative risk for 1-\>2.

- rr21:

  Relative risk for 2-\>1.

- rd13:

  Relative risk for death 1-\>3.

- rd23:

  Relative risk for death 2-\>3.

- rrc:

  Relative risk for censoring.

- gap.time:

  Gap time indicator. If true simulates gap-times with specified
  cumulative hazard.

- max.recurrent:

  Maximum recurrent events.

- dependence:

  Dependence structure (0-3).

- var.z:

  Variance of random effects.

- cor.mat:

  Correlation matrix.

- cens:

  Censoring rate.

- extend:

  Extend hazards.

- ...:

  Additional arguments.

## Value

Data frame with multi-state event history.

## Author

Thomas Scheike

## Examples

``` r
########################################
## getting some rates to mimick 
########################################
data(CPH_HPN_CRBSI)
dr <- CPH_HPN_CRBSI$terminal
base1 <- CPH_HPN_CRBSI$crbsi 
base4 <- CPH_HPN_CRBSI$mechanical
dr2 <- scalecumhaz(dr,1.5)
cens <- rbind(c(0,0),c(2000,0.5),c(5110,3))

iddata <- sim_multistate(100,base1,base1,dr,dr2,cens=cens)
dlist(iddata,.~id|id<3,n=0)
#> id: 1
#>   entry     time status rr death from to start     stop
#> 1     0 618.1512      3  1     1    1  3     0 618.1512
#> ------------------------------------------------------------ 
#> id: 2
#>   entry     time status rr death from to start     stop
#> 2     0 243.4911      3  1     1    1  3     0 243.4911
 
### estimating rates from simulated data  
c0 <- phreg(Surv(start,stop,status==0)~+1,iddata)
c3 <- phreg(Surv(start,stop,status==3)~+strata(from),iddata)
c1 <- phreg(Surv(start,stop,status==1)~+1,subset(iddata,from==2))
c2 <- phreg(Surv(start,stop,status==2)~+1,subset(iddata,from==1))
###
par(mfrow=c(2,3))
plot(c0)
lines(cens,col=2) 
plot(c3,main="rates 1-> 3 , 2->3")
lines(dr,col=1,lwd=2)
lines(dr2,col=2,lwd=2)
###
plot(c1,main="rate 1->2")
lines(base1,lwd=2)
###
plot(c2,main="rate 2->1")
lines(base1,lwd=2)
 
```
