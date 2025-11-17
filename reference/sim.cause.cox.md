# Simulation of cause specific from Cox models.

Simulates data that looks like fit from cause specific Cox models.
Censor data automatically. When censoring is given in the list of causes
this will give censoring that looks like the data. Covariates are drawn
from data-set with replacement. This gives covariates like the data.

## Usage

``` r
sim.cause.cox(coxs,n,data=NULL,cens=NULL,rrc=NULL,...)
```

## Arguments

- coxs:

  list of cox models.

- n:

  number of simulations.

- data:

  to extract covariates for simulations (draws from observed
  covariates).

- cens:

  specifies censoring model, if NULL then only censoring for each cause
  at end of last event of this type. if "is.matrix" then uses
  cumulative. hazard given, if "is.scalar" then uses rate for
  exponential, and if not given then takes average rate of in simulated
  data from cox model. But censoring can also be given as a cause.

- rrc:

  possible vector of relative risk for cox-type censoring.

- ...:

  arguments for rchaz, for example entry-time

## Author

Thomas Scheike

## Examples

``` r
nsim <- 100; data(bmt)

cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=bmt)
cox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=bmt)
coxs <- list(cox1,cox2)
dd <- sim.cause.cox(coxs,nsim,data=bmt)
scox1 <- phreg(Surv(time,status==1)~strata(tcell)+platelet,data=dd)
scox2 <- phreg(Surv(time,status==2)~tcell+strata(platelet),data=dd)
cbind(cox1$coef,scox1$coef)
#>                [,1]       [,2]
#> platelet -0.5658612 -0.3001387
cbind(cox2$coef,scox2$coef)
#>            [,1]      [,2]
#> tcell 0.4153706 0.1204181
par(mfrow=c(1,2))
plot(cox1); plot(scox1,add=TRUE); 
plot(cox2); plot(scox2,add=TRUE); 

```
