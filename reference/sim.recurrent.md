# Simulation of two-stage recurrent events data based on Cox/Cox or Cox/Ghosh-Lin structure

Simulation of two-stage recurrent events data based on Cox/Cox or
Cox/Ghosh-Lin structure

## Usage

``` r
sim.recurrent(cox1,coxd=NULL,coxc=NULL,n=1,data=NULL,
type=c("default","cox-cox","gl-cox"),id="id",
varz=1,share=1,cens=0.001,scale1=1,scaled=1,dependence=NULL,...)
```

## Arguments

- cox1:

  cox/ghosh-lin for recurrent events

- coxd:

  cox for terminal event (phreg)

- coxc:

  possible cox for censoring (phreg)

- n:

  number of id's

- data:

  on which the models are fitted (to draw covariates)

- type:

  to specify type of simulation, if not default

- id:

  name of id variable

- varz:

  dependence frailty

- share:

  to fit patly shared random effects model

- cens:

  censoring rate for exponential censoring

- scale1:

  to scale baseline of recurrent events model

- scaled:

  to scale baseline of terminal event

- dependence:

  if dependence different from NULL, then uses simRecurrentList based on
  models given

- ...:

  Additional arguments to simGLcox, nmin, nmax regulates linear
  approximation grid

## Details

Must specify two phreg objects, or a phreg and a recreg object, then
simulates data from two-stage model

## References

Scheike (2024), Twostage recurrent events models, under review.

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
hf <- hfactioncpx12
hf$x <- as.numeric(hf$treatment)
n <- 100
xr <- phreg(Surv(entry,time,status==1)~x+cluster(id),data=hf)
dr <- phreg(Surv(entry,time,status==2)~x+cluster(id),data=hf)
simcoxcox <- sim.recurrent(xr,dr,n=n,data=hf)
recGL <- recreg(Event(entry,time,status)~x+cluster(id),hf,death.code=2)
simglcox <- sim.recurrent(recGL,dr,n=n,data=hf)
```
