# Simulation of output from Cumulative incidence regression model

Simulates data that looks like fit from fitted cumulative incidence
model

## Usage

``` r
sim.cif(cif,n,data=NULL,Z=NULL,drawZ=TRUE,cens=NULL,rrc=NULL,
                cumstart=c(0,0),U=NULL,pU=NULL,...)
```

## Arguments

- cif:

  output form prop.odds.subdist or ccr (cmprsk), can also call
  invsubdist with with cumulative and linear predictor

- n:

  number of simulations.

- data:

  to extract covariates for simulations (draws from observed
  covariates).

- Z:

  to use these covariates for simulation rather than drawing new ones.

- drawZ:

  to random sample from Z or not

- cens:

  specifies censoring model, if "is.matrix" then uses cumulative hazard
  given, if "is.scalar" then uses rate for exponential, and if not given
  then takes average rate of in simulated data from cox model.

- rrc:

  possible vector of relative risk for cox-type censoring.

- cumstart:

  to start cumulatives at time 0 in 0.

- U:

  uniforms to use for drawing of timing for cumulative incidence.

- pU:

  uniforms to use for drawing event type (F1,F2,1-F1-F2).

- ...:

  arguments for simsubdist (for example Uniform variable for
  realizations)

## Author

Thomas Scheike

## Examples

``` r
data(bmt)

scif <-  cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,prop=NULL)
summary(scif)  
#> 
#>    n events
#>  408    161
#> 
#>  408 clusters
#> coeffients:
#>           Estimate      S.E.   dU^-1/2 P-value
#> tcell    -0.596990  0.270415  0.275780  0.0273
#> platelet -0.426702  0.180630  0.187719  0.0182
#> age       0.343724  0.080258  0.086277  0.0000
#> 
#> exp(coeffients):
#>          Estimate    2.5%  97.5%
#> tcell     0.55047 0.32401 0.9352
#> platelet  0.65266 0.45807 0.9299
#> age       1.41019 1.20493 1.6504
#> 
plot(scif)

################################################################
#  simulating several causes with specific cumulatives 
################################################################

cif1 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=1,prop=NULL)
cif2 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=2,prop=NULL)
# dd <- sim.cifsRestrict(list(cif1,cif2),200,data=bmt)
dd <- sim.cifs(list(cif1,cif2),200,data=bmt)
scif1 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=1)
scif2 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=2)
   
par(mfrow=c(1,2))   
plot(cif1); plot(scif1,add=TRUE,col=2)
plot(cif2); plot(scif2,add=TRUE,col=2)
```
