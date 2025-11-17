# Mediation analysis in survival context

Mediation analysis in survival context with robust standard errors
taking the weights into account via influence function computations.
Mediator and exposure must be factors. This is based on numerical
derivative wrt parameters for weighting. See vignette for more examples.

## Usage

``` r
mediatorSurv(
  survmodel,
  weightmodel,
  data = data,
  wdata = wdata,
  id = "id",
  silent = TRUE,
  ...
)
```

## Arguments

- survmodel:

  with mediation model (binreg, aalenMets, phreg)

- weightmodel:

  mediation model

- data:

  for computations

- wdata:

  weighted data expansion for computations

- id:

  name of id variable, important for SE computations

- silent:

  to be silent

- ...:

  Additional arguments to survival model

## Author

Thomas Scheike

## Examples

``` r
library(mets)
n <- 400
dat <- kumarsimRCT(n,rho1=0.5,rho2=0.5,rct=2,censpar=c(0,0,0,0),
          beta = c(-0.67, 0.59, 0.55, 0.25, 0.98, 0.18, 0.45, 0.31),
    treatmodel = c(-0.18, 0.56, 0.56, 0.54),restrict=1)
dfactor(dat) <- dnr.f~dnr
dfactor(dat) <- gp.f~gp
drename(dat) <- ttt24~"ttt24*"
dat$id <- 1:n
dat$ftime <- 1

weightmodel <- fit <- glm(gp.f~dnr.f+preauto+ttt24,data=dat,family=binomial)
wdata <- medweight(fit,data=dat)

### fitting models with and without mediator
aaMss2 <- binreg(Event(time,status)~gp+dnr+preauto+ttt24+cluster(id),data=dat,time=50,cause=2)
aaMss22 <- binreg(Event(time,status)~dnr+preauto+ttt24+cluster(id),data=dat,time=50,cause=2)

### estimating direct and indirect effects (under strong strong assumptions) 
aaMss <- binreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),
                data=wdata,time=50,weights=wdata$weights,cause=2)
## to compute standard errors , requires numDeriv
library(numDeriv)
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  800    364
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.277563  0.165585 -0.602103  0.046978  0.0937
#> dnr.f01     -0.013552  0.264909 -0.532765  0.505661  0.9592
#> dnr.f11      0.243449  0.081661  0.083396  0.403502  0.0029
#> preauto      0.317339  0.234407 -0.142091  0.776768  0.1758
#> ttt24        0.331677  0.264380 -0.186498  0.849853  0.2096
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.75763 0.54766 1.0481
#> dnr.f01      0.98654 0.58698 1.6581
#> dnr.f11      1.27564 1.08697 1.4971
#> preauto      1.37347 0.86754 2.1744
#> ttt24        1.39330 0.82986 2.3393
#> 
#> 
## not run bootstrap (to save time)
## bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=500)
```
