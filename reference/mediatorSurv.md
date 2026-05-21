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
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
#>    n events
#>  800    366
#> 
#>  400 clusters
#> coeffients:
#>              Estimate   Std.Err      2.5%     97.5% P-value
#> (Intercept) -0.354894  0.163603 -0.675550 -0.034239  0.0301
#> dnr.f01     -0.075193  0.253862 -0.572754  0.422368  0.7671
#> dnr.f11      0.256898  0.082388  0.095421  0.418376  0.0018
#> preauto      0.499033  0.232544  0.043255  0.954812  0.0319
#> ttt24        0.269618  0.261502 -0.242916  0.782153  0.3025
#> 
#> exp(coeffients):
#>             Estimate    2.5%  97.5%
#> (Intercept)  0.70125 0.50888 0.9663
#> dnr.f01      0.92756 0.56397 1.5256
#> dnr.f11      1.29291 1.10012 1.5195
#> preauto      1.64713 1.04420 2.5982
#> ttt24        1.30946 0.78434 2.1862
#> 
#> 
## not run bootstrap (to save time)
## bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=500)
```
