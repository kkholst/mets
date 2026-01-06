# Aalen frailty model

Additive hazards model with (gamma) frailty

## Usage

``` r
aalenfrailty(time, status, X, id, theta, B = NULL, ...)
```

## Arguments

- time:

  Time variable

- status:

  Status variable (0,1)

- X:

  Covariate design matrix

- id:

  cluster variable

- theta:

  list of thetas (returns score evaluated here), or starting point for
  optimization (defaults to magic number 0.1)

- B:

  (optional) Cumulative coefficients (update theta by fixing B)

- ...:

  Additional arguments to lower level functions

## Value

Parameter estimates

## Details

Aalen frailty model

## Author

Klaus K. Holst

## Examples

``` r
library("timereg")
#> Loading required package: survival
#> 
#> Attaching package: ‘timereg’
#> The following objects are masked from ‘package:mets’:
#> 
#>     Event, event.split, kmplot, plotConfregion
dd <- simAalenFrailty(5000)
f <- ~1##+x
X <- model.matrix(f,dd) ## design matrix for non-parametric terms
system.time(out<-timereg::aalen(update(f,Surv(time,status)~.),dd,n.sim=0,robust=0))
#>    user  system elapsed 
#>    0.01    0.00    0.01 
dix <- which(dd$status==1)
t1 <- system.time(bb <- .Call("Bhat",as.integer(dd$status),
                              X,0.2,as.integer(dd$id),NULL,NULL,
                              PACKAGE="mets"))
spec <- 1
##plot(out,spec=spec)
## plot(dd$time[dix],bb$B2[,spec],col="red",type="s",
##      ylim=c(0,max(dd$time)*c(beta0,beta)[spec]))
## abline(a=0,b=c(beta0,beta)[spec])
##'

if (FALSE) { # \dontrun{
thetas <- seq(0.1,2,length.out=10)
Us <- unlist(aalenfrailty(dd$time,dd$status,X,dd$id,as.list(thetas)))
##plot(thetas,Us,type="l",ylim=c(-.5,1)); abline(h=0,lty=2); abline(v=theta,lty=2)
op <- aalenfrailty(dd$time,dd$status,X,dd$id)
op
} # }
```
