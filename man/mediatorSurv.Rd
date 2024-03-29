% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediator-survival.R
\name{mediatorSurv}
\alias{mediatorSurv}
\alias{BootmediatorSurv}
\title{Mediation analysis in survival context}
\usage{
mediatorSurv(
  survmodel,
  weightmodel,
  data = data,
  wdata = wdata,
  id = "id",
  silent = TRUE,
  ...
)
}
\arguments{
\item{survmodel}{with mediation model (binreg, aalenMets, phreg)}

\item{weightmodel}{mediation model}

\item{data}{for computations}

\item{wdata}{weighted data expansion for computations}

\item{id}{name of id variable, important for SE computations}

\item{silent}{to be silent}

\item{...}{Additional arguments to survival model}
}
\description{
Mediation analysis in survival context  with robust standard errors taking the weights into account
via influence function computations. Mediator and exposure must be factors.  This is based on numerical
derivative wrt parameters for weighting.  See vignette for more examples.
}
\examples{

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
## not run bootstrap (to save time)
## bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=500)

}
\author{
Thomas Scheike
}
