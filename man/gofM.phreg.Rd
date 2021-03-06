% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gof-phreg.R
\name{gofM.phreg}
\alias{gofM.phreg}
\title{GOF for Cox covariates in  PH regression}
\usage{
gofM.phreg(
  formula,
  data,
  offset = NULL,
  weights = NULL,
  modelmatrix = NULL,
  n.sim = 1000,
  silent = 1,
  ...
)
}
\arguments{
\item{formula}{formula for cox regression}

\item{data}{data for model}

\item{offset}{offset}

\item{weights}{weights}

\item{modelmatrix}{matrix for cumulating residuals}

\item{n.sim}{number of simulations for score processes}

\item{silent}{to keep it absolutely silent, otherwise timing estimate will be prduced for longer jobs.}

\item{...}{Additional arguments to lower level funtions}
}
\description{
Cumulative residuals after model matrix for Cox PH regression
p-values based on Lin, Wei, Ying resampling.
}
\details{
That is, computes 
\deqn{
 U(t) = \int_0^t M^t d \hat M 
}
and resamples its asymptotic distribution. 

This will show if the residuals are consistent with the model. Typically,
M will be a design matrix for the continous covariates that gives for example
the quartiles, and then the plot will show if for the different quartiles of the covariate the risk
prediction is consistent over time  (time x covariate interaction).
}
\examples{
library(mets)
data(TRACE)
set.seed(1)
TRACEsam <- blocksample(TRACE,idvar="id",replace=FALSE,100)

dcut(TRACEsam)  <- ~. 
mm <- model.matrix(~-1+factor(wmicat.4),data=TRACEsam)
m1 <- gofM.phreg(Surv(time,status==9)~vf+chf+wmi,data=TRACEsam,modelmatrix=mm)
summary(m1)
if (interactive()) {
par(mfrow=c(2,2))
plot(m1)
}

m1 <- gofM.phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACEsam,modelmatrix=mm) 
summary(m1)

## cumulative sums in covariates, via design matrix mm 
mm <- cumContr(TRACEsam$wmi,breaks=10,equi=TRUE)
m1 <- gofM.phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACEsam,
		  modelmatrix=mm,silent=0)
summary(m1)

}
\author{
Thomas Scheike and Klaus K. Holst
}
