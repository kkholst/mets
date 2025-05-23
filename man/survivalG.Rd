% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phreg.R
\name{survivalG}
\alias{survivalG}
\alias{survivalGtime}
\title{G-estimator for Cox and Fine-Gray model}
\usage{
survivalG(
  x,
  data,
  time = NULL,
  Avalues = c(0, 1),
  varname = NULL,
  same.data = TRUE,
  id = NULL,
  subdata = NULL
)
}
\arguments{
\item{x}{phreg or cifreg object}

\item{data}{data frame for risk averaging}

\item{time}{for estimate}

\item{Avalues}{values to compare for first covariate A}

\item{varname}{if given then averages for this variable, default is first variable}

\item{same.data}{assumes that same data is used for fitting of survival model and averaging.}

\item{id}{might be given to link to data to iid decomposition of survival data, must be coded as 1,2,..,}

\item{subdata}{rows or TRUE/FALSE to select which part of the data that is used for the G-computation. Might be treated}
}
\description{
Computes G-estimator \deqn{ \hat S(t,A=a) = n^{-1} \sum_i \hat S(t,A=a,Z_i) }
for the Cox model based on phreg og the Fine-Gray model based on the
cifreg function. Gives influence functions of these risk estimates and SE's are 
based on  these.  If first covariate is a factor then all contrast are computed, 
and if continuous then considered covariate values are given by Avalues.
}
\examples{

data(bmt); bmt$time <- bmt$time+runif(408)*0.001
bmt$event <- (bmt$cause!=0)*1
dfactor(bmt) <- tcell.f~tcell

fg1 <- cifreg(Event(time,cause)~tcell.f+platelet+age,bmt,cause=1,
              cox.prep=TRUE,propodds=NULL)
summary(survivalG(fg1,bmt,50))

ss <- phreg(Surv(time,event)~tcell.f+platelet+age,bmt) 
summary(survivalG(ss,bmt,50))

ss <- phreg(Surv(time,event)~strata(tcell.f)+platelet+age,bmt) 
summary(survivalG(ss,bmt,50))

sst <- survivalGtime(ss,bmt,n=50)
plot(sst)

fg1t <- survivalGtime(fg1,bmt,n=50)
plot(fg1t)

}
\author{
Thomas Scheike
}
