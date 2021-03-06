% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recurrent.marginal.R
\name{simRecurrent}
\alias{simRecurrent}
\alias{showfitsim}
\alias{simRecurrentGamma}
\alias{covIntH1dM1IntH2dM2}
\alias{recurrentMarginalgam}
\alias{squareintHdM}
\alias{addCums}
\title{Simulation of recurrent events data based on cumulative hazards}
\usage{
simRecurrent(
  n,
  cumhaz,
  death.cumhaz = NULL,
  cumhaz2 = NULL,
  gap.time = FALSE,
  max.recurrent = 100,
  dhaz = NULL,
  haz2 = NULL,
  dependence = 0,
  var.z = 2,
  cor.mat = NULL,
  ...
)
}
\arguments{
\item{n}{number of id's}

\item{cumhaz}{cumulative hazard of recurrent events}

\item{death.cumhaz}{cumulative hazard of death}

\item{cumhaz2}{cumulative hazard of recurrent events  of type 2}

\item{gap.time}{if true simulates gap-times with specified cumulative hazard}

\item{max.recurrent}{limits number recurrent events to 100}

\item{dhaz}{rate for death hazard if it is extended to time-range of first event}

\item{haz2}{rate of second cause  if it is extended to time-range of first event}

\item{dependence}{=0 independence, =1 all share same random effect with variance var.z
=2 random effect exp(normal) with correlation structure from cor.mat,
first random effect is z1 and shared for a possible second cause,  second random effect is for death}

\item{var.z}{variance of random effects}

\item{cor.mat}{correlation matrix for var.z variance of random effects}

\item{...}{Additional arguments to lower level funtions}
}
\description{
Simulation of recurrent events data based on cumulative hazards
}
\details{
Must give hazard of death and recurrent events.  Possible with two
event types and their dependence can be specified but the two recurrent events need
to have the same random effect,  simRecurrentII more flexible !
}
\examples{
########################################
## getting some rates to mimick 
########################################

data(base1cumhaz)
data(base4cumhaz)
data(drcumhaz)
dr <- drcumhaz
base1 <- base1cumhaz
base4 <- base4cumhaz

 ######################################################################
 ### simulating simple model that mimicks data 
 ######################################################################
 rr <- simRecurrent(5,base1,death.cumhaz=dr)
 dlist(rr,.~id,n=0)

 rr <- simRecurrent(1000,base1,death.cumhaz=dr)
 par(mfrow=c(1,3))
 showfitsim(causes=1,rr,dr,base1,base1)

######################################################################
### simulating simple model that mimicks data 
### now with two event types and second type has same rate as death rate
######################################################################

 rr <- simRecurrent(1000,base1,death.cumhaz=dr,cumhaz2=base4)
 dtable(rr,~death+status)
 par(mfrow=c(2,2))
 showfitsim(causes=2,rr,dr,base1,base4)

######################################################################
### simulating simple model 
### random effect for all causes (Z shared for death and recurrent) 
######################################################################

 rr <- simRecurrent(1000,base1,death.cumhaz=dr,dependence=1,var.gamma=0.4)
 ### marginals do fit after input after integrating out
 par(mfrow=c(2,2))
 showfitsim(causes=1,rr,dr,base1,base1)

}
\author{
Thomas Scheike
}
