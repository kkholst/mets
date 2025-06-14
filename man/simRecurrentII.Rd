% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recurrent.marginal.R
\name{simRecurrentII}
\alias{simRecurrentII}
\alias{simRecurrent}
\alias{showfitsim}
\alias{covIntH1dM1IntH2dM2}
\alias{squareintHdM}
\alias{simRecurrentList}
\alias{showfitsimList}
\title{Simulation of recurrent events data based on cumulative hazards with two types of  recurrent events}
\usage{
simRecurrentII(
  n,
  cumhaz,
  cumhaz2,
  death.cumhaz = NULL,
  r1 = NULL,
  r2 = NULL,
  rd = NULL,
  rc = NULL,
  dependence = 0,
  var.z = 1,
  cor.mat = NULL,
  cens = NULL,
  gap.time = FALSE,
  max.recurrent = 100,
  ...
)
}
\arguments{
\item{n}{number of id's}

\item{cumhaz}{cumulative hazard of recurrent events}

\item{cumhaz2}{cumulative hazard of recurrent events  of type 2}

\item{death.cumhaz}{cumulative hazard of death}

\item{r1}{potential relative risk adjustment of rate}

\item{r2}{potential relative risk adjustment of rate}

\item{rd}{potential relative risk adjustment of rate}

\item{rc}{potential relative risk adjustment of rate}

\item{dependence}{0:independence; 1:all share same random effect with variance var.z; 2:random effect exp(normal) with correlation structure from cor.mat; 3:additive gamma distributed random effects, z1= (z11+ z12)/2 such that mean is 1 , z2= (z11^cor.mat(1,2)+ z13)/2, z3= (z12^(cor.mat(2,3)+z13^cor.mat(1,3))/2, with z11 z12 z13 are gamma with mean and variance 1 , first random effect is z1 and for N1 second random effect is z2 and for N2 third random effect is for death}

\item{var.z}{variance of random effects}

\item{cor.mat}{correlation matrix for var.z variance of random effects}

\item{cens}{rate of censoring exponential distribution}

\item{gap.time}{if true simulates gap-times with specified cumulative hazard}

\item{max.recurrent}{limits number recurrent events to 100}

\item{...}{Additional arguments to lower level funtions}
}
\description{
Simulation of recurrent events data based on cumulative hazards
}
\details{
Must give hazard of death and two recurrent events.  Possible with two
event types and their dependence can be specified but the two recurrent events need
to share random effect. Based on drawing the from cumhaz and cumhaz2 and 
taking the first event rather
the cumulative and then distributing it out. Key advantage of this is that 
there is  more flexibility wrt random effects
}
\examples{
########################################
## getting some rates to mimick 
########################################
data(CPH_HPN_CRBSI)
dr <- CPH_HPN_CRBSI$terminal
base1 <- CPH_HPN_CRBSI$crbsi 
base4 <- CPH_HPN_CRBSI$mechanical

######################################################################
### simulating simple model that mimicks data 
######################################################################
rr <- simRecurrent(5,base1)
dlist(rr,.~id,n=0)
rr <- simRecurrent(5,base1,death.cumhaz=dr)
dlist(rr,.~id,n=0)

rr <- simRecurrent(100,base1,death.cumhaz=dr)
par(mfrow=c(1,3))
showfitsim(causes=1,rr,dr,base1,base1)
######################################################################
### simulating simple model 
### random effect for all causes (Z shared for death and recurrent) 
######################################################################
rr <- simRecurrent(100,base1,death.cumhaz=dr,dependence=1,var.gamma=0.4)
dtable(rr,~death+status)

######################################################################
### simulating simple model that mimicks data 
### now with two event types and second type has same rate as death rate
######################################################################
set.seed(100)
rr <- simRecurrentII(100,base1,base4,death.cumhaz=dr)
dtable(rr,~death+status)
par(mfrow=c(2,2))
showfitsim(causes=2,rr,dr,base1,base4)

## general case, 3 recurrent events and 2 causes of death
set.seed(100)
cumhaz <- list(base1,base1,base4)
drl <- list(dr,base4)
rr <- simRecurrentList(100,cumhaz,death.cumhaz=drl,dependence=0)
dtable(rr,~death+status)
showfitsimList(rr,cumhaz,drl) 

}
\author{
Thomas Scheike
}
