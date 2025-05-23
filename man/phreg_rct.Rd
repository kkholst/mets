% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phreg_rct.R
\name{phreg_rct}
\alias{phreg_rct}
\title{Lu-Tsiatis More Efficient Log-Rank for Randomized studies with baseline covariates}
\usage{
phreg_rct(
  formula,
  data,
  cause = 1,
  cens.code = 0,
  typesR = c("R0", "R1", "R01"),
  typesC = c("C", "dynC"),
  weights = NULL,
  augmentR0 = NULL,
  augmentR1 = NULL,
  augmentC = NULL,
  treat.model = ~+1,
  RCT = TRUE,
  treat.var = NULL,
  km = TRUE,
  level = 0.95,
  cens.model = NULL,
  estpr = 1,
  pi0 = 0.5,
  base.augment = FALSE,
  return.augmentR0 = FALSE,
  mlogit = FALSE,
  ...
)
}
\arguments{
\item{formula}{formula with 'Surv' or 'Event' outcome (see \code{coxph}) and
treatment (randomization 0/1)}

\item{data}{data frame}

\item{cause}{to use for competing risks, recurrent events data}

\item{cens.code}{to use for competing risks, recurrent events data}

\item{typesR}{augmentations used for randomization}

\item{typesC}{augmentations used for censoring}

\item{weights}{weights for score equation}

\item{augmentR0}{formula for the randomization augmentation (~age+sex)}

\item{augmentR1}{formula for the randomization augmentation (~age+sex)}

\item{augmentC}{formula for the censoring augmentation (~age+sex)}

\item{treat.model}{propensity score model, default is ~+1, assuming an RCT study}

\item{RCT}{if false will use propensity score adjustment for marginal model}

\item{treat.var}{in case of twostage randomization, this variable is 1 for
the treatment times, if start,stop then default assumes that only one
treatment at first record}

\item{km}{use Kaplan-Meier for the censoring weights (stratified on
treatment)}

\item{level}{of confidence intervals}

\item{cens.model}{default is censoring model ~strata(treatment) but any
model can be used to make censoring martingales}

\item{estpr}{estimates propensity scores}

\item{pi0}{possible fixed propensity scores for randomizations}

\item{base.augment}{TRUE to covariate augment baselines (only for R0
augmentation)}

\item{return.augmentR0}{to return augmentation data}

\item{mlogit}{if TRUE then forces use of this function for propensity
scores, default for binary treatment is glm}

\item{...}{Additional arguments to phreg function}
}
\description{
Efficient implementation of the Lu-Tsiatis improvement using baseline covariates, extended to competing risks and recurrent events. Results
almost equivalent with the speffSurv function of the speff2trial function in the survival case. A dynamic 
censoring augmentation regression is also computed to gain even more from the censoring augmentation. Furhter, we also deal with twostage
randomizations. The function was implemented to deal with recurrent events (start,stop) + cluster, and  more examples in vignette.
}
\examples{
## Lu, Tsiatis simulation
data <- mets:::simLT(0.7,100)
dfactor(data) <- Z.f~Z

out <- phreg_rct(Surv(time,status)~Z.f,data=data,augmentR0=~X,augmentC=~factor(Z):X)
summary(out)
}
\references{
Lu, Tsiatis (2008), Improving the efficiency of the log-rank
  test using auxiliary covariates, Biometrika, 679--694

Scheike et al. (2024), WIP, Two-stage randomization for recurrent events,
}
\author{
Thomas Scheike
}
