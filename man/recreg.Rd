% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recreg.R
\name{recreg}
\alias{recreg}
\alias{strataAugment}
\alias{scalecumhaz}
\alias{GLprediid}
\alias{recregIPCW}
\alias{twostageREC}
\alias{IIDrecreg}
\alias{predicttime}
\title{Recurrent events regression with terminal event}
\usage{
recreg(
  formula,
  data,
  cause = 1,
  death.code = 2,
  cens.code = 0,
  cens.model = ~1,
  weights = NULL,
  offset = NULL,
  Gc = NULL,
  wcomp = NULL,
  marks = NULL,
  augmentation.type = c("lindyn.augment", "lin.augment"),
  ...
)
}
\arguments{
\item{formula}{formula with 'Event' outcome}

\item{data}{data frame}

\item{cause}{of interest (1 default)}

\item{death.code}{codes for death (terminating event, 2 default)}

\item{cens.code}{code of censoring (0 default)}

\item{cens.model}{for stratified Cox model without covariates}

\item{weights}{weights for score equations}

\item{offset}{offsets for model}

\item{Gc}{censoring weights for time argument, default is to calculate these with a Kaplan-Meier estimator, should then give G_c(T_i-)}

\item{wcomp}{weights for composite outcome, so when cause=c(1,3), we might have wcomp=c(1,2).}

\item{marks}{a mark value can be specified, this is vector from the data-frame where the mark value can be found at all events}

\item{augmentation.type}{of augmentation when augmentation model is given}

\item{...}{Additional arguments to lower level funtions}
}
\description{
Fits Ghosh-Lin IPCW Cox-type model
}
\details{
For Cox type model :
\deqn{
E(dN_1(t)|X) = \mu_0(t)dt exp(X^T \beta)
}
by solving Cox-type IPCW weighted score equations 
\deqn{
 \int (Z - E(t)) w(t) dN_1(t) 
}
where \deqn{w(t) = G(t) (I(T_i \wedge t < C_i)/G_c(T_i \wedge t))} and
\deqn{E(t) = S_1(t)/S_0(t)} and \deqn{S_j(t) = \sum X_i^j w_i(t) \exp(X_i^T \beta)}.


The iid decomposition of the beta's are on the form
\deqn{
\int (Z - E ) w(t) dM_1 + \int q(s)/p(s) dM_c
}
and returned as iid.

Events, deaths and censorings are specified via stop start structure and the Event call, that via a status vector 
and cause (code), censoring-codes (cens.code) and death-codes (death.code) indentifies these. See example and vignette.
}
\examples{
## data with no ties
library(mets)
data(hfactioncpx12)
hf <- hfactioncpx12
hf$x <- as.numeric(hf$treatment) 
dd <- data.frame(treatment=levels(hf$treatment),id=1)

gl <- recreg(Event(entry,time,status)~treatment+cluster(id),data=hf,cause=1,death.code=2)
summary(gl)
pgl <- predict(gl,dd,se=1); plot(pgl,se=1)

## censoring stratified after treatment 
gls <- recreg(Event(entry,time,status)~treatment+cluster(id),data=hf,
cause=1,death.code=2,cens.model=~strata(treatment))
summary(gls)

## IPCW at 2 years 
ll2 <- recregIPCW(Event(entry,time,status)~treatment+cluster(id),data=hf,
cause=1,death.code=2,time=2,cens.model=~strata(treatment))
summary(ll2)

}
\author{
Thomas Scheike
}
