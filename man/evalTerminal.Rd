% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA-estimands.R
\name{evalTerminal}
\alias{evalTerminal}
\title{Evaluates piece constant covariates at min(D,t) where D is a terminal event}
\usage{
evalTerminal(
  formula,
  data = data,
  death.code = 2,
  time = NULL,
  marks = NULL,
  mark.codes = NULL
)
}
\arguments{
\item{formula}{formula with 'Event' outcome and X to evaluate at min(D,t)}

\item{data}{data frame}

\item{death.code}{codes for death (terminating event, 2 default)}

\item{time}{for evaluation}

\item{marks}{for terminal events to add marks*I(D <=t ,epsilon "in" mark.codes)  to X(min(D,t))}

\item{mark.codes}{gives death codes for which to add mark value}
}
\description{
returns X(min(D,t)) and min(D,t) and their ratio. for censored observation 0. 
to use with the IPCW models implemented.
}
\author{
Thomas Scheike
}
