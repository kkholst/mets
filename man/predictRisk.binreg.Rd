% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/riskRegression.R
\name{predictRisk.binreg}
\alias{predictRisk.binreg}
\alias{predictRisk.phreg}
\alias{predictRisk.cifreg}
\alias{predictRisk.cifregFG}
\alias{predictRisk.recreg}
\title{Risk predictions to work with riskRegression package}
\usage{
\method{predictRisk}{binreg}(object, newdata, cause, times = NULL, ...)
}
\arguments{
\item{object}{phreg/binreg/cifreg object}

\item{newdata}{data.frame on which to make new predictions}

\item{cause}{cause (cif) to predict}

\item{times}{times for predictions}

\item{...}{additional arguments to lower level functions}
}
\description{
Risk predictions to work with riskRegression package
}
