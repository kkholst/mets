% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twostage.R
\name{twostageMLE}
\alias{twostageMLE}
\title{Twostage survival model fitted by pseudo MLE}
\usage{
twostageMLE(
  margsurv,
  data = parent.frame(),
  theta = NULL,
  theta.des = NULL,
  var.link = 0,
  method = "NR",
  no.opt = FALSE,
  weights = NULL,
  se.cluster = NULL,
  ...
)
}
\arguments{
\item{margsurv}{Marginal model from phreg}

\item{data}{data frame}

\item{theta}{Starting values for variance components}

\item{theta.des}{design for dependence parameters, when pairs are given this is could be a
(pairs) x (numer of parameters)  x (max number random effects) matrix}

\item{var.link}{Link function for variance  if 1 then uses exp link}

\item{method}{type of opitmizer, default is Newton-Raphson "NR"}

\item{no.opt}{to not optimize, for example to get score and iid for specific theta}

\item{weights}{cluster specific weights, but given with length equivalent to data-set, weights for score equations}

\item{se.cluster}{specifies how the influence functions are summed before squared when computing the variance. Note that the id from the marginal model is used to construct MLE, and then these scores can be summed with the se.cluster argument.}

\item{...}{arguments to be passed to  optimizer}
}
\description{
Fits Clayton-Oakes clustered  survival data
using marginals that are on Cox form in the likelihood for the dependence parameter
as in Glidden (2000). The dependence can be modelled via  a
\enumerate{
\item  Regression design on dependence parameter.
}

We allow a regression structure for the indenpendent gamma distributed
random effects  and their variances that may depend on cluster covariates. So
\deqn{
 \theta = h( z_j^T \alpha)
}
where \eqn{z} is specified by theta.des . The link function can be the exp when var.link=1
}
\examples{
data(diabetes)
dd <- phreg(Surv(time,status==1)~treat+cluster(id),diabetes)
oo <- twostageMLE(dd,data=diabetes)
summary(oo)

theta.des <- model.matrix(~-1+factor(adult),diabetes)

oo <-twostageMLE(dd,data=diabetes,theta.des=theta.des)
summary(oo)
}
\references{
Measuring early or late dependence for bivariate twin data
Scheike, Holst, Hjelmborg (2015), LIDA

Twostage modelling of additive gamma frailty models for survival data.
Scheike and Holst, working paper

Shih and Louis (1995) Inference on the association parameter in copula models for bivariate
survival data, Biometrics, (1995).

Glidden (2000), A Two-Stage estimator of the dependence
parameter for the Clayton Oakes model, LIDA, (2000).
}
\author{
Thomas Scheike
}
\keyword{survival}
