##' Restricted IPCW Mean for Censored Survival Data
##'
##' Provides a fast implementation of Inverse Probability of Censoring Weighting (IPCW) 
##' regression for a single time point. It fits the model:
##' \deqn{ E( \min(T, t) | X ) = \exp( X^T \beta) }
##' or, in the case of competing risks data:
##' \deqn{ E( I(\epsilon=1) (t - \min(T, t)) | X ) = \exp( X^T \beta) }
##' which represents the "Years Lost Due to Cause" (RMTL).
##'
##' The method solves the binomial regression IPCW response estimating equation:
##' \deqn{ X \left( \frac{\Delta(\min(T,t)) Y}{G_c(\min(T,t))} - \exp( X^T \beta) \right) = 0 }
##' where \eqn{\Delta(\min(T,t)) = I(\min(T,t) \leq C)} is the indicator of being uncensored 
##' at the time of interest.
##'
##' When the status variable is binary, the outcome is assumed to be \eqn{Y = \min(T,t)} (RMST).
##' If the status has more than two levels (competing risks), the outcome is 
##' \eqn{Y = (t - \min(T,t)) I(\text{status}=\text{cause})} (RMTL for a specific cause).
##'
##' The function supports:
##' \itemize{
##'   \item \strong{IPCW Adjustment}: Weights by the inverse of the censoring survival probability.
##'   \item \strong{Augmentation}: Can include an augmentation term (type="II" or "III") to improve efficiency 
##'     and robustness (Double Robust estimation).
##'   \item \strong{Variance Estimation}: Based on the influence function, including adjustments for 
##'     the estimation of the censoring model.
##' }
##'
##' @param formula Formula with an \code{Event} outcome (e.g., \code{Event(time, cause)}).
##' @param data Data frame containing the variables.
##' @param outcome Outcome type: \code{"rmst"} (Restricted Mean Survival Time) or 
##'   \code{"rmtl"} (Restricted Mean Time Lost).
##' @param ... Additional arguments passed to \code{binreg}, such as \code{time}, \code{cause}, 
##'   \code{cens.model}, \code{model}, \code{type}, etc.
##' @return An object of class \code{"binreg"} containing:
##'   \item{coef}{Coefficient estimates.}
##'   \item{se}{Standard errors.}
##'   \item{var}{Variance-covariance matrix.}
##'   \item{iid}{Influence function decomposition.}
##'   \item{naive.var}{Variance under known censoring model (if applicable).}
##'   \item{time}{Time point used.}
##'   \item{outcome}{Type of outcome analyzed.}
##' @author Thomas Scheike
##' @references 
##' Scheike, T. and Holst, K. K. (2024). Restricted mean time lost for survival and competing risks data using mets in R. WIP.
##' @seealso \code{\link{binreg}}, \code{\link{resmeanATE}}, \code{\link{rmstIPCW}}
##' @examples
##' data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001
##' 
##' # E( min(T;t) | X ) = exp( a+b X) with IPCW estimation 
##' out <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age, bmt,
##'                 time=50, cens.model=~strata(platelet), model="exp")
##' summary(out)
##' 
##' \donttest{ ## Reduce Ex.Timings
##' ## Weighted GLM version (RMST)
##' out2 <- logitIPCW(Event(time,cause!=0)~tcell+platelet+age, bmt,
##'             time=50, cens.model=~strata(platelet), model="exp", outcome="rmst")
##' summary(out2)
##' 
##' ### Time-lost (RMTL)
##' outtl <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age, bmt,
##'                 time=50, cens.model=~strata(platelet), model="exp", outcome="rmtl")
##' summary(outtl)
##' 
##' ### Same as Kaplan-Meier for full censoring model 
##' bmt$int <- with(bmt, strata(tcell, platelet))
##' out <- resmeanIPCW(Event(time,cause!=0)~-1+int, bmt, time=30,
##'                              cens.model=~strata(platelet, tcell), model="lin")
##' estimate(out)
##' out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet), data=bmt)
##' rm1 <- resmean_phreg(out1, times=30)
##' summary(rm1)
##' 
##' ### Years lost regression
##' outl <- resmeanIPCW(Event(time,cause!=0)~-1+int, bmt, time=30, outcome="years-lost",
##'                              cens.model=~strata(platelet, tcell), model="lin")
##' estimate(outl)
##' 
##' ## Competing risks years-lost for cause 1  
##' out <- resmeanIPCW(Event(time,cause)~-1+int, bmt, time=30, cause=1,
##'                             cens.model=~strata(platelet, tcell), model="lin")
##' estimate(out)
##' ## Same as integrated cumulative incidence 
##' rmc1 <- cif_yearslost(Event(time,cause)~strata(tcell,platelet), data=bmt, times=30)
##' summary(rmc1)
##' }
##' @aliases rmstIPCW 
##' @export
resmeanIPCW  <- function(formula,data,outcome=c("rmst","rmtl"),...)
{# {{{
   out <- binreg(formula,data,outcome=outcome[1],...)
   return(out)
}# }}}

##' @export
rmstIPCW <- function(formula,data,outcome=c("rmst","rmtl"),...)
{# {{{
   out <- binreg(formula,data,outcome=outcome[1],...)
   return(out)
}# }}}

preprrm <- function(cs,ss,X,times,data,model="exp") 
{# {{{

ttimes <- ss$cumhaz[,1]
ttimes <- ttimes[(ttimes< times)]
ttimes <- c(ttimes,times)
## cens model 
pcs <- predict(cs,data,se=FALSE,times=ttimes)
pcst0 <- c(tail(t(pcs$surv),n=1))
## survival model 
ps <- predict(ss,data,se=FALSE,times=ttimes)
###dLamt <- apply(cbind(0,ps$cumhaz),1,diff)
pst0 <- c(tail(t(ps$surv),n=1))
dtime <- diff(c(0,ttimes))
RRM <-  apply(dtime * t(ps$surv),2,cumsum)
RRMt0 <- c(tail(RRM,n=1))
dF <- -apply(cbind(1,ps$surv),1,diff)
I <- apply(ttimes^2/t(pcs$surv)*dF,2,sum)+times^2*pst0/pcst0 -RRMt0*apply(ttimes/t(pcs$surv)*dF,2,sum)-RRMt0*times*pst0/pcst0 
varY <- apply(ttimes^2*dF,2,sum)+times^2*pst0 - RRMt0^2

ctimes <- cs$cumhaz[,1]
ctimes <- ctimes[(ctimes< times)]
## cens model 
pcs <- predict(cs,data,se=FALSE,times=ctimes)
pcst0 <- c(tail(t(pcs$surv),n=1))
dLamc <- apply(cbind(0,pcs$cumhaz),1,diff)
## survival model 
ps <- predict(ss,data,se=FALSE,times=ctimes)
pst0 <- c(tail(t(ps$surv),n=1))
Fsst0 <- 1-pst0
dtime <- diff(c(0,ctimes))
###dS <- -apply(cbind(1,ps$surv),1,diff)
RRM <-  apply(dtime * t(ps$surv),2,cumsum)
###
RRMt0 <- c(tail(RRM,n=1))
dF <- dRRM <- t(RRMt0-t(RRM)) + ctimes*t(ps$surv)
II <- apply((dF^2/t(ps$surv*pcs$surv))*dLamc,2,sum)-RRMt0*apply((dF/t(pcs$surv))*dLamc,2,sum) 
if (model=="exp") DbetaF <- RRMt0 else DbetaF <- 1
h <- DbetaF/(I-II) 
###varYII <- 2*apply(ctimes*dtime*ps$surv,1,sum)-RRMt0^2

cmtimes <- matrix(ctimes,nrow(dLamc),ncol(dLamc))
tttimes <- matrix(pmin(data[,"time"],times),nrow(dLamc),ncol(dLamc),byrow=TRUE)
Augmentf <- (dF/t(ps$surv*pcs$surv))
AugmentC <- apply(Augmentf*dLamc*(cmtimes<=tttimes),2,sum)
###
n <- nrow(data)
AdN <- rep(0,n)
jc <- (1:nrow(data))[cs$ord][cs$jumps]
jc <-jc[1:length(ctimes)] 
AdN[jc] <- mdi(Augmentf,1:length(ctimes),jc)
Mc <- AdN- AugmentC

ph <- 1
if (model=="exp") ph <- RRMt0
Xaugment <- apply(X*Mc*ph,2,sum)
augment <- apply(X*Mc*h,2,sum)
hh <- (DbetaF/varY)
Faugment <- apply(X*hh*Mc,2,sum)
return(list(Mc=Mc,Xaugment=Xaugment,Faugment=Faugment,hXaugment=augment,h=h,hh=hh,varY=varY,RRMt0=RRMt0))
}# }}}


##' Average Treatment Effect for Restricted Mean Time
##'
##' Estimates the Average Treatment Effect (ATE) for Restricted Mean Survival Time (RMST) 
##' or Restricted Mean Time Lost (RMTL) in censored competing risks data using IPCW.
##' 
##' Under standard causal assumptions (Consistency, Ignorability, Positivity), the ATE is 
##' estimated as \eqn{E(Y(1) - Y(0))}, where \eqn{Y(a)} is the potential outcome under treatment \eqn{a}.
##' The method uses double robust estimating equations that are IPCW-adjusted for censoring.
##'
##' The first covariate in the formula must be the treatment effect (a factor). If the factor 
##' has more than two levels, multinomial logistic regression (mlogit) is used for propensity 
##' score modeling.
##'
##' @param formula Formula with an \code{Event} outcome. The first covariate must be the treatment factor.
##' @param data Data frame.
##' @param model Link function: \code{"exp"} (exponential) or \code{"lin"} (identity).
##' @param outcome Outcome type: \code{"rmst"} or \code{"rmtl"}.
##' @param ... Additional arguments passed to \code{binregATE}, such as \code{time}, \code{treat.model}, 
##'   \code{augmentR0}, \code{augmentC}, etc.
##' @return An object of class \code{"binregATE"} containing:
##'   \item{riskG}{Simple IPCW estimator results.}
##'   \item{riskDR}{Double Robust estimator results.}
##'   \item{riskG.iid, riskDR.iid}{Influence functions.}
##'   \item{coef}{Treatment effect estimates.}
##'   \item{se}{Standard errors.}
##' @author Thomas Scheike
##' @references 
##' Scheike, T. and Holst, K. K. (2024). Restricted mean time lost for survival and competing risks data using mets in R. WIP.
##' 
##' Scheike, T. and Tanaka, S. (2025). Restricted mean time lost ratio regression: Percentage of restricted mean time lost due to specific cause. WIP.
##' @seealso \code{\link{binregATE}}, \code{\link{ratioATE}}
##' @examples
##' data(bmt); bmt$event <- bmt$cause!=0; dfactor(bmt) <- tcell~tcell
##' 
##' out <- resmeanATE(Event(time,event)~tcell+platelet, data=bmt, time=40, 
##'                   treat.model=tcell~platelet, outcome="rmtl")
##' summary(out)
##' 
##' out1 <- resmeanATE(Event(time,cause)~tcell+platelet, data=bmt, cause=1, time=40,
##'                    treat.model=tcell~platelet, outcome="rmtl")
##' summary(out1)
##' 
##' @aliases rmstATE 
##' @export
resmeanATE <- function(formula,data,model="exp",outcome=c("rmst","rmtl"),...)
{# {{{
out <- 	binregATE(formula,data,outcome=outcome[1],model=model,...) 
return(out)
}# }}}

##' @export
rmstATE <- function(formula,data,model="exp",outcome=c("rmst","rmtl"),...)
{# {{{
out <- 	binregATE(formula,data,outcome=outcome[1],model=model,...) 
return(out)
}# }}}

##' Ratio of Average Treatment Effects
##'
##' Computes the ratio of two Average Treatment Effects (ATEs), typically comparing the 
##' ATE for a specific cause (e.g., RMTL due to cause 1) to the ATE for the total RMTL.
##' 
##' The function transforms the estimates (e.g., using log) to compute the ratio and 
##' its standard error using the delta method and the joint influence functions.
##'
##' @param rmtl Object containing the ATE for the total RMTL (from \code{resmeanATE}).
##' @param rmtl1 Object containing the ATE for the specific cause RMTL (from \code{resmeanATE}).
##' @param h Transformation function (e.g., \code{log}) applied to the ratio. Default is identity.
##' @param null Value under the null hypothesis for the ratio (default 1).
##' @return A list containing:
##'   \item{ratioG}{Ratio based on the simple IPCW estimator.}
##'   \item{ratioDR}{Ratio based on the double robust estimator.}
##' @author Thomas Scheike
##' @seealso \code{\link{resmeanATE}}
##' @examples
##' data(bmt); bmt$event <- bmt$cause!=0; dfactor(bmt) <- tcell~tcell
##' out <- resmeanATE(Event(time,event)~tcell+platelet, data=bmt, time=40, outcome="rmtl")
##' out1 <- resmeanATE(Event(time,cause)~tcell+platelet, data=bmt, cause=1, time=40, outcome="rmtl")
##' ratioATE(out, out1, h=log)
##' @export
ratioATE <- function(rmtl, rmtl1, h=NULL, null=1)
##' @export
ratioATE <- function(rmtl,rmtl1,h=NULL,null=1) { ## {{{

if (is.null(h)) h <- function(x) x

coefG <- c(rmtl$riskG,rmtl1$riskG)
iidG <- cbind( rmtl$riskG.iid, rmtl1$riskG.iid)
covG <- crossprod(iidG)
###
coefDR <- c(rmtl$riskDR,rmtl1$riskDR)
iidDR <- cbind(rmtl$riskDR.iid, rmtl1$riskDR.iid)
covDR <- crossprod(iidDR)

ratioG <- estimate(coef=coefG,vcov=covG,f=function(p) h(c(p[3:4]/p[1:2],p[1]*p[4]/(p[2]*p[3]),p[2]*p[3]/(p[1]*p[4]))),null=null)
ratioDR <- estimate(coef=coefDR,vcov=covDR,f=function(p) h(c(p[3:4]/p[1:2],p[1]*p[4]/(p[2]*p[3]),p[2]*p[3]/(p[1]*p[4]))),null=null)
out <- list(ratioG=ratioG,ratioDR=ratioDR)
return(out)

} ## }}}

