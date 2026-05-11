
##' Marginal mean estimation for recurrent events with a terminal event
##'
##' Estimates the marginal mean number of recurrent events over time in the
##' presence of a competing terminal event (e.g. death), using the nonparametric
##' estimator of Ghosh and Lin (2000). Two proportional hazards models are fitted
##' internally—one for the recurrent event rate and one for the terminal
##' event—and combined to form the estimator
##' \deqn{\mu(t) = \int_0^t S(u-)\,dR(u),}
##' where \eqn{S(u)} is the marginal survival probability at the baseline covariate
##' level and \eqn{dR(u)} is the baseline recurrent event rate among survivors.
##' Robust (sandwich) standard errors are computed via the influence-function
##' approach of Ghosh and Lin (2000).
##'
##' Jump times must be unique within each stratum. If ties are present, use
##' \code{\link{tie_breaker}} to resolve them before calling this function.
##'
##' @param formula A formula with an \code{\link{Event}} response on the left-hand
##'   side, specifying entry time, exit time, and event status. The right-hand side
##'   may include \code{cluster()} to identify subjects and \code{strata()} for a
##'   stratified analysis. A \code{cluster()} term is required.
##' @param data A data frame containing all variables in \code{formula}.
##' @param cause Integer code(s) for the recurrent event of interest. Default is
##'   \code{1}.
##' @param death.code Integer code(s) for the terminal event. Default is \code{2}.
##' @param test Logical. If \code{TRUE}, a logrank-type test comparing strata is
##'   computed and stored as an attribute of the result. Default is \code{FALSE}.
##' @param ... Further arguments passed to \code{\link{phreg}}.
##'
##' @return An object of class \code{"recurrent"} with the following components:
##'   \item{mu}{Estimated marginal mean \eqn{\mu(t)} at each jump time.}
##'   \item{se.mu}{Robust standard error of \code{mu}.}
##'   \item{times}{Jump times at which estimates are computed.}
##'   \item{St}{Marginal survival estimate \eqn{S(t)} at each jump time.}
##'   \item{cumhaz}{Two-column matrix of \code{(time, mu)}, suitable for plotting.}
##'   \item{se.cumhaz}{Two-column matrix of \code{(time, se.mu)}.}
##'   The object carries three attributes: \code{"logrank"} (the test result when
##'   \code{test = TRUE}, otherwise \code{NULL}), \code{"cause"}, and
##'   \code{"death.code"}.
##'
##' @author Thomas Scheike
##'
##' @seealso \code{\link{test_logrankRecurrent}}, \code{\link{tie_breaker}},
##'   \code{\link{prob_exceed_recurrent}}
##'
##' @references
##' Cook, R. J. and Lawless, J. F. (1997). Marginal analysis of recurrent events
##' and a terminating event. \emph{Statistics in Medicine}, 16, 911--924.
##'
##' Ghosh, D. and Lin, D. Y. (2000). Nonparametric analysis of recurrent events
##' and death. \emph{Biometrics}, 56, 554--562.
##'
##' @examples
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##' hf$x <- as.numeric(hf$treatment)
##'
##' ## Fit nonparametric baseline models for recurrent events and death
##' xr <- phreg(Surv(entry, time, status == 1) ~ cluster(id), data = hf)
##' dr <- phreg(Surv(entry, time, status == 2) ~ cluster(id), data = hf)
##'
##' par(mfrow = c(1, 3))
##' plot(dr, se = TRUE); title(main = "Death")
##' plot(xr, se = TRUE); title(main = "Recurrent events")
##'
##' ## Compare naive and robust standard errors for the recurrent event rate
##' rxr <- robust_phreg(xr, fixbeta = 1)
##' plot(rxr, se = TRUE, robust = TRUE, add = TRUE, col = 4)
##'
##' ## Marginal mean via formula interface
##' outN <- recurrent_marginal(Event(entry, time, status) ~ cluster(id),
##'                            data = hf, cause = 1, death.code = 2)
##' plot(outN, se = TRUE, col = 2, add = TRUE)
##' summary(outN, times = 1:5)
##'
##' ## Stratified analysis with logrank test
##' out <- recurrent_marginal(Event(entry, time, status) ~ strata(treatment) + cluster(id),
##'                           data = hf, cause = 1, death.code = 2, test = TRUE)
##' plot(out, se = TRUE, ylab = "Marginal mean", col = 1:2)
##' attr(out, "logrank")
##' summary(out, times = 1:5)
##'
##' ## Influence-function (iid) decomposition at a fixed time point
##' head(iid(outN, time = 3))
##'
##' @aliases recurrent_marginalAIPCW recurrentMarginalPhreg iidRecurrent
##' @export
recurrent_marginal <- function(formula,data,cause=1,...,death.code=2,test=FALSE)
{# {{{
  if (missing(formula)) { # Fall-back to recurrentMarginalPhreg for backward compatibility <= 1.3.5
    return(recurrentMarginalPhreg(...))
  }
  if (inherits(formula,"phreg")) { # Fall-back to recurrentMarginalPhreg
    if (inherits(data, "phreg")) {
        return(recurrentMarginalPhreg(formula, data, ...))
    } else {
        stop("You can give 2 phreg objects \n")
    }
  }

  ## {{{
  cl <- match.call()
    m <- match.call(expand.dots = TRUE)[1:3]
    des <- proc_design(
        formula,
        data = data,
        specials = c("offset", "weights", "cluster","strata","marks"),
        intercept = FALSE
    )
    Y <- des$y
    if (!inherits(Y, c("Event", "Surv"))) {
        stop("Expected a 'Surv' or 'Event'-object")
    }
    if (ncol(Y) == 2) {
        exit <- Y[, 1]
        entry <- rep(0, nrow(Y))
        status <- Y[, 2]
    } else {
        entry <- Y[, 1]
        exit <- Y[, 2]
        status <- Y[, 3]
    }
    X <- des$x
    strata <- des$strata
    specials = c("offset", "weights", "cluster","strata")
    Terms <- terms(formula, specials, data = data)
    ts <- survival::untangle.specials(Terms, "strata")
    if (!is.null(strata)) strata.name <- ts$vars else strata.name <- NULL
###    if (!is.null(strata)) strata.name <- names(des$xlevels) 
###    else strata.name <- NULL
    pos.strata <- NULL
    des.weights <- des$weights
    des.offset  <- des$offset
    id      <- des$cluster
    ## no use of 
    pos.cluster <- NULL
 ## }}}
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

 call.id <- id;
 conid <- construct_id(id,nrow(X),as.data=TRUE)
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid

  statusE <- 1*(status %in% cause)
###  if (is.null(death.code)) statusD <- 1*(!(status %in% cens.code)) else 
  statusD <- 1*(status %in% death.code)
  data$statusE__ <- statusE
  data$statusD__ <- statusD
  data$strata__  <- strata

  if (is.null(call.id))  stop("must specify 'cluster()'\n"); 
  ### setting up formulae for the two phreg (cause of interest and death)
###  if (is.null(call.id)) { 
###     stop("must give id\n"); 
###     formid <- update.formula(formula,~.+cluster(id)) 
###     data$id <- id
###     tt <- terms(formid)
###     tt <- delete.response(tt)
###     formid <- formula(tt)
###  }  else 
  formid <- formula
  tt <- terms(formid)
  tt <- delete.response(tt)
  formid <- formula(tt)

  data$exit__ <- exit
  if (ncol(Y)==3) {
     data$entry__ <- entry
     formE <- as.formula(paste("Surv(entry__,exit__,statusE__)~+1"))
     formD <- as.formula(paste("Surv(entry__,exit__,statusD__)~+1"))
  } else {
     formE <- as.formula(paste("Surv(exit__,statusE__)~+1"))
     formD <- as.formula(paste("Surv(exit__,statusD__)~+1"))
  }
 formE <- update.formula(formE,formid)
 formD <- update.formula(formD,formid)
 ## drop marks from terminal event formula
 formD <- drop.specials(formD,"marks")

  if (sum(statusE)==0) warning("No events of type 1\n"); 
  coxE <- phreg(formE,data=data,...)
  coxS <- phreg(formD,data=data,...)

  ### cif 
  if (sum(statusE)>0) meano <- recurrentMarginalPhreg(coxE,coxS) else meano <- coxE
  ## to work with predict function
  meano$no.opt <- TRUE

  if (test) {
	  logrank <- logrankRecurrentBase(coxE,coxS)
  } else logrank <- NULL

  attr(meano,"logrank") <- logrank
  attr(meano,"cause") <- cause
  attr(meano,"death.code") <- death.code
  return(meano)
}# }}}

##' @export
recurrentMarginal <- function(formula,data,...) recurrent_marginal(formula,data,...)

##' @export
recurrentMarginalPhreg <- function(recurrent,death,fixbeta=NULL,km=TRUE)
{# {{{
  xr <- recurrent
  dr <- death 

  ### sets fixbeta based on whether xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if ((xr$no.opt) | is.null(xr$coef)) fixbeta<- 1 else fixbeta <- 0

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## survival at t- to also work in competing risks situation
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
    St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazR <-  cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
  mu <- cumhazDR[,2]
# }}}

  ### robust standard errors 
  ### 1. sum_k ( int_0^t S(s)/S_0^r(s) dM_k.^r(s) )^2
 resIM1 <-  squareintHdM(xr,ft=St,fixbeta=fixbeta)
 ### 2. mu(t)^2 * sum_k ( int_0^t 1/S_0^d(s) dM_k.^d(s) )^2
 resIM2 <-  squareintHdM(dr,ft=NULL,fixbeta=fixbeta)
  ### 3. sum_k( int_0^t mu(s) /S_0^d(s) dM_k.^d(s))^2
 resIM3 <-  squareintHdM(dr,ft=mu,fixbeta=fixbeta)

 varA <-  resIM1$varInt+mu^2*resIM2$varInt+resIM3$varInt 

## covariances between different terms  13 23  12 12
## to allow different strata for xr and dr, but still nested strata
 if ((xr$nstrata>1 & dr$nstrata==1)) {
    cM1M3 <- covIntH1dM1IntH2dM2(resIM1,resIM3,fixbeta=fixbeta,mu=NULL)
    cM1M2 <- covIntH1dM1IntH2dM2(resIM1,resIM2,fixbeta=fixbeta,mu=mu)
 } else  {
    cM1M3 <- covIntH1dM1IntH2dM2(resIM3,resIM1,fixbeta=fixbeta,mu=NULL)
    cM1M2 <- covIntH1dM1IntH2dM2(resIM2,resIM1,fixbeta=fixbeta,mu=mu)
 }
 cM2M3 <- covIntH1dM1IntH2dM2(resIM2,resIM3,fixbeta=fixbeta,mu=mu)

 varA <- varA+2*cM1M3$cov12A-2*cM1M2$cov12A-2*cM2M3$cov12A 
### varA <- varA-2*cM1M3$cov12A+2*cM1M2$cov12A+2*cM2M3$cov12A 

 cov12aa <- cov13aa <- cov23aa <- 0

 if (fixbeta==0) {
    varA <-varA + cM2M3$covbeta - cM1M3$covbeta + cM1M2$covbeta 
 }

 varrs <- data.frame(mu=mu,cumhaz=mu,se.mu=varA^.5,time=xr$time,
		     se.cumhaz=varA^.5,strata=xr$strata,St=St)
 varrs <- varrs[c(xr$cox.prep$jumps)+1,]

 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,se.mu=varrs$se.mu,times=varrs$time,
     St=varrs$St,
     cumhaz=cbind(varrs$time,varrs$mu),se.cumhaz=cbind(varrs$time,varrs$se.mu),
     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
     strata.name=xr$strata.name,strata.level=recurrent$strata.level,call=xr$call,
     formula=xr$formula,no.opt=TRUE,exit=varrs$time,p=0)
 attr(out,"recurrent") <-  recurrent
 attr(out,"death") <- death
 class(out) <- rep("recurrent",2)
 return(out)
}# }}}

##' @export
plot.recurrent <- function(x,ylab=NULL,test=TRUE,...) {# {{{
 if (inherits(x,"recurrent") & is.null(ylab)) ylab <- "Mean events"
 baseplot(x,ylab=ylab,...)
 logrank <- attr(x,"logrank")
 if (!is.null(logrank) & test)
 legend("bottomright",paste("Logrank test, p=",round(logrank$compare$p.value,3)))
}# }}}

##' @export
summary.recurrent <- function(object,...) {# {{{
out <- summaryRecurrentobject(object,...)
return(out)
}# }}}

summaryRecurrentobject <- function(object,times=NULL,strata=NULL,estimates=FALSE,
			      name="mean",cumhaz="cumhaz",se.cumhaz="se.cumhaz",robust=FALSE,
	      conf.type=c("log","log-log","plain"),...) {# {{{
base <- basecumhaz(object,joint=1,robust=robust,cumhaz=cumhaz,se.cumhaz=se.cumhaz)
nstrata <- object$nstrata
stratobs <- attr(base,"stratobs")

baseci <- rep(list(NULL),nstrata)
if (length(stratobs)>0) 
 for (i in stratobs) {
   cumhaz <- base[[i+1]]$cumhaz
   if (nrow(cumhaz)>=1) {
     mu <- base[[i+1]]$cumhaz[,2]
     se.mu <- base[[i+1]]$cumhaz[,3]
	   if (length(mu)>=1) {
	   conf <- conftype(mu,se.mu,conf.type=conf.type[1],...)
	   out <- data.frame(times=cumhaz[,1],mu=mu,se.mu=se.mu,lower=conf$lower,upper=conf$upper,strata=i)
	   names(out) <- c("times",name,"se","CI-2.5%","CI-97.5%","strata")
	   baseci[[i+1]] <- out
	   } 
   } 
   }

i <- stratobs
pbaseci <- NULL
if (!is.null(times)) {
 if (length(stratobs)>0) 
 for (i in stratobs) 
if (!is.null(baseci[[i+1]])) pbaseci[[i+1]] <- predictCumhaz(rbind(0,baseci[[i+1]]),times) 
}
out <- list(baseci=baseci,pbaseci=pbaseci,times=times)

class(out) <- "summary.recurrent"
return(out)
}# }}}

##' @export
print.summary.recurrent  <- function(x,...) {# {{{
if (is.null(x$times)) print(x$baseci) else print(x$pbaseci)
} # }}}

iidRecurrent <- function(recurrent,death,wt=NULL,km=TRUE,start=0,time=NULL) { ## {{{ 
xr <- recurrent
dr <- death 
if (is.null(time)) time <- max(xr$jumptimes)
### marginal expected events  int_0^t w(s) S(s) \lambda_r(s) ds 
### weight is considered known 
## {{{
###
x <- dr
xx <- x$cox.prep
xstrata <- xx$strata
nstrata <- xx$nstrata

S0i2 <- S0i <- rep(0,length(xstrata))
S0i[xx$jumps+1] <-  1/x$S0
S0i2[xx$jumps+1] <- 1/x$S0^2
## survival at t- to also work in competing risks situation
if (!km) { 
cumhazD <- c(cumsumstratasum(S0i,xstrata,nstrata)$lagsum)
St      <- exp(-cumhazD)
} else St <- c(exp(cumsumstratasum(log(1-S0i),xstrata,nstrata)$lagsum))
x <- xr
xx <- x$cox.prep
S0i2 <- S0i <- rep(0,length(xstrata))
S0i[xx$jumps+1] <-  1/x$S0
S0i2[xx$jumps+1] <- 1/x$S0^2

btime <- (xx$time<=time  & xx$time>start)*1
wt.call <- wt
if (is.null(wt)) wt <- rep(1,length(xx$time))

cumhazR <-  cbind(xx$time,cumsumstrata(S0i*btime,xstrata,nstrata))
meanN <- cumhazDR <- cbind(xx$time,cumsumstrata(wt*St*S0i*btime,xstrata,nstrata))
mu <- cumhazDR[,2]
meanN <- meanN[xx$jumps+1,]

cumhazSR <- cumsumstrata(wt*St*S0i2*btime, xstrata, nstrata)
rr <- c(xx$sign * exp( xx$offset))
MSGt <- wt*St*S0i*btime - cumhazSR*rr*c(xx$weights)
id <- xx$id
MSGtiid <- apply(MSGt,2,sumstrata,id,max(id)+1)

x <- dr
xx <- x$cox.prep
rr <- c(xx$sign * exp( xx$offset))
S0i2 <- S0i <- rep(0,length(xx$strata))
S0i[xx$jumps+1] <-  1/x$S0
S0i2[xx$jumps+1] <- 1/x$S0^2

cumhazD <-  cumsumstrata(S0i2*btime,xstrata,nstrata)
cumhazmuD <-  cumsumstrata(mu*S0i2*btime,xstrata,nstrata)
tails <- tailstrata(xstrata,nstrata)
mut <- mu[tails][xstrata+1]

mutMDt <- mut*S0i*btime-mut*cumhazD*rr*c(xx$weights)
mutMDtiid <- apply(mutMDt,2,sumstrata,id,max(id)+1)

MmuDt <- mu*S0i*btime-cumhazmuD*rr*c(xx$weights)
MmuDtiid <- apply(MmuDt,2,sumstrata,id,max(id)+1)

MGAiid <- MSGtiid-mutMDtiid+MmuDtiid
mid <- max(id)+1

strataid <- xstrata[headstrata(xx$id,mid)]
id <- 1:nrow(MGAiid)

MGAiids <- matrix(0,nrow(MGAiid),nstrata)
cumhaz.time <- c()
sus <- sort(unique(strataid))
for (i in sus)  { 
 wi <- which(strataid==i)
 MGAiids[wi,i+1] <- MGAiid[wi,]
}
MGAiid <- MGAiids
if (is.matrix(MGAiid)) {
colnames(MGAiid) <- paste("strata",sus,sep="")
if (!is.null(x$call.id)) MGAiid <- nameme(MGAiid,x$name.id) 
}

var.iid <- crossprod(MGAiid)
se <- diag(var.iid)^.5
## }}} 

out <- list(iid=MGAiid,wt=wt.call,time=time,start=start,meanN=meanN,var=var.iid,se=diag(var.iid)^.5,mut=mu[tails])
return(out)
} ## }}}

##' @export
estimate.recurrent <- function(x,time=NULL,...) { ## {{{ 
   if (is.null(time)) stop("must give time for iid decomposition \n")
   ic <- IC(x,time=time)
   out <- estimate(coef=attr(ic,"coefs"),IC=ic,...)
return(out)
} ## }}}

##' @export
iid.recurrent <- function(x,...) { ## {{{ 
iid <- iidRecurrent(attr(x,"recurrent"),attr(x,"death"),...)
out <- iid$iid
attr(out,"coefs") <- iid$mut 
attr(out,"time") <-  iid$time
return(out)
} ## }}}

##' @export
IC.recurrent <- function(x,...) { ## {{{ 
iid <- iid(x,...)
ic <- iid*nrow(iid)
return(ic)
} ## }}}

##' Logrank-type test for comparing recurrent event marginal means between groups
##'
##' Tests whether the marginal mean number of recurrent events differs across
##' groups (strata), extending the classical logrank test to the setting of
##' recurrent events with a competing terminal event. The test statistic is
##' \deqn{z = \int_0^\tau w(s)\bigl[d\hat\mu_1(s) - d\hat\mu_2(s)\bigr],}
##' where \eqn{w(s)} is a weight function and \eqn{\hat\mu_j(s)} is the estimated
##' marginal mean for group \eqn{j}. Variance is estimated robustly via the
##' influence functions of Ghosh and Lin (2000).
##'
##' Three weight schemes are available:
##' \describe{
##'   \item{\code{"I"}}{(Default) \eqn{w(t) = R_1(t) R_2(t) / (R_1(t) + R_2(t))},
##'     where \eqn{R_j(t) = Y_j(t) / \hat S_j(t-)}. Analogous to the standard
##'     logrank weight.}
##'   \item{\code{"II"}}{\eqn{w(t) = Y_j(t)}, the raw risk-set size. Equivalent
##'     to using observed counts without survival adjustment.}
##'   \item{\code{"III"}}{A modified weight incorporating the cumulative incidence,
##'     analogous to Gray's test for competing risks.}
##' }
##'
##' @param recurrent Either a \code{"recurrent"} object returned by
##'   \code{\link{recurrent_marginal}}, or a \code{"phreg"} object for the
##'   recurrent event model (in which case \code{death} must also be supplied).
##' @param death A \code{"phreg"} object for the terminal event model. Required
##'   when \code{recurrent} is a \code{"phreg"} object; ignored otherwise.
##' @param weight Character string specifying the weight scheme: \code{"I"},
##'   \code{"II"}, or \code{"III"}. Default is \code{"I"}.
##' @param km Logical. If \code{TRUE} (default), the Kaplan-Meier estimator is
##'   used for the survival probability \eqn{S(t)}; otherwise the Nelson-Aalen
##'   estimator is used.
##' @param start Left truncation time for the integration. Default is \code{0}.
##' @param stop Right truncation time for the integration. Defaults to the last
##'   observed jump time.
##' @param at.risk Minimum combined risk-set size below which the weight is set
##'   to zero. Default is \code{5}.
##' @param cluster.id Optional vector of cluster identifiers for aggregating
##'   influence functions across clusters before forming the test statistic.
##' @param ... Currently unused.
##'
##' @return An object of class \code{"estimate"} (from the \pkg{lava} package)
##'   with the following components:
##'   \item{coef}{The weighted difference in marginal means between groups.}
##'   \item{se}{Robust standard error of the test statistic.}
##'   \item{lower, upper}{95\% confidence interval bounds.}
##'   \item{p.value}{Two-sided p-value for the null hypothesis of no difference.}
##'   The object also carries an \code{iid} attribute containing the
##'   subject-level influence function decomposition of the test statistic,
##'   which can be used for further inference or combination with other estimators.
##'
##' @author Thomas Scheike
##'
##' @seealso \code{\link{recurrent_marginal}}, \code{\link{logrankRecurrentBase}}
##'
##' @references
##' Ghosh, D. and Lin, D. Y. (2000). Nonparametric analysis of recurrent events
##' and death. \emph{Biometrics}, 56, 554--562.
##'
##' @examples
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##'
##' ## Test using two separate phreg models
##' xr <- phreg(Surv(entry, time, status == 1) ~ strata(treatment) + cluster(id), data = hf)
##' dr <- phreg(Surv(entry, time, status == 2) ~ strata(treatment) + cluster(id), data = hf)
##' out <- test_logrankRecurrent(xr, dr, stop = 5)
##' summary(out)
##'
##' ## Equivalently, using a recurrent_marginal object directly
##' outN <- recurrent_marginal(Event(entry, time, status) ~ strata(treatment) + cluster(id),
##'                            data = hf, cause = 1, death.code = 2)
##' test_logrankRecurrent(outN)
##'
##' @aliases logrankRecurrentBase
##' @export
test_logrankRecurrent <- function(recurrent,death,
              weight=c("I","II"),km=TRUE,start=0,stop=NULL,at.risk=5,cluster.id=NULL,...) { ## {{{ 
  if (inherits(recurrent,"phreg")) { # Fall-back to recurrentMarginalPhreg
    if (inherits(death, "phreg")) {
        return(logrankRecurrentBase(recurrent,death,
		weight=weight[1],km=km,start=start,stop=stop,at.risk=at.risk,
		cluster.id=cluster.id))
     } 
     return(NULL)
  } else if (inherits(recurrent,"recurrent")) { # Fall-back to recurrentMarginalPhreg
        return(logrankRecurrentBase(attr(recurrent,"recurrent"),attr(recurrent,"death"),
		weight=weight[1],km=km,start=start,stop=stop,at.risk=at.risk,cluster.id=cluster.id))
 } else stop("input either output from recurrentMarginal or two phreg's\n")

} ## }}}

##' @export
logrankRecurrentBase <- function(recurrent,death,weight=c("I","II","III"),km=TRUE,
       start=0,stop=NULL,at.risk=0,cluster.id=NULL) { ## {{{ 
  xr <- recurrent
  dr <- death 
  if (is.null(stop)) stop <- max(xr$jumptimes)

  nlev <- xr$nstrata
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  n <- dr$n

  Yrr <- revcumsumstrata(xx$sign,xx$strata,xx$nstrata)
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/Yrr[xx$jumps+1]

  Sta <- cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)
  St <- c(exp(Sta$lagsum))
  kmss <- exp(Sta$alllagsum)
  Yss <- .Call("_mets_riskstrataR",xx$sign,xx$strata,xx$nstrata)$risk

  if (weight[1]=="I") {
	  Rss <-  Yss/kmss
	  Rss[Yss==0] <- 0 
  }
  if (weight[1]=="III") {
     xx <- xr$cox.prep
     S0ri <- rep(0,length(xx$strata))
     S0ri[xx$jumps+1] <-  1/Yrr[xx$jumps+1]
     cifsm <- cumsumstratasum(St*S0ri,xx$strata,xx$nstrata)
     F1c <-  1-cifsm$alllagsum 
     Rss <-  F1c*Yss/kmss
     Rss[Yss==0] <- 0 
     F1ci <- rep(0,length(xx$strata))
     F1ci[xx$jumps+1] <-  1/(1-cifsm$lagsum[xx$jumps+1])
  }
  if (weight[1]=="II") {
     Rss <- Yss 
  }
  Rp <- apply(Rss,1,sum,na.rm=TRUE)
  Rssstrat <- mdi(Rss,1:length(xx$strata),xx$strata+1)

  contr <- contr.iid <- c()
  i <- 1
  for (i in 1:(nlev-1)) {
     wt <- rep(0,length(xx$strata))
     strati  <- which(xx$strata==(i-1))
     wt <- Rss[,i]*Rssstrat/Rp
     wt[strati] <- Rss[strati,i]*(Rp[strati]-Rss[strati,i])/Rp[strati]
     wt[is.na(wt)] <- 0
      if (at.risk>0) {
        minrisk <- pmin(Yss[,i],apply(Yss[,-i,drop=FALSE],1,sum))
        wt[minrisk<=at.risk] <- 0
      }
     ## adjust weight to get Grays test with [1-F_1(s-)]^-1
     if (weight[1]=="III") wt <- wt*F1ci

      lrt <- iidRecurrent(xr,dr,wt=wt)
      diff(lrt$mut)
      if (!is.null(cluster.id)) {
	      nid <- length(unique(cluster.id))
	      conid <- construct_id(cluster.id,nid)
	      lrt$iid <- apply(lrt$iid,2,sumstrata,conid$id,conid$nid)
      }
      diffi <- lrt$mut[i]-sum(lrt$mut[-i])
      contr <- c(contr,diffi)
      contr.iid <- cbind(contr.iid,lrt$iid[,i]-apply(lrt$iid[,-i,drop=FALSE],1,sum))
   }
   logrank <- estimate(coef=contr,IC=contr.iid*nrow(contr.iid),null=0)
   p <- length(contr)

   return(logrank)
} ## }}} 

##' @export
recurrent_marginalAIPCW <- function(formula,data=NULL,cause=1,cens.code=0,death.code=2,
	  cens.model=~1,km=TRUE,times=NULL,augment.model=~Nt,...)
{# {{{
    cl <- match.call()# {{{
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster","offset")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y,"Event")) stop("Expected a 'Event'-object")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- NULL ## rep(0,nrow(Y))
        status <- Y[,2]
    } else {
        entry <- Y[,1]
        exit <- Y[,2]
        status <- Y[,3]
    }
    id <- strata <- NULL
    if (!is.null(attributes(Terms)$specials$cluster)) {
        ts <- survival::untangle.specials(Terms, "cluster")
        pos.cluster <- ts$terms
        Terms  <- Terms[-ts$terms]
        id <- m[[ts$vars]]
    } else pos.cluster <- NULL
    if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
        ts <- survival::untangle.specials(Terms, "strata")
        pos.strata <- ts$terms
        Terms  <- Terms[-ts$terms]
        strata <- m[[ts$vars]]
        strata.name <- ts$vars
    }  else { strata.name <- NULL; pos.strata <- NULL}
    if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
        ts <- survival::untangle.specials(Terms, "offset")
        Terms  <- Terms[-ts$terms]
        offset <- m[[ts$vars]]
    }
    X <- model.matrix(Terms, m)
    if (!is.null(intpos  <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    ## }}}

   if (!is.null(id)) {
        ids <- unique(id)
        nid <- length(ids)
        if (is.numeric(id))
            id <-  fast.approx(ids,id)-1
        else  {
            id <- as.integer(factor(id,labels=seq(nid)))-1
        }
    } else { id <- as.integer(seq_along(entry))-1;  nid <- nrow(X); }
    ## orginal id coding into integers 1:...
    id.orig <- id+1;

   ### setting up with artificial names
 data$status__ <-  status 
 data$id__ <-  id
 ## lave Countcause
 data <- count_history(data,status="status__",id="id__",types=cause,multitype=TRUE)
 data$Count1__ <- data[,paste("Count",cause[1],sep="")]
 data$death__ <- (status %in% death.code)*1
 data$entry__ <- entry 
 data$exit__ <- exit 
 data$statusC__ <- (status %in% cens.code)*1
 data$status__cause <- (status %in% cause)*1

  xr <- phreg(Surv(entry__,exit__,status__cause)~Count1__+death__+cluster(id__),data=data,no.opt=TRUE,no.var=1)

  formC <- update.formula(cens.model,Surv(entry__,exit__,statusC__)~ . +cluster(id__))
  cr <- phreg(formC,data=data)
  whereC <- which(status %in% cens.code)

  if (length(whereC)>0) {# {{{
  ### censoring weights
  strat <- cr$strata[cr$jumps]
  x <- cr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## survival at t- to also work in competing risks situation
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
    St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  } else  St <- rep(1,length(xx$strata))
  Gc <- St
  ## }}}

  
formD <- as.formula(Surv(entry__,exit__,death__)~cluster(id__))
form1L <- as.formula(Surv(entry__,exit__,status__cause)~Count1__+death__+statusC__+cluster(id__))
form1 <- as.formula(Surv(entry__,exit__,status__cause)~cluster(id__))

 xr <- phreg(form1L,data=data,no.opt=TRUE,no.var=1)
 dr <- phreg(formD,data=data,no.opt=TRUE,no.var=1)

 ### augmenting partioned estimator computing \hat H_i(s,t) for fixed t
 data$Gctrr <- exp(-cpred(rbind(0,cr$cumhaz),exit)[,2])

 ### cook-lawless-ghosh-lin
 xr0 <- phreg(form1,data=data,no.opt=TRUE)
 clgl  <- recurrentMarginalPhreg(xr0,dr)
### bplot(clgl,se=1); print(cpred(clgl$cumhaz,times)); print(cpred(clgl$se.cumhaz,times)); 

  ####  First \mu_ipcw(t) \sum_i I(T_i /\ t \leq C_i)/G_c(T_i /\ t ) N_(T_i /\ t) {{{
  x <- xr
  xx <- xr$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## stay with N(D_i) when t is large so no -1 when death
  ## xx$X[,1] er Count1 dvs N(t-)
  Nt <- revcumsumstrata(xx$X[,1]*xx$sign,xx$strata,xx$nstrata)
  Nt <- Nt/Gc
  ## counting N(D) forward in time skal ikke checke ud når man dør N_(D_i) er i spil efter D_i
  NtD <- cumsumstrata(xx$X[,1]*(xx$X[,2]==1)*(xx$sign==1)/Gc,xx$strata,xx$nstrata)
  jump1 <- xx$jumps+1
  timeJ <- xx$time[jump1]
  avNtD <- (NtD+Nt)[jump1]/nid
  strataN1J <- xx$strata[jump1]

  varIPCW1 <- NULL
  seIPCW1 <- NULL

  ### IPCW estimator 
  cumhaz <- cbind(timeJ,avNtD)
# }}}

  ### Partitioned estimator , same as Lin, Lawless & Cook estimator {{{
  cumhazP <- c(cumsumstrata(1/Gc[jump1],strataN1J,xx$nstrata)/nid)
  cumhazP <- cbind(timeJ,cumhazP)

  ### variance of partitioned estimator 

  ### calculate E(H,s,t) = E(H,t) - E(H,s) 
  ### E(H,t) = 1/(S(t)*n) \sum_ \int Y_i(s)/G_c(s) dN_{1i} = 1/(S(t)) (\mu_p(t) - \mu_p(s)) 
  ### \hat H_i for all subjects, and look at together with Hst, eHst 
  ### for censoring martingale 
  data$Nt <- data$Count1__
  data$Nt2 <- data$Nt^2
  data$expNt <- exp(-data$Nt)
  data$NtexpNt <- data$Nt*exp(-data$Nt)

  Gcdata <- exp(-cpred(rbind(0,cr$cumhaz),exit)[,2])

  form <- as.formula(paste("Surv(entry__,exit__,statusC__)~+1"))
  desform <- update.formula(augment.model,~Hst + . + cluster(id__))
  form[[3]] <- desform[[2]]
  nterms <- length(all.vars(form[[3]]))-1

  if (!is.null(times)) {
       semuPA <-  muPA <- semuPA.times <-  muPA.times <- rep(0,length(times))
       ww <- fast.approx(timeJ,times,type="left")
       muP.times <- cumhazP[ww,2]
       semuP.times <- clgl$se.mu[ww]

  for (i in seq_along(times)) {
     timel <- times[i]
     data$Hst <- revcumsumstrata((exit<timel)*(status %in% cause)/Gcdata,id,nid)
     cr2 <- phreg(form,data=data,no.opt=TRUE,no.var=1)
     nterms <- cr2$p-1

     dhessian <- cr2$hessianttime
     dhessian <-  .Call("XXMatFULL",dhessian,cr2$p,PACKAGE="mets")$XXf
     ###  matrix(apply(dhessian,2,sum),3,3)
     timeb <- which(cr$cumhaz[,1]<timel)
     ### take relevant \sum H_i(s,t) (e_i - \bar e)
     covts <- dhessian[timeb,1+1:nterms,drop=FALSE]
     ### construct relevant \sum (e_i - \bar e)^2
     Pt <- dhessian[timeb,-c((1:(nterms+1)),(1:(nterms))*(nterms+1)+1),drop=FALSE]
     ###  matrix(apply(dhessian[,c(5,6,8,9)],2,sum),2,2)
     gammahat <- .Call("CubeVec",Pt,covts,1,PACKAGE="mets")$XXbeta
     S0 <- cr$S0[timeb]
     gammahat[is.na(gammahat)] <- 0
     gammahat[gammahat==Inf] <- 0
     Gctb <- Gc[cr$cox.prep$jumps+1][timeb]
     augment.times <- sum(apply(gammahat*cr2$U[timeb,1+1:nterms,drop=FALSE],1,sum))/nid
     mterms <- length(terms)
     mterms <- nterms
     ###
     varZ <- matrix(apply(Pt/Gctb^2,2,sum),mterms,mterms)
     gamma <- .Call("CubeVec",matrix(c(varZ),nrow=1),matrix(apply(covts/Gctb,2,sum),nrow=1),1,PACKAGE="mets")$XXbeta
     gamma <- c(gamma)
     gamma[is.na(gamma)] <- 0
     gamma[gamma==Inf] <- 0
     augment <- sum(apply(gamma*t(cr2$U[timeb,1+1:nterms,drop=FALSE])/Gctb,2,sum))/nid
     ###
     muPA[i] <- muP.times[i]+augment
     semuPA[i] <- (semuP.times[i]^2 +(gamma %*% varZ %*% gamma)/nid^2)^.5
     muPA.times[i] <- muP.times[i]+augment.times
     semuPA.times[i] <- (semuP.times[i]^2+sum(gammahat*.Call("CubeVec",Pt,gammahat,0,PACKAGE="mets")$XXbeta)/(nid^2))^.5
  }
  } else { Gctb <- NULL;gammahat <- NULL; 
  muPA.times <-semuPA.times <- muPA  <- semuPA <- NULL }
# }}}

  return(list(censoring.weights=Gctb,muP.all=cumhazP,Gcjump=Gc[jump1],gamma=gamma,gamma.time=gammahat,times=times,
  muP=muP.times,semuP=semuP.times, muPAt=muPA.times,semuPAt=semuPA.times, muPA=muPA,semuPA=semuPA))
}# }}}

##' @export
summaryTimeobject <-function(mutimes,mu,se.mu=NULL,times=NULL,type="log",level=0.95,...) {# {{{
 if (is.null(times)) times <- mutimes

 where <- fast.approx(c(0,mutimes),times,type="left")

 crit <- qnorm(1-(1-level)/2)
 ##  see if object is vector or matrix
 if (is.matrix(mu)) mu <- rbind(0,mu)[where,] else mu <- c(0,mu)[where]
 if (!is.null(se.mu)) {
     if (is.matrix(se.mu)) se.mu <- rbind(0,se.mu)[where,] else se.mu <- c(0,se.mu)[where]
 se.logmu=se.mu/mu
 if (type=="log") {
 lower <- exp(log(mu) - crit*se.logmu)
 upper <- exp(log(mu) + crit*se.logmu)
 } else {
 lower <- mu - crit*se.mu
 upper <- mu + crit*se.mu
 }
 } else {se.mu <- se.logmu <- lower <- upper <- NA}

 out <- data.frame(times=times,mu=mu,se.mu=se.mu,lower=lower,upper=upper)
 names(out) <- c("times","mean","se-mean","CI-2.5%","CI-97.5%")
 return(out)
}# }}}

recmarg <- function(recurrent,death,Xr=NULL,Xd=NULL,km=TRUE,...)
{# {{{
  xr <- recurrent
  dr <- death 

  if (!is.null(Xr)) rr <- exp(sum(xr$coef * Xr)) else rr <- 1
  if (!is.null(Xd)) rrd <- exp(sum(dr$coef * Xd)) else rrd <- 1

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD*rrd)
  } else St <- exp(rrd*c(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  ###
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
  mu <- rr*cumhazDR[,2]
# }}}

 varrs <- data.frame(mu=mu,time=xr$time,strata=xr$strata,St=St)
 varrs  <-  varrs[c(xr$cox.prep$jumps)+1,]

 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,time=varrs$time,
	     St=varrs$St,cumhaz=cbind(varrs$time,varrs$mu),
             strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
	     strata.name=xr$strata.name,strata.level=recurrent$strata.level)
 return(out)
}# }}}

squareintHdM <- function(phreg,ft=NULL,fixbeta=NULL,beta.iid=NULL,...)
{# {{{
###  sum_k ( int_0^t f(s)/S_0^r(s) dM_k.^r(s) )^2
###  strata "r" from object and "k" id from cluster 
  if (!inherits(phreg,"phreg")) stop("Must be phreg object\n"); 

  ### sets fixbeta based on whether xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if ((phreg$no.opt) | is.null(phreg$coef)) fixbeta<- 1 else fixbeta <- 0

  x <- phreg
  xx <- x$cox.prep
  ww <- xx$caseweights*xx$weights
 
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/(x$S0^2*ww[xx$jumps+1])
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
  U[xx$jumps+1,] <- x$U
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  if (is.null(ft))  ft <- rep(1,length(xx$time))
  cumS0i2 <- c(cumsumstrata(ft*S0i2,xx$strata,xx$nstrata))
  if (fixbeta==0) {
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  Ht <- apply(ft*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  } else Ht <- NULL
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
  else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  ### also weights 
  w <- c(xx$weights)
  xxx <- (ft*S0i-rr*cumS0i2*w)
  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
  ss <-  revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare*cumS0i2^2
  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*cumS0i2
  varA1 <- c(ssf+ss-2*covv)

  vbeta <- betaiidR <- NULL
  if (fixbeta==0) {# {{{
    if (!is.null(beta.iid))  betaiidR <- beta.iid else {
	     invhess <- -pinv(x$hessian)
	     MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-EdLam0)*rr*c(xx$weights)
	     UU <- apply(MGt,2,sumstrata,id,mid)
	     betaiidR <- UU %*% invhess
     }
     vbeta <- crossprod(betaiidR)
     varbetat <-   rowSums((Ht %*% vbeta)*Ht)
     ### writing each beta for all individuals 
     betakt <- betaiidR[id+1,,drop=FALSE]
     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- covk2*cumS0i2
     covv <- covk1-covk2
     varA1 <- varA1+varbetat-2*apply(covv*Ht,1,sum)
  }# }}}

  return(list(xx=xx,Ht=Ht,varInt=varA1,xxx=xxx,rr=rr,
	      cumhaz=cumhaz,cumS0i2=cumS0i2,mid=mid,id=id,
	      betaiid=betaiidR,vbeta=vbeta,covv=covv))
} # }}}

covIntH1dM1IntH2dM2 <- function(square1,square2,fixbeta=1,mu=NULL)
{# {{{

 ### strata and id same for two objects 
 xx <- square1$xx; xx2 <- square2$xx
 xxxR <- square1$xxx;     xxxD1 <- square2$xxx
 rrR  <- square1$rr;       rrD1 <- square2$rr
 id   <- id1 <- square1$id; id2 <- square2$id
 mid  <- square1$mid; w <- c(xx$weights)

 if (is.null(mu)) mu <- rep(1,length(xx$strata))

 cov12 <-  c(cumsumidstratasumCov(xxxR,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)
 cov122 <- c(revcumsumidstratasumCov(w*rrR,w*rrD1,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(square1$cumS0i2*square2$cumS0i2)
 cov123 <- covfridstrataCov(xxxR,w*rrR,xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(square2$cumS0i2)
 cov124 <- covfridstrataCov(xxxD1,w*rrD1,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(square1$cumS0i2)
 cov12A <- c(cov12+cov122+cov123+cov124)
 
 test <- 0
 if (test==1) {# {{{
	 print("________cov cov ___________________________"); 
	 print(summary(c(xxxR))); print(summary(c(xxxD1))); 
	 print(summary(c(rrD1))); print(summary(c(rrR)))
	 print(summary(c(square1$cumS0i2))); print(summary(c(square2$cumS0i2)));
	 print("-----------"); 
	 print(summary(cov12)); 
	 print(summary(cov122)); 
	 print(summary(c(cov123))); 
	 print(summary(c(cov124))); 
	 print("______________________________________"); 
	 jumps <- c(square1$xx$jumps,square2$xx$jumps)+1
	 print(summary(jumps))
         print(summary(cov12[jumps])); 
	 print(summary(cov122[jumps])); 
	 print(summary(c(cov123[jumps]))); 
	 print(summary(c(cov124[jumps]))); 
 }# }}}

 cov12aa <- 0
 if (fixbeta==0) {
 ### covariances between different terms and  beta's 
 # {{{
	 betaiidR <- square1$betaiid; betaiidD <- square2$betaiid
	 HtR <- square1$Ht; HtD <- square2$Ht
	 covbetaRD <- t(betaiidR) %*% betaiidD
	 covbeta <-   -1*rowSums((HtR %*% covbetaRD)*HtD)

	 ### cov12 wrt betaD and betaR
	 betakt <- betaiidD[id1+1,,drop=FALSE]
	 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id1,mid,xx$strata,xx$nstrata,type="sum")
	 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id1,mid,xx$strata,xx$nstrata,type="lagsum")
	 covk2 <- c(covk2)*c(square1$cumS0i2)
	 covRD12 <- apply((covk1-covk2)*HtD,1,sum)
	 betakt <- betaiidR[id2+1,,drop=FALSE]
	 covk1 <-apply(xxxD1*betakt,2,cumsumidstratasum,id2,mid,xx2$strata,xx$nstrata,type="sum")
	 covk2 <-apply(w*rrD1*betakt,2,revcumsumidstratasum,id2,mid,xx2$strata,xx$nstrata,type="lagsum")
	 covk2 <- c(covk2)*c(square2$cumS0i2)
	 covRD21 <- apply((covk1-covk2)*HtR,1,sum)
	 cov12aa <- 2*(covbeta + covRD12+covRD21)
	 test <- 0
	 if (test==1) {
		 print("--------------")
		 print(summary(2*mu*covbeta))
		 print(summary(2*mu*covRD12))
		 print(summary(2*mu*covRD21))
		 print("--------------")
	 }
 } # }}}

 cov12 <- (cov12A-cov12aa)*mu

 return(list(cov=cov12,cov12A=cov12A*mu,covbeta=cov12aa*mu))
} # }}}

##' Break ties in event times for recurrent event data
##'
##' Resolves tied event times in a counting-process dataset by adding a small
##' random perturbation to duplicated exit times of event rows. This is a
##' preprocessing step required by \code{\link{recurrent_marginal}} and related
##' functions, which assume unique jump times within each stratum.
##'
##' A tie is defined as an event exit time (rows where \code{status} is a cause
##' code) that coincides with another exit time in the dataset. When
##' \code{exit.unique = TRUE} (default), a tie is flagged whenever an event time
##' also appears among any other exit times (censored or event). When
##' \code{exit.unique = FALSE}, only exact ties between two event rows are
##' resolved.
##'
##' Tied event times are perturbed by adding \eqn{U \cdot \delta} where
##' \eqn{U \sim \text{Uniform}(0, 1)} and \eqn{\delta} is \code{ddt} (defaulting
##' to half the smallest observed positive gap between consecutive exit times).
##' When subject IDs are provided via \code{id}, the corresponding interval start
##' time of the immediately following row for the same subject is updated to
##' maintain a valid counting-process structure, and a logical \code{tiebreaker}
##' column is added to flag affected rows.
##'
##' @param data A data frame in counting-process format, sorted by subject and
##'   time.
##' @param stop Name of the column containing interval exit (stop) times.
##'   Default is \code{"time"}.
##' @param start Name of the column containing interval entry (start) times,
##'   used to update the following row when \code{id} is supplied. Default is
##'   \code{"entry"}.
##' @param status Name of the column containing event status codes. Default is
##'   \code{"status"}.
##' @param id Name of the column containing subject identifiers. If supplied,
##'   the start time of the next interval for the same subject is adjusted to
##'   match the perturbed stop time, preserving interval continuity. Default is
##'   \code{NULL} (no adjustment made).
##' @param ddt Maximum perturbation size. Tied event times are shifted by a
##'   uniform draw on \eqn{[0, \text{ddt}]}. If \code{NULL} (default), \code{ddt}
##'   is set to half the smallest positive gap between any two consecutive exit
##'   times in the data.
##' @param exit.unique Logical. If \code{TRUE} (default), an event time is
##'   considered tied whenever it coincides with \emph{any} exit time in the
##'   data (including censored rows). If \code{FALSE}, only ties among event rows
##'   are resolved.
##' @param cause Integer vector of status codes that identify events (non-censored
##'   rows). If \code{NULL} (default), all non-censoring status values are treated
##'   as events.
##' @param cens.code Integer code(s) for censoring. Rows with this status are
##'   never perturbed. Default is \code{0}.
##' @param seed Optional integer passed to \code{set.seed()} before drawing the
##'   perturbations, making results reproducible. Default is \code{NULL} (no seed
##'   set).
##'
##' @return The input data frame \code{data} with tied event exit times
##'   perturbed. If \code{id} is supplied, an additional logical column
##'   \code{tiebreaker} marks rows whose start time was adjusted as a consequence
##'   of a perturbation in the preceding row.
##'
##' @author Thomas Scheike
##'
##' @seealso \code{\link{recurrent_marginal}}, \code{\link{test_logrankRecurrent}}
##'
##' @examples
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##'
##' ## Check for ties in event exit times
##' ev <- hf[hf$status == 1, ]
##' any(duplicated(ev$time))
##'
##' ## Resolve ties before fitting the marginal mean model
##' hf_clean <- tie_breaker(hf, stop = "time", start = "entry",
##'                          status = "status", id = "id",
##'                          cause = 1, cens.code = 0)
##'
##' out <- recurrent_marginal(Event(entry, time, status) ~ cluster(id),
##'                            data = hf_clean, cause = 1, death.code = 2)
##' summary(out, times = 1:5)
##'
##' @export
tie_breaker <- function(data, stop = "time", start = "entry", status = "status",
                        id = NULL, cause = NULL, cens.code = 0,
                        exit.unique = TRUE, ddt = NULL, seed = NULL)
{# {{{
   id.col <- id
   if (!is.null(id)) id <- data[, id]
   ord  <- 1:nrow(data)
   stat <- data[, status]
   time <- data[, stop]

   ## determine event codes, excluding censoring
   if (is.null(cause)) cause <- unique(stat)
   type0 <- which(cause %in% cens.code)
   if (length(type0) > 0) cause <- cause[-type0]
   jumps <- stat %in% cause

   ## detect tied event times
   time1 <- data[jumps,  stop]
   time0 <- data[!jumps, stop]
   lt0   <- length(time0)
   ddp   <- duplicated(c(time0, time1))
   if (exit.unique)
     ties <- ddp[(lt0 + 1):nrow(data)]
   else
     ties <- duplicated(c(time1))

   nties <- sum(ties)

   ## initialise tiebreaker column unconditionally
   data[, "tiebreaker"] <- FALSE

   if (nties >= 1) {
     message(nties, " tied event time(s) found and perturbed.")

     ordties <- ord[jumps][ties]

     ## compute ddt from sorted unique event-time gaps, not row differences
     if (is.null(ddt)) {
       ut  <- sort(unique(time[jumps]))
       gaps <- diff(ut)
       ddt <- min(gaps[gaps > 0]) * 0.5
     }

     ## apply reproducible perturbation
     if (!is.null(seed)) set.seed(seed)
     time[ordties] <- time[ordties] + runif(nties) * ddt
     data[ordties, stop] <- time[ordties]

     ## mark perturbed rows
     data[ordties, "tiebreaker"] <- TRUE

     ## when id is available, propagate stop -> start for the next interval
     ## of the same subject to maintain a valid counting-process structure
     if (!is.null(id.col)) {
       ties_flag  <- ord %in% ordties
       lagties    <- dlag(ties_flag)
       change.start <- lagties == TRUE & id == dlag(id)
       change.start[is.na(change.start)] <- FALSE
       ocs <- ord[change.start]
       data[ocs, start] <- data[ocs - 1, stop]
     }
   }

   return(data)
 } # }}}


##' Simulate recurrent events with two event types and a terminal event
##'
##' Simulates recurrent event data with up to two distinct event types and an
##' optional terminal event (death), based on user-supplied cumulative hazard
##' functions. Dependence between processes can be introduced via shared or
##' correlated gamma-distributed frailties.
##'
##' The simulation proceeds by sequentially drawing the next event time from the
##' specified cumulative hazards, taking the minimum of the two recurrent event
##' times, and stopping each subject at death or administrative censoring.
##'
##' Dependence between processes is controlled by \code{dependence}:
##' \describe{
##'   \item{\code{0}}{Independence: all subjects have frailty fixed at 1.}
##'   \item{\code{1}}{Shared frailty: all processes share a single
##'     gamma-distributed random effect with mean 1 and variance \code{var.z}.}
##'   \item{\code{4}}{Recurrent-event frailty only: the two recurrent event
##'     processes share a gamma frailty but the terminal event is independent.}
##' }
##' For more complex correlation structures across two event types and death, use
##' \code{\link{sim_recurrentTS}}.
##'
##' @param n Number of subjects to simulate.
##' @param cumhaz Two-column matrix \code{(time, cumhaz)} giving the cumulative
##'   hazard of the first type of recurrent event.
##' @param cumhaz2 Two-column matrix \code{(time, cumhaz)} giving the cumulative
##'   hazard of the second type of recurrent event.
##' @param death.cumhaz Two-column matrix \code{(time, cumhaz)} giving the
##'   cumulative hazard of the terminal event. If \code{NULL}, no terminal event
##'   is simulated and follow-up ends at the end of \code{cumhaz}.
##' @param r1 Optional numeric vector of length \code{n} with subject-specific
##'   relative risk multipliers for the first event type.
##' @param r2 Optional numeric vector of length \code{n} with subject-specific
##'   relative risk multipliers for the second event type.
##' @param rd Optional numeric vector of length \code{n} with subject-specific
##'   relative risk multipliers for the terminal event.
##' @param rc Optional numeric vector of length \code{n} with subject-specific
##'   multipliers for the exponential censoring rate.
##' @param gap.time Logical. If \code{TRUE}, event times are drawn as gap times
##'   (time since the last event) rather than calendar times. Default is
##'   \code{FALSE}.
##' @param max.recurrent Maximum number of recurrent events allowed per subject.
##'   Default is \code{100}.
##' @param dependence Integer specifying the frailty structure. One of \code{0}
##'   (independence), \code{1} (shared gamma frailty), or \code{4} (shared frailty
##'   for recurrent events only). See Details.
##' @param var.z Variance of the gamma-distributed frailty. Default is \code{1}.
##' @param cor.mat Correlation matrix for the random effects. Used when
##'   \code{dependence = 2} (in \code{\link{sim_recurrent_list}}).
##' @param cens Rate of exponential censoring. If \code{NULL} (default), no
##'   additional censoring is applied.
##' @param ... Further arguments passed to \code{\link{sim_recurrent_list}}.
##'
##' @return A data frame in counting-process format (one row per event interval
##'   per subject) with columns:
##'   \item{id}{Subject identifier.}
##'   \item{start, entry}{Interval start time.}
##'   \item{stop, time}{Interval end time (event or censoring time).}
##'   \item{status}{Event type at \code{stop}: \code{1} or \code{2} for a
##'     recurrent event of the corresponding type, \code{0} for censoring.}
##'   \item{death}{Indicator for a terminal event (\code{1}) or censoring/survival
##'     (\code{0}).}
##'   Attributes \code{"cumhaz"}, \code{"death.cumhaz"}, \code{"rr"}, and
##'   \code{"rd"} store the inputs used for simulation.
##'
##' @author Thomas Scheike
##'
##' @seealso \code{\link{sim_recurrent}}, \code{\link{sim_recurrent_list}},
##'   \code{\link{sim_recurrentTS}}
##'
##' @examples
##' data(CPH_HPN_CRBSI)
##' dr    <- CPH_HPN_CRBSI$terminal
##' base1 <- CPH_HPN_CRBSI$crbsi
##' base4 <- CPH_HPN_CRBSI$mechanical
##'
##' ## Single recurrent event type, with and without terminal event
##' rr <- sim_recurrent(5, base1)
##' dlist(rr, . ~ id, n = 0)
##'
##' rr <- sim_recurrent(5, base1, death.cumhaz = dr)
##' dlist(rr, . ~ id, n = 0)
##'
##' ## Verify that estimated rates recover the true baselines (increase n for precision)
##' rr <- sim_recurrent(100, base1, death.cumhaz = dr)
##' par(mfrow = c(1, 3))
##' mets:::showfitsim(causes = 1, rr, dr, base1, base1)
##'
##' ## Shared frailty across all processes
##' rr <- sim_recurrent(100, base1, death.cumhaz = dr, dependence = 1, var.z = 0.4)
##' dtable(rr, ~death + status)
##'
##' ## Two event types; second type uses the mechanical complication rate
##' set.seed(100)
##' rr <- sim_recurrentII(100, base1, base4, death.cumhaz = dr)
##' dtable(rr, ~death + status)
##' par(mfrow = c(2, 2))
##' mets:::showfitsim(causes = 2, rr, dr, base1, base4)
##'
##' ## Three event types and two causes of death via sim_recurrent_list
##' set.seed(100)
##' cumhaz <- list(base1, base1, base4)
##' drl    <- list(dr, base4)
##' rr     <- sim_recurrent_list(100, cumhaz, death.cumhaz = drl, dependence = 0)
##' dtable(rr, ~death + status)
##' mets:::showfitsimList(rr, cumhaz, drl)
##'
##' @export
##' @aliases sim_recurrentII sim_recurrent_list
sim_recurrentII <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,r1=NULL,r2=NULL,rd=NULL,rc=NULL,dependence=0,var.z=1,
			   cor.mat=NULL,cens=NULL,gap.time=FALSE,max.recurrent=100,...) 
{# {{{
cumhazL <- list(cumhaz,cumhaz2)
rr <- cbind(r1,r2)
if (!is.null(death.cumhaz)) death.cumhaz <- list(death.cumhaz)
if (!is.null(r1)) {
	if (!is.null(r2)) r2 <- rep(1,length(r1))
	rr <- cbind(r1,r2)
}
data <-     sim_recurrent_list(n,cumhazL,death.cumhaz=death.cumhaz,rr=rr,
		     rd=rd,rc=rc,dependence=dependence,var.z=var.z,
		     cor.mat=cor.mat,cens=cens,gap.time=gap.time,
		     max.recurrent=max.recurrent,...)
return(data)
}# }}}

##' Simulate recurrent events with a single event type and a terminal event
##'
##' A convenience wrapper around \code{\link{sim_recurrentII}} for the common
##' case of a single recurrent event type. Frailty and censoring options are
##' passed through to \code{\link{sim_recurrent_list}}.
##'
##' @inherit sim_recurrentII examples author
##' @param n Number of subjects to simulate.
##' @param cumhaz Two-column matrix \code{(time, cumhaz)} giving the cumulative
##'   hazard of the recurrent event.
##' @param death.cumhaz Two-column matrix \code{(time, cumhaz)} giving the
##'   cumulative hazard of the terminal event. If \code{NULL}, no terminal event
##'   is included.
##' @param r1 Optional numeric vector of length \code{n} of subject-specific
##'   relative risks for the recurrent event.
##' @param rd Optional numeric vector of length \code{n} of subject-specific
##'   relative risks for the terminal event.
##' @param rc Optional numeric vector of length \code{n} of subject-specific
##'   multipliers for the exponential censoring rate.
##' @param ... Further arguments passed to \code{\link{sim_recurrent_list}},
##'   including \code{dependence}, \code{var.z}, \code{gap.time}, and
##'   \code{max.recurrent}.
##' @export
sim_recurrent <- function(n,cumhaz,death.cumhaz=NULL,r1=NULL,rd=NULL,rc=NULL,...) 
{# {{{
## wrapper for sim_recurrentII without type-2 events
if (!is.null(death.cumhaz)) death.cumhaz <- list(death.cumhaz)
if (!is.null(r1)) r1 <- as.matrix(r1,ncol=1)
if (!is.null(rd)) rd <- as.matrix(rd,ncol=1)

rr <- sim_recurrent_list(n,list(cumhaz),death.cumhaz=death.cumhaz,rr=r1,rd=rd,rc=rc,...)
return(rr)
}# }}}

##' @export
sim_recurrent_list <- function(n,cumhaz,death.cumhaz=NULL,rr=NULL,rd=NULL,rc=NULL,zzr=NULL,zzd=NULL,
  gap.time=FALSE,max.recurrent=100,dhaz=NULL,haz2=NULL,dependence=0,var.z=1,cor.mat=NULL,cens=NULL,extend=TRUE,...) 
{# {{{
  status <- fdeath <- dtime <- NULL # to avoid R-check 

  if (dependence==0) { # {{{
	  z <- z1 <- z2 <- zd <- rep(1,n) 
	  if (is.null(zzr)) zzr <- matrix(z,n,length(cumhaz)) 
	  if (!is.null(death.cumhaz))
	  if (is.null(zzd)) zzd <- matrix(zd,n,length(death.cumhaz))
     } else if (dependence==1) {
	  z <- rgamma(n,1/var.z[1])*var.z[1]
	  z1 <- z; z2 <- z; zd <- z
	  if (is.null(zzr)) zzr <- matrix(z,n,length(cumhaz)) 
	  if (!is.null(death.cumhaz))
	  if (is.null(zzd)) zzd <- matrix(zd,n,length(death.cumhaz))
     } else if (dependence==4) {
	      zz <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- zz; z2 <- zz; zd <- rep(1,n) 
	  if (is.null(zzr)) zzr <- matrix(zz,n,length(cumhaz)) 
	  if (!is.null(death.cumhaz))
	  if (is.null(zzd)) zzd <- matrix(zd,n,length(death.cumhaz))
     }    else stop("dependence 0,1,4"); # }}}

   if (is.null(rr)) rr <- matrix(1,n,length(cumhaz))
   if (!is.null(death.cumhaz))
   if (is.null(rd)) rd <- matrix(1,n,length(death.cumhaz))
   if (is.null(rc)) rc <- rep(1,n)

  ## extend cumulative for death to full range  of cause 1
  if (!is.null(death.cumhaz)) {
     out <- extendCums(c(cumhaz,death.cumhaz),NULL,extend=extend)
     l <- length(cumhaz)
     ld <- length(death.cumhaz)
     cumhaz <- out[1:l]
     cumhazd <- out[(l+1):(l+ld)]
  } else {
     if (length(cumhaz)>1) {
	     out <- extendCums(cumhaz,NULL,extend=extend)
	     l <- length(cumhaz)
	     cumhaz <- out[1:l]
      }
  }
  max.time <- tail(cumhaz[[1]][,1],1)

  rrz <- rr*zzr
  tall <- rchazl(cumhaz,rrz) 
  tall$id <- 1:n

### death time simulated
  if (!is.null(death.cumhaz)) {# {{{
          rrzd <- rd*zzd
	  timed   <- rchazl(cumhazd,rrzd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }
# }}}


### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  i <- 0; 
  while (any((tt$time<tt$dtime) & (tt$status!=0) & (i <= max.recurrent))) {
	  i <- i+1
	  still <- subset(tt,time<dtime & status!=0)
	  nn <- nrow(still)
	  tt <- rchazl(cumhaz,rrz[still$id,,drop=FALSE],entry=(1-gap.time)*still$time) 
	  if (i==max.recurrent) {
		  tt$time <- max.time; tt$status <- 0
	  }
	  if (gap.time) {
		  tt$entry <- still$time
		  tt$time  <- tt$time+still$time
		  if (i==max.recurrent) { tt$time <- max.time; tt$status <- 0}
	  }
          ###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt[1:nn,],row.names=NULL)
  } 

  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time
  tall$rr <- NULL

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"rr") <- rr
  attr(tall,"rd") <- rd
  attr(tall,"z") <- zzr
  attr(tall,"zd") <- zzd

  return(tall)
  }# }}}

showfitsimList <- function(rr,cumhaz,dr,...)  {# {{{
 colp <- max(length(cumhaz),length(dr))
 par(mfrow=c(2,colp))
 if (!is.null(cumhaz)) {
 for (i in 1:length(cumhaz)) {
	 pp <- phreg(Surv(entry,time,status==i)~+1,rr)
	 plot(pp,main=paste("status=",i),...); 
	 lines(cumhaz[[i]],col=2,lwd=2)
 }
 }
 if (!is.null(dr)) {
 for (i in 1:length(dr)) {
	 pp <- phreg(Surv(entry,time,death==i)~+1,rr)
	 plot(pp,main=paste("death=",i),...); 
	 lines(dr[[i]],col=2,lwd=2)
 }
 }
}
# }}}

showfitsim <- function(causes=2,rr,dr,base1,base4,which=1:3) 
{# {{{
if (1 %in% which) {
  drr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)
  basehazplot.phreg(drr,ylim=c(0,8))
  lines(dr,col=2)
}
###
if (2 %in% which) {
  xrr <- phreg(Surv(entry,time,status==1)~cluster(id),data=rr)
  basehazplot.phreg(xrr,add=TRUE)
###  basehazplot.phreg(xrr)
  lines(base1,col=2)
  if (causes>=2) {
	  xrr2 <- phreg(Surv(entry,time,status==2)~cluster(id),data=rr)
	  basehazplot.phreg(xrr2,add=TRUE)
	  lines(base4,col=2)
  }
  }
if (3 %in% which) {
  meanr1 <-   recurrentMarginalPhreg(xrr,drr)
  basehazplot.phreg(meanr1,se=TRUE)
  if (causes>=2) {
	  meanr2 <-   recurrentMarginalPhreg(xrr2,drr)
	  basehazplot.phreg(meanr2,se=TRUE,add=TRUE,col=2)
  }
}
}# }}}

##' Simulate recurrent events from a two-stage Cox or Ghosh-Lin model
##'
##' Simulates recurrent event data from a fitted two-stage model, where the
##' recurrent event process and the terminal event are each described by a
##' separate fitted model. The recurrent event model may be either a Cox
##' proportional hazards model (\code{phreg}) or a Ghosh-Lin marginal rate model
##' (\code{recreg}); the terminal event model must be a Cox model (\code{phreg}).
##'
##' Covariates are drawn by bootstrap from \code{data} (if supplied), and
##' subject-specific relative risks and strata are derived from the fitted
##' model objects. Stratified baselines are fully supported. When
##' \code{dependence} is \code{NULL} (default), the simulation uses the
##' two-stage structure from \code{\link{sim_GLcox}}; setting \code{dependence}
##' to an integer falls back to \code{\link{sim_recurrent_list}} with the
##' corresponding frailty model.
##'
##' @param cox1 A fitted \code{phreg} object for the recurrent event rate, or a
##'   fitted \code{recreg} (Ghosh-Lin) object. The model type is detected
##'   automatically from the class.
##' @param coxd A fitted \code{phreg} object for the terminal event. May be
##'   \code{NULL} if no terminal event is modelled.
##' @param n Number of subjects to simulate. Default is \code{1}.
##' @param data The data frame on which \code{cox1} and \code{coxd} were fitted,
##'   used to draw covariate values for the simulated subjects. If \code{NULL},
##'   covariates must be supplied via \code{r1}, \code{rd}, \code{strata1}, and
##'   \code{stratad}.
##' @param type Simulation type: \code{"default"} (auto-detected from class of
##'   \code{cox1}), \code{"cox-cox"}, or \code{"gl-cox"}.
##' @param id Name of the subject identifier variable in \code{data}. Default is
##'   \code{"id"}.
##' @param varz Variance of the frailty distribution in the two-stage model.
##'   Default is \code{1}.
##' @param share Proportion of the shared frailty assigned to the recurrent event
##'   process in the partial-sharing model. Default is \code{1}.
##' @param cens Rate of exponential censoring. Default is \code{0.001}.
##' @param scale1 Scalar multiplier for the baseline cumulative hazard of the
##'   recurrent event process. Default is \code{1}.
##' @param scaled Scalar multiplier for the baseline cumulative hazard of the
##'   terminal event. Default is \code{1}.
##' @param dependence If non-\code{NULL}, falls back to
##'   \code{\link{sim_recurrent_list}} with this frailty structure (see
##'   \code{\link{sim_recurrentII}} for valid values). Default is \code{NULL}.
##' @param r1 Optional numeric vector of length \code{n} of subject-specific
##'   relative risks for the recurrent event, used when \code{data = NULL}.
##' @param rd Optional numeric vector of length \code{n} of subject-specific
##'   relative risks for the terminal event, used when \code{data = NULL}.
##' @param rc Optional numeric vector of length \code{n} of subject-specific
##'   censoring rate multipliers.
##' @param strata1 Optional integer vector of length \code{n} specifying the
##'   stratum index (0-based) for the recurrent event model, used when
##'   \code{data = NULL}.
##' @param stratad Optional integer vector of length \code{n} specifying the
##'   stratum index (0-based) for the terminal event model, used when
##'   \code{data = NULL}.
##' @param death.code Integer status code used for the terminal event in the
##'   output \code{status} column. Default is \code{3}.
##' @param ... Further arguments passed to \code{\link{sim_GLcox}}, including
##'   \code{nmin} and \code{nmax} for the linear approximation grid.
##'
##' @return A data frame in counting-process format with one row per event
##'   interval per subject. Column names match those in the original model
##'   formula (entry, exit, and status variables). Additional columns include
##'   \code{id} and the covariates drawn from \code{data} (if supplied). The
##'   terminal event is coded as \code{death.code} in the status variable;
##'   recurrent events are coded as \code{1}.
##'
##' @author Thomas Scheike
##'
##' @seealso \code{\link{recurrent_marginal}}, \code{\link{sim_recurrent_list}},
##'   \code{\link{sim_GLcox}}
##'
##' @references
##' Scheike, T. H. (2026). Two-stage recurrent events random effects models.
##' \emph{Lifetime Data Analysis}.
##'
##' @examples
##' data(hfactioncpx12)
##' hf    <- hfactioncpx12
##' hf$x  <- as.numeric(hf$treatment)
##' n     <- 100
##'
##' ## Cox-Cox two-stage model
##' xr <- phreg(Surv(entry, time, status == 1) ~ x + cluster(id), data = hf)
##' dr <- phreg(Surv(entry, time, status == 2) ~ x + cluster(id), data = hf)
##' simcoxcox <- sim_recurrent_ts(xr, dr, n = n, data = hf, death.code = 2)
##'
##' ## Ghosh-Lin/Cox two-stage model
##' recGL  <- recreg(Event(entry, time, status) ~ x + cluster(id), hf, death.code = 2)
##' simglcox <- sim_recurrent_ts(recGL, dr, n = n, data = hf, death.code = 2)
##'
#' @export sim_recurrent_ts
sim_recurrent_ts <- function(cox1,coxd=NULL,
             n=1, data=NULL,type=c("default","cox-cox","gl-cox"),
             id="id",varz=1,share=1,cens=0.001,
      scale1=1,scaled=1,dependence=NULL,
 r1=NULL,rd=NULL,rc=NULL,strata1=NULL,stratad=NULL,death.code=3,
                          ...) {# {{{
## exp censoring default
statusD <- death <- NULL

if (type[1]=="default" & inherits(cox1,"recreg")) type <- "gl-cox" 
if (type[1]=="default" & inherits(cox1,"phreg")) type <- "cox-cox" 
if (type[1]=="cox-cox") type <- 3 else type <- 2

if (!is.null(data)) {
   coxs <- list(cox1,coxd)
   rrdata <- draw_phregs(coxs,n,data=data)
   rr1 <- rrdata$rr[,1]
   rstrata1 <- rrdata$strata[,1]

   if (!is.null(coxd)) {
      rrd <- rrdata$rr[,2]
      rstratad <- rrdata$strata[,2]
   }
} else { 
	data <- c()
	rr1 <- rep(1,n)
	rrd <- rep(1,n)
	rstrata1 <- rep(0,n)
	rstratad <- rep(0,n)
}

if (is.null(rc)) rc <- rep(1,n)

if (!is.null(coxd))    {
    coxd$cumhaz <- scalecumhaz(coxd$cumhaz,scaled); 
    LamD <- basecumhaz(coxd,only=1)
} else LamD <- NULL
cox1$cumhaz <-   scalecumhaz(cox1$cumhaz,scale1); 
Lam1 <- basecumhaz(cox1,only=1)

### when data not given then use r1,rd, strata1,stratad
if (is.null(r1)) r1 <- rr1
if (is.null(rd)) rd <- rrd
if (is.null(strata1)) strata1 <- rstrata1
if (is.null(stratad)) stratad <- rstratad


if (is.null(dependence) & (!is.null(LamD))) {
 ## go through combined strata for the two models 
strat1d <- mystrata(data.frame(strata1,stratad))
rrs <- c()
for (j in 1:attr(strat1d,"nlevel")) {
     r1i <- which(strat1d==j)
     strata1ss <- strata1[r1i[1]]+1
     stratadss <- stratad[r1i[1]]+1
     Lam1s <- Lam1[[strata1ss]]
     LamDs <- LamD[[stratadss]]
     ns <- length(r1i)
     rrss <- sim_GLcox(ns,Lam1s,LamDs,var.z=varz,r1=r1[r1i],rd=rd[r1i],rc=rc[r1i],
		model="twostage",cens=cens,type=type,share=share,...)
    rrss$ids <- rrss$id
    rrss$id <- r1i[rrss$id+1]
    rrs <- rbind(rrs,rrss)
}			  
rrs <- dtransform(rrs,statusD=death.code,statusD==3)
} else { 
if (is.null(dependence)) dependence <- 0
if (!is.null(LamD)) 
rrs <- sim_recurrent_list(n,Lam1,death.cumhaz=LamD,rr=matrix(r1,ncol=1),rd=matrix(rd,ncol=1),rc=rc,cens=cens,var.z=varz,dependence=dependence)
else rrs <- sim_recurrent_list(n,list(Lam1),rr=matrix(r1,ncol=1),rc=rc,cens=cens,var.z=varz,dependence=dependence)
rrs$Z <- attr(rrs,"z")[rrs$id+1]
rrs$statusD <- rrs$status
if (!is.null(LamD))  {
   rrs <- dtransform(rrs,statusD=death.code,death==1)
}
}

## add correct names to entry,time,status if possible
if (inherits(cox1,"phreg")) {
   varsY <- all.vars(update(drop.specials(cox1$formula,"cluster"),.~1)) 
   rrs[,varsY] <- cbind(rrs$start,rrs$stop,rrs$statusD)
   rrs <- dkeep(rrs,c(varsY,"id"))
}

## add covariates, 
if (!is.null(data)) rrs <- cbind(rrs,rrdata$data[rrs$id,])

return(rrs)
}
# }}}

sim_RecurrentIIHist <- function(n,cumhaz,death.cumhaz,cens=NULL,rr=NULL,rc=NULL,rd=NULL,
	    max.recurrent=100,dependence=0,var.z=0.22,cor.mat=NULL,
	    HistN1=~I(Nt^.5),HistD=~I(Nt^.5),HistN1.beta=c(1.0),HistD.beta=c(1.0),...) 
{# {{{

  ctime <- fdeath <- dtime <- NULL # to avoid R-check 
  status <- dhaz <- NULL; dhaz2 <- NULL

  if (dependence==0) { z <- z1 <- zc <- zd  <-  rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
###	      z <- exp(rnorm(n,1)*var.z[1]^.5)
	      z1 <- z; z2 <- z; zd <- z
	      if (!is.null(cor.mat)) { zd <- rep(1,n); }
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
###	      print(summary(z))
###	      print(cor(z))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
###	      print(summary(z))
###	      print(cor(z))
      } else stop("dependence 0-3"); # }}}

  ## covariate adjustment 
  if (is.null(rr))  rr <- z1; 
  if (is.null(rc))  rc <- zc; 
  if (is.null(rd))  rd <- zd; 
  if (length(rr)!=n) rr <- rep(rr[1],n)
  if (length(rc)!=n) rc <- rep(rc[1],n)
  if (length(rd)!=n) rd <- rep(rd[1],n)

  ll <- nrow(cumhaz)
  ### extend of cumulatives
  cumhaz <- rbind(c(0,0),cumhaz)
  death.cumhaz <- rbind(c(0,0),death.cumhaz)

  if (!is.null(cens)) {
	  if (is.matrix(cens))  {
             out <- extendCums(list(cumhaz,death.cumhaz,cens),NULL)
   	     cumcens <- out$cum3
	  } else {
             out <- extendCums(list(cumhaz,death.cumhaz),NULL)
	  }
  } else {
     out <- extendCums(list(cumhaz,death.cumhaz),NULL)
  }
  cumhaz <- out$cum1
  cumhazd <- out$cum2
  max.time <- tail(cumhaz[,1],1)

  ### draw censorting times 
  if (!is.null(cens)) { 
	  if (is.matrix(cens)) sdata <- rchaz(cens,rc) else 
          sdata <- data.frame(entry=0,time=pmin(max.time,rexp(n)/(c(rc)*cens)),
			      status=0,id=1:n)
  } else sdata <- data.frame(entry=0,time=rep(max.time,n),status=0,id=1:n)
  sdata$ctime <- sdata$time
  sdata$Nt <- 0

  i <- 0; 
  tall <- c()
  still <- tt1 <-  sdata
  ## start at 0
  tt1$time <- still$time <- 0
  while (any((still$time<still$ctime) &  (i < max.recurrent))) { ## {{{
	  still$Nt <- i
	  nn <- nrow(still)
	  z1r <- rr[still$id]
	  zdr <- rd[still$id]

	  m1 <- model.matrix(HistN1,still)[,-1,drop=FALSE]
	  md <- model.matrix(HistD,still)[,-1,drop=FALSE]
	  r1h <- c(exp(m1 %*% HistN1.beta))
	  rdh <- c(exp(md %*% HistN1.beta))

          tt1 <- rcrisk(cumhaz,cumhazd,z1r*r1h,zdr*rdh,entry=still$time)
	  tt1 <- cbind(tt1,dkeep(still,~id+ctime),row.names=NULL)
	  tt1 <- dtransform(tt1,status=0,time>ctime)
	  tt1 <- dtransform(tt1,time=ctime,time>ctime)
	  tt1$Nt <- i
          ## remove dead from tt1 to get those still there 
          deadid <- (tt1$status %in% c(0,2))
	  ### those that are still under risk 
	  still <- tt1[!deadid,,drop=FALSE]
	  ## also keep only those before max.time
          still <- subset(still,still$time<max.time)
	  tall <- rbind(tall,tt1,row.names=NULL)
	  i <- i+1
  }  # }}}

  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time
  tall <- dkeep(tall,~id+entry+time+status+ctime+Nt+start+stop)

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cens.cumhaz") <- cens

  return(tall)
  }# }}}

##' Simulate recurrent events from a two-stage model with structured gamma frailties
##'
##' Simulates recurrent event data with two event types and a terminal event,
##' using a parametric two-stage frailty model. The construction ensures that the
##' marginal rates are approximately correct: conditional on survival,
##' \eqn{E(dN_j \mid D > t) \approx} \code{cumhazj}, and the hazard of death
##' equals \code{death.cumhaz}.
##'
##' The frailty structure uses three gamma random variables \eqn{Z_{d1}},
##' \eqn{Z_{d2}}, \eqn{Z_{12}} to induce dependence:
##' \deqn{Z_\text{death} = Z_{d1} + Z_{d2}, \quad
##'       Z_1 = Z_{d1}^{\nu_1} Z_{12}, \quad
##'       Z_2 = Z_{d2}^{\nu_2} Z_{12}^{\nu_3}.}
##' The parameters \code{share1} and \code{vargamD} control how the death frailty
##' splits between the two components, and \code{vargam12} controls the shared
##' recurrent-event frailty. Setting \eqn{\nu = (1,1,1)} with \code{share1 = 0.5}
##' gives a symmetric structure; varying \eqn{\nu} allows asymmetric dependence.
##'
##' @param n Number of subjects to simulate.
##' @param cumhaz Two-column matrix \code{(time, cumhaz)} giving the target
##'   marginal cumulative rate of the first recurrent event type.
##' @param cumhaz2 Two-column matrix \code{(time, cumhaz)} giving the target
##'   marginal cumulative rate of the second recurrent event type.
##' @param death.cumhaz Two-column matrix \code{(time, cumhaz)} giving the
##'   cumulative hazard of the terminal event.
##' @param nu Numeric vector of length 3: the powers \eqn{(\nu_1, \nu_2, \nu_3)}
##'   applied to the frailty components (see Details). Must satisfy
##'   \eqn{\nu_j > -1/\text{shape}}. Default is \code{rep(1, 3)}.
##' @param share1 Proportion of the total death frailty variance assigned to the
##'   first component \eqn{Z_{d1}}. The remainder goes to \eqn{Z_{d2}}. Must be
##'   in \eqn{(0, 1)}. Default is \code{0.3}.
##' @param vargamD Total variance of the death frailty \eqn{Z_\text{death}}.
##'   Default is \code{2}.
##' @param vargam12 Variance of the shared recurrent-event frailty \eqn{Z_{12}}.
##'   Default is \code{0.5}.
##' @param gap.time Logical. If \code{TRUE}, event times are drawn as gap times
##'   rather than calendar times. Default is \code{FALSE}.
##' @param max.recurrent Maximum number of recurrent events per subject. Default
##'   is \code{100}.
##' @param cens Rate of exponential censoring. If \code{NULL} (default), no
##'   administrative censoring is applied.
##' @param ... Further arguments passed to lower-level simulation functions.
##'
##' @return A data frame in counting-process format (one row per event interval
##'   per subject) with columns \code{id}, \code{start}, \code{stop},
##'   \code{entry}, \code{time}, \code{status}, and \code{death}. Attributes
##'   store the (possibly adjusted) cumulative hazards used in simulation and
##'   the frailty parameters.
##'
##' @author Thomas Scheike
##'
##' @seealso \code{\link{sim_recurrentII}}, \code{\link{sim_recurrent_ts}}
##'
##' @examples
##' data(CPH_HPN_CRBSI)
##' dr    <- CPH_HPN_CRBSI$terminal
##' base1 <- CPH_HPN_CRBSI$crbsi
##' base4 <- CPH_HPN_CRBSI$mechanical
##'
##' rr <- sim_recurrentTS(1000, base1, base4, death.cumhaz = dr)
##' dtable(rr, ~death + status)
##' mets:::showfitsim(causes = 2, rr, dr, base1, base4)
##'
##' @export
sim_recurrentTS <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,
		    nu=rep(1,3),share1=0.3,vargamD=2,vargam12=0.5,
		    gap.time=FALSE,max.recurrent=100,cens=NULL,...) 
{# {{{

k <- 1
nu1 <- nu[1]; nu2 <- nu[2]; nu3 <- nu[3]
###nu1 <- 1; nu2 <- 1; nu3 <- 0.4
share2 <- (1-share1)
vargam <- vargamD
vargam12 <- 0.5
agam1 <- share1/vargam
agam2 <- share2/vargam
betagam=1/vargam
gamma1  <- rep(rgamma(n,agam1)*vargam,each=k)
gamma2  <- rep(rgamma(n,agam2)*vargam,each=k)
agam12 <- 1/vargam12
betagam12 <- 1/vargam12
gamma12 <- rep(rgamma(n,agam12)*vargam12,each=k)
agamD <- agam1+agam2
z1 <- (gamma1^nu1)*gamma12
z2 <- (gamma2^nu2)*gamma12^nu3
gamD <- gamma1+gamma2
zd <- gamD
egamma12nu3 <- (gamma(agam12+nu3)/gamma(agam12))*1/(betagam12)^nu3
zs <- cbind(z1,z2,zd)

 status <- fdeath <- dtime <- NULL # to avoid R-check 
 dhaz <- haz2 <- dhaz <- NULL

 ll <- nrow(cumhaz)
 max.time <- tail(cumhaz[,1],1)

 ################################################################
 ### approximate hazards to make marginals fit (approximately)
 ################################################################
 ## step-size set to one and range to that of base1
 base1 <- predictCumhaz(rbind(0,as.matrix(cumhaz)),1:round(tail(cumhaz[,1],1)) )
 if (!is.null(death.cumhaz)) 
     death.cumhaz <- predictCumhaz(rbind(0,as.matrix(death.cumhaz)),base1[,1])

 orig.death <- death.cumhaz
 dbase1 <- death.cumhaz
 gt <- exp(vargam*dbase1[,2]) 
 dtt <- diff(c(0,dbase1[,1]))
 lams <- (diff(c(0,dbase1[,2]))/dtt)*gt
 death.cumhaz <- cbind(dbase1[,1],cumsum(dtt*lams))

 dbase1 <- cpred(rbind(c(0,0),death.cumhaz),base1[,1])[,2]
 dtt <- diff(c(0,base1[,1]))
 gt <- (gamma(agam1+nu1)/gamma(agam1))*(1/(betagam+dbase1))^nu1
 lams <- (diff(c(0,base1[,2]))/dtt)*(1/gt)
 cumhaz <- cbind(base1[,1],cumsum(dtt*lams))

 base1 <- cumhaz2
 dbase1 <- cpred(rbind(c(0,0),death.cumhaz),base1[,1])[,2]
 dtt <- diff(c(0,base1[,1]))
 gt <- (gamma(agam2+nu2)/gamma(agam2))*(1/(betagam+dbase1))^nu2
 lams <-(1/egamma12nu3)*(diff(c(0,base1[,2]))/dtt)*(1/gt)
 cumhaz2 <- cbind(base1[,1],cumsum(dtt*lams))

 cumhaz <- rbind(c(0,0),cumhaz)
 cumhaz2 <- rbind(c(0,0),cumhaz2)
 death.cumhaz <- rbind(c(0,0),death.cumhaz)

## range max of cumhaz and cumhaz2 
  out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz),NULL)
  cumhaz <- out$cum1
  cumhaz2 <- out$cum2
  cumhazd <- out$cum3
  max.time <- tail(cumhaz[,1],1)

### recurrent first time
  tall1 <- rchaz(cumhaz,rr=z1)
  tall2 <- rchaz(cumhaz2,rr=z2)
  tall <- tall1 
  tall$status <- ifelse(tall1$time<tall2$time,tall1$status,2*tall2$status)
  tall$time <- ifelse(tall1$time<tall2$time,tall1$time,tall2$time)
  tall$id <- 1:n
  tall$rr2 <- tall2$rr
### death time simulated
  if (!is.null(death.cumhaz)) {# {{{
	  timed   <- rchaz(cumhazd,rr=zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }# }}}

  ### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  ### setting aside memory 
  tt1 <- tt2 <- tt
  i <- 1; 
  while (any((tt$time<tt$dtime) & (tt$status!=0) & (i < max.recurrent))) {
	  i <- i+1
	  still <- subset(tt,time<dtime & status!=0)
	  nn <- nrow(still)
          tt1 <- rchaz(cumhaz,rr=z1[still$id],entry=still$time)
          tt2 <- rchaz(cumhaz2,rr=z2[still$id],entry=still$time)
	  tt <- tt1
          tt$status <- ifelse(tt1$time<=tt2$time,tt1$status,2*tt2$status)
          tt$time <-   ifelse(tt1$time<=tt2$time,tt1$time,tt2$time)
	  tt$rr2 <- tt2$rr
          ###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt[1:nn,],row.names=NULL)
  }
  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"zs") <- zs

  attr(tall,"gamma.death") <- c(agam1,agam2,betagam,vargamD)
  attr(tall,"gamma.N12") <-   c(agam12,betagam12,vargam12)

  return(tall)
  }# }}}

##' Compute cumulative event counts as time-dependent covariates
##'
##' For each subject and each row of \code{data}, counts the number of prior
##' events of specified types in the recurrent event history. The resulting
##' count columns can be used as time-dependent covariates in subsequent models,
##' e.g. to capture event-history dependence in the recurrent event rate.
##'
##' When \code{lag = TRUE} (default), the count at each row reflects events that
##' occurred strictly before the current time point (i.e. \eqn{N(t-)}), making
##' it suitable as a left-continuous covariate in counting-process models.
##'
##' @param data A data frame in counting-process format, with one row per event
##'   interval per subject.
##' @param status Name of the column containing event status codes. Default is
##'   \code{"status"}.
##' @param id Name of the column containing subject identifiers. Default is
##'   \code{"id"}.
##' @param types Integer vector of status codes to count. Each value in
##'   \code{types} generates one new count column (when \code{multitype = FALSE})
##'   or contributes to a single combined count (when \code{multitype = TRUE}).
##'   Default is \code{1}.
##' @param names.count Prefix for the names of the new count columns. The
##'   status code is appended, e.g. \code{"Count1"}, \code{"Count2"}.
##'   Default is \code{"Count"}.
##' @param lag Logical. If \code{TRUE} (default), the count at each row is the
##'   number of events strictly before the current time (\eqn{N(t-)}). If
##'   \code{FALSE}, events at the current time are included (\eqn{N(t)}).
##' @param multitype Logical. If \code{TRUE}, events with status in \code{types}
##'   are aggregated into a single count column, optionally weighted by
##'   \code{marks}. If \code{FALSE} (default), a separate count column is created
##'   for each value in \code{types}.
##' @param marks Optional numeric vector of weights applied to events when
##'   \code{multitype = TRUE}. If \code{NULL} (default), each event has weight 1.
##'
##' @return The input data frame \code{data} with one or more new integer columns
##'   appended. With \code{multitype = FALSE}, columns are named
##'   \code{paste0(names.count, k)} for each \code{k} in \code{types}; with
##'   \code{multitype = TRUE}, a single column named
##'   \code{paste0(names.count, types[1])} is added. An internal bookkeeping
##'   column \code{lbnr__id} is also added.
##'
##' @author Thomas Scheike
##'
##' @examples
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##' dtable(hf, ~status)
##'
##' ## Separate counts for event types 1 and 2
##' rr <- count_history(hf, types = 1:2, id = "id", status = "status")
##' dtable(rr, ~"Count*" + status, level = 1)
##'
##' @export
count_history <- function(data,status="status",id="id",types=1,names.count="Count",lag=TRUE,multitype=FALSE,marks=NULL)
{# {{{
stat <- data[,status]

## also allowing marks  when multitype=TRUE
if (multitype & is.null(marks))  marks <- rep(1,nrow(data))

clusters <- data[,id]
if (is.numeric(clusters)) {
      clusters <- fast.approx(unique(clusters), clusters) - 1
      max.clust <- max(clusters)
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

data[,"lbnr__id"] <- cumsumstrata(rep(1,nrow(data)),clusters,max.clust+1) 
if (!multitype) {
for (i in types)  {
if (lag==TRUE)
data[,paste(names.count,i,sep="")] <- 
   cumsumidstratasum((stat==i),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
   else 
data[,paste(names.count,i,sep="")] <- 
   cumsumidstratasum((stat==i),rep(0,nrow(data)),1,clusters,max.clust+1)$sum 
}
} else {
if (lag==TRUE)
data[,paste(names.count,types[1],sep="")] <- 
   cumsumidstratasum((stat %in% types)*marks,rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
   else 
data[,paste(names.count,types[1],sep="")] <- 
   cumsumidstratasum((stat %in% types)*marks,rep(0,nrow(data)),1,clusters,max.clust+1)$sum 
}

return(data)
}# }}}

count_historyVar <- function(data,var="status",id="id",names.count="Count",lag=TRUE)
{# {{{
vvar <- data[,var]

clusters <- data[,id]
if (is.numeric(clusters)) {
      clusters <- fast.approx(unique(clusters), clusters) - 1
      max.clust <- max(clusters)
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

data[,"lbnr__id"] <- cumsumstrata(rep(1,nrow(data)),clusters,max.clust+1) 
if (lag==TRUE)
data[,names.count] <- 
   cumsumidstratasum(vvar,rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
   else 
data[,names.count] <- 
   cumsumidstratasum(vvar,rep(0,nrow(data)),1,clusters,max.clust+1)$sum 

return(data)
}# }}}

##' Estimate the probability of exceeding k recurrent events by time t
##'
##' Estimates \eqn{P(N(t) \geq k)} as a function of time \eqn{t}, for a range of
##' thresholds \eqn{k}, in the presence of a terminal event (death). The estimator
##' is based on the cumulative incidence of "reaching \eqn{k} events", treating
##' death as a competing risk. Confidence intervals are computed on the log or
##' plain scale.
##'
##' For each threshold \eqn{k} in \code{exceed}, the function identifies the
##' first time each subject reaches \eqn{k} events, then fits a competing risks
##' model (\code{\link{cif}}) with "reaching \eqn{k} events" as the event of
##' interest and death as the competing event. Strata are supported. When
##' \code{marks} is \code{NULL}, each event contributes equally; otherwise events
##' are weighted by their mark values before cumulative counts are formed.
##'
##' @param formula A formula with an \code{\link{Event}} response giving the
##'   exit time and status (and optionally entry time). The right-hand side may
##'   include \code{cluster()} and \code{strata()}.
##' @param data A data frame containing all variables in \code{formula}.
##' @param cause Integer code for the recurrent event of interest. Default is
##'   \code{1}.
##' @param death.code Integer code for the terminal event. Default is \code{2}.
##' @param cens.code Integer code for censoring. Default is \code{0}.
##' @param exceed Integer vector of thresholds \eqn{k} to evaluate. If
##'   \code{NULL} (default), all observed cumulative counts are used.
##' @param marks Optional numeric vector of event weights. If non-\code{NULL},
##'   cumulative counts are formed as weighted sums of events rather than simple
##'   counts. Must have the same length as \code{nrow(data)}.
##' @param all.cifs Logical. If \code{TRUE}, the fitted \code{cif} object for
##'   each threshold is returned in \code{cif.exceed}. Default is \code{FALSE}.
##' @param return.data Logical. If \code{TRUE}, the constructed dataset for each
##'   threshold is returned in \code{dataList}. Default is \code{FALSE}.
##' @param conf.type Type of confidence interval transformation: \code{"log"}
##'   (default) or \code{"plain"}.
##' @param level Confidence level. Default is \code{0.95}.
##' @param ... Further arguments passed to \code{\link{cif}}.
##'
##' @return An object of class \code{"exceed"} with the following components:
##'   \item{time}{Vector of evaluation time points.}
##'   \item{prob}{Array of dimension \code{(length(time), length(exceed) + 1,
##'     nstrata)} containing \eqn{P(N(t) \geq k)} for each threshold and stratum.
##'     The first column gives \eqn{P(N(t) < \text{exceed}[1])}.}
##'   \item{se.prob}{Standard errors of \code{prob}.}
##'   \item{lower, upper}{Pointwise confidence interval bounds.}
##'   \item{meanN}{Estimated mean number of events \eqn{E(N(t))} (single stratum
##'     only; \code{NULL} for stratified analyses).}
##'   \item{meanN2, varN}{Second moment and variance of \eqn{N(t)} (single
##'     stratum only).}
##'   \item{exceed}{Thresholds evaluated (excluding zero).}
##'   \item{cif.exceed}{List of fitted \code{cif} objects (if
##'     \code{all.cifs = TRUE}).}
##'   \item{dataList}{List of datasets for each threshold (if
##'     \code{return.data = TRUE}).}
##'   \item{nstrata, strata.levels, strata.name}{Stratification information.}
##'   Use \code{plot()} and \code{summary()} methods for visualisation and
##'   tabulation.
##'
##' @author Thomas Scheike
##'
##' @seealso \code{\link{recurrent_marginal}}, \code{\link{cif}}
##'
##' @references
##' Scheike, T. H., Eriksson, L., and Tribler, P. (2019). The mean, variance and
##' correlation for bivariate recurrent events with a terminal event.
##' \emph{Journal of the Royal Statistical Society, Series C}, 68(5).
##'
##' @examples
##' data(hfactioncpx12)
##' dtable(hfactioncpx12, ~status)
##'
##' oo <- prob_exceed_recurrent(Event(entry, time, status) ~ cluster(id),
##'                              hfactioncpx12, cause = 1, death.code = 2)
##' plot(oo)
##' summary(oo, times = c(1, 2, 5))
##' @export
prob_exceed_recurrent <- function(formula,data,cause=1,death.code=2,cens.code=0, exceed=NULL,marks=NULL,all.cifs=FALSE,
				  return.data=FALSE,conf.type=c("log","plain"),level=0.95,...)
{# {{{
    cifmets <- TRUE
    cl <- match.call()# {{{
    m <- match.call(expand.dots = TRUE)[1:3]
    des <- proc_design(
        formula,
        data = data,
        specials = c("marks","strata","offset", "weights", "cluster"),
        intercept = TRUE
    )
    Y <- des$y
    if (!inherits(Y, c("Event", "Surv"))) {
        stop("Expected a 'Surv' or 'Event'-object")
    }
    if (ncol(Y) == 2) {
        exit <- Y[, 1]
        entry <- rep(0, nrow(Y))
        status <- Y[, 2]
    } else {
        entry <- Y[, 1]
        exit <- Y[, 2]
        status <- Y[, 3]
    }
    X <- des$x
    des.weights <- des$weights
    des.offset  <- des$offset
    des.marks <- des$marks
    id      <- des$cluster
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    strata <- des$strata
    if (!is.null(strata))  {
      ns <- grep("strata",names(des$levels))
      strata.name  <-  names(des$levels)[ns]
    } else strata.name <- NULL
    ## no use of 
    pos.cluster <- pos.strata <- NULL

 ## take offset and weight first from formula, but then from arguments
  if (is.null(des.offset)) {
	  if (is.null(offset)) offset <- rep(0,length(exit)) 
  } else offset <- des.offset
  if (is.null(des.weights)) {
	  if (is.null(weights)) weights <- rep(1,length(exit)) 
  } else weights <- des.weights
  if (!is.null(des.marks) & is.null(marks))  marks <- des.marks
  marks.call <- marks
  if (is.null(marks)) marks <- rep(1,length(exit))
  ## }}}
  if (!is.null(strata)) {
	  nstrata <- nlevels(strata) 
          strata.levels <- levels(strata)
  } else { strata.levels <- NULL; nstrata <- 1}

 call.id <- id;
 conid <- construct_id(id,nrow(X),as.data=TRUE)
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
 data$id__  <- id

 times <- NULL
 statusD <- (status %in% death.code)*1
 statusE <- (status %in% cause)*1
 if (sum(statusE)==0) stop("none of type events")
 if (!is.null(strata) & !cifmets) stop("strata only for cifmets=TRUE\n")
 allvars <- all.vars(formula)
 if (is.null(times)) times <- sort(unique(exit[statusE==1]))

 countE <- cumsumstrata(statusE*marks,id,nid)
 if (is.null(exceed)) exceed <- sort(unique(countE))
 w0 <- which(exceed==0)
 if (length(w0)>=1) exceed <- exceed[-w0]

 if (is.null(marks.call)) {
    mc <- max(countE)+1
    idcount <- id*mc+countE
    idcount <- cumsumstrata(rep(1,length(idcount)),idcount,mc*(nid+1))
 } else {
    ex <- outer(countE,exceed,">=")[,1,]
    exC <- apply(ex,2,cumsumstrata,id,nid)
    idcount <- exC
 }

if (!cifmets) {
  ## can not handle start stop and id, so uses MG standard errors 
  pp <- as.formula(paste("Hist(entry=",allvars[1],",",allvars[2],",statN)~+1",sep=""))
  ###rhs <- update(formula,-1~.)
  ###form <- as.formula(update.formula(rhs,pp))
  form <- as.formula(pp)
} else {
###  pp <- as.formula(paste("Event(",allvars[1],",",allvars[2],",statN)~.",sep=""))
  pp <- as.formula(paste("Event(",allvars[2],",statN)~.",sep=""))
  rhs <- update(formula,-1~.)
  form <- as.formula(update.formula(rhs,pp))
}

cif.exceed <- NULL
dataList <- NULL
if (all.cifs) cif.exceed <- list() 
if (return.data) dataList <- list() 
se.probs <- probs <- lower <- upper <- array(0,c(length(times),length(exceed)+1,nstrata))
i <- 1
for (n1 in exceed) {# {{{
	i <- i+1
	### first time that get to n1
        if (is.null(marks.call)) keep <- (countE<n1 ) | (countE==n1 & idcount==1) else keep <- exC[,i-1]<=1 
###	keep <- (countE<n1 ) | (firsttoE==n1)

	### status, censoring, get to n1, or die
        statN <- rep(0,nrow(data))
        if (is.null(marks.call)) statN[countE==n1] <- 1 else statN[exC[,i-1]==1] <- 1
	statN[statusD==1] <- 2
	statN <- statN[keep]
	dataN <- data[keep,]
	dataN$statN <- statN

        idN <- revcumsumstrata(rep(1,nrow(dataN)),dataN$id__,nid)
        dataN <- subset(dataN,idN==1)
	if (return.data) dataList[[i-1]] <- dataN

   pN1 <-  suppressWarnings(cif(form,data=dataN))
   if (all.cifs) cif.exceed[[i-1]] <- pN1

###   lpN1 <- basecumhaz(pN1)
   where <- predictCumhaz(c(0,pN1$times),times)
###  where <- fast.approx(c(0,pN1$times),times,type="left")
###   print(c(length(pN1$mu),length(pN1$strata),pN1$nstrata))
   if (length(pN1$mu)>=1) {
   cifs <- vecAllStrata(pN1$mu,pN1$strata,pN1$nstrata)
   se.cifs <- vecAllStrata(pN1$se.mu,pN1$strata,pN1$nstrata)
   probs[,i,] <- rbind(0,cifs)[where,]
   se.probs[,i,] <- rbind(0,se.cifs)[where,]

   lu <- conftype(probs[,i,],se.probs[,i,],conf.type=conf.type[1],
         	 restrict=c("prob"),conf.int=level)
   lower[,i,] <- lu$lower
   upper[,i,] <- lu$upper
   } else {
   probs[,i,] <- se.probs[,i,] <- lower[,i,] <- upper[,i,] <- 0
   }

   ## surviving first level
   if (i==2) { probs[,1,]  <- 1-probs[,2,]; 
	    se.probs[,1,] <- se.probs[,2,]; 
	    lower[,1,] <- 1-lower[,2,]; 
	    upper[,1,] <- 1-upper[,2,]; 
   }

}# }}}

if (is.null(marks.call) & nstrata==1) {
dp <- -t(apply(cbind(probs[,-1,1],0),1,diff))
meanN <- apply(probs[,-1,1,drop=FALSE],1,sum)
meanN2 <- apply(t(exceed^2 * t(dp)),1,sum)
varN <- meanN2-meanN^2
} else dp <- meanN <- meanN2 <- varN <- NULL
 
colnames(probs) <- c(paste("N<",exceed[1],sep=""),paste("exceed>=",exceed,sep=""))
colnames(se.probs) <- c(paste("N<",exceed[1],sep=""),paste("exceed>=",exceed,sep=""))

out <- list(time=times,times=times,prob=probs,se.prob=se.probs,meanN=meanN,
	  lower=lower,upper=upper,meanN2=meanN2,varN=varN,exceed=exceed[-1],formula=form,
 cif.exceed=cif.exceed,dataList=dataList,
 nstrata=nstrata,strata.levels=strata.levels,strata.name=strata.name)
class(out) <- "exceed"
return(out)
}# }}}

##' @export
plot.exceed <- function(x,types=NULL,se=1,where=0.6,legend=NULL,strata=NULL,...) { ## {{{ 
if (is.null(types)) types <- 1:ncol(x$prob)
if (is.null(strata)) stratas <- seq(x$nstrata) else stratas <- strata
for (j in stratas) {
matplot(x$time,x$prob[,types,j],type="s",xlab="time",ylab="probabilty",xlim=range(c(0,x$times)),...)
if (x$nstrata>1) title(main=paste(x$strata.name,"=",x$strata.levels[j],sep=""))
if (se==1) 
for (i in types) plotConfRegion(x$time, cbind(x$lower[,i,j],x$upper[,i,j]) ,col=i)
if (!is.null(where)) {
if (is.null(legend)) legend <- colnames(x$prob)[types]
legend(0,0.6,legend=legend,lty=types,col=types)
}
}
} ## }}}

##' @export
summary.exceed <- function(object,times=NULL,types=NULL,strata=NULL,...) { ## {{{ 

if (is.null(types)) types <- 1:ncol(object$prob)
if (is.null(strata)) stratas <- seq(object$nstrata) else stratas <- strata
out <- list()

for (j in stratas) {
prob <- cbind(object$time,object$prob[,types,j])
se <- cbind(object$time,object$se.prob[,types,j])
lower <- cbind(object$time,object$lower[,types,j])
upper <- cbind(object$time,object$upper[,types,j])
if (is.null(times)) 
outl <- list(prob=prob,se=se,lower=lower,upper=upper)
else {
   rows <- fast.approx(object$time,times,type="left")
   probt <- cbind(times,prob[rows,])
   set <-   cbind(times,se[rows,])
   lowert <- cbind(times,lower[rows,])
   uppert <- cbind(times,upper[rows,])
   outl <- list(prob=probt,se=set,lower=lowert,upper=uppert)
}
if (length(stratas)>1) out[[j]] <- outl else out <- outl 
}
if (length(stratas)>1) names(out) <- object$strata.levels[stratas]
return(out)
} ## }}}

##' @export
print.exceed  <- function(x,...) summary.exceed(x,...)

