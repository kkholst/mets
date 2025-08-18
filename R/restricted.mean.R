##' Restricted IPCW mean for censored survival data 
##'
##' Simple and fast version for IPCW regression for just one time-point thus fitting the model 
##' \deqn{E( min(T, t) | X ) = exp( X^T beta) } or in the case of competing risks data
##' \deqn{E( I(epsilon=1) (t - min(T ,t)) | X ) = exp( X^T beta) } thus given years lost to 
##' cause, see \code{binreg} for the arguments. 
##'
##' When the status is binary assumes it is a survival setting and default is to consider outcome Y=min(T,t), 
##' if status has more than two levels, then computes years lost due to the specified cause, thus
##' using the response \deqn{ Y = (t-min(T,t)) I(status=cause) }
##'
##' Based on binomial regresion IPCW response estimating equation: 
##' \deqn{ X ( \Delta(min(T,t)) Y /G_c(min(T,t)) - exp( X^T beta)) = 0 }
##' for IPCW adjusted responses. Here \deqn{ \Delta(min(T,t)) = I ( min(T ,t) \leq C ) } is indicator of
##' being uncensored.  Concretely, the uncensored observations at time t will count those with an event (of any type) before t and those
##' with a censoring time at t or further out. One should therefore be a bit careful when data has been constructed such that
##' some of the event times T are equivalent to t. 
##'
##' Can also solve the binomial regresion IPCW response estimating equation: 
##' \deqn{ h(X) X ( \Delta(min(T,t)) Y /G_c(min(T,t)) - exp( X^T beta)) = 0 }
##' for IPCW adjusted responses where $h$ is given as an argument together with iid of censoring with h. 
##' 
##' By using appropriately  the h argument we can also do the efficient IPCW estimator estimator.
##' 
##' Variance is based on  \deqn{ \sum w_i^2 } also with IPCW adjustment, and naive.var is variance 
##' under known censoring model. 
##' 
##' When Ydirect is given it solves : 
##' \deqn{ X ( \Delta(min(T,t)) Ydirect /G_c(min(T,t)) - exp( X^T beta)) = 0 }
##' for IPCW adjusted responses. 
##'
##' The actual influence (type="II") function is based on augmenting with \deqn{ X \int_0^t E(Y | T>s) /G_c(s) dM_c(s) }
##' and alternatively just solved directly (type="I") without any additional terms. 
##'
##' Censoring model may depend on strata. 
##'
##' @param formula formula with outcome on Event form 
##' @param data data frame
##' @param outcome can do either rmst regression ('rmst') or years-lost regression  ('rmtl')
##' @param ... Additional arguments to binreg 
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001
##' # E( min(T;t) | X ) = exp( a+b X) with IPCW estimation 
##' out <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
##'                 time=50,cens.model=~strata(platelet),model="exp")
##' summary(out)
##' 
##' ## weighted GLM version   RMST
##' out2 <- logitIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
##'             time=50,cens.model=~strata(platelet),model="exp",outcome="rmst")
##' summary(out2)
##' 
##' ### time-lost
##' outtl <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
##'                 time=50,cens.model=~strata(platelet),model="exp",outcome="rmtl")
##' summary(outtl)
##' 
##' ### same as Kaplan-Meier for full censoring model 
##' bmt$int <- with(bmt,strata(tcell,platelet))
##' out <- resmeanIPCW(Event(time,cause!=0)~-1+int,bmt,time=30,
##'                              cens.model=~strata(platelet,tcell),model="lin")
##' estimate(out)
##' out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)
##' rm1 <- resmean.phreg(out1,times=30)
##' summary(rm1)
##' 
##' ### years lost regression
##' outl <- resmeanIPCW(Event(time,cause!=0)~-1+int,bmt,time=30,outcome="years-lost",
##'                              cens.model=~strata(platelet,tcell),model="lin")
##' estimate(outl)
##' 
##' ## competing risks years-lost for cause 1  
##' out <- resmeanIPCW(Event(time,cause)~-1+int,bmt,time=30,cause=1,
##'                             cens.model=~strata(platelet,tcell),model="lin")
##' estimate(out)
##' ## same as integrated cumulative incidence 
##' rmc1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=30)
##' summary(rmc1)
##' 
##' @export
##' @aliases rmstIPCW resmeanIPCWold 
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

##' Average Treatment effect for Restricted Mean for censored competing risks data using IPCW 
##'
##' Under the standard causal assumptions  we can estimate the average treatment effect E(Y(1) - Y(0)). We need Consistency, ignorability ( Y(1), Y(0) indep A given X), and positivity.
##'
##' The first covariate in the specification of the competing risks regression model must be the treatment effect that is a factor. If the factor has more than two levels
##' then it uses the mlogit for propensity score modelling.  We consider the outcome mint(T;tau) or  I(epsion==cause1)(t- min(T;t)) that gives years lost due to cause "cause" depending on 
##' the number of causes. The default model is the exp(X^ beta) and otherwise a linear model is used. 
##'
##' Estimates the ATE using the the standard binary double robust estimating equations that are IPCW censoring adjusted.
##'
##' @param formula formula with 'Event' outcome 
##' @param data data-frame 
##' @param model possible exp model for relevant mean model that is exp(X^t beta) 
##' @param outcome restricted mean time (rmst) or restricted mean time lost (rmtl)
##' @param ... Additional arguments to pass to binregATE 
##' @author Thomas Scheike
##' @examples
##' library(mets); data(bmt); bmt$event <- bmt$cause!=0; dfactor(bmt) <- tcell~tcell
##' out <- resmeanATE(Event(time,event)~tcell+platelet,data=bmt,time=40,treat.model=tcell~platelet)
##' summary(out)
##' 
##' out1 <- resmeanATE(Event(time,cause)~tcell+platelet,data=bmt,cause=1,time=40,
##'                    treat.model=tcell~platelet)
##' summary(out1)
##' 
##' ratioATE(out,out1,h=function(x) log(x))
##' @export
##' @aliases rmstATE ratioATE
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

##' @export
ratioATE <- function(rmtl,rmtl1,h=NULL) { ## {{{

if (is.null(h)) h <- function(x) x

coefG <- c(rmtl$riskG,rmtl1$riskG)
iidG <- cbind( rmtl$riskG.iid, rmtl1$riskG.iid)
covG <- crossprod(iidG)
###
coefDR <- c(rmtl$riskDR,rmtl1$riskDR)
iidDR <- cbind(rmtl$riskDR.iid, rmtl1$riskDR.iid)
covDR <- crossprod(iidDR)

ratioG <- estimate(coef=coefG,vcov=covG,f=function(p) h(c(p[3:4]/p[1:2],p[1]*p[4]/(p[2]*p[3]),p[2]*p[3]/(p[1]*p[4]))))
ratioDR <- estimate(coef=coefDR,vcov=covDR,f=function(p) h(c(p[3:4]/p[1:2],p[1]*p[4]/(p[2]*p[3]),p[2]*p[3]/(p[1]*p[4]))))
out <- list(ratioG=ratioG,ratioDR=ratioDR)
return(out)

} ## }}}

