##' Fast recurrent marginal mean when death is possible
##'
##' Fast Marginal means of recurrent events using the Lin and Ghosh (2000) 
##' standard errors.  
##' Fitting two models for death and recurent events these are
##' combined to prducte the estimator 
##' \deqn{ \int_0^t  S(u|x=0) dR(u|x=0) } the mean number of recurrent events, here
##' \deqn{ S(u|x=0) }  is the probability of survival for the baseline group, and 
##' \deqn{ dR(u|x=0) }  is the hazard rate of an event among survivors for the baseline. 
##' Here \deqn{ S(u|x=0) }  is estimated by \deqn{ exp(-\Lambda_d(u|x=0) }  with 
##'  \deqn{\Lambda_d(u|x=0) } being the cumulative baseline for death.
##' 
##' Assumes no ties in the sense that jump times needs to be unique, this is particularly so for the stratified version.
##'
##' @param formula with Event object
##' @param data data frame for computation
##' @param cause of interest (1 default)
##' @param ... Additional arguments to lower level funtions
##' @param death.code codes for death (terminating event, 2 default)
##' @author Thomas Scheike
##' 
##' @references 
##' Cook, R. J. and Lawless, J. F. (1997) Marginal analysis of recurrent events and a terminating event. Statist. Med., 16, 911–924.
##' Ghosh and Lin (2002) Nonparametric Analysis of Recurrent events and death, Biometrics, 554--562.
##' @examples
##' library(mets)
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##' hf$x <- as.numeric(hf$treatment) 
##'
##' ##  to fit non-parametric models with just a baseline 
##' xr <- phreg(Surv(entry,time,status==1)~cluster(id),data=hf)
##' dr <- phreg(Surv(entry,time,status==2)~cluster(id),data=hf)
##' par(mfrow=c(1,3))
##' plot(dr,se=TRUE)
##' title(main="death")
##' plot(xr,se=TRUE)
##' ### robust standard errors 
##' rxr <-  robust.phreg(xr,fixbeta=1)
##' plot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=4)
##' 
##' ## marginal mean of expected number of recurrent events 
##' ## out <- recurrentMarginalPhreg(xr,dr)
##' ## summary(out,times=1:5) 
##' 
##' ## marginal mean using formula  
##' outN <- recurrentMarginal(Event(entry,time,status)~cluster(id),hf,cause=1,death.code=2)
##' plot(outN,se=TRUE,col=2,add=TRUE)
##' summary(outN,times=1:5) 
##' 
##' ########################################################################
##' ###   with strata     ##################################################
##' ########################################################################
##' out <- recurrentMarginal(Event(entry,time,status)~strata(treatment)+cluster(id),
##'                          data=hf,cause=1,death.code=2)
##' plot(out,se=TRUE,ylab="marginal mean",col=1:2)
##' 
##' summary(out,times=1:5) 
##'
##' @aliases tie.breaker recmarg recurrentMarginalAIPCW  recurrentMarginalPhreg
##' @export
recurrentMarginal <- function(formula,data,cause=1,...,death.code=2)
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

  ### setting up formulae for the two phreg (cause of interest and death)
  if (is.null(call.id)) { 
     formid <- update.formula(formula,~.+cluster(id)) 
     data$id <- id
     tt <- terms(formid)
     tt <- delete.response(tt)
     formid <- formula(tt)
  }  else formid <- formula
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

  attr(meano,"cause") <- cause
  attr(meano,"death.code") <- death.code
  return(meano)
}# }}}

##' @export
recurrentMarginalPhreg <- function(recurrent,death,fixbeta=NULL,km=TRUE)
{# {{{
  xr <- recurrent
  dr <- death 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
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
 class(out) <- rep("recurrent",2)
 return(out)
}# }}}

##' @export
plot.recurrent <- function(x,ylab=NULL,...) {# {{{
 if (inherits(x,"recurrent") & is.null(ylab)) ylab <- "Mean events"
 baseplot(x,ylab=ylab,...)
}# }}}

##' @export
summary.recurrent <- function(object,...) {# {{{
out <- summaryRecurrentobject(object,...)
return(out)
}# }}}

##' @export
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

###summary.recurrent <- function(object,times=NULL,strata=NULL,estimates=FALSE,
###			      name="mean",robust=FALSE,
###	      conf.type=c("log","log-log","plain"),...) {# {{{
###base <- basecumhaz(object,joint=1,robust=robust)
###nstrata <- object$nstrata
###stratobs <- attr(base,"stratobs")
###
###baseci <- rep(list(NULL),nstrata)
###if (length(stratobs)>0) 
### for (i in stratobs) {
###   cumhaz <- base[[i+1]]$cumhaz
###   if (nrow(cumhaz)>=1) {
###     mu <- base[[i+1]]$cumhaz[,2]
###     se.mu <- base[[i+1]]$cumhaz[,3]
###	   if (length(mu)>=1) {
###	   conf <- conftype(mu,se.mu,conf.type=conf.type[1],...)
###	   out <- data.frame(times=cumhaz[,1],mu=mu,se.mu=se.mu,lower=conf$lower,upper=conf$upper,strata=i)
###	   names(out) <- c("times",name,"se","CI-2.5%","CI-97.5%","strata")
###	   baseci[[i+1]] <- out
###	   } 
###   } 
###   }
###
###pbaseci <- NULL
###if (!is.null(times)) {
### if (length(stratobs)>0) 
### for (i in stratobs) 
###if (!is.null(baseci[[i+1]])) pbaseci[[i+1]] <- predictCumhaz(rbind(0,baseci[[i+1]]),times) 
###}
###out <- list(baseci=baseci,pbaseci=pbaseci,times=times)
###
###class(out) <- "summary.recurrent"
###return(out)
###}# }}}

##' @export
print.summary.recurrent  <- function(x,...) {# {{{
if (is.null(x$times)) print(x$baseci) else print(x$pbaseci)
} # }}}

##' @export
recurrentMarginalAIPCW <- function(formula,data=data,cause=1,cens.code=0,death.code=2,
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
 data <- count.history(data,status="status__",id="id__",types=cause,multitype=TRUE)
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
  } else  St <- rep(1,nrow(xx$strata))
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
     gamma[gamma=Inf] <- 0
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
summaryTimeobject <-function(mutimes,mu,se.mu=NULL,times=NULL,type="log",...) {# {{{
 if (is.null(times)) times <- mutimes

 where <- fast.approx(c(0,mutimes),times,type="left")

 ##  see if object is vector or matrix
 if (is.matrix(mu)) mu <- rbind(0,mu)[where,] else mu <- c(0,mu)[where]
 if (!is.null(se.mu)) {
     if (is.matrix(se.mu)) se.mu <- rbind(0,se.mu)[where,] else se.mu <- c(0,se.mu)[where]
 se.logmu=se.mu/mu
 if (type=="log") {
 lower <- exp(log(mu) - 1.96*se.logmu)
 upper <- exp(log(mu) + 1.96*se.logmu)
 } else {
 lower <- mu - 1.96*se.mu
 upper <- mu + 1.96*se.mu
 }
 } else {se.mu <- se.logmu <- lower <- upper <- NA}


 out <- data.frame(times=times,mu=mu,se.mu=se.mu,lower=lower,upper=upper)
 names(out) <- c("times","mean","se-mean","CI-2.5%","CI-97.5%")
 return(out)
}# }}}

##' @export
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

##' @export
squareintHdM <- function(phreg,ft=NULL,fixbeta=NULL,beta.iid=NULL,...)
{# {{{
###  sum_k ( int_0^t f(s)/S_0^r(s) dM_k.^r(s) )^2
###  strata "r" from object and "k" id from cluster 
  if (!inherits(phreg,"phreg")) stop("Must be phreg object\n"); 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if ((phreg$no.opt) | is.null(phreg$coef)) fixbeta<- 1 else fixbeta <- 0

  x <- phreg
  xx <- x$cox.prep
  ww <- xx$caseweights*xx$weights
 
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/(x$S0^2*ww[xx$jump+1])
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
	     invhess <- -solve(x$hessian)
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

##' @export
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

##' @export
tie.breaker <- function(data,stop="time",start="entry",status="status",id=NULL,ddt=NULL,exit.unique=TRUE,cause=NULL,cens.code=0)
{# {{{

   if (!is.null(id)) id <- data[,id]
   ord <- 1:nrow(data)
   stat <- data[,status]
   time <- data[,stop]
   if (is.null(cause)) cause <- unique(stat)
   type0 <- which(cause %in% cens.code)
   if (length(type0)>0) cause <- cause[-type0]
   jumps <- stat %in% cause
   dupexit <- duplicated(time)
   time1 <- data[jumps,stop]
   time0 <- data[!jumps,stop]
   lt0 <- length(time0)
   ddp <- duplicated(c(time0,time1))
   if (exit.unique) ties <-ddp[(lt0+1):nrow(data)] else ties <- duplicated(c(time1))
   if (length(ties)>0) nties <- sum(ties) else nties <- 0
   if (nties>1) {
	   ordties <- ord[jumps][ties]
	   if (is.null(ddt)) {
		   abd <- abs(diff(data[,stop]))
		   abd <- min(abd[abd>0])
		   ddt <- abd*0.5
	   }
	   time[ordties] <- time[ordties]+runif(nties)*ddt

	   data[ordties,stop] <- time[ordties]
	   ties <- (ord %in% ordties)
	   if (!is.null(id)) {
	   lagties <- dlag(ties)
	   ### also move next start time if id the same 
	   change.start <- lagties==TRUE & id==dlag(id)
	   change.start[is.na(change.start)] <- FALSE
	   ocs <- ord[change.start]
	   data[ocs,start] <- data[ocs-1,stop]
	   data[,"tiebreaker"] <- FALSE
	   data[ocs,"tiebreaker"] <- TRUE
	   }
   }
   
   return(data)
 } # }}}

##' Simulation of recurrent events data based on cumulative hazards with two types of  recurrent events  
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Must give hazard of death and two recurrent events.  Possible with two
##' event types and their dependence can be specified but the two recurrent events need
##' to share random effect. Based on drawing the from cumhaz and cumhaz2 and 
##' taking the first event rather
##' the cumulative and then distributing it out. Key advantage of this is that 
##' there is  more flexibility wrt random effects 
##'
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param cumhaz2  cumulative hazard of recurrent events  of type 2
##' @param death.cumhaz cumulative hazard of death 
##' @param r1 potential relative risk adjustment of rate 
##' @param r2 potential relative risk adjustment of rate
##' @param rd potential relative risk adjustment of rate
##' @param rc potential relative risk adjustment of rate
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param dependence 0:independence; 1:all share same random effect with variance var.z; 2:random effect exp(normal) with correlation structure from cor.mat; 3:additive gamma distributed random effects, z1= (z11+ z12)/2 such that mean is 1 , z2= (z11^cor.mat(1,2)+ z13)/2, z3= (z12^(cor.mat(2,3)+z13^cor.mat(1,3))/2, with z11 z12 z13 are gamma with mean and variance 1 , first random effect is z1 and for N1 second random effect is z2 and for N2 third random effect is for death  
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param cens rate of censoring exponential distribution
##' @param ... Additional arguments to simRecurrentList
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##' library(mets)
##' data(CPH_HPN_CRBSI)
##' dr <- CPH_HPN_CRBSI$terminal
##' base1 <- CPH_HPN_CRBSI$crbsi 
##' base4 <- CPH_HPN_CRBSI$mechanical
##'
##'######################################################################
##' ### simulating simple model that mimicks data 
##'######################################################################
##' rr <- simRecurrent(5,base1)
##' dlist(rr,.~id,n=0)
##' rr <- simRecurrent(5,base1,death.cumhaz=dr)
##' dlist(rr,.~id,n=0)
##'
##' rr <- simRecurrent(100,base1,death.cumhaz=dr)
##' par(mfrow=c(1,3))
##' showfitsim(causes=1,rr,dr,base1,base1)
##' ######################################################################
##' ### simulating simple model 
##' ### random effect for all causes (Z shared for death and recurrent) 
##' ######################################################################
##' rr <- simRecurrent(100,base1,death.cumhaz=dr,dependence=1,var.z=0.4)
##' dtable(rr,~death+status)
##'
##' ######################################################################
##' ### simulating simple model that mimicks data 
##' ### now with two event types and second type has same rate as death rate
##' ######################################################################
##' set.seed(100)
##' rr <- simRecurrentII(100,base1,base4,death.cumhaz=dr)
##' dtable(rr,~death+status)
##' par(mfrow=c(2,2))
##' showfitsim(causes=2,rr,dr,base1,base4)
##'
##' ## general case, 3 recurrent events and 2 causes of death
##' set.seed(100)
##' cumhaz <- list(base1,base1,base4)
##' drl <- list(dr,base4)
##' rr <- simRecurrentList(100,cumhaz,death.cumhaz=drl,dependence=0)
##' dtable(rr,~death+status)
##' showfitsimList(rr,cumhaz,drl) 
##'
##' @export
##' @name simRecurrentII
##' @aliases showfitsim covIntH1dM1IntH2dM2 squareintHdM simRecurrentList showfitsimList
simRecurrentII <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,r1=NULL,r2=NULL,rd=NULL,rc=NULL,dependence=0,var.z=1,
			   cor.mat=NULL,cens=NULL,gap.time=FALSE,max.recurrent=100,...) 
{# {{{
cumhazL <- list(cumhaz,cumhaz2)
rr <- cbind(r1,r2)
if (!is.null(death.cumhaz)) death.cumhaz <- list(death.cumhaz)
if (!is.null(r1)) {
	if (!is.null(r2)) r2 <- rep(1,length(r1))
	rr <- cbind(r1,r2)
}
data <-     simRecurrentList(n,cumhazL,death.cumhaz=death.cumhaz,rr=rr,
		     rd=rd,rc=rc,dependence=dependence,var.z=var.z,
		     cor.mat=cor.mat,cens=cens,gap.time=gap.time,
		     max.recurrent=max.recurrent,...)
return(data)
}# }}}


##' @title Simulation of recurrent events data based on cumulative hazards for event and death process
##' @inheritParams simRecurrentII 
##' @inherit simRecurrentII examples author
##' @param ... Additional arguments to simRecurrentList
##' @export
simRecurrent <- function(n,cumhaz,death.cumhaz=NULL,r1=NULL,rd=NULL,rc=NULL,...) 
{# {{{
## wrapper for simRecurrentII without type-2 events
if (!is.null(death.cumhaz)) death.cumhaz <- list(death.cumhaz)
if (!is.null(r1)) r1 <- as.matrix(r1,ncol=1)
if (!is.null(rd)) rd <- as.matrix(rd,ncol=1)

rr <- simRecurrentList(n,list(cumhaz),death.cumhaz=death.cumhaz,rr=r1,rd=rd,rc=rc,...)
return(rr)
}# }}}

##' @export
simRecurrentList <- function(n,cumhaz,death.cumhaz=NULL,rr=NULL,rd=NULL,rc=NULL,zzr=NULL,zzd=NULL,
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

##' @export
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

##' @export
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

##' Simulation of two-stage recurrent events data based on Cox/Cox or Cox/Ghosh-Lin structure 
##'
##' Simulation of two-stage recurrent events data based on Cox/Cox or Cox/Ghosh-Lin structure 
##'
##' Must specify two phreg objects, or a phreg and a recreg object, then simulates data from two-stage model
##'
##' @param cox1 cox/ghosh-lin for recurrent events 
##' @param coxd cox for terminal event (phreg)
##' @param coxc possible cox for censoring (phreg)
##' @param n number of id's 
##' @param data on which the models are fitted (to draw covariates) 
##' @param type to specify type of simulation, if not default
##' @param id name of id variable
##' @param varz dependence frailty 
##' @param share to fit patly shared random effects model
##' @param cens censoring rate for exponential censoring
##' @param scale1 to scale baseline of recurrent events model
##' @param scaled to scale baseline of terminal event
##' @param dependence if dependence different from NULL, then uses simRecurrentList based on models given 
##' @param ... Additional arguments to simGLcox, nmin, nmax regulates linear approximation grid 
##' @author Thomas Scheike
##' @references 
##' Scheike (2024), Twostage recurrent events models, under review.
##' @examples
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##' hf$x <- as.numeric(hf$treatment)
##' n <- 100
##' xr <- phreg(Surv(entry,time,status==1)~x+cluster(id),data=hf)
##' dr <- phreg(Surv(entry,time,status==2)~x+cluster(id),data=hf)
##' simcoxcox <- sim.recurrent(xr,dr,n=n,data=hf)
##' recGL <- recreg(Event(entry,time,status)~x+cluster(id),hf,death.code=2)
##' simglcox <- sim.recurrent(recGL,dr,n=n,data=hf)
##'
#' @export sim.recurrent
#' @usage sim.recurrent(cox1,coxd=NULL,coxc=NULL,n=1,data=NULL,
#' type=c("default","cox-cox","gl-cox"),id="id",
#' varz=1,share=1,cens=0.001,scale1=1,scaled=1,dependence=NULL,...)
sim.recurrent <- function(cox1,coxd=NULL,coxc=NULL,
                          n=1,
                          data=NULL,
                          type=c("default","cox-cox","gl-cox"),
                          id="id",
                          varz=1,
                          share=1,
                          cens=0.001,
                          scale1=1,
                          scaled=1,
                          dependence=NULL,
                          ...) {# {{{
## exp censoring default
death <- NULL

if (type[1]=="default" & inherits(cox1,"recreg")) type <- "gl-cox" 
if (type[1]=="default" & inherits(cox1,"phreg")) type <- "cox-cox" 

scox1 <- draw.phreg(cox1,n,data=data)
if (!is.null(coxd)) scoxd <- draw.phreg(coxd,n,data=data,drawZ=FALSE,id=scox1$id)
if (!is.null(coxc)) scoxc <- draw.phreg(coxc,n,data=data,drawZ=FALSE,id=scox1$id)
if (type[1]=="cox-cox") type <- 3 else type <- 2
data <- scox1$data
ind <-  match(names(scox1$data), names(scoxd$data))
ind <- ind[!is.na(ind)]
if (length(ind)<ncol(scoxd$data))  data <- cbind(data,scoxd$data[,-ind])

Lam1 <- scalecumhaz(scox1$cumhaz,scale1); r1 <- scox1$rr
if (!is.null(coxc)) rc <-  scoxc$rr else rc <- rep(1,n)
if (!is.null(coxd))  {
LamD <- scalecumhaz(scoxd$cumhaz,scaled); rd <- scoxd$rr 
} else { LamD <- NULL; rd <- NULL; }
Lam2 <- scalecumhaz(scox1$cumhaz,0)
if (is.null(dependence) & (!is.null(LamD))) {
rrs <- simGLcox(n,Lam1,LamD,var.z=varz,r1=r1,rd=rd,rc=rc,model="twostage",cens=cens,type=type,share=share,...)
} else { 
if (is.null(dependence)) dependence <- 0
if (!is.null(LamD)) 
rrs <- simRecurrentList(n,list(Lam1),death.cumhaz=list(LamD),rr=matrix(r1,ncol=1),rd=matrix(rd,ncol=1),rc=rc,cens=cens,var.z=varz,dependence=dependence)
else rrs <- simRecurrentList(n,list(Lam1),rr=matrix(r1,ncol=1),rc=rc,cens=cens,var.z=varz,dependence=dependence)
rrs$Z <- attr(rrs,"z")[rrs$id]

rrs$statusD <- rrs$status
if (!is.null(LamD))  {
rrs <- dtransform(rrs,statusD=3,death==1)
}
rrs$id <- rrs$id-1
}

## add covariates 
rrs <- cbind(rrs,data[rrs$id+1,])

return(rrs)
}
# }}}

simRecurrentIIHist <- function(n,cumhaz,death.cumhaz,cens=NULL,rr=NULL,rc=NULL,rd=NULL,
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

##' Simulation of recurrent events data based on cumulative hazards: Two-stage model  
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Model is constructed such that marginals are on specified form by linear approximations
##' of cumulative hazards that are on a specific form to make them equivalent to marginals
##' after integrating out over survivors. Therefore E(dN_1 | D>t) = cumhaz, 
##' E(dN_2 | D>t) = cumhaz2,  and hazard of death is death.cumhazard 
##'
##' Must give hazard of death and two recurrent events.  Hazard of death is death.cumhazard  two
##' event types and their dependence can be specified but the two recurrent events need
##' to share random effect. 
##' 
##' Random effect for  death Z.death=(Zd1+Zd2), Z1=(Zd1^nu1) Z12,  Z2=(Zd2^nu2) Z12^nu3
##' \deqn{Z.death=Zd1+Zd2}  gamma distributions 
##' \deqn{Zdj}  gamma distribution  with mean parameters (sharej), vargamD,  share2=1-share1
##' \deqn{Z12}  gamma distribution with mean 1 and variance vargam12
##' 
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param cumhaz2  cumulative hazard of recurrent events  of type 2
##' @param death.cumhaz cumulative hazard of death 
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param nu powers of random effects where nu > -1/shape 
##' @param share1 how random effect for death splits into two parts 
##' @param vargamD variance of random effect  for death 
##' @param vargam12 shared random effect for N1 and N2 
##' @param cens rate of censoring exponential distribution
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##' data(CPH_HPN_CRBSI)
##' dr <- CPH_HPN_CRBSI$terminal
##' base1 <- CPH_HPN_CRBSI$crbsi 
##' base4 <- CPH_HPN_CRBSI$mechanical
##'
##' rr <- simRecurrentTS(1000,base1,base4,death.cumhaz=dr)
##' dtable(rr,~death+status)
##' showfitsim(causes=2,rr,dr,base1,base4)
##'
##' @export
simRecurrentTS <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,
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

##' Counts the number of previous events of two types for recurrent events processes
##'
##' Counts the number of previous events of two types for recurrent events processes
##'
##' @param data data-frame
##' @param status name of status 
##' @param id  id 
##' @param types types of the events (code) related to status (multiple values possible)
##' @param names.count name of Counts, for example Count1 Count2 when types=c(1,2)
##' @param lag if true counts previously observed, and if lag=FALSE counts up to know
##' @param multitype, if multitype is true then counts when status "in" types, otherwise counts for each value of type, types=c(1,2)
##' @param marks values related to status ("in" types), counts marks for types, only when multitype=TRUE
##' @author Thomas Scheike
##' @examples
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##' dtable(hf,~status)
##' rr <-  count.history(hf,types=1:2,id="id",status="status")
##' dtable(rr,~"Count*"+status,level=1)
##'
##' @aliases count.historyVar 
##' @export
count.history <- function(data,status="status",id="id",types=1,names.count="Count",lag=TRUE,multitype=FALSE,marks=NULL)
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

##' @export
count.historyVar <- function(data,var="status",id="id",names.count="Count",lag=TRUE)
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

##' Estimation of probability of more that k events for recurrent events process
##'
##' Estimation of probability of more that k events for recurrent events process
##' where there is terminal event, based on this also estimate of variance of recurrent events. The estimator is based on cumulative incidence of exceeding "k" events.
##' In contrast the probability of exceeding k events can also be computed as a counting process integral. 
##'
##' @param formula formula
##' @param data  data-frame 
##' @param cause of interest 
##' @param death.code for status 
##' @param cens.code censoring codes
##' @param exceed values (if not given then all observed values)
##' @param marks may be give for jump-times and then exceed values needs to be specified
##' @param all.cifs if true then returns list of all fitted objects in cif.exceed 
##' @param return.data if true then returns list of data for fitting the different excess thresholds 
##' @param conf.type  type of confidence interval c("log","plain")
##' @param  level of confidence intervals default is 0.95
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @references 
##'    Scheike, Eriksson, Tribler (2019), The mean, variance and correlation for bivariate 
##'                                        recurrent events with a terminal event, JRSS-C
##' @examples
##' library(mets)
##' data(hfactioncpx12)
##' dtable(hfactioncpx12,~status)
##' 
##' oo <- prob.exceed.recurrent(Event(entry,time,status)~cluster(id),
##'         hfactioncpx12,cause=1,death.code=2)
##' plot(oo)
##' summary(oo,times=c(1,2,5))
##' 
##' @export
##' @aliases summaryRecurrentobject summaryTimeobject
##' @export
prob.exceed.recurrent <- function(formula,data,cause=1,death.code=2,cens.code=0, exceed=NULL,marks=NULL,all.cifs=FALSE,
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
      strata.name  <-  names(des$levels)[1]
    } else strata.name <- NULL
    id      <- des$cluster
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
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

###prob.exceed.recurrent <- function(formula,data, cause=1, death.code=2, cens.code=0,
###                      exceed=NULL,marks=NULL,cifmets=TRUE,all.cifs=FALSE,return.data=FALSE,...)
###{# {{{
###    cl <- match.call()# {{{
###    m <- match.call(expand.dots = TRUE)[1:3]
###    des <- proc_design(
###        formula,
###        data = data,
###        specials = c("marks","strata","offset", "weights", "cluster"),
###        intercept = TRUE
###    )
###    Y <- des$y
###    if (!inherits(Y, c("Event", "Surv"))) {
###        stop("Expected a 'Surv' or 'Event'-object")
###    }
###    if (ncol(Y) == 2) {
###        exit <- Y[, 1]
###        entry <- rep(0, nrow(Y))
###        status <- Y[, 2]
###    } else {
###        entry <- Y[, 1]
###        exit <- Y[, 2]
###        status <- Y[, 3]
###    }
###    X <- des$x
###    des.weights <- des$weights
###    des.offset  <- des$offset
###    id      <- des$cluster
###    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
###    strata <- des$strata
###    if (!is.null(strata))  {
###      ns <- grep("strata",names(des$levels))
###      strata.name  <-  names(des$levels)[1]
###    } else strata.name <- NULL
###    des.marks <- des$marks
###    id      <- des$cluster
###    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
###    ## no use of 
###    pos.cluster <- pos.strata <- NULL
###
### ## take offset and weight first from formula, but then from arguments
###  if (is.null(des.offset)) {
###	  if (is.null(offset)) offset <- rep(0,length(exit)) 
###  } else offset <- des.offset
###  if (is.null(des.weights)) {
###	  if (is.null(weights)) weights <- rep(1,length(exit)) 
###  } else weights <- des.weights
###  if (!is.null(des.marks) & is.null(marks))  marks <- des.marks
###  marks.call <- marks
###  if (is.null(marks)) marks <- rep(1,length(statusE))
###  ## }}}
###
### call.id <- id;
### conid <- construct_id(id,nrow(X),as.data=TRUE)
### name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
###
### times <- NULL
### statusD <- (status %in% death.code)*1
### statusE <- (status %in% cause)*1
### if (sum(statusE)==0) stop("none of type events")
### if (!is.null(strata) & !cifmets) stop("strata only for cifmets=TRUE\n")
### allvars <- all.vars(formula)
### if (is.null(times)) times <- sort(unique(exit[statusE==1]))
###
### countE <- cumsumstrata(statusE*marks,id,nid)
### if (is.null(marks.call)) {
###    mc <- max(countE)+1
###    idcount <- id*mc+countE
###    idcount <- cumsumstrata(rep(1,length(idcount)),idcount,mc*(nid+1))
### } else {
###    ex <- outer(countE,exceed,">=")[,1,]
###    exC <- apply(ex,2,cumsumstrata,id,nid)
###    idcount <- exC
### }
###
### if (is.null(exceed)) exceed <- sort(unique(countE))
### w0 <- which(exceed==0)
### if (length(w0)>=1) exceed <- exceed[-w0]
###
###if (!cifmets) {
###  ## can not handle start stop and id, so uses MG standard errors 
###  pp <- as.formula(paste("Hist(entry=",allvars[1],",",allvars[2],",statN)~+1",sep=""))
###  ###rhs <- update(formula,-1~.)
###  ###form <- as.formula(update.formula(rhs,pp))
###  form <- as.formula(pp)
###} else {
###  pp <- as.formula(paste("Event(",allvars[1],",",allvars[2],",statN)~.",sep=""))
###  rhs <- update(formula,-1~.)
###  form <- as.formula(update.formula(rhs,pp))
###}
###
###cif.exceed <- NULL
###dataList <- NULL
###if (all.cifs) cif.exceed <- list() 
###if (return.data) dataList <- list() 
###se.probs <- probs <- matrix(0,length(times),length(exceed)+1)
###lower <-  matrix(0,length(times),length(exceed)+1)
###upper <-  matrix(0,length(times),length(exceed)+1)
###i <- 1
###for (n1 in exceed) {# {{{
###	i <- i+1
###	### first time that get to n1
###        if (is.null(marks.call)) keep <- (countE<n1 ) | (countE==n1 & idcount==1) else keep <- exC[,i-1]<=1 
######	keep <- (countE<n1 ) | (firsttoE==n1)
###
###	### status, censoring, get to n1, or die
###        statN <- rep(0,nrow(data))
###        if (is.null(marks.call)) statN[countE==n1] <- 1 else statN[exC[,i-1]==1] <- 1
###	statN[statusD==1] <- 2
###	statN <- statN[keep]
###	dataN <- data[keep,]
###	dataN$statN <- statN
###	if (return.data) dataList[[i-1]] <- dataN
###
###        if (!cifmets) pN1 <-  suppressWarnings(prodlim::prodlim(form,data=dataN)) else pN1 <-  suppressWarnings(cif(form,data=dataN))
###	if (all.cifs) cif.exceed[[i-1]] <- pN1
###
###	if ((sum(statN==1)==0) | !cifmets) {
###	lower[,i] <- upper[,i] <- se.probs[,i] <- probs[,i] <- rep(0,length(times)) } else  {
###
###        where <- fast.approx(c(0,pN1$times),times,type="left")
######     cifs <-cbind(pN1$times,vecAllStrata(pN1$cumhaz[,2],pN1$strata,pN1$nstrata))
######     se.cifs <-cbind(pN1$times,vecAllStrata(pN1$se.cumhaz[,2],pN1$strata,pN1$nstrata))
###	probs[,i] <- c(0,pN1$mu)[where]
###	se.probs[,i] <- c(0,pN1$se.mu)[where]
###	lower[,i] <- c(0,exp(log(pN1$mu)-1.96*pN1$se.mu/pN1$mu))[where]
###	upper[,i] <- c(0,exp(log(pN1$mu)+1.96*pN1$se.mu/pN1$mu))[where]
###	}
###	if (i==2) { probs[,1]    <- 1-probs[,2]; 
###                    se.probs[,1] <- se.probs[,2]; 
###                    lower[,1] <- 1-lower[,2]; 
###                    upper[,1] <- 1-upper[,2]; 
###	}
###}# }}}
###
###if (is.null(marks.call)) {
###dp <- -t(apply(cbind(probs[,-1],0),1,diff))
###meanN <- apply(probs[,-1,drop=FALSE],1,sum)
###meanN2 <- apply(t(exceed^2 * t(dp)),1,sum)
###varN <- meanN2-meanN^2
###} else dp <- meanN <- meanN2 <- varN <- NULL
### 
###colnames(probs) <- c(paste("N<",exceed[1],sep=""),paste("exceed>=",exceed,sep=""))
###colnames(se.probs) <- c(paste("N<",exceed[1],sep=""),paste("exceed>=",exceed,sep=""))
###
###out <- list(time=times,times=times,prob=probs,se.prob=se.probs,meanN=meanN,
###	    lower=lower,upper=upper,meanN2=meanN2,varN=varN,exceed=exceed[-1],formula=form,
###	    cif.exceed=cif.exceed,dataList=dataList)
###class(out) <- "exceed"
###return(out)
###}# }}}
###
###prob.exceed.recurrentO <- function(data,type,status="status",death="death",
### start="start",stop="stop",id="id",times=NULL,exceed=NULL,cifmets=TRUE,strata=NULL,all.cifs=FALSE,...)
###{# {{{
###### setting up data 
###stat <-     data[,status]
###dd   <-     data[,death]
###tstop <-    data[,stop]
###tstart <-   data[,start]
###clusters <- data[,id]
###
###if (sum(stat==type)==0) stop("none of type events")
###if (!is.null(strata) & !cifmets) stop("strata only for cifmets=TRUE\n")
###
###if (is.numeric(clusters)) {
###      clusters <- fast.approx(unique(clusters), clusters) - 1
###      max.clust <- max(clusters)
###}
###else {
###     max.clust <- length(unique(clusters))
###     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
###}
###
### count <- cumsumstrata((stat==type),clusters,max.clust+1)
###### count  <- cumsumidstratasum((stat==type),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
### mc <- max(count)+1
### idcount <- clusters*mc + count
### idcount <- cumsumstrata(rep(1,length(idcount)),idcount,mc*(max.clust+1))
###
###if (is.null(times)) times <- sort(unique(tstop[stat==type]))
###if (is.null(exceed)) exceed <- sort(unique(count))
###
###if (!cifmets) {
###   if (is.null(strata)) form <- as.formula(paste("Hist(entry=",start,",",stop,",statN)~+1",sep=""))
###   else form <- as.formula(paste("Hist(entry=",start,",",stop,",statN)~+",strata,sep="")) 
###}
###else {
###   if (is.null(strata)) form <- as.formula(paste("Event(",start,",",stop,",statN)~+1",sep=""))
###   else form <- as.formula(paste("Event(",start,",",stop,",statN)~strata(",strata,")",sep=""))
###}
###
###cif.exceed <- NULL
###if (all.cifs) cif.exceed <- list() 
###probs.orig <- se.probs <- probs <- matrix(0,length(times),length(exceed))
###se.lower <-  matrix(0,length(times),length(exceed))
###se.upper <-  matrix(0,length(times),length(exceed))
###i <- 1
###for (n1 in exceed[-1]) {# {{{
###	i <- i+1
###	### first time that get to n1
###	keep <- (count<n1 ) | (count==n1 & idcount==1)
###	### status, censoring, get to n1, or die
###        statN <- rep(0,nrow(data))
###	statN[count==n1] <- 1
###	statN[dd==1] <- 2
###	statN <- statN[keep]
###        if (!cifmets) 
###	pN1 <-  suppressWarnings(prodlim::prodlim(form,data=data[keep,]))
###        else pN1 <-  suppressWarnings(cif(form,data=data[keep,]))
###	if (all.cifs) cif.exceed[[i-1]] <- pN1
###
###	if (sum(statN)==0) {
###		se.lower[,i] <- se.upper[,i] <- se.probs[,i] <- probs[,i] <- rep(0,length(times)) } else  {
###
###                if (!cifmets) {
###			mps  <- summary(pN1,times=times,cause=1)
###			mps  <- suppressWarnings(summary(pN1,times=times,cause=1)$table)
###			if (is.list(mps)) mps <- mps$"1"
###			probs.orig[,i] <- ps <- mps[,5]
###			mm <- which.max(ps)
###			probs[,i] <- ps
###			probs[is.na(ps),i] <- ps[mm]
###			se.probs[,i] <- mps[,6]
###			se.probs[is.na(ps),i] <- se.probs[mm,i]
###			se.lower[,i] <- mps[,7] 
###			se.lower[is.na(ps),i] <- se.lower[mm,i]
###			se.upper[,i] <- mps[,8]
###			se.upper[is.na(ps),i] <- se.upper[mm,i]
###	        } else {
###			where <- fast.approx(c(0,pN1$times),times,type="left")
###	   	        probs[,i] <- c(0,pN1$mu)[where]
###			se.probs[,i] <- c(0,pN1$se.mu)[where]
###			se.lower[,i] <- probs[,i]-1.96*se.probs[,i] 
###			se.upper[,i] <- probs[,i]+1.96*se.probs[,i] 
###		}
###
###	}
###	if (i==2) { probs[,1]    <- 1-probs[,2]; 
###                    se.probs[,1] <- se.probs[,2]; 
###                    se.lower[,1] <- 1-se.lower[,2]; 
###                    se.upper[,1] <- 1-se.upper[,2]; 
###	}
###}# }}}
###
###dp <- -t(apply(cbind(probs[,-1],0),1,diff))
###meanN <- apply(probs[,-1,drop=FALSE],1,sum)
###meanN2 <- apply(t(exceed[-1]^2 * t(dp)),1,sum)
### 
###colnames(probs) <- c(paste("N=",exceed[1],sep=""),paste("exceed>=",exceed[-1],sep=""))
###colnames(se.probs) <- c(paste("N=",exceed[1],sep=""),paste("exceed>=",exceed[-1],sep=""))
###
###return(list(time=times,times=times,prob=probs,se.prob=se.probs,meanN=meanN,probs.orig=probs.orig[,-1],
###	    se.lower=se.lower,se.upper=se.upper,meanN2=meanN2,varN=meanN2-meanN^2,exceed=exceed[-1],formula=form,
###	    cif.exceed=cif.exceed))
###}# }}}
###
#####' @export
###plot.exceed <- function(x,types=NULL,se=1,where=0.6,legend=NULL,...) { ## {{{ 
###if (is.null(types)) types <- 1:ncol(x$prob)
###matplot(x$time,x$prob[,types],type="s",xlab="time",ylab="probabilty",xlim=range(c(0,x$times)),...)
###if (se==1) 
###for (i in types) plotConfRegion(x$time,cbind(x$lower[,i],x$upper[,i]),col=i)
###if (!is.null(where)) {
###if (is.null(legend)) legend <- colnames(x$prob)[types]
###legend(0,0.6,legend=legend,lty=types,col=types)
###}
###} ## }}}
###
#####' @export
###summary.exceed <- function(object,times=NULL,types=NULL,...) { ## {{{ 
###if (is.null(types)) types <- 1:ncol(object$prob)
###prob <- cbind(object$time,object$prob[,types])
###se <- cbind(object$time,object$se.prob[,types])
###lower <- cbind(object$time,object$lower[,types])
###upper <- cbind(object$time,object$upper[,types])
###if (is.null(times)) 
###out <- list(prob=prob,se=se,lower=lower,upper=upper)
###else {
###   rows <- fast.approx(object$time,times)
###   out <- list(prob=prob[rows,],se=se[rows,],lower=lower[rows,],upper=upper[rows,])
###}
###return(out)
###} ## }}}
###

