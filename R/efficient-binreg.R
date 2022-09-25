##' Efficient IPCW for binary data 
##'
##' Simple version of comp.risk function of timereg for just one time-point thus fitting the model 
##' \deqn{E(T \leq t | X ) = expit( X^T beta) }
##'
##' Based on binomial regresion IPCW response estimating equation: 
##' \deqn{ X  ( \Delta (T \leq t)/G_c(T_i-) - expit( X^T beta)) = 0 }
##' for IPCW adjusted responses. 
##'
##' Based on binomial regresion IPCW response estimating equation: 
##' \deqn{ h(X) X ( \Delta (T \leq t)/G_c(T_i-) - expit( X^T beta)) = 0 }
##' for IPCW adjusted responses where $h$ is given as an argument together with iid of censoring with h. 
##' By using appropriately the h argument we can also do the efficient IPCW estimator estimator this works 
##' the prepsurv and prepcif for survival or competing risks data. In this case also the censoring martingale 
##' should be given for variance calculation and this also comes out of the prepsurv or prepcif functions. 
##' (Experimental version at this stage).
##' 
##' Variance is based on  \deqn{ \sum w_i^2 } also with IPCW adjustment, and naive.var is variance 
##' under known censoring model. 
##'
##' Censoring model may depend on strata. 
##'
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param cause cause of interest
##' @param time  time of interest 
##' @param beta starting values 
##' @param offset offsets for partial likelihood 
##' @param weights for score equations 
##' @param cens.weights censoring weights 
##' @param cens.model only stratified cox model without covariates
##' @param se to compute se's  based on IPCW 
##' @param kaplan.meier uses Kaplan-Meier for IPCW in contrast to exp(-Baseline)
##' @param cens.code gives censoring code
##' @param no.opt to not optimize 
##' @param method for optimization 
##' @param model exp or linear 
##' @param augmentation to augment binomial regression 
##' @param h  h for estimating equation 
##' @param MCaugment iid of h and censoring model 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @export
Effbinreg <- function(formula,data,cause=1,time=NULL,beta=NULL,
   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",
   augmentation=NULL,h=NULL,MCaugment=NULL,...)
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
    stop("only right censored data, will not work for delayed entry\n"); 
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
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  ### possible handling of id to code from 0:(antid-1)
  if (!is.null(id)) {
          orig.id <- id
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else { orig.id <- NULL; nid <- nrow(X); id <- as.integer(seq_along(exit))-1; ids <- NULL}
  ### id from call coded as numeric 1 -> 
  id.orig <- id; 

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- (status==cens.code) 
  statusE <- (status==cause) & (exit<= time) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 
  kmt <- kaplan.meier

  statusC <- (status==cens.code) 
  data$id <- id
  data$exit <- exit
  data$statusC <- statusC 

  cens.strata <- cens.nstrata <- NULL 

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt)$surv
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  if (is.null(beta)) beta <- rep(0,ncol(X))
  p <- ncol(X)

  X <-  as.matrix(X)
  X2  <- .Call("vecMatMat",X,X)$vXZ
  Y <- c((status==cause)*(exit<=time)/cens.weights)
  nevent <- sum((status==cause)*(exit<=time))

  if (is.null(augmentation))  augmentation=rep(0,p)
  if (is.null(h))  h <- rep(1,length(exit))

obj <- function(pp,all=FALSE)
{ # {{{

lp <- c(X %*% pp+offset)
p <- expit(lp)
ploglik <- sum(weights*(Y-p)^2)

Dlogl <- weights*h*X*c(Y-p)
D2logl <- c(weights*h*p/(1+exp(lp)))*X2
D2log <- apply(D2logl,2,sum)
gradient <- apply(Dlogl,2,sum)+augmentation
hessian <- matrix(D2log,length(pp),length(pp))

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      beta.iid <-  apply(beta.iid,2,sumstrata,id,max(id)+1)
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
	 id=id,Dlogl=Dlogl,iid=beta.iid,robvar=robvar,var=robvar,se.robust=diag(robvar)^.5)
      return(val)
  }  
 structure(-ploglik,gradient=-gradient,hessian=hessian)
}# }}}

  p <- ncol(X)
  opt <- NULL
  if (p>0) {
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
	  tim <- system.time(opt <- lava::NR(beta,obj,...))
	  opt$timing <- tim
	  opt$estimate <- opt$par
      } else {
	  opt <- nlm(obj,beta,...)
	  opt$method <- "nlm"
      }
      cc <- opt$estimate; 
      if (!se) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }

  if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X)
  val <- c(val,list(time=time,formula=formula,formC=formC,
    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid))

 if (!is.null(MCaugment)) { se <- FALSE;}

 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    cens.weights <- cens.weights[ord]
    lp <- c(X %*% val$coef+offset)
    p <- expit(lp)
    Y <- c((status==cause)*weights*(exit<=time)/cens.weights)
    h <- h[ord]

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    ht  <-  apply(h*X*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply(ht*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*ht[xx$jumps+1,] /c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    MGCiid <- apply(MGt,2,sumstrata,xx$id,max(id)+1)
 
   }  else {
	  MGCiid <- 0
  }## }}}

  if (!is.null(MCaugment)) { MGCiid <- MCaugment*h*X}
  val$MGciid <- MGCiid
  val$MGtid <- id
  val$orig.id <- orig.id
  val$iid.origid <- ids 
  val$iid.naive <- val$iid 
  val$iid  <- beta.iid <- val$iid+(MGCiid %*% val$ihessian)
  beta.iid <- apply(beta.iid,2,sumstrata,id,max(id)+1)
  val$iid <- beta.iid
  val$naive.var <- val$var
  robvar <- crossprod(val$iid)
  val$var <-  val$robvar <- robvar
  val$se.robust <- diag(robvar)^.5
  val$se.coef <- diag(val$var)^.5

  class(val) <- "binreg"
  return(val)
}# }}}

prepsurv <- function(cs,ss,X,times,data) 
{# {{{
ttimes <- ss$cumhaz[,1]
ttimes <- ttimes[(ttimes< times)]
pcs <- predict(cs,data,se=FALSE,times=ttimes)
pcst0 <- c(tail(t(pcs$surv),n=1))
###
## survival model 
ps <- predict(ss,data,se=FALSE,times=ttimes)
dLamt <- apply(cbind(0,ps$cumhaz),1,diff)
pst0 <- c(tail(t(ps$surv),n=1))
Fsst0 <- 1-pst0
###
L <- apply((1/t(ps$surv*pcs$surv))*dLamt,2,sum)
L0 <- (1/pst0-1)

ctimes <- cs$cumhaz[,1]
ctimes <- ctimes[(ctimes< times)]
pcs <- predict(cs,data,se=FALSE,times=ctimes)
pcst0 <- c(tail(t(pcs$surv),n=1))
dLamc <- apply(cbind(0,pcs$cumhaz),1,diff)

###cctimes <- cs$cumhaz[,1]
###Lamcc <- cs$cumhaz[,2] %o%  exp(X[,2] * coef(cs)) 
###dLamcc <- diff(c(0,cs$cumhaz[cctimes<times,2])) %o%  exp(X[,2] * coef(cs)) 
###Lamcc[cctimes<times,]- t(pcs$cumhaz)
###matplot(ctimes, t(exp(-pcs$cumhaz)),type="s",lwd=0.5)

## survival model 
ps <- predict(ss,data,se=FALSE,times=ctimes)
pst0 <- c(tail(t(ps$surv),n=1))
Fsst0 <- 1-pst0
###L <- 1/(pst0*pcst0) - 1 - apply( (1/t(ps$surv*pcs$surv))*dLamc,2,sum)
h <- Fsst0/(pst0*L)

cmtimes <- matrix(ctimes,nrow(dLamc),ncol(dLamc))
tttimes <- matrix(pmin(data[,"time"],times),nrow(dLamc),ncol(dLamc),byrow=TRUE)
## augment
dF <- t(Fsst0-((1-ps$surv)))
Augmentf <- dF/t(ps$surv*pcs$surv)
AugmentC <- apply(Augmentf*dLamc*(cmtimes<tttimes),2,sum)
###
n <- nrow(data)
AdN <- rep(0,n)
jc <- (1:nrow(data))[cs$ord][cs$jumps]
jc <-jc[1:length(ctimes)] 
AdN[jc] <- mdi(Augmentf,1:length(ctimes),jc)
Augment <- AdN-AugmentC

saugment <- apply(X*Augment,2,sum)
augment <- apply(X*Augment*h,2,sum)
return(list(L=L,L0=L0,Mc=Augment,AugmentC=AugmentC,Xaugment=saugment,hXaugment=augment,h=h))
}# }}}

prepcif <- function(cs,ps,cif1,X,times,data) 
{# {{{

ctimes <- cs$cumhaz[,1]
ctimes <- ctimes[(ctimes< times)]
pcif1 <- predict(cif1,data,se=FALSE,times=ctimes)
pcift0 <- c(tail(t(pcif1$cif),n=1))
###
ctimes <- cs$cumhaz[,1]
ctimes <- ctimes[(ctimes< times)]
pcs <- predict(cs,data,se=FALSE,times=ctimes)
pcst0 <- c(tail(t(pcs$surv),n=1))
###
pps <- predict(ps,data,se=FALSE,times=ctimes)
dim(pps$surv)
pst0 <- c(tail(t(pps$surv),n=1))
psst0 <- 1-pst0
###
dLamc <- apply(cbind(0,pcs$cumhaz),1,diff)
Ia <- (1-pcift0)* ( pcift0/pcst0 - apply( t(pcif1$cif/pcs$surv)*dLamc,2,sum) )

times1 <- cif1$cumhaz[,1]
times1 <- times1[times1<times]
ppcs <- predict(cs,data,se=FALSE,times=times1)
ppcif1 <- predict(cif1,data,se=FALSE,times=times1)
dF1 <- apply(cbind(0,ppcif1$cif),1,diff)
I <- (1-pcift0)* apply(1/t(ppcs$surv)*dF1,2,sum)
###
ddF1 <-  t(((pcift0)-((pcif1$cif))))
II <- apply(ddF1^2/(t(pps$surv*pcs$surv))*dLamc,2,sum)-pcift0*apply(ddF1/(t(pcs$surv))*dLamc,2,sum) 
###
h <- pcift0*(1-pcift0)/(I-II) 
hc <- pcift0*(1-pcift0)/(Ia-II) 
###
cmtimes <- matrix(ctimes,nrow(dLamc),ncol(dLamc))
tttimes <- matrix(pmin(data[,"time"],times),nrow(dLamc),ncol(dLamc),byrow=TRUE)
med <- cmtimes<tttimes
## augment
dF <- t(((pcift0)-((pcif1$cif))))
Augmentf <- dF/t(pps$surv*pcs$surv)
AugmentC <- apply(Augmentf*dLamc*med,2,sum)
###
n <- nrow(data)
AdN <- rep(0,n)
dc <- length(ctimes)
jc <- (1:nrow(data))[cs$ord][cs$jumps]
jc <-jc[1:dc] 
AdN[jc] <- mdi(Augmentf,1:dc,jc)
Augment <- AdN- AugmentC
###
###X <- model.matrix(~formula,data)
saugment <- apply(X*Augment,2,sum)
augment <- apply(X*Augment*h,2,sum)

return(list(Mc=Augment,AugmentC=AugmentC,Xaugment=saugment,hXaugment=augment,h=h,hc=hc))
}# }}}

