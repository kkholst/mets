##' Restricted IPCW mean for censored survival data 
##'
##' Simple and fast version for IPCW regression for just one time-point thus fitting the model 
##' \deqn{E( min(T, t) | X ) = exp( X^T beta) } or in the case of competing risks data
##' \deqn{E( I(epsilon=1) (t - min(T ,t)) | X ) = exp( X^T beta) } thus given years lost to 
##' cause. 
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
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param cause cause of interest
##' @param time  time of interest 
##' @param type of estimator
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
##' @param Ydirect to bypass the construction of the response Y=min(T,tau) and use this instead
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##'
##' data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001
##' # E( min(T;t) | X ) = exp( a+b X) with IPCW estimation 
##' out <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
##'                 time=50,cens.model=~strata(platelet),model="exp")
##' summary(out)
##' 
##'  ### same as Kaplan-Meier for full censoring model 
##' bmt$int <- with(bmt,strata(tcell,platelet))
##' out <- resmeanIPCW(Event(time,cause!=0)~-1+int,bmt,time=30,
##'                              cens.model=~strata(platelet,tcell),model="lin")
##' estimate(out)
##' out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)
##' rm1 <- resmean.phreg(out1,times=30)
##' summary(rm1)
##' 
##' ## competing risks years-lost for cause 1  
##' out <- resmeanIPCW(Event(time,cause)~-1+int,bmt,time=30,cause=1,
##'                             cens.model=~strata(platelet,tcell),model="lin")
##' estimate(out)
##' ## same as integrated cumulative incidence 
##' rmc1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=30)
##' summary(rmc1)
##' @export
##' @aliases rmstIPCW 
resmeanIPCW  <- function(formula,data,cause=1,time=NULL,type=c("II","I"),
   beta=NULL,offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",model="exp",
   augmentation=NULL,h=NULL,MCaugment=NULL,Ydirect=NULL,...)
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
  statusC <- (status %in% cens.code) 
  statusE <- (status %in% cause) & (exit<= time) 
  if ((sum(statusE)==0) & is.null(Ydirect)) warning("No events of type 1 before time \n"); 
  kmt <- kaplan.meier

  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 
  data$id <- id
  data$exit <- exit
  data$statusC <- statusC 

  cens.strata <- cens.nstrata <- NULL 

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  if (is.null(beta)) beta <- rep(0,ncol(X))
  p <- ncol(X)

  X <-  as.matrix(X)
  X2  <- .Call("vecMatMat",X,X)$vXZ
  ## if event before time or alive, then uncensored, equality for both censored and events  
  obs <- (exit<=time & (status %in% Causes)) | (exit>=time)
  if (is.null(Ydirect))  {
	  if (!competing) Y <- c(pmin(exit,time)*obs)/cens.weights else 
	                  Y <- c((status==cause)*(time-pmin(exit,time))*obs)/cens.weights
  } else Y <- c(Ydirect*obs)/cens.weights

 if (is.null(augmentation))  augmentation=rep(0,p)
 nevent <- sum((status %in% cause)*(exit<=time))

 h.call <- h
 if (is.null(h))  h <- rep(1,length(exit))

 if (!is.null(MCaugment)) {se <- FALSE;}

 if (se) {## {{{ censoring adjustment of variance 

    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    X2 <- X2[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    if (!is.null(Ydirect)) Ydirect <- Ydirect[ord]
    cens.weights <- cens.weights[ord]
    h <- h[ord]
###    lp <- c(X %*% val$coef+offset)
###    p <- exp(lp)
    obs <- (exit<=time & status %in% cause) | (exit>=time)
    if (is.null(Ydirect))  {
	  if (!competing) Y <- c(pmin(exit,time)*obs)/cens.weights else 
	                  Y <- c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights
    } else Y <- c(Ydirect*obs)/cens.weights
    if (model=="exp" & is.null(h.call))  ph <- 1
    if (model=="exp" & !is.null(h.call)) ph <- h
    if (model!="exp" & is.null(h.call))  ph <- 1 
    if (model!="exp" & !is.null(h.call)) ph <- h 
    Xd <- ph*X

	    xx <- resC$cox.prep
	    S0i2 <- S0i <- rep(0,length(xx$strata))
	    S0i[xx$jumps+1]  <- 1/resC$S0
	    S0i2[xx$jumps+1] <- 1/resC$S0^2
	    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i) 
	    ## to make \int h(s)/Ys  dM_i^C(s) 
	    btime <- 1*(exit<time)

	    ht  <-  apply(Xd*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
	    ### Cens-Martingale as a function of time and for all subjects to handle strata 
	    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
	    IhdLam0 <- apply(ht*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
	    U <- matrix(0,nrow(xx$X),ncol(X))
	    U[xx$jumps+1,] <- (resC$jumptimes<=time)*ht[xx$jumps+1,]/c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

	    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
	    MGCiid <- apply(MGt,2,sumstrata,xx$id,max(id)+1)

           if (type[1]=="II") { ## psedo-value type augmentation
	    hYt  <-  revcumsumstrata(Y,xx$strata,xx$nstrata)
	    IhdLam0 <- cumsumstrata(hYt*S0i2*btime,xx$strata,xx$nstrata)
	    U <- rep(0,length(xx$strata))
	    U[xx$jumps+1] <- (resC$jumptimes<=time)*hYt[xx$jumps+1]/c(resC$S0)
	    MGt <- Xd*c(U-IhdLam0)*c(xx$weights)
	    MGtiid <- apply(MGt,2,sumstrata,xx$id,max(id)+1)

	    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
	    MGCiid <- MGtiid
	    augmentation  <-  apply(MGCiid,2,sum) + augmentation
           }
   }  else {
	  MGCiid <- 0
  }## }}}

 ## use data ordered by time (keeping track of id also)
 id <- xx$id

obj <- function(pp,all=FALSE)
{ # {{{

lp <- c(X %*% pp+offset)
if (model=="exp") p <- exp(lp) else p <- lp
ploglik <- sum(weights*(Y-p)^2)

if (model=="exp")  {
if (is.null(h.call)) ph <- 1 else ph  <- h
Dlogl <- weights*ph*X*c(Y-p)
D2logl <- c(weights*ph*p)*X2
} else {
if (is.null(h.call)) ph <- 1 else ph  <- h
Dlogl <- weights*ph*X*c(Y-p)
D2logl <- c(weights*ph)*X2
}
D2log <- apply(D2logl,2,sum)
gradient <- apply(Dlogl,2,sum)+augmentation
hessian <- matrix(D2log,length(pp),length(pp))

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      beta.iid <- apply(beta.iid,2,sumstrata,id,max(id)+1)
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
    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid,Y=Y))
  
  ph <- 1
  lp <- c(X %*% val$coe+offset)
  if (model=="exp") p <- exp(lp) else p <- lp
  if (!is.null(h.call)) ph<- h 
  if (model=="exp" & is.null(h.call)) ph<- 1
  if (model!="exp" & is.null(h.call)) ph<- 1
  if (!is.null(MCaugment)) MGCiid <- MCaugment*ph*X 
  val$MGciid <- MGCiid
  val$MGtid <- id
  val$orig.id <- orig.id
  val$iid.origid <- ids 
  val$iid.naive <- val$iid 
  val$iid  <- val$iid+(MGCiid %*% val$ihessian)
  val$naive.var <- val$var
  robvar <- crossprod(val$iid)
  val$var <-  val$robvar <- robvar
  val$se.robust <- diag(robvar)^.5
  val$se.coef <- diag(val$var)^.5
  val$cens.code <- cens.code
  val$model.type <- model
  val$type <- type[1]
  val$augmentation <- augmentation

  class(val) <- c("binreg","resmean")
  return(val)
}# }}}

##' @export
rmstIPCW <- function(formula,data,...)
{# {{{
   out <- resmeanIPCW(formula,data,...)
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
##' then it uses the mlogit for propensity score modelling.  We consider the outcome mint(T;tau) or
##' I(epsion==cause1)(t- min(T;t)) that gives years lost due to cause "cause".  
##* The default model is the exp(X^ \beta) 
##'
##' Estimates the ATE using the the standard binary double robust estimating equations that are IPCW censoring adjusted.
##'
##' @param formula formula with 'Event' outcome 
##' @param data data-frame 
##' @param outcome  "rmst"=E( min(T, t) | X) , or "rmst-cause"=E( I(epsilon==cause) ( t - mint(T,t)) ) | X) 
##' @param model possible exp model for relevant mean model that is exp(X^t beta) 
##' @param ... Additional arguments to pass to binregATE 
##' @author Thomas Scheike
##' @examples
##' library(mets); data(bmt); bmt$event <- bmt$cause!=0; dfactor(bmt) <- tcell~tcell
##' out <- resmeanATE(Event(time,event)~tcell+platelet,data=bmt,time=40,treat.model=tcell~platelet)
##' summary(out)
##' 
##' out1 <- resmeanATE(Event(time,cause)~tcell+platelet,data=bmt,cause=1,outcome="rmst-cause",
##'                    time=40,treat.model=tcell~platelet)
##' summary(out1)
##' 
##' @export
##' @aliases rmstATE
resmeanATE <- function(formula,data,outcome=c("rmst","rmst-cause"),model="exp",...)
{# {{{
out <- 	binregATE(formula,data,...,outcome=outcome,model=model) 
return(out)
}# }}}

##' @export
rmstATE <- function(formula,data,outcome=c("rmst","rmst-cause"),model="exp",...)
{# {{{
out <- 	resmeanATE(formula,data,...,outcome=outcome,model=model) 
return(out)
}# }}}

