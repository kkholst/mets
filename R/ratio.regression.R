##' Percentage of years lost due to cause regression 
##'
##' Estimates the percentage of the years lost that is due to a cause and how covariates affects this percentage by doing ICPW regression.
##'
##' Let the years lost be  \deqn{Y1= t- min(T ,) } and the years lost due to cause 1 \deqn{Y2= I(epsilon==1) ( t- min(T ,t) } , then
##' we model the ratio \deqn{logit( E(Y2 | X)/E(Y1 | X))  = X^T \beta }. Estimation is based on 
##' on binomial regresion IPCW response estimating equation: 
##' \deqn{ X ( \Delta^{ipcw}(t) Y2 expit(X^T \beta) -  Y1 ) = 0 }
##' where \deqn{\Delta^{ipcw}(t) = I((min(t,T)< C)/G_c(min(t,T)-)} is 
##' IPCW adjustment of the response \deqn{Y(t)= I(T \leq t, \epsilon=1 )}.  
##'
##' (type="I") sovlves this estimating equation using a stratified Kaplan-Meier for the
##' censoring distribution. For (type="II") the default an additional 
##' censoring augmentation term \deqn{X \int E(Y(t)| T>s)/G_c(s) d \hat M_c} is added.
##'
##' The variance is based on the squared influence functions that are also returned as the iid component. naive.var is variance 
##' under known censoring model. 
##'
##' Censoring model may depend on strata (cens.model=~strata(gX)). 
##'
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param cause cause of interest (numeric variable)
##' @param time  time of interest 
##' @param beta starting values 
##' @param type "II" adds augmentation term, and "I" classical outcome IPCW regression 
##' @param offset offsets for partial likelihood 
##' @param weights for score equations 
##' @param cens.weights censoring weights 
##' @param cens.model only stratified cox model without covariates
##' @param se to compute se's  based on IPCW 
##' @param kaplan.meier uses Kaplan-Meier for IPCW in contrast to exp(-Baseline)
##' @param cens.code gives censoring code
##' @param no.opt to not optimize 
##' @param method for optimization 
##' @param augmentation to augment binomial regression 
##' @param outcome  can do CIF regression "cif"=F(t|X), "rmtl"=E( t- min(T, t) | X)"
##' @param model  logit, exp or lin(ear) 
##' @param Ydirect use this Y instead of outcome constructed inside the program, should be a matrix with two column for numerator and denominator.
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data(bmt); bmt$time <- bmt$time+runif(408)*0.001
##' 
##' rmst30 <- rmstIPCW(Event(time,cause!=0)~platelet+tcell+age,bmt,time=30,cause=1)
##' rmst301 <- rmstIPCW(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=1)
##' rmst302 <- rmstIPCW(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=2)
##' 
##' estimate(rmst30)
##' estimate(rmst301)
##' estimate(rmst302)
##' 
##' ## percentage of total cumulative incidence due to cause 1
##' rmstratioI <- rmstRatio(Event(time,cause)~platelet+tcell+age,bmt,time=30,
##'                         cause=1,outcome="rmst")
##' summary(rmstratioI)
##' 
##' pp <- predict(rmstratioI,bmt)
##' ppb <- cbind(pp,bmt)
##' 
##' ## percentage of total cumulative incidence due to cause 1
##' cifratio <- binregRatio(Event(time,cause)~platelet+tcell+age,bmt,time=30,cause=1)
##' summary(cifratio)
##' pp <- predict(cifratio,bmt)
##' 
##' rmstratioI <- binregRatio(Event(time,cause)~platelet+tcell+age,bmt,
##'                                time=30,cause=1,outcome="rmst")
##' summary(rmstratioI)
##' 
##' pp <- predict(rmstratioI,bmt)
##' ppb <- cbind(pp,bmt)
##' @aliases rmstRatio
##' @export
binregRatio <- function(formula,data,cause=1,time=NULL,beta=NULL,type=c("II","I"),
	   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
	   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
	   outcome=c("cif","rmtl"),model=c("logit","exp","lin"),Ydirect=NULL,...)
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

 call.id <- id;
 conid <- construct_id(id,nrow(X))
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
 orig.id <- id

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- 1*(status %in% cens.code) 
  statusE <- (status %in% cause) & (exit<= time) 
  if (sum(statusE)==0) warning("No events of type 1 before time \n"); 
  kmt <- kaplan.meier

  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 
  data$id__ <- id
  data$exit <- exit
  data$statusC <- statusC 
  cens.strata <- cens.nstrata <- NULL 

 nevent <- sum((status %in% cause)*(exit<=time))
 ## if event before time or alive, then uncensored, equality for both censored and events  
 obs <- (exit<=time & (!statusC)) | (exit>=time)

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else { se <- FALSE; resC <- formC <- NULL}
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  p <- ncol(X)
  if (is.null(beta)) beta <- rep(0,p)
  if (is.null(augmentation))  augmentation=rep(0,p)
  X <-  as.matrix(X)
  X2  <- .Call("vecCPMat",X)$XX

 if (!is.null(Ydirect)) Y <-  Ydirect*obs/cens.weights else {
     if (outcome[1]=="cif") Y <- cbind( c((status %in% Causes)*(exit<=time)/cens.weights) , c((status %in% cause)*(exit<=time)/cens.weights) )
     else { 
            Y <- cbind(c((time-pmin(exit,time))*obs)/cens.weights, c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights)
     }
  }
  Yipcw <- Y

obj <- function(pp,all=FALSE)
{ # {{{
lp <- c(X %*% pp+offset)

if (model[1]=="logit") {
	p <- expit(lp)
        D2logl <- c(Y[,1]*weights*p/(1+exp(lp)))*X2
} else if (model[1]=="exp") {
	p <- exp(lp)
        D2logl <- c(Y[,1]*p*weights)*X2
} else {
	p <- lp
        D2logl <- c(Y[,1]*weights)*X2
}
Dlogl <- weights*X*c(Y[,1]*p-Y[,2])
D2log <- apply(D2logl,2,sum)
ploglik <- 0   ###sum(weights*(Y[,1]*p-Y[,2])^2)

gradient <- apply(Dlogl,2,sum)+augmentation
np <- length(pp)
hessian <- matrix(.Call("XXMatFULL",matrix(D2log,nrow=1),np,PACKAGE="mets")$XXf,np,np)

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      beta.iid <-  apply(beta.iid,2,sumstrata,id,max(id)+1)
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
	 id=id,Dlogl=Dlogl,iid=beta.iid,robvar=robvar,var=robvar,se.robust=diag(robvar)^.5)
      return(val)
  }  
 structure(0,gradient=-gradient/nid,hessian=-hessian/nid)
}# }}}

  ## first run without pseudo-value augmentation
  p <- ncol(X)
  opt <- NULL
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
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
  } else val <- c(list(coef=beta),obj(beta,all=TRUE))

 coefI <- val$coef
 gradientI <- val$gradient
 hessianI <- val$hessian
 ihessianI <- val$ihessian
 iidI <- val$iid

 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    X2 <-  X2[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    cens.weights <- cens.weights[ord]
    Y <- Y[ord,]
    lp <- c(X %*% val$coef+offset)
    p <- expit(lp)
    Yo <- Y[,1]*p-Y[,2]
    id <- id[ord]

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Yo,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    btime <- 1*(exit<time)
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply(h*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*h[xx$jumps+1,]/c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    mid <- max(id)
    MGCiidI <- MGCiid <- apply(MGt,2,sumstrata,xx$id,mid+1)

   if (type[1]=="II") { ##  pseudo-value type augmentation
    hYt  <-  revcumsumstrata(Yo,xx$strata,xx$nstrata)
    IhdLam0 <- cumsumstrata(hYt*S0i2*btime,xx$strata,xx$nstrata)
    U <- rep(0,length(xx$strata))
    U[xx$jumps+1] <- (resC$jumptimes<time)*hYt[xx$jumps+1]/c(resC$S0)
    MGt <- X*c(U-IhdLam0)*c(xx$weights)
    MGtiid <- apply(MGt,2,sumstrata,xx$id,mid+1)
    augmentation  <-  apply(MGtiid,2,sum) + augmentation
    ###
    EXt  <-  apply(X,2,revcumsumstrata,xx$strata,xx$nstrata)
    IEXhYtdLam0 <- apply(EXt*c(hYt)*S0i*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*hYt[xx$jumps+1]*EXt[xx$jumps+1,]/c(resC$S0)^2
    MGt2 <- (U[,drop=FALSE]-IEXhYtdLam0)*c(xx$weights)
    ###
    MGCiid2 <- apply(MGt2,2,sumstrata,xx$id,mid+1)
    ### Censoring Variance Adjustment 
    MGCiid <- MGCiid+(MGtiid-MGCiid2)
   }
   ## use data ordered by time (keeping track of id also)
   ## since X; Y and so forth are ordered in time
###   id <- xx$id
   }  else {
	 MGCiidI <-  MGCiid <- 0
  }## }}}

  ## first run without pseudo-value augmentation, then run with augmentation, if type="II"
 if (type[1]=="II") {
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
	  tim <- system.time(opt <- lava::NR(coefI,obj,...))
	  opt$timing <- tim
	  opt$estimate <- opt$par
      } else {
	  opt <- nlm(obj,beta,...)
	  opt$method <- "nlm"
      }
      cc <- opt$estimate; 
###	      if (!se) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
 }

 if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X)
  val <- c(val,list(time=time,formula=formula,formC=formC,
    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid))


  val$call <- cl
  val$MGciid <- MGCiid
  val$id <- orig.id
  val$call.id <- call.id
  val$name.id <- name.id
  val$nid <- nid
  val$iid.naive <- val$iid
  val$naive.var <- NULL 
  val$iidI <- iidI
  if (se) val$iid  <- val$iid+(MGCiid %*% val$ihessian)
  if (se) val$iidI  <- iidI+(MGCiidI %*% ihessianI)
  val$iid <- namesortme(val$iid,name.id)
  val$iidI <- namesortme(val$iidI,name.id)
  robvar <- crossprod(val$iid)
  val$var <-  val$robvar <- robvar
  val$se.robust <- diag(robvar)^.5
  val$se.coef <- diag(val$var)^.5
  val$cause <- cause
  val$cens.code <- cens.code 
  val$augmentation <- augmentation
  val$model <- model[1]
  val$outcome <- outcome[1]
  val$Yipcw <- Yipcw
  val$Causes <- Causes
  val$nevent <- nevent

  val$coefI <- coefI
  val$varI <- crossprod(val$iidI)
  val$se.coefI <- diag(val$varI)^.5
  val$gradientI <-  gradientI
  val$hessianI <- hessianI
  val$ihessianI <- ihessianI

  class(val) <- "binreg"
  return(val)
}# }}}

##' @export
rmstRatio <- function(formula,data,outcome=c("rmst","years-lost"),...)
{# {{{
   out <- binregRatio(formula,data,outcome=outcome[1],...)
   return(out)
}# }}}

