##' Binomial Regression for censored competing risks data 
##'
##' Simple version of comp.risk function of timereg for just one time-point thus fitting the model 
##' \deqn{P(T \leq t, \epsilon=1 | X ) = expit( X^T beta) }
##'
##' Based on binomial regresion IPCW response estimating equation: 
##' \deqn{ X ( \Delta^{ipcw}(t) I(T \leq t, \epsilon=1 ) - expit( X^T beta)) = 0 }
##' where \deqn{\Delta^{ipcw}(t) = I((min(t,T)< C)/G_c(min(t,T)-)} is 
##' IPCW adjustment of the response \deqn{Y(t)= I(T \leq t, \epsilon=1 )}.  
##'
##' (type="I") sovlves this estimating equation using a stratified Kaplan-Meier for the
##' censoring distribution. For (type="II") the default an additional 
##' censoring augmentation term \deqn{X \int E(Y(t)| T>s)/G_c(s) d \hat M_c} is added.
##'
##' logitIPCW instead considers 
##' \deqn{ X  I(min(T_i,t) < G_i)/G_c(min(T_i ,t)) ( I(T \leq t, \epsilon=1 ) - expit( X^T beta)) = 0 }
##' a standard logistic regression with weights that adjust for IPCW.
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
##' @param type "II" adds augmentation term, and "I" classic binomial regression 
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
##' @param outcome  can do CIF regression "cif"=F(t|X), "rmst"=E( min(T, t) | X) , or years-lost "rmtl"=E( I(epsilon==cause) ( t - mint(T,t)) ) | X) 
##' @param model link functions used, with defaults logit for cif, exp for rmst or rmtl, but can be logit, exp or lin (for identity link)
##' @param Ydirect use this Y instead of outcome constructed inside the program (e.g. I(T< t, epsilon=1)), then uses IPCW vesion of the Y, set outcome to "rmst" to fit using the model specified by model
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data(bmt); bmt$time <- bmt$time+runif(408)*0.001
##' # logistic regresion with IPCW binomial regression 
##' out <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50)
##' summary(out)
##' head(iid(out))
##'
##' predict(out,data.frame(tcell=c(0,1),platelet=c(1,1)),se=TRUE)
##'
##' outs <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50,cens.model=~strata(tcell,platelet))
##' summary(outs)
##'
##' ## glm with IPCW weights 
##' outl <- logitIPCW(Event(time,cause)~tcell+platelet,bmt,time=50)
##' summary(outl)
##'
##' ##########################################
##' ### risk-ratio of different causes #######
##' ##########################################
##' data(bmt)
##' bmt$id <- 1:nrow(bmt)
##' bmt$status <- bmt$cause
##' bmt$strata <- 1
##' bmtdob <- bmt
##' bmtdob$strata <-2
##' bmtdob <- dtransform(bmtdob,status=1,cause==2)
##' bmtdob <- dtransform(bmtdob,status=2,cause==1)
##' ###
##' bmtdob <- rbind(bmt,bmtdob)
##' dtable(bmtdob,cause+status~strata)
##' 
##' cif1 <- cif(Event(time,cause)~+1,bmt,cause=1)
##' cif2 <- cif(Event(time,cause)~+1,bmt,cause=2)
##' plot(cif1)
##' plot(cif2,add=TRUE,col=2)
##' 
##' cifs1 <- binreg(Event(time,cause)~tcell+platelet+age,bmt,cause=2,time=50)
##' cifs2 <- binreg(Event(time,cause)~tcell+platelet+age,bmt,cause=2,time=50)
##' summary(cifs1)
##' summary(cifs2)
##' 
##' cifdob <- binreg(Event(time,status)~-1+factor(strata)+
##' 	 tcell*factor(strata)+platelet*factor(strata)+age*factor(strata)
##' 	 +cluster(id),bmtdob,cause=1,time=50,cens.model=~strata(strata))
##' summary(cifdob)
##' 
##' riskratio <- function(p) {
##'   Z <- rbind(c(1,0,1,1,0,0,0,0), c(0,1,1,1,0,1,1,0))
##'   lp <- c(Z %*% p)
##'   p <- lava::expit(lp)
##'   return(p[1]/p[2])
##' }
##' 
##' lava::estimate(cifdob,f=riskratio)
##'
##' @aliases logitIPCW binregt
##' @export
binreg <- function(formula,data,cause=1,time=NULL,beta=NULL,type=c("II","I"),
	   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
	   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
	   outcome=c("cif","rmst","rmtl"),model=c("default","logit","exp","lin"),
	   Ydirect=NULL,...)
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
  call.id <- id 
  conid <- construct_id(id,nrow(X),namesX=rownames(X))
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
  ## id before time-sorting later 
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

  if (is.null(cens.weights))  { ## {{{
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else { se <- FALSE; resC <- formC <- NULL} ## }}}
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  p <- ncol(X)
  if (is.null(beta)) beta <- rep(0,p)
  if (is.null(augmentation))  augmentation=rep(0,p)
  X <-  as.matrix(X)
  X2  <- .Call("vecCPMat",X)$XX

 if (!is.null(Ydirect)) Y <-  Ydirect*obs/cens.weights else {
     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time)/cens.weights)
     else { if (!competing) {
	     if (outcome[1]=="rmst")
	     Y <-  c(pmin(exit,time)*obs)/cens.weights 
             else Y <-  c((time-pmin(exit,time))*obs)/cens.weights 
            } else Y <- c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights
     }
  }
  Yipcw <- Y

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
###    lp <- c(X %*% val$coef+offset)
###    p <- expit(lp)
###    Y <- c((status %in% cause)*weights*(exit<=time)/cens.weights)
    Y <- Y[ord]

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    btime <- 1*(exit<time)
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply(h*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*h[xx$jumps+1,]/c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    mid <- max(id)
    MGCiid <- apply(MGt,2,sumstrata,xx$id,mid+1)

   if (type[1]=="II") { ##  pseudo-value type augmentation
    hYt  <-  revcumsumstrata(Y,xx$strata,xx$nstrata)
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
   id <- xx$id
   }  else {
	  MGCiid <- 0
  }## }}}

 ## default links
 if (model[1]=="default") {
	 if (outcome[1]=="cif") model <- "logit"
	 if (outcome[1]=="rmst") model <- "exp"
	 if (outcome[1]=="rmtl") model <- "exp"
	 if (outcome[1]=="years-lost") model <- "exp"
 }

obj <- function(pp,all=FALSE)
{ # {{{
lp <- c(X %*% pp+offset)

    if (model[1]=="exp") {
	 p <- exp(lp) 
         D2logl <- c(weights*p)*X2 
     } else if (model[1]=="lin") {
	 p <- lp
         D2logl <- c(weights)*X2
       }
     else if (model[1]=="logit") {
	p <- expit(lp)
        D2logl <- c(weights*p/(1+exp(lp)))*X2
} else stop("link functions must be logit,exp,lin\n") 
ploglik <- sum(weights*(Y-p)^2)
## to avoid using ploglik for fitting, with exp only solve estimating equation
if (model[1]=="exp") ploglik <- 0

Dlogl <- weights*X*c(Y-p)
D2log <- apply(D2logl,2,sum)
gradient <- apply(Dlogl,2,sum)+augmentation
np <- length(pp)
hessian <- matrix(.Call("XXMatFULL",matrix(D2log,nrow=1),np,PACKAGE="mets")$XXf,np,np)

  if (all) {
      ploglik <- sum(weights*(Y-p)^2)
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

## setting default for NR 
dots <- list(...)
if (length(dots)==0) {
   if (model[1]=="exp") control <- list(tol=1e-10,stepsize=0.5)  
   else control <- NULL
} else control <- dots[[1]]

  p <- ncol(X)
  opt <- NULL
  if (p>0) {
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
          tim <- system.time(opt <- lava::NR(beta,obj,control=control))
	  opt$timing <- tim
	  opt$estimate <- opt$par
      } else {
	  opt <- nlm(obj,beta,...)
	  opt$method <- "nlm"
      }
      cc <- opt$estimate; 
###	      if (!se) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else {
	      val <- c(list(coef=beta),obj(beta,all=TRUE))
  }
  } else {
      val <- obj(0,all=TRUE)
  }

  if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X)
  val <- c(val,list(time=time,formula=formula,formC=formC,
    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid))

  val$call <- cl
  val$MGciid <- MGCiid
  val$call.id <- call.id
  val$id <- orig.id
  val$name.id <- name.id
  val$nid <- nid
  val$iid.naive <- val$iid
  val$naive.var <- NULL 
  if (se)  val$iid  <- val$iid+(MGCiid %*% val$ihessian)
  if (!is.null(call.id)) val$iid <- namesortme(val$iid,name.id)
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

  class(val) <- "binreg"
  return(val)
}# }}}

##' @export
IC.binreg  <- function(x,...) {# {{{
  x$iid*NROW(x$iid)
}# }}}

##' @export
print.binreg  <- function(x,...) {# {{{
  print(summary(x),...)
}# }}}

##' @export
summary.binreg <- function(object,...) {# {{{

if (!is.null(object$gradient)) { ## write warning if gradient not small 
gradient <- max(abs(object$gradient))
if (gradient > 0.000001) { cat("gradient:\n"); print(object$gradient) }
}

cc  <- estimate(coef=object$coef,vcov=object$var)$coefmat
V=object$var

res <- list(coef=cc,n=object$n,nevent=object$nevent,strata=NULL,
	    ncluster=object$ncluster,var=V,model=object$model[1])

if (object$model[1]=="exp" | object$model[1]=="logit") {
expC <- exp(lava::estimate(coef=coef(object),vcov=object$var)$coefmat[,c(1,3,4),drop=FALSE])
res <- c(res,list(exp.coef=expC))
}

## to add marginal estimates for binregATE 
if (!is.null(object$riskDR))  {
    marginalDR <- estimate(coef=object$riskDR,vcov=object$var.riskDR)$coefmat
    difmarginalDR <- estimate(coef=object$difriskDR,vcov=as.matrix(object$var.difriskDR))$coefmat
###    rownames(difmarginalDR) <- "differenceDR"
    marginalDR <- rbind(marginalDR,difmarginalDR)

    marginalG <- estimate(coef=object$riskG,vcov=object$var.riskG)$coefmat
    difG <- estimate(coef=object$difriskG,vcov=as.matrix(object$var.difriskG))$coefmat
###    rownames(difG) <- "differenceG"
    marginalG <- rbind(marginalG,difG) 
    res <- c(res,list(ateG=marginalG,ateDR=marginalDR))

    if (!is.null(object$attc)) {
    attc <- estimate(coef=object$attc,vcov=object$var.attc)$coefmat
    res <- c(res,list(attc=attc))
    }
}

class(res) <- "summary.binreg"
return(res)
}# }}}

##' @export
print.summary.binreg <- function(x,max.strata=5,...) { ## {{{

  nn <- cbind(x$n, x$nevent)
  rownames(nn) <- levels(x$strata); colnames(nn) <- c("n","events")
  if (is.null(rownames(nn))) rownames(nn) <- rep("",NROW(nn))
  if (length(x$strata)>max.strata) {
      nn <- rbind(c(colSums(nn),length(x$strata)));
      colnames(nn) <- c("n","events","stratas")
      rownames(nn) <- ""
  } 
  print(nn,quote=FALSE)  
  if (!is.null(x$ncluster)) cat("\n ", x$ncluster, " clusters\n",sep="")
  if (!is.null(x$coef)) {
    cat("coeffients:\n")
    printCoefmat(x$coef,...)
    cat("\n")
    if (x$model[1]=="exp" | x$model[1]=="logit") {
    cat("exp(coeffients):\n")
    printCoefmat(x$exp.coef,...)
    }
  }
  cat("\n")

 ## for binreg ATE
 if (!is.null(x$ateDR)) {
    cat("Average Treatment effects (G-formula) :\n")
    printCoefmat(x$ateG,...)
    cat("\n")

    cat("Average Treatment effects (double robust) :\n")
    printCoefmat(x$ateDR,...)
    cat("\n")

###    if (!is.null(x$attc)) {
###    cat("Average Treatment effects on Treated/Non-Treated (DR) :\n")
###    printCoefmat(x$attc,...)
###    cat("\n")
###    }

  }
  cat("\n")

} ## }}}

###}}} 

##' @export
vcov.binreg <- function(object,...) {# {{{
	return(object$var)
}# }}}

##' @export
coef.binreg <- function(object,...) {# {{{
	return(object$coef)
}# }}}

##' @export
predict.binreg <- function(object,newdata,se=TRUE,iid=FALSE,...)
{# {{{

  xlev <- lapply(object$model.frame,levels)
  ff <- unlist(lapply(object$model.frame,is.factor))
  upf <- update(object$formula,~.)
  tt <- terms(upf)
  tt <- delete.response(tt)
  Z <- model.matrix(tt,data=newdata,xlev=xlev)
  ## assuming that cluster comes after Z's 
  Z <- as.matrix(Z)[,1:length(object$coef),drop=FALSE]
  clusterTerm<- grep("^cluster[(][A-z0-9._:]*",colnames(object$model.frame),perl=TRUE)
  if (length(clusterTerm)>=1) 
	  if (ncol(object$model.frame)!=clusterTerm) stop("cluster term must be last\n")

  if (object$model[1]=="logit") {
  lp <- c(Z %*% object$coef)
  p <- expit(lp)
  preds <- p

  if (se) {
     if (is.null(object$var)) covv <- vcov(object)  else covv <- object$var
     Dpv <- Z*exp(-lp)*p^2
     se <- apply((Dpv %*% covv)* Dpv,1,sum)^.5
     cmat <- data.frame(pred=p,se=se,lower=p-1.96*se,upper=p+1.96*se)
     names(cmat)[1:4] <- c("pred","se","lower","upper")
     preds <- cmat
  }
  if (iid) {
      Piid <- object$iid  %*% t(Dpv)
      preds <- list(pred=preds,iid=Piid)
  }

  } else {

  lp <- c(Z %*% object$coef)
  if (object$model[1]=="exp") p <- exp(lp) else p <- lp
  preds <- p

  if (se) {
     if (is.null(object$var)) covv <- vcov(object)  else covv <- object$var
###     if (object$model[1]=="dexp") Dpv <- Z*p^2 else 
     if (object$model[1]=="exp") Dpv <- Z*p else Dpv <- Z 
     se <- apply((Dpv %*% covv)* Dpv,1,sum)^.5
     cmat <- data.frame(pred=p,se=se,lower=p-1.96*se,upper=p+1.96*se)
     names(cmat)[1:4] <- c("pred","se","lower","upper")
     preds <- cmat
  }
  if (iid) {
      Piid <- object$iid  %*% t(Dpv)
      preds <- list(pred=preds,iid=Piid)
  }

  }
return(preds)
} # }}}

##' @export
binregt <- function(formula,data,cause=1,time=NULL,beta=NULL,
	   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
	   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,...)
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

  call.id <- id
  conid <- construct_id(id,nrow(X),namesX=rownames(X))
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
  ## id before time-sorting
  orig.id <- id

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- (status %in% cens.code) 
  statusE <- (status %in% cause)*outer(exit,time,"<=") 
  if (any(apply(statusE,2,sum))==0) stop("No events of type 1 before time \n"); 
  kmt <- kaplan.meier

  statusC <- (status %in% cens.code) 
  data$id__ <- id
  data$exit <- exit
  data$statusC <- statusC 
  cens.strata <- cens.nstrata <- NULL 

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      cens.weights <- c()
      for (tt in time)  {
      exittime <- pmin(exit,tt)
      cens.weights1 <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      cens.weights <- cbind(cens.weights,cens.weights1)
      }
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  if (length(time)>1) {
      Yt <- outer(exit,time,"<=") 
      nevent <- apply(c(status %in% cause)*Yt,2,sum)
      Y <- c((status %in%cause)/cens.weights)*Yt
  } else {
      Y <- c((status %in%cause)*(exit<=time)/cens.weights)
      nevent <- sum((status %in%cause)*(exit<=time))
  }

  ppt <- length(time)
  p <- pd <- ncol(X)+(length(time)-1)
  if (is.null(beta)) beta <- rep(0,ncol(X)+ppt-1)
  X <-  as.matrix(X)
  X2o  <- .Call("vecMatMat",X[,-1,drop=FALSE],X[,-1,drop=FALSE])$vXZ
###  X2  <- .Call("vecMatMat",X,X)$vXZ

  if (is.null(augmentation))  augmentation=rep(0,p)

obj <- function(pp,all=FALSE)
{ # {{{

lp <- c(X[,-1,drop=FALSE] %*% pp[-seq(1,ppt)])
lp <- outer(lp,pp[1:ppt],"+")
p <- expit(lp)
ploglik <- sum(weights*(Y-p)^2)

Dloglt <- weights*(Y-p)
Dloglo <- weights*X[,-1]*apply((Y-p),1,sum)
Dlogl <- cbind(Dloglt,Dloglo)

wwp <- (weights*p/(1+exp(lp)))
wwpt <- apply(wwp,2,sum)
wwpo <- apply(wwp,1,sum)
D2loglo <- c(wwpo)*X2o
hessian <- matrix(0,pd,pd)
hessian[seq(1,ppt),seq(1,ppt)] <- diag(wwpt,ppt,ppt)
for (i in 1:ppt) hessian[i,(ppt+1):pd] <- apply(wwp[,i]*X[,-1,drop=FALSE],2,sum)
hessian[(ppt+1):pd,1:ppt]  <- t(hessian[1:ppt,(ppt+1):pd])
hessian[(ppt+1):pd,(ppt+1):pd] <- c(apply(D2loglo,2,sum))

gradient <- apply(Dlogl,2,sum)+augmentation

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      beta.iid <-  apply(beta.iid,2,sumstrata,id,max(id)+1)
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
	 id=id,Dlogl=Dlogl,iid=beta.iid,robvar=robvar,var=robvar,se.robust=diag(robvar)^.5)
      return(val)
  }  
 structure(-ploglik/nid,gradient=-gradient/nid,hessian=hessian/nid)
}# }}}

	  p <- pd 
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

	  if (length(time)==1) {
	  if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X) }
	  else names(val$coef) <- c(paste("Intercept",time,sep=""),colnames(X)[-1])
	  val <- c(val,list(time=time,formula=formula,formC=formC,
	    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
	    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid))
	  

 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    cens.weights <- cens.weights[ord]
    if (length(time)>1) {
         Yt <- outer(exit,time,"<=") 
         Y <- c((status%in% cause)/cens.weights)*Yt
	 Ys <- apply(Y,1,sum)
     } else {
         Ys <- Y <- c((status%in%cause)*(exit<=time)/cens.weights)
     }

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(cbind(Y,X[,-1]*Ys),2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
###    btime <- 1*outer(exit,time,"<") 
    IhdLam0 <- apply(h*c(S0i2),2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(X),pd)
    jtime <- outer(resC$jumptimes,time,"<")*1
    U[xx$jumps+1,] <- (h[xx$jumps+1,])/c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    MGCiid <- apply(MGt,2,sumstrata,xx$id,max(id)+1)
 
   }  else {
	  MGCiid <- 0
  }## }}}

  val$call <- cl
  val$MGciid <- MGCiid
  val$call.id <- call.id
  val$name.id <- name.id
  val$.id <- orig.id
  val$nid <- nid
  val$iid.naive <- val$iid 
  if (se) val$iid  <- val$iid+(MGCiid %*% val$ihessian) else val$iid  <- val$iid
  val$naive.var <- val$var
  robvar <- crossprod(val$iid)
  val$var <-  val$robvar <- robvar
  val$se.robust <- diag(robvar)^.5
  val$se.coef <- diag(val$var)^.5
  val$cause <- cause
  val$cens.code <- cens.code 
  val$model <- "logit"


  class(val) <- "binreg"
  return(val)
}# }}}

##' @export
logitIPCW <- function(formula,data,cause=1,time=NULL,beta=NULL,
	   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
	   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
	   outcome=c("cif","rmst","rmtl"),model=c("default","logit","exp","lin"),Ydirect=NULL,...)
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

  call.id <- id
  conid <- construct_id(id,nrow(X),namesX=rownames(X))
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
  orig.id <- id

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- (status%in% cens.code) 
  statusE <- (status%in%cause) & (exit<= time) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 
  kmt <- kaplan.meier

  statusC <- (status %in%cens.code) 
  data$id__ <- id
  data$exit <- exit
  data$statusC <- statusC 
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 

  cens.strata <- cens.nstrata <- NULL 

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  if (is.null(beta)) beta <- rep(0,ncol(X))
  p <- ncol(X)

  X <-  as.matrix(X)
  X2  <- .Call("vecCPMat",X)$XX
  obs <- (exit<=time & status!=cens.code) | (exit>=time)
  weights <- obs*weights/c(cens.weights)
  cens.weights <- c(cens.weights)

  if (!is.null(Ydirect)) Y <-  Ydirect else {
     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time))
     else { if (!competing) {
	     if (outcome[1]=="rmst")
	     Y <-  c(pmin(exit,time)) 
             else Y <-  c((time-pmin(exit,time)))
            } else Y <- c((status %in% cause)*(time-pmin(exit,time)))
     }
  }
  Yipcw <- Y
 
 if (is.null(augmentation))  augmentation=rep(0,p)
 nevent <- sum(Y)

 ## default links
 if (model[1]=="default") {
	 if (outcome[1]=="cif") model <- "logit"
	 if (outcome[1]=="rmst") model <- "exp"
	 if (outcome[1]=="rmtl") model <- "exp"
	 if (outcome[1]=="years-lost") model <- "exp"
 }


obj <- function(pp,all=FALSE)
{ # {{{

lp <- c(X %*% pp+offset)

    if (model[1]=="exp") {
	 p <- exp(lp) 
         D2logl <- c(weights*p)*X2 
     } else if (model[1]=="lin") {
	 p <- lp
         D2logl <- c(weights)*X2
     } else if (model[1]=="logit") {
	p <- expit(lp)
        D2logl <- c(weights*p/(1+exp(lp)))*X2
     } else stop("link functions must be logit,exp,lin\n") 
ploglik <- sum(weights*(Y-p)^2)

if (model[1]=="exp") ploglik <- 0
Dlogl <- weights*X*c(Y-p)
D2log <- apply(D2logl,2,sum)
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
 structure(-ploglik,gradient=-gradient,hessian=hessian)
}# }}}

## setting default for NR 
dots <- list(...)
if (length(dots)==0) {
   if (model[1]=="exp") control <- list(tol=1e-10,stepsize=0.5)  
   else control <- NULL
} else control <- dots[[1]]

  p <- ncol(X)
  opt <- NULL
  if (p>0) {
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
	  tim <- system.time(opt <- lava::NR(beta,obj,control=control))
	  opt$timing <- tim
	  opt$estimate <- opt$par
      } else {
	  opt <- nlm(obj,beta,...)
	  opt$method <- "nlm"
      }
      cc <- opt$estimate; 
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }

  if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X)
  val <- c(val,list(time=time,formula=formula,formC=formC,
    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid,weights=weights))
  
 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    lp <- c(X %*% val$coef+offset)
    p <- expit(lp)
    Yglm <- weights*c(Y[ord]-p) # *(exit<=time)

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Yglm,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply((exit<time)*h*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*h[xx$jumps+1,] /c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    MGCiid <- apply(MGt,2,sumstrata,xx$id,max(id)+1)
  }  else {
	  MGCiid <- 0
  }## }}}

    val$Y <- Yipcw
    val$model <- model[1]
    val$outcome <- outcome[1]
    val$call <- cl
    val$MGciid <- MGCiid
    val$id <- orig.id
    val$call.id <- call.id
    val$name.id <- name.id
    val$nid <- nid
    val$iid.naive <- val$iid 
    val$iid  <- val$iid+(MGCiid %*% val$ihessian)
    if (!is.null(call.id)) val$iid <- namesortme(val$iid,name.id)
    val$naive.var <- val$var
    robvar <- crossprod(val$iid)
    val$var <-  val$robvar <- robvar
    val$se.robust <- diag(robvar)^.5
    val$se.coef <- diag(val$var)^.5
    val$cause <- cause
    val$cens.code <- cens.code 
    val$augmentation <- augmentation

  class(val) <- "binreg"
  return(val)
}# }}}

logitIPCW <- function(formula,data,cause=1,time=NULL,beta=NULL,
	   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
	   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
	   outcome=c("cif","rmst","rmtl"),model=c("default","logit","exp","lin"),Ydirect=NULL,...)
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

  call.id <- id
  conid <- construct_id(id,nrow(X),namesX=rownames(X))
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
  orig.id <- id

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- (status%in% cens.code) 
  statusE <- (status%in%cause) & (exit<= time) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 
  kmt <- kaplan.meier

  statusC <- (status %in%cens.code) 
  data$id__ <- id
  data$exit <- exit
  data$statusC <- statusC 
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 

  cens.strata <- cens.nstrata <- NULL 

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  if (is.null(beta)) beta <- rep(0,ncol(X))
  p <- ncol(X)

  X <-  as.matrix(X)
  X2  <- .Call("vecCPMat",X)$XX
  obs <- (exit<=time & status!=cens.code) | (exit>=time)
  weights <- obs*weights/c(cens.weights)
  cens.weights <- c(cens.weights)

  if (!is.null(Ydirect)) Y <-  Ydirect else {
     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time))
     else { if (!competing) {
	     if (outcome[1]=="rmst")
	     Y <-  c(pmin(exit,time)) 
             else Y <-  c((time-pmin(exit,time)))
            } else Y <- c((status %in% cause)*(time-pmin(exit,time)))
     }
  }
  Yipcw <- Y
 
 if (is.null(augmentation))  augmentation=rep(0,p)
 nevent <- sum(Y)

 ## default links
 if (model[1]=="default") {
	 if (outcome[1]=="cif") model <- "logit"
	 if (outcome[1]=="rmst") model <- "exp"
	 if (outcome[1]=="rmtl") model <- "exp"
	 if (outcome[1]=="years-lost") model <- "exp"
 }


obj <- function(pp,all=FALSE)
{ # {{{

lp <- c(X %*% pp+offset)

    if (model[1]=="exp") {
	 p <- exp(lp) 
         D2logl <- c(weights*p)*X2 
     } else if (model[1]=="lin") {
	 p <- lp
         D2logl <- c(weights)*X2
     } else if (model[1]=="logit") {
	p <- expit(lp)
        D2logl <- c(weights*p/(1+exp(lp)))*X2
     } else stop("link functions must be logit,exp,lin\n") 
ploglik <- sum(weights*(Y-p)^2)

if (model[1]=="exp") ploglik <- 0
Dlogl <- weights*X*c(Y-p)
D2log <- apply(D2logl,2,sum)
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
 structure(-ploglik,gradient=-gradient,hessian=hessian)
}# }}}

## setting default for NR 
dots <- list(...)
if (length(dots)==0) {
   if (model[1]=="exp") control <- list(tol=1e-10,stepsize=0.5)  
   else control <- NULL
} else control <- dots[[1]]

  p <- ncol(X)
  opt <- NULL
  if (p>0) {
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
	  tim <- system.time(opt <- lava::NR(beta,obj,control=control))
	  opt$timing <- tim
	  opt$estimate <- opt$par
      } else {
	  opt <- nlm(obj,beta,...)
	  opt$method <- "nlm"
      }
      cc <- opt$estimate; 
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }

  if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X)
  val <- c(val,list(time=time,formula=formula,formC=formC,
    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid,weights=weights))
  
 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    lp <- c(X %*% val$coef+offset)
    if (model[1]=="exp") { p <- exp(lp) } else if (model[1]=="lin") { p <- lp } else if (model[1]=="logit") {p <- expit(lp) } 
    Yglm <- weights*c(Y[ord]-p) # *(exit<=time)

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Yglm,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply((exit<time)*h*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*h[xx$jumps+1,] /c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    MGCiid <- apply(MGt,2,sumstrata,xx$id,max(id)+1)
  }  else {
	  MGCiid <- 0
  }## }}}

    val$Y <- Yipcw
    val$model <- model[1]
    val$outcome <- outcome[1]
    val$call <- cl
    val$MGciid <- MGCiid
    val$id <- orig.id
    val$call.id <- call.id
    val$name.id <- name.id
    val$nid <- nid
    val$iid.naive <- val$iid 
    val$iid  <- val$iid+(MGCiid %*% val$ihessian)
    if (!is.null(call.id)) val$iid <- namesortme(val$iid,name.id)
    val$naive.var <- val$var
    robvar <- crossprod(val$iid)
    val$var <-  val$robvar <- robvar
    val$se.robust <- diag(robvar)^.5
    val$se.coef <- diag(val$var)^.5
    val$cause <- cause
    val$cens.code <- cens.code 
    val$augmentation <- augmentation

  class(val) <- "binreg"
  return(val)
}# }}}

##' Average Treatment effect for censored competing risks data using Binomial Regression 
##'
##' Under the standard causal assumptions  we can estimate the average treatment effect E(Y(1) - Y(0)). We need Consistency, ignorability ( Y(1), Y(0) indep A given X), and positivity.
##'
##' The first covariate in the specification of the competing risks regression model must be the treatment variable that should be coded as a factor. 
##' If the factor has more than two levels
##' then it uses the mlogit for propensity score modelling. If there are no censorings this is the same as ordinary logistic regression modelling. 
##'
##' Estimates the ATE using the the standard binary double robust estimating equations that are IPCW censoring adjusted.
##' Rather than binomial regression we also consider a IPCW weighted version of standard logistic regression logitIPCWATE. 
##'
##' typeATE="II" will augment the estimating equation with \deqn{ (A/\pi(X)) \int E( O(t) | T \geq t, S(X))/ G_c(t,S(X)) d \hat M_c(s) } when estimating the mean
##' outcome for the treated. 
##'
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param cause cause of interest
##' @param time  time of interest 
##' @param beta starting values 
##' @param treat.model logistic treatment model given covariates 
##' @param cens.model only stratified cox model without covariates
##' @param offset offsets for partial likelihood 
##' @param weights for score equations 
##' @param cens.weights censoring weights 
##' @param se to compute se's with IPCW  adjustment, otherwise assumes that IPCW weights are known
##' @param type "II" adds augmentation term, and "I" classic binomial regression 
##' @param kaplan.meier uses Kaplan-Meier for IPCW in contrast to exp(-Baseline)
##' @param cens.code gives censoring code
##' @param no.opt to not optimize 
##' @param method for optimization 
##' @param augmentation for augment binomial regression 
##' @param outcome can do CIF regression "cif"=F(t|X), "rmst"=E( min(T, t) | X) , or E( I(epsilon==cause) ( t - mint(T,t)) ) | X) depending on the number of the number of causes. 
##' @param model exp or linear model for E( min(T, t) | X)=exp(X^t beta), or E( I(epsilon==cause) ( t - mint(T,t)) ) | X)=exp(X^t beta) 
##' @param Ydirect use this outcome Y with IPCW vesion
##' @param typeATE "II" to censor augment  the estimating equation
##' @param ... Additional arguments to lower level funtions (binreg that fits outcome model)
##' @author Thomas Scheike
##' @examples
##' library(mets); data(bmt)
##' dfactor(bmt)  <-  ~.
##'
##' brs <- binregATE(Event(time,cause)~tcell.f+platelet+age,bmt,time=50,cause=1,
##'	  treat.model=tcell.f~platelet+age)
##' summary(brs)
##'
##' brsi <- binregATE(Event(time,cause)~tcell.f+tcell.f*platelet+tcell.f*age,bmt,time=50,cause=1,
##'	  treat.model=tcell.f~platelet+age)
##' summary(brsi)
##'
##'
##' @aliases logitIPCWATE logitATE normalATE kumarsim kumarsimRCT 
##' @export
binregATE <- function(formula,data,cause=1,time=NULL,beta=NULL,treat.model=~+1,cens.model=~+1,
   offset=NULL,weights=NULL,cens.weights=NULL,se=TRUE,type=c("II","I"),
   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
   outcome=c("cif","rmst","rmtl"),model=c("default","logit","exp","lin"),Ydirect=NULL,typeATE="II",...)
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
 conid <- construct_id(id,nrow(X),namesX=rownames(X))
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
 idclust <- id; nclust <- nid
 if (nid!=length(id)) { nid <- length(id); id <- 0:(nid-1); }

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  ## {{{ setting up status, response, censoring model 
  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- (status %in%cens.code) 
  statusE <- (status %in% cause) & (exit<= time) 
  if (sum(statusE)==0) warning("No events of type 1 before time \n"); 
  kmt <- kaplan.meier
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 

  statusC <- (status %in%cens.code) 
  data$id__ <- id 
  data$exit <- exit
  data$statusC <- statusC 
  n <- length(exit)
  cens.strata <- cens.nstrata <- NULL 

  call.cw <- cens.weights
  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  p <- ncol(X)
  if (is.null(beta)) beta <- rep(0,p)
  X <-  as.matrix(X)
  ###  X2  <- .Call("vecMatMat",X,X)$vXZ
  X2  <- .Call("vecCPMat",X)$XX
  ###
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)>=1) Causes <- ucauses[-ccc] else Causes <- ucauses
  obs <- (exit<=time & (status %in% Causes)) | (exit>=time)

 if (!is.null(Ydirect)) Y <-  Ydirect*obs/cens.weights else {
     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time)/cens.weights)
     else { if (!competing) {
	     if (outcome[1]=="rmst") Y <-  c(pmin(exit,time)*obs)/cens.weights 
	     else Y <-  c((time-pmin(exit,time))*obs)/cens.weights 
            } else Y <- c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights
     }
  }
 nevent <- sum((status %in% cause)*(exit<=time))

 ## }}}

## change id from call to id__
formula <- drop.specials(formula,"cluster")
formula <- update(formula, .~.+cluster(id__))

  val <- binreg(formula,data,cause=cause,time=time,beta=beta,type=type,
	cens.model=cens.model,se=se,kaplan.meier=kaplan.meier,
        offset=offset,weights=weights,cens.weights=call.cw,
        cens.code=cens.code,no.opt=no.opt,method=method,augmentation=augmentation,
        outcome=outcome[1],model=model[1],Ydirect=Ydirect,...)

model <- val$model[1]

# {{{ computation of ate, att, atc and their influence functions

### treatment is rhs of treat.model 
if (length(all.vars(treat.model))==0) 
{ respformula  <-  as.formula( paste(all.vars(formula)[3],"~."))
  treat.model  <- update.formula(treat.model,respformula)
###  cat("Using propensity score model :\n "); print(treat.model)
}
treat.name <-  all.vars(treat.model)[1]
treatvar <- data[,treat.name]
if (!is.factor(treatvar)) stop(paste("treatment=",treat.name," must be coded as factor \n",sep="")); 
## treatvar, 1,2,...,nlev or 1,2
nlev <- nlevels(treatvar)
nlevs <- levels(treatvar)
###treatvar <- as.numeric(treatvar)
ntreatvar <- as.numeric(treatvar)
ytreat <- ntreatvar-1

## dropping cluster here 
if (nlev==2) {
   treat.model <- drop.specials(treat.model,"cluster")
   treat <- glm(treat.model,data,family="binomial")
   iidalpha <- iid(treat,id=id)
   lpa <- treat$linear.predictors 
   pal <- expit(lpa)
   pal <-cbind(1-pal,pal)
   ppp <- (pal/pal[,1])
   spp <- 1/pal[,1]
} else {  
   treat.modelid <- update.formula(treat.model,.~.+cluster(id__))
   treat <- mlogit(treat.modelid,data)
   iidalpha <- iid(treat)
   pal <- predict(treat,data,se=0,response=FALSE)
   ppp <- (pal/pal[,1])
   spp <- 1/pal[,1]
}

   ###########################################################
   ### computes derivative of D (1/Pa) propensity score 
   ###########################################################
   Xtreat <- model.matrix(treat.model,data)
   tvg2 <- 1*(ntreatvar>=2)
   pA <- c(mdi(pal, 1:length(treatvar), ntreatvar))
   pppy <- c(mdi(ppp,1:length(treatvar), ntreatvar))
   Dppy <-  (spp*tvg2-pppy) 
   Dp <- c()
   for (i in seq(nlev-1)) Dp <- cbind(Dp,Xtreat*ppp[,i+1]*Dppy/spp^2);  
   DPai <- -1*Dp/pA^2
   p1lp <-   X %*% val$coef+offset
   if (model[1]=="logit") { p1 <- expit(p1lp) } else {
    if (model[1]=="exp") { p1 <- exp(p1lp); } else { p1 <- p1lp;}
   }


k <- 0
DePsia <- DariskG <- DaPsia <- list(); 
pia <- Ya <- riskG <- riska <- c(); 
datA <- dkeep(data,x=all.vars(formula))
xlev <- lapply(datA,levels)
formulanc <- drop.specials(formula,"cluster")
a <- nlevs[1]
for (a in nlevs) {# {{{
	k <- k+1
	datA[,treat.name] <- a
	Xa <- model.matrix(formulanc[-2],datA,xlev=xlev)
        lpa <- Xa %*% val$coef+offset
	if (model[1]=="logit") {
	   ma <- expit(lpa); Dma  <-  Xa*c(ma/(1+exp(lpa)))
	} else {
	    if (model[1]=="exp") { ma <- exp(lpa);  Dma<-c(ma)*Xa; } else { ma <- lpa; Dma <- Xa }
	}
	paka <- pal[,k]
	riska <- cbind(riska,((treatvar==a)/paka)*(Y-ma)+ma)
	riskG <- cbind(riskG,ma)
	Ya <- cbind(Ya,Y*((treatvar==a)/paka))
	pia <- cbind(pia,((treatvar==a)/paka))
        DariskG[[k]] <- apply(Dma,2,sum)
        DePsia[[k]] <-  apply(Dma*(1-(treatvar==a)/paka),2,sum)
        DaPsia[[k]] <-  apply(DPai*(treatvar==a)*c(Y-p1),2,sum)
}# }}}

augmentationATE <- 0
 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    Xtreat <- Xtreat[ord,,drop=FALSE]
    ytreat <- ytreat[ord]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    Ya <- Ya[ord,]
    pal <- pal[ord]
    piaC <- pia[ord,]
    cens.weights <- cens.weights[ord]
    lp <- c(X %*% val$coef+offset)
    obs <- obs[ord]
    p <- expit(lp)

  Y <- Y[ord]
  Yc <- matrix(Y,ncol=1)

 Y <- Y*weights 
 xx <- resC$cox.prep
 S0i2 <- S0i <- rep(0,length(xx$strata))
 S0i[xx$jumps+1]  <- 1/resC$S0
 S0i2[xx$jumps+1] <- 1/resC$S0^2
 ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
 ## to make \int h(s)/Ys  dM_i^C(s) 
 has  <-  apply(Ya,2,revcumsumstrata,xx$strata,xx$nstrata)
###    hattc  <-  apply(cbind(ytreat-pal*(1-ytreat)/(1-pal),-(1-ytreat)+(1-pal)*ytreat/pal)*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    btime <- 1*(exit<time)
    IhdLamhas <- apply(has*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    Uas <- matrix(0,nrow(xx$X),ncol(has))
    Uas[xx$jumps+1,] <- (resC$jumptimes<time)*has[xx$jumps+1,] /c(resC$S0)
    MGtas <- (Uas[,drop=FALSE]-IhdLamhas)*c(xx$weights)

    if (typeATE=="II") {
            h  <-  apply(Yc,2,revcumsumstrata,xx$strata,xx$nstrata)
	    IhdLam <- apply(h*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
	    U <- matrix(0,nrow(xx$X),1)
	    U[xx$jumps+1,] <- (resC$jumptimes<time)*h[xx$jumps+1,]/c(resC$S0)
	    MGateII <- piaC*c(U[,drop=FALSE]-IhdLam)*c(xx$weights)
	    augmentationATE <- apply(MGateII,2,sum)

	    hapi  <-  apply(piaC,2,revcumsumstrata,xx$strata,xx$nstrata)
	    IhaidLam <- apply(hapi*c(h)*S0i*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
	    Uai <- matrix(0,nrow(xx$X),ncol(hapi))
	    Uai[xx$jumps+1,] <- (resC$jumptimes<time)*c(h[xx$jumps+1,])*hapi[xx$jumps+1,]/c(resC$S0)^2
	    MGait <- (Uai[,drop=FALSE]-IhaidLam)*c(xx$weights)
	    MGtII <- MGateII - MGait
    }
    ### Censoring Variance Adjustment \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    mid <- max(xx$id)+1
    MGCiidas <- apply(MGtas,2,sumstrata,xx$id,mid)
    if (typeATE=="II") {
        MGtII <-  apply(MGtII,2,sumstrata,xx$id,mid)
        MGCiidas  <-  MGCiidas+MGtII
    } else MGtII <- 0
 
  }  else { MGCiid <- MGCiidas <- 0 }
## }}}

val$call.id <- call.id
val$name.id  <- name.id
val$id <- id
val$nid  <- nid
val$se.coef <- diag(val$var)^.5

################################
### estimates risk, att, atc
################################
###val$riskDR <- apply(riska,2,mean)
val$riskDR <- apply(riska,2,mean)+augmentationATE/nid
val$riskG<-  apply(riskG,2,mean)
names(val$riskDR) <- paste("treat",nlevs,sep="")
names(val$riskG) <- paste("treat",nlevs,sep="")

################################
## iid's of marginal risk estimates 
################################

k <- 1
iidrisk <- c()
riskG.iid <- c()
for (a in nlevs) {
	iidbasea <- c(sumstrata(riska[,k]-val$riskDR[k],id,nid))
	iidcifa <- c(DePsia[[k]] %*% t(val$iid))
	iidpala <- c(DaPsia[[k]] %*% t(iidalpha))
	if (se)  iidGca <- MGCiidas[,k] else iidGca<-0 
        ###
	iidriska <- (iidbasea+iidcifa+iidpala+iidGca)/nid
        iidrisk <- cbind(iidrisk,iidriska)
	iidriskG <- c(sumstrata(riskG[,k]-val$riskG[k],id,nid))
	riskGa.iid <- c(iidriskG)/nid+c(DariskG[[k]] %*% t(val$iid))/nid
        riskG.iid <- cbind(riskG.iid,riskGa.iid)
	k <- k+1
}
# }}}

# {{{ output variances and se for ate; cluster correction
val$id <- idclust
val$nid  <- nclust

val$iid <- apply(val$iid,2,sumstrata,idclust,nclust)
robvar <- crossprod(val$iid)
val$var <-  val$robvar <- robvar
val$se.robust <- diag(robvar)^.5
val$se.coef <- diag(val$var)^.5
val$MGciid <- NULL
val$iid.naive <- NULL

## outcome model from binreg
iidrisk <- apply(iidrisk,2,sumstrata,idclust,nclust)
val$var.riskDR <- crossprod(iidrisk); 
val$se.riskDR <- diag(val$var.riskDR)^.5
val$riskDR.iid <- iidrisk

pdiff <- function(x) lava::contr(lapply(seq(x-1), function(z) seq(z,x)))
contrast <- -pdiff(length(nlevs))
nncont <- c()
for (k in seq_along(nlevs[-length(nlevs)])) nncont <-c(nncont, paste("treat:",nlevs[-seq(k)],"-",nlevs[k],sep="")) 
rownames(contrast) <- nncont

mm <- estimate(coef=val$riskDR,vcov=val$var.riskDR,contrast=contrast)
val$difriskDR <- mm$coef 
names(val$difriskDR) <-  rownames(contrast) 
val$var.difriskDR <- mm$vcov 
val$se.difriskDR <- diag(val$var.difriskDR)^.5

riskG.iid <- apply(riskG.iid,2,sumstrata,idclust,nclust)
###
val$riskG.iid <- riskG.iid
val$var.riskG <- crossprod(val$riskG.iid)
val$se.riskG <- diag(val$var.riskG)^.5
###
mm <- estimate(coef=val$riskG,vcov=val$var.riskG,contrast=contrast)
val$difriskG <- mm$coef 
names(val$difriskG) <-  rownames(contrast) 
val$var.difriskG <- mm$vcov 
val$se.difriskG <- diag(val$var.difriskG)^.5

### DR-estimator ; G -estimator sort after namid; also outcome model
if (!is.null(call.id)) {
	iidrisk <- namesortme(iidrisk,name.id)
	riskG.iid <-  namesortme(riskG.iid,name.id)
	val$iid <-  namesortme(val$iid,name.id)
}
## }}}

  class(val) <- "binreg"
  return(val)
}# }}}

binregATEold <- function(formula,data,cause=1,time=NULL,beta=NULL,treat.model=~+1,cens.model=~+1,
   offset=NULL,weights=NULL,cens.weights=NULL,se=TRUE,type=c("II","I"),
   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
   outcome=c("cif","rmst","rmtl"),model=c("default","logit","exp","lin"),Ydirect=NULL,...)
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
 conid <- construct_id(id,nrow(X),namesX=rownames(X))
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
 idclust <- id; nclust <- nid
 if (nid!=length(id)) { nid <- length(id); id <- 0:(nid-1); }

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  ## {{{ setting up status, response, censoring model 
  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- (status %in%cens.code) 
  statusE <- (status %in% cause) & (exit<= time) 
  if (sum(statusE)==0) warning("No events of type 1 before time \n"); 
  kmt <- kaplan.meier
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 

  statusC <- (status %in%cens.code) 
  data$id__ <- id 
  data$exit <- exit
  data$statusC <- statusC 
  n <- length(exit)
  cens.strata <- cens.nstrata <- NULL 

  call.cw <- cens.weights
  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  p <- ncol(X)
  if (is.null(beta)) beta <- rep(0,p)
  X <-  as.matrix(X)
  ###  X2  <- .Call("vecMatMat",X,X)$vXZ
  X2  <- .Call("vecCPMat",X)$XX
  ###
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)>=1) Causes <- ucauses[-ccc] else Causes <- ucauses
  obs <- (exit<=time & (status %in% Causes)) | (exit>=time)

 if (!is.null(Ydirect)) Y <-  Ydirect*obs/cens.weights else {
     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time)/cens.weights)
     else { if (!competing) {
	     if (outcome[1]=="rmst") Y <-  c(pmin(exit,time)*obs)/cens.weights 
	     else Y <-  c((time-pmin(exit,time))*obs)/cens.weights 
            } else Y <- c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights
     }
  }
 nevent <- sum((status %in% cause)*(exit<=time))

 ## }}}

## change id from call to id__
formula <- drop.specials(formula,"cluster")
formula <- update(formula, .~.+cluster(id__))

  val <- binreg(formula,data,cause=cause,time=time,beta=beta,type=type,
	cens.model=cens.model,se=se,kaplan.meier=kaplan.meier,
        offset=offset,weights=weights,cens.weights=call.cw,
        cens.code=cens.code,no.opt=no.opt,method=method,augmentation=augmentation,
        outcome=outcome[1],model=model[1],Ydirect=Ydirect,...)

model <- val$model[1]

# {{{ computation of ate, att, atc and their influence functions

### treatment is rhs of treat.model 
if (length(all.vars(treat.model))==0) 
{ respformula  <-  as.formula( paste(all.vars(formula)[3],"~."))
  treat.model  <- update.formula(treat.model,respformula)
###  cat("Using propensity score model :\n "); print(treat.model)
}
treat.name <-  all.vars(treat.model)[1]
treatvar <- data[,treat.name]
if (!is.factor(treatvar)) stop(paste("treatment=",treat.name," must be coded as factor \n",sep="")); 
## treatvar, 1,2,...,nlev or 1,2
nlev <- nlevels(treatvar)
nlevs <- levels(treatvar)
###treatvar <- as.numeric(treatvar)
ntreatvar <- as.numeric(treatvar)
ytreat <- ntreatvar-1

## dropping cluster here 
if (nlev==2) {
   treat.model <- drop.specials(treat.model,"cluster")
   treat <- glm(treat.model,data,family="binomial")
   iidalpha <- iid(treat,id=id)
   lpa <- treat$linear.predictors 
   pal <- expit(lpa)
   pal <-cbind(1-pal,pal)
   ppp <- (pal/pal[,1])
   spp <- 1/pal[,1]
} else {  
   treat.modelid <- update.formula(treat.model,.~.+cluster(id__))
   treat <- mlogit(treat.modelid,data)
   iidalpha <- iid(treat)
   pal <- predict(treat,data,se=0,response=FALSE)
   ppp <- (pal/pal[,1])
   spp <- 1/pal[,1]
}


   ###########################################################
   ### computes derivative of D (1/Pa) propensity score 
   ###########################################################
   Xtreat <- model.matrix(treat.model,data)
   tvg2 <- 1*(ntreatvar>=2)
   pA <- c(mdi(pal, 1:length(treatvar), ntreatvar))
   pppy <- c(mdi(ppp,1:length(treatvar), ntreatvar))
   Dppy <-  (spp*tvg2-pppy) 
   Dp <- c()
   for (i in seq(nlev-1)) Dp <- cbind(Dp,Xtreat*ppp[,i+1]*Dppy/spp^2);  
   DPai <- -1*Dp/pA^2
   p1lp <-   X %*% val$coef+offset
   if (model[1]=="logit") { p1 <- expit(p1lp) } else {
    if (model[1]=="exp") { p1 <- exp(p1lp); } else { p1 <- p1lp;}
   }


k <- 0
DePsia <- DariskG <- DaPsia <- list(); 
Ya <- riskG <- riska <- c(); 
datA <- dkeep(data,x=all.vars(formula))
xlev <- lapply(datA,levels)
formulanc <- drop.specials(formula,"cluster")
a <- nlevs[1]
for (a in nlevs) {# {{{
	k <- k+1
	datA[,treat.name] <- a
	Xa <- model.matrix(formulanc[-2],datA,xlev=xlev)
        lpa <- Xa %*% val$coef+offset
	if (model[1]=="logit") {
	   ma <- expit(lpa); Dma  <-  Xa*c(ma/(1+exp(lpa)))
	} else {
	    if (model[1]=="exp") { ma <- exp(lpa);  Dma<-c(ma)*Xa; } else { ma <- lpa; Dma <- Xa }
	}
	paka <- pal[,k]
	riska <- cbind(riska,((treatvar==a)/paka)*(Y-ma)+ma)
	riskG <- cbind(riskG,ma)
	Ya <- cbind(Ya,Y*((treatvar==a)/paka))
        DariskG[[k]] <- apply(Dma,2,sum)
        DePsia[[k]] <-  apply(Dma*(1-(treatvar==a)/paka),2,sum)
        DaPsia[[k]] <-  apply(DPai*(treatvar==a)*c(Y-p1),2,sum)
}# }}}

 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    Xtreat <- Xtreat[ord,,drop=FALSE]
    ytreat <- ytreat[ord]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    Ya <- Ya[ord,]
    pal <- pal[ord]
    cens.weights <- cens.weights[ord]
    lp <- c(X %*% val$coef+offset)
    obs <- obs[ord]
    p <- expit(lp)

 if (!is.null(Ydirect)) Y <-  Ydirect[ord]*obs/cens.weights else {
     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time)/cens.weights)
     else { if (!competing) {
	     if (outcome[1]=="rmst") Y <-  c(pmin(exit,time)*obs)/cens.weights 
	     else Y <-  c((time-pmin(exit,time))*obs)/cens.weights 
            } else Y <- c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights
     }
  }

    Y <- Y*weights 
    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    has  <-  apply(Ya,2,revcumsumstrata,xx$strata,xx$nstrata)
###    hattc  <-  apply(cbind(ytreat-pal*(1-ytreat)/(1-pal),-(1-ytreat)+(1-pal)*ytreat/pal)*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    btime <- 1*(exit<time)
    IhdLam0 <- apply(h*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*h[xx$jumps+1,]/c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    IhdLamhas <- apply(has*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    Uas <- matrix(0,nrow(xx$X),ncol(has))
    Uas[xx$jumps+1,] <- (resC$jumptimes<time)*has[xx$jumps+1,] /c(resC$S0)
    MGtas <- (Uas[,drop=FALSE]-IhdLamhas)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    mid <- max(xx$id)+1
    MGCiid <- apply(MGt,2,sumstrata,xx$id,mid)
    MGCiidas <- apply(MGtas,2,sumstrata,xx$id,mid)
 
  }  else { MGCiid <- MGCiidas <- 0 }
## }}}

val$call.id <- call.id
val$name.id  <- name.id
val$id <- id
val$nid  <- nid
val$se.coef <- diag(val$var)^.5

################################
### estimates risk, att, atc
################################
val$riskDR <- apply(riska,2,mean)
val$riskG<-  apply(riskG,2,mean)
names(val$riskDR) <- paste("treat",nlevs,sep="")
names(val$riskG) <- paste("treat",nlevs,sep="")

################################
## iid's of marginal risk estimates 
################################

k <- 1
iidrisk <- c()
riskG.iid <- c()
for (a in nlevs) {
	iidbasea <- c(sumstrata(riska[,k]-val$riskDR[k],id,nid))
	iidcifa <- c(DePsia[[k]] %*% t(val$iid))
	iidpala <- c(DaPsia[[k]] %*% t(iidalpha))
	if (se)  iidGca <- MGCiidas[,k] else iidGca<-0 
        ###
	iidriska <- (iidbasea+iidcifa+iidpala+iidGca)/nid
        iidrisk <- cbind(iidrisk,iidriska)
	iidriskG <- c(sumstrata(riskG[,k]-val$riskG[k],id,nid))
	riskGa.iid <- c(iidriskG)/nid+c(DariskG[[k]] %*% t(val$iid))/nid
        riskG.iid <- cbind(riskG.iid,riskGa.iid)
	k <- k+1
}
# }}}

# {{{ output variances and se for ate; cluster correction
val$id <- idclust
val$nid  <- nclust

val$iid <- apply(val$iid,2,sumstrata,idclust,nclust)
robvar <- crossprod(val$iid)
val$var <-  val$robvar <- robvar
val$se.robust <- diag(robvar)^.5
val$se.coef <- diag(val$var)^.5
val$MGciid <- NULL
val$iid.naive <- NULL

## outcome model from binreg
iidrisk <- apply(iidrisk,2,sumstrata,idclust,nclust)
val$var.riskDR <- crossprod(iidrisk); 
val$se.riskDR <- diag(val$var.riskDR)^.5
val$riskDR.iid <- iidrisk

pdiff <- function(x) lava::contr(lapply(seq(x-1), function(z) seq(z,x)))
contrast <- -pdiff(length(nlevs))
nncont <- c()
for (k in seq_along(nlevs[-length(nlevs)])) nncont <-c(nncont, paste("treat:",nlevs[-seq(k)],"-",nlevs[k],sep="")) 
rownames(contrast) <- nncont

mm <- estimate(coef=val$riskDR,vcov=val$var.riskDR,contrast=contrast)
val$difriskDR <- mm$coef 
names(val$difriskDR) <-  rownames(contrast) 
val$var.difriskDR <- mm$vcov 
val$se.difriskDR <- diag(val$var.difriskDR)^.5

riskG.iid <- apply(riskG.iid,2,sumstrata,idclust,nclust)
###
val$riskG.iid <- riskG.iid
val$var.riskG <- crossprod(val$riskG.iid)
val$se.riskG <- diag(val$var.riskG)^.5
###
mm <- estimate(coef=val$riskG,vcov=val$var.riskG,contrast=contrast)
val$difriskG <- mm$coef 
names(val$difriskG) <-  rownames(contrast) 
val$var.difriskG <- mm$vcov 
val$se.difriskG <- diag(val$var.difriskG)^.5

### DR-estimator ; G -estimator sort after namid; also outcome model
if (!is.null(call.id)) {
	iidrisk <- namesortme(iidrisk,name.id)
	riskG.iid <-  namesortme(riskG.iid,name.id)
	val$iid <-  namesortme(val$iid,name.id)
}
## }}}

  class(val) <- "binreg"
  return(val)
}# }}}

##' G-estimator for binomial regression model (Standardized estimates) 
##'
##' Computes G-estimator \deqn{ \hat F(t,A=a) = n^{-1} \sum_i \hat F(t,A=a,Z_i) }.
##' Assumes that the first covariate is $A$.
##' Gives influence functions of these risk estimates and SE's are based on these.  
##' If first covariate is a factor then all contrast are computed, and if continuous 
##' then considered covariate values are given by Avalues.
##'
##' @param x binreg object
##' @param data data frame for risk averaging
##' @param Avalues values to compare for first covariate A
##' @param varname if given then averages for this variable, default is first variable
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data(bmt); bmt$time <- bmt$time+runif(408)*0.001
##' bmt$event <- (bmt$cause!=0)*1
##'
##' b1 <- binreg(Event(time,cause)~age+tcell+platelet,bmt,cause=1,time=50)
##' sb1 <- binregG(b1,bmt,Avalues=c(0,1,2))
##' summary(sb1)
##'
##' @export
binregG <- function(x,data,Avalues=c(0,1),varname=NULL)
{# {{{

if (is.null(varname))  {
treat.name <- all.vars(update.formula(x$formula,1~.))[1]
} else treat.name <- varname
treatvar <- data[,treat.name]

if (is.factor(treatvar)) {
   nlevs <- levels(treatvar)
} else {
   nlevs <- Avalues
}

formulaX <- update.formula(x$formula,.~.)
formulanc <- drop.specials(formulaX,"cluster")

datA <- dkeep(data,x=all.vars(formulaX))
xlev <- lapply(datA,levels)

DariskG <- list()
risks <- c()
k <- 0
for (a in nlevs) {# {{{
k <- k+1
datA[,treat.name] <- a
Xa <- model.matrix(formulanc[-2],datA,xlev=xlev)
lpa <- Xa %*% x$coef
## only for logit link so far 
if (x$model[1]=="logit") {
   pa <- expit(lpa)
   Dma  <-  Xa*c(pa/(1+exp(lpa)))
} else if (x$model[1]=="exp") {
   pa <- exp(lpa)
   Dma  <-  Xa*c(pa) 
} else {
   pa <- lpa
   Dma  <-  Xa
}
risks <- cbind(risks,pa)
DariskG[[k]] <- apply(Dma,2,sum)
}# }}}

Grisk <- apply(risks,2,mean)
Gest <- list(Gest=Grisk,iid=t(t(risks)-Grisk))

nid <- max(x$id)+1
n <- length(x$id)
risk.iid <- apply(Gest$iid,2,sumstrata,x$id,nid)/n 
if (!is.null(x$call.id)) risk.iid <- namesortme(risk.iid,x$name.id)
DGiid <- risk.iid

for (a in seq_along(nlevs)) {
	DGiid[,1] <- c(x$iid  %*% DariskG[[a]])/n
	risk.iid[,a] <- risk.iid[,a]+ c(x$iid  %*% DariskG[[a]])/n
}
vv <- crossprod(risk.iid)

Gout <- estimate(coef=Gest$Gest,vcov=vv,labels=paste("risk",nlevs,sep=""))
ed <-  estimate(coef=Gest$Gest,vcov=vv,f=function(p) p[-1]-p[1])
rd <- estimate(coef=Gest$Gest,vcov=vv,f=function(p) log(p[-1]/p[1]),null=0)
out <- list(risk.iid=risk.iid,risk=Gout,difference=ed,ratio=rd,vcov=vv,model=x$model[1])
class(out) <- "survivalG"
return(out)
} ## }}}

##' @export
logitIPCWATE <- function(formula,data,cause=1,time=NULL,beta=NULL,treat.model=~+1, cens.model=~+1,
	   offset=NULL,weights=NULL,cens.weights=NULL,se=TRUE,
	   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
           outcome=c("cif","rmst","rmtl"),model=c("default","logit","exp","lin"),
	   Ydirect=NULL,...)
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

### call.id <- id;
### conid <- construct_id(id,nrow(X),namesX=rownames(X))
### name.id <- conid$name.id; id <- conid$id; nid <- conid$nid

 call.id <- id;
 conid <- construct_id(id,nrow(X),namesX=rownames(X))
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
 idclust <- id; nclust <- nid
 if (nid!=length(id)) { nid <- length(id); id <- 0:(nid-1); }

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- (status %in% cens.code) 
  statusE <- (status %in% cause) & (exit<= time) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 
  kmt <- kaplan.meier
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 

  statusC <- (status %in% cens.code) 
  data$id__ <- id
  data$exit <- exit
  data$statusC <- statusC 
  n <- length(exit)

  cens.strata <- cens.nstrata <- NULL 

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
###      cens.weights <- predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt)$surv
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL

  if (is.null(beta)) beta <- rep(0,ncol(X))
  p <- ncol(X)

  X <-  as.matrix(X)
  X2  <- .Call("vecMatMat",X,X)$vXZ
  obs <- (exit<=time & status!=cens.code) | (exit>=time)
  weights <- obs*weights/c(cens.weights)
  cens.weights <- c(cens.weights)

  Y <- c((status==cause)*(exit<=time))

  if (!is.null(Ydirect)) Y <-  Ydirect else {
     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time))
     else { if (!competing) {
	     if (outcome[1]=="rmst") Y <-  c(pmin(exit,time))
	     else Y <-  c((time-pmin(exit,time)))
            } else Y <- c((status %in% cause)*(time-pmin(exit,time)))
     }
  }
 Yout <- Y
 if (is.null(augmentation))  augmentation=rep(0,p)
 nevent <- sum(Y)

 ## default links
if (model[1]=="default") {
	 if (outcome[1]=="cif") model <- "logit"
	 if (outcome[1]=="rmst") model <- "exp"
	 if (outcome[1]=="rmtl") model <- "exp"
 }

obj <- function(pp,all=FALSE)
{ # {{{

lp <- c(X %*% pp+offset)
if (model[1]=="exp") {
    p <- exp(lp) 
    D2logl <- c(weights*p)*X2 
} else if (model[1]=="lin") {
    p <- lp
    D2logl <- c(weights)*X2
} else if (model[1]=="logit") {
   p <- expit(lp)
   D2logl <- c(weights*p/(1+exp(lp)))*X2
} 

ploglik <- sum(weights*(Y-p)^2)
Dlogl <- weights*X*c(Y-p)
D2log <- apply(D2logl,2,sum)
gradient <- apply(Dlogl,2,sum)+augmentation
hessian <- matrix(D2log,length(pp),length(pp))

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      beta.iid <-  apply(beta.iid,2,sumstrata,id,max(id)+1)
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
	 id=id,Dlogl=Dlogl,
	 iid=beta.iid,robvar=robvar,var=robvar,se.robust=diag(robvar)^.5)
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
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }

  if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X)
  val <- c(val,list(time=time,formula=formula,formC=formC,
    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
    model.frame=m,n=length(exit),nevent=nevent,ncluster=nid,weights=weights))

  
# {{{ computation of ate, att, atc and their influence functions
if (length(all.vars(treat.model))==0) 
{ respformula  <-  as.formula( paste(all.vars(formula)[3],"~."))
  treat.model  <- update.formula(treat.model,respformula)
###  cat("Using propensity score model :\n "); print(treat.model)
}

## dropping cluster here 
treat.model <- drop.specials(treat.model,"cluster")
treat <- glm(treat.model,data,family="binomial")
Xtreat <- model.matrix(treat$formula,data)
ytreat <- treat$y

lpa <- treat$linear.predictors 
pal <- expit(lpa)
iidalpha <- iid(treat,id=id)

### treatment is rhs of treat.model 
treat.name <-  all.vars(treat.model)[1]
dat0 <- data
dat0[,treat.name] <- 0 ## treat.contrast[1]
dat1 <- data
dat1[,treat.name] <- 1 ## treat.contrast[2]

formulanc <- drop.specials(formula,"cluster")

X1 <- model.matrix(formulanc[-2],dat1)
X0 <- model.matrix(formulanc[-2],dat0)

p11lp <- X1 %*% val$coef+offset
p10lp <- X0 %*% val$coef+offset
p1lp <-   X %*% val$coef+offset

if (model[1]=="exp") {
   p1 <- exp(p1lp); p10 <- exp(p10lp); p11 <- exp(p11lp);
} else if (model[1]=="lin") {
   p1 <- p1lp; p10 <- p10lp; p11 <- p11lp; 
} else if (model[1]=="logit") {
    p1 <- expit(p1lp); p10 <- expit(p10lp); p11 <- expit(p11lp);
} else stop("link functions must be logit,exp,lin\n") 

Yorig <- Y
Y <- c(Y*obs/cens.weights)
risk1 <- ytreat*(Y-p11)/pal+p11
risk0 <- (1-ytreat)*(Y-p10)/(1-pal)+p10
riskG1 <- p11
riskG0 <- p10

ntreat <- sum(ytreat)
###att <- ytreat*Y-(pal*(1-ytreat)*Y + (ytreat - pal)* p10)/(1-pal)
###atc <- ((1-pal)*ytreat*Y - ((ytreat - pal)* p11)/pal)-(1-ytreat)*Y

if (model[1]=="logit") {
	Dp1 <- X * c(p1/(1+exp(p1lp)))
	Dp11 <- X1 * c(p11/(1+exp(p11lp)))
	Dp10 <- X0 * c(p10/(1+exp(p10lp)))
} else if (model[1]=="lin") {
      Dp1 <- X; Dp11 <- X1; Dp10 <- X0;
} else {
      Dp1 <- c(p1)*X; Dp11 <- c(p11)*X1; Dp10 <- c(p10)*X0;
}

Dpai <-  - Xtreat * exp(-lpa)
D1mpai <-   Xtreat * exp(lpa)

DaPsi1 <-  apply( Dpai * ytreat * c( Y - p1),2,sum)
DaPsi0 <-  apply( D1mpai * (1-ytreat) * c( Y - p1),2,sum)
DePsi1 <-  apply( Dp11 * ( 1- ytreat/pal),2,sum)
DePsi0 <-  apply( Dp10 * ( 1- (1-ytreat)/(1-pal)),2,sum)

###DePsiatt <- - apply( Dp10* (ytreat-pal)/(1-pal),2,sum)
###DaPsiatt <- apply(c((1-ytreat)*(Y-p10))*D1mpai,2,sum)
###
###DePsiatc <- -apply( Dp11* (ytreat-pal)/pal,2,sum)
###DaPsiatc <- apply(c(ytreat*(Y-p11))*Dpai,2,sum)

DriskG1 <- apply(Dp11,2,sum)
DriskG0 <- apply(Dp10,2,sum)
DdifriskG <- DriskG1-DriskG0

 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    Xtreat <- Xtreat[ord,,drop=FALSE]
    ytreat <- ytreat[ord]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    cens.weights <- cens.weights[ord]
    lp <- c(X %*% val$coef+offset)
    if (model[1]=="exp") { p <- exp(lp) } else if (model[1]=="lin") { p <- lp } else if (model[1]=="logit") {p <- expit(lp) } 
    Yglm <- weights*(Yout[ord]-p)
    Y <- Yout[ord]

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Yglm,2,revcumsumstrata,xx$strata,xx$nstrata)
    h10  <-  apply(cbind((ytreat/pal)*Y,I(ytreat==0)/(1-pal))*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
###  hattc  <-  apply(cbind(ytreat-pal*(1-ytreat)/(1-pal),-(1-ytreat)+(1-pal)*ytreat/pal)*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply((exit<=time)*h*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<=time)*h[xx$jumps+1,] /c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    IhdLamh10 <- apply(h10*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)
    U10 <- matrix(0,nrow(xx$X),2)
    U10[xx$jumps+1,] <- h10[xx$jumps+1,] /c(resC$S0)
    MGt10 <- (U10[,drop=FALSE]-IhdLamh10)*c(xx$weights)

###    IhdLamhattc <- apply(hattc*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)
###    Uattc <- matrix(0,nrow(xx$X),2)
###    Uattc[xx$jumps+1,] <- hattc[xx$jumps+1,]/c(resC$S0)
###    MGtattc <- (Uattc[,drop=FALSE]-IhdLamhattc)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    mid <- max(id)+1
    MGCiid <- apply(MGt,2,sumstrata,xx$id,mid)
    MGCiid10 <- apply(MGt10,2,sumstrata,xx$id,mid)
###    MGCiidattc <- apply(MGtattc,2,sumstrata,xx$id,mid)
  }  else { MGCiidattc <- MGCiid <- 0; MGCiid10 <- 0 }
## }}}

val$MGciid <- MGCiid
val$MGciid10 <- MGCiid10
val$id <- id
val$call.id <- call.id
val$name.id <- name.id
val$nid <- nid
val$iid.naive <- val$iid 
if (se) val$iid  <- val$iid+(MGCiid %*% val$ihessian)
val$naive.var <- val$var
robvar <- crossprod(val$iid)
val$var <-  val$robvar <- robvar
val$se.robust <- diag(robvar)^.5
val$se.coef <- diag(val$var)^.5
val$cause <- cause
val$cens.code <- cens.code 

### estimates risk, att, atc
val$riskDR <- c(mean(risk1),mean(risk0))
val$riskG<- c(mean(riskG1),mean(riskG0))

###val$att <- sum(att)/ntreat
###val$atc <- sum(atc)/(n-ntreat)
###val$attc <- c(val$att,val$atc)
###names(val$attc) <- c("ATT","ATC")

names(val$riskDR) <- paste("treat",1:0,sep="-")
names(val$riskG) <- paste("treat",1:0,sep="-")

## iid's of marginal risk estimates 
iidbase1 <- c(risk1-val$riskDR[1])
iidcif1 <- c(c(DePsi1) %*% t(val$iid))
iidpal1 <- c(c(DaPsi1) %*% t(iidalpha))
if (se)  {
iidGc1 <- MGCiid10[,1]; iidGc0 <- MGCiid10[,2]
###iidGatt <-  MGCiidattc[,1]; iidGatc <-  MGCiidattc[,2]
}  else { iidGc1 <- iidGatt  <- iidGatc  <- iidGc0  <- 0 } 

iidbase0 <- c(risk0-val$riskDR[2])
iidcif0 <- c(c(DePsi0) %*% t(val$iid))
iidpal0 <- c(c(DaPsi0) %*% t(iidalpha))

iidrisk1 <- iidbase1+iidcif1+iidpal1+iidGc1
iidrisk0 <- iidbase0+iidcif0+iidpal0+iidGc0
difriskiid <- (iidrisk1-iidrisk0)/n

iidrisk <- cbind(iidrisk1,iidrisk0)/n

###iidatt <- att-ytreat*val$att
###iidatc <- atc-(1-ytreat)*val$atc
###iidcifatt <- c(c(DePsiatt) %*% t(val$iid))
###iidpalatt <- c(c(DaPsiatt) %*% t(iidalpha))  
###iidcifatc <- c(c(DePsiatc) %*% t(val$iid))
###iidpalatc <- c(c(DaPsiatc) %*% t(iidalpha)) 
###iidatt <- iidatt+iidcifatt+iidpalatt+iidGatt
###iidatc <- iidatc+iidcifatc+iidpalatc+iidGatc

# }}}


val$call.id <- call.id
val$name.id  <- name.id
val$id <- id
val$nid  <- nid

# {{{ output variances and se for ate; cluster correction
val$id <- idclust
val$nid  <- nclust

## make iid form before combining into clusters 
riskG.iid <- cbind(c(p11-val$riskG[1]) + c(DriskG1 %*% t(val$iid)),
	          c(p10-val$riskG[2]) + c(DriskG0 %*% t(val$iid)))/n

## outcome model from IPCW logistic regression 
val$iid <- apply(val$iid,2,sumstrata,idclust,nclust)
robvar <- crossprod(val$iid)
val$var <-  val$robvar <- robvar
val$se.robust <- diag(robvar)^.5
val$se.coef <- diag(val$var)^.5
val$MGciid <- NULL
val$iid.naive <- NULL

## DR estimates 
iidrisk <- apply(iidrisk,2,sumstrata,idclust,nclust)
val$var.riskDR <- crossprod(iidrisk); 
val$se.riskDR <- diag(val$var.riskDR)^.5
val$riskDR.iid <- iidrisk

###nlevs <- 2
###pdiff <- function(x) lava::contr(lapply(seq(x-1), function(z) seq(z,x)))
###contrast <- -pdiff(length(nlevs))
###nncont <- c()
###for (k in seq_along(nlevs[-length(nlevs)])) nncont <-c(nncont, paste("treat:",nlevs[-seq(k)],"-",nlevs[k],sep="")) 
###rownames(contrast) <- nncont
contrast <- rbind(c(1,-1))

mm <- estimate(coef=val$riskDR,vcov=val$var.riskDR,contrast=rbind(c(1,-1)))
val$difriskDR <- mm$coef 
names(val$difriskDR) <-  rownames(contrast) 
val$var.difriskDR <- mm$vcov 
val$se.difriskDR <- diag(val$var.difriskDR)^.5

riskG.iid <- apply(riskG.iid,2,sumstrata,idclust,nclust)
###
val$riskG.iid <- riskG.iid
val$var.riskG <- crossprod(riskG.iid)
val$se.riskG <- diag(val$var.riskG)^.5
###
mm <- estimate(coef=val$riskG,vcov=val$var.riskG,contrast=contrast)
val$difriskG <- mm$coef 
names(val$difriskG) <-  rownames(contrast) 
val$var.difriskG <- mm$vcov 
val$se.difriskG <- diag(val$var.difriskG)^.5

### DR-estimator ; G -estimator sort after namid; also outcome model
if (!is.null(call.id)) {
	iidrisk <- namesortme(iidrisk,name.id)
	riskG.iid <-  namesortme(riskG.iid,name.id)
	val$iid <-  namesortme(val$iid,name.id)
}
## }}}

  val$model <- model[1]
  val$outcome <- outcome[1]
  val$call.id <- call.id
  class(val) <- "binreg"
  return(val)
}# }}}


##' @export
kumarsim <- function (n,rho1=0.71,rho2=0.40,rate = c(6.11,24.2),
		      beta=c(-0.67,0.59,-0.55,0.25,0.68,0.18,0.45,0.31),
		      labels= c("gp","dnr","preauto","ttt24(24,300]"),
		      depcens=0,type = c("logistic", "cloglog"),restrict=0 )
{# {{{
    p = length(beta)/2
    tt <- seq(0, 150, by = 1)
    if (length(rate) == 1)
    rate <- rep(rate, 2)
    Lam1 <- rho1 * (1 - exp(-tt/rate[1]))
    Lam2 <- rho2 * (1 - exp(-tt/rate[2]))
    
    ## fully saturated model for kumar covariates 
    Zdist <- c(0.21064815,0.02083333,0.05555556,0.01504630,
	       0.13888889,0.15393519, 0.02662037, 0.04398148, 
	       0.04745370, 0.02430556, 0.02199074, 0.03935185,
               0.02430556, 0.05555556, 0.01273148, 0.10879630)
    Zs <- expand.grid(gp=c(0,1),dnr=c(0,1),preauto=c(0,1),ttt24=c(0,1))
    Zs <- dsort(Zs,~ttt24+gp+dnr+preauto)

    samn <- sample(1:16,n,replace=TRUE,prob=c(Zdist))
    Z <- Zs[samn,]

    colnames(Z) <- labels
    cif1 <- setup.cif(cbind(tt, Lam1), beta[1:4], Znames = colnames(Z),
        type = type[1])
    cif2 <- setup.cif(cbind(tt, Lam2), beta[5:8], Znames = colnames(Z),
        type = type[1])
    if (restrict==0) 
    data <- sim.cifs(list(cif1, cif2), n, Z = Z)
    else {
    ## keep model 2 on logistic form
    data <- sim.cifsRestrict(list(cif2, cif1), n, Z = Z)
    data$status21 <- data$status
    data$status21[data$status==1] <- 2
    data$status21[data$status==2] <- 1
    data$status <- data$status21
    }

    ## kumar censoring, cox model 
    c0 <- list()     
    c0$cumhaz <- cbind(c(0,20,60,90,160),
		       c(0, 0.07797931, 0.28512764, 0.76116180, 1.95720759))
    c0$coef <- c(1.8503113,-0.6976226,0.5828763,-0.2000003)

    if (depcens == 0)
        censor  <- rchaz(c0$cumhaz,n=n)
    else {
        rrc <- exp(as.matrix(Z) %*% c0$coef)
        censor  <- rchaz(c0$cumhaz,rrc)
    }
    status = data$status * (data$time <= censor$time)
    time = pmin(data$time, censor$time)
    data <- data.frame(time = time, status = status)
    return(cbind(data, Z))
}# }}}

##' @export
kumarsimRCT <- function (n,rho1=0.71,rho2=0.40,rate = c(6.11,24.2),
		      beta=c(-0.67,0.59,-0.55,0.25,0.68,0.18,0.45,0.31),
		      labels= c("gp","dnr","preauto","ttt24(24,300]"),
		      nocens=0,addcens=1,rct=1,
		      type = c("logistic", "cloglog"),restrict=1,
		      censpar=c(1,1,1,1), F1par=c(1,1,1,1), F2par=c(1,1,1,1),
		      treatmodel=c(-0.18,-0.16,0.06,0.24) 
)
{# {{{
    p = length(beta)/2
    tt <- seq(0, 150, by = 1)
    if (length(rate) == 1)
    rate <- rep(rate, 2)
    Lam1 <- rho1 * (1 - exp(-tt/rate[1]))
    Lam2 <- rho2 * (1 - exp(-tt/rate[2]))
    
    ## fully saturated model for kumar covariates 
    Zdist <- c(0.21064815,0.02083333,0.05555556,0.01504630,
	       0.13888889,0.15393519, 0.02662037, 0.04398148, 
	       0.04745370, 0.02430556, 0.02199074, 0.03935185,
               0.02430556, 0.05555556, 0.01273148, 0.10879630)
    Zs <- expand.grid(gp=c(0,1),dnr=c(0,1),preauto=c(0,1),ttt24=c(0,1))
    Zs <- dsort(Zs,~ttt24+gp+dnr+preauto)

    samn <- sample(1:16,n,replace=TRUE,prob=c(Zdist))
    Z <- Zs[samn,]
    ## randomized gp instead
    if (rct==1) Z[,1] <- rbinom(n,1,0.5)
    ## randomized gp given other covariates 
    if (rct==2) {
	    lp <- as.matrix(cbind(1,Z[,-1])) %*% treatmodel 
	    Z[,1] <- rbinom(n,1, expit(lp))
    }
    colnames(Z) <- labels
    cif1 <- setup.cif(cbind(tt, Lam1), F1par*beta[1:4], Znames = colnames(Z), type = type[1])
    cif2 <- setup.cif(cbind(tt, Lam2), F2par*beta[5:8], Znames = colnames(Z), type = type[1])
    if (restrict==0) 
    data <- sim.cifs(list(cif1, cif2), n, Z = Z)
    else {
    ## keep model 2 on logistic form
    data <- sim.cifsRestrict(list(cif2, cif1), n, Z = Z)
    data$status21 <- data$status
    data$status21[data$status==1] <- 2
    data$status21[data$status==2] <- 1
    data$status <- data$status21
    }

    if (nocens==0) {
    ## kumar censoring, cox model 
    c0 <- list()     
    c0$cumhaz <- cbind(c(0,20,60,90,160),
		       c(0, 0.07797931, 0.28512764, 0.76116180, 1.95720759))
    c0$coef <-censpar* c(1.8503113,-0.6976226,0.5828763,-0.2000003)

	if (addcens)  {
	  ## draw from cox gp model  
         rrc <- exp(as.matrix(Z) %*% c(c(1,0,0,0)*c0$coef))
         censorgp  <- rchaz(c0$cumhaz,rrc)
	 ## draw from cox other covariates   
         rrc <- exp(as.matrix(Z) %*% c(c(0,1,1,1) * c0$coef))
         censoroc  <- rchaz(c0$cumhaz,rrc)
	 ## miniumum of the two censoring times 
	 censor <- censorgp
	 censor$time <- pmin(censoroc$time,censorgp$time)
	 censor$status <- pmax(censoroc$status,censorgp$status)
	} else {
        rrc <- exp(as.matrix(Z) %*% c0$coef)
        censor  <- rchaz(c0$cumhaz,rrc)
	}
    status = data$status * (data$time <= censor$time)
    time = pmin(data$time, censor$time)
    data <- data.frame(time = time, status = status)
    }
    return(cbind(data, Z))
}# }}}

##' @export
logitATE <- function(formula,data,binreg=TRUE,...)
{# {{{
   ## use IPCW machine in no-censoring case
    cl <- match.call()
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster", "offset")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (inherits(Y,"Event")) {
    if (binreg) out <- binregATE(formula,data,...) else out <- logitIPCWATE(formula,data,...)
    } else {
      response <- all.vars(formula)[1]
      Ydirect <-  as.numeric(data[,response]) 
      data$time <- 2
      data$event <- 1
      time <- 2
      Survform <-  update.formula(formula,Event(time,event)~.)
      n <- nrow(data)
      if (binreg)
      out <- binregATE(Survform,data,se=0,cens.weights=rep(1,n),time=time,Ydirect=Ydirect,...)
      else
      out <- logitIPCWATE(Survform,data,se=0,cens.weights=rep(1,n),time=time,Ydirect=Ydirect,...)
    }

   return(out)
}# }}}

##' @export
normalATE <- function(formula,data,binreg=TRUE,model="lin",...)
{# {{{
   if (binreg) 
   out <- logitATE(formula,data,binreg=binreg,model=model,outcome="rmst",...)
   else 
   out <- logitATE(formula,data,binreg=binreg,...)
   return(out)
}# }}}

##' Augmentation for Binomial regression based on stratified NPMLE Cif (Aalen-Johansen) 
##'
##' Computes  the augmentation term for each individual as well as the sum
##' \deqn{
##' A = \int_0^t H(u,X) \frac{1}{S^*(u,s)} \frac{1}{G_c(u)} dM_c(u)
##' }
##' with 
##' \deqn{
##' H(u,X) = F_1^*(t,s) - F_1^*(u,s)
##' }
##' using a KM for \deqn{G_c(t)} and a working model for cumulative baseline
##' related to \deqn{F_1^*(t,s)} and \deqn{s} is strata, \deqn{S^*(t,s) = 1 - F_1^*(t,s) - F_2^*(t,s)}. 
##'
##' Standard errors computed under assumption of correct \deqn{G_c(s)} model.
##' 
##' @param formula formula with 'Event', strata model for CIF given by strata, and strataC specifies censoring strata
##' @param data data frame
##' @param offset offsets for cox model
##' @param data data frame
##' @param cause of interest 
##' @param cens.code code of censoring 
##' @param km to use Kaplan-Meier
##' @param time of interest 
##' @param weights weights for estimating equations
##' @param offset offsets for logistic regression
##' @param ... Additional arguments to binreg function.
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data(bmt)
##' dcut(bmt,breaks=2) <- ~age 
##' out1<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
##'			  strata(platelet,agecat.2),data=bmt,cause=1,time=40)
##' summary(out1)
##' 
##' out2<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
##'     strata(platelet,agecat.2)+strataC(platelet),data=bmt,cause=1,time=40)
##' summary(out2)
##' @export
BinAugmentCifstrata <- function(formula,data=data,cause=1,cens.code=0,km=TRUE,time=NULL,weights=NULL,offset=NULL,...)
{# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset","strataC")
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
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else strata.name <- NULL
  if (!is.null(stratapos <- attributes(Terms)$specials$strataC)) {
    ts <- survival::untangle.specials(Terms, "strataC")
    Terms  <- Terms[-ts$terms]
    strataC <- as.numeric(m[[ts$vars]])-1
    strataC.name <- ts$vars
  }  else { strataC <- NULL; strataC.name <- NULL}

  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms  <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  id.orig <- id; 
  if (!is.null(id)) {
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else { id <- as.integer(seq_along(exit))-1;  nid <- length(exit); }


 p <- ncol(X)
 beta <- NULL
  if (is.null(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(exit)))

  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }

if (is.null(strataC)) { strataC <- rep(0,length(exit)); nstrataC <- 1; strataC.level <- NULL; } else {
	  strataC.level <- levels(strataC)
	  ustrataC <- sort(unique(strataC))
	  nstrataC <- length(ustrataC)
	  strataC.values <- ustrataC
      if (is.numeric(strataC)) strataC <-  fast.approx(ustrataC,strataC)-1 else  {
      strataC <- as.integer(factor(strataC,labels=seq(nstrataC)))-1
    }
  }

  cens.strata <- strataC
  cens.nstrata <- nstrataC 

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
  strata.call <- strata
  Z <- NULL
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z

  ## possible casewights to use for bootstrapping and other things
  case.weights <- NULL 
  if (is.null(case.weights)) case.weights <- rep(1,length(exit)) 

  trunc <- (!is.null(entry))
  if (!trunc) entry <- rep(0,length(exit))
  call.id <- id

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1; 
   ## orginal id coding into integers 
   id.orig <- id+1; 

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 

  statusC <- (status%in%cens.code) 
  statusE <- (status%in%cause) & (exit<= time) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 

  Zcall <- cbind(status,strata)
  dd <- .Call("FastCoxPrepStrata",entry,exit,statusC,X,id, 
	     trunc,strataC,weights,offset,Zcall,case.weights,PACKAGE="mets")

  jumps <- dd$jumps+1
  xxstrataC <- c(dd$strata)
  xxstatus  <- dd$Z[,1]
  xxstrata  <- dd$Z[,2]
  jumpsD <- which(xxstatus!=cens.code)
  jumps1 <- which(xxstatus==cause)
  rr <- c(dd$sign*exp(dd$offset))
  S0 = c(revcumsumstrata(rr,xxstrata,nstrata))
  ## S0 after strataC
  S00C = c(revcumsumstrata(rr,xxstrataC,nstrataC))

  S0C <- rep(0,length(dd$strata))

  ## censoring MG, strataC
  stratJumps <- dd$strata[jumps]
  S00i <- rep(0,length(dd$strata))
  S00i[jumps] <-  1/S00C[jumps]

  ## cif calculation, uses strata {{{
  S0Di <- S02i <- S01i <- rep(0,length(dd$strata))
  S01i[jumps1] <-  1/S0[jumps1]
  S0Di[jumpsD] <-  1/S0[jumpsD]

  ## strata-Cif G_T(t)
  if (!km) { 
     cumhazD <- cumsumstratasum(S0Di,xxstrata,nstrata)
     Stm <- exp(-cumhazD$lagsum)
     St <- exp(-cumhazD$sum)
  } else { 
     StA <- cumsumstratasum(log(1-S0Di),xxstrata,nstrata)
     Stm <- exp(StA$lagsum)
     St <- exp(StA$sum)
  }

  ## G_c(t-) 
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S00i,xxstrataC,nstrataC)$lagsum)
    Gc      <- exp(-cumhazD)
  } else Gc <- c(exp(cumsumstratasum(log(1-S00i),xxstrataC,nstrataC)$lagsum))

# }}}
  
 btime <- c(1*(dd$time<=time))
 cif1 <- cumsumstrata(Stm*S01i*btime,xxstrata,nstrata)
 ### final F_1^s(time) 
 cif1time <- cif1[tailstrata(xxstrata,nstrata)]
 ciftt <- cif1time[xxstrata+1]
 ft <-  (ciftt-cif1)/(Gc*St)
 ft[Gc<0.00001] <- 0
 ft[St<0.00001] <- 0

 S1ft <- c(revcumsumstrata(rr*ft,xxstrata,nstrata))
 Eft <- S1ft/S0

 Z <- dd$X
 U1 <- matrix(0,nrow(Z),1)
 U1[jumps,] <- ft[jumps]*btime[jumps]
 E1dLam0 <- cumsum2strata(ft,S00i*btime,xxstrata,nstrata,xxstrataC,nstrataC,cif1time)$res

 ### Martingale  as a function of time and for all subjects to handle strata 
 mm <- cbind(U1,E1dLam0)
 MGt <- Z*c(U1-E1dLam0)*rr*c(dd$weights)
 MGiid <- apply(MGt,2,sumstrata,dd$id,max(id)+1)
 augment <- apply(MGt,2,sum)

 ## drop strata's from formula and run with augmention term
 formulans <- drop.strata(formula)
 ## censoring weights not used 
 if (nstrataC==1) cens.model <- ~+1 else cens.model <- ~strata(strataCC)
 data$strataCC <- strataC

 bra <- binreg(formulans,data=data,cause=cause,augmentation=augment,time=time,
	       cens.model=cens.model,...)

 ## adjust SE and var based on augmentation term
 ## only report SE based on iid 
 bra$var.orig <- bra$var
 bra$augment <- augment
 ## with correct augmentation term, things cancel out 
 bra$iid <- bra$iid.naive + MGiid %*%  bra$ihessian
 bra$var <- crossprod(bra$iid)
 bra$se.coef <-  diag(bra$var)^.5
 bra$robvar <- bra$var
 bra$se.robust <-bra$se.coef
 bra$MGciid <- MGiid
 bra$model <- "logit"

 allAugment <- list(MGiid=MGiid,augment=augment,id=id,id.orig=id.orig,
	       cif=cif1,St=St,Gc=Gc,strata=xxstrata,strataC=xxstrataC,time=dd$time)
 bra$allAugment <- allAugment

 return(bra)
}# }}})

