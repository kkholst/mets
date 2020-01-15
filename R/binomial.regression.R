##' Binomial Regression for censored competing risks data 
##'
##' Simple version of comp.risk function of timereg for just one time-point thus fitting the model 
##' \deqn{P(T \leq t, \epsilon=1 | X ) = expit( X^T beta) }
##'
##' Based on binomial regresion IPCW response estimating equation: 
##' \deqn{ X ( \Delta I(T \leq t, \epsilon=1 )/G_c(T_i-) - expit( X^T beta)) = 0 }
##' for IPCW adjusted responses. 
##'
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param cause cause of interest
##' @param time  time of interest 
##' @param beta starting values 
##' @param offset offsets for partial likelihood 
##' @param weights for score equations 
##' @param cens.weights censoring weights 
##' @param cens.model stratified cox model 
##' @param se to compute se's  based on IPCW 
##' @param kaplan.meier uses Kaplan-Meier for baseline than standard Cox 
##' @param cens.code gives censoring code
##' @param no.opt to not optimize 
##' @param method for optimization 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##'
##' data(bmt)
##' # logistic regresion with IPCW binomial regression 
##' out <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50)
##' summary(out)
##' predict(out,data.frame(tcell=c(0,1),platelet=c(1,1)),se=TRUE)
##'
##' outs <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50,cens.model=~strata(tcell,platelet))
##' summary(outs)
##' 
##' @export
binreg <- function(formula,data,cause=1,time=NULL,beta=NULL,
		   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,
		   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",...)
{# {{{

  cl <- match.call()# {{{
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (class(Y)!="Event") stop("Expected a 'Event'-object")
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
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(exit))-1; 
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
      cens.weights <- predict(resC,data,times=exit,tminus=TRUE,individual.time=TRUE,se=FALSE,km=kmt)$surv
      ## strata from original data 
      cens.strata <- resC$strata[resC$ord]
      cens.nstrata <- resC$nstrata
  }
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  if (is.null(beta)) beta <- rep(0,ncol(X))
  p <- ncol(X)

  X <-  as.matrix(X)
  X2  <- .Call("vecMatMat",X,X)$vXZ
###mm <-  .Call("CubeVec",D2logl,Dlogl)

obj <- function(pp,all=FALSE)
{ # {{{

lp <- c(X %*% pp)
p <- expit(lp)
###
Y <- c((status==cause)*(exit<=time)/cens.weights)
ploglik <- sum(weights*(Y-p)^2)

Dlogl <- weights*X*c(Y-p)
D2logl <- c(weights*p/(1+exp(lp)))*X2
D2log <- apply(D2logl,2,sum)
###
gradient <- apply(Dlogl,2,sum)
hessian <- matrix(D2log,length(pp),length(pp))

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      beta.iid <-  apply(beta.iid,2,sumstrata,id,max(id)+1)
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
	 id=id,Dlogl=Dlogl,
	 iid=beta.iid,robvar=robvar,var=robvar,
         se=diag(robvar)^.5,se.robust=diag(robvar)^.5)
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
	  val <- c(val,list(time=time,formula=formula,
	    exit=exit,
	    cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata,
			    model.frame=m))

 if (se) {## {{{ censoring adjustment of variance 
	ord <- resC$cox.prep$ord+1
        X <-  X[ord,,drop=FALSE]
        status <- status[ord]
	exit <- exit[ord]
        cens.weights <- cens.weights[ord]
	lp <- c(X %*% val$coef)
	p <- expit(lp)
        Y <- c((status==cause)*(exit<=time)/cens.weights)

       Dlogl <- weights*X*c(Y-p)
       hessian <- val$hessian 
       xx <- resC$cox.prep
       S0i2 <- S0i <- rep(0,length(xx$strata))
       S0i[xx$jumps+1]  <- 1/resC$S0
       S0i2[xx$jumps+1] <- 1/resC$S0^2
       U <- E <- matrix(0,nrow(xx$X),ncol(X))
###       Ys <- revcumsumstrata(xx$sign,xx$strata,xx$nstrata)
	
       ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
       ## to make \int h(s)/Ys  dM_i^C(s) 
       h  <-  apply(X*Y*(exit<=time),2,revcumsumstrata,xx$strata,xx$nstrata)
       h2  <- .Call("vecMatMat",h,h)$vXZ
       U[xx$jumps+1,] <- h[xx$jumps+1]/resC$S0
       hdLam0 <- apply(h*S0i^2,2,cumsumstrata,xx$strata,xx$nstrata)
       ### Cens-Martingale as a function of time and for all subjects to handle strata 
       MGt <- (U[,drop=FALSE]-hdLam0)*c(xx$weights)
       ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s)
       ### estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
       h2dLam0 <- apply(h2*S0i2,2,sum)
       varadjC <- matrix(h2dLam0,length(val$coef),length(val$coef))
       MGCiid <- apply(MGt,2,sumstrata,id,max(id)+1)
       val$varadjC <- val$ihessian %*% varadjC %*% val$ihessian
 
       val$nc.iid <- val$iid 
       beta.iid <- val$iid+(MGCiid %*% val$ihessian)
       val$iid  <- beta.iid
       val$naive.var <- val$var
       val$var <- val$var - val$varadjC
       robvar <- crossprod(beta.iid)
       val$robvar <- robvar
       val$se.robust <- diag(robvar)^.5
       val$se <- diag(val$var)^.5
  } ## }}}

  class(val) <- "binreg"
  return(val)
}# }}}

##' @export
iid.binreg  <- function(x,...) {# {{{
  x$iid
}# }}}

##' @export
print.binreg  <- function(x,...) {# {{{
  print(summary(x),...)
}# }}}

##' @export
summary.binreg <- function(object,or=TRUE,...) {# {{{
if (or)  {
cat("OR estimates \n"); 
estimate(coef=object$coef,vcov=object$var,f=function(p) exp(p))
} else  {
	cat("log-OR estimates \n"); 
estimate(coef=object$coef,vcov=object$var)
}
}# }}}

##' @export
vcov.binreg <- function(object,...) {# {{{
	return(object$var)
}# }}}

##' @export
predict.binreg <- function(object,newdata,se=TRUE,...)
{# {{{

  xlev <- lapply(object$model.frame,levels)
  ff <- unlist(lapply(object$model.frame,is.factor))
  upf <- update(object$formula,~.)
  tt <- terms(upf)
  tt <- delete.response(tt)
  Z <- model.matrix(tt,data=newdata,xlev=xlev)
  Z <- as.matrix(Z)
  expit  <- function(z) 1/(1+exp(-z)) ## expit
  lp <- c(Z %*% object$coef)
  p <- expit(lp)
  preds <- p

  if (se) {
  preds <- c()
  for (i in 1:length(lp)) {
     if (is.null(object$var)) covv <- vcov(object)  else covv <- object$var
     Dp <- Z[i,]*exp(-lp[i])*p[i]^2
     se <- (Dp %*% covv %*% Dp)^.5
     cmat <- data.frame(pred=p[i],se=se,lower=p[i]-1.96*se,upper=p[i]+1.96*se)
     names(cmat)[1:4] <- c("pred","se","lower","upper")
     preds <- rbind(preds,cmat)
  } 
  }

return(preds)
} # }}}

###predict.binreg <- function(x,newdata,se=TRUE,...)
###{# {{{
###  Z <- as.matrix(model.matrix(x$formula,newdata))
###  expit  <- function(z) 1/(1+exp(-z)) ## expit
###  lp <- c(Z %*% x$coef)
###  p <- expit(lp)
###  preds <- p
###
###  if (se) {
###    Ft <- function(p,lpi=1)
###    {# {{{
###	 p <- expit(lpi)
###	 return(p)
###    }# }}}
###
###  preds <- c()
###  for (i in 1:length(lp)) {
###     if (is.null(x$var)) covv <- vcov(x)  else covv <- x$var
###     eud <- estimate(coef=x$coef,vcov=covv,f=function(p) Ft(p,lpi=lp[i]))
###     cmat <- data.frame(eud$coefmat)
###     names(cmat)[1:4] <- c("pred","se","lower","upper")
###     preds <- rbind(preds,cmat)
###  } 
###  }
###
###return(preds)
###} # }}}
###

