##' Discrete time to event interval censored data 
##'
##' We consider the cumulative odds model 
##' \deqn{
##'    P(T \leq t | x) =  \frac{G(t) \exp(x \beta) }{1 + G(t) exp( x \beta) }
##' }
##' or equivalently 
##' \deqn{
##'    logit(P(T \leq t | x)) = log(G(t)) + x \beta
##' }
##' and we can thus also compute the probability of surviving 
##' \deqn{
##'    P(T >t | x) =  \frac{1}{1 + G(t) exp( x \beta) }
##' }
##' 
##' The baseline \eqn{G(t)} is written as \eqn{cumsum(exp(\alpha))} and this is not the standard
##' parametrization that takes log of \eqn{G(t)} as the parameters. Note that the regression 
##' coefficients are describing the probability of dying before or at time t. 
##' 
##' Input are intervals given by ]t_l,t_r] where t_r can be infinity for right-censored intervals 
##' When truly discrete ]0,1] will be an observation at 1, and  ]j,j+1] will be an observation at j+1.
##' Can be used for fitting the usual ordinal regression model (with logit link) that in contrast, however, 
##' describes the probibility of surviving time t (thus leads to -beta).
##' 
##' Likelihood is maximized:
##' \deqn{
##'  \prod  P(T_i >t_{il} | x) - P(T_i> t_{ir}| x) 
##' }
##' 
##' @param formula  formula
##' @param data  data 
##' @param beta starting values 
##' @param no.opt optimization TRUE/FALSE 
##' @param method NR, nlm 
##' @param stderr to return only estimate 
##' @param weights weights following id for GLM 
##' @param offsets following id  for GLM
##' @param exp.link parametrize increments exp(alpha) > 0
##' @param increment using increments dG(t)=exp(alpha) as parameters
##' @param ... Additional arguments to lower level funtions lava::NR  optimizer or nlm
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data(ttpd) 
##' dtable(ttpd,~entry+time2)
##' 
##' out <- interval.logitsurv.discrete(Interval(entry,time2)~X1+X2+X3+X4,ttpd)
##' summary(out)
##' head(iid(out)) 
##' 
##' pred <- predictlogitSurvd(out,se=FALSE)
##' plotSurvd(pred)
##' 
##' ttpd <- dfactor(ttpd,fentry~entry)
##' out <- cumoddsreg(fentry~X1+X2+X3+X4,ttpd)
##' summary(out)
##' head(iid(out)) 
##' 
##' @aliases Interval dInterval simlogitSurvd predictlogitSurvd cumoddsreg simTTP predictSurvd plotSurvd 
##' @export
interval.logitsurv.discrete <- function (formula,data,beta=NULL,no.opt=FALSE,method="NR",
	   stderr=TRUE,weights=NULL,offsets=NULL,exp.link=1,increment=1,...)
{ ## {{{ 

  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")

 if (ncol(Y)==2) {
	time2 <- eventtime <- Y[,2]
	entrytime <- Y[,1]
	left <- 0
    } else {
	time2 <- eventtime <- Y[,2]
	status <- delta  <- Y[,3]
	entrytime <- Y[,1]
	left <- 1
	if (max(entrytime)==0) left <- 0
    }

  if (any(entrytime==time2)) stop("left==right not possible for discrete data")

  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    pos.cluster <- ts$terms
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  } else pos.cluster <- NULL
  if (!is.null(attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}

  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
  X <- X[,-intpos,drop=FALSE]
  if (ncol(X)>0) X.names <- colnames(X) else X.names <- NULL

###  if (!is.null(id)) {
###	  ids <- unique(id)
###	  nid <- length(ids)
###      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
###      id <- as.integer(factor(id,labels=seq(nid)))-1
###     }
###   } else { id <- as.integer(seq_along(time2))-1; nid <- length(time2) }
###   ## orginal id coding into integers 
###   id.orig <- id+1; 

  call.id <- id
  conid <- construct_id(id,nrow(X))
  id <- conid$id; nid <- conid$nid; name.id <- conid$name.id

  ## times 1 -> mutimes , 0 til start
  utimes <- sort(unique(c(time2,entrytime)))
  mutimes <- max(utimes[utimes<Inf])
  n <- length(time2)
  
  ## setting up designs for t_l and t_r
  tR <- tL <- matrix(0,n,mutimes)

  ## design for ]t_l,t_r], for t_l=0 row is  0
  if (increment==0) {
	  tL <- matdoubleindex(tL,1:n,entrytime,rep(1,n))
	  tR <- matdoubleindex(tL,1:n,time2,rep(1,n))
  } else {
     for (i in 1:mutimes) {
        tL[,i] <- (i <= entrytime) 
        tR[,i] <- (i <= time2) 
     }
  }

  ## computing X^2 
  if (ncol(X)>0) {
  X2  <- .Call("vecMatMat",X,X)$vXZ
  XtL  <- .Call("vecMatMat",X,tL)$vXZ
  XtR  <- .Call("vecMatMat",X,tR)$vXZ
  }
  tL2  <- .Call("vecMatMat",tL,tL)$vXZ
  tR2  <- .Call("vecMatMat",tR,tR)$vXZ

  ## weights/offets will follow id 
  if (is.null(weights))  weights <- rep(1,n); #  else wiid <- weights
  if (is.null(offsets))  offsets <- rep(0,n); # else offsets <- offsets

  if (is.null(beta)) {
     beta <- rep(0,ncol(X)+mutimes)
     Set <- 1-cumsum(table(time2))/n
     dHt <- log(diff(c(0,(1/Set-1))))
     beta[1:mutimes] <- dHt[1:mutimes]
  }


obj <- function(pp,all=FALSE)
{ # {{{

if (exp.link==1) theta <- exp(pp[1:mutimes]) else theta  <-  pp[1:mutimes]

if (ncol(X)>0) {
betal   <- pp[-c(1:mutimes)]
Zbeta <- c(X %*% betal+offsets)
} else {Zbeta <- offsets }
xltheta <-  c(tL %*% theta)
xrtheta <-  c(tR %*% theta)
EZbeta <- exp(Zbeta)
GEl <- xltheta * EZbeta
GEr <- xrtheta * EZbeta
Stl <- 1/(1+GEl) 
Str <- 1/(1+GEr) 

############################################
## likelihood
############################################
# {{{
p <- Stl-Str*(time2<Inf)
logp <- log(p)
logl <- weights*logp 
# }}}
############################################
## Derivative 
############################################
# {{{
if (ncol(X)>0) {
DbetaStl <- -X*GEl*Stl^2
DbetaStr <- -(time2<Inf)*X*GEr*Str^2
Dbetalogp <- (DbetaStl-DbetaStr)/p
} 
DgStl <- -tL*EZbeta*Stl^2
DgStr <- -(time2<Inf)*tR*EZbeta*Str^2
Dglogp <- (DgStl-DgStr)/p

if (ncol(X)>0) {
Dlogliid  <- cbind(Dglogp*weights,Dbetalogp*weights)
Dlogl  <- apply(Dlogliid,2,sum)
} else {
	Dlogliid  <- Dglogp*weights
        Dlogl  <- apply(Dlogliid,2,sum)
}
if (exp.link==1) {
Dlogliid[,1:mutimes] <- t( t(Dlogliid[,1:mutimes])*theta)
Dlogl[1:mutimes] <- Dlogl[1:mutimes]*theta
}
# }}}
############################################
## 2nd Derivative 
############################################
# {{{
if (ncol(X)>0) {
D2betaStl <- X2*(GEl^2-GEl)*Stl^3
D2betaStr <- X2*(time2<Inf)*(GEr^2-GEr)*Str^3
DbetagStl <- XtL*            (EZbeta*(GEl-1))*Stl^3
DbetagStr <- XtR*(time2<Inf)*(EZbeta*(GEr-1))*Str^3
}
D2gStl <-             2*tL2*EZbeta^2*Stl^3
D2gStr <- 2*tR2*(time2<Inf)*EZbeta^2*Str^3

DglpDglp       <- .Call("vecMatMat",Dglogp,Dglogp)$vXZ
D2glogp <-    (D2gStl-D2gStr)/p  - DglpDglp
D2g <-    apply(weights*D2glogp,2,sum)

## cross products of derivatives
if (ncol(X)>0) {
DbetalpDbetalp <- .Call("vecMatMat",Dbetalogp,Dbetalogp)$vXZ
DbetalpDglp    <- .Call("vecMatMat",Dbetalogp,Dglogp)$vXZ

D2betalogp <- (D2betaStl-D2betaStr)/p - DbetalpDbetalp
Dbetaglogp <- (DbetagStl-DbetagStr)/p  - DbetalpDglp
D2beta <- apply(weights*D2betalogp,2,sum)
Dbetag <- apply(weights*Dbetaglogp,2,sum)
}

D2log <- matrix(0,length(pp),length(pp))
D2log[1:mutimes,1:mutimes] <- D2g

if (ncol(X)>0) {
###D2log[1:mutimes,(mutimes+1):length(pp)] <- Dbetag
D2log[(mutimes+1):length(pp),1:mutimes] <- Dbetag
D2log[1:mutimes,(mutimes+1):length(pp)] <- t(D2log[(mutimes+1):length(pp),1:mutimes])
D2log[(mutimes+1):length(pp),(mutimes+1):length(pp)] <- D2beta
}

diagM <- function(vec)
{# {{{
	p <- length(vec)
	out <- matrix(0,p,p)
	out <- matdoubleindex(out,1:p,1:p,vec)
	return(out)
}
# }}}

if (exp.link==1) {
D2log[1:mutimes,1:mutimes] <- diagM(Dlogl[1:mutimes]) + D2log[1:mutimes,1:mutimes]*(theta %o% theta)  
if (ncol(X)>0) {
D2log[(mutimes+1):length(pp),1:mutimes] <- Dbetag*rep(theta,each=length(betal))
###D2log[(mutimes+1):length(pp),1:mutimes] <- Dbetag*rep(theta,length(betal))
D2log[1:mutimes,(mutimes+1):length(pp)] <- t(D2log[(mutimes+1):length(pp),1:mutimes])
}
}

# }}}
############################################

ploglik <- sum(logl)
gradient <- Dlogl 
hessian <- D2log

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogliid %*% ihess ## %*% t(Dlogl) 
      beta.iid <- apply(beta.iid,2,sumstrata,id,nid)
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
		  iid=beta.iid,robvar=robvar,var=robvar,ihessian=ihess,id=id,
		  se=diag(robvar)^.5,coef=pp,se.coef=diag(robvar)^.5)
      return(val)
  }  
 structure(-ploglik/nid,gradient=-gradient/nid,hessian=-hessian/nid)
}# }}}

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
      if (!stderr) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
  } else val <- c(list(coef=beta),obj(beta,all=TRUE))

  uu <- utimes[utimes<Inf]
  Xnames <- c(paste("time",uu[-1],sep=""),X.names)
  if (length(val$coef)==length(Xnames)) names(val$coef) <- Xnames

  if (is.matrix(val$iid)) 
	  if (length(name.id)==nrow(val$iid)) rownames(val$iid) <- name.id

  val <- c(list(increment=increment,exp.link=exp.link,ntimes=mutimes,utimes=utimes,
		name.id=name.id,call.id=call.id,nid=nid),val)

  class(val) <- c("cumoddsreg")
  return(val)
} ## }}} 

##' @export
IC.cumoddsreg <- function(x,...) { x$iid*NROW(x$iid) }

##' @export
cumoddsreg <- function (formula,data,...)
{ ## {{{ 
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  nt <- nlevels(Y)
  time2__ <- as.numeric(Y)
  entrytime__ <- time2__-1
  time2__[time2__==nt] <- Inf

  xf <- update.formula(formula,Interval(entrytime__,time2__)~.)
  data$entrytime__ <- entrytime__
  data$time2__ <- time2__
  out <- interval.logitsurv.discrete(xf,data,...)

 return(out)
} ## }}}

##' @export
Interval <- function (time, time2 , ...)
{# {{{
    out <- cbind(time, time2 )
    colnames(out) <- c("entry", "time2")
    class(out) <- "Interval"
    return(out)
}# }}}

##' @export
dInterval <- function (time, time2 ,cuts=NULL,cut.first=0,show=FALSE, ...) 
{# {{{
if (is.null(cuts)) cuts <-  sort(unique(time,time2))
if (min(cuts)> 0)  cuts <- c(cut.first,cuts)
###
lleft <- fast.approx(cuts,time,type="left")
rright <- fast.approx(cuts,time2,type="right")

out <- data.frame(time=time,time2=time2,left=lleft-1,right=rright-1)
attr(out,"cuts") <- cuts 
out$leftd <- cuts[lleft]
out$rightd <- cuts[rright]

if (max(cuts)==Inf) out$right[out$rightd==Inf] <- Inf

if (any(out$time<out$leftd)) warning("left not in [leftd,rightd]\n")
if (any(out$time2>out$rightd)) warning("right  not in [leftd,rightd]\n") 

if (show) {
	mtime <- mmtime <- max(cuts[cuts<Inf])
	n <- length(time)
	if (max(cuts)==Inf) mtime <- mmtime+1
	plot(0,0,xlim=c(0,mtime),ylim=c(0,n),type="n")
	abline(v=cuts,col=3)
	lines(c(min(cuts),mmtime),c(-1,-1),lwd=2)
	if (max(cuts)==Inf) lines(c(mtime,mmtime),c(-1,-1),lwd=2,col=2)
	for (i in 1:n)
	{
	lines(c(out$time[i],min(out$time2[i],mtime)),rep(i,2),col=1)
	lines(c(out$leftd[i],min(out$rightd[i],mtime)),rep(i+0.5,2),col=2)
	}
}

return(out)
}# }}}

##' @export
simlogitSurvd <- function(coef,Z,n=NULL,times=0:6)
{# {{{
  if (missing(Z)) Z <- NULL
  if (!is.null(Z)) n <- nrow(Z) 
  Z <- as.matrix(Z)

  g <- coef[times[-1]]
  beta <- coef[-(times[-1])]
  Gt <- cumsum(c(0,exp(g)))
  if (!is.null(Z)) EZbeta <-  c(exp(Z %*% beta)) else EZbeta <- rep(1,n)
  EZbeta <- rep(EZbeta,each=length(times))
  Gti <- rep(Gt,n)
  GE <- Gti*EZbeta
  Ft <- GE/(1+GE)

  ru <- runif(n)
  out <- data.frame(.Call("wherestrataR",ru,Ft,rep(0:(n-1),each=length(times)),n))
  status <- rep(1,n)
  names(out)[1] <- "time"
  status[out$time==(length(times)-1)] <- 0
  out$status <- status
  out <- cbind(out,Z)

return(out)
}# }}}

##' @export
predictlogitSurvd <- function(x,Z,n=NULL,times=NULL,se=FALSE,type="prob")
{# {{{

if (missing(Z)) Z <- NULL
  if (!is.null(Z))  {
	  n <- nrow(Z)
	  Z <- as.matrix(Z)
  } else { if (is.null(n)) n <- 1; }
  ccc <- x$coef

  ntimes <- x$ntimes
  g <- ccc[1:ntimes]
  beta <- ccc[-(1:ntimes)]
  if (is.null(times)) times <- 0:ntimes 

  Gt <- cumsum(c(0,exp(g)))
  if (!is.null(Z)) EZbeta <-  c(exp(Z %*% beta)) else EZbeta <- rep(1,n)

  if (!se) {# {{{
	  EZbeta <- rep(EZbeta,each=length(times))
	  id <- rep(1:n,each=length(times))
	  Gti <- rep(Gt,n)
	  GE <- Gti*EZbeta
	  Gt <- 1/(1+GE)
	  preds <- data.frame(pred=Gt,id=id,times=rep(times,n))
# }}}
  } else {# {{{

    Ft <- function(p,Zi=1)
    {# {{{
         g <- p[1:ntimes]
         beta <- p[-(1:ntimes)]
         Gt <- cumsum(c(0,exp(g)))
	 if (length(beta)>0) EZ <- exp(sum(Zi*beta)) else EZ <- 1
	 pred <- 1/(1+Gt*EZ)
	 return(pred)
    }# }}}

  preds <- c()
  for (i in 1:length(EZbeta)) {
     if (is.null(x$var)) covv <- vcov(x)  else covv <- x$var
     eud <- estimate(coef=x$coef,vcov=covv,f=function(p) Ft(p,Zi=Z[i,]))
     cmat <- data.frame(eud$coefmat)
     cmat$id <- i
     cmat$times <- times
     names(cmat)[1:4] <- c("pred","se","lower","upper")
     preds <- rbind(preds,cmat)
  } 

  }# }}}

return(preds)
} ## }}} 

##' @export
simTTP <- function(coef=NULL,n=100,Xglm=NULL,times=NULL)
{# {{{
	  
  Z <- Xglm  
  if (!is.null(Z)) n <- nrow(Z) 

  if (!is.null(Z)) data <- Z else data <- data.frame(id=1:n)

  if (!is.null(times)) {
     timesf <- data.frame(times=rep(times,n),id=rep(1:n,each=length(times)))
     data <- merge(data,timesf,by.x="id",by.y="id")
     mt <- model.matrix(~factor(times),data)
     nm <- match(c("id","times"),names(data))
     Z <- cbind(mt,data[,-nm])
  }

  p <- c(expit(as.matrix(Z) %*% coef))
  y <- rbinom(length(p),1,p)

  data <- cbind(y,data)
  data <- count.history(data,status="y",id="id",types=1)
  data <- subset(data,data$Count1<=0)

  attr(data,"coef") <- beta
  return(data)
 }# }}}

##' @export
summary.cumoddsreg <- function(object,...) { ## {{{ 
	ntimes <- object$ntimes
	coefb <- coef(object)[1:ntimes]
	vcovb <- object$var[1:object$ntimes,1:object$ntimes]
	outb <- lava::estimate(coef=coefb,vcov=vcovb,...)

	if (length(object$coef)>ntimes) {
		ll <- length(object$coef)
		takex <- (ntimes+1):ll
		coefx <- object$coef[takex]
		vcovx <- object$var[takex,takex]
     	        outx <- lava::estimate(coef=coefx,vcov=vcovx,...)
	        eoutx <- exp(outx$coefmat[,c(1,3,4)])
	        out <- list(baseline=outb,logor=outx,or=eoutx)
	} else out <- list(baseline=outb)
	return(out)
} ## }}} 


##' @export
print.cumoddsreg <- function(x,...) summary(x,...)

##' @export
vcov.cumoddsreg <- function(object,...) return(object$var) 

##' @export
coef.cumoddsreg <- function(object,...) return(object$coef)

##' @export
predictSurvd <- function(ds,Z,times=1:6,se=FALSE,type="prob")
{# {{{
  if (!is.null(Z)) n <- nrow(Z) 
  if (!is.null(Z)) data <- Z else data <- data.frame(id=1:n)
  Z <- data.frame(Z)
  Z$id <- 1:n
  ccc <- ds$coef

  if (!se) {# {{{{{{
	  data <- Z
	  if (!is.null(times)) {
	     timesf <- data.frame(times=rep(times,n),id=rep(1:n,each=length(times)))
	     data <- merge(data,timesf,by.x="id",by.y="id")
	     mt <- model.matrix(~factor(times),data)
	     nm <- match(c("id","times"),names(data))
	     Z <- cbind(mt,data[,-nm])
	  }
	  if (ncol(Z)!=length(c(ccc))) {
		  print(head(Z))
		  print(ccc)
		  stop("dimension of Z not consistent with length of coefficients"); 
	  }

	  p <- c(expit(as.matrix(Z) %*% ccc))

	  preds <- data.frame(p=p,id=data$id,times=data$times)
	  survt <- exp(cumsumstrata(log(1-preds$p),data$id-1,6))
	  if (type=="prob") pred <- 1-survt
	  if (type=="surv") pred <- survt
	  if (type=="hazard") pred <- p
	  if (type=="rrm") { ## restricted residual mean 
		  ll <- length(survt)
	        pred <- cumsum(c(1,survt[-ll]))
	  }
	  preds <- cbind(preds,pred)
# }}}
  } else {# {{{

    Ft <- function(p)
    {
	   xp <- as.matrix(Zi) %*% p
	   lam <- expit(xp)
	   st <- cumprod(1-lam)
	   if (type=="prob") st <- 1-st 
	   if (type=="surv") st <- st 
	   if (type=="hazard") st <- lam
           if (type=="rrm") { ## restricted residual mean 
		ll <- length(st)
	        st <- cumsum(c(1,st[-ll]))
	   }
	   return(st)
    }

  preds <- c()
  for (i in 1:nrow(Z)) {
     Zi <- data.frame(Z[i,,drop=FALSE])
     data <- Zi
     if (!is.null(times)) {
        timesf <- data.frame(times=rep(times,n),id=rep(1:n,each=length(times)))
        data <- merge(data,timesf,by.x="id",by.y="id")
        mt <- model.matrix(~factor(times),data)
        nm <- match(c("id","times"),names(data))
        Zi <- cbind(mt,data[,-nm])
     }
     if (is.null(ds$var)) covv <- vcov(ds)  else covv <- ds$var
     eud <- estimate(coef=ds$coef,vcov=covv,f=function(p) Ft(p))
     cmat <- data.frame(eud$coefmat)
     cmat$id <- i
     cmat$times <- times
     names(cmat)[1:4] <- c("pred","se","lower","upper")
     preds <- rbind(preds,cmat)
  } 

  }# }}}

return(preds)
}# }}}

## }}} 

##' @export
plotSurvd <- function(ds,ids=NULL,add=FALSE,se=FALSE,cols=NULL,ltys=NULL,...)
{# {{{

 if (is.null(ids)) ids <- unique(ds$id)
 if (is.null(cols)) cols <- 1:length(ids)
 if (is.null(ltys)) ltys <- 1:length(ids)

  k <- 1
  fplot <- 0
  for (i in ids) {
	  timei <- ds$time[ds$id==i]
	  predi <- ds$pred[ds$id==i]

 	  if (fplot==0) {
	  if (!add) plot(timei,predi,type="s",col=cols[k],lty=ltys[k],...)
	  if (add) lines(timei,predi,type="s",col=cols[k],lty=ltys[k],...)
	  fplot <- 1
	  } else lines(timei,predi,type="s",col=cols[k],lty=ltys[k],...)

          if (se) {
	  loweri <- ds$lower[ds$id==i]
	  upperi <- ds$upper[ds$id==i]
	  plotConfRegion(timei,cbind(loweri,upperi),col=cols[k])
	  }
	  k <- k+1
  }

} ## }}} 

