###{{{ phreg0 
###
######phreg0 <- function(X,entry,exit,status,id=NULL,strata=NULL,beta,stderr=TRUE,method="NR",...) {# {{{ 
######  p <- ncol(X) 
######  if (missing(beta)) beta <- rep(0,p)
######  if (p==0) X <- cbind(rep(0,length(exit)))
######  if (!is.null(strata)) { # {{{
######    stratalev <- levels(strata)
######    strataidx <- lapply(stratalev,function(x) which(strata==x))
######    if (!all(unlist(lapply(strataidx,function(x) length(x)>0))))
######      stop("Strata without any observation")
######    dd <- lapply(strataidx, function(ii) {
######        entryi <- entry[ii]
######        trunc <- !is.null(entryi)
######        if (!trunc) entryi <- rep(0,length(exit[ii]))
######                 .Call("FastCoxPrep",
######                       entryi,exit[ii],status[ii],
######                       as.matrix(X)[ii,,drop=FALSE],
######                       id[ii],
######                       trunc,
######                       PACKAGE="mets")
######                 })
######    if (!is.null(id))
######      id <- unlist(lapply(dd,function(x) x$id[x$jumps+1]))
######      obj <- function(pp,U=FALSE,all=FALSE) {
######      val <- lapply(dd,function(d)
######                    with(d,
######                         .Call("FastCoxPL",pp,X,XX,sign,jumps,PACKAGE="mets")))
######      ploglik <-     Reduce("+",lapply(val,function(x) x$ploglik))
######      gradient <-    Reduce("+",lapply(val,function(x) x$gradient))
######      hessian <-     Reduce("+",lapply(val,function(x) x$hessian))
######      if (all) {
######        U <- do.call("rbind",lapply(val,function(x) x$U))
######        hessiantime <- do.call("rbind",lapply(val,function(x) x$hessianttime))
######        time <- lapply(dd,function(x) x$time[x$ord+1])
######        ord <- lapply(dd,function(x) x$ord+1)
######        jumps <- lapply(dd,function(x) x$jumps+1)
######        jumptimes <- lapply(dd,function(x) x$time[x$ord+1][x$jumps+1])
######        S0 <- lapply(val,function(x) x$S0)
######        nevent  <- unlist(lapply(S0,length))
######        return(list(ploglik=ploglik,gradient=gradient,hessian=hessian,
######                    U=U,S0=S0,nevent=nevent,hessianttime=hessiantime,
######                    ord=ord,time=time,jumps=jumps,jumptimes=jumptimes))
######      }
######      structure(-ploglik,gradient=-gradient,hessian=-hessian)
######    }# }}}
######  } else { # {{{
######      trunc <- !is.null(entry)
######      if (!trunc) entry <- rep(0,length(exit))
######      system.time(dd <- .Call("FastCoxPrep",
######                              entry,exit,status,X,
######                              as.integer(seq_along(entry)),
######                              !is.null(entry),
######                              PACKAGE="mets"))
######
######      if (!is.null(id))
######          id <- dd$id[dd$jumps+1]
######      obj <- function(pp,U=FALSE,all=FALSE) {
######          val <- with(dd,
######                      .Call("FastCoxPL",pp,X,XX,sign,jumps,PACKAGE="mets"))
######          if (all) {
######              val$time <- dd$time[dd$ord+1]
######              val$ord <- dd$ord+1
######              val$jumps <- dd$jumps+1
######              val$jumptimes <- val$time[val$jumps]
######              val$nevent <- length(val$S0)
######              return(val)
######          }
######          with(val, structure(-ploglik, gradient=-gradient, hessian=-hessian))
######      }
######  }# }}}
######
######  opt <- NULL
######  if (p>0) {
######      if (tolower(method)=="nr") {
######          opt <- lava::NR(beta,obj,...)
######          opt$estimate <- opt$par
######      } else {
######          opt <- nlm(obj,beta,...)
######          opt$method <- "nlm"
######      }
######      cc <- opt$estimate;  names(cc) <- colnames(X)
######      if (!stderr) return(cc)
######      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
######  } else {
######      val <- obj(0,all=TRUE)
######      val[c("ploglik","gradient","hessian","U")] <- NULL
######  }
######
######  ### computes Breslow estimator 
######  cumhaz <- NULL
######
######  res <- c(val,
######           list(strata=strata,
######                entry=entry,
######                exit=exit,
######                status=status,                
######                p=p,
######                X=X,
######                id=id, opt=opt,cum=cumhaz))
######  class(res) <- "phreg"
######  res
######} # }}}
######
###
###}}} phreg0

###{{{ phreg01

phreg01 <- function(X,entry,exit,status,id=NULL,strata=NULL,
	   offset=NULL,weights=NULL,strata.name=NULL,cumhaz=TRUE,
             beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
	     case.weights=NULL,no.var=0,...) {
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
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
  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z

  ## possible casewights to use for bootstrapping and other things
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

   dd <- .Call("FastCoxPrepStrata", entry,exit,status,X, id, 
	     trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")

   dd$nstrata <- nstrata
   obj <- function(pp,U=FALSE,all=FALSE) {# {{{
      if (length(dd$jumps) >0) {
	      if (is.null(propodds) & is.null(AddGam)) {
		val <- with(dd, .Call("FastCoxPLstrata",pp,X,XX,sign,jumps, strata,nstrata,weights,offset,ZX,caseweights,PACKAGE="mets"))
	      } else if (is.null(AddGam))
		val <- with(dd, .Call("FastCoxPLstrataPO",pp,X,XX,sign,jumps, strata,nstrata,weights,offset,ZX,propodds,PACKAGE="mets"))
	      else val <- with(dd, .Call("FastCoxPLstrataAddGam",pp,X,XX,sign,jumps, strata,nstrata,weights,offset,ZX,
					 AddGam$theta,AddGam$dimthetades,AddGam$thetades,AddGam$ags,AddGam$varlink,AddGam$dimjumprv,AddGam$jumprv,
					 AddGam$JumpsCauses,PACKAGE="mets"))
      } else {
	      val <- list(jumps=NULL,ploglik=NULL,U=NULL,gradient=NULL,hessian=NULL,hessiantime=NULL,S2S0=NULL,E=NULL,S0=NULL)
      }
	 

	  if (all) {
	      val$time <- dd$time
	      ### plus 1 for R index
	      val$ord <- dd$ord+1
	      val$jumps <- dd$jumps+1
	      val$jumptimes <- val$time[val$jumps]
	      val$weightsJ <- dd$weights[val$jumps]
	      val$case.weights <- dd$case.weights[val$jumps]
	      val$nevent <- length(val$S0)
	      val$nstrata <- dd$nstrata
	      val$strata <- dd$strata
	      val$strata.jumps <- val$strata[val$jumps]
	      return(val)
	  }
          n <- length(dd$time)
	 with(val,structure(-ploglik/n,gradient=-gradient/n,hessian=-hessian/n))
	}# }}}

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
      cc <- opt$estimate;  names(cc) <- colnames(X)
      if (!stderr) return(cc)
      val <- c(list(coef=cc), obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta), obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }

  se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL
  II <- NULL
  if (no.opt==FALSE & p!=0) {
         II <- - tryCatch(solve(val$hessian),error=
             function(e) matrix(0,nrow(val$hessian),ncol(val$hessian)) )
  } else II <- matrix(0,p,p)

  ### computes Breslow estimator 
  if (cumhaz==TRUE & (length(val$jumps)>0)) { # {{{
	 strata <- val$strata[val$jumps]
	 nstrata <- val$nstrata
	 jumptimes <- val$jumptimes

	 ## Brewslow estimator
	 cumhaz <- cbind(jumptimes,cumsumstrata(1/val$S0,strata,nstrata))
	 if ((no.opt==FALSE & p!=0)) { 
	     DLambeta.t <- apply(val$E/c(val$S0),2,cumsumstrata,strata,nstrata)
	     varbetat <-   rowSums((DLambeta.t %*% II)*DLambeta.t)
	 ### covv <-  apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
	 } else varbetat <- 0
	 var.cumhaz <- cumsumstrata(1/val$S0^2,strata,nstrata)+varbetat
	 se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)

	 colnames(cumhaz)    <- c("time","cumhaz")
	 colnames(se.cumhaz) <- c("time","se.cumhaz")
 } # }}} 
 else {cumhaz <- se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL}

  res <- c(val,
           list(cox.prep=dd,
                strata.call=strata.call, strata.level=strata.level,
                entry=entry,
                exit=exit,
                status=status,                
                p=p,
                X=X,
                offsets=offset,
                weights=weights,
                id=id.orig, call.id=call.id,
                opt=opt,
                no.opt=no.opt,
                cumhaz=cumhaz, se.cumhaz=se.cumhaz,
                lcumhaz=lcumhaz, lse.cumhaz=lse.cumhaz,
                ihessian=II,
                II=II,strata.name=strata.name,propodds=propodds))
  class(res) <- "phreg"

  ## also computing robust variance 
  if (p>0 & no.var==0) {
  ii <- iid(res)
  phvar <- crossprod(ii)
  colnames(phvar) <- rownames(phvar) <- names(res$coef)
  res$var <- phvar
  } else res$var <- 0

  return(res)
}

###}}} phreg0

###{{{ phreg

##' Fast Cox PH regression
##'
##' Fast Cox PH regression
##' Robust variance is default variance with the summary. 
##'
##' influence functions (iid) will follow numerical order of given cluster variable
##' so ordering after $id will give iid in order of data-set.
##'
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param offset offsets for cox model
##' @param weights weights for Cox score equations
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases phreg phreg.par robust.phreg readPhreg IIDbaseline.phreg
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' out2 <- phreg(Event(time,status)~vf+chf+strata(wmicat.4),data=TRACE)
##' ## tracesim <- timereg::sim.cox(out1,1000)
##' ## sout1 <- phreg(Surv(time,status==1)~vf+chf+strata(wmicat.4),data=tracesim)
##' ## robust standard errors default 
##' summary(out1)
##' out1 <- phreg(Surv(time,status!=0)~vf+chf+strata(wmicat.4),data=TRACE)
##' summary(out2)

##' 
##' par(mfrow=c(1,2))
##' bplot(out1)
##' ## bplot(sout1,se=TRUE)
##' 
##' ## computing robust variance for baseline
##' rob1 <- robust.phreg(out1)
##' bplot(rob1,se=TRUE,robust=TRUE)
##' 
##' ## making iid decomposition of regression parameters
##' betaiiid <- lava::iid(out1)
##' 
##' ## making iid decomposition of baseline at a specific time-point
##' Aiiid <- mets:::IIDbaseline.phreg(out1,time=30)
##' 
##' @export
phreg <- function(formula,data,offset=NULL,weights=NULL,...) {# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
###  if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
  if (!inherits(Y,c("Event","Surv"))) stop("Expected a 'Surv' or 'Event'-object")
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
  if (!is.null(attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}
###  if (!is.null(attributes(Terms)$specials$offset)) {
###    ts <- survival::untangle.specials(Terms, "offset")
###    pos.offset <- ts$terms
###    Terms  <- Terms[-ts$terms]
###    offset <- m[[ts$vars]]
###  }  else pos.offset <- NULL
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  res <- c(phreg01(X,entry,exit,status,id,strata,offset,weights,strata.name,...),
   list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,
	cluster.pos=pos.cluster,n=length(exit)))
  class(res) <- "phreg"
  
  res
}# }}}

##' @export
readPhreg <- function (object, newdata, nr=TRUE, ...)
{# {{{
     exit <- entry <- status  <- clusters <- NULL
     if (missing(newdata)) { # {{{
         X <- object$X
         strataNew <- object$strata
	 if (!nr) { 
	 exit <- object$exit; entry <- object$entry; status <- object$status; 
	 clusters <- object$id
	 }
     } else { ## make design for newdata
       xlev <- lapply(object$model.frame,levels)
       ff <- unlist(lapply(object$model.frame,is.factor))
       upf <- update(object$formula,~.)
       tt <- terms(upf)
       if (nr) tt <- delete.response(tt)
       X <- model.matrix(tt,data=newdata,xlev=xlev)[,-1,drop=FALSE]
       if (!nr) {
         allvar <- all.vars(tt)
	 pr <- length(allvar)-ncol(object$model.frame)+1
	 if (pr==2) { 
	    exit <- newdata[,allvar[1]]
	    status <- newdata[,allvar[2]]
	 } else {
	    entry <- newdata[,allvar[1]]
	    exit <- newdata[,allvar[2]]
	    status <- newdata[,allvar[3]]
	 }
       }
       clusterTerm<- grep("^cluster[(][A-z0-9._:]*",colnames(X),perl=TRUE)
       ## remove clusterTerm from design
       if (length(clusterTerm)==1) { 
	       clusters <- X[,clusterTerm]
	       X <- X[,-clusterTerm,drop=FALSE]
	       id <- clusters
               id.orig <- id; 
               if (!is.null(id)) {
	          ids <- unique(id)
	          nid <- length(ids)
                  if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
                        id <- as.integer(factor(id,labels=seq(nid)))-1
                  }
		  clusters <- id
           } else clusters <- (1:nrow(newdata))-1
       } 
       strataTerm<- grep("^strata[(][A-z0-9._:]*",colnames(X),perl=TRUE)
       ## remove strataTerm from design, and construct numeric version of strata
       if (length(strataTerm)>=1) { 
	       strataNew <- X[,strataTerm,drop=FALSE]
               whichstrata  <-  paste(object$strata.name,object$strata.level,sep="")
	       if (length(strataTerm)>=1) { ## construct strata levels numeric
               mm <- match(colnames(X)[strataTerm], whichstrata)-1
               strataNew <- c(strataNew %*% mm )
               }
	       X <- X[,-strataTerm,drop=FALSE]
       } else strataNew <- rep(0,nrow(X))
     }# }}}
return(list(X=X,strata=strataNew,entry=entry,exit=exit,status=status,clusters=clusters))
}# }}}

###}}} phreg

###{{{ phregR


##' Fast Cox PH regression and calculations done in R to make play and adjustments easy 
##'
##' Fast Cox PH regression with R implementation to play and adjust in R function: FastCoxPLstrataR
##'
##' Robust variance is default variance with the summary. 
##'
##' influence functions (iid) will follow numerical order of given cluster variable
##' so ordering after $id will give iid in order of data-set.
##'
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param offset offsets for cox model
##' @param weights weights for Cox score equations
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases  FastCoxPLstrataR
##' 
##' @export
phregR <- function(formula,data,offset=NULL,weights=NULL,...) {# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
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
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  res <- c(phreg01R(X,entry,exit,status,id,strata,offset,weights,strata.name,...),
   list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster,
	n=length(exit)))
  class(res) <- "phreg"
  
  res
}# }}}

phreg01R <- function(X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
             strata.name=NULL,cumhaz=TRUE,
             beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
	     case.weights=NULL,...) {# {{{
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
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
  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z

  ## possible casewights to use for bootstrapping and other things
  if (is.null(case.weights)) case.weights <- rep(1,length(exit)) 

  trunc <- (!is.null(entry))
  if (!trunc) entry <- rep(0,length(exit))

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1; 
   ## orginal id coding into integers 
   id.orig <- id+1; 

   dd <- .Call("FastCoxPrepStrata", entry,exit,status,X, id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")

   dd$nstrata <- nstrata

	obj <- function(pp,U=FALSE,all=FALSE) {# {{{
	if (is.null(propodds) & is.null(AddGam)) 
		cc <- system.time(
	  val <- with(dd, FastCoxPLstrataR(pp,X,XX,sign,jumps, strata,nstrata,weights,offset,ZX,caseweights))
	  )
         else if (is.null(AddGam)) 
		 val <- with(dd, .Call("FastCoxPLstrataPO",pp,X,XX,sign,jumps, strata,nstrata,weights,offset,ZX,propodds,PACKAGE="mets"))
	 else val <- with(dd, .Call("FastCoxPLstrataAddGam",pp,X,XX,sign,jumps, strata,nstrata,weights,offset,ZX,
		    AddGam$theta,AddGam$dimthetades,AddGam$thetades,AddGam$ags,AddGam$varlink,AddGam$dimjumprv,AddGam$jumprv,AddGam$JumpsCauses,PACKAGE="mets"))
	 
	  if (all) {
	      val$time <- dd$time
	      val$ord <- dd$ord+1
	      val$jumps <- dd$jumps+1
	      val$jumptimes <- val$time[val$jumps]
	      val$weightsJ <- dd$weights[val$jumps]
	      val$case.weights <- dd$case.weights[val$jumps]
	      val$strata.jumps <- val$strata[val$jumps]
	      val$nevent <- length(val$S0)
	      val$nstrata <- dd$nstrata
	      val$strata <- dd$strata
	      return(val)
	  } 
	 with(val,structure(-ploglik,gradient=-gradient,hessian=-hessian))
	}# }}}

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
      cc <- opt$estimate;  names(cc) <- colnames(X)
      if (!stderr) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }

  se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL
  II <- NULL
  ### computes Breslow estimator 
  if (cumhaz==TRUE) { # {{{
	 if (no.opt==FALSE & p!=0) {
               II <- - tryCatch(solve(val$hessian),error=
	              function(e) matrix(0,nrow(val$hessian),ncol(val$hessian)) )
	 } else II <- matrix(0,p,p)
	 strata <- val$strata[val$jumps]
	 nstrata <- val$nstrata
	 jumptimes <- val$jumptimes

	 ## Brewslow estimator
	 cumhaz <- cbind(jumptimes,cumsumstrata(1/val$S0,strata,nstrata))
	 if ((no.opt==FALSE & p!=0)) { 
	     DLambeta.t <- apply(val$E/c(val$S0),2,cumsumstrata,strata,nstrata)
	     varbetat <-   rowSums((DLambeta.t %*% II)*DLambeta.t)
	 ### covv <-  apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
	 } else varbetat <- 0
	 var.cumhaz <- cumsumstrata(1/val$S0^2,strata,nstrata)+varbetat
	 se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)

	 colnames(cumhaz)    <- c("time","cumhaz")
	 colnames(se.cumhaz) <- c("time","se.cumhaz")
 } # }}} 
 else {cumhaz <- se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL}

  res <- c(val,
           list(cox.prep=dd,
		strata.call=strata.call, strata.level=strata.level,
                entry=entry,
                exit=exit,
                status=status,                
                p=p,
                X=X,
		offsets=offset,
		weights=weights,
                id=id.orig, 
		opt=opt, 
		cumhaz=cumhaz, se.cumhaz=se.cumhaz,
		lcumhaz=lcumhaz, lse.cumhaz=lse.cumhaz,
		II=II,strata.name=strata.name,propodds=propodds))
  class(res) <- "phreg"

  ## also computing robust variance 
  if (p>0) {
  phvar <- crossprod(iid(res))
  colnames(phvar) <- rownames(phvar) <- names(res$coef)
  res$var <- phvar
  } else res$var <- 0

  return(res)
}# }}}

##' @export
FastCoxPLstrataR <- function(beta, X, XX, Sign, Jumps, strata, nstrata, weights, offsets, ZX, caseweights) 
{# {{{
	p=length(beta)
	strata=c(strata)
	Xb = c(X %*% beta+offsets)
	eXb = c(exp(Xb)*weights);
	if (nrow((Sign))==length(eXb)) { ## Truncation
		eXb = c(Sign)*eXb;
	}

	S0 = c(revcumsumstrata(eXb,strata,nstrata))
	E=apply(eXb*X,2,revcumsumstrata,strata,nstrata)/S0; 
	Jumps=Jumps+1

	E = E[Jumps,];
        E2=.Call("vecMatMat",E,E)$vXZ;  

	XX2=apply(XX*eXb,2,revcumsumstrata,strata,nstrata)/S0; 
##	mat XX2=revcumsumstrataMatCols(XX,eXb,S0,strata,nstrata); 

	XX2 = XX2[Jumps,];
	weightsJ=weights[Jumps];  
	caseweightsJ=caseweights[Jumps];  
	S0 = S0[Jumps];
	grad = (X[Jumps,]-E);                         ## Score
	val =  (Xb[Jumps]-log(S0));                   ## Partial log-likelihood

	##  colvec iweightsJ=1/weightsJ; 
	S02 = S0/(caseweightsJ*weightsJ);             ## S0 with weights to estimate baseline 
	grad2= grad*(caseweightsJ*weightsJ);          ## score  with weights
	gradient=apply(grad2,2,sum)
	val2 = caseweightsJ*weightsJ*val;             ## Partial log-likelihood with weights

	hesst = -(XX2-E2);                            ## hessian contributions in jump times 
	hesst2 = hesst*(caseweightsJ*weightsJ);       ## hessian over time with weights 
	hess2 = matrix(apply(hesst2,2,sum),p,p);      ## hessian with weights 


	out=list(jumps=Jumps, ploglik=sum(val2),U=grad2, 
		 gradient=matrix(gradient,1,p), hessian=hess2, hessianttime=hesst2, S2S0=XX2, E=E, S0=S02 )
	return(out)
}# }}}

###}}} 

###{{{ simcox

##' @export
simCox <- function(n=1000, seed=1, beta=c(1,1), entry=TRUE) {
  if (!is.null(seed))
      set.seed(seed)
  m <- lava::lvm()
  lava::regression(m,T~X1+X2) <- beta
  lava::distribution(m,~T+C) <- lava::coxWeibull.lvm(scale=1/100)
  lava::distribution(m,~entry) <- lava::coxWeibull.lvm(scale=1/10)
  m <- lava::eventTime(m,time~min(T,C=0),"status")
  d <- lava::sim(m,n);
  if (!entry) d$entry <- 0
  else d <- subset(d, time>entry,select=-c(T,C))
  return(d)
}

###}}} simcox

###{{{ vcov

##' @export
vcov.phreg  <- function(object,...) {    
 if ((length(class(object))==1) & inherits(object,"phreg")) {
  res <- as.matrix(object$var)  ### objectcrossprod(ii <- iid(object,...))
###  attributes(res)$ncluster <- attributes(ii)$ncluster
###  attributes(res)$invhess <- attributes(ii)$invhess
  colnames(res) <- rownames(res) <- names(coef(object))
} else { ##if ((length(class(object))==2) & class(object)[2]=="cifreg") {
  res <- as.matrix(object$var)
  colnames(res) <- rownames(res) <- names(coef(object))
}
  res
}

###}}} vcov

###{{{ coef

##' @export
coef.phreg  <- function(object,...) {
  object$coef
}

###}}} coef

###{{{ iid & Robust variances 

##' @export
IC.phreg  <- function(x,type="robust",all=FALSE,baseline=FALSE,...) {# {{{
  if (baseline) {
    if (inherits(x, "cifreg")) {
      res <- IIDbaseline.cifreg(x, ...)$base.iid
    } else {
      res <- IIDbaseline.phreg(x, ...)$base.iid
    }
    return(res*NROW(res))
  }
  classes1 <- "mlogit"
  if ((length(class(x))==1) || inherits(x, classes1)) {
    invhess <- -solve(x$hessian)
    orig.order <- FALSE

if (is.null(x$propodds)) {
  if (type=="robust") {	# {{{ cox model 
	  xx <- x$cox.prep
	  ii <- invhess 
	  S0i <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/x$S0
	  Z <- xx$X
	  U <- E <- matrix(0,nrow(xx$X),x$p)
	  E[xx$jumps+1,] <- x$E
	  U[xx$jumps+1,] <- x$U
	  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
	  ### Martingale  as a function of time and for all subjects to handle strata 
	  MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-EdLam0)*rr*c(xx$weights)
	  if (orig.order) {
	     oo <- (1:nrow(xx$X))[xx$ord+1]
	     oo <- order(oo)
	     ### back to order of iid variable 
	     MGt <- MGt[oo,,drop=FALSE]   ## sum after id later so not needed
	     id <- xx$id[oo]
	  } else id <-  xx$id
  } else  { 
     MGt <- x$U; MG.base <- 1/x$S0; 
  }# }}}
} else { # {{{  prop-odds model logitSurv
    xx <- x$cox.prep
    ii <- invhess 
    S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1] <- 1/x$S0
    Z <- xx$X
    U <- E <- matrix(0,nrow(xx$X),x$p)
    E[xx$jumps+1,] <- x$E
    U[xx$jumps+1,] <- x$U
    cumhazA <- cumsumstratasum(S0i,xx$strata,xx$nstrata,type="all")
    cumhaz <- c(cumhazA$sum)
    cumhazm <- cumhazA$lagsum
    ###
    EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
    rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
    rro <- c(exp(Z %*% coef(x) + xx$offset))
    S0star <- revcumsumstrata(rr/(1+rro*cumhazm),xx$strata,xx$nstrata)
    S0 <- revcumsumstrata(rr,xx$strata,xx$nstrata)
    S1 <- apply(Z*rr,2,revcumsumstrata,xx$strata,xx$nstrata)
    Et <- S1/c(S0)
    lt <- apply((Z-Et)*c(rr*rro/(1+rro*cumhazm)),2,revcumsumstrata,xx$strata,xx$nstrata)
    Estar <- S0star/S0
    EstardLam <- cumsumstrata(Estar*S0i,xx$strata,xx$nstrata)
    k <- exp(-EstardLam)
    basecor <- apply(lt*c(k*S0i),2,revcumsumstrata,xx$strata,xx$nstrata)
    basecor <- basecor/c(k*S0)
    www <- x$propoddsW*x$weightsJ
    U[xx$jumps+1,] <- U[xx$jumps+1,]-c(www)*basecor[xx$jumps+1,]

    MGt <- U; 
    if (type=="robust") {	
       baseDLam0 <- apply(basecor*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
       ### Martingale  as a function of time and for all subjects to handle strata 
       MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0-baseDLam0)*rr*c(xx$weights) ### -baseIII
     }  else MGt <- MGt[xx$jumps+1,]
     MG.base <- 1/x$S0; 
     id <-  xx$id
  }# }}}

  ncluster <- NULL
  if (type=="robust" & (!is.null(x$id) | any(x$entry>0))) {
    if (type=="martingale") id <- x$id[x$jumps]
    ###  ii <- mets::cluster.index(id)
    UU <- apply(MGt,2,sumstrata,id,max(id)+1)
    ### names of clusters given in call 
    ncluster <- nrow(UU)
  } else {
     UU <- MGt
  }
 res <-  structure(UU%*%invhess,invhess=invhess,ncluster=ncluster)
 res <- res*NROW(res)
 return(res)
} else if (inherits(x,c("cifreg","recreg"))) {
  res <- x$iid*NROW(x$iid)
  return(res)
}
} # }}}

##' @export IIDbaseline.phreg 
IIDbaseline.phreg <- function(x,time=NULL,ft=NULL,fixbeta=NULL,...)
{# {{{
###  sum_i int_0^t f(s)/S_0(s) dM_{ki}(s) - P(t) \beta_k
###  with possible strata and cluster "k", and i in clusters 
  if (!inherits(x,"phreg")) stop("Must be phreg object\n"); 
  if (is.null(time)) stop("Must give time for iid of baseline")

  if (!is.null(x$propodds))  stop("Only for Cox model") 
  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if ((x$no.opt) | is.null(x$coef)) fixbeta<- 1 else fixbeta <- 0

  xx <- x$cox.prep
  btimexx <- c(1*(xx$time < time))

  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <- 1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
  U[xx$jumps+1,] <- x$U
  ###    
  cumhaz <- c(cumsumstrata(S0i,xx$strata,xx$nstrata))
  ### only up to time t
  if (is.null(ft))  ft <- rep(1,length(xx$time))
  cumS0i2 <- c(cumsumstrata(ft*S0i2*btimexx,xx$strata,xx$nstrata))
  if (fixbeta==0) {
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  Ht <- apply(ft*E*S0i*btimexx,2,cumsumstrata,xx$strata,xx$nstrata)
  } else Ht <- NULL
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset)) else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  MGAiid <- ft*S0i*btimexx-cumS0i2*rr*c(xx$weights)

  MGtiid <- NULL
  if (fixbeta==0) {# {{{
     invhess <- -solve(x$hessian)
     MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx$weights)
     MGt <- MGt %*% invhess
     MGtiid <- apply(MGt,2,sumstrata,id,mid)
     ## Ht efter strata
     Htlast <- tailstrata(xx$strata,xx$nstrata)
     HtS <- Ht[Htlast,,drop=FALSE]
  }# }}}

 ### \hat beta - \beta = \sum_i \beta_i  (iid) 
 ### iid after baseline:
 ### \hat A_s-A_s=\sum_{i clusters} \sum_{j: i(j)=i, s(j)=s} \int_0^t 1/S_0 dM^s_{i,j} - P^s(t) \sum_i \beta_i
 ### = \sum_{i clusters} ( \sum_{j \in i(j)=i, s(j)=s} \int_0^t 1/S_0 dM^s_j - P^s(t) \beta_i ) 

 ## sum after id's within strata and order 
 MGAiids <- c()
 sus <- sort(unique(xx$strata))
 for (i in sus)  { 
	 wi <- which(xx$strata==i)
         MGAiidl <- sumstrata(MGAiid[xx$strata==i],xx$id[xx$strata==i],mid)

        if (fixbeta==0) {
           UU <-  apply(HtS[i+1,]*t(MGtiid),2,sum)
           MGAiidl <- MGAiidl - UU
         }
         MGAiids <- cbind(MGAiids,MGAiidl)
 }
 MGAiid <- MGAiids
 colnames(MGAiid) <- paste("strata",sus,sep="")

 return(list(time=time,base.iid=MGAiid,strata=xx$strata,nstrata=xx$nstrata,
	     beta.id=id,beta.iid=MGtiid,model.frame=x$model.frame,formula=x$formula))
} # }}}


##' @export
residuals.phreg  <- function(object,cumsum=FALSE,...) {# {{{
  orig.order <- FALSE
  x <- object

if (is.null(x$propodds)) { # {{{ cox model 
	  xx <- x$cox.prep
	  dN <- S0i <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/x$S0
	  dN[xx$jumps+1] <- 1 
	  cumhaz <- cumsumstrata(S0i,xx$strata,xx$nstrata)
	  Z <- xx$X
	  if (is.null(coef(x))) 
		  rr <- c(xx$sign*exp(xx$offset))
          else 
	  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
	  ### dMartingale  as a function of time and for all subjects to handle strata 
	  Lamt <- (cumhaz)*rr*c(xx$weights)
	  dMGt <- dN-Lamt
	  if (orig.order) {
	     oo <- (1:nrow(xx$X))[xx$ord+1]
	     oo <- order(oo)
	     ### back to order of iid variable 
	     dMGt <- dMGt[oo,,drop=FALSE]   ## sum after id later so not needed
	     id <- xx$id[oo]
	  } else id <-  xx$id
} else { ## }}}
	# {{{  prop-odds model logitSurv
	  xx <- x$cox.prep
	  dN <- S0i <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/x$S0
	  dN[xx$jumps+1] <- 1 
	  U <- E <- matrix(0,nrow(xx$X),x$p)
	  E[xx$jumps+1,] <- x$E
	  U[xx$jumps+1,] <- x$U
	  S0i <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/x$S0
	  cumhazA <- cumsumstratasum(S0i,xx$strata,xx$nstrata,type="all")
	  cumhaz <- c(cumhazA$sum)
	  cumhazm <- cumhazA$lagsum
          ###
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
	  rro <- c(exp(Z %*% coef(x) + xx$offset))
          S0star <- revcumsumstrata(rr/(1+rro*cumhazm),xx$strata,xx$nstrata)
          S0 <- revcumsumstrata(rr,xx$strata,xx$nstrata)
          S1 <- apply(Z*rr,2,revcumsumstrata,xx$strata,xx$nstrata)
	  Et <- S1/c(S0)
          lt <- apply((Z-Et)*c(rr*rro/(1+rro*cumhazm)),2,revcumsumstrata,xx$strata,xx$nstrata)
	  Estar <- S0star/S0
	  EstardLam <- cumsumstrata(Estar*S0i,xx$strata,xx$nstrata)
	  k <- exp(-EstardLam)
	  basecor <- apply(lt*c(k*S0i),2,revcumsumstrata,xx$strata,xx$nstrata)
	  basecor <- basecor/c(k*S0)
	  www <- x$propoddsW*x$weightsJ
	  U[xx$jumps+1,] <- U[xx$jumps+1,] - c(www)* basecor[xx$jumps+1,]
	  baseDLam0 <- apply(basecor*S0i,2,cumsumstrata,xx$strata,xx$nstrata)

	  ### Martingale  as a function of time and for all subjects to handle strata 
	  MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0-baseDLam0)*rr*c(xx$weights)

	  if (orig.order) {
	     oo <- (1:nrow(xx$X))[xx$ord+1]
	     oo <- order(oo)
	     ### back to order of data-set
	     MGt <- MGt[oo,,drop=FALSE]  
	     id <- xx$id[oo]
	  } else id <-  xx$id
  }# }}}

  ncluster <- NULL
  if (!cumsum) {
	  mid <- max(id)+1
    Mt <- sumstrata(dMGt,id,mid)
    cumhaz <- sumstrata(Lamt,id,mid)
    out <- list(residuals=c(Mt),cumhaz=c(cumhaz))
  } else {
     mid <- max(id)+1
     out <- data.frame(dMGt=c(dMGt),Lamt=Lamt,status=dN,
		time=c(xx$time),id=c(xx$id)+1,sign=xx$sign)
     dsort(out) <- ~time+status+id-sign
     Mt <- cumsumstrata(out$dMGt,out$id-1,mid)
     cumhaz <- cumsumstrata(out$Lamt,out$id-1,mid)
    out <- cbind(out,Mt,cumhaz)
    out <- subset(out,sign==1)
    ddrop(out)  <- ~sign+dMGt+Lamt 
  }
  
  return(out)
} # }}}

##' @export
robust.basehaz.phreg  <- function(x,type="robust",fixbeta=NULL,...) {# {{{

  IsdM <- squareintHdM(x,ft=NULL,fixbeta=fixbeta,...)
  varA <-   IsdM$varInt[x$jumps]
  strata <- x$strata[x$jumps]
  cumhaz <- x$cumhaz
  se.cumhaz <- cbind(cumhaz[,1],varA^.5)
  colnames(se.cumhaz) <- c("time","se.cumhaz")
  
  return(list(cumhaz=cumhaz,se.cumhaz=se.cumhaz,strata=strata))
} # }}}

##' @export robust.phreg
robust.phreg  <- function(x,fixbeta=NULL,...) {

  if (is.null(fixbeta)) 
  if ((x$no.opt) | is.null(x$coef)) fixbeta<- 1 else fixbeta <- 0

 if (fixbeta==0)  {
    beta.iid <- iid(x)
    robvar <- crossprod(beta.iid)
 } else robvar <- beta.iid <- NULL
 baseline <- robust.basehaz.phreg(x,fixbeta=fixbeta,...); 
 ## add arguments so that we can call basehazplot.phreg
 res <- c(x,list(beta.iid=beta.iid,robvar=robvar,robse.cumhaz=baseline$se.cumhaz))
 class(res) <- "phreg"
 return(res)
}

###}}}

###{{{ summary

##' @export
summary.phreg <- function(object,type=c("robust","martingale"),augment.type=c("var.augment.times","var.augment"),...) {
  expC <- cc <- ncluster <- V <- NULL

   if (length(object$p)>0 & object$p>0 & (!object$no.opt)) {
    I <- -solve(object$hessian)
    if ( (length(class(object))==2) && ( inherits(object,c("cifreg","recreg")))) {
	    V <- object$var
	    ncluster <- object$ncluster ## nrow(object$Uiid)
            if (!is.null(object$augmentation)) { V <- object[[augment.type[1]]]; 
	    }
    } else  {  ## phreg
	    V <- vcov(object,type=type[1])
            ncluster <- object$n
    }
    cc <- cbind(coef(object),diag(V)^0.5,diag(I)^0.5)
    cc  <- cbind(cc,2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
    colnames(cc) <- c("Estimate","S.E.","dU^-1/2","P-value")
    if (length(class(object))==1) if (!is.null(ncluster <- attributes(V)$ncluster))
    rownames(cc) <- names(coef(object))
    expC <- exp(lava::estimate(coef=coef(object),vcov=V)$coefmat[,c(1,3,4),drop=FALSE])
  } 
  Strata <- levels(object$strata)
  n <- object$n
  res <- list(coef=cc,n=n,nevent=object$nevent,strata=Strata,ncluster=ncluster,var=V,exp.coef=expC)
  class(res) <- "summary.phreg"
  res
}

###}}} summary

###{{{ print.summary

##' @export
print.summary.phreg  <- function(x,max.strata=5,...) {

  if (length(class(x))==2 & inherits(x,"cifreg")) cat("Competing risks regression \n"); 
  if (!is.null(x$propodds)) { 
       cat("Proportional odds model, log-OR regression \n"); 
  } else cat("\n")

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
    cat("exp(coeffients):\n")
    printCoefmat(x$exp.coef,...)
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

}

###}}} print.summary

## {{{ Utility functions, cumsumstrata ....

##' @export
tailstrata <- function(strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("tailstrataR",length(strata),strata,nstrata,PACKAGE="mets")
if (any(res$found<0.5))  { warning("Not all strata found");  cat((1:nstrata)[res$found>0.5]); }
return(res$where)
}# }}}

##' @export
headstrata <- function(strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("headstrataR",length(strata),strata,nstrata,PACKAGE="mets")
if (any(res$found<0.5))  { warning("Not all strata found");  cat((1:nstrata)[res$found>0.5]); }
return(res$where)
}# }}}

##' @export
sumstrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (length(x)!=length(strata)) stop("length of x and strata must be same\n"); 
res <- .Call("sumstrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}

##' @export
cumsumstrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (length(x)!=length(strata)) stop("length of x and strata must be same\n"); 
res <- .Call("cumsumstrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}

##' @export
diffstrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (length(x)!=length(strata)) stop("length of x and strata must be same\n"); 
res <- .Call("diffstrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}


##' @export
revcumsumstrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (length(x)!=length(strata)) stop("length of x and strata must be same\n"); 
res <- .Call("revcumsumstrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}

##' @export
revcumsum2strata <- function(x,strata,nstrata,strata2,nstrata2,lag=FALSE)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (any(strata2<0) | any(strata2>nstrata2-1)) stop("strata2 index not ok\n"); 
if (length(x)!=length(strata))  stop("length of x and strata must be same\n"); 
if (length(x)!=length(strata2)) stop("length of x and strata2 must be same\n"); 
res <- .Call("revcumsum2strataR",as.double(x),strata,nstrata,strata2,nstrata2,PACKAGE="mets")
return(res)
}# }}}

##' @export
revcumsum2stratafdN <- function(x,y,strata,nstrata,strata2,nstrata2,startx)
{# {{{
if (any(strata<0) | any(strata>nstrata-1))    stop("strata index not ok\n"); 
if (any(strata2<0) | any(strata2>nstrata2-1)) stop("strata2 index not ok\n"); 
if (length(x)!=length(strata))  stop("length of x and strata must be same\n"); 
if (length(x)!=length(strata2)) stop("length of x and strata2 must be same\n"); 
if (length(x)!=length(y)) stop("length of x and y must be same\n"); 
res <- .Call("revcumsum2stratafdNR",as.double(x),as.double(y),strata,nstrata,strata2,nstrata2,as.double(startx),PACKAGE="mets")
return(res)
}# }}}


##' @export
cumsum2strata <- function(x,y,strata,nstrata,strata2,nstrata2,startx)
{# {{{
if (any(strata<0) | any(strata>nstrata-1))    stop("strata index not ok\n"); 
if (any(strata2<0) | any(strata2>nstrata2-1)) stop("strata2 index not ok\n"); 
if (length(x)!=length(strata))  stop("length of x and strata must be same\n"); 
if (length(x)!=length(strata2)) stop("length of x and strata2 must be same\n"); 
res <- .Call("cumsum2strataR",as.double(x),as.double(y),strata,nstrata,strata2,nstrata2,as.double(startx),PACKAGE="mets")
return(res)
}# }}}


##' @export
vecAllStrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (length(x)!=length(strata))  stop("length of x and strata must be same\n"); 
res <- .Call("vecAllStrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}



##' @export
revcumsum <- function(x)
{# {{{
res <- .Call("revcumsumR",x,PACKAGE="mets")$res
return(res)
}# }}}

##' @export
revcumsumstratasum <- function(x,strata,nstrata,type="all")
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (length(x)!=length(strata)) stop("length of x and strata must be same\n"); 
if (type=="sum")    res <- .Call("revcumsumstratasumR",x,strata,nstrata)$sum
if (type=="lagsum") res <- .Call("revcumsumstratasumR",x,strata,nstrata)$lagsum
if (type=="all")    res <- .Call("revcumsumstratasumR",x,strata,nstrata)
return(res)
}# }}}

##' @export
cumsumstratasum <- function(x,strata,nstrata,type="all")
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (length(x)!=length(strata)) stop("length of x and strata must be same\n"); 
if (type=="sum")    res <- .Call("cumsumstratasumR",x,strata,nstrata)$sum
if (type=="lagsum") res <- .Call("cumsumstratasumR",x,strata,nstrata)$lagsum
if (type=="all")    res <- .Call("cumsumstratasumR",x,strata,nstrata)
return(res)
}# }}}

##' @export
matdoubleindex <- function(x,rows,cols,xvec=NULL)
{# {{{
if (!is.matrix(x)) stop("x must be matrix")
ncols <- ncol(x); nrows <- nrow(x)
## to avoid warnings when going to C, get rid of Inf 
cols[cols==Inf] <- ncol(x)+1; 
rows[rows==Inf] <- nrow(x)+1; 
if (length(rows)==1) rows <- rep(rows,length(cols))
if (length(cols)==1) cols <- rep(cols,length(rows))
if (length(cols)!=length(rows)) stop("rows and cols different lengths\n"); 
if (is.null(xvec)) { assign <- 0; xvec <- 1} else { 
	assign <- 1; 
        if (length(cols)!=length(xvec)) stop("rows and cols and xvec differ \n"); 
}
res <- .Call("Matdoubleindex",x,rows-1,cols-1,length(cols),assign,xvec)$mat
return(res)
}# }}}

##' @export
mdi <- function(x,...) matdoubleindex(x,...)

##' @export
covfr  <- function(x,y,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("covrfR",x,y,strata,nstrata)
return(res)
}# }}}

##' @export
revcumsumidstratasum <- function(x,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")    res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)$sum
if (type=="lagsum")    res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)$lagsum
if (type=="lagsumsquare") res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)$lagsumsquare
if (type=="all")    res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
revcumsumidstratasumCov <- function(x,y,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")    res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$sum
if (type=="lagsum")    res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$lagsum
if (type=="lagsumsquare") res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$lagsumsquare
if (type=="all")    res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
cumsumidstratasumCov <- function(x,y,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")   res <- .Call("cumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$sum
else res <- .Call("cumsumidstratasumCovR",x,y,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
cumsumidstratasum <- function(x,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")   res <- .Call("cumsumidstratasumR",x,id,nid,strata,nstrata)$sum
else res <- .Call("cumsumidstratasumR",x,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
covfridstrata  <- function(x,y,id,nid,strata,nstrata)
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("covrfstrataR",x,y,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
covfridstrataCov  <- function(x,y,x1,y1,id,nid,strata,nstrata)
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("covrfstrataCovR",x,y,x1,y1,id,nid,strata,nstrata)
return(res)
}# }}}


## }}}

# {{{ Restricted mean for stratified Kaplan-Meier with martingale standard errors 

##' Restricted mean for stratified Kaplan-Meier or Cox model with martingale standard errors 
##' 
##' Restricted mean for stratified Kaplan-Meier or stratified Cox with martingale 
##' standard error. Standard error is computed using linear interpolation between 
##' standard errors at jump-times. Plots gives restricted mean at all times. 
##' Years lost can be computed based on this and decomposed into years lost for
##' different causes using the cif.yearslost function that is based on  
##' integrating the cumulative incidence functions.  
##' One particular feature of these functions are that the restricted mean and years-lost are 
##' computed for all event times as functions and can be viewed.  When times are given and beyond
##' the last event time withn a strata the curves are extrapolated using the estimates of 
##' cumulative incidence. 
##' 
##' @param x phreg object 
##' @param times possible times for which to report restricted mean 
##' @param covs possible covariate for Cox model 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' data(bmt); bmt$time <- bmt$time+runif(408)*0.001
##' out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)
##' 
##' rm1 <- resmean.phreg(out1,times=10*(1:6))
##' summary(rm1)
##' par(mfrow=c(1,2))
##' plot(rm1,se=1)
##' plot(rm1,years.lost=TRUE,se=1)
##' 
##' ## years.lost decomposed into causes
##' drm1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=10*(1:6))
##' summary(drm1)
##' @export
##' @aliases cif.yearslost  rmst.phreg
resmean.phreg <- function(x,times=NULL,covs=NULL,...) 
{# {{{
  ii <- invhess <- x$II

 if (!is.null(times)) {
	   tt <- expand.grid(times,0:(x$nstrata-1))
           mm <- cbind(x$jumptimes,x$strata.jumps,1)
           mm <- rbind(mm,cbind(tt[,1],tt[,2],0))
           ord <- order(mm[,1])
           mm <- mm[ord,]
	   phd <- which(mm[,3]==1)
	   newtimes <- which(mm[,3]==0)
	   S0i <- S0i2 <- rep(0,nrow(mm))
	   S0i[phd] <- c(1/x$S0)
	   S0ii2 <- c(1/(x$S0*(x$S0-1)))
	   S0ii2[x$S0==1] <- 0
	   S0i2[phd] <- S0ii2
           strata.jumps <- mm[,2]
	   nstrata <- x$nstrata; 
	   jumptimes <- mm[,1]
	   E <- matrix(0,nrow(mm),x$p)
	   E[phd,] <- x$E
    }   else { 
	   strata.jumps <- x$strata.jumps; 
	   nstrata <- x$nstrata; 
	   jumptimes <- x$jumptimes; 
	   S0i <- c(1/x$S0)
	   S0i2 <- c(1/(x$S0*(x$S0-1)))
	   S0i2[x$S0==1] <- 0
	   E <- x$E
 }

  cumhaz <- cumsumstrata(S0i,strata.jumps,nstrata)
  var.cumhazMG <- cumsumstrata(S0i2,strata.jumps,nstrata)

  ## baseline survival 
  km <- exp(cumsumstratasum(log(1-S0i),strata.jumps,nstrata,type="lagsum"))

  ## If covariates compute rr = exp( X beta)  and adjust km 
  if (!is.null(covs)) { rr <- exp(sum(covs * coef(x))) } else rr <- 1 
  km <- km^rr

  ## start integral in 0 
  dtime <- c(diffstrata(jumptimes,strata.jumps,nstrata))
  intkm <-  cumsumstrata( c(km)*dtime,strata.jumps,nstrata)

  ### variance of baseline term 
  var.intkmcumhaz <- cumsumstrata(intkm*S0i2,strata.jumps,nstrata)
  var.intkm2cumhaz <- cumsumstrata(intkm^2*S0i2,strata.jumps,nstrata)
  var.resmean <- intkm^2*var.cumhazMG+var.intkm2cumhaz-2*intkm*var.intkmcumhaz

   if (!is.null(x$coef)) {
      intp <- apply(E*S0i,2,cumsumstratasum,strata.jumps,nstrata,type="lagsum")
      intpS <-  apply(intp*c(km)*dtime,2,cumsumstrata,strata.jumps,nstrata)
      Dbeta <- -intpS 
      if (!is.null(covs)) {
	      intLam <- cumsumstratasum(S0i,strata.jumps,nstrata,type="lagsum")
	      intLamS <-  cumsumstrata( c(intLam)*km*dtime,strata.jumps,nstrata)
	      Dbeta <- Dbeta + matrix(covs,nrow=nrow(intLam),ncol=length(covs),byrow=TRUE)*c(intLamS)
      }
      varbetat <-   rowSums((Dbeta %*% x$II)*Dbeta)
      var.resmean <- rr^2*(var.resmean+varbetat)
  } else {
      varbetat <- 0
  }

  time<- jumptimes
  meanm <- cbind(time,intkm)
  se.resmean <- var.resmean^.5
  se.resmean[var.resmean<0] <- 0
  se.mm <- cbind(time,se.resmean)

  ### make output at specified times
  if (!is.null(times)) {
    intkmtimes <- meanm[newtimes,]
    se.intkmtimes <- se.resmean[newtimes]
    skmtimes <- mm[newtimes,2]
    years.lost <- intkmtimes[,1]-intkmtimes[,2]
    intkmtimes <- cbind(skmtimes,intkmtimes,se.intkmtimes,years.lost)
    colnames(intkmtimes) <- c("strata","times","rmean","se.rmean","years.lost")
    intkmtimes=data.frame(intkmtimes)
###    logintkmtimes=cbind(intkmtimes[,1:2],log(intkmtimes[,3]),se.intkmtimes/intkmtimes[,3])
###    colnames(logintkmtimes) <- c("strata","times","log-rmean","log-se.rmean")
  } else intkmtimes <- se.intkmtimes <- NULL


 out <- list(cumhaz=meanm,se.cumhaz=se.mm,covs=covs,
       time=time, strata=strata.jumps,nstrata=nstrata,
       jumps=1:length(km),strata.name=x$strata.name,
       strata.level=x$strata.level,
       intkmtimes=intkmtimes
       )
class(out) <- c("resmean_phreg")
return(out)
}# }}}

##' @export
rmst.phreg <- function(x,times=NULL,covs=NULL,...) 
{# {{{
out <- resmean.phreg(x,times=NULL,covs=NULL,...) 
return(out)
}# }}}

##' @export
cif.yearslost <- function(formula,data=data,cens.code=0,times=NULL,...)
{# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
###  if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
  if (!inherits(Y,c("Event","Surv"))) stop("Expected a 'Surv' or 'Event'-object")
  if (ncol(Y)==2) {
    exit <- Y[,1]
    entry <- NULL ## rep(0,nrow(Y))
    status <- Y[,2]
  } else {
    entry <- Y[,1]
    exit <- Y[,2]
    status <- Y[,3]
  }

  x <- phreg(formula,data=data,no.opt=TRUE,no.var=1,Z=as.matrix(status,ncol=1))
  causes <- sort(unique(x$cox.prep$Z[,1]))
  ccc <- which(causes %in% cens.code)
  if (length(ccc)>=1) causes <- causes[-ccc] 

 if (!is.null(times)) {# {{{
	   tt <- expand.grid(times,0:(x$nstrata-1))
###           cause.jumps <- x$U+x$E
           cause.jumps <- x$cox.prep$Z[x$cox.prep$jumps+1,1]
           mm <- cbind(x$jumptimes,x$strata.jumps,1,cause.jumps)
           mm <- rbind(mm,cbind(tt[,1],tt[,2],0,cens.code))
           ord <- order(mm[,1])
           mm <- mm[ord,]
	   phd <- which(mm[,3]==1)
	   newtimes <- which(mm[,3]==0)
	   S0i <- S0i2 <- rep(0,nrow(mm))
	   S0i[phd] <- c(1/x$S0)
	   S0ii2 <- c(1/(x$S0*(x$S0-1)))
	   S0ii2[x$S0==1] <- 0
	   S0i2[phd] <- S0ii2
           strata.jumps <- mm[,2]
	   nstrata <- x$nstrata; 
	   jumptimes <- mm[,1]
	   cause.jumps <- mm[,4]
    }   else { 
	   strata.jumps <- x$strata.jumps; 
	   nstrata <- x$nstrata; 
	   jumptimes <- x$jumptimes; 
           cause.jumps <- x$U+x$E
	   S0i <- c(1/x$S0)
	   S0i2 <- c(1/(x$S0*(x$S0-1)))
	   S0i2[x$S0==1] <- 0
 }# }}}

 ### formula from Pepe-Mori: SIM 93, 737-
 pepemori <- 0
 ### formula from PKA SIM 2013
 pka <- 0

  years.lostF1 <- se.years.lostF1m  <- se.years.lostF1 <- se.years.lostF1pm <- c()
  for (i in seq_along(causes)) {
	  jumpsi <- (cause.jumps==causes[i])*1
	  jumpsni <- (cause.jumps %in% causes[-i])*1
	  cumhaz <- cumsumstrata(S0i*jumpsi,strata.jumps,nstrata)
	  var.cumhazMG <- cumsumstrata(S0i2*jumpsi,strata.jumps,nstrata)
	  var.cumhazniMG <- cumsumstrata(S0i2*jumpsni,strata.jumps,nstrata)

	  ## baseline survival 
	  km <- exp(cumsumstratasum(log(1-S0i),strata.jumps,nstrata,type="lagsum"))
	  F1 <- cumsumstratasum(jumpsi*km*S0i,strata.jumps,nstrata,type="lagsum")
	  ### 1 - F2  = S + F1
	  F1n1 <- km+F1

	  ## start integral in 0 
	  dtime <- c(diffstrata(jumptimes,strata.jumps,nstrata))
	  intkm <-  cumsumstrata( c(km)*dtime,strata.jumps,nstrata)
	  intF1 <- cumsumstrata( c(F1)*dtime,strata.jumps,nstrata)

	  years.lostF1 <- cbind(years.lostF1,intF1)

	  if (pka==1) {
	  ### variance of baseline term  \int_0^t 1/Y^2 \int_s^t S(u) du  dN^i{{{
	  var.intkmcumhaz <- cumsumstrata(jumpsi*intkm*S0i2,strata.jumps,nstrata)
	  var.intkm2cumhaz <- cumsumstrata(jumpsi*intkm^2*S0i2,strata.jumps,nstrata)
	  var.resmean1 <- intkm^2*var.cumhazMG+var.intkm2cumhaz-2*intkm*var.intkmcumhaz

	  ### variance of baseline term  \int_0^t 1/Y^2 \int_s^t F_i(u) du  dN^-i
	  var.intF1cumhaz <- cumsumstrata(jumpsni*intF1*S0i2,strata.jumps,nstrata)
	  var.intF12cumhaz <- cumsumstrata(jumpsni*intF1^2*S0i2,strata.jumps,nstrata)
	  var.resmean2 <- intF1^2*var.cumhazniMG+var.intF12cumhaz-2*intF1*var.intF1cumhaz

	  var.intF1 <- (var.resmean1+var.resmean2)

          ## possibly negative due to rounding errors
	  var.intF1[var.intF1<0] <- 0
	  se.intF1 <- var.intF1^.5
	  se.years.lostF1 <- cbind(se.years.lostF1,se.intF1) # }}}
	  } 

          ### variance of baseline term  \int_0^t 1/Y^2 ( IF(t)-IF(s) - (t-s) (1-Fn1(s))) dN_1^i{{{
	  vars.1 <- cumsumstrata(jumpsi*intF1^2*S0i2,strata.jumps,nstrata)
	  vars.2 <- cumsumstrata(jumpsi*intF1*S0i2,strata.jumps,nstrata)
	  vars.11 <- intF1^2*var.cumhazMG+vars.1-2*intF1*vars.2

	  vars.3 <- cumsumstrata(jumpsi*F1n1^2*S0i2,strata.jumps,nstrata)
	  vars.4 <- cumsumstrata(jumpsi*F1n1^2*jumptimes^2*S0i2,strata.jumps,nstrata)
	  vars.5 <- cumsumstrata(jumpsi*F1n1^2*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.21 <- jumptimes^2*vars.3+vars.4-2*jumptimes*vars.5

          vars.6 <- cumsumstrata(jumpsi*F1n1*S0i2,strata.jumps,nstrata)
	  vars.7 <- cumsumstrata(jumpsi*F1n1*intF1*S0i2,strata.jumps,nstrata)
	  vars.8 <- cumsumstrata(jumpsi*F1n1*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.9 <- cumsumstrata(jumpsi*F1n1*intF1*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.31 <- jumptimes*intF1*vars.6+vars.9-intF1*vars.8-jumptimes*vars.7
 	  varsj1 <- vars.11+vars.21-2*vars.31

          ### variance of baseline term  \int_0^t 1/Y^2 ( IF(t)-IF(s) - (t-s) F1(s)) dN_n1^i
	  vars.1 <- cumsumstrata(jumpsni*intF1^2*S0i2,strata.jumps,nstrata)
	  vars.2 <- cumsumstrata(jumpsni*intF1*S0i2,strata.jumps,nstrata)
	  vars.11 <- intF1^2*var.cumhazniMG+vars.1-2*intF1*vars.2

	  vars.3 <- cumsumstrata(jumpsni*F1^2*S0i2,strata.jumps,nstrata)
	  vars.4 <- cumsumstrata(jumpsni*F1^2*jumptimes^2*S0i2,strata.jumps,nstrata)
	  vars.5 <- cumsumstrata(jumpsni*F1^2*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.21 <- jumptimes^2*vars.3+vars.4-2*jumptimes*vars.5

          vars.6 <- cumsumstrata(jumpsni*F1*S0i2,strata.jumps,nstrata)
	  vars.7 <- cumsumstrata(jumpsni*F1*intF1*S0i2,strata.jumps,nstrata)
	  vars.8 <- cumsumstrata(jumpsni*F1*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.9 <- cumsumstrata(jumpsni*F1*intF1*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.31 <- jumptimes*intF1*vars.6+vars.9-intF1*vars.8-jumptimes*vars.7
	  varsnj1 <- vars.11+vars.21-2*vars.31

	  varss <-  varsj1+varsnj1
          varss[varss<0] <- 0
	  se.intF1m <- varss^.5
	  se.years.lostF1m <- cbind(se.years.lostF1m,se.intF1m)
# }}}

	  if (pepemori==1) {
          ### variance Pepe-Mori {{{
          ### variance of dN_1 term  \int_0^t 1/Y^2   dN_n1^i
	  vars.1 <- cumsumstrata(jumpsi*intF1^2*S0i2,strata.jumps,nstrata)
	  vars.2 <- cumsumstrata(jumpsi*intF1*S0i2,strata.jumps,nstrata)
	  vars.11 <- intF1^2*var.cumhazMG+vars.1-2*intF1*vars.2

	  vars.4 <- cumsumstrata(jumpsi*jumptimes^2*S0i2,strata.jumps,nstrata)
	  vars.5 <- cumsumstrata(jumpsi*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.21 <- F1n1^2*(jumptimes^2*var.cumhazMG+vars.4-2*jumptimes*vars.5)

          vars.6 <- cumsumstrata(jumpsi*S0i2,strata.jumps,nstrata)
	  vars.7 <- cumsumstrata(jumpsi*intF1*S0i2,strata.jumps,nstrata)
	  vars.8 <- cumsumstrata(jumpsi*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.9 <- cumsumstrata(jumpsi*intF1*jumptimes*S0i2,strata.jumps,nstrata)
	  vars.31 <- F1n1*(jumptimes*intF1*vars.6+vars.9-intF1*vars.8-jumptimes*vars.7)
 	  varsj1 <- vars.11+vars.21-2*vars.31

          ### variance of baseline term  \int_0^t 1/Y^2 ( IF(t)-IF(s)) dN_n1^i
	  vars.1 <- cumsumstrata(jumpsni*intF1^2*S0i2,strata.jumps,nstrata)
	  vars.2 <- cumsumstrata(jumpsni*intF1*S0i2,strata.jumps,nstrata)
	  vars.11 <- intF1^2*var.cumhazniMG+vars.1-2*intF1*vars.2

	  varspm <- vars.11+varsj1

          varspm[varspm<0] <- 0
	  se.intF1pm <- varspm^.5
	  se.years.lostF1pm <- cbind(se.years.lostF1pm,se.intF1pm)
# }}}
	  }
 }
 years.lostF1 <- cbind(jumptimes,years.lostF1)
 colnames(years.lostF1) <- c("time",paste("intF_",causes,sep=""))


  ### make output at specified times
  if (!is.null(times)) {
    intF1times <- years.lostF1[newtimes,]
    years.lost <- apply(intF1times[,-1,drop=FALSE],1,sum)
    if (pka==1) se.intF1times <- se.years.lostF1[newtimes,-1]
    se.intF1mtimes <- se.years.lostF1m[newtimes,]
    stratatimes <- strata.jumps[newtimes]
    intF1mtimes <- cbind(stratatimes,intF1times,se.intF1mtimes,years.lost)
    colnames(intF1mtimes) <- c("strata","times",
			      paste0("intF1",causes),paste0("se.intF1",causes),"total-years-lost")
    intF1mtimes=data.frame(intF1mtimes)
    if (pka==1) { 
	    intF1pkatimes <- cbind(stratatimes,intF1times,se.intF1times)
            colnames(intF1pkatimes) <- c("strata","times",
			      paste0("intF1",causes),paste0("se.intF1",causes))
            intF1pkatimes=data.frame(intF1times)
    } else intF1pkatimes <- NULL 
    if (pepemori==1) {
          se.intF1pmtimes <- se.years.lostF1pm[newtimes,]
          intF1pmtimes <- cbind(stratatimes,intF1times,se.intF1pmtimes)
          colnames(intF1pmtimes) <- c("strata","times",
			      paste0("intF1",causes),paste0("se.intF1",causes))
         intF1pmtimes=data.frame(intF1pmtimes)
    } else intF1pmtimes <- NULL
  } else intF1pmtimes <-intF1mtimes <- intF1pkatimes <- NULL

 se.years.lostF1 <- cbind(jumptimes,se.years.lostF1m)
 colnames(se.years.lostF1) <- c("time",paste("intF_",causes,sep=""))


 ## name things to make use of other programs
 out <- list(cumhaz=years.lostF1,se.cumhaz=se.years.lostF1,
       time=jumptimes, strata=strata.jumps,nstrata=nstrata,
       jumps=1:length(km),strata.name=x$strata.name,
       strata.level=x$strata.level,
       intkmtimes=intF1mtimes, intF1times=intF1mtimes,
       intF1pmtimes=intF1pmtimes, intF1pkatimes=intF1pkatimes
       )
class(out) <- c("resmean_phreg")
return(out)
}# }}}

##' @export
summary.resmean_phreg <- function(object,...)
{# {{{
if (is.null(object$intkmtimes)) return(cbind(object$cumhaz,object$se.cumhaz[,2])) else return(object$intkmtimes)
}# }}}

##' @export
print.resmean_phreg <- function(x,...)
{# {{{
print(summary.resmean_phreg(x,...))
}# }}}

##' @export
plot.resmean_phreg <- function(x, se=FALSE,time=NULL,add=FALSE,ylim=NULL,xlim=NULL,
    lty=NULL,col=NULL,lwd=NULL,legend=TRUE,ylab=NULL,xlab=NULL,
    polygon=TRUE,level=0.95,stratas=NULL,robust=FALSE,years.lost=FALSE,cause=1,...) {# {{{

	if (inherits(x,"phreg") & is.null(ylab)) ylab <- "Cumulative hazard"
	if (inherits(x,"km") & is.null(ylab)) ylab <- "Survival probability"
	if (inherits(x,"cif") & is.null(ylab)) ylab <- "Probability"
	if (inherits(x,"resmean_phreg") & is.null(ylab)) ylab <- "Restricted residual mean life"
	if (years.lost) ylab <- "Years lost up to t: t - E(min(T,t))"
	if (is.null(xlab)) xlab <- "time"
   level <- -qnorm((1-level)/2)
   if (years.lost) rr <- range(x$cumhaz[,1]-x$cumhaz[,cause+1])  else rr <- range(x$cumhaz[,cause+1]) 
   strat <- x$strata[x$jumps]
   ylimo <- ylim
   if (is.null(ylim)) ylim <- rr
   if (is.null(xlim)) xlim <- range(x$cumhaz[,1])
   if (se==TRUE) {
      if (is.null(x$se.cumhaz) & is.null(x$robse.cumhaz) ) 
		   stop("phreg must be with cumhazard=TRUE\n"); 
   if (years.lost) rrse <- range(c(x$cumhaz[,1]-x$cumhaz[,cause+1]+level*x$se.cumhaz[,cause+1])) else 
       rrse <- range(c(x$cumhaz[,cause+1]+level*x$se.cumhaz[,cause+1])) 
       if (is.null(ylimo)) ylim <- rrse
   }

   ## all strata
   if (is.null(stratas)) stratas <- 0:(x$nstrata-1) 
   ltys <- lty
   cols <- col
   lwds <- lwd

   if (length(stratas)>0 & x$nstrata>1) { ## with strata
   lstrata <- x$strata.level[(stratas+1)]
   stratn <-  substring(x$strata.name,8,nchar(x$strata.name)-1)
   stratnames <- paste(stratn,lstrata,sep=":")
   
      if (!is.matrix(lty)) {
         if (is.null(lty)) ltys <- 1:length(stratas) else if (length(lty)!=length(stratas)) ltys <- rep(lty[1],length(stratas))
      } else ltys <- lty
      if (!is.matrix(col)) {
         if (is.null(col)) cols <- 1:length(stratas) else 
		 if (length(col)!=length(stratas)) cols <- rep(col[1],length(stratas))
      } else cols <- col
      if (!is.matrix(lwd)) {
         if (is.null(lwd)) lwds <- rep(1,length(stratas)) else 
		 if (length(lwd)!=length(stratas)) lwds <- rep(lwd[1],length(stratas))
      } else lwds <- lwd
   } else { 
     stratnames <- "Baseline" 
     if (is.matrix(col))  cols <- col
     if (is.null(col)) cols <- 1  else cols <- col[1]
     if (is.matrix(lty))  ltys <- lty
     if (is.null(lty)) ltys <- 1  else ltys <- lty[1]
     if (is.matrix(lwd))  lwds <- lwd
     if (is.null(lwd)) lwds <- 1  else lwds <- lwd[1]
   }

  if (!is.matrix(ltys))  ltys <- cbind(ltys,ltys,ltys)
  if (!is.matrix(cols))  cols <- cbind(cols,cols,cols)
  if (!is.matrix(lwds))  lwds <- cbind(lwds,lwds,lwds)

  if (years.lost) x$cumhaz[,cause+1] <- x$cumhaz[,1]-x$cumhaz[,cause+1]
  first <- 0
  for (i in seq(stratas)) {
      j <- stratas[i]
        cumhazard <- x$cumhaz[strat==j,c(1,cause+1),drop=FALSE]
        if (!is.null(cumhazard)) {
	if (nrow(cumhazard)>1) {
        if (add | first==1) 
        lines(cumhazard,type="l",lty=ltys[i,1],col=cols[i,1],lwd=lwds[i,1])   
       else {
	  first <- 1
          plot(cumhazard,type="l",lty=ltys[i,1],col=cols[i,1],lwd=lwds[i,1],ylim=ylim,ylab=ylab,xlab=xlab,xlim=xlim,...)
       }
       if (se==TRUE) {
	    if (robust==TRUE) secumhazard  <- x$robse.cumhaz[strat==j,c(1,cause+1),drop=FALSE]
	    else secumhazard <- x$se.cumhaz[strat==j,c(1,cause+1),drop=FALSE]
		 ul <-cbind(cumhazard[,1],cumhazard[,2]+level*secumhazard[,2])
		 nl <-cbind(cumhazard[,1],cumhazard[,2]-level*secumhazard[,2])
		 if (inherits(x,"km")) { ul[,2] <- x$upper[x$strata==j]; 
		                          nl[,2] <- x$lower[x$strata==j];
		                          wna <- which(is.na(ul[,2]))
		                          ul[wna,2] <- 0
		                          nl[wna,2] <- 0

		 }
	      if (!polygon) {
		  lines(nl,type="l",lty=ltys[i,2],col=cols[i,2],lwd=lwds[i,2])
		  lines(ul,type="l",lty=ltys[i,3],col=cols[i,3],lwd=lwds[i,3])
	      } else {
                 ## type="l" confidence regions
		 ll <- length(nl[,1])
		 timess <- nl[,1]
		 tt <- c(timess,rev(timess))
		 yy <- c(nl[,2],rev(ul[,2]))
		 col.alpha<-0.1
		 col.ci<-cols[i,]
		 col.trans <- sapply(col.ci, FUN=function(x) 
		   do.call(grDevices::rgb,as.list(c(grDevices::col2rgb(x)/255,col.alpha))))
		 polygon(tt,yy,lty=0,col=col.trans)
	      }
        }
     }
     }
    }

    where <- "topleft"; 
    if (inherits(x,"km")) where <-  "topright"
    if (legend & (!add)) 
    graphics::legend(where,legend=stratnames,col=cols[,1],lty=ltys[,1])

}# }}}

# }}}

##' G-estimator for Cox and Fine-Gray model 
##'
##' Computes G-estimator \deqn{ \hat S(t,A=a) = n^{-1} \sum_i \hat S(t,A=a,Z_i) }
##' for the Cox model based on phreg og the Fine-Gray model based on the
##' cifreg function. Gives influence functions of these risk estimates and SE's are 
##' based on  these.  If first covariate is a factor then all contrast are computed, 
##' and if continuous then considered covariate values are given by Avalues.
##'
##' @param x phreg or cifreg object
##' @param data data frame for risk averaging
##' @param time for estimate
##' @param Avalues values to compare for first covariate A
##' @param varname if given then averages for this variable, default is first variable
##' @author Thomas Scheike
##' @examples
##' 
##' data(bmt); bmt$time <- bmt$time+runif(408)*0.001
##' bmt$event <- (bmt$cause!=0)*1
##' dfactor(bmt) <- tcell.f~tcell
##'
##' fg1 <- cifreg(Event(time,cause)~tcell.f+platelet+age,bmt,cause=1,
##'               cox.prep=TRUE,propodds=NULL)
##' summary(survivalG(fg1,bmt,50))
##'
##' ss <- phreg(Surv(time,event)~tcell.f+platelet+age,bmt) 
##' summary(survivalG(ss,bmt,50))
##' @export
survivalG <- function(x,data,time=NULL,Avalues=c(0,1),varname=NULL)
{# {{{

if (is.null(time)) stop("Give time for estimation of survival/cumulative incidence\n")

if (inherits(x,"cifreg"))
Aiid <- IIDbaseline.cifreg(x,time=time) else Aiid <- IIDbaseline.phreg(x,time=time)

### dealing with first variable that is a factor 

if (is.null(varname))  {
treat.name <- all.vars(update.formula(x$formula,1~.))[1]
} else treat.name <- varname
treatvar <- data[,treat.name]

if (is.factor(treatvar)) {
###if (!is.factor(treatvar)) stop(paste("treatment=",treat.name,"must be coded as factor \n",sep="")); 
## treatvar, 1,2,...,nlev or 1,2
nlev <- nlevels(treatvar)
nlevs <- levels(treatvar)
###treatvar <- as.numeric(treatvar)
ntreatvar <- as.numeric(treatvar)
ytreat <- ntreatvar-1
} else {
   nlevs <- Avalues
}

formulaX <- update.formula(x$formula,.~.)
formulaX <- drop.specials(formulaX,"cluster")
formulaX[-2]
datA <- dkeep(data,x=all.vars(formulaX))
xlev <- lapply(datA,levels)

###
###Gf <- function(p,ic=0) {# {{{
###   risks <- c()
###   a <- nlevs[1]
###   for (a in nlevs) {# {{{
###      datA[,treat.name] <- a
###      Xa <- model.matrix(formulaX[-2],datA,xlev=xlev)[,-1]
###      rra <- exp(Xa %*% p[-1])
###
###      if (inherits(x,"cifreg")) { ps0 <- 1- exp(-p[1]*rra) } 
###      else if (inherits(x,"recreg")) { ps0 <- p[1]*rra }
###      else { ps0 <-  exp(-p[1]*rra) } 
###      risks <- cbind(risks,ps0)
###    }# }}}
###
###Gest <- apply(risks,2,mean)
###if (ic==1) Gest <- list(Gest=Gest,iid=t(t(risks)-Gest))
###return(Gest)
###}
#### }}}
###

cumhaz.time <- Cpred(x$cumhaz,time)[-1]
k <- 1; risks <- c(); DariskG <- list()
for (a in nlevs) { ## {{{
 datA[,treat.name] <- a
 Xa <- model.matrix(formulaX[-2],datA,xlev=xlev)[,-1]
 rra <- c(exp(Xa %*% x$coef))
 if (inherits(x,"phreg"))  { ps0 <- exp(-cumhaz.time*rra); Dma  <- -cbind(rra*ps0,ps0*cumhaz.time*rra*Xa);    }
 if (inherits(x,"cifreg")) { ps0 <- 1- exp(-cumhaz.time*rra); Dma  <- cbind(rra*(1-ps0),(1-ps0)*cumhaz.time*rra*Xa);    }
  else if (inherits(x,"recreg")) { ps0 <- cumhaz.time*rra ; Dma <-  cbind(rra,cumhaz.time*rra*Xa) }
 risks <- cbind(risks,ps0)
 DariskG[[k]] <- apply(Dma,2,sum)
 k <- k+1
} ## }}}

Grisk <- apply(risks,2,mean)
risk.iid  <- t(t(risks)-Grisk)
###icf <- Gf(theta,ic=1)
######
###DG <- numDeriv::jacobian(Gf,theta,ic=0)
nid <- max(x$id)
risk.iid <- apply(risk.iid,2,sumstrata,x$id-1,nid)/nid 
for (a in seq_along(nlevs)) risk.iid[,a] <- risk.iid[,a]+ cbind(Aiid$base.iid,Aiid$beta.iid)%*% DariskG[[a]]/nid
vv <- crossprod(risk.iid)

###estimate(lava::estimate(coef=theta,vcov=vv,f=function(p) Gf(p,ic=0))
out <- estimate(coef=Grisk,vcov=vv,labels=paste("risk",nlevs,sep=""))
ed <- estimate(coef=Grisk,vcov=vv,out,function(p) p[-1]-p[1])
rd <- estimate(coef=Grisk,vcov=vv,out,function(p) p[-1]/p[1],null=1)
out <- list(risk.iid=risk.iid,risk=out,difference=ed,ratio=rd,vcov=vv)
class(out) <- "survivalG"
return(out)
} ## }}}


###{{{ summary 

##' @export
summary.survivalG <- function(object,...) {
  res <- list(risk=object$risk,difference=object$difference,ratio=object$ratio)
  class(res) <- "summary.survivalG"
  res
}

##' @export
print.summary.survivalG  <- function(x,...) {
    cat("risk:\n")
    print(x$risk,...)
    cat("\n")

    cat("Average Treatment effects (G-estimator) :\n")
    print(x$difference,...)
    cat("\n")

    cat("Average Treatment effect ratio (G-estimator) :\n")
    print(x$ratio$coefmat,...)
    cat("\n")

}

###}}} summary 
# }}}



##' Fast additive hazards model with robust standard errors 
##'
##' Fast Lin-Ying additive hazards model with a possibly stratified baseline. 
##' Robust variance is default variance with the summary. 
##'
##' influence functions (iid) will follow numerical order of given cluster variable
##' so ordering after $id will give iid in order of data-set.
##'
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param ... Additional arguments to phreg 
##' @author Thomas Scheike
##' @examples
##' 
##' data(bmt); bmt$time <- bmt$time+runif(408)*0.001
##' out <- aalenMets(Surv(time,cause==1)~tcell+platelet+age,data=bmt)
##' summary(out)
##' 
##' ## out2 <- timereg::aalen(Surv(time,cause==1)~const(tcell)+const(platelet)+const(age),data=bmt)
##' ## summary(out2)
##' 
##' @export
aalenMets <- function(formula,data=data,...)
{# {{{
formula.call <- formula
x <- phreg(formula,data=data,no.opt=TRUE,...)

xx <- x$cox.prep
###ii <- solve(x$hessian)
###S0i <- rep(0,length(xx$strata))
###S0i[xx$jumps+1] <- 1/x$S0
###Z <- xx$X
###U <- E <- matrix(0,nrow(xx$X),x$p)
###E[xx$jumps+1,] <- x$E
###U[xx$jumps+1,] <- x$U
###cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
###EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)

### computation of intZHZ matrix, and gamma and baseline ##################
# {{{
eXb <- c(xx$sign) * c(xx$weights)
S0 = c(revcumsumstrata(eXb,xx$strata,xx$nstrata))
E=apply(eXb*xx$X,2,revcumsumstrata,xx$strata,xx$nstrata)/S0; 
S2=apply(eXb*xx$XX,2,revcumsumstrata,xx$strata,xx$nstrata)
###
E2=.Call("vecMatMat",E,E)$vXZ;  
###
dts <- c(diffstrata(xx$time,xx$strata,xx$nstrata))
dt <- diff(c(0,xx$time))
###
intZbar <- apply(E*dts,2,cumsumstrata,xx$strata,xx$nstrata) 
intZHZt <- apply((S2-S0*E2)*dts,2,cumsum) 
intZHZ <- matrix(tail(intZHZt,1),ncol(E), ncol(E))
IintZHZ  <-  solve(intZHZ)
intZHdN <- matrix(x$gradient,ncol(E),1)
gamma <- IintZHZ %*% intZHdN
rownames(gamma)  <-  colnames(x$X)
x$cumhaz[,2] <- x$cumhaz[,2]- intZbar[xx$jumps+1,] %*% gamma
# }}}

### iid gamma #########################################
# {{{
mm <-  xx$X * c(xx$X %*% gamma) * c(xx$time)+ apply( E* c(E %*% gamma)*dts,2,cumsumstrata,xx$strata,xx$nstrata) -  xx$X* c(intZbar %*% gamma) - c(xx$X %*% gamma)* intZbar 
mm <- c(xx$weights) * mm

id <- xx$id
mm <- apply(mm,2,sumstrata,id,max(id)+1)
MGt <- t(x$hessian %*% t(iid(x))) + mm
iid <-  MGt %*% IintZHZ
# }}}

coef <- c(gamma)
names(coef) <- rownames(gamma)
x$coef <- coef
x$var <- crossprod(iid) 
x$iid <- iid
x$intZHdN <- intZHdN
x$intZHZ  <-  intZHZ
x$formula <- formula.call
x$ihessian <- IintZHZ
## score of gamma 
x$gradient <- c(intZHdN-intZHZ %*% gamma)

x$no.opt <- FALSE
class(x) <- c(class(x),"aalenMets")

return(x)
}# }}}


##' Kaplan-Meier with robust standard errors 
##'
##' Kaplan-Meier with robust standard errors 
##' Robust variance is default variance with the summary. 
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param conf.type transformation 
##' @param conf.int level of confidence intervals 
##' @param robust for robust standard errors based on martingales 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' data(TRACE)
##' TRACE$cluster <- sample(1:100,1878,replace=TRUE)
##' out1 <- km(Surv(time,status==9)~strata(vf,chf),data=TRACE)
##' out2 <- km(Surv(time,status==9)~strata(vf,chf)+cluster(cluster),data=TRACE)
##' 
##' par(mfrow=c(1,2))
##' bplot(out1,se=TRUE)
##' bplot(out2,se=TRUE)
##' @export
km <- function(formula,data=data,conf.type="log",conf.int=0.95,robust=TRUE,...)
{# {{{
 coxo <- phreg(formula,data=data)
 coxo <- robust.phreg(coxo)

 chaz <-     coxo$cumhaz[,2]
 time <-     coxo$cumhaz[,1]
 if (robust) std.err <-  coxo$robse.cumhaz[,2]
 else std.err <-  coxo$se.cumhaz[,2]
 strat <-    coxo$strata[coxo$jumps]

 S0i  <-  1/coxo$S0
 kmt <- exp(cumsumstrata(log(1-S0i),strat,coxo$nstrata))
 temp <- list(surv=kmt)

 zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)

  ### different conf-types
 if (conf.type == "plain") {# {{{
    temp1 <- temp$surv + zval * std.err * temp$surv
    temp2 <- temp$surv - zval * std.err * temp$surv
    temp <- c(temp, list(upper = pmin(temp1, 1), lower = pmax(temp2,
	0), conf.type = "plain", conf.int = conf.int))
 }
 if (conf.type == "log") {
    xx <- ifelse(temp$surv == 0, 1, temp$surv)
    temp1 <- ifelse(temp$surv == 0, NA, exp(log(xx) + zval * std.err))
    temp2 <- ifelse(temp$surv == 0, NA, exp(log(xx) - zval * std.err))
    temp <- c(temp, list(upper = pmin(temp1, 1), lower = temp2,
	conf.type = "log", conf.int = conf.int))
 }
 if (conf.type == "log-log") {
    who <- (temp$surv == 0 | temp$surv == 1)
    temp3 <- ifelse(temp$surv == 0, NA, 1)
    xx <- ifelse(who, 0.1, temp$surv)
    temp1 <- exp(-exp(log(-log(xx)) + zval * std.err/log(xx)))
    temp1 <- ifelse(who, temp3, temp1)
    temp2 <- exp(-exp(log(-log(xx)) - zval * std.err/log(xx)))
    temp2 <- ifelse(who, temp3, temp2)
    temp <- c(temp, list(upper = temp1, lower = temp2,
	conf.type = "log-log", conf.int = conf.int))
 }# }}}

 ### to use basehazplot.phreg
 temp <- c(temp,
 list(cumhaz=cbind(time,kmt),se.cumhaz=cbind(time,kmt*std.err),time=time,
	strata=strat,nstrata=coxo$nstrata,
        jumps=1:length(kmt),strata.name=coxo$strata.name,strata.level=coxo$strata.level))
 class(temp) <- c("km","phreg")
 return(temp)
}# }}}

##' Cumulative incidence with robust standard errors 
##'
##' Cumulative incidence with robust standard errors 
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param cause NULL looks at all, otherwise specify which cause to consider
##' @param cens.code censoring code "0" is default
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' data(TRACE)
##' TRACE$cluster <- sample(1:100,1878,replace=TRUE)
##' out1 <- cif(Event(time,status)~+1,data=TRACE,cause=9)
##' out2 <- cif(Event(time,status)~+1+cluster(cluster),data=TRACE,cause=9)
##' 
##' out1 <- cif(Event(time,status)~strata(vf,chf),data=TRACE,cause=9)
##' out2 <- cif(Event(time,status)~strata(vf,chf)+cluster(cluster),data=TRACE,cause=9)
##' 
##' par(mfrow=c(1,2))
##' bplot(out1,se=TRUE)
##' bplot(out2,se=TRUE)
##' @export
cif <- function(formula,data=data,cause=1,cens.code=0,...)
{# {{{

  cl <- match.call()
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
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else strata.name <- NULL
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms  <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  id.orig <- id; 
  if (!is.null(id)) {
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(exit))-1; 


  statusE <- 1*(status==cause)
  statusD <- 1*(status!=cens.code)
  if (ncol(Y)==3) {
	  if (!is.null(strata)) {
  formE <- as.formula(paste("Surv(entry,exit,statusE)~strata(strata)+cluster(id_1_)",sep=""))
  formD <- as.formula(paste("Surv(entry,exit,statusD)~strata(strata)+cluster(id_1_)",sep=""))
	  } else {
  formE <- as.formula(paste("Surv(entry,exit,statusE)~1+cluster(id_1_)",sep=""))
  formD <- as.formula(paste("Surv(entry,exit,statusD)~1+cluster(id_1_)",sep=""))

	  }
  } else {
	  if (!is.null(strata)) {
  formE <- as.formula(paste("Surv(exit,statusE)~strata(strata)+cluster(id_1_)",sep=""))
  formD <- as.formula(paste("Surv(exit,statusD)~strata(strata)+cluster(id_1_)",sep=""))
	  } else {
  formE <- as.formula(paste("Surv(exit,statusE)~cluster(id_1_)",sep=""))
  formD <- as.formula(paste("Surv(exit,statusD)~cluster(id_1_)",sep=""))
	  } 
  }

  data$id_1_ <- id

  if (sum(statusE)==0) warning("No events of type 1\n"); 

  coxE <- phreg(formE,data=data,...)
  coxS <- phreg(formD,data=data,...)

  ### cif 
  cifo <- recurrentMarginal(coxE,coxS)

  ### to use basehazplot.phreg
  class(cifo) <- c("cif","phreg")
  return(cifo)
}# }}}


##' Proportional odds survival model
##'
##' Semiparametric Proportional odds model, that has the advantage that 
##' \deqn{
##' logit(S(t|x)) = \log(\Lambda(t)) + x \beta
##' }
##' so covariate effects give OR of survival. 
##'
##' This is equivalent to using a hazards model 
##' \deqn{
##'   Z \lambda(t) \exp(x \beta)
##' }
##' where Z is gamma distributed with mean and variance 1.
##'
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param offset offsets for exp(x beta) terms 
##' @param weights weights for score equations
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @references
##'
##' The proportional odds cumulative incidence model for competing risks,
##' Eriksson, Frank and Li, Jianing and Scheike, Thomas and Zhang, Mei-Jie,
##' Biometrics, 2015, 3, 687--695, 71,
##'
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- logitSurv(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' summary(out1)
##' gof(out1)
##' plot(out1)
##' @export
logitSurv <- function(formula,data,offset=NULL,weights=NULL,...)
{# {{{
   out <- phreg(formula,data,offset=offset,weights=weights,propodds=1,...)
   return(out)
}# }}}

###{{{ predict with se for baseline

predictPhreg <- function(x,jumptimes,S0,beta,time=NULL,X=NULL,surv=FALSE,band=FALSE,...) {
    strata <- x$strata[x$jumps]# {{{
    nstrata <- x$nstrata
    
    ## Brewslow estimator
    if (is.null(x$cumhaz)) {
        ##II <- x$II
        ##x$jumptimes
        II <- -solve(x$hessian)
        chaz <- cbind(jumptimes,cumsumstrata(1/S0,strata,nstrata))
        DLambeta.t <- apply(x$E/c(x$S0),2,cumsumstrata,strata,nstrata)
        varbetat <- apply((DLambeta.t %*%  II)*DLambeta.t,1,sum)
        se.chaz <- cbind(jumptimes,(cumsumstrata(1/S0^2,strata,nstrata)+varbetat)^.5)
    } else {
        chaz <- x$cumhaz
        se.chaz <- x$se.cumhaz
    }

    if (!is.null(time)) {
	### do within strata
        chaz <- cpred(chaz,time)
        se.chaz <- cpred(se.chaz,time)
    }
    colnames(chaz) <- c("time","chaz")
    colnames(se.chaz) <- c("time","se.chaz")

    if (band==TRUE) { ## on log-scale  for one strata# {{{
      ii <- -solve(x$hessian)
      Ubeta <- x$U
      betaiid <- t(ii %*% t(Ubeta))
      cumhaz <-  x$cumhaz[,1,drop=FALSE]
      se.chaz <- x$se.cumhaz[,1]
      ###
      rr <- c( exp(sum(c(X) %*% x$coef)))
      Pt      <- outer(cumhaz[,2],c(X))
      DLambeta.t <- apply(x$E/c(x$S0),2,cumsumstrata,strata,nstrata)
      Pt <- DLambeta.t  - Pt
      ### se of cumulaive hazard for this covariate , can use different versions of variance for beta
      varbetat <- rowSums((Pt %*% ii)*Pt)
      se.chazexb <- cbind(jumptimes,rr*(cumsumstrata(1/S0^2,strata,nstrata)+varbetat)^.5)
    }# }}}

    if (!is.null(X)) {
      H <- exp(X%*%beta)
      if (nrow(chaz)==length(H)) {
        chaz[,2] <- chaz[,2]*H
      } else {
        chaz2 <- c()
        X <- rbind(X)
        for (i in seq(nrow(X)))
          chaz2 <- rbind(chaz2,
                         cbind(chaz[,1],chaz[,2]*H[i],
                               rep(1,nrow(chaz))%x%X[i,,drop=FALSE]))
        chaz <- chaz2;
        nn <- c("time","chaz",names(beta))
        colnames(chaz) <- nn
      }
    }
    if (surv) {    
      chaz[,2] <- exp(-chaz[,2])
      colnames(chaz)[2] <- "surv"
    }
    return(chaz)
}# }}}

##' Predictions from proportional hazards model
##'
##' @param object phreg object
##' @param newdata data.frame
##' @param times Time where to predict variable, default is all time-points from the object sorted
##' @param individual.time to use one (individual) time per subject, and then newdata and times have same length and makes only predictions for these individual times.
##' @param tminus to make predictions in T- that is strictly before given times, useful for IPCW techniques
##' @param se with standard errors and upper and lower confidence intervals.
##' @param robust to get robust se's. 
##' @param conf.type transformation for suvival estimates, default is log
##' @param conf.int significance level
##' @param km to use Kaplan-Meier product-limit for baseline \deqn{S_{s0}(t)= (1 - dA_{s0}(t))}, otherwise take exp of cumulative baseline.
##' @param ... Additional arguments to plot functions
##' @aliases headstrata tailstrata revcumsumstrata revcumsumstratasum cumsumstrata sumstrata covfr covfridstrata covfridstrataCov cumsumidstratasum cumsumidstratasumCov cumsumstratasum revcumsum revcumsumidstratasum revcumsumidstratasumCov robust.basehaz.phreg matdoubleindex mdi cumsum2strata revcumsum2strata revcumsum2stratafdN
##' @export
predict.phreg <- function(object,newdata,times=NULL,individual.time=FALSE,tminus=FALSE,se=TRUE,robust=FALSE,conf.type="log",conf.int=0.95,km=FALSE,...) 
{# {{{ default is all time-points from the object

   ### take baseline and strata from object# {{{
   jumptimes <- object$cumhaz[,1]
   chaz <- object$cumhaz[,2]
   if (is.null(object$nstrata)) {  ## try to make more robust
	nstrata <- 1; 
	strata <- rep(1,length(jumptimes))
   } else {
      nstrata <- object$nstrata
      strata <- object$strata[object$jumps]
   }
   if (length(jumptimes)==0) se <- FALSE
   if (se) {
   if (!robust) { 
	   se.chaz <- object$se.cumhaz[,2] 
	   varbeta <- object$ihessian  
	   if (!object$no.opt) Pt <- apply(object$E/c(object$S0),2,cumsumstrata,strata,nstrata)
           else Pt <- 0
   } else {
          if ((object$no.opt) | is.null(object$coef)) fixbeta<- 1 else fixbeta <- 0
          IsdM <- squareintHdM(object,ft=NULL,fixbeta=fixbeta,...)
          ###
          se.chaz <-   IsdM$varInt[object$jumps]^.5
	  covv <- IsdM$covv[object$jumps,,drop=FALSE]
	  varbeta <- IsdM$vbeta
	  Pt <- IsdM$Ht[object$jumps,,drop=FALSE]
   }
   } # }}}

   
### setting up newdata with factors and strata 
desX <- readPhreg(object,newdata) 
X <- desX$X
strataNew <- desX$strata

if (is.null(times)) times <- sort(unique(c(object$exit)))
if (individual.time & is.null(times)) times <- c(object$exit)
if (individual.time & length(times)==1) times <- rep(times,length(object$exit)) 

    se.cumhaz <- NULL
    if (!individual.time) {
       pcumhaz <- surv <- matrix(0,nrow(X),length(times))
       if (se) se.cumhaz <- matrix(0,nrow(X),length(times))
    } else { 
        pcumhaz <- surv <- matrix(0,nrow(X),1)
        if (se) se.cumhaz <- matrix(0,nrow(X),1)
    }
    hazt <- length(times)

    for (j in unique(strataNew)) {
###        where <- sindex.prodlim(c(0,jumptimes[strata==j]),times,strict=tminus)
       where <- predictCumhaz(c(0,jumptimes[strata==j]),times,type="left",tminus=tminus)

###	if (sum(abs(whereO-where))>=1) {
###	print(c(0,jumptimes[strata==j]))
###	print(times)
###	print(cbind(where,whereO,where-whereO,times,c(0,jumptimes[strata==j])[where],c(0,jumptimes[strata==j])[whereO]))
###	print("sindex.prodlim - phreg-predict"); 
###	}
	plhazt <- hazt <- c(0,chaz[strata==j])
	if (km) { plhazt <- suppressWarnings(c(1,exp(cumsum(log(1-diff(hazt))))));  plhazt[is.na(hazt)] <- 0 }
	if (se) se.hazt <- c(0,se.chaz[strata==j])
	Xs <- X[strataNew==j,,drop=FALSE]
        if (object$p==0) RR <- rep(1,nrow(Xs)) else RR <- c(exp( Xs %*% coef(object)))
	if (se)  { # {{{ based on Hazard's 
		if (object$p>0) {
			Ps <- Pt[strata==j,,drop=FALSE]
			Ps <- rbind(0,Ps)[where,,drop=FALSE]
			Xbeta <- Xs %*% varbeta
			seXbeta <- rowSums(Xbeta*Xs)^.5
			cov2 <- cov1 <- Xbeta %*% t(Ps*hazt[where])
		        if (robust)	{
			   covvs <- covv[strata==j,,drop=FALSE]
			   covvs <- rbind(0,covvs)[where,,drop=FALSE]
                           covv1 <- Xs %*% t((covvs*hazt[where]))
			   cov1 <- cov1-covv1
			}
		} else cov1 <- 0 
	}# }}}
   haztw <- hazt[where] 
   if (se) se.haztw <- se.hazt[where] 
	if (is.null(object$propodds)) {
           plhaztw <- plhazt[where] 
 	   if (!individual.time) pcumhaz[strataNew==j,]  <- RR%o%haztw else pcumhaz[strataNew==j,] <- RR*haztw[strataNew==j]
           if (!km) {
	     if (!individual.time) surv[strataNew==j,]  <- exp(- RR%o%haztw)
	     else surv[strataNew==j,]  <- exp(-RR*haztw[strataNew==j])
	   } else {
             if (!individual.time) surv[strataNew==j,]  <- exp( RR %o% log(plhaztw))
	     else surv[strataNew==j,]  <- plhaztw[strataNew==j]^RR
	   }
	} else {
	  if (!individual.time) surv[strataNew==j,]  <- 1/(1+RR%o%haztw)
          else surv[strataNew==j,]  <- 1/(1+RR*haztw[strataNew==j])
	}
	if (se) {# {{{
	    if (object$p>0)  {
	       if (!individual.time) se.cumhaz[strataNew==j,]  <- 
		     ((RR %o% se.haztw)^2+(c(RR*seXbeta) %o% haztw)^2-2*RR^2*cov1)^.5
	        else se.cumhaz[strataNew==j,]  <- RR* (se.haztw^2+(c(seXbeta)*haztw)^2-2*diag(cov1))^.5
	    } else {
	       if (!individual.time) se.cumhaz[strataNew==j,]  <- RR %o% (se.haztw)
	        else se.cumhaz[strataNew==j,]  <- RR* se.haztw[strataNew==j]  
	    }
	}# }}}
    }


    zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)
    std.err <- se.cumhaz
    cisurv  <- list()
    cisurv$upper <- NULL
    cisurv$lower <- NULL

### different conf-types for surv
 if (se) {# {{{
 if (conf.type == "plain") {# {{{
    temp1 <- surv + zval * std.err * surv
    temp2 <- surv - zval * std.err * surv
    cisurv <- list(upper = pmin(temp1, 1), lower = pmax(temp2,
	0), conf.type = "plain", conf.int = conf.int)
 }
 if (conf.type == "log") {
    xx <- ifelse(surv == 0, 1, surv)
    temp1 <- ifelse(surv == 0, NA, exp(log(xx) + zval * std.err))
    temp2 <- ifelse(surv == 0, NA, exp(log(xx) - zval * std.err))
    cisurv <- list(upper = pmin(temp1, 1), lower = temp2,
	conf.type = "log", conf.int = conf.int)
 }
 if (conf.type == "log-log") {
    who <- (surv == 0 | surv == 1)
    temp3 <- ifelse(surv == 0, NA, 1)
    xx <- ifelse(who, 0.1, surv)
    temp1 <- exp(-exp(log(-log(xx)) + zval * std.err/log(xx)))
    temp1 <- ifelse(who, temp3, temp1)
    temp2 <- exp(-exp(log(-log(xx)) - zval * std.err/log(xx)))
    temp2 <- ifelse(who, temp3, temp2)
    cisurv <- list(upper = temp1, lower = temp2,
	conf.type = "log-log", conf.int = conf.int)
 }# }}}
 }# }}}

 if (object$p>0) RR <-  exp(X %*% coef(object)) else RR <- rep(1,nrow(X))

### non-cox setting
if (!is.null(object$propodds)) pcumhaz <- -log(surv)

 out <- list(surv=surv,times=times,
	      surv.upper=cisurv$upper,surv.lower=cisurv$lower,cumhaz=pcumhaz,se.cumhaz=se.cumhaz,strata=strataNew,X=X, RR=RR)
 if (length(class(object))==2 && substr(class(object)[2],1,3)=="cif") {
	 out <- c(out,list(cif=1-out$surv,cif.lower=1-out$surv.upper, cif.upper=1-out$surv.lower))
 }
 class(out) <- c("predictphreg")
 if (length(class(object))==2) class(out) <- c("predictphreg",class(object)[2])
 return(out)
}# }}}


##' @export
print.predictphreg  <- function(x,se=FALSE,...) {# {{{

   if (is.null(x$se.cumhaz) & se==TRUE)  {
       warning("predict.phreg must be with se=TRUE\n"); 
       se <- FALSE
   }
  type <- "surv"
  if ((length(class(x))==2) && (substr(class(x)[2],1,3)=="cif")) type <- "cif"

  if (type[1]=="surv") { 
	  xx <- x$surv 
  } else if (type[1]=="cif") {
     xx <- x$cif 
  } else { xx <- x$cumhaz; names(xx) <- "cumhaz"}
  ## colnames(xx) <- type

  if (se) {
     if (type[1]=="surv") {
	     upper <- x$surv.upper; lower <- x$surv.lower
     } else if (type[1]=="cif") {
	     lower <- x$cif.lower; upper <- x$cif.upper
     } else { upper <- NA; lower <- NA;}
     ci <- cbind(lower,upper)
     ## colnames(ci) <- c("lower","upper")
     xx <- cbind(xx,ci)
  }

  print(xx)
  invisible(xx)
}# }}}


##' @export
plot.predictphreg  <- function(x,se=FALSE,add=FALSE,ylim=NULL,xlim=NULL,lty=NULL,col=NULL,type=c("surv","cumhaz","cif"),ylab=NULL,xlab=NULL,
    polygon=TRUE,level=0.95,whichx=NULL,robust=FALSE,...) {# {{{
   if (type[1]=="surv" & is.null(ylab)) ylab <- "Survival probability"
   if (type[1]=="cif" & is.null(ylab)) ylab <- "Cumulative probability"
   if (type[1]=="surv" & length(class(x))==2) ylab <- "Cumulative probability"
   if (type[1]=="cumhaz" & is.null(ylab)) ylab <- "Cumulative  hazard"
   if (is.null(xlab)) xlab <- "time"
   level <- -qnorm((1-level)/2)
   if (type[1]=="surv") rr <- c(0,1) 
   if (type[1]=="cumhaz") rr <- range(c(0,x$cumhaz))
   ylimo <- ylim
   if (is.null(ylim)) ylim <- rr
   if (is.null(xlim)) xlim <- range(x$times)
   if (is.null(x$se.cumhaz) & se==TRUE)  {
	  warning("predict.phreg must be with se=TRUE\n"); 
          se <- FALSE
   }
   if (se==TRUE) {
      if (is.null(x$se.cumhaz)) stop("predict.phreg must be with se=TRUE\n"); 
   if (type[1]=="surv") rrse <- range(c(x$surv.upper,x$surv.lower)) 
   if (type[1]=="cumhaz") {
	   cumhaz.upper <- x$cumhaz+level*x$se.cumhaz
	   cumhaz.lower <- x$cumhaz-level*x$se.cumhaz
       rrse <- range(c(cumhaz.upper,cumhaz.lower)) 
   }
   if (type[1]=="surv") rrse <- c(0,1)
   if (type[1]=="cumhaz") rrse <- c(max(0,rrse[1]),rrse[2])
   if (is.null(ylimo)) ylim <- rrse
   }

   ## all covriates 
   nx <- nrow(x$surv)
   if (is.null(whichx)) whichx <- 1:nx
   stratas <- whichx

   ltys <- lty
   cols <- col

   if (length(whichx)>0 ) { ## with X 

      if (!is.matrix(lty)) {
         if (is.null(lty)) ltys <- 1:length(whichx) else 
		 if (length(lty)!=length(whichx)) ltys <- rep(lty[1],length(whichx)) else ltys <- lty
      } else ltys <- lty

      if (!is.matrix(col)) {
         if (is.null(col)) cols <- 1:length(stratas) else 
		 if (length(col)!=length(stratas)) cols <- rep(col[1],length(stratas))
      } else cols <- col
   } else { 
     if (is.matrix(col))  cols <- col
     if (is.null(col)) cols <- 1  else cols <- col[1]
     if (is.matrix(lty))  ltys <- lty
     if (is.null(lty)) ltys <- 1  else ltys <- lty[1]
   }

  if (!is.matrix(ltys))  ltys <- cbind(ltys,ltys,ltys)
  if (!is.matrix(cols))  cols <- cbind(cols,cols,cols)

  i <- 1
  j <- whichx[i]

  cifreg <- FALSE
  if ((length(class(x))==2) && (substr(class(x)[2],1,3)=="cif")) cifreg <- TRUE
  if (type[1]=="surv") { 
	  xx <- x$surv 
	  if (cifreg)  xx <- x$cif 
  } else if (type[1]=="cif") {
	  xx <- 1-x$surv 
  } else xx <- x$cumhaz
  if (se) {
     if (type[1]=="surv") {
	     upper <- x$surv.upper; lower <- x$surv.lower 
	     if (cifreg) { upper <- x$cif.upper; lower <- x$cif.lower} 
     } else if (type[1]=="cif") {
	     upper <- 1-x$surv.lower; lower <- 1-x$surv.upper

     } else { upper <- cumhaz.upper; lower <- cumhaz.lower}
  }

  if (!add) 
  plot(x$times,xx[j,],type="s",lty=ltys[i,1],col=cols[i,1],ylim=ylim,ylab=ylab,xlim=xlim,xlab=xlab,...)
  else lines(x$times,xx[j,],type="s", lty=ltys[i,1],col=cols[i,1],...)
  if (length(whichx)>1) 
  for (i in seq(2,length(whichx))) lines(x$times,xx[whichx[i],],type="s",lty=ltys[i,1],col=cols[i,1],...)

    if (se==TRUE) {
    for (i in seq(1,length(whichx))) {
      j <- whichx[i]
      ul <- upper[j,]; nl <- lower[j,]; 
      if (!polygon) {
      lines(x$times,nl,type="s",lty=ltys[i,2],col=cols[i,2])
      lines(x$times,ul,type="s",lty=ltys[i,3],col=cols[i,3])
      } else {
	 ll <- length(x$times)
         tt <- c(x$times,rev(x$times))
         yy <- c(nl,rev(ul))
         ttp <- c(x$times[1],rep(x$times[-c(1,ll)],each=2),x$times[ll])
         tt <- c(ttp,rev(ttp))
         yy <- c(rep(nl[-ll],each=rep(2)),rep(rev(ul[-ll]),each=2))
         col.alpha<-0.1
         col.ci<-cols[i,1]
         col.trans <- sapply(col.ci, FUN=function(x) 
         do.call(grDevices::rgb,as.list(c(grDevices::col2rgb(x)/255,col.alpha))))
	 polygon(tt,yy,lty=0,col=col.trans)
     }
    }
  }

  where <- "topleft"; 
  if (type[1]=="surv") where <-  "topright"

}# }}}

###}}} predict

###{{{ plot

##' Plotting the baslines of stratified Cox 
##'
##' Plotting the baslines of stratified Cox 
##' @param x phreg object
##' @param se to include standard errors
##' @param time to plot for specific time variables
##' @param add to add to previous plot 
##' @param ylim to give ylim 
##' @param xlim to give xlim 
##' @param lty to specify lty of components
##' @param col to specify col of components
##' @param lwd to specify lwd of components
##' @param legend to specify col of components
##' @param ylab to specify ylab 
##' @param xlab to specify xlab 
##' @param polygon to get standard error in shaded form
##' @param level of standard errors
##' @param stratas wich strata to plot 
##' @param robust to use robust standard errors if possible
##' @param conf.type "plain" or "log" transformed 
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases basehazplot.phreg  bplot  basecumhaz plotConfRegion  plotConfRegionSE plotstrata
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' 
##' par(mfrow=c(2,2))
##' bplot(out1)
##' bplot(out1,stratas=c(0,3))
##' bplot(out1,stratas=c(0,3),col=2:3,lty=1:2,se=TRUE)
##' bplot(out1,stratas=c(0),col=2,lty=2,se=TRUE,polygon=FALSE)
##' bplot(out1,stratas=c(0),col=matrix(c(2,1,3),1,3),lty=matrix(c(1,2,3),1,3),se=TRUE,polygon=FALSE)
##' @export
basehazplot.phreg  <- function(x,se=FALSE,time=NULL,add=FALSE,ylim=NULL,xlim=NULL,
    lty=NULL,col=NULL,lwd=NULL,legend=TRUE,ylab=NULL,xlab=NULL,
    polygon=TRUE,level=0.95,stratas=NULL,robust=FALSE,conf.type=c("plain","log"),...) {# {{{
	if (inherits(x,"phreg") & is.null(ylab)) ylab <- "Cumulative hazard"
	if (inherits(x,"km") & is.null(ylab)) ylab <- "Survival probability"
	if (inherits(x,"cif") & is.null(ylab)) ylab <- "Probability"
	if (is.null(xlab)) xlab <- "time"
   level <- -qnorm((1-level)/2)
   rr <- range(x$cumhaz[,-1]) 
   strat <- x$strata[x$jumps]
   ylimo <- ylim
   if (is.null(ylim)) ylim <- rr
   if (is.null(xlim)) xlim <- range(x$cumhaz[,1])
   if (se==TRUE) {
	   if (is.null(x$se.cumhaz) & is.null(x$robse.cumhaz) ) 
		   stop("phreg must be with cumhazard=TRUE\n"); 
       if (conf.type[1]=="plain")
       rrse <- range(c(x$cumhaz[,-1]+level*x$se.cumhaz[,-1])) 
       else {
	       relse <- exp(level*x$se.cumhaz[,-1]/x$cumhaz[,-1])
	       rrse <- range(x$cumhaz[,-1]/relse,x$cumhaz[,-1]*relse) 
       }
       if (inherits(x,"km")) rrse <- c(min(x$lower,na.rm=TRUE),1)
       if (is.null(ylimo)) ylim <- rrse
   }

   ## all strata
   if (is.null(stratas)) stratas <- 0:(x$nstrata-1) 
   ltys <- lty
   cols <- col
   lwds <- lwd

   if (length(stratas)>0 & x$nstrata>1) { ## with strata
   lstrata <- x$strata.level[(stratas+1)]
   stratn <-  substring(x$strata.name,8,nchar(x$strata.name)-1)
   stratnames <- paste(stratn,lstrata,sep=":")
   
      if (!is.matrix(lty)) {
         if (is.null(lty)) lty <- ltys <- 1:length(stratas) 
         if (length(lty)!=length(stratas)) ltys <- rep(lty[1],length(stratas)) else ltys <- lty
      } else ltys <- lty
      if (!is.matrix(col)) {
         if (is.null(col)) col <- cols <- 1:length(stratas)  
	 if (length(col)!=length(stratas)) cols <- rep(col[1],length(stratas)) else cols <- col
      } else cols <- col
      if (!is.matrix(lwd)) {
         if (is.null(lwd)) lwd <- lwds <- rep(1,length(stratas))  
	 if (length(lwd)!=length(stratas)) lwds <- rep(lwd[1],length(stratas)) else lwds <- lwd
      } else lwds <- lwd
   } else { 
     stratnames <- "Baseline" 
     if (is.matrix(col))  cols <- col
     if (is.null(col)) cols <- 1  else cols <- col[1]
     if (is.matrix(lty))  ltys <- lty
     if (is.null(lty)) ltys <- 1  else ltys <- lty[1]
     if (is.matrix(lwd))  lwds <- lwd
     if (is.null(lwd)) lwds <- 1  else lwds <- lwd[1]
   }

  if (!is.matrix(ltys))  ltys <- cbind(ltys,ltys,ltys)
  if (!is.matrix(cols))  cols <- cbind(cols,cols,cols)
  if (!is.matrix(lwds))  lwds <- cbind(lwds,lwds,lwds)


  first <- 0
  for (i in seq(stratas)) {
      j <- stratas[i]
        cumhazard <- x$cumhaz[strat==j,,drop=FALSE]
        if (!is.null(cumhazard)) {
	if (nrow(cumhazard)>1) {
        if (add | first==1) 
        lines(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1],lwd=lwds[i,1])   
       else {
	  first <- 1
          plot(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1],lwd=lwds[i,1],ylim=ylim,ylab=ylab,xlab=xlab,xlim=xlim,...)
       }
       if (se==TRUE) {
	    if (robust==TRUE) secumhazard  <- x$robse.cumhaz[strat==j,,drop=FALSE]
	    else secumhazard <- x$se.cumhaz[strat==j,,drop=FALSE]
	    xx <- cumhazard[,2]
	    std.err <- secumhazard[,2]
	 if (conf.type[1] == "log") {
	    temp1 <-  exp(log(xx) + level * std.err/xx)
	    temp2 <-  exp(log(xx) - level * std.err/xx)
	    ul = cbind(cumhazard[,1],temp1); 
	    nl <- cbind(cumhazard[,1],temp2);
	 } 
	 if (conf.type[1] == "plain") {
		 ul <-cbind(cumhazard[,1],cumhazard[,2]+level*secumhazard[,2])
		 nl <-cbind(cumhazard[,1],cumhazard[,2]-level*secumhazard[,2])
	 }
	 if (inherits(x,"km")) { ul[,2] <- x$upper[x$strata==j]; 
	                          nl[,2] <- x$lower[x$strata==j];
	                          wna <- which(is.na(ul[,2]))
	                          ul[wna,2] <- 0
	                          nl[wna,2] <- 0

	 }
	      if (!polygon) {
		  lines(nl,type="s",lty=ltys[i,2],col=cols[i,2],lwd=lwds[i,2],...)
		  lines(ul,type="s",lty=ltys[i,3],col=cols[i,3],lwd=lwds[i,3],...)
	      } else {
                 ## type="s" confidence regions
		 ll <- length(nl[,1])
		 timess <- nl[,1]
		 ttp <- c(timess[1],rep(timess[-c(1,ll)],each=2),timess[ll])
		 tt <- c(ttp,rev(ttp))
		 yy <- c(rep(nl[-ll,2],each=rep(2)),rep(rev(ul[-ll,2]),each=2))
		 col.alpha<-0.1
		 col.ci<-cols[j+1]
		 col.trans <- sapply(col.ci, FUN=function(x) 
			   do.call(grDevices::rgb,as.list(c(grDevices::col2rgb(x)/255,col.alpha))))
		 polygon(tt,yy,lty=0,col=col.trans)
	      }
        }
     }
     }
    }


    where <- "topleft"; 
    if (inherits(x,"km")) where <-  "topright"
    if (legend & (!add)) 
    graphics::legend(where,legend=stratnames,col=cols[,1],lty=ltys[,1])

}# }}}


##' @export
plotConfRegion <- function(x,band,add=TRUE,polygon=TRUE,col=1,type="s",...)
{# {{{
nl <- cbind(x,band[,1])
ul <- cbind(x,band[,2])

  if (!polygon) {
      lines(nl,type=type,...)
      lines(ul,type=type,...)
      } else {
	 ll <- length(nl[,1])
         timess <- nl[,1]
         ttp <- c(timess[1],rep(timess[-c(1,ll)],each=2),timess[ll])
         tt <- c(ttp,rev(ttp))
         yy <- c(rep(nl[-ll,2],each=rep(2)),rep(rev(ul[-ll,2]),each=2))
         col.alpha<-0.1
         col.ci<-col[1]
         col.trans <- sapply(col.ci, FUN=function(x) 
           do.call(grDevices::rgb,as.list(c(grDevices::col2rgb(x)/255,col.alpha))))
	 polygon(tt,yy,lty=0,col=col.trans,...)
      }
}# }}}

##' @export
plotConfRegionSE <- function(x,est,se,...)
{# {{{
ul <- est+1.96*se; nl <- est-1.96*se
plotConfRegion(x,cbind(nl,ul),...)
}# }}}


##' @export
bplot <- function(x,...) basehazplot.phreg(x,...)

##' @export
basecumhaz <- function(x,type="matrix",robust=FALSE,...) {# {{{
   ## all strata
   strat <- x$strata[x$jumps]
   stratas <- 0:(x$nstrata-1) 

   ###   se.cum <- cum <- x$cumhaz
   se.cum <- cum <- c()
   strata <- rep(0,nrow(x$cumhaz))
   if (type=="matrix") { se.cum <- cum <- x$cumhazard }
   if (robust==TRUE) secum <- x$robse.cumhaz else secum  <- x$se.cumhaz
   if (is.null(secum)) nose <- TRUE else nose <- FALSE

   start <- 1
   for (i in stratas) {
	   cumhazard <- x$cumhaz[strat==i,,drop=FALSE]
	   if (!is.null(cumhazard)) {
		   nr <- nrow(cumhazard)
		   if (nr>=1) {
		   slut <- start-1+nr
		   cum <- rbind(cum,cumhazard)
		   if (!nose) se.cum <- rbind(se.cum,secum[strat==i,])
		   strata[start:slut] <- i
		   start <- slut+1
	      }
	   }
   }

   list(cumhaz=cum,se.cumhaz=se.cum,strata=strata)
}# }}}

##' @export
lines.phreg <- function(x,...,add=TRUE) plot(x,...,add=add)

###}}} plot

###{{{ plot

##' @export
plot.phreg  <- function(x,...) {
bplot(x,...)
}

##' @export
plotstrata <- function(x,y,strata,add=FALSE,where="topright",legend=TRUE,...) 
{ # {{{
  ## all strata
  stratas <- sort(unique(strata))
  cols <- 1:length(stratas)
  ltys <- 1:length(stratas)
  xlim <- range(x)
  ylim <- range(y)

  cumhaz <- cbind(x,y)

  k <- first <- 0
  for (i in stratas) {
    k <- 1+k
    cumhazs <- cumhaz[strata==i,,drop=FALSE]
    if (!is.null(cumhazs)) {
       if (nrow(cumhazs)>1) {
           if (add | first==1) 
           lines(cumhazs,type="s",lty=ltys[k],col=cols[k])   
	 else { first <- 1
         plot(cumhazs,type="s",lty=ltys[k],col=cols[k],ylim=ylim,xlim=xlim,...)
         } 
    } 
  }
  }

  if (legend & (!add)) 
  graphics::legend(where,legend=paste(stratas),col=cols,lty=ltys)

}  # }}}


##' @export
lines.phreg <- function(x,...,add=TRUE) plot(x,...,add=add)

###}}} plot

###{{{ print
##' @export
print.phreg  <- function(x,...) {
  cat("Call:\n")
  dput(x$call)
  print(summary(x),...)
}
###}}} print

