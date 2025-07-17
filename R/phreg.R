
##' @export
namesortme <- function(iid,name.id,sort=TRUE) { ## {{{
if (is.matrix(iid))  
	if (nrow(iid)==length(name.id)) {
		rownames(iid) <- name.id
		if (sort)  {
		oid <- order(name.id)
		iid <- iid[oid,,drop=FALSE]
		}
}
return(iid)
} ## }}}

##' @export
construct_id <- function(id,nid,namesX=NULL,as.data=FALSE) { ## {{{ 
  call.id <- id

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
     order.ids <- order(ids)
     id.name <- ids[order.ids]
   } else { 
	   id <- 1:nid-1;  
	   ids <- id+1
	   order.ids <- ids
###	   if (!is.null(namesX)) id.name <- namesX else id.name <- ids
           id.name <- ids
   } 
   ## orginal id coding into integers 
   ## id from 0,1,...,nid-1

    if (as.data) {
        id  <- (0:(nid - 1))[order(ids)][id +1]
        id.name <- ids
    }

  return(list(call.id=call.id,id=id,nid=nid,unique.id=ids,name.id=id.name))
} ## }}} 

###{{{ phreg01
phreg01 <- function(X,entry,exit,status,id=NULL,strata=NULL, offset=NULL,weights=NULL,strata.name=NULL,cumhaz=TRUE,
             beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
	     case.weights=NULL,no.var=0,augmentation=0,...) {
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

  call.id <- id; 
  conid <- construct_id(id,nrow(X),namesX=rownames(X))
  id <- conid$id; nid <- conid$nid; name.id <- conid$name.id

  dd <- .Call("FastCoxPrepStrata", entry,exit,status,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")

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
	      val$gradient <- val$gradient+augmentation
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
         II <- - tryCatch(solve(val$hessian),error=function(e) matrix(0,nrow(val$hessian),ncol(val$hessian)) )
  } else II <- matrix(0,p,p)

  ## Brewslow estimator, to handle also possible weights, caseweights that are 0
  ww <- val$caseweightsJ * val$weightsJ
  val$S0[abs(ww)<.00000001] <- 0 

  ### computes Breslow estimator 
  if (cumhaz==TRUE & (length(val$jumps)>0)) { # {{{
	 strata <- val$strata[val$jumps]
	 nstrata <- val$nstrata
	 jumptimes <- val$jumptimes

	 ## Brewslow estimator, to handle also possible weights/case-weights, and caseweights that are 0
	 S0i2 <- S0i <- rep(0,length(val$S0))
         wwJ <- val$caseweightsJ*val$weightsJ
	 S0i[val$S0>0] <- 1/val$S0[val$S0>0]
	 S0i2[val$S0>0] <- 1/(val$S0[val$S0>0]^2*wwJ[val$S0>0])
	 cumhaz <- cbind(jumptimes,cumsumstrata(S0i,strata,nstrata))
	 if ((no.opt==FALSE & p!=0)) { 
	     DLambeta.t <- apply(val$E*S0i,2,cumsumstrata,strata,nstrata)
	     varbetat <-   rowSums((DLambeta.t %*% II)*DLambeta.t)
	 ### covv <-  apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
	 } else varbetat <- 0
	 var.cumhaz <- cumsumstrata(S0i2,strata,nstrata)+varbetat
	 se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)

	 colnames(cumhaz)    <- c("time","cumhaz")
	 colnames(se.cumhaz) <- c("time","se.cumhaz")
 } # }}} 
 else {cumhaz <- se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL}

  res <- c(val,
           list(cox.prep=dd,
                strata.call=strata.call, strata.level=strata.level,
                entry=entry, exit=exit, status=status, X=X,
                p=p, offsets=offset, weights=weights,
                id=id,call.id=call.id,nid=nid,name.id=name.id,
                opt=opt, no.opt=no.opt, cumhaz=cumhaz, se.cumhaz=se.cumhaz,
                lcumhaz=lcumhaz, lse.cumhaz=lse.cumhaz, ihessian=II,
                II=II,strata.name=strata.name,propodds=propodds))
  class(res) <- "phreg"

  ## also computing robust variance 
  if (p>0 & no.var==0) {
  beta.iid <- iid(res)
  phvar <- crossprod(beta.iid)
###  estimate.beta <- estimate(coef=val$coef,IC=nrow(beta.iid)*beta.iid,id=name.id)
###  res$estimate.beta <- estimate.beta
  colnames(phvar) <- rownames(phvar) <- names(res$coef)
  res$var <- phvar
  res$beta.iid <- beta.iid
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
##' @param offset offsets for Cox model
##' @param weights weights for Cox score equations
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases phreg phreg.par robust.phreg readPhreg conftype plotO.predictphreg plotpredictphreg predictO.phreg predictrecreg summarybase.phreg namesortme construct_id
##' @examples
##' library(mets)
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' summary(out1)
##' 
##' par(mfrow=c(1,2))
##' plot(out1)
##' 
##' ## computing robust variance for baseline
##' rob1 <- robust.phreg(out1)
##' plot(rob1,se=TRUE,robust=TRUE)
##' 
##' ## iid decomposition, with scaled influence functions
##' ## for regression parameters
##' head(iid(out1))
##' ## making iid decomposition of baseline at a specific time-point
##' Aiiid <- iid(out1,time=30)
##' head(Aiiid)
##' ## both iid decompositions
##' dd <- iidBaseline(out1,time=30)
##' head(dd$beta.iid)
##' head(dd$base.iid)
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

###{{{ iid & Robust variances

##' @export

##' Influence functions or IID decomposition of baseine for recrec/phreg/cifregFG
##'
##' @title Influence functions or IID decomposition of baseine for recrec/phreg/cifregFG
##' @param object phreg/recreg/cifregFG object
##' @param time for baseline IID 
##' @param ft function to compute IID of baseline integrated against f(t) 
##' @param fixbeta to fix the coefficients 
##' @param beta.iid to use these iid of beta 
##' @param tminus to get predictions in t-  
##' @param sort to sort after object$name.id when returning iid decomposition
##' @param ... additional arguments to lower level functions
##' @author Thomas Scheike
##' @aliases iidBaseline
##' @export
iidBaseline <- function(object,time=NULL,ft=NULL,fixbeta=NULL,beta.iid=NULL,tminus=FALSE,sort=TRUE,...) UseMethod("iidBaseline")
##' @export 
iidBaseline.phreg <- function(object,time=NULL,ft=NULL,fixbeta=NULL,beta.iid=NULL,tminus=FALSE,sort=FALSE,...)
{# {{{
###  sum_i int_0^t f(s)/S_0(s) dM_{ki}(s) - P(t) \beta_k
###  with possible strata and cluster "k", and i in clusters 
  x <- object
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
  ww <- xx$caseweights*xx$weights
  S0i2[xx$jumps+1] <- 1/(x$S0^2*ww[xx$jump+1])
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
  U[xx$jumps+1,] <- x$U
  ### only up to time t for A
  if (is.null(ft))  ft <- rep(1,length(xx$time))
  if (is.matrix(ft)) ft <- c(ft)
  cumhaz <- c(cumsumstrata(ft*S0i,xx$strata,xx$nstrata))
  cumS0i2 <- c(cumsumstrata(ft*S0i2*btimexx,xx$strata,xx$nstrata))
  if (fixbeta==0) {
	  EdLam0 <- apply(ft*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  Ht <- apply(ft*E*S0i*btimexx,2,cumsumstrata,xx$strata,xx$nstrata)
  } else Ht <- NULL
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset)) else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  MGAiid <- ft*S0i*btimexx-cumS0i2*rr*c(xx$weights)

  MGtiid <- NULL
  if (fixbeta==0) {# {{{
     if (!is.null(beta.iid)) MGtiid <- beta.iid else {
     invhess <- -solve(x$hessian)
     MGt <- ft*U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx$weights)
     MGt <- MGt %*% invhess
     MGtiid <- apply(MGt,2,sumstrata,id,mid)
    }
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
 if (is.matrix(MGAiid)) {
    colnames(MGAiid) <- paste("strata",sus,sep="")
    MGAiid <- namesortme(MGAiid,x$name.id,sort=sort)
    ###    if (length(x$name.id)==nrow(MGAiid)) rownames(MGAiid) <- x$name.id
 }
 MGtiid <- namesortme(MGtiid,x$name.id,sort=sort)

 return(list(time=time,base.iid=MGAiid,strata=xx$strata,nstrata=xx$nstrata,
	     id=id,beta.iid=MGtiid,model.frame=x$model.frame,formula=x$formula))
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

###{{{ phreg: coef vcov print summary print.summary  plot lines 

##' @export
coef.phreg  <- function(object,...) { ## {{{
  object$coef
} ## }}}

##' @export
vcov.phreg  <- function(object,...) {     ## {{{
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
} ## }}}


##' @export
summary.phreg <- function(object,type=c("robust","martingale"),augment.type=c("var.augment.times","var.augment"),...) { ## {{{
  expC <- cc <- ncluster <- V <- NULL

   if (length(object$p)>0 & object$p>0 & (!object$no.opt)) {
    I <- -solve(object$hessian)
    if ( (length(class(object))==2) && ( inherits(object,c("cifreg","recreg")))) {
	    V <- object$var
	    ncluster <- object$ncluster ## nrow(object$Uiid)
            if (!is.null(object$augmentation)) { V <- object[[augment.type[1]]]; }
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
} ## }}}

##' @export
summarybase.phreg <- function(object,robust=FALSE,...) { ## {{{
  out <- summaryRecurrentobject(object,robust=robust,...)
  class(out) <- "summary.recurrent"
  return(out)
}# }}}


##' @export
print.phreg  <- function(x,...) { ## {{{
  cat("Call:\n")
  dput(x$call)
  print(summary(x),...)
}  ## }}}

##' Plotting the baselines of stratified Cox 
##'
##' Plotting the baselines of stratified Cox
##' @param x phreg object
##' @param ... Additional arguments to baseplot funtion
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases basehazplotO.phreg  baseplot bplot  basecumhaz plotConfRegion  plotConfRegionSE plotstrata kmplot plotConfregion
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' 
##' par(mfrow=c(2,2))
##' plot(out1)
##' plot(out1,stratas=c(0,3))
##' plot(out1,stratas=c(0,3),col=2:3,lty=1:2,se=TRUE)
##' plot(out1,stratas=c(0),col=2,lty=2,se=TRUE,polygon=FALSE)
##' plot(out1,stratas=c(0),col=matrix(c(2,1,3),1,3),lty=matrix(c(1,2,3),1,3),se=TRUE,polygon=FALSE)
##' @aliases plotConfRegion  plotConfRegionSE plotstrata kmplot plotConfregion
##' @export
plot.phreg  <- function(x,...)  baseplot(x,...) 

##' @export
lines.phreg <- function(x,...,add=TRUE) plot(x,...,add=add)

##' @export
bplot <- function(x,...,add=TRUE) plot(x,...,add=add)

##' @export
print.summary.phreg  <- function(x,max.strata=5,...) { ## {{{

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

} ## }}}

###}}} 

##' @export
conftype <- function(x,std.err,conf.type=c("log","plain"),restrict=c("positive","prob","none"),conf.int=0.95)
{ ## {{{
zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)
 if (conf.type[1] == "plain") {# {{{
    temp1 <- x + zval * std.err 
    temp2 <- x - zval * std.err 
    if (restrict[1]=="prob") cix <- list(upper = pmin(temp1, 1), lower = pmax(temp2, 0), conf.type = "plain", conf.int = conf.int)
    if (restrict[1]=="positive") cix <- list(upper = pmax(temp1, 0), lower = pmax(temp2, 0), conf.type = "plain", conf.int = conf.int)
    if (restrict[1]=="none") cix <- list(upper = temp1, lower = temp2, conf.type = "plain", conf.int = conf.int)
 }
 if (conf.type[1] == "log") {
    xx <- ifelse(x <= 0, NA, x)
    temp1 <- ifelse(x <= 0, NA, exp(log(xx) + zval * std.err/xx))
    temp2 <- ifelse(x <= 0, NA, exp(log(xx) - zval * std.err/xx))
    if (restrict[1]=="prob") cix <- list(upper = pmin(temp1, 1), lower = temp2, conf.type = "log", conf.int = conf.int)
    else cix <- list(upper = temp1, lower = temp2, conf.type = "log", conf.int = conf.int)
 }
### if (conf.type[1] == "log-log") {
###    who <- (x == 0 | x == 1)
###    xx <- ifelse(who, NA, x)
###    temp1 <- exp(-exp(log(-log(xx)) + zval * std.err/log(xx)))
###    temp2 <- exp(-exp(log(-log(xx)) - zval * std.err/log(xx)))
###    if (restrict[1]=="prob") {
###	    temp3 <- ifelse(x == 0, NA, 1)
###	    temp1 <- ifelse(who, temp3, temp1)
###	    temp2 <- ifelse(who, temp3, temp2)
###    } else { 
###    cix <- list(upper = temp1, lower = temp2, conf.type = "log-log", conf.int = conf.int)
###    }
### }# }}}
 return(cix)
}# }}}

##' @export
basecumhaz <- function(x,type=c("list"),only=0,joint=0,cumhaz="cumhaz",
		       se.cumhaz="se.cumhaz",robust=FALSE,...) {# {{{
   ## all strata
   stratobs <- x$strata[x$jumps]
   ustratobs <- sort(unique(stratobs))
   stratas <- 0:(x$nstrata-1) 

   se.cum <- cum <- c()
   se.cum <- cum <- x[[cumhaz]]
   if (robust==TRUE) {
	   ## take robse if there otherwise stay with se.cumhaz
	   secum <- x$robse.cumhaz 
	   if (is.null(secum)) secum <- x[[se.cumhaz]]
   } else secum  <- x[[se.cumhaz]]
   if (is.null(secum)) nose <- TRUE else nose <- FALSE

   out <- rep(list(NULL),x$nstrata)
   if (length(ustratobs)>0)
   for (i in ustratobs) {
	   if (!is.null(x[[cumhaz]])) {
	   cumhazard <- x[[cumhaz]][stratobs==i,,drop=FALSE]
           nr <- nrow(cumhazard)
	   if (nr>=1) {
		   if (!nose) se.cum <- secum[stratobs==i,,drop=FALSE] else se.cum <- NULL
		   if (only==0) {
			  if (joint==0) out[[i+1]] <- list(cumhaz=cumhazard,se.cumhaz=se.cum,strata=i) else 
		                        out[[i+1]] <- list(cumhaz=cbind(cumhazard,se.cum[,2],i)) 
	           } else out[[i+1]] <- cumhazard
         } 
	 } 
   }

   attr(out,"stratobs") <- ustratobs
   return(out) 
}# }}}

###
###basecumhaz <- function(x,type=c("list"),only=0,joint=0,robust=FALSE,...) {# {{{
###   ## all strata
###   stratobs <- x$strata[x$jumps]
###   ustratobs <- sort(unique(stratobs))
###   stratas <- 0:(x$nstrata-1) 
###
###   se.cum <- cum <- c()
###   se.cum <- cum <- x$cumhaz 
###   if (robust==TRUE) {
###	   ## take robse if there otherwise stay with se.cumhaz
###	   secum <- x$robse.cumhaz 
###	   if (is.null(secum)) secum <- x$se.cumhaz
###   } else secum  <- x$se.cumhaz
###   if (is.null(secum)) nose <- TRUE else nose <- FALSE
###
###   out <- rep(list(NULL),x$nstrata)
###   if (length(ustratobs)>0)
###   for (i in ustratobs) {
###	   if (!is.null(x$cumhaz)) {
###	   cumhazard <- x$cumhaz[stratobs==i,,drop=FALSE]
###           nr <- nrow(cumhazard)
###	   if (nr>=1) {
###		   if (!nose) se.cum <- secum[stratobs==i,,drop=FALSE] else se.cum <- NULL
###		   if (only==0) {
###			  if (joint==0) out[[i+1]] <- list(cumhaz=cumhazard,se.cumhaz=se.cum,strata=i) else 
###		                        out[[i+1]] <- list(cumhaz=cbind(cumhazard,se.cum[,2],i)) 
###	           } else out[[i+1]] <- cumhazard
###         } 
###	 } 
###   }
###
###   attr(out,"stratobs") <- ustratobs
###   return(out) 
###}# }}}
###

##' @export
baseplot  <- function(x,se=FALSE,time=NULL,add=FALSE,ylim=NULL,xlim=NULL,
	 lty=NULL,col=NULL,lwd=NULL,legend=TRUE,ylab="Cumulative hazard",xlab="time",
	 polygon=TRUE,level=0.95,stratas=NULL,robust=FALSE,cumhaz="cumhaz",se.cumhaz="se.cumhaz",
 conf.type=c("log","plain"),restrict = c("positive","prob", "none"),...) {# {{{

   base <- basecumhaz(x,joint=1,robust=robust,cumhaz=cumhaz,se.cumhaz=se.cumhaz)
   nstrata <- x$nstrata
   stratobs <- attr(base,"stratobs")
   ###
   if (is.null(stratas)) stratas <- stratobs else { 
         wm <-   match(stratas,stratobs)
         wm <- wm[!is.na(wm)]
         if (length(wm)==0) stratas <- rep(0,0) else stratas <- stratobs[wm]
   }

   ## all confidence intervals for relevant strata  to make ylim
   if (length(stratobs)>0) {
   xx <- conftype(x$cumhaz[,2],x$se.cumhaz[,2],conf.type=conf.type[1],restrict=restrict[1],conf.int=level)

   if (is.null(xlim)) {
	   xlim <-  range(x$cumhaz[,1],na.rm=TRUE) 
   }
   if (se) {
       if (is.null(x$se.cumhaz) & is.null(x$robse.cumhaz)) stop("phreg must be with cumhazard=TRUE\n"); 
       rrse <- range(c(xx$lower,xx$upper),na.rm=TRUE) 
       if (is.null(ylim))
       { 
        if (Inf %in% abs(rrse)) ylim <- range(x$cumhaz[,2],na.rm=TRUE) else ylim <- rrse
       } 
   } else {
       if (is.null(ylim)) ylim <- range(x$cumhaz[,2],na.rm=TRUE)
   }

   ## all strata
   ## {{{
   ltys <- lty; cols <- col; lwds <- lwd

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
  ## }}}

  first <- 1
  j <- 1
  for (i in stratas)  {
	if (!is.null(base[[i+1]])) {
	cumhazard <- base[[i+1]]$cumhaz
	if (nrow(cumhazard)>1) {
        if (add | first==0) lines(cumhazard,type="s",lty=ltys[j,1],col=cols[j,1],lwd=lwds[j,1])   
        else {
	  first <- 0
          plot(cumhazard,type="s",lty=ltys[j,1],col=cols[j,1],lwd=lwds[j,1],ylim=ylim,ylab=ylab,xlab=xlab,xlim=xlim,...)
        }
        if (se==TRUE) {
            xx <- conftype(cumhazard[,2],cumhazard[,3],conf.type=conf.type[1],restrict=restrict[1],conf.int=level)
            nl <- cbind(cumhazard[,1],xx$lower); ul <- cbind(cumhazard[,1],xx$upper)
	      if (!polygon) {
		  lines(nl,type="s",lty=ltys[j,2],col=cols[j,2],lwd=lwds[i,2],...)
		  lines(ul,type="s",lty=ltys[j,3],col=cols[j,3],lwd=lwds[j,3],...)
	      } else plotConfRegion(nl[,1],cbind(nl[,2],ul[,2]),col=cols[j,1])
        }
     }
     }
    j <- j+1
  }

    where <- "topleft"; 
    if (inherits(x,"km")) where <-  "topright"
    if (legend & (!add)) 
    graphics::legend(where,legend=stratnames,col=cols[,1],lty=ltys[,1])
   } else {
	   if (!add) plot(0,0,ylim=ylim,ylab=ylab,xlab=xlab,xlim=xlim,...)
   }
}# }}}

###{{{ phreg: basehazplot.phreg predictO.phreg plotO.predictphreg (old)

##' Plotting the baselines of stratified Cox 
##'
##' Plotting the baselines of stratified Cox
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
##' @aliases basehazplot.phreg  
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' 
##' par(mfrow=c(2,2))
##' plot(out1)
##' plot(out1,stratas=c(0,3))
##' plot(out1,stratas=c(0,3),col=2:3,lty=1:2,se=TRUE)
##' plot(out1,stratas=c(0),col=2,lty=2,se=TRUE,polygon=FALSE)
##' plot(out1,stratas=c(0),col=matrix(c(2,1,3),1,3),lty=matrix(c(1,2,3),1,3),se=TRUE,polygon=FALSE)
##' @export
basehazplot.phreg  <- function(x,se=FALSE,time=NULL,add=FALSE,ylim=NULL,xlim=NULL,lty=NULL,col=NULL,
			       lwd=NULL,legend=TRUE,ylab=NULL,xlab=NULL,polygon=TRUE,level=0.95,
			       stratas=NULL,robust=FALSE,conf.type=c("plain","log"),...) {# {{{
   if (inherits(x,"phreg") & is.null(ylab)) ylab <- "Cumulative hazard"
   if (is.null(xlab)) xlab <- "time"
   level <- -qnorm((1-level)/2)
   rr <- range(x$cumhaz[,-1]) 
   strat <- x$strata[x$jumps]
   ylimo <- ylim
   if (is.null(ylim)) ylim <- rr
   if (is.null(xlim)) xlim <- range(x$cumhaz[,1])
   if (se==TRUE) {
	   if (is.null(x$se.cumhaz) & is.null(x$robse.cumhaz)) stop("phreg must be with cumhazard=TRUE\n"); 
       if (conf.type[1]=="plain")
       rrse <- range(c(x$cumhaz[,-1]+level*x$se.cumhaz[,-1])) 
       else { ## type="log"
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
	      } else plotConfRegion(nl[,1],cbind(nl[,2],ul[,2]),col=cols[j+1])
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
predictO.phreg <- function(object,newdata,times=NULL,individual.time=FALSE,tminus=FALSE,se=TRUE,robust=FALSE,conf.type="log",conf.int=0.95,km=FALSE,...) 
{# {{{ default is all time-points from the object

   ### take baseline and strata from object# {{{
   ocumhaz <- object$cumhaz
   jumptimes <- ocumhaz[,1]
   chaz <- ocumhaz[,2]
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
 if (length(class(object))==2 && ( substr(class(object)[2],1,3)=="cif" | substr(class(object)[1],1,3)=="cif")) {
	 out <- c(out,list(cif=1-out$surv,cif.lower=1-out$surv.upper, cif.upper=1-out$surv.lower))
 }
if (length(class(object))==2 && ( substr(class(object)[2],1,3)=="rec" | substr(class(object)[1],1,3)=="rec")) {
	 out <- c(out,list(mean=-log(out$surv),mean.lower=-log(out$surv.upper),mean.upper=-log(out$surv.lower)))
 }

   class(out) <- c("predictphreg",class(object)[1])
### class(out) <- c("predictphreg")
### if (length(class(object))==2) class(out) <- c("predictphreg",class(object)[1])
 return(out)
}# }}}

##' @export
plotO.predictphreg  <- function(x,se=FALSE,add=FALSE,ylim=NULL,xlim=NULL,lty=NULL,col=NULL,type=c("default","surv","cumhaz","cif"),ylab=NULL,xlab=NULL,
    polygon=TRUE,level=0.95,whichx=NULL,robust=FALSE,...) {# {{{
   if (type[1]=="cumhaz" & is.null(ylab)) ylab <- "Cumulative hazard"
   if (type[1]=="cif" & is.null(ylab)) ylab <- "Cumulative probability"
   if (type[1]=="surv" & is.null(ylab)) ylab <- "Surival probability"
   if (type[1]=="default" & (class(x)[2]=="phreg"))  { if (is.null(ylab)) ylab <- "Survival probability";   type <- "surv"}
   if (type[1]=="default" & (class(x)[2]=="cifreg"))  {if (is.null(ylab))  ylab <- "Cumulative probability"; type <- "cif"}
   if (type[1]=="default" & (class(x)[2]=="recreg"))  {if (is.null(ylab))  ylab <- "Cumulative mean";        type <- "cumhaz"}

   if (is.null(xlab)) xlab <- "time"
   level <- -qnorm((1-level)/2)
   if (type[1]=="cif") rr <- c(0,1) 
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
   if (type[1]=="surv") rrse <- range(c(x$surv.upper,x$surv.lower)) 
   if (type[1]=="cumhaz") {
	   cumhaz.upper <- x$cumhaz+level*x$se.cumhaz
	   cumhaz.lower <- x$cumhaz-level*x$se.cumhaz
       rrse <- range(c(cumhaz.upper,cumhaz.lower)) 
   }
   if (type[1]=="surv") rrse <- c(0,1)
   if (type[1]=="cif") rrse <- c(0,1)
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

  recreg <- cifreg <- FALSE
  if ((length(class(x))==2) && (substr(class(x)[2],1,3)=="cif")) cifreg <- TRUE
###  if ((length(class(x))==2) && (substr(class(x)[2],1,3)=="rec")) recreg <- TRUE
  if (type[1]=="surv") { 
	  xx <- x$surv 
	  if (cifreg)  xx <- x$cif 
	  if (recreg)  xx <- x$mean 
  } else if (type[1]=="cif") {
	  xx <- 1-x$surv 
  } else xx <- x$cumhaz
  if (se) {
     if (type[1]=="surv") {
	     upper <- x$surv.upper; lower <- x$surv.lower 
	     if (cifreg) { upper <- x$cif.upper; lower <- x$cif.lower} 
	     if (recreg) { upper <- x$mean.upper; lower <- x$mean.lower} 
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
      } else plotConfRegion(x$times,cbind(nl,ul),col=cols[i,1])
    }
  }

  where <- "topleft"; 
  if (type[1]=="surv") where <-  "topright"

}# }}}

###}}} plot

##' @export
plotpredictphreg <- function(x,se=FALSE,add=FALSE,ylim=NULL,xlim=NULL,lty=NULL,col=NULL,
       type=c("surv","cumhaz","cif"),ylab=NULL,xlab=NULL,
       polygon=TRUE,whichx=NULL,robust=FALSE,...) {# {{{
   if (type[1]=="surv" & is.null(ylab)) ylab <- "Surival probability"
   if (type[1]=="cumhaz" & is.null(ylab)) ylab <- "Cumulative hazard"
   if (type[1]=="cif" & is.null(ylab)) ylab <- "Cumulative probability"

   if (is.null(xlim)) xlim <- range(x$times)
   if (is.null(xlab)) xlab <- "time"
   if (is.null(x$se.cumhaz) & se==TRUE)  {
	  warning("predict.phreg must be with se=TRUE\n"); 
          se <- FALSE
   }
   ylimo <- ylim
   nlower <- paste(type[1],".lower",sep="")
   nupper <- paste(type[1],".upper",sep="")
   xx <- x[[type[1]]]

   if (se==TRUE) {
      upper <- x[[nupper]]
      lower <- x[[nlower]]
      if (is.null(x$se.cumhaz)) stop("predict.phreg/cifreg/recreg must be with se=TRUE\n"); 
      rrse <- range(c(0,x[[nlower]],x[[nupper]]),na.rm=TRUE) 
      if (is.null(ylimo)) ylim <- rrse
   } else {
      rrse <- range(c(xx),na.rm=TRUE) 
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
      } else plotConfRegion(x$times,cbind(nl,ul),col=cols[i,1])
    }
  }

  where <- "topleft"; 
  if (type[1]=="surv") where <-  "topright"

}# }}}

###{{{ predict.phreg print.predict plot.predict.phreg

##' Predictions from proportional hazards model
##'
##' @param object phreg object
##' @param newdata data.frame
##' @param times Time where to predict variable, default is all time-points from the object sorted
##' @param individual.time to use one (individual) time per subject, and then newdata and times have same length and makes only predictions for these individual times.
##' @param tminus to make predictions in T- that is strictly before given times, useful for IPCW techniques
##' @param se with standard errors and upper and lower confidence intervals.
##' @param robust to get robust se's also default for most functions (uses robse.cumhaz otherwise se.cumhaz). 
##' @param conf.type transformation for suvival estimates, default is log
##' @param conf.int significance level
##' @param km to use Kaplan-Meier product-limit for baseline \deqn{S_{s0}(t)= (1 - dA_{s0}(t))}, otherwise take exp of cumulative baseline.
##' @param ... Additional arguments to plot functions
##' @aliases headstrata tailstrata revcumsumstrata revcumsumstratasum cumsumstrata sumstrata covfr covfridstrata covfridstrataCov cumsumidstratasum cumsumidstratasumCov cumsumstratasum revcumsum revcumsumidstratasum revcumsumidstratasumCov robust.basehaz.phreg matdoubleindex mdi cumsum2strata revcumsum2strata revcumsum2stratafdN
##' @export
predict.phreg <- function(object,newdata,times=NULL,individual.time=FALSE,tminus=FALSE,se=TRUE,robust=FALSE,conf.type="log",conf.int=0.95,km=FALSE,...) 
{# {{{ default is all time-points from the object
### take baseline and strata from object# {{{
ocumhaz <- object$cumhaz
jumptimes <- ocumhaz[,1]
chaz <- ocumhaz[,2]
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
} else { ## only relevant for phreg objects 
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
        where <- predictCumhaz(c(0,jumptimes[strata==j]),times,type="left",tminus=tminus)
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

  if (se)  {
     cibase <- conftype(pcumhaz,se.cumhaz,conf.type=conf.type[1],restrict="positive",conf.int=conf.int)
     cisurv <- conftype(surv,surv*se.cumhaz,conf.type=conf.type[1],restrict="prob",conf.int=conf.int)
  } else {
     cibase <- list(upper=NULL,lower=NULL,conf.type=conf.type[1],restrict="positive",conf.int=conf.int)
     cisurv <- list(upper=NULL,lower=NULL,conf.type=conf.type[1],restrict="positive",conf.int=conf.int)
  }

 if (object$p>0) RR <-  exp(X %*% coef(object)) else RR <- rep(1,nrow(X))

### non-cox setting
if (!is.null(object$propodds)) pcumhaz <- -log(surv)

 out <- list(surv=surv,times=times,
    surv.upper=cisurv$upper,surv.lower=cisurv$lower,cumhaz=pcumhaz,se.cumhaz=se.cumhaz,
    se.surv=surv*se.cumhaz,se.cif=surv*se.cumhaz,
    cif=1-surv,cif.lower=1-cisurv$upper,cif.upper=1-cisurv$lower,
    cumhaz.upper=cibase$upper,cumhaz.lower=cibase$lower,strata=strataNew,X=X,RR=RR)

class(out) <- c("predictphreg",class(object)[1])
return(out)
}# }}}

##' @export
summary.predictphreg <- function(object,times=NULL,type=c("cif","cumhaz","surv")[3],np=10,...) {# {{{
ret <- summary.predictrecreg(object,times=times,type=type[1],np=np,...)
return(ret)
}# }}}


##' @export
plot.predictphreg  <- function(x,se=FALSE,add=FALSE,ylim=NULL,xlim=NULL,lty=NULL,col=NULL,type=c("surv","cumhaz","cif"),ylab=NULL,xlab=NULL,
    polygon=TRUE,level=0.95,whichx=NULL,robust=FALSE,...) {# {{{
   if (type[1]=="cumhaz" & is.null(ylab)) ylab <- "Cumulative hazard"
   if (type[1]=="cif" & is.null(ylab)) ylab <- "Cumulative probability"
   if (type[1]=="surv" & is.null(ylab)) ylab <- "Surival probability"

   if (is.null(xlab)) xlab <- "time"
   if (type[1]=="cif") rr <- c(0,1) 
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
   if (type[1]=="surv") rrse <- range(c(x$surv.upper,x$surv.lower,na.rm=TRUE)) 
   if (type[1]=="surv") rrse <- range(c(x$surv.upper,x$surv.lower,na.rm=TRUE)) 
   if (type[1]=="cumhaz") {
	   cumhaz.upper <- x$cumhaz+level*x$se.cumhaz
	   cumhaz.lower <- x$cumhaz-level*x$se.cumhaz
       rrse <- range(c(cumhaz.upper,cumhaz.lower),na.rm=TRUE) 
   }
   if (type[1]=="surv") rrse <- c(0,1)
   if (type[1]=="cif") rrse <- c(0,1)
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

  recreg <- cifreg <- FALSE
  if ((length(class(x))==2) && (substr(class(x)[2],1,3)=="cif")) cifreg <- TRUE
###  if ((length(class(x))==2) && (substr(class(x)[2],1,3)=="rec")) recreg <- TRUE
  if (type[1]=="surv") { 
	  xx <- x$surv 
	  if (cifreg)  xx <- x$cif 
	  if (recreg)  xx <- x$mean 
  } else if (type[1]=="cif") {
	  xx <- 1-x$surv 
  } else xx <- x$cumhaz
  if (se) {
     if (type[1]=="surv") {
	     upper <- x$surv.upper; lower <- x$surv.lower 
	     if (cifreg) { upper <- x$cif.upper; lower <- x$cif.lower} 
	     if (recreg) { upper <- x$mean.upper; lower <- x$mean.lower} 
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
      } else plotConfRegion(x$times,cbind(nl,ul),col=cols[i,1])
    }
  }

  where <- "topleft"; 
  if (type[1]=="surv") where <-  "topright"

}# }}}

###}}} predict

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
##' computed for all event times as functions and can be plotted/viewed.  When times are given and beyond
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
##' par(mfrow=c(1,2)); plot(drm1,cause=1,se=1); plot(drm1,cause=2,se=1);
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
    intkmtimes <- meanm[newtimes,,drop=FALSE]
    se.intkmtimes <- se.resmean[newtimes]
    skmtimes <- mm[newtimes,2]
    years.lost <- intkmtimes[,1]-intkmtimes[,2]
    intkmtimes <- cbind(skmtimes,intkmtimes,se.intkmtimes,years.lost)
    colnames(intkmtimes) <- c("strata","times","rmean","se.rmean","years.lost")
    rownames(intkmtimes) <- rep(x$strata.level,length(times)); 
    intkmtimes=data.frame(intkmtimes)
###    logintkmtimes=cbind(intkmtimes[,1:2],log(intkmtimes[,3]),se.intkmtimes/intkmtimes[,3])
###    colnames(logintkmtimes) <- c("strata","times","log-rmean","log-se.rmean")
  } else intkmtimes <- se.intkmtimes <- NULL

  rmst <- cbind(meanm[,1],vecAllStrata(meanm[,2],strata.jumps,nstrata))
  se.rmst <- cbind(se.mm[,1],vecAllStrata(se.mm[,2],strata.jumps,nstrata))

 out <- list(cumhaz=meanm,se.cumhaz=se.mm,covs=covs,
       time=time,strata=strata.jumps,nstrata=nstrata,
       jumps=1:length(km),strata.name=x$strata.name,
       strata.level=x$strata.level,
       intkmtimes=intkmtimes,
       rmst=rmst,se.rmst=se.rmst,
       coef=intkmtimes[,3],var=diag(intkmtimes[,4]^2))
class(out) <- c("resmean_phreg")
return(out)
}# }}}

##' @export
vcov.resmean_phreg <- function(object,cause=1,...) 
{# {{{
if (is.na(match("rmean",names(object$intkmtimes)))) name <- paste("se.intF1",cause,sep="") else name <- "se.rmean"
return(diag(object$intkmtimes[,name]^2))
}# }}}

##' @export
coef.resmean_phreg <- function(object,cause=1,...) 
{# {{{
if (is.na(match("rmean",names(object$intkmtimes)))) name <- paste("intF1",cause,sep="") else name <- "rmean"
return(object$intkmtimes[,name])
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
    rownames(intF1mtimes) <- rep(x$strata.level,length(times)); 
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
	if (inherits(x,"resmean_phreg") & !is.null(x$intF1times)) ylab <- "Years lost up to t: t - E(min(T,t))"
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

##' IPTW Cox, Inverse Probaibilty of Treatment Weighted Cox regression 
##'
##' Fits Cox model with treatment weights \deqn{ w(A)= \sum_a I(A=a)/\pi(a|X)}, where
##'  \deqn{\pi(a|X)=P(A=a|X)}. Computes
##' standard errors via influence functions that are returned as the IID argument. 
##' Propensity scores are fitted using either logistic regression (glm) or the multinomial model (mlogit) when there are 
##' than treatment categories. The treatment needs to be a factor and is identified on the rhs
##' of the "treat.model". Recurrent events can be considered with start,stop structure and then cluster(id) must be
##' specified. Robust standard errors are computed in all cases. 
##'
##' Time-dependent propensity score weights can also be computed when treat.var is used, it must be 1 at the time
##' of first (A_0) and 2nd treatment (A_1), then uses weights \deqn{w_0(A_0) * w_1(A_1)^{t>T_r}} where \deqn{T_r} is
##' time of 2nd randomization.
##'
##' @param formula for phreg 
##' @param data data frame for risk averaging
##' @param treat.model propensity score model (binary or multinomial) 
##' @param treat.var a 1/0 variable that indicates when treatment is given and the propensity score is computed 
##' @param weights may be given, and then uses weights*w(A) as the weights
##' @param estpr (=1, default) to estimate propensity scores and get infuence function contribution to uncertainty
##' @param pi0 fixed simple weights 
##' @param se.cluster to compute GEE type standard errors when additional cluster structure is present
##' @param ...  arguments for phreg call
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data <- mets:::simLT(0.7,100,beta=0.3,betac=0,ce=1,betao=0.3)
##' dfactor(data) <- Z.f~Z
##' out <- phreg_IPTW(Surv(time,status)~Z.f,data=data,treat.model=Z.f~X)
##' summary(out)
##' @export
phreg_IPTW <- function (formula, data,treat.model = NULL, treat.var = NULL,weights = NULL, estpr = 1, pi0 = 0.5,se.cluster=NULL,...)
{# {{{
    cl <- match.call()
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster", "offset")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, c("Event", "Surv")))
        stop("Expected a 'Surv' or 'Event'-object")
    if (ncol(Y) == 2) {
        exit <- Y[, 1]
        entry <- NULL
        status <- Y[, 2]
    }
    else {
        entry <- Y[, 1]
        exit <- Y[, 2]
        status <- Y[, 3]
    }
    id <- strata <- NULL
    if (!is.null(attributes(Terms)$specials$cluster)) {
        ts <- survival::untangle.specials(Terms, "cluster")
        pos.cluster <- ts$terms
        Terms <- Terms[-ts$terms]
        id <- m[[ts$vars]]
    } else pos.cluster <- NULL

 call.id <- id;
 conid <- construct_id(id,length(exit))
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid

    if (!is.null(se.cluster)) {
        sds <- sort(unique(se.cluster))
        sid <- length(sds)
        if (is.numeric(se.cluster))
            se.cluster <- fast.approx(sds, se.cluster) - 1
        else se.cluster <- as.integer(factor(se.cluster, labels = seq(sid))) - 1
    }

    id <- id + 1
    nid <- length(unique(id))
    data$id__ <- id
    data$cid__ <- cumsumstrata(rep(1, length(id)), id - 1, nid)
    treats <- function(treatvar) {
        treatvar <- droplevels(treatvar)
        nlev <- nlevels(treatvar)
        nlevs <- levels(treatvar)
        ntreatvar <- as.numeric(treatvar)
        return(list(nlev = nlev, nlevs = nlevs, ntreatvar = ntreatvar))
    }
    fittreat <- function(treat.model, data, id, ntreatvar, nlev) {
        if (nlev == 2) {
            treat.model <- drop.specials(treat.model, "cluster")
            treat <- glm(treat.model, data, family = "binomial")
            iidalpha <- lava::iid(treat, id = id)
            lpa <- treat$linear.predictors
            pal <- expit(lpa)
            pal <- cbind(1 - pal, pal)
            ppp <- (pal/pal[, 1])
            spp <- 1/pal[, 1]
        }
        else {
            treat.modelid <- update.formula(treat.model, . ~ . + cluster(id__))
            treat <- mlogit(treat.modelid, data)
            iidalpha <- lava::iid(treat)
            pal <- predict(treat, data, se = 0, response = FALSE)
            ppp <- (pal/pal[, 1])
            spp <- 1/pal[, 1]
        }
        Xtreat <- model.matrix(treat.model, data)
        tvg2 <- 1 * (ntreatvar >= 2)
        pA <- c(mdi(pal, 1:length(ntreatvar), ntreatvar))
        pppy <- c(mdi(ppp, 1:length(ntreatvar), ntreatvar))
        Dppy <- (spp * tvg2 - pppy)
        Dp <- c()
        for (i in seq(nlev - 1)) Dp <- cbind(Dp, Xtreat * ppp[, i + 1] * Dppy/spp^2)
        DPai <- -1 * Dp/pA^2
        out <- list(iidalpha = iidalpha, pA = pA, Dp = Dp, pal = pal,
            ppp = ppp, spp = spp, id = id, DPai = DPai)
        return(out)
    }
    expit <- function(x) 1/(1 + exp(-x))
    if (!is.null(treat.var)) {
        weightWT <- data[, treat.var]
        whereW <- which(weightWT == 1)
        CountW <- cumsumstrata(weightWT, id - 1, nid)
        dataW <- data[whereW, ]
        idW <- id[whereW]
    } else {
        whereW <- 1:nrow(data)
        dataW <- data
        idW <- id
        CountW <- cumsumstrata(rep(1,nrow(data)), id-1,nid)
    }
    treat.name <- all.vars(treat.model)[1]
    treatvar <- dataW[, treat.name]
    if (!is.factor(treatvar))
        stop(paste("treatment=", treat.name, " must be coded as factor \n",
            sep = ""))
    treats <- treats(treatvar)
    wlPA <- ww <- rep(1, nrow(data))
    idWW <- mystrata2index(cbind(id, CountW))
    if (estpr[1] == 1) {
        fitt <- fittreat(treat.model, dataW, idW, treats$ntreatvar, treats$nlev)
        iidalpha0 <- fitt$iidalpha
        wPA <- c(fitt$pA)
        DPai <- fitt$DPai
    } else {
        if (length(pi0)==length(treats$ntreatvar)) wPA <- pi0 
	else wPA <- ifelse(treats$ntreatvar==2,pi0[1],1-pi0[1])
        ## pi0 <- rep(pi0, treats$nlev)
        DPai <- matrix(0, nrow(data), 1)
    }
    ww <- rep(1, nrow(data))
    ww[whereW] <- wPA
    wlPA <- exp(cumsumstrata(log(ww), idWW - 1, attr(idWW, "nlevel")))
    wwt <- c(exp(cumsumstrata(log(ww), id - 1, nid)))
    ## P(t) = P_0 * P_1^(I(t>T1)), time-dependent weights
    ## DP = P(t) \sum_j P_j I(t> TJ) (-DP_j/P_j^2)
    if (estpr[1] == 1) {
        DPait <- matrix(0, nrow(data), ncol(DPai))
        DPait[whereW, ] <- DPai
        DPait <- apply(DPait * c(wlPA), 2, cumsumstrata, id - 1, nid)/wwt
    } else DPait <- NULL
    if (is.null(weights)) ww <- 1/wwt else ww <- weights/wwt
    phw <- phreg(formula, data, weights = ww, Z = DPait, ...)

###    check.derivative <- 0
###	if (check.derivative == 1) { ## {{{ 
###	### for checking derivative 
###	fpar <- glm(treat.model,dataw,family=binomial)
###	mm <- model.matrix(treat.model,dataw)
###	cpar <- coef(fpar)
###   	 ### library(numderiv)
###
###	ff <- function(par,base=0) {
###	pa <-        expit(mm %*% par)
###	www <- ifelse(dataw[,treat.name] == "1", pa, 1 - pa)
###	ww <- rep(1, nrow(data))
###	ww[wherew] <- www
###	wlpa <- exp(cumsumstrata(log(ww), idww - 1, attr(idww, "nlevel")))
###	wwwt <- exp(cumsumstrata(log(ww), id - 1, nid))
###
###	pp <- phreg(formula,data,weights=1/wwwt,no.opt=true,beta=coef(phw))
###	if (base==1) po <- c(pp$cumhaz[,2]) else po <- pp$gradient
###	return(po)
###	}
###
###	print(ff(cpar))
###	gf <- jacobian(ff,cpar)
###	print(t(gf))
###        ###
###        print(ff(cpar,base=1))
###	gf <- jacobian(ff,cpar,base=1)
###	print(gf)
###	print("___________________________________"); 
###	} ## }}}

if (estpr[1] == 1) { 
	xx <- phw$cox.prep
	nid <- max(xx$id)
	S0i <- rep(0, length(xx$strata))
	wPAJ <- xx$weights[xx$jumps+1]*xx$caseweights[xx$jumps+1]
	Xt <- xx$X
	S0 <- phw$S0 * wPAJ
	S0i[xx$jumps + 1] <- 1/S0
	U <- E <- matrix(0, nrow(xx$X), phw$p)
	U[xx$jumps + 1, ] <- phw$U/wPAJ
	E[xx$jumps + 1, ] <- phw$E


        if (phw$p>0) {
	rr <- c(xx$sign * exp(Xt %*% coef(phw) + xx$offset) * xx$weights)
	rrnw <- c(xx$sign * exp(Xt %*% coef(phw) + xx$offset))
        } else {
	rr <- c(xx$sign * exp(xx$offset) * xx$weights)
	rrnw <- c(xx$sign * exp(xx$offset))
        }
	DWX = .Call("vecMatMat", xx$Z, Xt)$vXZ
	S1 = apply(Xt * rr, 2, revcumsumstrata, xx$strata, xx$nstrata)
	###S00 = revcumsumstrata( rr, xx$strata, xx$nstrata)
	DS1 = apply(DWX * rrnw, 2, revcumsumstrata, xx$strata, xx$nstrata)
	DS0 = apply(xx$Z * rrnw, 2, revcumsumstrata, xx$strata, xx$nstrata)
	if (phw$p>0) {
		DS0S1 = .Call("vecMatMat", DS0[xx$jumps + 1, , drop = FALSE],S1[xx$jumps + 1, , drop = FALSE])$vXZ
		DUa2 <- apply(wPAJ * DS1[xx$jumps + 1, , drop = FALSE]/c(S0),2, sum) - apply(wPAJ * DS0S1/c(S0^2), 2, sum)
		DUa2 <- matrix(DUa2, ncol(fitt$DPai), phw$p)
		DUa1 <- t(xx$Z[xx$jumps + 1, ]) %*% (phw$U/wPAJ)
		DUa <- DUa1 - DUa2
		iidpal <- iidalpha0 %*% DUa
		iidbeta <- lava::iid(phw) + iidpal %*% phw$ihess
		phw$DUa <- DUa
		phw$IID <- iidbeta
		phw$naive.var <- phw$var
		phw$var <- crossprod(iidbeta)
	} else iidbeta <- NULL

        ###
	DAw <- apply(xx$caseweights[xx$jumps+1]*xx$Z[xx$jumps + 1, ,drop=FALSE]/c(S0),2,cumsumstrata,phw$strata.jumps,xx$nstrata)
	DA2 <- apply(wPAJ*DS0[xx$jumps + 1, , drop = FALSE]/c(S0)^2,2,cumsumstrata,phw$strata.jumps,xx$nstrata)
	DAt <- 1*(DAw-DA2)
	phw <- robust.phreg(phw,beta.iid=iidbeta)
	varA <- phw$robse.cumhaz[,2]^2
	phw$naive.se.cumhaz <- phw$robse.cumhaz 
	vtheta <- crossprod(iidalpha0)
        varthetat <-   rowSums((DAt %*% vtheta)*DAt)
        ###
        ww <- xx$caseweights*xx$weights
	S0i2 <- rep(0,length(xx$strata))
        S0i2[xx$jumps+1] <- 1/(phw$S0^2*ww[xx$jumps+1])
        cumS0i2 <- c(cumsumstrata(S0i2,xx$strata,xx$nstrata))
        xxx <- (S0i-rr*cumS0i2)
	id <- xx$id
	mid <- max(xx$id)+1
        thetat <- iidalpha0[id+1,,drop=FALSE]
        covk1 <- apply(S0i*thetat,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
        covk2 <- apply(rr*thetat,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
        ###
        covk1 <- apply(S0i*thetat,2,cumsumstrata,xx$strata,xx$nstrata)
        covk2 <- apply(rr*thetat,2,revcumsumstrata,xx$strata,xx$nstrata)
        ###
        covk2 <- covk2*cumS0i2
        covv <- (covk1-covk2)[xx$jumps+1,,drop=FALSE]
        covvt <- 2*apply(covv*DAt,1,sum)
        varA <- varA+varthetat+2*apply(covv*DAt,1,sum)

	Ht <- NULL
	if (phw$p>0) {
		covbetatheta <- t(iidalpha0) %*%  iidbeta 
		Ht <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
		Ht <- Ht[xx$jumps+1,,drop=FALSE]
		covvBT <- 2* rowSums((Ht %*% t(covbetatheta))*DAt)
		varA <- varA + covvBT 
	}
	phw$robse.cumhaz[,2] <- varA^.5
	phw$se.cumhaz[,2] <- varA^.5
	phw$var.cumhaz <- cbind(phw$se.cumhaz[,1],varA)
	resAiid <- list(DAt=DAt,Ht=Ht,iidalpha0=iidalpha0,iidbeta=phw$IID)
	phw$resAiid <- resAiid
	phw$iptw <- ww
} else {
   phw$iptw <- ww
   phw <- robust.phreg(phw)
   phw$naive.se.cumhaz <- phw$se.cumhaz 
   phw$se.cumhaz <- phw$robse.cumhaz
   phw$naive.var <- phw$var
} 


if (is.matrix(phw$IID)) 
	if (nrow(phw$IID)==length(name.id)) rownames(phw$IID) <- name.id

if (!is.null(se.cluster)) {
	phw$IID.simple <- phw$IID
	phw$var.simple <- phw$var
	phw$IID <- apply(phw$IID,2,sumstrata,se.cluster,sid)
	phw$var <- crossprod(phw$IID)
	phw$se.cluster <- se.cluster
}
class(phw) <- c("phreg","IPTW")

return(phw)
}# }}}

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
##' @param same.data assumes that same data is used for fitting of survival model and averaging. 
##' @param id might be given to link to data to iid decomposition of survival data, must be coded as 1,2,..,  
##' @param subdata rows or TRUE/FALSE to select which part of the data that is used for the G-computation. Might be treated 
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
##'
##' ss <- phreg(Surv(time,event)~strata(tcell.f)+platelet+age,bmt) 
##' summary(survivalG(ss,bmt,50))
##'
##' sst <- survivalGtime(ss,bmt,n=50)
##' plot(sst)
##'
##' fg1t <- survivalGtime(fg1,bmt,n=50)
##' plot(fg1t)
##' 
##' @export
##' @aliases survivalGtime
survivalG <- function(x,data,time=NULL,Avalues=c(0,1),varname=NULL,same.data=TRUE,id=NULL,subdata=NULL)
{# {{{

if (is.null(time)) stop("Give time for estimation of survival/cumulative incidence\n")

if (inherits(x,c("cifreg","phreg","recreg")))
Aiid <- iidBaseline(x,time=time) else stop("Must be cifregFG/phreg/recreg object \n"); 

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

if (is.null(subdata))  subdata <- 1:length(x$id)
## for cluster case take first record for each subject
cid <- countID(data.frame(id=x$id[subdata]))
FirstId <- which(cid$Countid==1)
datA <- data[subdata,][FirstId,]
id.data <- x$id[subdata][FirstId]
formulaX <- update.formula(x$formula,.~.)
formulaX <- drop.specials(formulaX,"cluster")
datA <- dkeep(datA,x=all.vars(formulaX))
xlev <- lapply(datA,levels)

idname <- attr(formulaX,"variables")$cluster
##### to work with predict in case id used in object
datA[,idname] <- 1

cumhaz.time <- cpred(x$cumhaz,time)[-1,]
k <- 1; risks <- c(); DariskG <- list()
for (a in nlevs) { ## {{{
 datA[,treat.name] <- a

 pp <- predict(x,datA,time,se=0)
 Xbase <- 1*outer(pp$strata,0:(x$nstrata-1),"==")

 if (inherits(x,"phreg"))  { ps0 <- c(pp$surv); Dma  <- -cbind(c(pp$RR)*ps0*Xbase,c(ps0*pp$cumhaz)*pp$X);    }
 if (inherits(x,"cifreg")) { ps0 <- c(pp$cif); Dma  <- cbind(c(pp$RR)*c(1-ps0)*Xbase,c((1-ps0)*pp$cumhaz)*pp$X); }
  else if (inherits(x,"recreg")) { ps0 <- pp$cumhaz; Dma <-  cbind(c(pp$RR)*Xbase,c(pp$cumhaz)*pp$X) }

 risks <- cbind(risks,ps0)
 DariskG[[k]] <- apply(Dma,2,sum)
 k <- k+1
} ## }}}

predictAiid <- NULL
###
Grisk <- apply(risks,2,mean)
risk.iid  <- t(t(risks)-Grisk)
###
nid <- x$nid
ndata <- length(unique(x$id[subdata]))
risk.iid <- apply(risk.iid,2,sumstrata,id.data,nid)/ndata
coxiid <- cbind(Aiid$base.iid,Aiid$beta.iid)

if (same.data) {
   for (a in seq_along(nlevs)) risk.iid[,a] <- risk.iid[,a]+ coxiid %*% DariskG[[a]]/ndata
   vv <- crossprod(risk.iid)
} else {
   predictAiid <- matrix(0,nid,ncol(risks))
   for (a in seq_along(nlevs))  {
	risk.iid[,a] <- risk.iid[,a] 
	predictAiid[,a] <- coxiid %*% DariskG[[a]]/ndata
   }
   vv <- crossprod(risk.iid)+crossprod(predictAiid)
}

###estimate(lava::estimate(coef=theta,vcov=vv,f=function(p) Gf(p,ic=0))

if (inherits(x,"phreg"))  { 
out <- estimate(coef=1-Grisk,vcov=vv,labels=paste("risk",nlevs,sep=""))
sout <- estimate(coef=Grisk,vcov=vv,labels=paste("risk",nlevs,sep=""))
ed <- estimate(coef=Grisk,vcov=vv,f=function(p) (1-p[-1])-(1-p[1]))
rd <- estimate(coef=Grisk,vcov=vv,f=function(p) (1-p[-1])/(1-p[1]),null=1)
srd <- estimate(coef=Grisk,vcov=vv,f=function(p) (p[-1])/(p[1]),null=1)
out <- list(risk.iid=risk.iid,survivalG=sout,risk=out,difference=ed,ratio=rd,survival.ratio=srd,vcov=vv)
} 
if (inherits(x,"cifreg") | inherits(x,"recreg")) { 
out <- estimate(coef=Grisk,vcov=vv,labels=paste("risk",nlevs,sep=""))
ed <- estimate(coef=Grisk,vcov=vv,f=function(p) (p[-1])-(p[1]))
rd <- estimate(coef=Grisk,vcov=vv,f=function(p) (p[-1])/(p[1]),null=1)
out <- list(risk.iid=risk.iid,risk=out,difference=ed,ratio=rd,vcov=vv)
}

class(out) <- "survivalG"
return(out)
} ## }}}

##' @export
survivalGtime <- function(x,data,time=NULL,n=100,...)
{# {{{

if (!is.null(time)) if (time=="all") time <- x$cumhaz[,1] 

if (is.null(time)) {
       rr <- range(x$cumhaz[,1])
       time <- seq(rr[1],rr[2],length=n)
}

survivalG <- risk <- difference <- ratio <- survival.ratio <- c()

for (tt in time) {
  Gt <- survivalG(x,data,time=tt,...)
  strata <- strata(rownames(Gt$risk$coefmat))
  survivalG <- rbind(survivalG,cbind(tt,Gt$survivalG$coefmat))
  risk <- rbind(risk,cbind(tt,Gt$risk$coefmat))
  difference <- rbind(difference,cbind(tt,Gt$difference$coefmat))
  ratio <- rbind(ratio,cbind(tt,Gt$ratio$coefmat))
  survival.ratio <- rbind(survival.ratio,cbind(tt,Gt$survival.ratio$coefmat))
}
  colnames(survivalG)[1] <- "time"
  colnames(risk)[1] <- "time"
  colnames(difference)[1] <- "time"
  colnames(ratio)[1] <- "time"
  colnames(survival.ratio)[1] <- "time"
  strata <- strata(rownames(risk))
out <- list(time=time,survivalG=survivalG,risk=risk,difference=difference,
	    ratio=ratio,survival.ratio=survival.ratio,strata=strata)

class(out) <- "survivalGtime"
return(out)
} ## }}}

##' @export
plot.survivalGtime <- function(x,type=c("survival","risk","survival.ratio","difference","ratio"),ylim=NULL,legend=NULL,...) {# {{{

  ## to deal with fine-gray based things, that do not contain survival estimates
  if ((ncol(x$survivalG)==1) & type[1]=="survival") type <- "risk"
  us <- unique(x$strata)
  cols <- 1:length(us)
  ltys <- cols
  ss0 <- x$strata==us[1]
 
  if (type[1]=="survival") {
  if (is.null(ylim))  ylim  <- range(x$survivalG[,c(4,5)])
  plot(x$time,x$survivalG[ss0,2],type="s",xlab="time",ylab="Survival",col=cols[1],lty=ltys[1],ylim=ylim,...)
  plotConfRegion(x$time,x$survivalG[ss0,c(4,5)],col=cols[1])

  k <- 2
for (ss in us[-1]) {
   ss0 <- x$strata==ss
   lines(x$time,x$survivalG[ss0,2],type="s",ylim=c(0,1),col=cols[k],lty=ltys[k])
   plotConfRegion(x$time,x$survivalG[ss0,c(4,5)],col=cols[k])
   k <- k+1
}
  }

  if (type[1]=="risk") {
   if (is.null(ylim))  ylim  <- range(x$risk[,c(4,5)]) 
   plot(x$time,x$risk[ss0,2],type="s",xlab="time",ylab="risk",col=cols[1],lty=ltys[1],ylim=ylim,...)
   plotConfRegion(x$time,x$risk[ss0,c(4,5)],col=cols[1])
   k <- 2
   for (ss in us[-1]) {
   ss0 <- x$strata==ss
   lines(x$time,x$risk[ss0,2],type="s",ylim=c(0,1),col=cols[k],lty=ltys[k])
   plotConfRegion(x$time,x$risk[ss0,c(4,5)],col=cols[k])
   k <- k+1
   }
  }

  if (type[1]=="difference")  {
	  if (is.null(ylim))  ylim  <- range(x$risk[,c(4,5)]) 
  plot(x$time,x$difference[,2],type="s",xlab="time",ylab="difference in risk",ylim=range(x$difference[,2]),col=cols[1],lty=ltys[1],...)
  plotConfRegion(x$time,x$difference[,c(4,5)],col=cols[1])
  }
  if (type[1]=="ratio")  {
  if (is.null(ylim))  ylim  <- range(x$ratio[,c(4,5)]) 
  plot(x$time,x$ratio[,2],type="s",ylim=ylim,xlab="time",ylab="ratio",col=cols[1],lty=ltys[1],...)
  plotConfRegion(x$time,x$ratio[,c(4,5)],col=cols[1])
}
 if (type[1]=="survival.ratio")  {
   if (is.null(ylim))  ylim  <- range(x$survival.ratio[,c(4,5)]) 
  plot(x$time,x$survival.ratio[,2],type="s",xlab="time",ylab="risk",col=cols[1],lty=ltys[1],ylim=ylim,...)
  plotConfRegion(x$time,x$survival.ratio[,c(4,5)],col=cols[1])
}

 if (is.null(legend)) legend("topleft",levels(x$strata),col=cols,lty=ltys)

}
# }}}

###{{{ summary 

##' @export
print.survivalG <- function(x,...) {
  print(summary(x,...))
}

##' @export
summary.survivalG <- function(object,...) {
  res <- list(risk=object$risk,difference=object$difference,ratio=object$ratio,survival.ratio=object$survival.ratio)
  class(res) <- "summary.survivalG"
  res
}

##' @export
print.summary.survivalG  <- function(x,...) {
    cat("risk:\n")
    print(x[["risk"]],...)
    cat("\n")

    cat("Average Treatment effects (G-estimator) :\n")
    print(x$difference,...)
    cat("\n")

    cat("Average Treatment effect risk-ratio (G-estimator) :\n")
    print(x$ratio$coefmat,...)
    cat("\n")

    cat("Average Treatment effect (1-risk=survival)-ratio (G-estimator) :\n")
    print(x$survival.ratio$coefmat,...)
    cat("\n")

}

###}}} summary 

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
##' @param no.baseline to fit model without baseline hazard
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
aalenMets <- function(formula,data=data,no.baseline=FALSE,...)
{# {{{
formula.call <- formula
x <- phreg(formula,data=data,no.opt=TRUE,...)
xx <- x$cox.prep

### computation of intZHZ matrix, and gamma and baseline ##################
# {{{
eXb <- c(xx$sign) * c(xx$weights)
S0 = c(revcumsumstrata(eXb,xx$strata,xx$nstrata))
S0[S0==0] <- 1
E=apply(eXb*xx$X,2,revcumsumstrata,xx$strata,xx$nstrata)/S0; 
if (no.baseline) E <- 0*E
###if (!no.int) E <- cbind(1,E)
###if (no.int) X <-  xx$X else X <- cbind(1,xx$X)
###if (no.int) XX <-  xx$XX else { XX <- .Call("vecCPMat",X)$XX }
X <- xx$X; XX <- xx$XX
###
S2=apply(eXb*XX,2,revcumsumstrata,xx$strata,xx$nstrata)
E2  <- .Call("vecCPMat",E)$XX
###
dts <- c(diffstrata(xx$time,xx$strata,xx$nstrata))
dt <- diff(c(0,xx$time))
###
intZbar <- apply(E*dts,2,cumsumstrata,xx$strata,xx$nstrata) 
intZHZt <- apply((S2-S0*E2)*dts,2,cumsum) 
p <- ncol(E)
intZHZ <-  matrix(.Call("XXMatFULL",tail(intZHZt,1),p,PACKAGE="mets")$XXf,p,p)

IintZHZ  <-  solve(intZHZ)
intZHdN <- matrix(x$gradient,ncol(E),1)
XJ <- X[xx$jumps+1,]
if (no.baseline) intZHdN <- matrix(apply(XJ,2,sum),ncol(E),1)
gamma <- IintZHZ %*% intZHdN
###if (!no.int) nn <- c("int",colnames(x$X)) else 
nn <- colnames(x$X)
rownames(gamma)  <-  nn
x$cumhaz[,2] <- x$cumhaz[,2]- intZbar[xx$jumps+1,] %*% gamma
# }}}

### iid gamma #########################################
# {{{
id <- xx$id
if (no.baseline==FALSE) {
  mm <-  X * c(X %*% gamma) * c(xx$time)+ apply( E* c(E %*% gamma)*dts,2,cumsumstrata,xx$strata,xx$nstrata) -  X* c(intZbar %*% gamma) - c(X %*% gamma)* intZbar 
  mm <- c(xx$weights*xx$sign) * mm
  mm <- apply(mm,2,sumstrata,id,max(id)+1)
  ###
  MGt <- t(x$hessian %*% t(lava::iid(x))) + mm 
  XJ2l  <- .Call("vecCPMat",x$U)$XX
  XJ2l <- matrix(apply(XJ2l,2,sum),nrow=1)
  varmg <-  matrix(.Call("XXMatFULL",XJ2l,p,PACKAGE="mets")$XXf,p,p)
  varmg <- IintZHZ %*% varmg %*% IintZHZ
} else {
  mm <-  X * c(X %*% gamma) * c(xx$time)
  mm <- c(xx$weights*xx$sign) * mm
  mm <- apply(mm,2,sumstrata,id,max(id)+1)
  XJ2l  <- .Call("vecCPMat",XJ)$XX
  XJ2l <- matrix(apply(XJ2l,2,sum),nrow=1)
  varmg <-  matrix(.Call("XXMatFULL",XJ2l,p,PACKAGE="mets")$XXf,p,p)
  varmg <- IintZHZ %*% varmg %*% IintZHZ
  XdN  <- apply(XJ,2,sumstrata,id[xx$jumps+1],max(id)+1)
  MGt <-  XdN-mm
}
iid <-  MGt %*% IintZHZ
# }}}

### output  #########################################
# {{{
coef <- c(gamma)
names(coef) <- nn
x$coef <- coef
x$var <- crossprod(iid) 
x$varmg <- varmg
x$iid <- iid
x$intZHdN <- intZHdN
x$intZHZ  <-  intZHZ
x$formula <- formula.call
x$ihessian <- IintZHZ
## score of gamma 
x$gradient <- c(intZHdN-intZHZ %*% gamma)

x$no.opt <- FALSE
class(x) <- c(class(x),"aalenMets")
# }}}

return(x)
}# }}}

aalenMetsOld <- function(formula,data=data,...)
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
###E2=.Call("vecMatMat",E,E)$vXZ;  
E2  <- .Call("vecCPMat",E)$XX
###
dts <- c(diffstrata(xx$time,xx$strata,xx$nstrata))
dt <- diff(c(0,xx$time))
###
intZbar <- apply(E*dts,2,cumsumstrata,xx$strata,xx$nstrata) 
intZHZt <- apply((S2-S0*E2)*dts,2,cumsum) 
intZHZ <- matrix(0,ncol(E),ncol(E));
intZHZ[lower.tri(intZHZ,diag=TRUE)] <- tail(intZHZt,1)
intZHZ<- intZHZ+t(intZHZ)
diag(intZHZ) <- diag(intZHZ)/2
###intZHZ <- matrix(tail(intZHZt,1),ncol(E), ncol(E))
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
##' Robust variance is default variance and obtained from the predict call 
##' @param formula formula with 'Surv' 'Event' outcome 
##' @param data data frame
##' @param ... Additional arguments to phreg 
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' data(sTRACE)
##' sTRACE$cluster <- sample(1:100,500,replace=TRUE)
##' out1 <- km(Surv(time,status==9)~strata(vf,chf),data=sTRACE)
##' out2 <- km(Surv(time,status==9)~strata(vf,chf)+cluster(cluster),data=sTRACE)
##' 
##' summary(out1,times=1:3)
##' summary(out2,times=1:3)
##' 
##' par(mfrow=c(1,2))
##' plot(out1,se=TRUE)
##' plot(out2,se=TRUE)
##' @export
km <- function(formula,data=data,...)
{# {{{
### coxo <- phreg(formula,data=data,...)
### coxo <- robust.phreg(coxo)
### chaz <-     coxo$cumhaz[,2]
### time <-     coxo$cumhaz[,1]
### if (robust) std.err <-  coxo$robse.cumhaz[,2] else std.err <-  coxo$se.cumhaz[,2]
### strat <-    coxo$strata[coxo$jumps]
### S0i <- rep(0,length(coxo$S0))
### S0i[coxo$S0>0] <- 1/coxo$S0[coxo$S0>0]
### kmt <- exp(cumsumstrata(log(1-S0i),strat,coxo$nstrata))
### ### to use basehazplot.phreg
### res <- list(cumhaz=cbind(time,kmt),se.cumhaz=cbind(time,kmt*std.err),time=time,surv=kmt,
### strata=strat,nstrata=coxo$nstrata,jumps=1:length(kmt),strata.name=coxo$strata.name,strata.level=coxo$strata.level,
### model.frame=coxo$model.frame,formula=coxo$formula)
 
res <- phreg(formula,data=data,...)

rhs <- update(formula,-1~.)
varss <- all.vars(rhs)
## find all strata 
first <- c(headstrata(res$strata.call,res$nstrata))

ddf <- data[first,varss]
pres <- predict(res,ddf,robust=TRUE,...)
pres$formula <- formula
pres$call <- match.call()

attr(pres,"data") <- ddf
class(pres) <- c("km","predictphreg")
return(pres)
}# }}}

##' @export
summary.km <- function(object,times=NULL,type=c("cif","cumhaz","surv")[3],...) { ## {{{
   out <- summary.predictrecreg(object,times=times,type=type[1],...)
   return(out)
} ## }}}

##' @export
plot.km <- function(x,...) { ## {{{
   plot.predictphreg(x,...)
}# }}}

###predict.km <- function(object,newdata,...) { ## {{{
## take strata after readPhreg
### take relevant parts of prediction that is only for each strata
## out <- predict.phreg(object,newdata,se=se,conf.type=conf.type[1],...)
###}# }}}

##' Cumulative incidence with robust standard errors 
##'
##' Cumulative incidence with robust standard errors 
##' @param formula formula with 'Event' outcome and strata (only!)
##' @param data data frame
##' @param cause NULL looks at all, otherwise specify which cause to consider
##' @param cens.code censoring code "0" is default, and death is cens.code!=0
##' @param death.code alternative to cens.code give codes of death 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' data(bmt)
##' bmt$cluster <- sample(1:100,408,replace=TRUE)
##' out1 <- cif(Event(time,cause)~+1,data=bmt,cause=1)
##' out2 <- cif(Event(time,cause)~+1+cluster(cluster),data=bmt,cause=1)
##' 
##' par(mfrow=c(1,2))
##' plot(out1,se=TRUE)
##' plot(out2,se=TRUE)
##' @export
cif <- function(formula,data=data,cause=1,cens.code=0,death.code=NULL,...)
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

  statusE <- 1*(status %in% cause)
  if (is.null(death.code)) statusD <- 1*(!(status %in% cens.code)) else statusD <- 1*(status %in% death.code)
  data$statusE__ <- statusE
  data$statusD__ <- statusD
  data$strata__  <- strata

  ### setting up formulae for the two phreg (cause of interest and death)
  if (is.null(id.orig)) { 
     formid <- update.formula(formula,~.+cluster(id)) 
     data$id <- id
  } else formid <- formula
  tt <- terms(formid)
  tt <- delete.response(tt)
  formid <- formula(tt)

  if (ncol(Y)==3) {
     formE <- as.formula(paste("Surv(entry,exit,statusE__)~+1"))
     formD <- as.formula(paste("Surv(entry,exit,statusD__)~+1"))
  } else {
     formE <- as.formula(paste("Surv(exit,statusE__)~+1"))
     formD <- as.formula(paste("Surv(exit,statusD__)~+1"))
  }
 formE <- update.formula(formE,formid)
 formD <- update.formula(formD,formid)

  if (sum(statusE)==0) warning("No events of type 1\n"); 
  coxE <- phreg(formE,data=data,...)
  coxS <- phreg(formD,data=data,...)

  ### cif 
  if (sum(statusE)>0) cifo <- recurrentMarginalPhreg(coxE,coxS) else cifo <- coxE
  ## to work with predict function
  ##  cifo$no.opt <- TRUE

  ### to use basehazplot.phreg
  class(cifo) <- c("cif","phreg")
  attr(cifo,"cause") <- cause
  attr(cifo,"cens.code") <- cens.code
  attr(cifo,"death.code") <- death.code
  return(cifo)
}# }}}

##' @export
plot.cif <- function(x,se=FALSE,ylab=NULL,ylim=c(0,1),conf.type=c("log","plain"),...) { ## {{{
   if (inherits(x,"cif") & is.null(ylab)) ylab <- "Probability"
   baseplot(x,se=se,ylab=ylab,ylim=ylim,restrict="prob",conf.type=conf.type[1],...)
} ## }}}

##' @export
summary.cif <- function(object,se=FALSE,ylab=NULL,times=NULL,conf.type=c("log","plain"),...) { ## {{{
   if (inherits(object,"cif") & is.null(ylab)) ylab <- "Probability"
   out <- summaryRecurrentobject(object,name="probability",times=times,conf.type=conf.type[1],restrict="prob",...)
   return(out)
} ## }}}

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
   out$se.cumhaz <- NULL
   return(out)
}# }}}

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

