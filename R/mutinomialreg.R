##' Multinomial regression based on phreg regression
##'
##' Fits multinomial regression model 
##' \deqn{ P_i = \frac{ \exp( X^\beta_i ) }{ \sum_{j=1}^K \exp( X^\beta_j ) }} 
##' for \deqn{i=1,..,K}
##' where \deqn{\beta_1 = 0}, such that \deqn{\sum_j P_j = 1} using phreg function. 
##' Thefore the ratio \deqn{\frac{P_i}{P_1} = \exp( X^\beta_i )}
##'
##' Coefficients give log-Relative-Risk relative to baseline group (first level of factor, so that it can reset by relevel command).  
##' Standard errors computed based on sandwhich form \deqn{ DU^-1  \sum U_i^2 DU^-1}.  
##'
##' Can also get influence functions (possibly robust) via iid() function, response should be a factor. 
##'
##' Can fit cumulative odds model as a special case of interval.logitsurv.discrete
##'
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param weights for score equations 
##' @param offset offsets for partial likelihood 
##' @param fix.X to have same coefficients for all categories
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##'
##' data(bmt)
##' dfactor(bmt) <- cause1f~cause
##' drelevel(bmt,ref=3) <- cause3f~cause
##' dlevels(bmt)
##'
##' mreg <- mlogit(cause1f~+1,bmt)
##' summary(mreg)
##'
##' mreg <- mlogit(cause1f~tcell+platelet,bmt)
##' summary(mreg)
##' 
##' mreg3 <- mlogit(cause3f~tcell+platelet,bmt)
##' summary(mreg3)
##' 
##' ## inverse information standard errors 
##' estimate(coef=mreg3$coef,vcov=mreg3$II)
##' 
##' 
##' @export
mlogit <- function(formula,data,offset=NULL,weights=NULL,fix.X=FALSE,...)
{# {{{

  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
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
###  if (!is.null(intpos  <- attributes(Terms)$intercept))
###    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  res <- mlogit01(X,Y,id=id,strata=strata,offset=offset,weights=weights,strata.name=strata.name,fix.X=fix.X,...) ###,
  return(res)
}# }}}


mlogit01 <- function(X,Y,id=NULL,strata=NULL,offset=NULL,weights=NULL,
       strata.name=NULL,cumhaz=FALSE,
       beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,
       propodds=NULL,AddGam=NULL,case.weights=NULL,fix.X=FALSE,...) {# {{{
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(Y)))
  if (is.null(offset)) offset <- rep(0,length(Y)) 
  if (is.null(weights)) weights <- rep(1,length(Y)) 
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z
  if (is.null(case.weights)) case.weights <- rep(1,length(Y)) 

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(Y))-1; 
   ## orginal id coding into integers 
   id.orig <- id+1; 

  types <- unique(as.numeric(Y))
  nlev <- length(types)

  nX <- nrow(X)
  idrow <- rep(1:nX,each=nlev)
  X <- X[idrow,,drop=FALSE]
  px <- ncol(X)
  Y <- Y[idrow]
  id <- id[idrow]
  status <- rep(0,nrow(X))
  nY <- as.numeric(Y)

  refg <- 1  ### else refg <- match(ref,types)
  nrefs <- (1:nlev)[-refg]
  for (i in 1:nlev) status[nY==i] <- rep(((1:nlev)==i),sum(nY==i)/nlev)
  time <- id
  strat <- rep(1:nlev,nX)
  XX <- c()
  nn <- c()
  namesX <- colnames(X); 
  if (!fix.X) {
     for (i in nrefs) { XX <- cbind(XX,X*(strat==i)); 
                      nn<-c(nn, paste(namesX,i,sep=".")) 
     }  
     colnames(XX) <- nn; 
  } else {
     ### different intercepts
     for (i in nrefs) { XX <- cbind(XX,cbind(X[,1]*(strat==i))); 
                        nn<-c(nn, paste(namesX[1],i,sep=".")) 
     }  
     ### same covariate effects 
     XX<- cbind(XX,X[,-1]*(strat!=1))
     nn <- c(nn,colnames(X[,-1]))
     colnames(XX) <- nn; 
  }
  rownames(XX) <- NULL

  datph=data.frame(time=time,status=status,XX=XX,id=id,idrow=idrow)
  loffset <- offset[idrow]
  lweights<- weights[idrow]

###  print(list(...))
  res <- phreg(Surv(time,status)~XX+strata(idrow)+cluster(id),datph,weights=lweights,offset=loffset,...)

  res$px <- px
  res$nlev <- nlev
  return(res)
}# }}}


###predmlogit(mreg,bmt[1:2,])
###X <- cbind(1,0,0)
###predmlogit <- function(object,X)
###{# {{{
###
###  expit  <- function(z) 1/(1+exp(-z)) ## expit
###
###  refg <- 1  ### else refg <- match(ref,types)
###  nrefs <- (1:(object$nlev-1))
###  px <- ncol(X)
###  Xbeta <- c()
###  k <- 1
###  for (i in nrefs) { Xbeta <- cbind(Xbeta,X %*% object$coef[(1:px)+px*(i-1)]);  }
###  for (i in nrefs)  print(object$coef[(1:px)+px*(i-1)]); 
###  head(X)
###  Xbeta
###
###  return(res)
###}# }}}


