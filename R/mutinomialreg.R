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
##' library(mets)
##' data(bmt)
##' bmt$id <- sample(200,408,replace=TRUE)
##' dfactor(bmt) <- cause1f~cause
##' drelevel(bmt,ref=3) <- cause3f~cause
##' dlevels(bmt)
##'
##' mreg <- mlogit(cause1f~+1+cluster(id),bmt)
##' summary(mreg)
##' head(iid(mreg))
##' dim(iid(mreg))
##'
##' mreg <- mlogit(cause1f~tcell+platelet,bmt)
##' summary(mreg)
##' head(iid(mreg))
##' dim(iid(mreg))
##' 
##' mreg3 <- mlogit(cause3f~tcell+platelet,bmt)
##' summary(mreg3)
##' 
##' ## inverse information standard errors 
##' lava::estimate(coef=mreg3$coef,vcov=mreg3$II)
##' 
##' ## predictions based on seen response or not 
##' ## all probabilities
##' head(predict(mreg,response=FALSE))
##' head(predict(mreg))
##' ## using newdata 
##' newdata <- data.frame(tcell=c(1,1,1),platelet=c(0,1,1),cause1f=c("2","2","0"))
##' ## only probability of seen response 
##' predict(mreg,newdata)
##' ## without response
##' predict(mreg,newdata,response=FALSE)
##' ## given indexx of P(Y=j)
##' predict(mreg,newdata,Y=c(1,2,3))
##' ##  reponse not given 
##' newdata <- data.frame(tcell=c(1,1,1),platelet=c(0,1,1))
##' predict(mreg,newdata)
##' @export
##' @aliases predict 
mlogit <- function(formula,data,offset=NULL,weights=NULL,fix.X=FALSE,...)
{# {{{
  cl <- match.call()
    m <- match.call(expand.dots = TRUE)[1:3]
    des <- proc_design(
        formula,
        data = data,
        specials = c("offset", "weights", "cluster"),
        intercept = TRUE
    )
    Y <- des$y
    X <- des$x
    des.weights <- des$weights
    des.offset  <- des$offset
    id      <- des$cluster

 ## take offset and weight first from formula, but then from arguments
  if (is.null(des.offset)) {
	  if (is.null(offset)) offset <- rep(0,nrow(X)) 
  } else offset <- des.offset
  if (is.null(des.weights)) {
	  if (is.null(weights)) weights <- rep(1,nrow(X)) 
  } else weights <- des.weights

   if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
   strata.name <- NULL

  res <- mlogit01(X,Y,id=id,strata=strata,offset=offset,weights=weights,strata.name=strata.name,
		  fix.X=fix.X,formula.call=formula,...) 
  res$design <- des
  return(res)
}# }}}

mlogit01 <- function(X,Y,id=NULL,strata=NULL,offset=NULL,weights=NULL,
       strata.name=NULL,cumhaz=FALSE,
       beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,
       propodds=NULL,AddGam=NULL,case.weights=NULL,fix.X=FALSE,formula.call=NULL,
       X.call=NULL,Y.call=NULL,...) {# {{{
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(Y)))
  if (is.null(offset)) offset <- rep(0,length(Y)) 
  if (is.null(weights)) weights <- rep(1,length(Y)) 
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z
  if (is.null(case.weights)) case.weights <- rep(1,length(Y)) 

  call.id <- id
  if (is.null(id)) id <- 1:nrow(X)

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

  res <- phreg(Surv(time,status)~XX+strata(idrow)+cluster(id),datph,weights=lweights,offset=loffset,...)
  res$formula <- formula.call
  res$px <- px
  res$nlev <- nlev
  class(res) <- rev(c("phreg","mlogit"))
  return(res)
}# }}}

##' @export
predict.mlogit <- function (object, newdata , se = TRUE, response=TRUE , Y=NULL,level=0.95,...)
{# {{{

  ylev <- levels(object$design$y)
  if (missing(newdata)) {
     X <- object$design$x
     ## take response from newdata if it is there and it is not given  
     if (response) Y <- as.numeric(object$design$y) 
  }  else {
     respindata <- length((grep(all.vars(object$formula)[1],names(newdata))))
     if (respindata>0 & response) x <- update_design(object$design,data = newdata,response=TRUE)
     else x <- update_design(object$design,data = newdata)
     X <- x$x
     ## take response from newdata if it is there and it is not given  
     if (is.null(Y))  { 
         if (response & !is.null(x$y))   
	   Y <- as.numeric(factor(x$y,levels=ylev))
     } else {
	     if (!is.numeric(Y)) 
		     Y <- as.numeric(factor(Y,levels=ylev))
     }
  }

  refg <- 1  ### else refg <- match(ref,types)
  nrefs <- (1:(object$nlev-1))
  px <- ncol(X)
  Xbeta <- c()
  for (i in nrefs) { Xbeta <- cbind(Xbeta,X %*% object$coef[(1:px)+px*(i-1)]);  }

  ppp <- cbind(1,exp(Xbeta))
  spp <- apply(ppp,1,sum)
  pp <- ppp/spp
  colnames(pp) <- ylev
  alpha <- 1-level
  crit <- -qnorm(1-alpha/2)

  if (!is.null(Y)) {
	  Yg2 <- 1*(Y>=2)
          pp <- p <- c(mdi(pp,1:length(Y),Y)) 
          pppy <- c(mdi(ppp,1:length(Y),Y)) 
     if (se) {
	     Dppy <-  (spp*Yg2-pppy) 
	     Dp <- c()
             for (i in nrefs) Dp <- cbind(Dp,X*ppp[,i+1]*Dppy/spp^2);  
             if (is.null(object$var)) covv <- vcov(object) else covv <- object$var
	     se <-  apply((Dp %*% covv) * Dp,1,sum)^.5
	     cmat <- data.frame(pred = p, se = se, lower = p - crit * se, upper = p + crit * se)
             names(cmat)[1:4] <- c("pred", "se", "lower", "upper")
             pp <- cmat
     }
 }

  return(pp)
}# }}}

