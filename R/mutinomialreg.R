##' Multinomial regression based on phreg regression
##'
##' returns also influence functions (possibly robust) 
##'
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##'
##' data(bmt)
##' mreg <- mlogit(cause~tcell+platlet,bmt)
##' summary(mreg,type="martingale")
##' 
##' 
##' @export
mlogit <- function(formula,data,...)
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
###  if (ncol(Y)==2) {
###    exit <- Y[,1]
###    entry <- NULL ## rep(0,nrow(Y))
###    status <- Y[,2]
###  } else {
###    entry <- Y[,1]
###    exit <- Y[,2]
###    status <- Y[,3]
###  }
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
###  if (!is.null(intpos  <- attributes(Terms)$intercept))
###    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  res <- mymlogit01(X,Y,id,strata,offset,weights,strata.name,...) ###,
###	   list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster))
  return(res)
}# }}}

mlogit01 <- function(X,Y,id=NULL,strata=NULL,offset=NULL,weights=NULL,
             strata.name=NULL,cumhaz=FALSE,
             beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
	     case.weights=NULL,...) {# {{{
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(Y)))
  if (is.null(strata)) { strata <- rep(0,length(Y)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }
  if (is.null(offset)) offset <- rep(0,length(Y)) 
  if (is.null(weights)) weights <- rep(1,length(Y)) 
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z

  ## possible casewights to use for bootstrapping and other things
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

  nlev <- length(unique(levels(Y)))

  nX <- nrow(X)
  id <- rep(1:nX,each=nlev)
  X <- X[id,,drop=FALSE]
  Y <- Y[id]
  status <- rep(0,nrow(X))
  nY <- as.numeric(Y)
  for (i in 1:nlev) status[nY==i] <- rep(((1:nlev)==i),sum(nY==i)/nlev)
  time <- id
  strat <- rep(1:nlev,nX)
  XX <- c()
  for (i in 1:(nlev-1)) XX <- cbind(XX,X*(strat==i))
  rownames(XX) <- NULL
  data=data.frame(time=time,status=status,XX=XX,id=id)
  system.time(
     res <- phreg(Surv(time,status)~XX+strata(id),data)
   )
  return(res)
}# }}}

