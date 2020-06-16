##' CIF regression
##'
##' CIF logistic for propodds=1 default
##' CIF Fine-Gray (cloglog) regression for propodds=NULL
##'
##' For FG model: 
##' \deqn{
##' \int (X - E ) Y_1(t) w(t) dM_1
##' }
##' is computed and summed over clusters  and returned multiplied with inverse 
##' of second derivative as iid.naive. Where \deqn{w(t) = G(t) (I(T_i \wedge t < C_i)/G_c(T_i \wedge t))} and
##' \deqn{E(t) = S_1(t)/S_0(t)} and \deqn{S_j(t) = \sum X_i^j Y_{i1}(t) w_i(t) \exp(X_i^T \beta)}
##'
##'
##' The iid decomposition of the beta's, however, also have a censoring term that is also
##' is computed and added to UUiid (still scaled with inverse second derivative)
##' \deqn{
##' \int (X - E ) Y_1(t) w(t) dM_1 + \int q(s)/p(s) dM_c 
##' }
##' and returned as iid 
##'
##' For logistic link standard errors are slightly to small since uncertainty from recursive baseline is not considered, so for smaller
##' data-sets it is recommended to use the prop.odds.subdist of timereg that is also more efficient due to use of different weights for 
##' the estimating equations. Alternatively, one can also bootstrap the standard errors. 
##'
##' @param formula formula with 'Event' outcome 
##' @param data data frame
##' @param cause of interest 
##' @param cens.code code of censoring 
##' @param cens.model for stratified Cox model without covariates 
##' @param offset offsets for cox model
##' @param weights weights for Cox score equations
##' @param Gc censoring weights for time argument, default is to calculate these with a Kaplan-Meier estimator, should then give G_c(T_i-) 
##' @param propodds 1 is logistic model, NULL is fine-gray model 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ## data with no ties
##' data(bmt,package="timereg")
##' bmt$time <- bmt$time+runif(nrow(bmt))*0.01
##' bmt$id <- 1:nrow(bmt)
##' 
##' ## logistic link  OR interpretation
##' ll=cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
##' plot(ll)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pll <- predict(ll,nd)
##' plot(pll)
##' 
##' ## Fine-Gray model
##' llfg=cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,propodds=NULL)
##' plot(llfg)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pll <- predict(ll,nd)
##' plot(pll)
##' 
##' sllfg=cifreg(Event(time,cause)~strata(tcell)+platelet+age,data=bmt,cause=1,propodds=NULL)
##' plot(sllfg)
##' @aliases revcumsum2strata vecAllStrata
##' @export
cifreg <- function(formula,data=data,cause=1,cens.code=0,cens.model=~1,
			weights=NULL,offset=NULL,Gc=NULL,propodds=1,...)
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
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

###  if (!is.null(id)) {
###	  ids <- sort(unique(id))
###	  nid <- length(ids)
###      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
###      id <- as.integer(factor(id,labels=seq(nid)))-1
###     }
###   } else id <- as.integer(seq_along(exit))-1; 
###  ### id from call coded as numeric 1 -> 
###  id.orig <- id; 

# }}}

 res <- c(cifreg01(data,X,exit,status,id,strata,offset,weights,strata.name,
		   cens.model=cens.model,
	   cause=cause,cens.code=cens.code,Gc=Gc,propodds=propodds,...),
   list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster)
   )

  class(res) <- c("phreg","cif.reg")
  return(res)
}# }}}

cifreg01 <- function(data,X,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
             strata.name=NULL,beta,stderr=TRUE,method="NR",no.opt=FALSE,propodds=1,profile=0,
   case.weights=NULL,cause=1,cens.code=0,Gc=NULL,cens.model=~+1,augmentation=0,...) {# {{{

##  setting up weights, strata, beta and so forth before the action starts# {{{
 p <- ncol(X)
 if (missing(beta)) beta <- rep(0,p)
 if (p==0) X <- cbind(rep(0,length(exit)))

 cause.jumps <- which(status==cause)
 max.jump <- max(exit[cause.jumps])
 other <- which((!(status %in% c(cens.code,cause)) ) & (exit< max.jump))

 entry <- NULL
 n <- length(exit)
 trunc <- (!is.null(entry))
  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }

 if (!trunc) entry <- rep(0,length(exit))
 if (is.null(offset)) offset <- rep(0,length(exit)) 
 if (is.null(weights)) weights <- rep(1,length(exit)) 
 if (is.null(case.weights)) case.weights <- rep(1,length(exit)) 
 strata.call <- strata
 Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1; 
   ## orginal id coding into integers 1:...
   id.orig <- id+1; 

# }}}


 ### censoring weights constructed
 whereC <- which(status==cens.code)
 time <- exit
 statusC <- (status==cens.code) 
 data$id <- id
 data$exit <- exit
 data$statusC <- statusC 
 cens.strata <- cens.nstrata <- NULL 

 if (is.null(Gc)) {
      kmt <- TRUE
      if (class(cens.model)[1]=="formula")
      {
         formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id))
         cens.model <- phreg(formC,data)
      }
      if (cens.model$p>0) kmt <- FALSE
      Pcens.model <- predict(cens.model,data,times=exit,tminus=TRUE,individual.time=TRUE,se=FALSE,km=kmt)
      Stime <- Pcens.model$surv <- c(Pcens.model$surv)
      ## strata from original data 
      nCstrata <- cens.model$nstrata
  } else  {
      formC <- NULL
      Stime <- Gc
      Pcens.model <- list(time=exit,surv=Gc,strata=0)
      nCstrata <- 1
 }

 ## setting up all jumps of type "cause", need S0, S1, S2 at jumps of "cause"
 stat1 <- 1*(status==cause)
 trunc <- FALSE
 xx2 <- .Call("FastCoxPrepStrata",entry,exit,stat1,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")
 xx2$nstrata <- nstrata
 jumps <- xx2$jumps+1
 jumptimes <- xx2$time[jumps]
 strata1jumptimes <- xx2$strata[jumps]
 Xj <- xx2$X[jumps,,drop=FALSE]

 ## G(T_j-) at jumps of type "cause"
 if (length(whereC)>0) {
    if (is.null(Gc)) {
        whereaJ <- fast.approx(c(0,cens.model$cumhaz[,1]),jumptimes,type="left")
        Gts <- vecAllStrata(cens.model$cumhaz[,2],cens.model$strata.jump,cens.model$nstrata)
	### back to km product-limit form
        Gts <- apply(rbind(0,Gts),2,diff)
	### back to km 
	GtsAll <- Gts <- apply(Gts,2,function(x) exp(cumsum(log(1-x))))
	Gts <- rbind(1,Gts)[whereaJ,]
	Gts[is.na(Gts)] <- 0
	Gjumps <- Gts
    } else Gts <- Gjumps <- c(1,Pcens.model$surv)[fast.approx(c(0,Pcens.model$time),jumptimes)]
 } else {
       Gts <-   Gjumps <- rep(1,length(jumptimes))
 }

 ## computing terms for those experiencing another cause, need S0, S1, S2
 if (length(other)>=1) {# {{{
	 trunc <- TRUE
	 weightso <- weights[other]/Stime[other]
	 timeoo <- rep(max(exit)+1,length(other))
	 statuso <- rep(1,length(other))
	 Xo <- X[other,,drop=FALSE]
	 offseto <- offset[other]
	 entryo <- exit[other]
	 ido <- id[other]
	 stratao <- strata[other]
	 ###
	 if (nCstrata>1) {
	    Cstratao <- Pcens.model$strata[other]
	    Zcall <- matrix(Cstratao,length(other),1)  
	 }  else { 
	    Cstratao <- rep(0,length(other))
	    Zcall <- matrix(0,1,1); 
	 }
	 xx <- .Call("FastCoxPrepStrata",entryo,timeoo,statuso,Xo,
	        ido,trunc,stratao,weightso,offseto,Zcall,case.weights[other],PACKAGE="mets")
	 xx$nstrata <- nstrata

	 timeo  <- xx$time
	 if (nCstrata>1) xxCstrata <- c(xx$Z) else xxCstrata <- rep(0,length(timeo))
	 ## use right because we want S_0(T_jump) 
	 where <- indexstrata(timeo,xx$strata,jumptimes,strata1jumptimes,nstrata,type="right")
 }# }}}

obj <- function(pp,all=FALSE) {# {{{

if (length(other)>=1)  {
if (nCstrata==1) {# {{{
	rr <- c(xx$sign*exp(xx$X %*% pp + xx$offset)*xx$weights)
	S0no <- revcumsumstrata(rr,xx$strata,xx$nstrata)
	S1no  <- apply(xx$X*rr,2,revcumsumstrata,xx$strata,xx$nstrata); 
        S2no  <- apply(xx$XX*rr,2,revcumsumstrata,xx$strata,xx$nstrata); 
	Gjumps <- c(Gjumps)

	S0no <- Gjumps*S0no[where]
	S1no <- Gjumps*S1no[where,,drop=FALSE]
	S2no <- Gjumps*S2no[where,,drop=FALSE]
# }}}
}  else {# {{{

	ffp <- function(x,strata,nstrata,strata2,nstrata2)
	{# {{{
	 x <- revcumsum2strata(x,strata,nstrata,strata2,nstrata2)$mres
	 ### take relevant S0sc (s=strata,c=cstrata) at jumptimes so that strata=s also match
	 print(round(cbind(x,xx$time,strata,nstrata,strata2,nstrata2),2))
	 x <- rbind(0,x)[where,]
	 print(round(x,2))
###	 print("S0no")
###	 print(cbind(jumptimes,x,strata1jumptimes))
	 ### average over Gc
	 x <- apply(x*Gts,1,sum)
	 return(x)
	}# }}}

	ff <- function(x,strata,nstrata,strata2,nstrata2)
	{# {{{
	 x <- revcumsum2strata(x,strata,nstrata,strata2,nstrata2)$mres
	 ### take relevant S0sc (s=strata,c=cstrata) at jumptimes so that strata=s also match
	 x <- x[where,]
	 x <- apply(x*Gts,1,sum)
	 return(x)
	}# }}}

	rr <- c(xx$sign*exp(xx$X %*% pp + xx$offset)*xx$weights)
	S0no  <- ff(rr,xx$strata,xx$nstrata,xxCstrata,nCstrata)
	S1no  <- apply(xx$X*rr,2,ff,xx$strata,xx$nstrata,xxCstrata,nCstrata); 
	S2no  <- apply(xx$XX*rr,2,ff,xx$strata,xx$nstrata,xxCstrata,nCstrata); 

}# }}}
} else { Gjumps <- S0no <- S1no <-  S2no <- 0} 

rr2 <- c(xx2$sign*exp(xx2$X %*% pp + xx2$offset)*xx2$weights)
rr2now <- c(xx2$sign*exp(xx2$X %*% pp + xx2$offset))
S0oo <- revcumsumstrata(rr2,xx2$strata,xx2$nstrata)
S1oo  <- apply(xx2$X*rr2,2,revcumsumstrata,xx2$strata,xx2$nstrata); 
S2oo  <- apply(xx2$XX*rr2,2,revcumsumstrata,xx2$strata,xx2$nstrata); 
S0oo <- S0oo[jumps,]
S1oo <- S1oo[jumps,,drop=FALSE]
S2oo <- S2oo[jumps,,drop=FALSE]

S0 <- c(S0oo+S0no)
S1 <- S1oo+S1no 

E <- S1/S0
weightsJ <- xx2$weights[jumps]
caseweightsJ <- xx2$caseweights[jumps]
strataJ <- xx2$strata[jumps]
rr2now <- rr2now[jumps]
U <- (Xj-E)
ploglik <- (log(rr2now)-log(S0))*weightsJ*caseweightsJ; 

###for (i in (0:(xx2$nstrata-1))) {
###print("S0 ========================================")
###print(i)
###print(S0oo[strataJ==i])
###print(i)
###print(S0no[strataJ==i])
###}
###for (i in (0:(xx2$nstrata-1))) {
###print("S1 ========================================")
###print(i)
###print(S1oo[strataJ==i,])
###print(S1no[strataJ==i,])
###}
###for (i in (0:(xx2$nstrata-1))) {
###print("S2 ========================================")
###print(i)
###print(S2oo[strataJ==i,])
###print(S2no[strataJ==i,])
###ll <- matrix(apply(S2oo,2,sum),p,p)
###print(ll)
###ll <- matrix(apply(S2no,2,sum),p,p)
###print(ll)
###}
###for (i in (0:(xx2$nstrata-1))) {
###print("E ========================================")
###print(i)
###print(E[strataJ==i,])
###}

if (!is.null(propodds)) {
   pow <- c(.Call("cumsumstrataPOR",weightsJ,S0,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$pow); 
   DLam <-.Call("DLambetaR",weightsJ,S0,E,Xj,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$res; 
   Dwbeta <- DLam*rr2now+(pow-1)*Xj
   DUadj  <- .Call("vecMatMat",Dwbeta,U,PACKAGE="mets")$vXZ 
}

Ut <- caseweightsJ*weightsJ*U
## E^2, as n x (pxp)
Et2 <-  .Call("vecMatMat",E,E,PACKAGE="mets")$vXZ
S2S0 <-  (S2oo+S2no)/S0
DUt <-  -(S2S0-Et2)

if (!is.null(propodds)) {
	Ut  <- pow*Ut
	S0 <- S0/pow
	DUt <- pow*DUt
	DUt <- DUt+DUadj
	if (profile==1) {
		Ut <- Ut+c(ploglik)*Dwbeta
		## not implemented 
		DUt <- DUt 
	}
	ploglik <- pow*ploglik
}

U  <- apply(Ut,2,sum)
DUt <- caseweightsJ*weightsJ*DUt
DU <- -matrix(apply(DUt,2,sum),p,p)
###if (!is.null(Saug)) {
###   S1ss <-  .Call("vecMatMat",S1a,S1aug,PACKAGE="mets")$vXZ
###   corDU <- matrix(apply(caseweightsJ*weightsJ*S1ss/S0^2,2,sum),p,p)
###   DU <- DU-corDU
###}

ploglik <- sum(ploglik)

U <- U+augmentation

out <- list(ploglik=ploglik,gradient=U,hessian=-DU,cox.prep=xx2,
	    hessiantime=DUt,weightsJ=weightsJ,caseweightsJ=caseweightsJ,
	   jumptimes=jumptimes,strata=strataJ,nstrata=nstrata,
	   time=jumptimes,S0=S0/(caseweightsJ*weightsJ),S2S0=S2S0,E=E,U=Ut,X=Xj,Gjumps=Gjumps)


if (all) return(out) else with(out,structure(-ploglik, gradient=-gradient, hessian=-hessian))
}# }}}

 if (length(jumps)==0) no.opt <- TRUE

 opt <- NULL
  if (p>0) {# {{{
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
          opt <- lava::NR(beta,obj,...)
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
  }# }}}

### opt <- lava::NR(beta,obj); beta.s <- opt$par 
 beta.s <- val$coef 
 ## getting final S's 
 opt <-  obj(beta.s,all=TRUE)

 ### iid version 
# {{{
 ###  iid.phreg
 ##iid robust phreg
 S0i <- rep(0,length(xx2$strata))
 S0i[jumps] <- 1/opt$S0
 Z <- xx2$X
 U <- E <- matrix(0,nrow(Z),p)
 E[jumps,] <- opt$E
 U[jumps,] <- opt$U
 cumhazA <- cumsumstratasum(S0i,xx2$strata,xx2$nstrata,type="all")
 cumhaz <- c(cumhazA$sum)
 rr <- c(xx2$sign*exp(Z %*% beta.s + xx2$offset))
 if (!is.null(propodds)) {
    cumhazm <- c(cumhazA$lagsum)
    S0star <- cumsumstrata(rr/(1+rr*cumhazm),xx2$strata,xx2$nstrata)
 }
 EdLam0 <- apply(E*S0i,2,cumsumstrata,xx2$strata,xx2$nstrata)

 ### Martingale  as a function of time and for all subjects to handle strata 
 MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx2$weights)
 mid <- max(xx2$id)
 UU <- apply(MGt,2,sumstrata,xx2$id,mid+1)

### for (i in 0:(xx2$nstrata-1)) {
### print(paste("====================",i))
### print(cbind(U,cumhaz,rr,xx2$weights,xx2$id)[xx2$strata==i,])
### print(cbind(Z,EdLam0,rr,xx2$weights,xx2$id)[xx2$strata==i,])
### print(MGt[xx2$strata==i,])
### print(UU[xx2$strata==i,])
### print(xx2$id)
### }

 if (length(other)>=1) {
   ### T_j for jumps of other type 
   Gt <- matrix(0,length(xx2$strata),nCstrata)
   Gt[jumps,] <- Gts
   #### matdoubleindex(Gts,1:nrow(Gts),xxCstrata)
   MGtGtx   <- matrix(0,length(other),ncol(E))
   for (i in 0:(nCstrata-1)) {
   rcumhazGt <- apply(S0i*Gt[,i+1,drop=FALSE],2,revcumsumstrata,xx2$strata,xx2$nstrata)
   rEdLam0Gt <- apply(E*S0i*Gt[,i+1],2,revcumsumstrata,xx2$strata,xx2$nstrata)
   ###
   ### type right since all jumps after entryo should be used 
   where <- indexstrata(c(rep(0,xx2$nstrata),xx2$time), c(0:(xx2$nstrata-1),xx2$strata),
			entryo[Cstratao==i],stratao[Cstratao==i],xx2$nstrata,type="right")
   rcumhazGtx <- c(0,rcumhazGt)[where-xx2$nstrata+1]
   rEdLam0Gtx <- rbind(0,rEdLam0Gt)[where-xx2$nstrata+1,]  
###   print(dim(rEdLam0Gtx))
###   print(table(Cstratao))
   rrx <- c(exp(Xo[Cstratao==i,,drop=FALSE] %*% beta.s + offseto[Cstratao==i])*weightso[Cstratao==i])
###   print(rrx)
   if (!is.null(propodds)) { 
   }
   ###	
   MGtGtx[Cstratao==i,] <- -(Xo[Cstratao==i,]*rcumhazGtx-rEdLam0Gtx)*rrx
   }

### for (i in 0:(xx2$nstrata-1)) 
###  print(cbind(MGtGtx,ido,stratao)[stratao==i,])
   UU2 <- apply(MGtGtx,2,sumstrata,ido,mid+1)

### for (i in 0:(xx2$nstrata-1)) {
### print(paste("====================",i))
###### print(cbind(U,cumhaz,rr,xx2$weights,xx2$id)[xx2$strata==i,])
###### print(cbind(Z,EdLam0,rr,xx2$weights,xx2$id)[xx2$strata==i,])
### print(cbind(Xo[Cstratao==i,],rcumhazGtx,rrx,ido)[stratao==i,])
### print(cbind(rEdLam0Gtx,rrx,ido)[stratao==i,])
### print(MGtGtx[stratao==i,])
### print(UU2[stratao==i,])
###### print(xx2$id)
### }

   UU  <-  UU+UU2
   }
# }}}

if ((length(other)>=1) & (length(whereC)>0) & (nCstrata==1)) {
 ### Censoring adjustment for jumps of other type but only for KM-case {{{
 where <- fast.approx(xx2$time,entryo,type="right")
 rrrx <- rep(0,length(xx2$strata))
 rrrx[where] <- rrx
 Xos <- matrix(0,length(xx2$time),ncol(Xo));
 Xos[where,] <- Xo
 ###
 Xos <- apply(Xos,2,cumsum)
 rro <- cumsum(rrrx)
 ###
 q <- -(Xos*c(rcumhazGt)-rEdLam0Gt*rro)

 idloc__ <- id
 cens.mgs = phreg(Surv(exit,status==0)~+cluster(idloc__),data=data,no.opt=TRUE)
 cxx <- cens.mgs$cox.prep
 ###
 Gt <- S0i <-  S0i2 <- rep(0,length(cxx$strata))
 S0i[cxx$jumps+1] <- 1/cens.mgs$S0
 S0i2[cxx$jumps+1] <- 1/cens.mgs$S0^2
 qc <- matrix(0,nrow(q),ncol(q)) 
 ## sort q after censoring times
 qc[cxx$jumps+1,] <- q[cxx$jumps+1]
 ###
 EdLam0q <- apply(qc*S0i2,2,cumsumstrata,cxx$strata,cxx$nstrata)
 ### Martingale  as a function of time and for all subjects to handle strata 
 MGc <- qc[,drop=FALSE]*S0i-EdLam0q
 MGc <- apply(MGc,2,sumstrata,cxx$id,mid+1)
### print(crossprod(MGc))
### print(crossprod(UU))
### print(t(UU) %*% MGc )
# }}}
} else MGc <- 0

 iH <- - tryCatch(solve(opt$hessian),error=function(e) matrix(0,nrow(opt$hessian),ncol(opt$hessian)) )
 Uiid <-  (UU+MGc) %*% iH
 UUiid <- UU %*% iH
 var1 <-  crossprod(UUiid) 
 varm <-  crossprod(Uiid) 
 strata <- xx2$strata[jumps]
 cumhaz <- cbind(opt$time,cumsumstrata(1/opt$S0,strata,nstrata))
 colnames(cumhaz)    <- c("time","cumhaz")

 ## SE of estimator ignoring some censoring terms 
 if (no.opt==FALSE & p!=0) { 
     DLambeta.t <- apply(opt$E/c(opt$S0),2,cumsumstrata,strata,nstrata)
     varbetat <-   rowSums((DLambeta.t %*% iH)*DLambeta.t)
     ### covariance is 0 for cox model 
     ### covv <-  apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
 } else varbetat <- 0
 var.cumhaz <- cumsumstrata(1/opt$S0^2,strata,nstrata)+varbetat
 se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)
 colnames(se.cumhaz) <- c("time","se.cumhaz")


out <- list(coef=beta.s,var=varm,se.coef=diag(varm)^.5,iid.naive=UUiid,
	    iid=Uiid,
	    ihessian=iH,hessian=opt$hessian,var1=var1,se1.coef=diag(var1)^.5,
	    ploglik=opt$ploglik,gradient=opt$gradient,
	    cumhaz=cumhaz, se.cumhaz=se.cumhaz,MGciid=MGc,
	    strata=xx2$strata,nstrata=nstrata,strata.name=strata.name,
	    strata.level=strata.level,propodds=propodds,
	    S0=opt$S0,E=opt$E,S2S0=opt$S2S0,time=opt$time,Ut=opt$U,
            jumps=jumps,II=iH,exit=exit,p=p,opt=opt,n=nrow(X),nevent=length(jumps),
	    Pcens.model=Pcens.model,Gjumps=Gjumps
	    )

return(out)
}# }}}

##' @export
S0_FG_Gct <- function(S0,Gct,strata,nstrata,strata2,nstrata2)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (any(strata2<0) | any(strata2>nstrata2-1)) stop("strata2 index not ok\n"); 
if (length(S0)!=length(strata))  stop("length of x and strata must be same\n"); 
if (length(S0)!=length(strata2)) stop("length of x and strata2 must be same\n"); 
if (length(Gct)!=length(S0)) stop("length of S0 and Gct must be same\n"); 
res <- .Call("S0_FG_GcR",as.double(S0),as.double(Gct),
	     strata,nstrata,strata2,nstrata2,PACKAGE="mets")$res
return(res)
}# }}}


