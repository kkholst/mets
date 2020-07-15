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
##' @param offset offsets for FG  model
##' @param weights weights for FG score equations
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
##' summary(ll)
##' plot(ll)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pll <- predict(ll,nd)
##' plot(pll)
##' 
##' ## Fine-Gray model
##' fg=cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,propodds=NULL)
##' summary(fg)
##' plot(fg)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pfg <- predict(fg,nd)
##' plot(pfg)
##' 
##' sfg=cifreg(Event(time,cause)~strata(tcell)+platelet+age,data=bmt,cause=1,propodds=NULL)
##' summary(sfg)
##' plot(sfg)
##' @aliases vecAllStrata diffstrata
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
   list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,
	cluster.pos=pos.cluster,n=nrow(X),nevent=sum(status==cause))
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
   } else { id <- as.integer(seq_along(entry))-1;  nid <- nrow(X); }
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
	GtsAl<- Gts <- apply(Gts,2,function(x) exp(cumsum(log(1-x))))
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
	    iid=Uiid,ncluster=nid,
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

##' @export
strataC <- survival:::strata

##' Augmentation for Fine-Gray model based on stratified NPMLE Cif (Aalen-Johansen) 
##'
##' Computes  the augmentation term for each individual as well as the sum
##' \deqn{
##' A(\beta) = \int H(t,X,\beta) \frac{F_2^*(t,s)}{S^*(t,s)} \frac{1}{G_c(t)} dM_c
##' }
##' with 
##' \deqn{
##' H(t,X,\beta) = \int_t^\infty (X - E(\beta,t) ) G_c(t) d\Lambda_1^*i(t,s)
##' }
##' using a KM for \deqn{G_c(t)} and a working model for cumulative baseline
##' related to \deqn{F_1^*(t,s)} and \deqn{s} is strata, 
##' \deqn{S^*(t,s) = 1 - F_1^*(t,s) - F_2^*(t,s)}, and
##' \deqn{E(\beta^p,t)} is given.
##'
##' After a couple of iterations we end up with a solution of 
##' \deqn{
##' \int (X - E(\beta) ) Y_1(t) w(t) dM_1 + A(\beta)
##' }
##' the augmented FG-score. 
##'
##' Standard errors computed under assumption of correct \deqn{G_c} model.
##'
##' @param formula formula with 'Event', strata model for CIF given by strata, and strataC specifies censoring strata
##' @param data data frame
##' @param offset offsets for cox model
##' @param data data frame
##' @param E from FG-model 
##' @param cause of interest 
##' @param cens.code code of censoring 
##' @param km to use Kaplan-Meier
##' @param case.weights weights for FG score equations (that follow dN_1) 
##' @param weights weights for FG score equations
##' @param offset offsets for FG   model
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' rho1 <- 0.2; rho2 <- 10
##' n <- 400
##' beta=c(0.0,-0.1,-0.5,0.3)
##' dats <- simul.cifs(n,rho1,rho2,beta,rc=0.2)
##' dtable(dats,~status)
##' dsort(dats) <- ~time
##' fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
##' summary(fg)
##' 
##' fgaugS <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fg$E)
##' summary(fgaugS)
##' fgaugS2 <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fgaugS$E)
##' 
##' @aliases strataC  simul.cifs setup.cif
##' @export
FG_AugmentCifstrata <- function(formula,data=data,E=NULL,cause=NULL,cens.code=0,km=TRUE,case.weights=NULL,weights=NULL,offset=NULL,...)
{# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset","strataC")
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


  statusC <- (status==cens.code) 
  statusE <- (status==cause) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 

  ## sorting after time and statusC, but event times unique after order
  Zcall <- cbind(status,strata)
  dd <- .Call("FastCoxPrepStrata",entry,exit,statusC,X,id, 
	     trunc,strataC,weights,offset,Zcall,case.weights,PACKAGE="mets")

  Z <- dd$X
  jumps <- dd$jumps+1
  xxstrataC <- c(dd$strata)
  xxstatus  <- dd$Z[,1]
  xxstrata  <- dd$Z[,2]
  other <- which((!(xxstatus %in% c(cens.code,cause)) ) )               
  jumps1 <- which(xxstatus==cause)
  jumpsD <- which(xxstatus!=cens.code)
  rr <- c(dd$sign*exp(dd$offset))
  ## S0 after strata
  S0 = c(revcumsumstrata(rr,strata,nstrata))
  ## S0 after strataC
  S00C = c(revcumsumstrata(rr,xxstrataC,nstrataC))

  ## censoring MG, strataC
  stratJumps <- dd$strata[jumps]
  S00i <- rep(0,length(dd$strata))
  S00i[jumps] <-  1/S00C[jumps]

  ## cif calculation, uses strata {{{
  S0Di <- S02i <- S01i <- rep(0,length(dd$strata))
  S01i[jumps1] <-  1/S0[jumps1]
  S02i[other] <-  1/S0[other]
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
 cif1 <- cumsumstrata(Stm*S01i,xxstrata,nstrata)
 cif2 <- cumsumstrata(Stm*S02i,xxstrata,nstrata)
# }}}

 Et <- matrix(0,nrow(Z),ncol(Z))
 Et[jumps1,] <- E

 Lam1fg <- -log(1-cif1)
 laststrata <- tailstrata(xxstrata,nstrata)
 gtstart <- Lam1fg[laststrata]
 dLam1fg <- c(diffstrata(Lam1fg,xxstrata,nstrata))

 ## Gc~strataC, dLam1fg~strata
 tailcstrata <- tailstrata(xxstrataC,nstrataC)
 Gcstart <- Gc[tailcstrata]

 dstrata <- mystrata(data.frame(cbind(xxstrataC,xxstrata)))
 ndstrata <- attr(dstrata,"nlevel")
 lastt <- tailstrata(dstrata-1,ndstrata)

 ### ## \int_t^\infty G_c^j(t) d\Lambda_1^k(t)
 G0start <- rep(1,nstrataC)
 cLam1fg  <- cumsum2strata(Gc,dLam1fg,xxstrataC,nstrataC,xxstrata,nstrata,G0start)$res
 lastt <- tailstrata(dstrata-1,ndstrata)
 RLam1fg <- cLam1fg[lastt][dstrata]-cLam1fg

 ## E(s) from FG without strata
 ## \int_0^t  G_c^j(s) E(s) d\Lambda_1^k(s)
 fff <- function(x) {
   cx  <- cumsum2strata(Gc,x*dLam1fg,xxstrataC,nstrataC,xxstrata,nstrata,G0start)$res
    return(cx[lastt][dstrata]-cx)
 }
 ERLam1fg0  <- apply(Et,2,fff)

 gt <-  RLam1fg*c(cif2/(Gc*St))
 gt[gt==Inf] <- 0
 gt[is.na(gt)] <- 0
 ERLam1fg  <- ERLam1fg0*c(cif2/(Gc*St))
 ERLam1fg[ERLam1fg==Inf] <- 0
 ERLam1fg[is.na(ERLam1fg)] <- 0

 sss <- headstrata(dstrata-1,ndstrata)
 gtstart <- gt[sss]
 E1dLam0 <- cumsum2strata(gt,S00i,dstrata-1,ndstrata,xxstrataC,nstrataC,gtstart)$res

 fff <- function(x) {
    gtstart <- x[sss]
    cx  <- cumsum2strata(x,S00i,dstrata-1,ndstrata,xxstrataC,nstrataC,gtstart)$res
    return(cx)
 }
 E2dLam0 <- apply(ERLam1fg,2,fff)

 U1 <- matrix(0,nrow(Z),1)
 U2 <- matrix(0,nrow(Z),ncol(Z))
 U1[jumps,] <- gt[jumps]
 U2[jumps,] <- ERLam1fg[jumps,]

 ### Martingale  as a function of time and for all subjects to handle strata 
 MG1t <- Z*c(U1[,,drop=FALSE]-E1dLam0)*rr*c(dd$weights)
 MG2t <- (U2[,,drop=FALSE]-E2dLam0)*rr*c(dd$weights)
 MGt <- MG1t-MG2t
 MGiid <- apply(MGt,2,sumstrata,dd$id,max(dd$id)+1)
 augment <- apply(MGt,2,sum)

 augment <- list(MGiid=MGiid,augment=augment,id=id,id.orig=id.orig,
###  dd=dd,div=div,Et=Et,E=E, gt=gt, ERLam1fg=ERLam1fg,RLam1fg=RLam1fg,ERLam1fg0=ERLam1fg0,
     jumps1=jumps1,jumps=jumps,other=other,
     nstrata=nstrata,nstrataC=nstrataC,dstrata=dstrata,ndstrata=ndstrata,
     cif1=cif1,cif2=cif2,St=St,Gc=Gc,strata=xxstrata,strataC=xxstrataC,time=dd$time)

 ## drop strata's from formula and run wiht augmention term
 # {{{
 drop.strata <- function(x) {
   mm <- unlist(Specials(x,"strata"))
 print(mm)
   for (i in mm) x <- update(x, paste(".~.-strata(",i,")"))
 print(x)
   mm <- unlist(Specials(x,"strataC"))
 print(mm)
   for (i in mm) x <- update(x, paste(".~.-strataC(",i,")"))
 print(x)
   return(x)
 }

### formula=Event(time,status)~Z1+Z2+strata(Z1,Z2)+strataC(Z1)
 print(formula)
 formulans <- drop.strata(formula)
 print(formulans)
# }}}

  if (nstrataC==1) cens.model <- ~+1 else cens.model <- ~strata(strataCC)
  data$strataCC <- cens.strata

 fga <- cifreg(formulans,data=data,cause=cause,
     propodds=NULL,augmentation=augment$augment,cens.model=cens.model,...)

 ## adjust SE and var based on augmentation term
 fga$var.orig <- fga$var
 fga$augment <- augment$augment
 fga$iid <- fga$iid + MGiid %*% fga$ihessian 
 fga$var <- crossprod(fga$iid)
 fga$se.coef <-  diag(fga$var)^.5
 fga$MGciid <- MGiid

 return(fga)
}# }}})

##' @export
simul.cifs <- function(n,rho1,rho2,beta,rc=0.5,depcens=0,rcZ=0.5,bin=1,type=
		       c("cloglog","logistic"),rate=1) {# {{{
p=length(beta)/2
tt <- seq(0,6,by=0.1)
if (length(rate)==1) rate <- rep(rate,2)
Lam1 <- rho1*(1-exp(-tt/rate[1]))
Lam2 <- rho2*(1-exp(-tt/rate[2]))

if (length(bin)==1) bin <- rep(bin,2)

Z=cbind((bin[1]==1)*(2*rbinom(n,1,1/2)-1)+(bin[1]==0)*rnorm(n),(bin[2]==1)*(rbinom(n,1,1/2))+(bin[2]==0)*rnorm(n))
colnames(Z) <- paste("Z",1:2,sep="")

cif1 <- setup.cif(cbind(tt,Lam1),beta[1:2],Znames=colnames(Z),type=type[1])
cif2 <- setup.cif(cbind(tt,Lam2),beta[3:4],Znames=colnames(Z),type=type[1])
###
data <- sim.cifsRestrict(list(cif1,cif2),n,Z=Z)

if (depcens==0) censor=pmin(rexp(n,1)*(1/rc),6) else censor=pmin(rexp(n,1)*(1/(rc*exp(rcZ*Z[,1]))),6)

status=data$status*(data$time<=censor)
time=pmin(data$time,censor)
data <- data.frame(time=time,status=status)
return(cbind(data,Z))

}# }}}

simul.mod <- function(n,rho1,rho2,beta,rc=0.5,k=1,depcens=0) {# {{{
p=length(beta)/2
tt <- seq(0,6,by=0.1)
Lam1 <- rho1*(1-exp(-tt))
Lam2 <- rho2*(1-exp(-tt))

Z=cbind(2*rbinom(n,1,1/2)-1,rnorm(n))
colnames(Z) <- paste("Z",1:2,sep="")
cif1 <- setup.cif(cbind(tt,Lam1),beta[1:2],Znames=colnames(Z),type="cloglog")
cif2 <- setup.cif(cbind(tt,Lam2),beta[3:4],Znames=colnames(Z),type="cloglog")
###
data <- sim.cifsRestrict(list(cif1,cif2),n,Z=Z)

censhaz  <-  cbind(tt,k*tt)
if (depcens==1) {
datc <- rchaz(censhaz,exp(Z[,1]*rc))
} else datc <- rchaz(censhaz,n=n)

data$time <- pmin(data$time,datc$time)
data$status <- ifelse(data$time<datc$time,data$status,0)
dsort(data) <- ~time

return(data)
}# }}}

#' @export 
setup.cif  <- function(cumhazard,coef,Znames=NULL,type="logistic")
{# {{{
cif <- list()
cif$cumhaz <- cumhazard
cif$coef <- coef
cif$model <- type
class(cif) <- "defined"
attr(cif,"znames") <- Znames
return(cif)
}# }}}

