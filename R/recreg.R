##' Recurrent events regression with terminal event 
##'
##' Fits Ghosh-Lin IPCW Cox-type model
##'
##' For Cox type model :
##' \deqn{
##' E(dN_1(t)|X) = \mu_0(t)dt exp(X^T \beta)
##' }
##' by solving Cox-type IPCW weighted score equations 
##' \deqn{
##'  \int (Z - E(t)) w(t) dN_1(t) 
##' }
##' where \deqn{w(t) = G(t) (I(T_i \wedge t < C_i)/G_c(T_i \wedge t))} and
##' \deqn{E(t) = S_1(t)/S_0(t)} and \deqn{S_j(t) = \sum X_i^j w_i(t) \exp(X_i^T \beta)}.
##'
##'
##' The iid decomposition of the beta's are on the form
##' \deqn{
##' \int (Z - E ) w(t) dM_1 + \int q(s)/p(s) dM_c
##' }
##' and returned as iid.
##'
##' Events, deaths and censorings are specified via stop start structure and the Event call, that via a status vector 
##' and cause (code), censoring-codes (cens.code) and death-codes (death.code) indentifies these. See example and vignette. 
##'
##' @param formula formula with 'Event' outcome
##' @param data data frame
##' @param cause of interest (1 default)
##' @param death.code codes for death (terminating event, 2 default)
##' @param cens.code code of censoring (0 default)
##' @param cens.model for stratified Cox model without covariates
##' @param weights weights for score equations
##' @param offset offsets for model
##' @param Gc censoring weights for time argument, default is to calculate these with a Kaplan-Meier estimator, should then give G_c(T_i-)
##' @param wcomp weights for composite outcome, so when cause=c(1,3), we might have wcomp=c(1,2).
##' @param augmentation.type of augmentation when augmentation model is given 
##' @param marks  a mark value can be specified, this is vector from the data-frame where the mark value can be found at all events
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ## data with no ties
##' library(mets)
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##' hf$x <- as.numeric(hf$treatment) 
##' dd <- data.frame(treatment=levels(hf$treatment),id=1)
##'
##' gl <- recreg(Event(entry,time,status)~treatment+cluster(id),data=hf,cause=1,death.code=2)
##' summary(gl)
##' head(iid(gl))
##' pgl <- predict(gl,dd,se=1); plot(pgl,se=1)
##' 
##' ## censoring stratified after treatment 
##' gls <- recreg(Event(entry,time,status)~treatment+cluster(id),data=hf,
##' cause=1,death.code=2,cens.model=~strata(treatment))
##' summary(gls)
##' 
##' glss <- recreg(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,
##' cause=1,death.code=2,cens.model=~strata(treatment))
##' summary(glss)
##' plot(glss)
##' 
##' ## IPCW at 2 years 
##' ll2 <- recregIPCW(Event(entry,time,status)~treatment+cluster(id),data=hf,
##' cause=1,death.code=2,time=2,cens.model=~strata(treatment))
##' summary(ll2)
##' 
##' ll2i <- recregIPCW(Event(entry,time,status)~-1+treatment+cluster(id),data=hf,
##' cause=1,death.code=2,time=2,cens.model=~strata(treatment))
##' summary(ll2i)
##' @aliases marks strataAugment scalecumhaz GLprediid recregIPCW IIDrecreg predicttime
##' @export
recreg <- function(formula,data,cause=1,death.code=2,cens.code=0,cens.model=~1,weights=NULL,offset=NULL,Gc=NULL,wcomp=NULL,marks=NULL,augmentation.type=c("lindyn.augment","lin.augment"),...)
{ ## {{{
	cl <- match.call()
	outi  <- recregBN(formula,data,cause=cause,death.code=death.code,cens.code=cens.code,cens.model=cens.model,weights=weights,offset=offset,Gc=Gc,wcomp=wcomp,marks=marks,...)

	if (!is.null(outi$lindyn.augment)) {
		outA  <- recregBN(formula,data,cause=cause,death.code=death.code,cens.code=cens.code,cens.model=cens.model,weights=weights,offset=offset,Gc=Gc,wcomp=wcomp,marks=marks,augmentation=outi[[augmentation.type[1]]],...)
		outi <- outA
	}

	outi$call <- cl
	return(outi)
} ## }}}

##' @export
marks <- function(x) x

recregBN <- function(formula,data=data,cause=c(1),death.code=c(2),cens.code=c(0),cens.model=~1,weights=NULL,offset=NULL,Gc=NULL,wcomp=NULL,marks=NULL,...)
{ ## {{{
cl <- match.call()
    m <- match.call(expand.dots = TRUE)[1:3]
    des <- proc_design(
        formula,
        data = data,
        specials = c("offset", "weights", "cluster","strata","marks"),
        intercept = FALSE
    )
    Y <- des$y
    if (!inherits(Y, c("Event"))) {
        stop("Expected an 'Event'-object")
    }
    if (ncol(Y) == 2) {
        exit <- Y[, 1]
        entry <- rep(0, nrow(Y))
        status <- Y[, 2]
    } else {
        entry <- Y[, 1]
        exit <- Y[, 2]
        status <- Y[, 3]
    }
    X <- des$x
    strata <- des$strata
    if (!is.null(strata))  {
      ns <- grep("strata",names(des$levels))
      strata.name  <-  names(des$levels)[1]
    } else strata.name <- NULL
    des.weights <- des$weights
    des.offset  <- des$offset
    des.marks <- des$marks
    id      <- des$cluster
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    ## no use of 
    pos.cluster <- pos.strata <- NULL

 ## take offset and weight first from formula, but then from arguments
  if (is.null(des.offset)) {
	  if (is.null(offset)) offset <- rep(0,length(exit)) 
  } else offset <- des.offset
  if (is.null(des.weights)) {
	  if (is.null(weights)) weights <- rep(1,length(exit)) 
  } else weights <- des.weights
  if (!is.null(des.marks) & is.null(marks))  marks <- des.marks

 res <- c(recregN01(data,X,entry,exit,status,id=id,strata=strata,offset=offset,weights=weights,strata.name=strata.name,cens.model=cens.model,cause=cause,
		   death.code=death.code,cens.code=cens.code,Gc=Gc,wcomp=wcomp,
		   case.weights=marks,...),
 list(call=cl,formula=formula,strata.pos=pos.strata,
      cluster.pos=pos.cluster,n=length(status),nevent=sum(status %in% cause)))
	colnames(res$iid) <- names(res$coef)
  res$design <- des

	class(res) <- c("recreg","phreg")
	return(res)
} ## }}}

recregN01 <- function(data,X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
		      strata.name=NULL,beta,stderr=1,method="NR",no.opt=FALSE,propodds=NULL,zero.remove=1,
		      case.weights=NULL,cause=1,death.code=2,cens.code=0,Gc=NULL,cens.model=~+1,augmentation=NULL,
		      cox.prep=TRUE,wcomp=NULL,augment.model=NULL,adm.cens.time=NULL,ftime.augment=NULL,twostage=FALSE,...) { ## {{{
	p <- ncol(X)  # setting up weights, strata, beta and so forth before the action starts
	if (missing(beta)) beta <- rep(0,p)
	if (p==0) X <- cbind(rep(0,length(exit)))
	augmentation.call <- augmentation
	if (is.null(augmentation)) augmentation <- 0
	cause.jumps <- which(status %in% cause)
	if (length(cause.jumps)>0) {
		max.jump <- max(exit[cause.jumps])
		other <- which((status %in% death.code ) & (exit< max.jump))
	} else  {
		###	warning("no jumps of cause type\n"); 
		max.jump <- max(exit)
		other <- which((status %in% death.code ) )
	}

	n <- length(exit)
	if (is.null(strata)) {
		strata <- rep(0,length(exit))
		nstrata <- 1
		strata.level <- NULL
	} else {
		strata.level <- levels(strata)
		ustrata <- sort(unique(strata))
		nstrata <- length(ustrata)
		strata.values <- ustrata
		if (is.numeric(strata))
			strata <-  fast.approx(ustrata,strata)-1
		else  {
			strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
		}
	}

	if (is.null(entry)) entry <- rep(0,length(exit))
	trunc <- (any(entry>0))
	if (is.null(offset)) offset <- rep(0,length(exit))
	if (is.null(weights)) weights <- rep(1,length(exit))
	if (is.null(case.weights)) case.weights <- rep(1,length(exit))
	if (!is.null(wcomp))  {
		if (length(wcomp)!=length(cause)) stop("weights follow the causes, and length must be the same\n"); 
		wwcomp <- rep(1,length(exit)); 
		k <- 1
		for (i in cause) { wwcomp[status==i] <- wcomp[k];k <- k+1} 
		case.weights <- case.weights* wwcomp
	}
	strata.call <- strata

	call.id <- id 
	conid <- construct_id(id,nrow(X))
	name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
	orig.id <- id

	### censoring weights constructed
	whereC <- which(status %in% cens.code)
	time <- exit
	statusC <- c(status %in% cens.code)
	data$id__ <- id  
	data$exit__ <- exit
	data$entry__ <- entry
	data$statusC <- statusC
	data$status__ <- (status %in% cause)*1
	cens.strata <- cens.nstrata <- NULL
	## lag-count to use for augment.model=~Nt+X1+X2
	data <- count.history(data,status="status__",id="id__",types=cause,multitype=TRUE)
	data$Nt <- data[,paste("Count",cause[1],sep="")]

	## augmentation model and remove intercept
	if (!is.null(augment.model)) { XXA <- model.matrix(augment.model,data)[,-1,drop=FALSE]; namesXXA <- colnames(XXA); } else XXA <- NULL

	if ((length(whereC)>0) ) { ## {{{
		if (is.null(Gc)) {
			kmt <- TRUE
			if (class(cens.model)[1]=="formula") {
				formC <- update.formula(cens.model,Surv(entry__,exit__,statusC)~ . +cluster(id__))
				cens.model <- phreg(formC,data)
			}
			if (cens.model$p>0) kmt <- FALSE
			###        Pcens.model <- predict(cens.model,data,times=exit,tminus=TRUE,individual.time=TRUE,se=FALSE,km=kmt)
			Pcens.model <- predict(cens.model,data,times=exit,individual.time=TRUE,se=FALSE,km=kmt)
			Stime <- Pcens.model$surv <- c(Pcens.model$surv)
			## strata from original data
			nCstrata <- cens.model$nstrata
			cens.strata <- Pcens.model$strata
		} else {
			formC <- NULL
			Stime <- Gc
			Pcens.model <- list(time=exit,surv=Gc,strata=0)
			nCstrata <- 1
			cens.strata <- rep(0,length(exit))
		}
	} else { 
		formC <- NULL
		Stime <- Gc  <- rep(1,length(exit))
		Pcens.model <- list(time=exit,surv=Gc,strata=0)
		nCstrata <- 1
		cens.strata <- rep(0,length(exit))
	} ## }}}

	Zcall <- cbind(status,cens.strata,Stime,strata,strata,1) ## to keep track of status and Censoring strata

	trunc <- TRUE
	if (twostage) { ## to use old structure for twostageREC
		xxg <- .Call("FastCoxPrepStrata",entry,exit,(status %in% cause)*1,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")
		xxg$nstrata <- nstrata
	}

	## computing terms for those experiencing another cause, need S0, S1, S2
	if ((length(other)>=1)) { ## {{{ 
		weightso <- weights[other]/Stime[other]
		timeoo <- rep(max(exit)+1,length(other))
		if (is.null(adm.cens.time))
			timeoo <- rep(max(exit)+1,length(other)) else timeoo <- adm.cens.time[other] 
		statuso <- rep(cens.code[1],length(other))
		Xo <- X[other,,drop=FALSE]
		offseto <- offset[other]
		entryo <- exit[other]
		ido <- id[other]
		stratao <- strata[other]
		type <- rep(2,length(other))
		###

		entry <- c(entry,entryo)
		exit <- c(exit,timeoo)
		status <- c(status,statuso)
		X <- rbind(X,Xo)
		id <- c(id,ido)
		strata <- c(strata,stratao)
		weights <- c(weights,weightso)
		offset <- c(offset,offseto)
		case.weights <- c(case.weights,case.weights[other])

		Zcallo <-  Zcall[other,]
		Zcallo[,6] <- 2
		Zcall <- rbind(Zcall,Zcallo)
		Zcall <- cbind(Zcall,rbind(XXA,XXA[other,,drop=FALSE]))
	} ## }}}

	stat1 <- 1*(status %in% cause)
	xx2 <- .Call("FastCoxPrepStrata",entry,exit,stat1,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")

	### remove all the initial zero's to save space/time, that is, the first streak of sign=-1, 
	### that is not important for risk and events
	if (zero.remove==1) { ## {{{
		first <- which(xx2$sign==1)[1]-1
		if (first >1) {
			fentry <- (1:first)
			xx2$id     <- xx2$id[-fentry]
			xx2$time   <- xx2$time[-fentry]
			xx2$status <- xx2$status[-fentry]
			xx2$sign   <- xx2$sign[-fentry]
			xx2$X      <- xx2$X[-fentry,,drop=FALSE]
			xx2$XX     <- xx2$XX[-fentry,,drop=FALSE]
			if (nrow(xx2$ZX)==nrow(xx2$X)) xx$ZX <- xx2$ZX[-fentry,,drop=FALSE]
			xx2$Z      <- xx2$Z[-fentry,,drop=FALSE]
			xx2$offset <-xx2$offset[-fentry]
			xx2$weights <-xx2$weights[-fentry]
			xx2$caseweights <-xx2$caseweights[-fentry]
			xx2$strata <-xx2$strata[-fentry]
			xx2$jumps <- xx2$jumps-first
		}
	} ## }}}

	jumps <- xx2$jumps+1
	typexx2 <- xx2$Z[,6]
	Xj <- xx2$X[jumps,,drop=FALSE]
	xx2$nstrata <- nstrata
	jumptimes <- xx2$time[jumps]
	strata1jumptimes <- xx2$strata[jumps]
	if ((length(whereC)>0)) {
		###
		rr0 <- xx2$sign*(typexx2==1)
		jumpsC <- which((xx2$Z[,1] %in% cens.code) & xx2$sign==1 & typexx2==1)
		strataCxx2 <- xx2$Z[,2]
		S0iC2  <-  S0iC <- rep(0,length(xx2$status))
		S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
		S0iC[jumpsC] <- 1/S0rrr[jumpsC]
		S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
		## Gc(t) computed  along all times of combined data-set: data + [D,\infty] 
		Gcxx2 <- exp(cumsumstrata(log(1-S0iC),strataCxx2,nCstrata))
		Gstart <- rep(1,nCstrata)
		Gjumps <- Gcxx2[jumps,]

		njumps <- length(jumps)
		GtsS0ooN <-   .Call("_mets_GcjumpsR",Gcxx2,c(xx2$status),strataCxx2,nCstrata,Gstart,njumps)
		Gts <- GtsS0ooN$Gcjumps
	} else {
		Gcxx2 <- rep(1,length(xx2$sign))
		strataCxx2 <- rep(0,length(xx2$sign))
		nCstrata <- 1
		Gstart <- 1
		Gjumps <- 1
		Gts <- matrix(1,length(jumps),1)
	} 

	obj <- function(pp,all=FALSE) { ## {{{

		rr2 <- c(xx2$sign*exp(xx2$X %*% pp + xx2$offset)*xx2$weights)
		rr2now <- c(xx2$sign*exp(xx2$X %*% pp + xx2$offset))

		###     S0ooN <-   .Call("_mets_S0_FG_GcR",rr2,Gcxx2,typexx2-1,c(xx2$status),xx2$strata,xx2$nstrata,strataCxx2,nCstrata,Gstart)

		###     S0ooA <-   .Call("_mets_S0_FGRN",rr2,typexx2-1,c(xx2$status),xx2$strata,xx2$nstrata,strataCxx2,nCstrata,Gts)
		S0oo <-   .Call("_mets_S0_FGRN",rr2,typexx2-1,c(xx2$status),xx2$strata,xx2$nstrata,strataCxx2,nCstrata,Gts)$S0

		f  <-  function(x) {
			ll <-   .Call("_mets_S0_FGRN",x,typexx2-1,c(xx2$status),xx2$strata,xx2$nstrata,strataCxx2,nCstrata,Gts)$S0
		}
		###	f  <-  function(x) {
		###           ll <-   .Call("_mets_S0_FG_GcR",x,Gcxx2,typexx2-1,c(xx2$status),xx2$strata,xx2$nstrata,strataCxx2,nCstrata,Gstart)$S0
		###	}

		S1oo  <- apply(xx2$X*rr2,2,f)
		S2oo  <- apply(xx2$XX*rr2,2,f)

		S0 <- c(S0oo) ## [jumps,]
		S1 <- S1oo ## [jumps,,drop=FALSE]
		S2 <- S2oo ## [jumps,,drop=FALSE]

		E <- S1/S0
		weightsJ <- xx2$weights[jumps]
		caseweightsJ <- xx2$caseweights[jumps]
		strataJ <- xx2$strata[jumps]
		rr2now <- rr2now[jumps]
		U <- (Xj-E)
		ploglik <- (log(rr2now)-log(S0))*weightsJ*caseweightsJ;

		if (!is.null(propodds)) {
			pow <- c(.Call("cumsumstrataPOR",weightsJ,S0,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$pow);
			DLam <-.Call("DLambetaR",weightsJ,S0,E,Xj,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$res;
			Dwbeta <- DLam*rr2now+(pow-1)*Xj
			DUadj  <- .Call("vecMatMat",Dwbeta,U,PACKAGE="mets")$vXZ
		}

		Ut <- caseweightsJ*weightsJ*U
		## E^2, as n x (pxp)
		Et2 <-  .Call("vecCPMat",E,PACKAGE="mets")$XX
		S2S0 <-  S2/S0
		DUt <-  -(S2S0-Et2)
		np <- length(pp)

		if (!is.null(propodds)) {
			Ut  <- pow*Ut
			S0 <- S0/pow
			DUt <- pow*DUt
			DUt <- .Call("XXMatFULL",DUt,np,PACKAGE="mets")$XXf
			if (ncol(DUt)>0) DUt <- DUt+DUadj 
			ploglik <- pow*ploglik
		}

		U  <- apply(Ut,2,sum)
		DUt <- caseweightsJ*weightsJ*DUt
		DU <- -apply(DUt,2,sum)
		np <- length(pp)
		if (ncol(DUt)!=p*p) {
			DU <- matrix(.Call("XXMatFULL",matrix(DU,nrow=1),np,PACKAGE="mets")$XXf,np,np)
		} else  DU <- matrix(DU,p,p)
		ploglik <- sum(ploglik)
		U <- U+augmentation

		out <- list(ploglik=ploglik,gradient=U,hessian=-DU,cox.prep=xx2,
			    hessiantime=DUt,weightsJ=weightsJ,caseweightsJ=caseweightsJ,
			    jumptimes=jumptimes,strata=strataJ,nstrata=nstrata,S0s=S0,
			    time=jumptimes,S0=S0/(caseweightsJ*weightsJ),S2S0=S2S0,E=E,U=Ut,X=Xj,Gjumps=Gjumps)

		nn <- length(xx2$time)
		if (all)
			return(out)
		else
			with(out,structure(-ploglik/nn, gradient=-gradient/nn, hessian=-hessian/nn))
	} ## }}}

	if (length(jumps)==0) no.opt <- TRUE
	opt <- NULL
	if (p>0) {
		if (no.opt==FALSE) {
			if (tolower(method)=="nr") {
				opt <- lava::NR(beta,obj,...)
				opt$estimate <- opt$par
			} else {
				opt <- nlm(obj,beta,...)
				opt$method <- "nlm"
			}
			cc <- opt$estimate;  names(cc) <- colnames(X)
			if (stderr==2) return(cc)
			val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
		} else val <- c(list(coef=beta),obj(beta,all=TRUE))
	} else {
		no.opt <- TRUE
		val <- obj(0,all=TRUE)
	}

	beta.s <- val$coef
	if (is.null(beta.s)) beta.s <- 0
	## getting final S's
	opt <-  val ## obj(beta.s,all=TRUE)

	if (p>0) {
		iH <- - tryCatch(solve(opt$hessian),error=function(e) matrix(0,nrow(opt$hessian),ncol(opt$hessian)) )
		opt$ihessian <- iH
		opt$no.opt <- FALSE
		dd <- IIDrecreg(xx2,opt,cause=cause,cens.code=cens.code,death.code=death.code,adm.cens=adm.cens.time) 

		Uiid <-  dd$beta.iid.naive 
		UUiid <- dd$beta.iid
		UU <- dd$MGt
		MGc <- dd$MGc

		Uiid <-  (UU+MGc) %*% iH
		UUiid <- UU %*% iH
		var1 <-  crossprod(UUiid)
		varmc <-  crossprod(Uiid)

		## compute regression augmentation for censoring martingale 
		if ((!is.null(augment.model)) & (length(whereC)>0) ) {# ## {{{ 

			CovZXstrata <- function(X,Ej,Z,Sign,strata,nstrata,jumps) 
			{
				strata  <- c(strata); Sign <- c(Sign)
				###	Ej <- Ej[jumps,,drop=FALSE]; Ej <- Ej
				ZE <- apply(Z*Sign,2,revcumsumstrata,strata,nstrata)[jumps,,drop=FALSE]; 
				XZ  <- .Call("vecMatMat",X,Z)$vXZ;  
				XZ <- apply(XZ*Sign,2,revcumsumstrata,strata,nstrata)[jumps,,drop=FALSE]; 
				EXZ  <- .Call("vecMatMat",Ej,ZE)$vXZ;  
				out <- XZ-EXZ
				return(out)
			}

			## regress U(s)=\int_s^\infty (Z-E) w(s) dM(s) on agument-model among survivors 
			## U(s) = U(\infty) - \int_0^s (Z-E) w(s)  dM(s)
			## sum (e_i - \bar e) U(s) Y_i(s)

			# construct censoring weights going along with all data, with added [D,\infty], start stop
			rr0 <- c(xx2$sign)*(typexx2==1)
			jumpsC <- which((xx2$Z[,1] %in% cens.code) & xx2$sign==1 & typexx2==1)
			strataCxx2 <- xx2$Z[,2]
			S0iC2  <-  S0iC <- rep(0,length(xx2$status))
			S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
			S0iC[jumpsC] <- 1/S0rrr[jumpsC]
			S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
			S0c <- c(S0rrr[jumpsC])
			## Gc(t) computed  as exp(- Cumhaz) to avoid some "0"s
			Gcj <- Gcxx2 <- exp(-cumsumstrata(S0iC,strataCxx2,nCstrata))[jumpsC]
			#

			XXA <- xx2$Z[,-(1:6),drop=FALSE]
			EXXA <- apply(XXA*c(rr0),2,revcumsumstrata,strataCxx2,nCstrata)
			EA <- EXXA[jumpsC,,drop=FALSE]/S0rrr[jumpsC]
			UA <- (XXA[jumpsC,,drop=FALSE]-EA)

			###
			E2A <- .Call("vecMatMat",EA,EA)$vXZ;  
			XX2A <- .Call("vecMatMat",XXA,XXA)$vXZ;  
			S2A <- apply(XX2A*c(rr0),2,revcumsumstrata,strataCxx2,nCstrata)
			###
			hessiant <- -(S2A[jumpsC,,drop=FALSE]/S0c-E2A)
			hesst <- hessiant

			### X fra GL + tail-death 
			rr <- c(exp(xx2$X %*% beta.s+ xx2$offset)*xx2$weights)*(typexx2==1)
			Zrr <- xx2$X*rr
			ZEdN <- apply(dd$Ut,2,revcumsumstrata,xx2$id,nid)

			covXsZ <-   CovZXstrata(XXA,EA,Zrr,rr0,strataCxx2,nCstrata,jumpsC) 
			covXsrr <-  CovZXstrata(XXA,EA,as.matrix(rr,ncol=1),rr0,strataCxx2,nCstrata,jumpsC) 
			covXsUs3 <- .Call("vecMatMat",covXsrr,dd$EdLam0[jumpsC,,drop=FALSE])$vXZ;  
			covXsUs2 <- covXsZ*dd$cumhaz[jumpsC]-covXsUs3 
			### U(infty)= UU
			Uinfiid <- UU[xx2$id+1,,drop=FALSE]
			fid <- headstrata(xx2$id,nid)
			cZEdN <- ZEdN[fid,,drop=FALSE][xx2$id+1,,drop=FALSE]-ZEdN
			Us1 <- Uinfiid-cZEdN
			covXsUs1 <- CovZXstrata(XXA,EA,Us1,rr0,strataCxx2,nCstrata,jumpsC) 
			## scale with Y_(s) because hessiantime is also scaled with this 
			covXsYs <- (covXsUs1+covXsUs2)/S0c; ## /c(cr2$S0)

			pXXA <- ncol(XXA)
			gammat <-  .Call("CubeMattime",hesst,covXsYs,pXXA,pXXA,pXXA,p,1,0,0,PACKAGE="mets")$XXX
			gammat[is.na(gammat)] <- 0
			gammat[gammat==Inf] <- 0
			namesG <- c(); for (i in 1:p) namesG <- c(namesG,paste(namesXXA,"-e",i,sep=""))
			colnames(gammat) <- namesG
			augmentt <- .Call("CubeMattime",gammat,UA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
			augment.times <- -apply(augmentt,2,sum)
			gain.times <- .Call("CubeMattime",covXsYs,gammat,pXXA,p,pXXA,p,0,1,0,PACKAGE="mets")$XXX
			gain.times <- matrix(apply(gain.times,2,sum),p,p)
			var.augment.times <-  gain.times 

			###
			time.gammat <- timeC <- xx2$time[jumpsC]
			if (is.null(ftime.augment)) {
				### simple default parabola
				maxt <- max(timeC)
				ftime <- timeC*(timeC-maxt)/maxt^2
			} else { 
				if (is.list(ftime.augment)) ftime <- ftime.augment[[1]](timeC) else ftime <- ftime.augment(timeC)
				if (length(ftime.augment)==2) {
					timepar <- ftime.augment[[2]](timeC)
					parap <- lm(gammat~-1+timepar)
					gammat <- parap$fitted.values
				}
			}
			ftime.gamma <- ftime
			varZdN <- matrix(apply(ftime^2*hesst/c(Gcj^2),2,sum),pXXA,pXXA)
			covXYdN <- matrix(apply(ftime*covXsYs/c(Gcj),2,sum),p,pXXA,byrow=TRUE) 
			gamma <- -1*.Call("CubeMattime",matrix(varZdN,nrow=1),matrix(covXYdN,nrow=1),pXXA,pXXA,p,pXXA,1,0,1,PACKAGE="mets")$XXX
			gamma <- matrix(gamma,p,pXXA,byrow=TRUE)
			gamma[is.na(gamma)] <- 0; gamma[gamma==Inf] <- 0
			colnames(gamma) <- namesXXA
			augment <- c(gamma %*% apply(ftime*UA/c(Gcj),2,sum))
			var.augment  <-  gamma %*% t(covXYdN) ###  /(nid^2)

			if (!is.null(augmentation.call)) { ## update variance when called with augmentation
				#### iid magic  for censoring augmentation martingale
				### int_0^infty gamma (e_i - ebar(s)) 1/G_c(s) dM_i^c
				S0iG <- S0i <- rep(0,length(xx2$strata))
				S0iG[jumpsC] <- ftime/(S0rrr[jumpsC]*c(Gcj))
				S0i[jumpsC] <- c(1/S0rrr[jumpsC])
				U <- E <- matrix(0,nrow(xx2$X),pXXA)
				E[jumpsC,] <- EA; 
				U[jumpsC,] <- ftime*UA/c(Gcj)
				cumhaz <- cumsumstrata(S0iG,strataCxx2,nCstrata)
				EdLam0 <- apply(E*S0iG,2,cumsumstrata,strataCxx2,nCstrata)
				MGCt <- U[,drop=FALSE]-(XXA*c(cumhaz)-EdLam0)*c(rr0)
				MGCtiid <- apply(MGCt,2,sumstrata,xx2$id,nid)
				iid.augment <-  (MGCtiid %*% t(gamma)) %*% iH

				gammasEs <- .Call("CubeMattime",gammat,EA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
				gammasE <- matrix(0,nrow(XXA),p)
				gammatt  <-    matrix(0,nrow(XXA),pXXA*p)
				gammasE[jumpsC,] <- gammasEs
				gammatt[jumpsC,] <- gammat
				gammaEsdLam0 <- apply(gammasE*S0i,2,cumsumstrata,strataCxx2,nCstrata)
				gammadLam0 <-   apply(gammatt*S0i,2,cumsumstrata,strataCxx2,nCstrata)
				XgammadLam0 <- .Call("CubeMattime",gammadLam0,XXA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
				Ut <- Et <- matrix(0,nrow(XXA),p)
				Ut[jumpsC,] <- augmentt
				MGCtt <- Ut[,drop=FALSE]-(XgammadLam0-gammaEsdLam0)*c(rr0)
				MGCttiid <- apply(MGCtt,2,sumstrata,xx2$id,nid)
				iid.augment.times <-  MGCttiid %*% iH
				Uiid.augment <- Uiid-iid.augment
				Uiid.augment.times <- Uiid-iid.augment.times
				## so that iid corresponds to dynamic censorig augmentation 
				Uiid <- Uiid.augment.times
				## scale with G_c(t) to compare with gamma
				gammat <- gammat * c(Gcj)

				var.augment <-  varmc  -  iH %*% var.augment %*% iH
				var.augment.times <-  varmc  +  iH %*% var.augment.times %*% iH
				varmc <- var.augment.times
				var.augment.iid <-  crossprod(Uiid.augment) 
				var.augment.times.iid <- crossprod(Uiid.augment.times) 
			} else {
				var.augment <-  var.augment.times <-  var.augment.iid <-  var.augment.times.iid <- NULL
				Uiid.augment.times <- Uiid.augment <- NULL
			}
		} else {
			iid.augment <- iid.augment.times <- augment <- augment.times <- NULL 
			var.augment.times <- var.augment <- NULL
			var.augment.times.iid <- var.augment.iid <- NULL
			Uiid.augment.times <- Uiid.augment <- NULL
			time.gammat <- gamma <- gammat <- NULL
			ftime.gamma <- NULL
			Gcj <- NULL
		}  ## }}}
		### end if (p>0)
	} else {
		iid.augment <- iid.augment.times <- augment <- augment.times <- NULL 
		var.augment.times <- var.augment <- NULL
		var.augment.times.iid <- var.augment.iid <- NULL
		Uiid.augment.times <- Uiid.augment <- NULL
		time.gammat <- gamma <- gammat <- NULL
		ftime.gamma <- NULL
		Gcj <- NULL
		varmc <- var1 <- 0; MGc <- iH <- UUiid <- Uiid <- NULL
	}
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

	if (!is.null(call.id)) {
		MGc <-  nameme(MGc,name.id)
		Uiid <- nameme(Uiid,name.id)
		UUiid <- nameme(UUiid,name.id)
		Uiid.augment <- nameme(Uiid.augment,name.id)
		Uiid.augment.times  <- nameme(Uiid.augment.times,name.id)
	}

	out <- list(coef=beta.s,var=varmc,se.coef=diag(varmc)^.5,iid.naive=UUiid,
   iid=Uiid,ncluster=nid,ihessian=iH,hessian=opt$hessian,var1=var1,se1.coef=diag(var1)^.5,
		    hessianttime=opt$hessianttime,
		    ploglik=opt$ploglik,gradient=opt$gradient,
		    cumhaz=cumhaz,se.cumhaz=se.cumhaz,MGciid=MGc,
		    id=orig.id,call.id=call.id,name.id=name.id,nid=nid,
		    strata.jumps=opt$strata,strata=xx2$strata,
		    nstrata=nstrata,strata.name=strata.name,strata.level=strata.level,
		    propodds=propodds,
		    S0=opt$S0,E=opt$E,S2S0=opt$S2S0,time=opt$time,Ut=opt$U,
		    jumps=jumps,exit=exit,
		    p=p,S0s=val$S0s,
		    no.opt=no.opt,##n=nrow(X),nevent=length(jumps),
		    Pcens.model=Pcens.model,Gjumps=Gjumps,
		    cens.code=cens.code,death.code=death.code,cause=cause,
		    augmentation=augmentation.call,
		    var.augment=var.augment,var.augment.times=var.augment.times,
		    var.augment.iid=var.augment.iid,var.augment.times.iid=var.augment.times.iid,
		    lin.augment=c(augment),lindyn.augment=c(augment.times),
		    iid.augment=Uiid.augment,iid.augment.times=Uiid.augment.times,
		    gamma=gamma, gamma.times=gammat, time.gammat=time.gammat,ftime.gamma=ftime.gamma,Gcj=Gcj,
		    adm.cens.time=adm.cens.time,cens.weights=Stime
		    )

	if (cox.prep) out <- c(out,list(cox.prep=xx2))
	if (twostage) out <- c(out,list(cox.prep.twostage=xxg))

	return(out)
} ## }}}

##' @export
IIDrecreg <- function(coxprep,x,time=NULL,cause=1,cens.code=0,death.code=2,fixbeta=NULL,beta.iid=NULL,adm.cens=NULL,tminus=FALSE)
{ ## {{{
	if (is.null(fixbeta)) 
	if ((x$no.opt) | is.null(x$coef)) fixbeta<- 1 else fixbeta <- 0

	xx2 <- coxprep
	status <- xx2$Z[,1]
	cause.jumps <- xx2$jumps+1 
	exit <- xx2$time
	max.jump <- max(exit[cause.jumps])+1
	other <- (length(which(status %in% death.code)) > 0)
	whereC <- which( (status %in% cens.code) & xx2$sign==1)
        anyC <- (length(whereC)>0)

	# construct censoring weights going along with all data, with added [D,\infty], start stop
	jumps <- xx2$jumps+1
	jumptimes <- xx2$time[jumps]
	strata1jumptimes <- xx2$strata[jumps]
	Xj <- xx2$X[jumps,,drop=FALSE]
	###
	typexx2 <- xx2$Z[,6]
	rr0 <- xx2$sign*(typexx2==1)
	jumpsC <- which((xx2$Z[,1] %in% cens.code) & xx2$sign==1 & typexx2==1)

	if (anyC) {
		strataCxx2 <- xx2$Z[,2]
		S0iC2  <-  S0iC <- rep(0,length(xx2$status))
		nCstrata <- max(strataCxx2)+1
		S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
		if (length(jumpsC)>0) {
			S0iC[jumpsC] <- 1/S0rrr[jumpsC]
			S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
		}
		## Gc(t) computed  along all times of combined data-set: data + [D,\infty] 
		Gcxx2 <- exp(cumsumstrata(log(1-S0iC),strataCxx2,nCstrata))
		Gstart <- rep(1,nCstrata)
		Gjumps <- Gcxx2[jumps,]
	} else  Gcxx2 <- rep(1,length(xx2$stata))
	### all administrative censoring, or no censoring at all 
	if (!anyC) typexx2 <- 1

	### iid version given G_c when covariates are there, iid robust 
	S0i <- rep(0,length(xx2$strata))
	S0i[jumps] <- 1/x$S0
	Z <- xx2$X
	p <- ncol(x$E)
	if ( (!is.null(beta.iid)) | fixbeta==0) {
		U <- E <- matrix(0,nrow(Z),p)
		E[jumps,] <- x$E
		U[jumps,] <- x$U
		EdLam0 <- apply(E*S0i,2,cumsumstrata,xx2$strata,xx2$nstrata)
	} else U <- NULL
	cumhazA <- cumsumstratasum(S0i,xx2$strata,xx2$nstrata,type="all")
	cumhaz <- c(cumhazA$sum)

	## make indicator for times of interest for integrals up to time 
	if (!is.null(time))  btimexx <- (xx2$time<time) else btimexx <- rep(1,length(xx2$time))

	if (fixbeta==0) {
		rr <- c(xx2$sign*exp(Z %*% x$coef + xx2$offset))
		Ht <- apply(E*S0i*btimexx,2,cumsumstrata,xx2$strata,xx2$nstrata); 
	} else { Ht <- NULL; rr <- c(xx2$sign*exp(xx2$offset)) }
	rrw <- rr*c(xx2$weights)

	mid <- max(xx2$id)
	if ( (!is.null(beta.iid)) | fixbeta==0) {
		### Martingale  as a function of time and for all subjects to handle strata
		MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rrw*(typexx2==1)
		UU <- apply(MGt,2,sumstrata,xx2$id,mid+1)
	} else UU <- 0

	if (!is.null(time)) {
		## baseline
		MGAiid <- NULL
		S0i2 <- rep(0,length(xx2$strata))
		ww <- xx2$caseweights*xx2$weights
		S0i2[jumps] <- 1/(x$S0^2*ww[jumps])
		MGAiid <- matrix(0,length(S0i2),1)
		MGAiid2 <- matrix(0,length(S0i2),1)
		cumhazAA <- cumsumstrata(S0i2*btimexx,xx2$strata,xx2$nstrata)
		MGAiid <- S0i*btimexx-cumhazAA*rrw*(typexx2==1)
	} else MGAiid <- NULL

	if (other & anyC) { ## "martingale part" for type-2 after T ## {{{
	   ## tail part with \int (Z_i-E) w_i(t) dM_i = \int_D_i^\tau (Z_i-E) Gc(t) dM_i/Gc(D_i) 
	   rrw2 <- rrw*(typexx2==2)
	   GdL <- c(cumsum2strata(Gcxx2,S0i,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)$res)
	   fff <- function(x) {
		cx  <- cumsum2strata(Gcxx2,x*S0i,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)$res
		return(cx)
	   }
           if ((!is.null(beta.iid)) | fixbeta==0) EGdL  <- apply(E,2,fff)

           if ( ((!is.null(beta.iid)) | fixbeta==0)  ) {
		MGt2  <- -(Z*GdL-EGdL)*rrw2
		UU2 <- apply(MGt2,2,sumstrata,xx2$id,mid+1)
		UU  <-  UU+UU2
	   }

	   dstrata <- mystrata(data.frame(strataCxx2,xx2$strata))
	   ndstrata <- attr(dstrata,"nlevel")
	   lastt <- tailstrata(dstrata-1,ndstrata)

	   if (!is.null(time) ) {
	   HBt <- cumsum2strata(Gcxx2,S0i2*btimexx,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)
	   HBtinf <- HBt$res[lastt][dstrata]-HBt$res
###	   if ( is.null(adm.cens)) {
		## baseline
		MGAiid2 <- -HBt$res*c(rrw2)
		MGAiid <- MGAiid+MGAiid2
###		MGAiid <- apply(MGAiid,2,sumstrata,xx2$id,mid+1)
###	   }
	   }
	   if ( ((!is.null(beta.iid)) | fixbeta==0) ) {
		Htinf <- GdL[lastt][dstrata]-GdL
		ff <- function(x) x[lastt][dstrata]-x
		EHtinf  <- apply(EGdL,2,ff)
	   }
	} # ## }}}

	if (!is.null(time) )  MGAiid <- apply(MGAiid,2,sumstrata,xx2$id,mid+1)

	## censoring terms for influene functions 
	if (other & anyC) { ##  ## {{{
		### Censoring adjustment for jumps of other type but only for KM-case 
		### first time we see them with type2 event 

		if (is.null(adm.cens)) {
                   rrw2j <- -c(rrw2*(xx2$sign==-1))
		   Xos <- Z*rrw2j
		   rrsx <- cumsumstrata(rrw2j,strataCxx2,nCstrata)
		   Xos <- apply(Xos,2,cumsumstrata,strataCxx2,nCstrata)

		   if ( (!is.null(beta.iid)) | fixbeta==0) q2 <- (Xos*c(Htinf)-EHtinf*c(rrsx))
		   if (!is.null(time))  qB2 <- rrsx*c(HBtinf) 
		} else {
                   lastid <- tailstrata(xx2$id,mid+1)	
	 	   GdLast <- GdL[lastid][xx2$id+1]
		   EGdLast <- EGdL[lastid,][xx2$id+1,]
                        
		   GadmXE2 <- apply(GdLast*Z*rrw2,2,revcumsumstrata,strataCxx2,nCstrata)
		   XE2 <- apply(Z*rrw2,2,revcumsumstrata,strataCxx2,nCstrata)
		   EGadmRR2 <- apply(EGdLast*rrw2,2,revcumsumstrata,strataCxx2,nCstrata)
		   RR2 <- revcumsumstrata(rrw2,strataCxx2,nCstrata)

		   q11 <- GadmXE2-XE2*GdL
		   q21 <- EGadmRR2-c(RR2)*EGdL
                   ###
		   qq <- q11-q21

                   if ( (!is.null(beta.iid)) | fixbeta==0) q2 <- qq
                   if (!is.null(time))  {
		       HBlast <- HBt$res[lastid][xx2$id+1]
	               HBadmRR2 <- revcumsumstrata(HBlast*rrw2,strataCxx2,nCstrata)
		       qB2 <- (HBadmRR2-HBt$res*RR2)
                   }

              }

		sss <- headstrata(dstrata-1,ndstrata)
		fff <- function(x) {
		 gtstart <- x[sss]
		 cx  <- cumsum2strata(x,S0iC2,dstrata-1,ndstrata,strataCxx2,nCstrata,gtstart)$res
		 return(cx)
		}

		### Martingale  as a function of time and for all subjects to handle strata
		if ( (!is.null(beta.iid)) | fixbeta==0) {
			EdLam0q2 <- apply(q2,2,fff)
			MGc <- q2*S0iC-EdLam0q2*c(xx2$sign)*(typexx2==1)
			MGc <- apply(MGc,2,sumstrata,xx2$id,mid+1)
		}
		if (!is.null(time) ) {
			EBdLam0q2 <- apply(qB2,2,fff)
			MGBc <- qB2*S0iC-EBdLam0q2*c(xx2$sign)*(typexx2==1)
			MGBc <- apply(MGBc,2,sumstrata,xx2$id,mid+1)
		}

	} else { MGc <- 0; MGBc <- 0}   ## }}}

	if (!is.null(time) & anyC) { MGAiid <- MGAiid+MGBc }  

	if ( (!is.null(beta.iid)) | fixbeta==0) {
		Uiid <-  (UU+MGc) %*% x$ihessian
		Uiid.naive <-  (UU) %*% x$ihessian
	} else {
		Uiid <- beta.iid
		Uiid.naive <- NULL
	} 

	if ( (!is.null(beta.iid)) | fixbeta==0) { #
		Htlast <- tailstrata(xx2$strata,xx2$nstrata)
		HtS <- Ht[Htlast,,drop=FALSE]
	} #

	## sum after id's within strata and order 
	if (!is.null(time))  {
		MGAiids <- c()
		cumhaz.time <- c()
		sus <- sort(unique(xx2$strata))
		fid <- headstrata(xx2$id,mid+1)
		wis <- xx2$strata[fid]

		for (i in sus)  { 
			ws <- 1*(wis==i)
			##
			cumhazs <- rbind(0,x$cumhaz[x$strata[x$jumps]==i,])
			cumhaz.time <- c(cumhaz.time,cpred(cumhazs,time,tminus=tminus)[,2])

			if (fixbeta==0) {
				UU <-  apply(HtS[i+1,]*t(Uiid),2,sum)
				MGAiidl <- ws*MGAiid - UU
			} else MGAiidl <- ws*MGAiid 
			MGAiids <- cbind(MGAiids,MGAiidl)
		}
		colnames(MGAiids) <- paste("strata",sus,sep="")
		names(cumhaz.time) <- paste("strata",sus,sep="")
	} else { sus <- MGAiids <- cumhaz.time <- NULL }

	if (!is.null(x$call.id)) {
		MGAiids <- nameme(MGAiids,x$name.id)
		Uiid    <- nameme(Uiid,x$name.id)
	}

	if (inherits(x,c("cifreg","recreg"))) {
		out <- list(time=time,base.iid=MGAiids,nstrata=xx2$nstrata, beta.iid=Uiid,
			    strata.call=x$strata.call,id=xx2$id,call.id=x$call.id,
			    coef=coef(x),cumhaz=x$cumhaz,cumhaz.strata=x$strata[x$jumps],
			    cumhaz.time=cumhaz.time,strata.time=sus,
			    nstrata=x$nstrata,strata.name=x$strata.name,strata.level=x$strata.level,
			    formula=x$formula,Ut=U)
	} else {
		out <- list(time=time,base.iid=MGAiid,id=xx2$id,beta.iid=Uiid,beta.iid.naive=Uiid.naive, MGt=UU,MGc=MGc,Ut=U,EdLam0=EdLam0,cumhaz=cumhaz)
	}
	return(out)
}  ## }}}

##' @export
iidBaseline.recreg <- function(object,time=NULL,ft=NULL,fixbeta=NULL,beta.iid=object$iid,tminus=FALSE,...)
{ ## {{{
	if (is.null(object$cox.prep)) stop("must call cifreg/recreg with cox.prep=TRUE\n")
	out <- IIDrecreg(object$cox.prep,object,time=time,fixbeta=fixbeta,beta.iid=beta.iid,adm.cens=object$adm.cens,tminus=tminus,...)
	out$design <- object$design
	return(out)
}  ## }}}

##' @export
GLprediid <- function(...)
{
	out <- FGprediid(...,model="GL")
	return(out)
}

##' @export
IC.recreg <- function(x,time=NULL,...) {
	if (!is.null(time)) {
		res <- iidBaseline(x,time=time,...)$base.iid
		return(res*NROW(res))
	}
	res <- with(x, iid * NROW(iid))
	return(res)
}


##' @export
plot.recreg <- function(x,se=FALSE,ylab=NULL,...) { #
	if (inherits(x,"recreg") & is.null(ylab)) ylab <- "Mean number"
	if (!se) baseplot(x,se=se,ylab=ylab,...)
	else {
		warning("Standard errors approximative (but too small), use predict and type='cumhaz' \n")
		baseplot(x,se=se,ylab=ylab,...)
	}
} #

##' @export
predict.recreg <- function(object,newdata,se=FALSE,times=NULL,np=50,...) { #
	if (!se) out <- predict.phreg(object,newdata,se=se,times=times,...)
	else {
		out <- predictrecreg(object,newdata,times=times,np=np,...)
	}
	class(out) <- c("predictrecreg",class(object)[1])
	return(out)
} #

##' @export
summary.predictrecreg <- function(object,times=NULL,strata=NULL,type=c("cif","cumhaz","surv")[2],np=10,...) { ## {{{
	if (!is.null(times)) {
		indexcol <- predictCumhaz(c(0,object$times),times,return.index=TRUE)
	} else {
		## all predictions
		nl <- length(object$times)+1 
		indexcol  <- seq(1,nl,length=np)
		times <- c(0,object$times)[indexcol]
		print("summary.predictrecreg")
		print(times)
	}

	out <- object[[type[1]]]
	nlower <- paste(type[1],".lower",sep="")
	nupper <- paste(type[1],".upper",sep="")
	lower <- object[[nlower]]
	upper <- object[[nupper]]
	nse <-  paste("se.",type[1],sep="")
	se.out  <- object[[paste("se.",type[1],sep="")]]
	if (type[1]=="surv") {
		out <- cbind(1,out) 
		if (!is.null(se.out)) se.out <- cbind(0,se.out)
		if (!is.null(lower)) lower <- cbind(1,lower) 
		if (!is.null(upper)) upper <- cbind(1,upper) 
	} else { 
		out <- cbind(0,out)
		if (!is.null(se.out)) se.out <- cbind(0,se.out)
		if (!is.null(lower)) lower <- cbind(0,lower) 
		if (!is.null(upper)) upper <- cbind(0,upper) 
	}
	if (length(lower)>1) { se <- 1; } else  { se <- 0; lower <- upper <- NULL}

	if (!is.null(lower)) ret <- list(pred=out[,indexcol],se.pred=se.out[,indexcol],lower=lower[,indexcol],upper=upper[,indexcol],times=times)
	else  ret <- list(pred=out[,indexcol],times=times)
	rownames(ret) <- NULL
	###ret$strata <- object$strata; ret$X <- object$X; ret$RR <- object$RR
	class(ret) <- "summarypredictrecreg"
	return(ret)
} ## }}}


##' @export
plot.predictrecreg <- function(x,se=FALSE,ylab=NULL,type="cumhaz",...)
{ ## {{{
	if (inherits(x,"predictrecreg") & is.null(ylab)) ylab <- "Mean number"
	plotpredictphreg(x,se=se,ylab=ylab,type=type[1],...)
} ## }}}

##' @export
predictrecreg <- function(x,newdata,times=NULL,individual.time=FALSE,tminus=FALSE,conf.type="log",conf.int=0.95,np=50,...)
{ ## {{{
	if (!inherits(x,c("cifreg","recreg","phreg")))  stop("only for phreg/recreg/cifreg models\n")
	se <- TRUE

	if (is.null(times))  {
		if (is.null(np)) times <- x$cumhaz[,1] else 
			times <- quantile(x$cumhaz[,1],probs=seq(0,1,length=np))
	} 
	des <- readPhreg(x,newdata)

	if (x$p>0)  {
		RRj <- RR <- c(exp(des$X %*% x$coef))
		Xj <- X <- RR*des$X
	} else { RRj <- RR <- rep(1,length(des$strata)); Xj <- NULL}

	surv <- surv.upper <- surv.lower <- cift <- cif.upper <- cif.lower <- cumhaz <- se.cumhaz <- base.lower <- base.upper <- se.cif <- se.surv <- c()
	j <- 1
	for (tt in times) {
		bt <- iidBaseline(x,time=tt,tminus=tminus)
		covv <- crossprod(with(bt,cbind(base.iid,beta.iid)))
		if (individual.time)  {
			baset <- bt$cumhaz.time[des$strata[j]+1]
			Xbase <- 1*outer(des$strata[j],0:(bt$nstrata-1),"==")
			RRj <- RR[j]; 
			if (x$p>0) Xj <- X[j,,drop=FALSE]*baset else Xj <- NULL
		} else {  
			baset <- bt$cumhaz.time[des$strata+1]
			Xbase <- 1*outer(des$strata,0:(bt$nstrata-1),"==")
			if (x$p>0) Xj <- X*baset else Xj <- NULL
			RRj <- RR
		}
		j <- j+1
		Xall <- cbind(RRj*Xbase,Xj)
		###   
		seLamt <- apply((Xall %*% covv)* Xall,1,sum)^.5
		Lamt <- baset*RRj
		F1t <- 1-exp(-Lamt)
		St <- exp(-Lamt)
		seF1t <- St*seLamt
		seSt <- St*seLamt
		###   
		ciLam <- conftype(Lamt,seLamt,conf.type=conf.type[1],restrict="positive",conf.int=conf.int)
		ciSt <- conftype(St,seSt,conf.type=conf.type[1],restrict="prob",conf.int=conf.int)
		ciF1 <- conftype(F1t,seF1t,conf.type=conf.type[1],restrict="prob",conf.int=conf.int)

		surv <- cbind(surv,St)
		se.surv <- cbind(se.surv,seSt)
		se.cif <- cbind(se.cif,seF1t)
		surv.upper <- cbind(surv.upper,ciSt$upper)
		surv.lower <- cbind(surv.lower,ciSt$lower)
		cift <- cbind(cift,F1t)
		cif.upper <- cbind(cif.upper,ciF1$upper)
		cif.lower <- cbind(cif.lower,ciF1$lower)
		cumhaz <- cbind(cumhaz,Lamt)
		se.cumhaz <- cbind(se.cumhaz,seLamt)
		base.upper  <- cbind(base.upper,ciLam$upper)
		base.lower  <- cbind(base.lower,ciLam$lower)
	}

	out <- list(times=times,surv=surv,surv.upper=surv.upper,surv.lower=surv.lower,
		    cumhaz=cumhaz,se.cumhaz=se.cumhaz,cif=cift,cif.upper=cif.upper,cif.lower=cif.lower,
		    cumhaz.upper=base.upper,cumhaz.lower=base.lower,strata=des$strata,X=des$X,RR=RR,
		    se.cif=se.cif,se.surv=se.surv)

	return(out)
} ## }}}

##' @export
recregIPCW <- function(formula,data=data,cause=1,cens.code=0,death.code=2,
		       cens.model=~1,km=TRUE,times=NULL,beta=NULL,offset=NULL,type=c("II","I"),
		       marks=NULL,weights=NULL,model=c("exp","lin"),no.opt=FALSE,augmentation=NULL,method="nr",se=TRUE,...)
{ ## {{{
   ## method=c("incIPCW","IPCW","rate")
   estimator=c("incIPCW")
   cl <- match.call()
   ## {{{ reading design 
    m <- match.call(expand.dots = TRUE)[1:3]
    des <- proc_design(
        formula,
        data = data,
        specials = c("offset", "weights", "cluster","marks"),
        intercept = TRUE
    )
    Y <- des$y
    if (!inherits(Y, c("Event", "Surv"))) {
        stop("Expected a 'Surv' or 'Event'-object")
    }
    if (ncol(Y) == 2) {
        exit <- Y[, 1]
        entry <- rep(0, nrow(Y))
        status <- Y[, 2]
    } else {
        entry <- Y[, 1]
        exit <- Y[, 2]
        status <- Y[, 3]
    }
    X <- des$x
    des.weights <- des$weights
    des.offset  <- des$offset
    des.marks <- des$marks
    id      <- des$cluster

  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  ### possible handling of id to code from 0:(antid-1)
  call.id <- id 
  conid <- construct_id(id,nrow(X),namesX=rownames(X))
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
  ## id before time-sorting later 
  orig.id <- id

  ## take offset and weight first from formula, but then from arguments
  if (is.null(des.offset)) {
	  if (is.null(offset)) offset <- rep(0,length(exit)) 
  } else offset <- des.offset
  if (is.null(des.weights)) {
	  if (is.null(weights)) weights <- rep(1,length(exit)) 
  } else weights <- des.weights
  if (!is.null(des.marks) & is.null(marks))  marks <- des.marks
  if (is.null(marks)) marks <- rep(1,length(id))
# }}}

	### setting up with artificial names
	data$status__ <-  status 
	data$id__ <-  id
	## lave Countcause
	data <- count.history(data,status="status__",id="id__",types=cause,multitype=TRUE)
	data$Count1__ <- data[,paste("Count",cause[1],sep="")]
	data$death__ <- (status %in% death.code)*1
	data$entry__ <- entry 
	data$exit__ <- exit 
	data$marks__ <- marks 
	statusC <- data$statusC__ <- (status %in% cens.code)*1
	data$status__cause <- (status %in% cause)*1
	data$rid__ <- revcumsumstrata(rep(1,length(entry)),id,nid)
	dexit <- exit
	dstatus <- status
	## to define properly 
	Dtime <- NULL

	formC <- update.formula(cens.model,Surv(entry__,exit__,statusC__)~ .+cluster(id__))
	cr <- phreg(formC,data=data,no.opt=TRUE,no.var=1)
	whereC <- which(status %in% cens.code)

	if (length(whereC)>0) {
		### censoring weights
		strat <- cr$strata[cr$jumps]
		x <- cr
		xx <- x$cox.prep
		S0i2 <- S0i <- rep(0,length(xx$strata))
		S0i[xx$jumps+1] <-  1/x$S0
		S0i2[xx$jumps+1] <- 1/x$S0^2
		## survival at t- to also work in competing risks situation
		if (!km) { 
			cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
			St      <- exp(-cumhazD)
		} else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
	} else  St <- rep(1,length(exit))
	Gc <- St

	###formD <- as.formula(Surv(entry__,exit__,death__)~cluster(id__))
	form1L <- as.formula(Surv(entry__,exit__,status__cause)~Count1__+death__+statusC__+cluster(id__))
	xr <- phreg(form1L,data=data,no.opt=TRUE,no.var=1,Z=matrix(data$marks__,ncol=1))
	###  xr0 <- phreg(form1,data=data,no.opt=TRUE)
	###  dr <- phreg(formD,data=data,no.opt=TRUE,no.var=1)
	###  clgl  <- recurrentMarginal(xr0,dr)
	###  plot(clgl)

	####  First partitioned estimator everywhere n^-1 sum_i \int_0^t m_i(s)/G_c(s) Y_i(s) I(D_i > s) dN_i(s) 
	x <- xr
	xx <- xr$cox.prep
	marksxx <- xx$Z[,1]
	jump1 <- xx$jumps+1
	timeJ <- xx$time[jump1]
	strataN1J <- xx$strata[jump1]
	### Partitioned estimator , same as Ghosh-Lin+Lawless-Cook estimator
	cumhazP <- c(cumsum(marksxx[jump1]/Gc[jump1])/nid)
	cumhazP <- cbind(timeJ,cumhazP)

	if (is.null(times)) stop("time for recurrent events regression must be given\n")

	### setting up regression setting with Y(t) =\int_0^t 1/G(s) dN_i(s)
	if (estimator[1]=="incIPCW") 
		Ydata <- Y <- sumstrata(marksxx*(xx$status!=0)*(xx$time<times)/Gc,xx$id,nid)
	else if (estimator[1]=="IPCW")  {
		obs <- (exit<=time & (!statusC)) | (exit>=time)/Gc
		Ydata <- Y <- sumstrata(marksxx*(xx$status!=0)*(xx$time<times),xx$id,nid)*obs
	} else {
		obs <- (exit<=time & (!statusC)) | (exit>=time)/Gc
		NtD <- sumstrata(marksxx*(xx$status!=0)*(xx$time<times),xx$id,nid)
		Ydata <- Y <- obs*NtD/pmin(times,Dtime)
	}
	nevent <- sum((xx$status!=0)*(xx$time<times))

	lastrecord <- which(data$rid__==1)
	##  1 record per subject 
	Ydata <- Y 
	###
	if (is.null(offset)) offset <- rep(0, length(exit))
	if (is.null(weights)) weights <- rep(1, length(exit))
	###  
	Xorig <- X <- as.matrix(X)
	px <- ncol(X)
	if (is.null(augmentation))  augmentation=rep(0,px)

	## order after id, similar to Y
	oid <- order(id[lastrecord])
	Xdata <- X <- X[lastrecord,,drop=FALSE][oid,,drop=FALSE]
	offset <- offset[lastrecord][oid]
	weights <- weights[lastrecord][oid]
	status <- status[lastrecord][oid]
	exit <- exit[lastrecord][oid]
	idR <- id[lastrecord][oid]
	X2 <- .Call("vecMatMat", X, X)$vXZ
	ph <- 1
	if (is.null(beta)) beta <- rep(0,ncol(X))
	## take iid vession of data 
	dataiid <- data[lastrecord,][oid,]

	if (type[1]=="II") { 
		Gcdata <- suppressWarnings(predict(cr,data,times=dexit,individual.time=TRUE,se=FALSE,km=km,tminus=TRUE)$surv)
		Gcdata[Gcdata<0.000001] <- 0.00001
		Hst <- Y[data$id__ +1]-cumsumstratasum((dexit<times)*marks*(dstatus %in% cause)/Gcdata,data$id__,nid)$lagsum
		#data$Hst <- cumsumstratasum(marks*(dstatus %in% cause)/Gcdata,data$id__,nid)$lagsum
		data$Hst <- Hst
		formC <- as.formula(paste("Surv(entry__,exit__,statusC__)~+1"))
		desform <- update.formula(cens.model,as.formula(paste("~ Hst+ . + cluster(id__)")))
		formC[[3]] <- desform[[2]]

		resC <- phreg(formC,data=data,no.opt=TRUE,no.var=1,Z=Xorig)
		xx <- resC$cox.prep
		Xt <- xx$Z
		S0i2 <- S0i <- rep(0, length(xx$strata))
		S0i[xx$jumps+1] <- 1/resC$S0
		S0i2[xx$jumps+1] <- 1/resC$S0^2
		E <- matrix(0, nrow(xx$X), 1)
		E[xx$jumps+1,] <- resC$E*((resC$jumptimes<times)*1)
		btime <- c(1 * (xx$time < times))
		EdLam0 <- apply(E*c(S0i)*btime,2,cumsumstrata,xx$strata,xx$nstrata)

		mid <- max(idR)
		MGt <- Xt*c(E[, drop = FALSE] - EdLam0 * c(xx$sign))*c(xx$weights)
		MGtiid <- apply(MGt,2,sumstrata,xx$id,mid+1)
		augmentation <- augmentation+apply(MGtiid,2,sum)
		###        
		EXt  <-  apply(Xt*c(xx$sign),2,revcumsumstrata,xx$strata,xx$nstrata)
		IEXhYtdLam0 <- apply(EXt*c(E)*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
		U <- matrix(0,nrow(xx$X),ncol(X))
		U[xx$jumps+1,] <- (resC$jumptimes<times)*E[xx$jumps+1]*EXt[xx$jumps+1,]/c(resC$S0)
		MGt2 <- (U[,drop=FALSE]-IEXhYtdLam0*c(xx$sign))*c(xx$weights)
		###
		MGCiid2 <- apply(MGt2,2,sumstrata,xx$id,mid+1)
		MGCiid2 <- MGtiid-MGCiid2
	}  else MGCiid2 <- 0   
	#

	obj <- function(pp, all = FALSE) { ## {{{
		lp <- c(X %*% pp + offset)
		if (model[1] == "exp") p <- exp(lp) else p <- lp
		if (model[1] == "dexp") p <- exp(lp) 
		if (model[1] == "dexp") {
			Dlogl <- weights * X *p * c(Y - p)
			D2logl <- c(weights) *p^2* X2
		} else {
			Dlogl <- weights *  X * c(Y - p)
			if (model[1]=="exp") D2logl <- c(weights) * c(p)* X2 else D2logl <- c(weights) * X2
		}
		if (model[1] == "exp") ploglik <- 0 else ploglik <- sum(weights * (Y - p)^2)
		D2log <- apply(D2logl, 2, sum)
		gradient <- apply(Dlogl, 2, sum) + augmentation
		hessian <- matrix(D2log, length(pp), length(pp))

		if (all) {
			ploglik <- sum(weights * (Y - p)^2)
			ihess <- solve(hessian)
			beta.iid <- Dlogl %*% ihess
			beta.iid <- apply(beta.iid, 2, sumstrata, idR, max(id) + 1)
			robvar <- crossprod(beta.iid)
			val <- list(par = pp, ploglik = ploglik, gradient = gradient,
				    hessian = hessian, ihessian = ihess, id = idR,
				    Dlogl = Dlogl, iid = beta.iid, robvar = robvar,
				    var = robvar, se.robust = diag(robvar)^0.5)
			return(val)
		}
		structure(-ploglik, gradient = -gradient, hessian = hessian)
	} ## }}}

	## setting default for NR 
	dots <- list(...)
	if (length(dots)==0) {
		if (model[1]=="exp") control <- list(tol=1e-10,stepsize=0.5)  
		else control <- NULL
	} else control <- dots[[1]]

	p <- ncol(X)
	opt <- NULL
	if (p > 0) {
		if (no.opt == FALSE) {
			if (tolower(method) == "nr") {
				tim <- system.time(opt <- lava::NR(beta, obj,control=control))
				opt$timing <- tim
				opt$estimate <- opt$par
			}
			else {
				opt <- nlm(obj, beta, ...)
				opt$method <- "nlm"
			}
			cc <- opt$estimate
			###            if (!se) return(cc)
			val <- c(list(coef = cc), obj(opt$estimate, all = TRUE))
		}
		else val <- c(list(coef = beta), obj(beta, all = TRUE))
	}
	else {
		val <- obj(0, all = TRUE)
	}
	if (length(val$coef) == length(colnames(X))) names(val$coef) <- colnames(X)

	val <- c(val, list(times = times, Y=Y, ncluster=nid, nevent=nevent,n=length(exit),X=X))


	if (se) {
		Gcdata <- suppressWarnings(predict(cr,data,times=dexit,individual.time=TRUE,se=FALSE,km=km,tminus=TRUE)$surv)
		Gcdata[Gcdata<0.000001] <- 0.00001
		## check data sorted in dexit 
		if (type[1]!="II") Hst <- Y[data$id__ +1]-cumsumstratasum((dexit<times)*marks*(dstatus %in% cause)/Gcdata,data$id__,nid)$lagsum
		data$Hst <- Hst
		if (model[1]=="dexp") HstX <-c(exp(as.matrix(Xorig) %*% val$coef))*Xorig*c(data$Hst) else HstX <- Xorig*c(data$Hst) 
		ccn <- paste("nn__nn",1:ncol(Xorig),sep="")
		colnames(HstX) <- ccn
		nncovs <- c()
		for (i in 1:ncol(Xorig)) nncovs <- c(paste(nncovs,paste("+",ccn[i],sep="")))
		formC <- as.formula(paste("Surv(entry__,exit__,statusC__)~+1"))
		desform <- update.formula(cens.model,as.formula(paste("~",nncovs,"+ . + cluster(id__)")))
		formC[[3]] <- desform[[2]]

		data <- cbind(data,HstX)
		resC <- phreg(formC,data=data,no.opt=TRUE,no.var=1)
		xx <- resC$cox.prep
		S0i2 <- S0i <- rep(0, length(xx$strata))
		S0i[xx$jumps + 1] <- 1/resC$S0
		S0i2[xx$jumps + 1] <- 1/resC$S0^2
		E <- U <- matrix(0, nrow(xx$X), ncol(X))
		E[xx$jumps + 1, ] <- resC$E*((resC$jumptimes < times)*1)
		btime <- c(1 * (xx$time < times))
		EdLam0 <- apply(E*c(S0i)*btime,2,cumsumstrata,xx$strata,xx$nstrata)
		MGt <- (E[, drop = FALSE] - EdLam0 * c(xx$sign) )*c(xx$weights)
		MGCiid <- apply(MGt, 2, sumstrata, xx$id, max(id) + 1)
		#MGCiid <- MGCiid-MGCiid2
		MGCiid <- MGCiid+MGCiid2
	} else  MGCiid <- 0 


	val$call <- cl
	val$formula <- formula
	val$model <- model[1]
	val$model.type <- model[1]
	val$Y <- Ydata
	val$X <- Xdata
	val$id <- orig.id
	val$call.id <- call.id
	val$nid <- nid
	val$name.id <- name.id
	val$iid.naive <- val$iid
	val$naive.var <- val$var
	if (se) { 
		val$MGCiid <- MGCiid
		MGCiid <- MGCiid %*% val$ihessian 
		val$iid <- val$iid + MGCiid
	} else val$MGCiid <- MGCiid 
	if (is.matrix(val$iid)) 
		if (length(name.id)==nrow(val$iid)) {
			rownames(val$iid) <- name.id
			oid <- order(name.id)  
			val$iid <- val$iid[oid,,drop=FALSE]
		}
	if (is.matrix(val$MGCiid))  {
		if (length(name.id)==nrow(val$MGCiid)) rownames(val$MGCiid) <- name.id
		val$MGCiid <- val$MGCiid[oid,,drop=FALSE]
	}

	robvar <- crossprod(val$iid)
	val$var <- val$robvar <- robvar
	val$se.robust <- diag(robvar)^0.5
	val$se.coef <- diag(val$var)^0.5
	val$cens.code <- cens.code
	val$cause <- cause
	val$death.code <- death.code
	val$model.estimator <- estimator
	val$augmentation <- augmentation
	val$type <- type
	val$cumhazP <- cumhazP
	val$nevent <- nevent
	val$design <- des
	class(val) <- c("binreg", "resmean")
	return(val)
} ## }}}

strataAugment <- survival::strata

##' Simulation of two-stage recurrent events data based on Cox/Cox or Cox/Ghosh-Lin structure 
##'
##' Simulation of two-stage recurrent events data based on Cox/Cox or Cox/Ghosh-Lin structure. type=3 will generate
##' Cox/Cox twostage mode, type=2 will generate Ghosh-Lin/Cox model. 
##' If the variance is var.z=0, then generates data without any dependence or frailty. If model="twostage" then default is to generate data from Ghosh-Lin/Cox model, and
##' if type=3 then will generate data with marginal Cox models (Cox/Cox). 
##' Simulation based on linear aproximation of hazard for two-stage models based on grid on time-scale. Must be sufficientyly fine. 
##'
##' Must specify baselines of recurrent events and terminal event and possible covariate effects.
##'
##' @param n number of id's 
##' @param base1 baseline for cox/ghosh-lin models
##' @param drcumhaz baseline for terminal event
##' @param var.z variance of gamma frailty 
##' @param r1 relative risk term for baseline 
##' @param rd relative risk term for terminal event 
##' @param rc relative risk term for censorings
##' @param fz possible transformation (function) of frailty term 
##' @param fdz possible transformation (function) of frailty term for death 
##' @param model twostage, frailty, shared (partly shared two-stage model)
##' @param type type of simulation, default is decided based on model
##' @param cens to right censor
##' @param share to fit patly shared random effects model
##' @param cens censoring rate for exponential censoring
##' @param nmin default 100, at least nmin or number of rows of the two-baselines max(nmin,nrow(base1),nrow(drcumhaz)) points in time-grid from 0 to maximum time for base1
##' @param nmax default 1000, at most nmax points in time-grid 
##' @references 
##' Scheike (2025), Two-stage recurrent events random effects models, LIDA, to appear
##' @export
simGLcox <- function(n,base1,drcumhaz,var.z=0,r1=NULL,rd=NULL,rc=NULL,fz=NULL,fdz=NULL,
		     model=c("twostage","frailty","shared"),type=NULL,share=1,cens=NULL,nmin=100,nmax=1000)
{ ## {{{
	## setting up baselines for simulations 
	maxt <- tail(base1[,1],1)
	base1 <- as.matrix(base1); drcumhaz <- as.matrix(drcumhaz)
	nmin <- max(nrow(base1),nrow(drcumhaz),nmin)
	nmin <- min(nmax,nmin)
	seqt <- seq(from=0,to=maxt,length.out=nmin)
	if (base1[1,1]!=0) base1 <- rbind(0,base1) 
	if (drcumhaz[1,1]!=0) drcumhaz <- rbind(0,drcumhaz) 
	base1 <- cbind(seqt, lin.approx(seqt,base1))
	cumD <- cbind(seqt, lin.approx(seqt,drcumhaz))
	###
	St <- exp(-cumD[,2])
	Stm <- cbind(base1[,1],St)
	###
	dbase1 <- diff(c(0,base1[,2]))
	dcum <- cbind(base1[,1],dbase1)
	maxtime <- tail(base1[,1],1)

	if (is.null(r1)) r1 <- rep(1,n)
	if (is.null(rd)) rd <- rep(1,n)
	if (is.null(rc)) rc <- rep(1,n)

	fz.orig <- fz
	if (is.null(fz)) fz <- function(x) x

	if (var.z[1]>0) {
		z1 <- z <- rgamma(n,share/var.z[1])*var.z[1] 
		if (share<1) { 
			z2 <- rgamma(n,(1-share)/var.z[1])*var.z[1] 
			z <- z+z2
		} 
		fzz <- fz(z1)
		if (share<1) fzz <- fzz/share; 
		mza <- mean(fzz)
		if (n<10000 & (!is.null(fz.orig))) {
			zl <- rgamma(100000,share/var.z[1])*var.z[1] 
			fzl <- fz(z)
			mza <- mean(fzl)
		} 
	}  else fzz <- z <- z1 <- rep(1,n)

	if (var.z[1]==0) model <- "frailty"
	if (is.null(type)) 
		if (model[1]=="twostage") type <- 2 else type <- 1
	## for frailty setting we also consider any function of z 
	if (!is.null(fdz)) { fdzz <- fdz(z); rd <- rd*fdzz; z <- rep(1,n);}

	## survival censoring given X, Z, either twostage or frailty-model 
	if (type>=2) stype <- 2 else stype <- 1
	if (var.z[1]==0) stype <- 1
	### dd <- .Call("_mets_simSurvZ",as.matrix(rbind(c(0,1),Stm)),rd,z,var.z[1],stype)
	dd <- .Call("_mets_simSurvZ",as.matrix(rbind(Stm)),rd,z,var.z[1],stype)
	dd <- data.frame(time=dd[,1],status=(dd[,1]<maxtime))
	if (!is.null(cens)) cens <- rexp(n)/(rc*cens) else cens <- rep(maxtime,n)
	dd$status <- ifelse(dd$time<cens,dd$status,0)
	dd$time <- pmin(dd$time,cens)

	## to avoid R check error
	reverseCountid  <-  death  <- NULL

	if (model[1]=="multiplicative") {
		## other random effect 
		z2 <- rgamma(n,share/var.z[2])*var.z[2] 
		fzz <- z1*z2
		type <- 3
	}
	## type=2 draw recurrent process given X,Z with rate:
	##  Z exp(X^t beta_1) d \Lambda_1(t)/S(t|X,Z) 
	## such that GL model holds with exp(X^t beta_1) \Lambda_1(t)
	## type=3, observed hazards on Cox form among survivors
	## twostage shared<1: W_1 ~ N1, W_1+W_2 ~ D observed hazards on Cox form among survivors
	## twostage share=1: or W_1 ~ N1, W_1~ D   observed hazards on Cox form among survivors
	## multiplicatve:       W_2 * W_1 ~ N1, W_1~ D   observed hazards on Cox form among survivors
	dcum <- cbind(base1[,1],dbase1)
	### ll <- .Call("_mets_simGL",as.matrix(rbind(0,dcum)),c(1,St),r1,rd,z1,fzz,dd$time,type,var.z[1],nmax,1)
	ll <- .Call("_mets_simGL",as.matrix(dcum),c(St),r1,rd,z1,fzz,dd$time,type,var.z[1],nmax,1)
	colnames(ll) <- c("id","start","stop","death")
	ll <- data.frame(ll)
	ll$death <- dd$status[ll$id+1]
	## add frailty to data for possible validation
	ll$z <- z1[ll$id+1]
	ll$fz <- fzz[ll$id+1]
	## add counts of id
	ids <- countID(ll)
	ll <- cbind(ll,ids[,c(2,4,5)]); 
	ll$status <- 1; 
	ll <- dtransform(ll,status=0,reverseCountid==1)
	ll$statusD <- ll$status
	ll <- dtransform(ll,statusD=3,reverseCountid==1 & death==1)

	attr(ll,"base1events") <- base1
	attr(ll,"deathcumbase") <- cumD
	return(ll)
} ## }}}


simGLcoxC <- function(n,base1,drcumhaz,var.z=0,r1=NULL,rd=NULL,rc=NULL,fz=NULL,fdz=NULL,
		      model=c("twostage","frailty","shared","multiplicative"),type=NULL,share=1,cens=NULL,nmax=200,by=1)
{ ## {{{
	## setting up baselines for simulations 
	maxt <- tail(base1[,1],1)
	seqt <- seq(by,maxt,by=by)
	base1 <- predictCumhaz(rbind(0, as.matrix(base1)), seqt)
	cumD <- predictCumhaz(rbind(0, as.matrix(drcumhaz)), seqt)
	###
	St <- 1-cumD[,2]
	Stm <- cbind(base1[,1],St)
	###
	dbase1 <- diff(c(0,base1[,2]))
	dcum <- cbind(base1[,1],dbase1)
	maxtime <- tail(base1[,1],1)

	if (is.null(r1)) r1 <- rep(1,n)
	if (is.null(rd)) rd <- rep(1,n)
	if (is.null(rc)) rc <- rep(1,n)

	fz.orig <- fz
	if (is.null(fz)) fz <- function(x) x

	if (var.z[1]>0) {
		z1 <- z <- rgamma(n,share/var.z[1])*var.z[1] 
		if (share<1) { 
			z2 <- rgamma(n,(1-share)/var.z[1])*var.z[1] 
			z <- z+z2
		} 
		fzz <- fz(z1)
		if (share<1) fzz <- fzz/share; 
		mza <- mean(fzz)
		if (n<10000 & (!is.null(fz.orig))) {
			zl <- rgamma(100000,share/var.z[1])*var.z[1] 
			fzl <- fz(z)
			mza <- mean(fzl)
		} 
	}  else fzz <- z <- z1 <- rep(1,n)

	if (var.z[1]==0) model <- "frailty"
	if (is.null(type)) 
		if (model[1]=="twostage") type <- 2 else type <- 1
	## for frailty setting we also consider any function of z 
	if (!is.null(fdz)) { fdzz <- fdz(z); rd <- rd*fdzz; z <- rep(1,n);}

	## survival censoring given X, Z, either twostage or frailty-model 
	if (type>=2) stype <- 2 else stype <- 1
	cumHazS <- -log(St)
	dd <- rchaz(cumHazS,rd)
	dd <- data.frame(time=dd[,1],status=(dd[,1]<maxtime))
	if (!is.null(cens)) cens <- rexp(n)/(rc*cens) else cens <- rep(maxtime,n)
	dd$status <- ifelse(dd$time<cens,dd$status,0)
	dd$time <- pmin(dd$time,cens)

	## to avoid R check error
	reverseCountid  <-  death  <- NULL
	## type=2 draw recurrent process given X,Z with rate:
	##  Z exp(X^t beta_1) d \Lambda_1(t)/S(t|X,Z) 
	## such that GL model holds with exp(X^t beta_1) \Lambda_1(t)
	## type=3, observed hazards on Cox form among survivors
	## W_1 ~ N1, W_1+W_2 ~ D observed hazards on Cox form among survivors
	## or W_1 ~ N1, W_1~ D   observed hazards on Cox form among survivors
	## or W_2 * W_1 ~ N1, W_1~ D   observed hazards on Cox form among survivors
	dcum <- cbind(base1[,1],dbase1)
	ll <- .Call("_mets_simGL",as.matrix(rbind(0,dcum)),c(1,St),r1,rd,z1,fzz,dd$time,type,var.z[1],nmax,1)
	colnames(ll) <- c("id","start","stop","death")
	ll <- data.frame(ll)
	ll$death <- dd$status[ll$id+1]
	## add frailty to data for possible validation
	ll$z <- z1[ll$id+1]
	ll$fz <- fzz[ll$id+1]
	## add counts of id
	ids <- countID(ll)
	ll <- cbind(ll,ids[,c(2,4,5)]); 
	ll$status <- 1; 
	ll <- dtransform(ll,status=0,reverseCountid==1)
	ll$statusD <- ll$status
	ll <- dtransform(ll,statusD=3,reverseCountid==1 & death==1)

	return(ll)
} ## }}}


simRecurrentCox <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,X=NULL,r1=NULL,r2=NULL,rd=NULL,rc=NULL, 
			    model=c("not-random","random"),frailty=TRUE,var.z=0.5,death.code=3,alpha=1,...)
{ ## {{{
	if (is.null(r1)) r1 <- rep(1,n)
	if (is.null(r2)) r2 <- rep(1,n)
	if (is.null(rd)) rd <- rep(1,n)
	if (is.null(rc)) rc <- rep(1,n)

	## to avoid error 
	ctime <- death <- NULL

	################################################################
	### approximate hazards to make marginals fit (approximately)
	################################################################
	###laplace and inverse laplace of gamma
	lap<-function(theta,t) { return( (1+t/theta)^(-theta)) }
	ilap<-function(theta,t) { itheta<-1/theta; return((t^(-itheta)-1)/(itheta)) }

	## addapt to make recurrent mean on cox form with this baseline
	base1 <- cumhaz
	if (is.null(death.cumhaz)) stop("Modification for death in this function otherwise just use simRecurrentII\n")
	if (is.null(X)) stop("X must be given to link with simulated data\n"); 

	### Cox baseline 
	St <- exp(-cpred(rbind(c(0,0),death.cumhaz),base1[,1])[,2])

	## unique hazard combinations Death Relative Risk
	rds <- unique(rd)
	if (model[1]=="random") {
		z <- rgamma(n,1/var.z)*var.z 
		if (n<10000) {
			zl <- rgamma(100000,1/var.z)*var.z 
			mza <- mean(zl^alpha)
		} else mza <- mean(z^alpha)
		rds <- rd 
	}  else z <- rep(1,n)
	dtt <- diff(c(0,base1[,1]))
	dbase1 <- diff(c(0,base1[,2]))
	data <- c()

	XX <- c()
	nks <- 0
	k <- 1
	if (model[1]!="random") {
		for (rdss in rds) {
			lam1ms <- (dbase1)/St^rdss
			where <- which(rdss==rd)
			nk <- length(where)
			cumhaz1 <- cbind(base1[,1],cumsum(lam1ms))
			LamDr <- scalecumhaz(death.cumhaz,rdss) 

			datss <- simRecurrentII(nk,cumhaz1,cumhaz2,death.cumhaz=LamDr,
						r1=r1[where],r2=NULL,rd=NULL,rc=rc[where],...)
			Xs <- X[where,,drop=FALSE][datss$id,,drop=FALSE]
			XX <- rbind(XX,Xs)
			datss$id <- datss$id+nks
			nks <- nks+nk
			data <- rbind(data,as.matrix(datss))
		}
	}

	if (model[1]=="random") {
		for (k in 1:n) {
			rdss <- rd[k]
			if (frailty==TRUE) lam1ms <- (z[k]^alpha/mza)*(dbase1)/St^(z[k]*rdss) else {
				Stt <-  exp(-z[k]*ilap(1/var.z,St^rdss))
				lam1ms <- (z[k]^alpha/mza)*(dbase1)/Stt
			}
			where <- k
			nk <- 1
			cumhaz1 <- cbind(base1[,1],cumsum(lam1ms))

			if (!frailty) LamDr <- cbind(base1[,1],-log(Stt)) 
			if (frailty) LamDr <- scalecumhaz(death.cumhaz,z[k]*rdss) 

			datss <- simRecurrentII(nk,cumhaz1,cumhaz2,death.cumhaz=LamDr,
						r1=r1[where],r2=NULL,rd=NULL,rc=rc[where],...)
			Xs <- X[where,,drop=FALSE][datss$id,,drop=FALSE]
			XX <- rbind(XX,Xs)
			datss$id <- datss$id+nks
			nks <- nks+nk
			data <- rbind(data,as.matrix(datss))
		}
	}

	rownames(XX) <- NULL
	rownames(data) <- NULL
	colnames(XX) <- paste("X",1:ncol(X),sep="")
	data <- cbind(data,XX)
	data <- as.data.frame(data)

	dsort(data) <- ~id+entry+time
	data$revnr <- revcumsumstrata(rep(1,nrow(data)),data$id-1,n)
	data$statusD <- data$status
	data <- dtransform(data,statusD=death.code,death==1)

	return(data)
} ## }}}


simMarginalMeanCox <- function(n,cens=3/5000,k1=0.1,k2=0,bin=1,Lam1=NULL,Lam2=NULL,LamD=NULL,
			       beta1=rep(0,2),betad=rep(0,2),betac=rep(0,2),X=NULL,...)
{ ## {{{

	### to avoid R-check error
	revnr <- death <- status <- NULL
	p <- length(beta1)

	if (is.null(X))  {
		if (bin==1) X <- matrix(rbinom(n*p,1,0.5),n,p) else  X <- matrix(rnorm(n*p),n,p)
		colnames(X) <- paste("X",1:p,sep="")
	}
	r1 <- exp( X %*% beta1)
	rd <- exp( X %*% betad)
	rc <- exp( X %*% betac)

	if (is.null(Lam2)) Lam2 <- Lam1; 

	rr <- simRecurrentCox(n,scalecumhaz(Lam1,k1),cumhaz2=scalecumhaz(Lam1,k2),death.cumhaz=LamD,X=X,cens=cens,r1=r1,rd=rd,rc=rc,...)

	if (bin==0) dcut(rr,breaks=4) <- X1g~X1 else rr$X1g <- rr$X1
	if (bin==0) dcut(rr,breaks=4) <- X2g~X2 else rr$X2g <- rr$X2
	return(rr)
}

##' @export
scalecumhaz <- function(cumt,k)
{
	return( t(t(cumt)*c(1,k)))
}

##' @export
GLprediid <- function(...)
{
	out <- FGprediid(...,model="GL")
	return(out)
} ## }}}


boottwostageREC <- function(margsurv,recurrent,data,bootstrap=100,id="id",stepsize=0.5,...) 
{ ## {{{
	n <- max(margsurv$id)
	K <- bootstrap
	formid <- as.formula(paste("~",id))
	rrb <- blocksample(data, size = n*K, formid)
	rrb$strata <- floor((rrb[,id]-0.01)/n)

	outb <- outd <- outr <- c()
	for (i in 1:K)
	{
		rrbs <- subset(rrb,strata==i-1)
		drb <- phreg(margsurv$formula,data=rrbs)
		if (inherits(recurrent,"recreg")) {
			xrb <- recreg(recurrent$formula,data=rrbs,
				      cause=recurrent$cause,death.code=recurrent$death.code,cens.code=recurrent$cens.code,twostage=TRUE)
		} else xrb <- phreg(recurrent$formula,data=rrbs)
		outbl <- tryCatch(twostageREC(drb,xrb,rrbs,numderiv=0,control=list(stepsize=stepsize,iter.max=10),...),error=function(x) NULL)
     if (!is.null(outbl)) outb <- rbind(outb,outbl$coef)
     outd <- rbind(outd,coef(drb))
     outr <- rbind(outr,coef(xrb))
  }
  varb <- cov(outb)
  vard <- cov(outd)
  varr <- cov(outr)

  list(outb=outb,var=varb,se=diag(varb)^.5,se.coxD=diag(vard)^.5,se.coxR=diag(varr)^.5)
} ## }}}

##' Fittting of Two-stage recurrent events random effects model based on Cox/Cox or Cox/Ghosh-Lin marginals 
##'
##' Fittting of Two-stage recurrent events random effects model based on Cox/Cox or Cox/Ghosh-Lin marginals. Random
##' effects model fore recurrent events with terminal  event. Marginal models fitted first and given to twostageREC function.
##'
##' @param margsurv marginal model for terminal event 
##' @param recurrent marginal model for recurrent events
##' @param data used for fitting
##' @param theta starting value for total variance of gamma frailty
##' @param model can fully shared "full", partly shared "shared", or non-shared where the random effect acts only on recurrent events
##' @param ghosh.lin to force use ghosh.lin marginals based on recurrent (taking baseline and coefficients) 
##' @param theta.des regression design for variance
##' @param var.link possible link  function 1 is exponential link
##' @param method NR
##' @param no.opt to not optimize
##' @param weights possible weights
##' @param se.cluster  to combine influence functions for naive variance based on these clusters GEE style
##' @param fnu a function to make transformation for nu (amount shared)
##' @param nufix to fix the amount shared
##' @param nu starting value for amount shared
##' @param numderiv uses numerical derivatives for some derivatives
##' @param derivmethod method for numerical derivative
##' @param ... arguments for 
##' @references 
##' Scheike (2025), Two-stage recurrent events random effects models, LIDA, to appear
##' @export
twostageREC  <-  function (margsurv,recurrent, data = parent.frame(), theta = NULL, model=c("full","shared","non-shared"),ghosh.lin=NULL,
  theta.des = NULL, var.link = 0, method = "NR", no.opt = FALSE, weights = NULL, se.cluster = NULL, 
  fnu=NULL,nufix=0,nu=NULL,numderiv=1,derivmethod=c("simple","Richardson"),...)
{ ## {{{
    newrec <- 1; at.risk <- 1
    if (!inherits(margsurv, "phreg")) stop("Must use phreg for death model\n")
    if (!inherits(recurrent, "phreg")) stop("Must use phreg for recurrent model\n")
    if (is.null(recurrent$cox.prep) & newrec==0) stop("recreg must be called with cox.prep=TRUE\n")
    if (inherits(recurrent, "recreg") & is.null(ghosh.lin))  ghosh.lin <- 1 
    if ((!inherits(recurrent, "recreg")) & is.null(ghosh.lin))  ghosh.lin <- 0 
    if (is.null(recurrent$cox.prep.twostage) & newrec==1 & (ghosh.lin==1)) stop("recreg must be called with twostage=TRUE\n")
    clusters <- margsurv$cox.prep$id
    n <-  max(clusters)+1
    if (is.null(theta.des) == TRUE) ptheta <- 1
    if (is.null(theta.des) == TRUE) theta.des <- matrix(1, n, ptheta) else theta.des <- as.matrix(theta.des)
    ptheta <- ncol(theta.des)
    if (nrow(theta.des) != n) stop("Theta design does not have correct dim")
    if (is.null(theta) == TRUE) {
        if (var.link == 1) theta <- rep(0, ptheta)
        if (var.link == 0) theta <- rep(1, ptheta)
    }
    if (is.null(nu)  & (nufix==0)) { if (is.null(fnu))  nu <-  0.5 else nu <- fnu[[1]](0.5) }
    if (length(theta) != ptheta) theta <- rep(theta[1], ptheta)
    if (length(nu) != ptheta) nu <- rep(nu,ptheta)
    theta.score <- rep(0, ptheta)
    Stheta <- var.theta <- matrix(0, ptheta, ptheta)
    max.clust <- length(unique(clusters))
    theta.iid <- matrix(0, max.clust, ptheta)
    xx <- margsurv$cox.prep
    xr <- recurrent$cox.prep
    if (ghosh.lin & newrec==1) {
	    xr <- recurrent$cox.prep.twostage
    } 
    nn <- length(xx$strata)
    if (is.null(weights)) weights <- rep(1, nn)
    if (length(weights) != nn) stop("Weights do not have right length")
    statusxx <- rep(0, length(xx$strata))
    statusxx[xx$jumps + 1] <- 1
    xx$status <- statusxx
    mid <- max(xx$id) + 1
    Nsum <- cumsumstratasum(statusxx, xx$id, mid, type = "all")
    Ni.tau <- sumstrata(statusxx, xx$id, mid)
    S0i2 <- S0i <- rep(0, length(xx$strata))
    S0i[xx$jumps + 1] <- 1/margsurv$S0
    cumhazD <- cumsumstrata(S0i, xx$strata, xx$nstrata)
    if (!is.null(margsurv$coef)) RR <- exp(xx$X %*% margsurv$coef) else RR <- rep(1, nn)
    HD <- c(cumhazD * RR)
    cumDYt <- sumstrata(HD*xx$sign,xx$id,mid)
    statusx1 <- rep(0, length(xx$strata))
    statusx1[xr$jumps + 1] <- 1
    xr$status <- statusx1
    N1sum <- cumsumstratasum(statusx1, xr$id, mid, type = "all")
    N1i.tau <- sumstrata(statusx1, xr$id, mid)
    S01i2 <- S01i <- rep(0, length(xr$strata))
    S01i[xr$jumps + 1] <- 1/recurrent$S0
    cumhaz1 <- cumsumstrata(S01i, xr$strata, xr$nstrata)
    if (!is.null(recurrent$coef)) RR1 <- exp(xr$X %*% recurrent$coef) else RR1 <- rep(1, nn)
    H1 <- c(cumhaz1 * RR1)
    ###
    ## designs of fixed time covariates and weights
    cc <- cluster.index(xx$id)
    firstid <- cc$firstclustid + 1
    if (max(cc$cluster.size) == 1) stop("No clusters !, maxclust size=1\n")
    ###
    theta.des <- theta.des[xx$id+ 1, , drop = FALSE]
    theta.des <- theta.des[xx$ord + 1, , drop = FALSE]
    weightsid <- weights <- weights[xx$ord + 1]
    weights <- weights[firstid]
    thetaX <- as.matrix(theta.des[firstid, , drop = FALSE])
    Xdeath <- as.matrix(xx$X[firstid,,drop = FALSE])
    Xrecurrent <- as.matrix(xr$X[firstid,,drop = FALSE])
    statusxb <- statusxx+statusx1
    ###
    Nsumb <- Nsum$lagsum+N1sum$lagsum
    rd <- RR[firstid]
    r1 <- RR1[firstid]
    lastid <- tailstrata(xx$id,mid)
    ###
    cumDL <- HD[lastid]
    nuX <- thetaX
    ###
    idD <- xx$id[statusxx==1]

    obj <- function(par, all = FALSE) {
        if (var.link == 1) epar <- c(exp(c(par))) else epar <- c(par)
        thetav <- c(as.matrix(theta.des) %*% c(epar))
        thetai <- thetav[firstid]
	###
###        tildeL <- .Call("_mets_tildeLambda1",S01i,cumhazD,r1,rd,thetai,xx$id)
###	if (at.risk==1)  
###		tildeL <- apply(tildeL*c(xr$sign),2,cumsumstrata,xr$id,mid)
        tildeL <- .Call("_mets_tildeLambda1R",S01i,cumhazD,r1,rd,thetai+1*ghosh.lin,xr$id,xr$sign)
        tildeLast <- tildeL[lastid,]
	Ht <- thetav*tildeL[,1]+exp(thetav*HD)
	Hr <- thetai*tildeLast[,1]+exp(thetai*cumDL)
	DHt <- tildeL[,1]+thetav*tildeL[,2]+HD*exp(thetav*HD)
	DHr <- tildeLast[,1]+thetai*tildeLast[,2]+cumDL*exp(thetai*cumDL)
        ###
	D2Ht <- 2*tildeL[,2]+thetav*tildeL[,3]+HD^2*exp(thetav*HD)
	D2Hr <- 2*tildeLast[,2]+thetai*tildeLast[,3]+cumDL^2*exp(thetai*cumDL)
	aHt <- abs(Ht)
	aDHt <- abs(DHt)
	aD2Ht <- abs(D2Ht)
        ###
        l1 <- sumstrata(log(1 + thetav * N1sum$lagsum) * statusxb, xx$id, mid)
        l2 <- sumstrata(statusxb * HD, xx$id, mid)
	if (at.risk==0) l3 <- -(1/thetai + N1i.tau) * log(Hr) 
	if (at.risk==1) l3 <- - sumstrata(c(xr$sign)*(1/thetav+N1sum$sum)*log(aHt),xr$id,mid) 
	HtD <- Ht[statusxx==1]
	DHtD <- DHt[statusxx==1]
	D2HtD <- D2Ht[statusxx==1]
        l4 <- -sumstrata(log(HtD),idD,mid)
        logliid <- (l1 + thetai * l2 + l3 + l4) * c(weights)
        logl <- sum(logliid)
        ploglik <- logl
        ###
        l1s <- sumstrata(N1sum$lagsum/(1 + thetav * N1sum$lagsum) * statusxb, xx$id, mid)
	if (at.risk==0) l3s <- -(1/thetai + N1i.tau) * DHr/Hr + log(Hr)/thetai^2;
	if (at.risk==1) l3s <- sumstrata(c(xr$sign)*(-(1/thetav+N1sum$sum)*DHt/Ht+log(aHt)/thetav^2),xr$id,mid);
	l4s <-  -sumstrata((DHtD/HtD),idD, mid) 
	Dltheta <- (l1s+l2+l3s+l4s)*c(weights) 
	scoreiid <- thetaX * c(Dltheta)
        D2N <- -sumstrata(N1sum$lagsum^2/(1 + thetav * N1sum$lagsum)^2 * statusxb, xx$id, mid)
	if (at.risk==0) 
        Dhes <- (2/thetai^2) * DHr/Hr -(1/thetai + N1i.tau)*(D2Hr*Hr-DHr^2)/Hr^2 - (2/thetai^3) * log(Hr) 
	if (at.risk==1) 
        Dhes <- sumstrata(c(xr$sign)* (
	(2/thetav^2) * aDHt/aHt -(1/thetav + N1sum$sum)*(aD2Ht*Ht-aDHt^2)/aHt^2 -(2/thetav^3) * log(aHt)) ,xr$id,mid)
	D2l4 <- -sumstrata(((D2HtD*HtD-DHtD^2)/HtD^2), idD, mid)
        Dhes <- c(Dhes+D2N+D2l4) * c(weights)
        if (var.link == 1) {
            scoreiid <- scoreiid * c(thetai)
            Dhes <- Dhes * thetai^2 + thetai * Dltheta
        }
        gradient <- apply(scoreiid, 2, sum)
        hessian <- crossprod(thetaX, thetaX * c(Dhes))
        hess2 <- crossprod(scoreiid)
        val <- list(id = xx$id, score.iid = scoreiid, logl.iid = logliid,
            ploglik = ploglik, gradient = gradient, hessian = hessian,
            hess2 = hess2)
        if (all)
            return(val)
        with(val, structure(-ploglik, gradient = -gradient, hessian = -hessian))
    }


    objNonShared <- function(par, all = FALSE) {
        if (var.link == 1) epar <- c(exp(c(par))) else epar <- c(par)
        thetav <- c(as.matrix(theta.des) %*% c(epar))
        thetai <- thetav[firstid]
	###
	Ht <- 1+thetav*H1 
	Hr <- Ht[lastid]
	DHt <- H1
	DHr <- H1[lastid]
        ###
	D2Ht <- 0
	D2Hr <- 0
	aHt <- abs(Ht)
	aDHt <- abs(DHt)
	aHr <- abs(Hr)
        ###
        l1 <- sumstrata(log(1 + thetav * N1sum$lagsum) * statusx1, xx$id, mid)
	if (at.risk==0) l3 <- -(1/thetai + N1i.tau) * log(Hr) 
	if (at.risk==1) l3 <- - sumstrata(c(xr$sign)*(1/thetav+N1sum$sum)*log(aHt),xr$id,mid) 
        logliid <- (l1 + l3) * c(weights)
        logl <- sum(logliid)
        ploglik <- logl
        ###
        l1s <- sumstrata(N1sum$lagsum/(1 + thetav * N1sum$lagsum) * statusx1, xx$id, mid)
	if (at.risk==0) l3s <- -(1/thetai + N1i.tau) * DHr/Hr + log(Hr)/thetai^2;
	if (at.risk==1) l3s <- sumstrata(c(xr$sign)*(-(1/thetav+N1sum$sum)*DHt/Ht+log(aHt)/thetav^2),xr$id,mid);
	Dltheta <- (l1s+l3s)*c(weights) 
	scoreiid <- thetaX * c(Dltheta)
        D2N <- -sumstrata(N1sum$lagsum^2/(1 + thetav * N1sum$lagsum)^2 * statusx1, xx$id, mid)
	if (at.risk==0) 
        Dhes <- (2/thetai^2) * DHr/Hr -(1/thetai + N1i.tau)*(-DHr^2)/Hr^2 - (2/thetai^3) * log(Hr) 
	if (at.risk==1) 
        Dhes <- sumstrata(c(xr$sign)* (
	(2/thetav^2) * aDHt/aHt -(1/thetav + N1sum$sum)*(-aDHt^2)/aHt^2 -(2/thetav^3) * log(aHt)) ,xr$id,mid)
        Dhes <- c(Dhes+D2N) * c(weights)
        if (var.link == 1) {
            scoreiid <- scoreiid * c(thetai)
            Dhes <- Dhes * thetai^2 + thetai * Dltheta
        }
        gradient <- apply(scoreiid, 2, sum)
        hessian <- crossprod(thetaX, thetaX * c(Dhes))
        hess2 <- crossprod(scoreiid)
        val <- list(id = xx$id, score.iid = scoreiid, logl.iid = logliid,
            ploglik = ploglik, gradient = gradient, hessian = hessian,
            hess2 = hess2)
        if (all)
            return(val)
        with(val, structure(-ploglik, gradient = -gradient, hessian = -hessian))
    }


    ## default is simple identity 
###    if (is.null(fnu)) { fw <- function(x) 1/(1+exp(x)); Dfw <- function(x) -exp(x)/(1+exp(x))^2;} else { 
   if (is.null(fnu)) { fw <- function(x) x; Dfw <- function(x) 1;} else { fw <- fnu[[1]]; Dfw <- fnu[[2]]; } 
   nudes <- theta.des  
   p <- ncol(theta.des)

    objShared <- function(par, all = FALSE) {
        if (var.link == 1) epar <- c(exp(par[1:p])) else epar <- c(par[1:p])
        thetav <- c(as.matrix(theta.des) %*% epar)
	if (nufix==1) nu1 <- c(as.matrix(nudes) %*% nu)
	else nu1 <- c(as.matrix(nudes) %*% par[(p+1):2*p])
        nu1i <-nu1[firstid]
	tbeta1 <- fw(nu1); tbeta2 <- 1-tbeta1; 
        thetai <- thetav[firstid]; tbeta1i <- tbeta1[firstid]; tbeta2i <- tbeta2[firstid]
	###
	R <-  exp( - thetav*HD);  DR <- -HD*exp( - thetav* HD); D2R <-  HD^2*exp( - thetav* HD) 
        tildeL <- .Call("_mets_tildeLambda1R",S01i,cumhazD,r1,rd,thetai+1*ghosh.lin,xr$id,xr$sign)
	tildeLast <- tildeL[lastid,]
	Ht <- (thetav/tbeta1)*tildeL[,1]+exp(thetav*HD)
	Hr <- (thetai/tbeta1i)*tildeLast[,1]+exp(thetai*cumDL)
	DHt <- (tildeL[,1]+thetav*tildeL[,2])/tbeta1+HD*exp(thetav*HD)
	DHr <- (tildeLast[,1]+thetai*tildeLast[,2])/tbeta1i+cumDL*exp(thetai*cumDL)
	DHtn <- -(thetav/tbeta1^2)*tildeL[,1]
	DHrn <- -(thetai/tbeta1i^2)*tildeLast[,1]
        ###
	D2Ht <- (2*tildeL[,2]+thetav*tildeL[,3])/tbeta1+HD^2*exp(thetav*HD)
	D2Hr <- (2*tildeLast[,2]+thetai*tildeLast[,3])/tbeta1i+cumDL^2*exp(thetai*cumDL)
        ###
        N <- (tbeta1 + thetav * N1sum$lagsum)/Ht + tbeta2*R
        l1d <- sumstrata(log(N) * statusxx, xx$id, mid)
        l11 <- sumstrata(log(tbeta1 + thetav * N1sum$lagsum) * statusx1, xx$id, mid)
	l1 <- l1d+l11
        l2 <- sumstrata(statusxb * HD, xx$id, mid) 
	l22 <-  -log(tbeta1i)*N1i.tau 
        if (at.risk==0) l3 <- -(tbeta1i/thetai + N1i.tau)*log(Hr) 
        if (at.risk==1) l3 <- - sumstrata(c(xr$sign)*(tbeta1/thetav+N1sum$sum)*log(Ht),xr$id,mid) 
        l4 <-  (tbeta1i)*cumDYt
        logliid <- (l1 + thetai*l2+ l22 + l3 + l4) * c(weights)
        logl <- sum(logliid)
        ploglik <- logl
        ###
	DN1tt <- (Ht*N1sum$lagsum-DHt*(tbeta1+thetav*N1sum$lagsum))
	DN1t <-  DN1tt/Ht^2 
	DNt <- c(DN1t+tbeta2*DR)
        l1ds <- sumstrata((DNt/N)*statusxx,xx$id,mid)
        l11s <- sumstrata((N1sum$lagsum/(tbeta1 + thetav * N1sum$lagsum)) * statusx1, xx$id, mid)
	l1s <- l1ds+l11s
        if (at.risk==0) 
        l3s <- -(tbeta1i/thetai + N1i.tau) * (DHr/Hr) + (tbeta1i/thetai^{2}) * log(Hr) 
        if (at.risk==1) 
	l3s <- sumstrata(c(xr$sign)*(-(tbeta1/thetav+N1sum$sum)*DHt/Ht+log(Ht)*tbeta1/thetav^2),xr$id,mid);
        Dltheta <- (l1s+l2+l3s)*c(weights)
        ### 
	DNtn <- (Ht-DHtn*(tbeta1+thetav*N1sum$lagsum))
	DNn <- DNtn/Ht^2 
        l1dn <- sumstrata( c(DNn-R)/N * statusxx, xx$id, mid)
        l11n <- sumstrata(1/(tbeta1 + thetav * N1sum$lagsum) * statusx1, xx$id, mid)
	l1n <- l1dn+l11n
	l22n <-  -N1i.tau/tbeta1i 
	if (at.risk==0) 
		l3n <- -(tbeta1i/thetai+N1i.tau)*(DHrn/Hr)-(1/thetai)*log(Hr)
	if (at.risk==1) 
        l3n <- -sumstrata(c(xr$sign)*(+(tbeta1/thetav+N1sum$sum)*(DHtn/Ht)+(1/thetav)*log(Ht)), xr$id,mid);
        Dlnu <- (l1n+l22n+l3n+cumDYt)*c(weights)
        ###
	if (nufix==1)
        scoreiid <- thetaX * c(Dltheta)
        else  scoreiid <- cbind(thetaX * c(Dltheta),nuX*c(Dlnu)*Dfw(nu1i))
        ###   ###
        Dl11s <- -sumstrata(N1sum$lagsum^2/(tbeta1 + thetav * N1sum$lagsum)^2 * statusxb, xx$id, mid)
        Dl3s <- (2*tbeta1i/thetai^2) * DHr/Hr -(tbeta1i/thetai+ N1i.tau)*(D2Hr*Hr-DHr^2)/Hr^2 - (2*tbeta1i/thetai^3) * log(Hr) 
        ###
	D2N1t <- (Ht^2*(DHt*N1sum$lagsum-D2Ht*(tbeta1+thetav*N1sum$lagsum)-DHt*N1sum$lagsum)-2*DHt*DN1tt)/Ht^4 
	D2Nt <-  D2N1t+tbeta2*D2R
	Dl1ds <- -sumstrata(((D2Nt*N-DNt^2)/N^2)*statusxx, xx$id, mid)
        Dhes <- c(Dl1ds+Dl11s+Dl3s) * c(weights)
        if (var.link == 1) {
            scoreiid[,1:p] <- scoreiid[,1:p] * c(thetai)
            Dhes <- Dhes * thetai^2 + thetai * Dltheta
        }
        gradient <- apply(scoreiid, 2, sum)
        hessian <- crossprod(thetaX, thetaX * c(Dhes))
	if (nufix==0) {
###	    hessiann <- crossprod(nuX, nuX* c(D2nu))
###	    hessiannp <- crossprod(thetaX,nuX*Dtn)
###	    hessian <- cbind(hessian,hessianp)
###	    hessian <- rbind(hessian,cbind(t(hessianp),hessiann))
	}
        hess2 <- crossprod(scoreiid)
        val <- list(id = xx$id, score.iid = scoreiid, logl.iid = logliid, ploglik = ploglik, gradient = gradient, hessian = -hess2, hess2 = hess2)
        if (all)
            return(val)
###        with(val, structure(-ploglik, gradient = -gradient, hessian = -hessian))
        with(val, structure(-ploglik, gradient = -gradient, hessian = hess2))
    }


   if (model[1]=="shared") obj <- objShared
   if (model[1]=="non-shared") obj <- objNonShared
   if (nufix==0 & model[1]=="shared") par <- c(theta,nu) else par <- theta

    opt <- NULL
    if (no.opt == FALSE) {
        if (tolower(method) == "nr") {
            opt <- lava::NR(par, obj,...)
            opt$estimate <- opt$par
        }
        else {
            opt <- nlm(obj, par,...)
            opt$method <- "nlm"
        }
        cc <- opt$estimate
        val <- c(list(coef = cc), obj(opt$estimate, all = TRUE))
    }
    else val <- c(list(coef = par), obj(par, all = TRUE))
    val$score <- val$gradient
    theta <- matrix(c(val$coef), length(c(val$coef)), 1)

    if (!is.null(colnames(theta.des)))
        thetanames <- colnames(theta.des)
    else thetanames <- paste("dependence", 1:length(c(theta)), sep = "")
    if (nufix==0 & model[1]=="shared") thetanames <-
     thetanames <- c(paste("dependence", 1:p, sep = ""),paste("share", 1:p,sep = ""))

    if (length(thetanames) == length(c(theta))) {
        rownames(theta) <- thetanames
        names(val$coef)  <- thetanames
    }

    if (numderiv==1 & model[1]=="shared") {
	    dobj <- function(p) {
		    oo <- obj(p)
		    return(attr(oo,"gradient"))
	    }
	    hessian <- numDeriv::jacobian(dobj,val$coef,method=derivmethod[1])
	    val$hessian <- -hessian
    }

    hessianI <- solve(val$hessian)
    val$theta.iid.naive <- val$score.iid %*% hessianI

    if (!is.null(se.cluster))
        if (length(se.cluster) != length(clusters))
            stop("Length of seclusters and clusters must be same\n")
    if (!is.null(se.cluster)) {
        iids <- unique(se.cluster)
        nseclust <- length(iids)
        if (is.numeric(se.cluster))
            se.cluster <- fast.approx(iids, se.cluster) - 1
        else se.cluster <- as.integer(factor(se.cluster, labels = seq(nseclust))) - 1
        val$theta.iid <- apply(val$theta.iid,2,sumstrata,se.cluster, nseclust)
        val$theta.iid.naive <- apply(val$theta.iid.naive,2,sumstrata,se.cluster, nseclust)
    }
    var <- robvar.theta <- var.theta <- crossprod(val$theta.iid)
    naive.var <- crossprod(val$theta.iid.naive)
    val <- c(val, list(theta = theta, var.theta = var,n=mid,p=ncol(thetaX),var.link=var.link,
                       robvar.theta = var, var = var, thetanames = thetanames,
                       model = model[1], se = diag(var)^0.5), 
	     var.naive = naive.var,no.opt=no.opt,ghosh.lin=ghosh.lin)
    class(val) <- "twostageREC"
    attr(val, "clusters") <- clusters
    attr(val, "fnu") <- fw
    attr(val, "Dfnu") <- Dfw
    attr(val, "secluster") <- c(se.cluster)
    attr(val, "var.link") <- var.link
    attr(val, "ptheta") <- ptheta
    attr(val, "n") <- n
    attr(val, "response") <- "survival"
    attr(val, "additive-gamma") <- 0
    attr(val, "twostage") <- "two.stage"
    return(val)
} ## }}}

##' @export
summary.twostageREC <- function(object,vcov=NULL,delta=0,...) { ## {{{
    I <- -solve(object$hessian)
    if (!is.null(vcov)) V <- vcov else V <- object$var
    ncluster <- object$n
    cc <- estimate(coef=object$coef,vcov=V)$coefmat
    pd <- object$p
    if (object$var.link==1 & object$model=="full") f <- function(p) exp(p)
    if (object$var.link==0 & object$model=="full") f <- function(p) p
    if (object$var.link==1 & object$model=="non-shared") f <- function(p) exp(p)
    if (object$var.link==0 & object$model=="non-shared") f <- function(p) p

    if (object$var.link==1 & object$model=="shared") f <- function(p) c(exp(p[1]),attr(object,"fnu")(p[2]))
    if (object$var.link==0 & object$model=="shared") f <- function(p) c(p[1],attr(object,"fnu")(p[2]))
    if (delta==1 | object$model=="full" | object$model=="non-shared") 
    expC <- lava::estimate(coef=object$coef,vcov=V,f=f)$coefmat ##[,c(1,3,4),drop=FALSE]
    else expC <- apply(cc[,c(1,3,4),drop=FALSE],2,f) 
  n <- object$n
  res <- list(coef=cc,n=n,nevent=object$nevent,ncluster=ncluster,var=V,exp.coef=expC,var.link=object$var.link,
	      ghosh.lin=object$ghosh.lin)
  class(res) <- "summary.twostageREC"
  res
}  ## }}}

##' @export
print.summary.twostageREC  <- function(x,max.strata=5,...) { ## {{{
  if (x$ghosh.lin==0) cat("Cox(recurrent)-Cox(terminal) intensity model\n"); 
  if (x$ghosh.lin==1) cat("Ghosh-Lin(recurrent)-Cox(terminal) mean model\n"); 
  if (!is.null(x$ncluster)) cat("\n ", x$ncluster, " clusters\n",sep="")
  if (!is.null(x$coef)) {
    cat("coeffients:\n")
    printCoefmat(x$coef,...)
    cat("\n")
    if (x$var.link==1) cat("var=exp(coeffients),shared:\n")
    if (x$var.link==0) cat("var,shared:\n")
    printCoefmat(x$exp.coef,...)
  }
  cat("\n")
}  ## }}}


##' @export
print.twostageREC  <- function(x,...) {
  print(summary(x),...)
} 

