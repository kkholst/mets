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
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ## data with no ties
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' Lam1 <- base1cumhaz;  Lam2 <- base4cumhaz;  LamD <- drcumhaz
##' ## simulates recurrent events of types 1 and 2 and with terminal event D and censoring
##' rr <- simRecurrentII(100,Lam1,cumhaz2=Lam2,death.cumhaz=LamD,cens=3/5000)
##' rr <- count.history(rr)
##' rr$cens <- 0
##' nid <- max(rr$id)
##' rr$revnr <- revcumsumstrata(rep(1,nrow(rr)),rr$id-1,nid)
##' rr$x <- rnorm(nid)[rr$id]
##' rr$statusG <- rr$status
##' rr <- dtransform(rr,statusG=3,death==1)
##' dtable(rr,~statusG+status+death)
##' dcut(rr) <- gx~x
##'
##' ll <- recreg(Event(start,stop,statusG)~x+cluster(id),data=rr,cause=1,death.code=3)
##' summary(ll)
##' 
##' ## censoring stratified after quartiles of x
##' lls <- recreg(Event(start, stop, statusG)~x+cluster(id),data=rr,cause=1,
##'               death.code=3,cens.model=~strata(gx))
##' summary(lls)
##' 
##' @aliases IIDbaseline.recreg strataAugment scalecumhaz GLprediid recregIPCW twostageREC simGLcox
##' @export
recreg <- function(formula,data,cause=1,death.code=c(2),cens.code=0,cens.model=~1,weights=NULL,offset=NULL,Gc=NULL,wcomp=NULL,
		   augmentation.type=c("lindyn.augment","lin.augment"),...)
{# {{{
outi  <- recregB(formula,data,cause=cause,death.code=death.code,cens.code=cens.code,cens.model=cens.model,weights=weights,offset=offset,Gc=Gc,wcomp=wcomp,...)

if (!is.null(outi$lindyn.augment)) {
outA  <- recregB(formula,data,cause=cause,death.code=death.code,cens.code=cens.code,cens.model=cens.model,weights=weights,offset=offset,Gc=Gc,wcomp=wcomp,
		 augmentation=outi[[augmentation.type[1]]],...)
outi <- outA
}
return(outi)
}# }}}

recregB <- function(formula,data=data,cause=1,death.code=c(2),cens.code=0,cens.model=~1,weights=NULL,offset=NULL,Gc=NULL,wcomp=NULL,...)
{# {{{
    cl <- match.call()# {{{
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster","offset","strataAugment")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (inherits(Y,"EventCens")) stop("Change to Event call, see example, EventCens disabled")
    if (!inherits(Y,"Event")) stop("Expected a 'Event'-object")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- rep(0,nrow(Y))
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

###   if (!is.null(stratapos <- attributes(Terms)$specials$strataAugment)) {
###    ts <- survival::untangle.specials(Terms, "strataAugment")
###    Terms  <- Terms[-ts$terms]
###    strataAugment <- as.numeric(m[[ts$vars]])-1
###    strataAugment.name <- ts$vars
###  }  else { strataAugment <- NULL; strataAugment.name <- NULL}

    if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
        ts <- survival::untangle.specials(Terms, "offset")
        Terms  <- Terms[-ts$terms]
        offset <- m[[ts$vars]]
    }
    X <- model.matrix(Terms, m)
    if (!is.null(intpos  <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

    ## }}}

    res <- c(
        recreg01(data, X, entry, exit, status,
            id = id, strata = strata, offset = offset, weights = weights,
            cens.model = cens.model, cause = cause, strata.name = strata.name, strataA = NULL, ## strataAugment,
            death.code = death.code, cens.code = cens.code, Gc = Gc, wcomp = wcomp, ...
        ),
        list(call = cl, model.frame = m, formula = formula, strata.pos = pos.strata, cluster.pos = pos.cluster, n = nrow(X), nevent = sum(status %in% cause))
    )
    colnames(res$iid) <- names(res$coef)

    class(res) <- c("recreg", "phreg")
    return(res)
}# }}}

##' @export IIDbaseline.recreg 
IIDbaseline.recreg <- function(x,time=NULL,fixbeta=NULL,...)
{# {{{
   return(IIDbaseline.cifreg(x,time=time,fixbeta=fixbeta,...))
} # }}}


##' @export
IC.recreg <- function(x, ...) {
  res <- with(x, iid * NROW(iid))
  return(res)
}
recreg01 <- function(data,X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,strataA=NULL,
          strata.name=NULL,beta,stderr=1,method="NR",no.opt=FALSE, propodds=NULL,profile=0,
          case.weights=NULL,cause=1,death.code=2,cens.code=0,Gc=NULL,cens.model=~+1,augmentation=NULL,
	  cox.prep=FALSE,wcomp=NULL,augment.model=NULL,ftime.augment=NULL,
	  adm.cens.code=NULL,adm.cens.time=NULL,...) { # {{{
# {{{ setting up weights, strata, beta and so forth before the action starts
    p <- ncol(X)
    if (missing(beta)) beta <- rep(0,p)
    if (p==0) X <- cbind(rep(0,length(exit)))
    augmentation.call <- augmentation
    if (is.null(augmentation)) augmentation <- 0

    cause.jumps <- which(status %in% cause)
    max.jump <- max(exit[cause.jumps])
###    other <- which((status %in% death.code ) )
    other <- which((status %in% death.code ) & (exit< max.jump))

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

   orig.strataA <- strataA
   if (is.null(strataA)) {  strataA <- rep(0,length(exit)); 
              nstrataA <- 1; strataA.level <- NULL; 
   } else {
	  strataA.level <- levels(strataA)
	  ustrataA <- sort(unique(strataA))
	  nstrataA <- length(ustrataA)
	  strataA.values <- ustrataA
      if (is.numeric(strataA)) strataA <-  fast.approx(ustrataA,strataA)-1 else  {
      strataA <- as.integer(factor(strataA,labels=seq(nstrataA)))-1
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

    if (!is.null(id)) {
        ids <- unique(id)
        nid <- length(ids)
        if (is.numeric(id))
            id <-  fast.approx(ids,id)-1
        else  {
            id <- as.integer(factor(id,labels=seq(nid)))-1
        }
    } else { id <- as.integer(seq_along(entry))-1;  nid <- nrow(X); }
    ## orginal id coding into integers 1:...
    id.orig <- id+1;

    ## }}}

   ### censoring weights constructed
    whereC <- which( status %in% cens.code )
    time <- exit
    cens <- statusC <- c(status %in% cens.code )
    data$id <- id
    data$exit__ <- exit
    data$entry__ <- entry
    data$statusC <- statusC
    data$status__ <- (status %in% cause)*1
    cens.strata <- cens.nstrata <- NULL
    data <- count.history(data,status="status__",id="id",types=cause,multitype=TRUE)
    data$Nt <- data[,paste("Count",cause[1],sep="")]
 
    ## augmentation model and remove intercept
    if (!is.null(augment.model)) { XXA <- model.matrix(augment.model,data)[,-1,drop=FALSE]; namesXXA <- colnames(XXA); } else XXA <- NULL

    if (length(whereC)>0) {# {{{
    if (is.null(Gc)) {
        kmt <- TRUE
        if (inherits(cens.model,"formula")) {
            formC <- update.formula(cens.model,Surv(entry__,exit__,statusC)~ . +cluster(id))
            cens.model <- phreg(formC,data)
        }
        if (cens.model$p>0) kmt <- FALSE
        Pcens.model <- suppressWarnings(predict(cens.model,data,times=exit,individual.time=TRUE,se=FALSE,km=kmt))
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
    }# }}}

    trunc <- TRUE
    Zcall <- cbind(status,cens.strata,Stime,cens,strata,strataA,XXA) ## to keep track of status and Censoring strata
    ## setting up all jumps of type "cause", need S0, S1, S2 at jumps of "cause"
    stat1 <- 1*(status %in% cause)
    xx2 <- .Call("FastCoxPrepStrata",entry,exit,stat1,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")
    xx2$nstrata <- nstrata
    jumps <- xx2$jumps+1
    jumptimes <- xx2$time[jumps]
    strata1jumptimes <- xx2$strata[jumps]
    Xj <- xx2$X[jumps,,drop=FALSE]
    mdif <- min(diff(c(0,jumptimes)))

    ## G(T_j-) at jumps of type "cause"
    if (length(whereC)>0) {
        if (is.null(Gc)) {
            whereaJ <- fast.approx(c(0,cens.model$cumhaz[,1]),jumptimes,type="left")
            Gts <- vecAllStrata(cens.model$cumhaz[,2],cens.model$strata.jump,cens.model$nstrata)
            ### back to km product-limit form
            Gts <- apply(rbind(0,Gts),2,diff)
            GtsAl<- Gts <- suppressWarnings(apply(Gts,2,function(x) exp(cumsum(log(1-x)))))
            Gts <- rbind(1,Gts)[whereaJ,]
            Gts[is.na(Gts)] <- 0
            Gjumps <- Gts
        } else Gts <- Gjumps <- c(1,Pcens.model$surv)[fast.approx(c(0,Pcens.model$time),jumptimes,type="left")]
    } else {
        Gts <-   Gjumps <- rep(1,length(jumptimes))
    }

    ## computing terms for those experiencing another cause, need S0, S1, S2
    if (length(other)>=1) {# {{{
        trunc <- TRUE
        weightso <- weights[other]/Stime[other]
        if (is.null(adm.cens.time))
        timeoo <- rep(max(exit)+1,length(other))
        else timeoo <- adm.cens.time[other] 
        statuso <- rep(1,length(other))
        Xo <- X[other,,drop=FALSE]
        offseto <- offset[other]
	## mdif to avoid double counting for composite where death.code and cause share types
        entryo <- exit[other]+mdif/10
        ido <- id[other]
        stratao <- strata[other]
        if (nCstrata>1) {
	    Cstratao <- cens.strata[other]
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
	### gives index of timeo related to jumptimes and same strata
	### the value 0 means that jumptime has no point in time0, thus S0other=0
        where <- indexstratarightR(timeo,xx$strata,jumptimes,strata1jumptimes,nstrata)
    }# }}}

    obj <- function(pp,all=FALSE) {# {{{

        if (length(other)>=1)  {
            if (nCstrata==1) {# {{{
                rr <- c(xx$sign*exp(xx$X %*% pp + xx$offset)*xx$weights)
                S0no <- c(0,revcumsumstrata(rr,xx$strata,xx$nstrata))
                S1no  <- rbind(0,apply(xx$X*rr,2,revcumsumstrata,xx$strata,xx$nstrata))
                S2no  <- rbind(0,apply(xx$XX*rr,2,revcumsumstrata,xx$strata,xx$nstrata));
                Gjumps <- c(Gjumps)

                S0no <- Gjumps*S0no[where+1]
                S1no <- Gjumps*S1no[where+1,,drop=FALSE]
                S2no <- Gjumps*S2no[where+1,,drop=FALSE]
                ## }}}
            }  else {# {{{

                ff <- function(x,strata,nstrata,strata2,nstrata2)
                {# {{{
                    x <- rbind(0,revcumsum2strata(x,strata,nstrata,strata2,nstrata2)$mres)
                    ### take relevant S0sc (s=strata,c=cstrata) at jumptimes so that strata=s also match
                    x <- x[where+1,]
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

        if (!is.null(propodds)) {
            pow <- c(.Call("cumsumstrataPOR",weightsJ,S0,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$pow);
            DLam <-.Call("DLambetaR",weightsJ,S0,E,Xj,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$res;
            Dwbeta <- DLam*rr2now+(pow-1)*Xj
            DUadj  <- .Call("vecMatMat",Dwbeta,U,PACKAGE="mets")$vXZ
        }

        Ut <- caseweightsJ*weightsJ*U
        ## E^2, as n x (pxp)
###        Et2 <-  .Call("vecMatMat",E,E,PACKAGE="mets")$vXZ
        Et2  <- .Call("vecCPMat",E)$XX
        S2S0 <-  (S2oo+S2no)/S0
        DUt <-  -(S2S0-Et2)

        if (!is.null(propodds)) {
            Ut  <- pow*Ut
            S0 <- S0/pow
            DUt <- pow*DUt
            DUt <- DUt+DUadj
            DUt <-  .Call("XXMatFULL",DUt,p,PACKAGE="mets")$XXf
            if (profile==1) {
		Ut <- Ut+c(ploglik)*Dwbeta
		## not implemented
		DUt <- DUt
            }
            ploglik <- pow*ploglik
        }

        U  <- apply(Ut,2,sum)
        DUt <- caseweightsJ*weightsJ*DUt
###        DU <- -matrix(apply(DUt,2,sum),p,p)
	if (ncol(DUt)!=p*p) {
        DU <- matrix(0,p,p);
        DU[lower.tri(DU,diag=TRUE)] <- -apply(DUt,2,sum)
        DU<- DU+t(DU)
        diag(DU) <- diag(DU)/2
	} else  DU <- -matrix(apply(DUt,2,sum),p,p)
        ploglik <- sum(ploglik)
        U <- U+augmentation

        out <- list(ploglik=ploglik,gradient=U,hessian=-DU,cox.prep=xx2,
                    hessianttime=DUt,weightsJ=weightsJ,caseweightsJ=caseweightsJ,
                    jumptimes=jumptimes,strata=strataJ,nstrata=nstrata,S0s=cbind(S0oo,S0no),
                    time=jumptimes,S0=S0/(caseweightsJ*weightsJ),S2S0=S2S0,E=E,U=Ut,X=Xj,Gjumps=Gjumps)

        if (all)
            return(out)
        else
            with(out,structure(-ploglik, gradient=-gradient, hessian=-hessian))
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
            if (stderr==2) return(cc)
            val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
        } else val <- c(list(coef=beta),obj(beta,all=TRUE))
    } else {
	no.opt <- TRUE
        val <- obj(0,all=TRUE)
    }# }}}

    beta.s <- val$coef
    if (is.null(beta.s)) beta.s <- 0
    ## getting final S's
    opt <-  val ## obj(beta.s,all=TRUE)

    if (p>0) {
    iH <- - tryCatch(solve(opt$hessian),error=function(e) matrix(0,nrow(opt$hessian),ncol(opt$hessian)) )

    ### iid version given G_c
    ## {{{
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

    if (length(other)>=1) { ## martingale part for type-2 after T
        ### xx2 all data with start stop structure, takes position of death times
        otherxx2 <- which((xx2$Z[,1] %in% death.code) & xx2$sign==1)
        statusxx2 <- xx2$Z[,1]
        rr0 <- xx2$sign
        jumpsC <- which((xx2$Z[,1] %in% cens.code) & xx2$sign==1)
        strataCxx2 <- xx2$Z[,2]
        S0iC2  <-  S0iC <- rep(0,length(xx2$status))
        S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
        S0iC[jumpsC] <- 1/S0rrr[jumpsC]
        S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
	## Gc(t) computed 
        Gcxx2 <- exp(cumsumstrata(log(1-S0iC),strataCxx2,nCstrata))
        Gstart <- rep(1,nCstrata)

        dstrata <- mystrata(data.frame(strataCxx2,xx2$strata))
        ndstrata <- attr(dstrata,"nlevel")
        lastt <- tailstrata(dstrata-1,ndstrata)
        ll <-  cumsum2strata(Gcxx2,S0i,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)
        Htsj <- ll$res[lastt][dstrata]-ll$res

        fff <- function(x) {
            cx  <- cumsum2strata(Gcxx2,x*S0i,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)
            cx <- cx$res[lastt][dstrata]-cx$res
            return(cx)
        }
        EHtsj  <- apply(E,2,fff)
        rrx2 <- rr[otherxx2]*xx2$weights[otherxx2]/xx2$Z[otherxx2,3]
        MGt2  <- -(Z[otherxx2,,drop=FALSE]*Htsj[otherxx2,]-EHtsj[otherxx2,,drop=FALSE])*rrx2
        UU2 <- apply(MGt2,2,sumstrata,xx2$id[otherxx2],mid+1)
        UU  <-  UU+UU2
    }
    ## }}}

    if ((stderr==1) & (length(other)>=1) & (length(whereC)>0)) { ### Censoring adjustment for jumps of other type but only for KM-case {{{

        Xos <- matrix(0,nrow(Z),ncol(Z));
        Xos[otherxx2,] <- Z[otherxx2,]*rrx2
        rrx <- rep(0,nrow(Z))
        rrx[otherxx2] <- rrx2
        rrsx <- cumsumstrata(rrx,strataCxx2,nCstrata)
        Xos <- apply(Xos,2,cumsumstrata,strataCxx2,nCstrata)
        q2 <- (Xos*c(Htsj)-EHtsj*c(rrsx))

        sss <- headstrata(dstrata-1,ndstrata)
        fff <- function(x) {
            gtstart <- x[sss]
            cx  <- cumsum2strata(x,S0iC2,dstrata-1,ndstrata,strataCxx2,nCstrata,gtstart)$res
            return(cx)
        }
        EdLam0q2 <- apply(q2,2,fff)

        ### Martingale  as a function of time and for all subjects to handle strata
        MGc <- q2*S0iC-EdLam0q2*c(xx2$sign)
        MGc <- apply(MGc,2,sumstrata,xx2$id,mid+1)
        
    } else MGc <- 0 ## }}}

    Uiid <-  (UU+MGc) %*% iH
    UUiid <- UU %*% iH
    var1 <-  crossprod(UUiid)
    varmc <-  crossprod(Uiid)

    ## compute regression augmentation for censoring martingale 
    if ((!is.null(augment.model)) & (length(whereC)>0)) {## {{{

	CovZXstrata <- function(X,Ej,Z,Sign,strata,nstrata,jumps) 
	{# {{{
		strata  <- c(strata); Sign <- c(Sign)
         ###	Ej <- Ej[jumps,,drop=FALSE]; Ej <- Ej
		ZE <- apply(Z*Sign,2,revcumsumstrata,strata,nstrata)[jumps,,drop=FALSE]; 
		XZ  <- .Call("vecMatMat",X,Z)$vXZ;  
		XZ <- apply(XZ*Sign,2,revcumsumstrata,strata,nstrata)[jumps,,drop=FALSE]; 
		EXZ  <- .Call("vecMatMat",Ej,ZE)$vXZ;  
		out <- XZ-EXZ
		return(out)
	}# }}}

       ## regress U(s)=\int_s^\infty (Z-E) w(s) dM(s) on agument-model among survivors 
       ## U(s) = U(\infty) - \int_0^s (Z-E) w(s)  dM(s)
       ## sum (e_i - \bar e) U(s) Y_i(s)

       ## setup design for augmentation regression 
	###    desform <- update.formula(augment.model,~+ . + cluster(id))
	###    formC[[3]] <- desform[[2]]
	###    ## set-up censoring martingale 
	###    cr2 <- phreg(formC,data=data,no.opt=TRUE,no.var=1)
	###    xxc <- cr2$cox.prep
	###    jumpsC <- xxc$jumps+1
	###    Gcj <- exp(-cr2$cumhaz[,2])

    rr0 <- c(xx2$sign)
    jumpsC <- which((xx2$Z[,1] %in% cens.code) & xx2$sign==1)
    strataCxx2 <- xx2$Z[,2]
    S0iC2  <-  S0iC <- rep(0,length(xx2$status))
    S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
    S0iC[jumpsC] <- 1/S0rrr[jumpsC]
    S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
    S0c <- c(S0rrr[jumpsC])
    ## Gc(t) computed  as exp(- Cumhaz) to avoid some "0"s
    Gcj <- Gcxx2 <- exp(-cumsumstrata(S0iC,strataCxx2,nCstrata))[jumpsC]

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
    rr <- c(exp(xx2$X %*% beta.s+ xx2$offset)*xx2$weights)
    Zrr <- xx2$X*rr
    ZEdN <- apply(U,2,revcumsumstrata,xx2$id,nid)

    covXsZ <-   CovZXstrata(XXA,EA,Zrr,rr0,strataCxx2,nCstrata,jumpsC) 
    covXsrr <-  CovZXstrata(XXA,EA,as.matrix(rr,ncol=1),rr0,strataCxx2,nCstrata,jumpsC) 
    covXsUs3 <- .Call("vecMatMat",covXsrr,EdLam0[jumpsC,,drop=FALSE])$vXZ;  
    covXsUs2 <- covXsZ*cumhaz[jumpsC]-covXsUs3 
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
### solve(matrix(hesst[1,],3,3)) %*% matrix(covXsYs[1,],3,2)
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
   #### iid magic  for censoring augmentation martingale{{{
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
   ## scale with G_c(t) to compare with gamma
   gammat <- gammat * c(Gcj)
# }}}

   var.augment <-  varmc  -  iH %*% var.augment %*% iH
   var.augment.times <-  varmc  +  iH %*% var.augment.times %*% iH
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
   } ## }}}

    if (!is.null(orig.strataA)) { ## compute augmentation term based on this beta{{{
	xx2strataA <- xx2$Z[,6]
        Xos <- matrix(0,nrow(Z),ncol(Z));
        Xos[otherxx2,] <- Z[otherxx2,]*rrx2
        rrx <- rep(0,nrow(Z))
        rrx[otherxx2] <- rrx2
        S0A <- revcumsumstrata(rr0,xx2strataA,nstrataA)
	S0A[S0A==0] <- 1
        rrsx <- cumsumstrata(rrx,xx2strataA,nstrataA)/S0A
        Xos <- apply(Xos,2,cumsumstrata,xx2strataA,nstrataA)/c(S0A)
        q2 <- (Xos*c(Htsj)-EHtsj*c(rrsx))

        dstrataA  <- mystrata(data.frame(strataCxx2,xx2strataA))
        ndstrataA <- attr(dstrataA,"nlevel")

        sss <- headstrata(dstrataA-1,ndstrataA)
        fff <- function(x) {
            gtstart <- x[sss]
            cx  <- cumsum2strata(x,S0iC,dstrataA-1,ndstrataA,strataCxx2,nCstrata,gtstart)$res
            return(cx)
        }
        EdLam0qA2 <- apply(q2,2,fff)

        ### Martingale  as a function of time and for all subjects to handle strata
        MGAc <- q2*(S0iC!=0)-EdLam0q2*c(xx2$sign)
        MGAc <- apply(MGAc,2,sumstrata,xx2$id,mid+1)
	augment.new <- apply(MGAc,2,sum)
    } else { augment.new <- MGAc <- NULL}
# }}}

    ### end if (p>0)
    } else {
          iid.augment <- iid.augment.times <- augment <- augment.times <- NULL 
          var.augment.times <- var.augment <- NULL
          var.augment.times.iid <- var.augment.iid <- NULL
          Uiid.augment.times <- Uiid.augment <- NULL
          time.gammat <- gamma <- gammat <- NULL
          ftime.gamma <- NULL
          Gcj <- NULL
	    varmc <- var1 <- 0; augment.new <- MGAc <- MGc <- iH <- UUiid <- Uiid <- NULL}
    strata <- xx2$strata[jumps]
    cumhaz <- cbind(opt$time,cumsumstrata(1/opt$S0,strata,nstrata))
    colnames(cumhaz)    <- c("time","cumhaz")

## SE of estimator ignoring some censoring terms
if (no.opt==FALSE & p!=0) {
DLambeta.t <- apply(opt$E/c(opt$S0),2,cumsumstrata,strata,nstrata)
varbetat <-   rowSums((DLambeta.t %*% iH)*DLambeta.t)
} else varbetat <- 0
wwJ <- opt$caseweightsJ*opt$weightsJ
var.cumhaz <- cumsumstrata(1/(opt$S0^2*wwJ),strata,nstrata)+varbetat
se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)
colnames(se.cumhaz) <- c("time","se.cumhaz")

out <- list(coef=beta.s,var=varmc,se.coef=diag(varmc)^.5,iid.naive=UUiid,
	iid=Uiid,ncluster=nid,
	ihessian=iH,hessian=opt$hessian,
	hessianttime=opt$hessianttime,var1=var1,se1.coef=diag(var1)^.5,
	ploglik=opt$ploglik,gradient=opt$gradient,
	cumhaz=cumhaz, se.cumhaz=se.cumhaz,MGciid=MGc,
	id=id.orig,call.id=call.id,
	strata.jumps=opt$strata.jumps,
	strata=xx2$strata,
	nstrata=nstrata,strata.name=strata.name,strata.level=strata.level,
	propodds=propodds,
	S0=opt$S0,E=opt$E,S2S0=opt$S2S0,time=opt$time,Ut=opt$U,
	jumps=jumps,exit=exit,p=p,S0s=val$S0s,
	no.opt=no.opt,##n=nrow(X),nevent=length(jumps),
	Pcens.model=Pcens.model,Gjumps=Gjumps,cens.code=cens.code,
	death.code=death.code, cause=cause, 
	strataA=strataA,nstrataA=nstrataA,augmentA=augment.new,MGAc=MGAc,
	augmentation=augmentation.call,
	var.augment=var.augment,var.augment.times=var.augment.times,
	var.augment.iid=var.augment.iid,var.augment.times.iid=var.augment.times.iid,
	lin.augment=c(augment),lindyn.augment=c(augment.times),
	iid.augment=Uiid.augment,iid.augment.times=Uiid.augment.times,
	gamma=gamma, gamma.times=gammat, time.gammat=time.gammat,ftime.gamma=ftime.gamma,Gcj=Gcj
	)

if (cox.prep) out <- c(out,list(cox.prep=xx2))

    return(out)
}# }}}

##' @export
recregIPCW <- function(formula,data=data,cause=1,cens.code=0,death.code=2,
	       cens.model=~1,km=TRUE,times=NULL,beta=NULL,offset=NULL,type=c("incIPCW"),
	       weights=NULL,model="exp",no.opt=FALSE,method="nr",augment.model=~+1,se=TRUE,...)
{# {{{
   ## type=c("incIPCW","IPCW","rate")
    cl <- match.call()# {{{
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
###    if (!is.null(intpos  <- attributes(Terms)$intercept))
###        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    ## }}}

   if (!is.null(id)) {
        ids <- unique(id)
        nid <- length(ids)
        if (is.numeric(id))
            id <-  fast.approx(ids,id)-1
        else  {
            id <- as.integer(factor(id,labels=seq(nid)))-1
        }
    } else { id <- as.integer(seq_along(entry))-1;  nid <- nrow(X); }
    ## orginal id coding into integers 1:...
    orig.id <- id.orig <- id+1;
    nid <- length(unique(id))

    ## to avoid Rcheck wawning

 ### setting up with artificial names
 data$status__ <-  status 
 data$id__ <-  id
 ## lave Countcause
 data <- count.history(data,status="status__",id="id__",types=cause,multitype=TRUE)
 data$Count1__ <- data[,paste("Count",cause[1],sep="")]
 data$death__ <- (status %in% death.code)*1
 data$entry__ <- entry 
 data$exit__ <- exit 
 statusC <- data$statusC__ <- (status %in% cens.code)*1
 data$status__cause <- (status %in% cause)*1
 data$rid__ <- revcumsumstrata(rep(1,length(entry)),id,nid)
 dexit <- exit
 dstatus <- status
 ## to define properly 
 Dtime <- NULL

  xr <- phreg(Surv(entry__,exit__,status__cause)~Count1__+death__+cluster(id__),data=data,no.opt=TRUE,no.var=1)
  formC <- update.formula(cens.model,Surv(entry__,exit__,statusC__)~ . +cluster(id__))
  cr <- phreg(formC,data=data)
  whereC <- which(status %in% cens.code)

###  xr <- phreg(Surv(entry__,exit__,status__ %in% cause)~cens__+cluster(id__),data=data,no.opt=TRUE,no.var=1)

  formC <- update.formula(cens.model,Surv(entry__,exit__,statusC__)~ .+cluster(id__))
  cr <- phreg(formC,data=data,no.opt=TRUE,no.var=1)
  whereC <- which(status %in% cens.code)

  if (length(whereC)>0) {# {{{
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
  } else  St <- rep(1,nrow(xx$strata))
  Gc <- St
  ## }}}

  formD <- as.formula(Surv(entry__,exit__,death__)~cluster(id__))
  form1L <- as.formula(Surv(entry__,exit__,status__cause)~Count1__+death__+statusC__+cluster(id__))
  form1 <- as.formula(Surv(entry__,exit__,status__cause)~cluster(id__))
  xr <- phreg(form1L,data=data,no.opt=TRUE,no.var=1)
###  dr <- phreg(formD,data=data,no.opt=TRUE,no.var=1)
###   ### cook-lawless ghosh-lin
###   xr0 <- phreg(form1,data=data,no.opt=TRUE)
###   clgl  <- recurrentMarginal(xr0,dr)
###   plot(clgl)
 
  ####  First \mu_ipcw(t) \sum_i I(T_i /\ t \leq C_i)/G_c(T_i /\ t ) N_(T_i /\ t) {{{
  x <- xr
  xx <- xr$cox.prep
  jump1 <- xx$jumps+1
  timeJ <- xx$time[jump1]
  strataN1J <- xx$strata[jump1]
 ### Partitioned estimator , same as Ghosh-Lin+Lawless-Cook estimator
 cumhazP <- c(cumsum(1/Gc[jump1])/nid)
 cumhazP <- cbind(timeJ,cumhazP)
# }}}

  if (is.null(times)) stop("time for recurrent events regression must be given\n")

  ### setting up regression setting with Y(t) =\int_0^t 1/G(s) dN_i(s)
 if (type[1]=="incIPCW") 
  Ydata <- Y <- sumstrata((xx$status!=0)*(xx$time<times)/Gc,xx$id,nid)
 else if (type[1]=="IPCW")  {
     obs <- (exit<=time & (!statusC)) | (exit>=time)/cens.weights
  Ydata <- Y <- sumstrata((xx$status!=0)*(xx$time<times),xx$id,nid)*obs
  } else {
     obs <- (exit<=time & (!statusC)) | (exit>=time)/cens.weights
  NtD <- sumstrata((xx$status!=0)*(xx$time<times),xx$id,nid)
  Ydata <- Y <- obs*NtD/pmin(times,Dtime)
 }

  ## back to ordering in data
  Ydata <- Y <- Y[id[data$rid__==1]+1]
  ###
  if (is.null(offset)) offset <- rep(0, length(exit))
  if (is.null(weights)) weights <- rep(1, length(exit))
  ###  
  Xorig <- X <- as.matrix(X)
  Xdata <- X <- X[data$rid__==1,,drop=FALSE]
  offset <- offset[data$rid__==1]
  weights <- weights[data$rid__==1]
  status <- status[data$rid__==1]
  exit <- exit[data$rid__==1]
  id <- id[data$rid__==1]
  X2 <- .Call("vecMatMat", X, X)$vXZ
  ph <- 1
  if (is.null(beta)) beta <- rep(0,ncol(X))
  ## take iid vession of data 
  dataiid <- data[data$rid__==1,]

  nevent <- sum(status %in% cause)

    obj <- function(pp, all = FALSE) {# {{{
        lp <- c(X %*% pp + offset)
        if (model == "exp") p <- exp(lp) else p <- lp
        if (model == "dexp") p <- exp(lp) 
        ploglik <- sum(weights * (Y - p)^2)
        if (model == "dexp") {
            Dlogl <- weights * X *p * c(Y - p)
            D2logl <- c(weights) *p^2* X2
        } else {
            Dlogl <- weights *  X * c(Y - p)
	if (model=="exp") D2logl <- c(weights) * c(p)* X2 else D2logl <- c(weights) * X2
	}
        D2log <- apply(D2logl, 2, sum)
        gradient <- apply(Dlogl, 2, sum) 
        hessian <- matrix(D2log, length(pp), length(pp))

        if (all) {
            ihess <- solve(hessian)
            beta.iid <- Dlogl %*% ihess
            beta.iid <- apply(beta.iid, 2, sumstrata, id, max(id) + 1)
            robvar <- crossprod(beta.iid)
            val <- list(par = pp, ploglik = ploglik, gradient = gradient,
                hessian = hessian, ihessian = ihess, id = id,
                Dlogl = Dlogl, iid = beta.iid, robvar = robvar,
                var = robvar, se.robust = diag(robvar)^0.5)
            return(val)
        }
        structure(-ploglik, gradient = -gradient, hessian = hessian)
    }# }}}

    p <- ncol(X)
    opt <- NULL
    if (p > 0) {
        if (no.opt == FALSE) {
            if (tolower(method) == "nr") {
                tim <- system.time(opt <- lava::NR(beta, obj, ...))
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

   val <- c(val, list(times = times, Y=Y, ncluster=nid, nevent=nevent, model.frame=m, n=length(exit),X=X))

    if (se & type[1]=="incIPCW") {# {{{

       Gcdata <- suppressWarnings(predict(cr,data,times=dexit,individual.time=TRUE,se=FALSE,km=km,tminus=TRUE)$surv)
       Gcdata[Gcdata<0.000001] <- 0.00001
       data$Hst <- revcumsumstrata((dexit<times)*(dstatus %in% cause)/Gcdata,data$id__,nid)
       if (model=="dexp") HstX <-c(exp(as.matrix(Xorig) %*% val$coef))*Xorig*c(data$Hst) else HstX <- Xorig*c(data$Hst) 
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
       E[xx$jumps + 1, ] <- resC$E
       btime <- c(1 * (xx$time < times))
       EdLam0 <- apply(E*c(S0i)*btime,2,cumsumstrata,xx$strata,xx$nstrata)
       MGt <- E[, drop = FALSE] - EdLam0 * c(xx$sign) 
       MGCiid <- apply(MGt, 2, sumstrata, xx$id, max(id) + 1)
    } else if (se & type[1]!="incIPCW") {

    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    cens.weights <- cens.weights[ord]
    lp <- c(X %*% val$coef+offset)
    p <- expit(lp)
    Y <- c((status %in% cause)*weights*(exit<=time)/cens.weights)

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply(h*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- h[xx$jumps+1,] /c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    MGCiid <- apply(MGt,2,sumstrata,xx$id,max(id)+1)
   }  else { 
	  MGCiid <- 0
  }## }}}

    if (se) val$MGciid <- MGCiid %*% val$ihessian else val$MGciid <- MGCiid 
    val$id <- id
    val$Y <- Ydata
    val$X <- Xdata
    val$orig.id <- orig.id
    val$iid.origid <- ids
    val$iid.naive <- val$iid
    if (se) val$iid <- val$iid + val$MGciid
    val$naive.var <- val$var
    robvar <- crossprod(val$iid)
    val$var <- val$robvar <- robvar
    val$se.robust <- diag(robvar)^0.5
    val$se.coef <- diag(val$var)^0.5
    val$cens.code <- cens.code
    val$cause <- cause
    val$death.code <- death.code
    val$model.type <- model
    val$cumhazP <- cumhazP
    class(val) <- c("binreg", "resmean")
    return(val)
} # }}}

##' @export
strataAugment <- survival:::strata

##' @export
simGLcox <- function(n,base1,drcumhaz,var.z=0,r1=NULL,rd=NULL,rc=NULL,fz=NULL,fdz=NULL,
   model=c("twostage","frailty","shared","multiplicative"),type=NULL,share=1,cens=NULL,nmin=100,nmax=1000)
{# {{{
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
 ## W_1 ~ N1, W_1+W_2 ~ D observed hazards on Cox form among survivors
 ## or W_1 ~ N1, W_1~ D   observed hazards on Cox form among survivors
 ## or W_2 * W_1 ~ N1, W_1~ D   observed hazards on Cox form among survivors
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

 return(ll)
}# }}}

simGLcoxC <- function(n,base1,drcumhaz,var.z=0,r1=NULL,rd=NULL,rc=NULL,fz=NULL,fdz=NULL,
	     model=c("twostage","frailty","shared","multiplicative"),type=NULL,share=1,cens=NULL,nmax=200,by=1)
{# {{{
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
}# }}}

simRecurrentCox <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,X=NULL,r1=NULL,r2=NULL,rd=NULL,rc=NULL, 
                     model=c("not-random","random"),frailty=TRUE,var.z=0.5,death.code=3,alpha=1,...)
{# {{{
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
 if (model[1]!="random") {# {{{
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
 }# }}}

 if (model[1]=="random") {# {{{
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
 }# }}}

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
}# }}}

simMarginalMeanCox <- function(n,cens=3/5000,k1=0.1,k2=0,bin=1,Lam1=NULL,Lam2=NULL,LamD=NULL,
			       beta1=rep(0,2),betad=rep(0,2),betac=rep(0,2),X=NULL,...)
{# {{{
###

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
}# }}}

##' @export
scalecumhaz <- function(cumt,k)
{# {{{
	return( t(t(cumt)*c(1,k)))
}# }}}

##' @export
GLprediid <- function(...)
{# {{{
out <- FGprediid(...,model="GL")
return(out)
}# }}}

boottwostageREC <- function(margsurv,recurrent,data,bootstrap=100,id="id",stepsize=0.5,...) 
{# {{{
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
      cause=recurrent$cause,death.code=recurrent$death.code,cens.code=recurrent$cens.code,cox.prep=TRUE)
    } else xrb <- phreg(recurrent$formula,data=rrbs)
     outbl <- tryCatch(twostageREC(drb,xrb,rrbs,numderiv=0,control=list(stepsize=stepsize),...),error=function(x) NULL)
     if (!is.null(outbl)) outb <- rbind(outb,outbl$coef)
     outd <- rbind(outd,coef(drb))
     outr <- rbind(outr,coef(xrb))
  }
  varb <- cov(outb)
  vard <- cov(outd)
  varr <- cov(outr)

  list(outb=outb,var=varb,se=diag(varb)^.5,se.coxD=diag(vard)^.5,se.coxR=diag(varr)^.5)
}# }}}

##' @export
twostageREC  <-  function (margsurv,recurrent, data = parent.frame(), theta = NULL, model=c("full","shared","non-shared"),ghosh.lin=NULL,
  theta.des = NULL, var.link = 0, method = "NR", no.opt = FALSE, weights = NULL, se.cluster = NULL, 
  fnu=NULL,nufix=0,nu=NULL,at.risk=1,numderiv=1,derivmethod=c("simple","Richardson"),...)
{# {{{
    if (!inherits(margsurv, "phreg")) stop("Must use phreg for death model\n")
    if (!inherits(recurrent, "phreg")) stop("Must use phreg for recurrent model\n")
    if (is.null(recurrent$cox.prep)) stop("recreg must be called with cox.prep=TRUE\n")
    if (inherits(recurrent, "recreg") & is.null(ghosh.lin))  ghosh.lin <- 1 
    if ((!inherits(recurrent, "recreg")) & is.null(ghosh.lin))  ghosh.lin <- 0 
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

    obj <- function(par, all = FALSE) {# {{{
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
# }}}

    objNonShared <- function(par, all = FALSE) {# {{{
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
# }}}

    ## default is simple identity 
###    if (is.null(fnu)) { fw <- function(x) 1/(1+exp(x)); Dfw <- function(x) -exp(x)/(1+exp(x))^2;} else { 
   if (is.null(fnu)) { fw <- function(x) x; Dfw <- function(x) 1;} else { fw <- fnu[[1]]; Dfw <- fnu[[2]]; } 
   nudes <- theta.des  
   p <- ncol(theta.des)

    objShared <- function(par, all = FALSE) {# {{{
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
# }}}

   if (model[1]=="shared") obj <- objShared
   if (model[1]=="non-shared") obj <- objNonShared
   if (nufix==0 & model[1]=="shared") par <- c(theta,nu) else par <- theta

    opt <- NULL
    if (no.opt == FALSE) {
        if (tolower(method) == "nr") {
            opt <- lava::NR(par, obj, ...)
            opt$estimate <- opt$par
        }
        else {
            opt <- nlm(obj, par, ...)
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
}
# }}}

##' @export
summary.twostageREC <- function(object,vcov=NULL,delta=0,...) {# {{{
    I <- -solve(object$hessian)
    if (!is.null(vcov)) V <- vcov else V <- object$var
    ncluster <- object$n
    cc <- estimate(coef=object$coef,vcov=V)$coefmat
    pd <- object$p
    if (object$var.link==1 & object$model=="full") f <- function(p) exp(p)
    if (object$var.link==0 & object$model=="full") f <- function(p) p
    if (object$var.link==1 & object$model=="shared") f <- function(p) c(exp(p[1]),attr(object,"fnu")(p[2]))
    if (object$var.link==0 & object$model=="shared") f <- function(p) c(p[1],attr(object,"fnu")(p[2]))
    if (delta==1 | object$model=="full") 
    expC <- lava::estimate(coef=object$coef,vcov=V,f=f)$coefmat ##[,c(1,3,4),drop=FALSE]
    else expC <- apply(cc[,c(1,3,4),drop=FALSE],2,f) 
  n <- object$n
  res <- list(coef=cc,n=n,nevent=object$nevent,ncluster=ncluster,var=V,exp.coef=expC,var.link=object$var.link,
	      ghosh.lin=object$ghosh.lin)
  class(res) <- "summary.twostageREC"
  res
} # }}}

##' @export
print.summary.twostageREC  <- function(x,max.strata=5,...) {# {{{
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
} # }}}

##' @export
print.twostageREC  <- function(x,...) {# {{{
  print(summary(x),...)
} # }}}

##' Evaluates piece constant covariates at min(D,t) where D is a terminal event
##'
##' returns X(min(D,t)) and min(D,t) and their ratio. for censored observation 0. 
##' to use with the IPCW models implemented. 
##'
##' @param formula formula with 'Event' outcome and X to evaluate at min(D,t)
##' @param data data frame
##' @param death.code codes for death (terminating event, 2 default)
##' @param time for evaluation 
##' @author Thomas Scheike
##' @export
evalTerminal <- function(formula,data=data,death.code=2,time=NULL)
{# {{{
    cl <- match.call()# {{{
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
        pos.cluster <- ts$terms
        Terms  <- Terms[-ts$terms]
        id <- m[[ts$vars]]
    } else pos.cluster <- NULL
###    if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
###        ts <- survival::untangle.specials(Terms, "strata")
###        pos.strata <- ts$terms
###        Terms  <- Terms[-ts$terms]
###        strata <- m[[ts$vars]]
###        strata.name <- ts$vars
###    }  else { strata.name <- NULL; pos.strata <- NULL}
###    if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
###        ts <- survival::untangle.specials(Terms, "offset")
###        Terms  <- Terms[-ts$terms]
###        offset <- m[[ts$vars]]
###    }
    X <- model.matrix(Terms, m)
    if (!is.null(intpos  <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    ## }}}

   if (!is.null(id)) {
	   call.id <- id
        ids <- unique(id)
        nid <- length(ids)
        if (is.numeric(id))
            id <-  fast.approx(ids,id)-1
        else  {
            id <- as.integer(factor(id,labels=seq(nid)))-1
        }
    } else { call.id <- id <- as.integer(seq_along(entry))-1;  nid <- nrow(X); }
    ## orginal id coding into integers 1:...
   id <- id+1

   dd <- data.frame(id=id)
   dd <- countID(dd,sorted=TRUE)
   ## new id 1,2,.... and so on, referring to rows of data
   id <- dd$indexid+1

 ###
 indexD  <- which(exit <= time & (status %in% death.code))
 ### rr[indexD,c("time1","statusD")]
 idD <- id[indexD]
 Dmintid <- rep(0,nid)
 Dmintid[idD] <- exit[indexD]
 ### 
 indexA<- which(entry <= time & time <= exit)
 idA <- id[indexA]
 Dmintid[idA] <- time
 Dmint <- Dmintid[id]
 ###
 XminDtid <- matrix(0,nid,ncol(X))
 colnames(XminDtid) <- colnames(X)
 XminDtid[idD,] <- X[indexD,]
 XminDtid[idA,]  <- X[indexA,]
 XminDt <- XminDtid[id,,drop=FALSE]
 ratio <- XminDt/Dmint
 ratio[is.na(ratio)] <- 0

 nX <- colnames(XminDt) <- paste(colnames(X),"minDt",sep="") 
 dd <- data.frame(cbind(XminDt,Dmint,id,call.id,ratio))
 colnames(dd) <- c(nX,"minDt","nid","call.id","ratio")

 return(dd)
} # }}}


dynCensAugOld <- function(formC,data,augmentC=~+1,response="Yipcw",time=NULL,Z=NULL) {# {{{ 
if (is.null(time)) stop("must give time of response \n")
   data$Y__ <- data[,response]
   varsC <- c("Y__",attr(terms(augmentC), "term.labels"))
   formCC <- update(formC, reformulate(c(".", varsC)))
###	www <- rep(1,nrow(data))
###	if (treat.specific.cens==1)  www <-data$WW1__;  
   cr2 <- phreg(formCC, data = data, no.opt = TRUE, no.var = 1,Z=Z)
   xx <- cr2$cox.prep
   icoxS0 <- rep(0,length(cr2$S0))
   icoxS0[cr2$S0>0] <- 1/cr2$S0[cr2$S0>0]
   S0i <- rep(0, length(xx$strata))
   S0i[xx$jumps + 1] <- icoxS0
   km <- TRUE
   if (!km) {
        cumhazD <- c(cumsumstratasum(S0i, xx$strata, xx$nstrata)$lagsum)
        St <- exp(-cumhazD)
   } else St <- c(exp(cumsumstratasum(log(1 - S0i), xx$strata, xx$nstrata)$lagsum))

   nterms <- cr2$p-1
   dhessian <- cr2$hessianttime
   dhessian <-  .Call("XXMatFULL",dhessian,cr2$p,PACKAGE="mets")$XXf
   ###  matrix(apply(dhessian,2,sum),3,3)
   timeb <- which(cr2$cumhaz[,1]<time)
   ### take relevant \sum H_i(s,t) (e_i - \bar e)
   covts <- dhessian[timeb,1+1:nterms,drop=FALSE]
   ### construct relevant \sum (e_i - \bar e)^2
   Pt <- dhessian[timeb,-c((1:(nterms+1)),(1:(nterms))*(nterms+1)+1),drop=FALSE]
   ###  matrix(apply(dhessian[,c(5,6,8,9)],2,sum),2,2)
   gammatt <- .Call("CubeVec",Pt,covts,1,PACKAGE="mets")$XXbeta
   S0 <- cr2$S0[timeb]
   gammatt[is.na(gammatt)] <- 0
   gammatt[gammatt==Inf] <- 0
   Gctb <- St[cr2$cox.prep$jumps+1][timeb]
   augmentt.times <- apply(gammatt*cr2$U[timeb,1+1:nterms,drop=FALSE],1,sum)
   augment.times <- sum(augmentt.times)
   if (!is.null(Z)) {
             Zj <- cr2$cox.prep$Z[cr2$cox.prep$jumps+1,][timeb]
             Xaugment.times <- apply( augmentt.times*Zj,2,sum)
   }

   p <- 1
   #### iid magic  for censoring augmentation martingale ## {{{
   ### int_0^infty gamma(t) (e_i - ebar(s)) 1/G_c(s) dM_i^c
   xx <- cr2$cox.prep
   nid <- max(xx$id)+1
   jumpsC <- xx$jumps+1
   rr0 <- xx$sign
   S0i <- rep(0,length(xx$strata))
   S0i[jumpsC] <- c(1/(icoxS0*St[jumpsC]))
   S0i[jumpsC] <- icoxS0

   pXXA <- ncol(cr2$E)-1
   EA <- cr2$E[timeb,-1,drop=FALSE]
   gammasEs <- .Call("CubeMattime",gammatt,EA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
   gammasE <- matrix(0,length(xx$strata),1)
   gammattt  <-    matrix(0,length(xx$strata),pXXA*1)
   jumpsCt <- jumpsC[timeb]
   gammasE[jumpsCt,] <- gammasEs
   gammattt[jumpsCt,] <- gammatt
   gammaEsdLam0 <- apply(gammasE*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
   gammadLam0 <-   apply(gammattt*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
   XgammadLam0 <- .Call("CubeMattime",gammadLam0,xx$X[,-1],pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
   Ut <- Et <- matrix(0,length(xx$strata),1)
   Ut[jumpsCt,] <- augmentt.times
   MGCt <- Ut[,drop=FALSE]-(XgammadLam0-gammaEsdLam0)*c(rr0)
   MGCiid <- apply(MGCt,2,sumstrata,xx$id,nid)
   MGCiid <- MGCiid/nid
   ## }}}

   nid <- max(cr2$id)
   ids <- headstrata(cr2$id-1,nid)
   ids <- cr2$call.id[ids]

   res <- list(MGCiid=MGCiid,gammat=gammatt,augment=augment.times, id=ids,n=nid)
} # }}}
