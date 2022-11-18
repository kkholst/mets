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
##' @param formula formula with 'EventCens' outcome
##' @param data data frame
##' @param cause of interest
##' @param death.code codes for death (terminating event)
##' @param cens.code code of censoring (1 default)
##' @param cens.model for stratified Cox model without covariates
##' @param weights weights for score equations
##' @param offset offsets for model
##' @param Gc censoring weights for time argument, default is to calculate these with a Kaplan-Meier estimator, should then give G_c(T_i-)
##' @param wcomp weights for composite outcome, so when cause=c(1,3), we might have wcomp=c(1,2).
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ## data with no ties
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' Lam1 <- base1cumhaz;  Lam2 <- base4cumhaz;  LamD <- drcumhaz
##' ## simulates recurrent events of types 1 and 2 and with terminal event D and censoring
##' rr <- simRecurrentII(1000,Lam1,cumhaz2=Lam2,death.cumhaz=LamD,cens=3/5000)
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
##' @aliases strataAugment scalecumhaz GLprediid recregIPCW
##' @export
recreg <- function(formula,data=data,cause=1,death.code=c(2),cens.code=0,cens.model=~1,weights=NULL,offset=NULL,Gc=NULL,
		   wcomp=NULL,...)
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

    res <- c(recreg01(data,X,entry,exit,status,id=id,strata=strata,offset=offset,weights=weights,
		      cens.model=cens.model, cause=cause, strata.name=strata.name, strataA=NULL,## strataAugment,
		      death.code=death.code,cens.code=cens.code,Gc=Gc,wcomp=wcomp,...),
             list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster,n=nrow(X),nevent=sum(status %in%cause))
             )

    class(res) <- c("phreg","recreg")
    return(res)
}# }}}

recreg01 <- function(data,X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,strataA=NULL,
          strata.name=NULL,beta,stderr=1,method="NR",no.opt=FALSE, propodds=NULL,profile=0,
          case.weights=NULL,cause=1,death.code=2,cens.code=0,Gc=NULL,cens.model=~+1,augmentation=0,
	  cox.prep=FALSE,wcomp=NULL,...) { # {{{
# {{{ setting up weights, strata, beta and so forth before the action starts
    p <- ncol(X)
    if (missing(beta)) beta <- rep(0,p)
    if (p==0) X <- cbind(rep(0,length(exit)))

    cause.jumps <- which(status %in% cause)
    max.jump <- max(exit[cause.jumps])
    other <- which((status %in% death.code ) )
###    other <- which((status %in% death.code ) & (exit< max.jump))

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
    whereC <- which( status %in% cens.code)
    time <- exit
    cens <- statusC <- c(status %in% cens.code)
    data$id <- id
    data$exit__ <- exit
    data$entry__ <- entry
    data$statusC <- statusC
    cens.strata <- cens.nstrata <- NULL

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
    Zcall <- cbind(status,cens.strata,Stime,cens,strata,strataA) ## to keep track of status and Censoring strata
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
        timeoo <- rep(max(exit)+1,length(other))
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
        ploglik <- sum(ploglik)
        U <- U+augmentation

        out <- list(ploglik=ploglik,gradient=U,hessian=-DU,cox.prep=xx2,
                    hessiantime=DUt,weightsJ=weightsJ,caseweightsJ=caseweightsJ,
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

        ### xx2 data all data with start stop structure, takes position of death times
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
        ## }}}
    } else MGc <- 0

    if (!is.null(orig.strataA)) { ### compute augmentation term based on this beta
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


    iH <- - tryCatch(solve(opt$hessian),error=function(e) matrix(0,nrow(opt$hessian),ncol(opt$hessian)) )
    Uiid <-  (UU+MGc) %*% iH
    UUiid <- UU %*% iH
    var1 <-  crossprod(UUiid)
    varmc <-  crossprod(Uiid)
    ### end if (p>0)
    } else {varmc <- var1 <- 0; augment.new <- MGAc <- MGc <- iH <- UUiid <- Uiid <- NULL}
    strata <- xx2$strata[jumps]
    cumhaz <- cbind(opt$time,cumsumstrata(1/opt$S0,strata,nstrata))
    colnames(cumhaz)    <- c("time","cumhaz")

## SE of estimator ignoring some censoring terms
if (no.opt==FALSE & p!=0) {
DLambeta.t <- apply(opt$E/c(opt$S0),2,cumsumstrata,strata,nstrata)
varbetat <-   rowSums((DLambeta.t %*% iH)*DLambeta.t)
} else varbetat <- 0
var.cumhaz <- cumsumstrata(1/opt$S0^2,strata,nstrata)+varbetat
se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)
colnames(se.cumhaz) <- c("time","se.cumhaz")


out <- list(coef=beta.s,var=varmc,se.coef=diag(varmc)^.5,iid.naive=UUiid,
	iid=Uiid,ncluster=nid,
	ihessian=iH,hessian=opt$hessian,var1=var1,se1.coef=diag(var1)^.5,
	ploglik=opt$ploglik,gradient=opt$gradient,
	cumhaz=cumhaz, se.cumhaz=se.cumhaz,MGciid=MGc,
	strata=xx2$strata,
	nstrata=nstrata,strata.name=strata.name,strata.level=strata.level,
	propodds=propodds,
	S0=opt$S0,E=opt$E,S2S0=opt$S2S0,time=opt$time,Ut=opt$U,
	jumps=jumps,exit=exit,p=p,S0s=val$S0s,
	no.opt=no.opt,##n=nrow(X),nevent=length(jumps),
	Pcens.model=Pcens.model,Gjumps=Gjumps,cens.code=cens.code,
	death.code=death.code, cause=cause, 
	strataA=strataA,nstrataA=nstrataA,augment=augment.new,MGAc=MGAc
	)

if (cox.prep) out <- c(out,list(cox.prep=xx2))

    return(out)
}# }}}

##' @export
recregIPCW <- function(formula,data=data,cause=1,cens.code=0,death.code=2,
	       cens.model=~1,km=TRUE,times=NULL,beta=NULL,offset=NULL,
	       weights=NULL,model="exp",no.opt=FALSE,method="nr",augment.model=~+1,se=TRUE,...)
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

 ### setting up with artificial names
 data$status__ <-  status 
 data$id__ <-  id
 ## lave Countcause
 data <- count.history(data,status="status__",id="id__",types=cause,multitype=TRUE)
 data$Count1__ <- data[,paste("Count",cause[1],sep="")]
 data$death__ <- (status %in% death.code)*1
 data$entry__ <- entry 
 data$exit__ <- exit 
 data$statusC__ <- (status %in% cens.code)*1
 data$status__cause <- (status %in% cause)*1
 data$rid__ <- revcumsumstrata(rep(1,length(entry)),id,nid)
 dexit <- exit
 dstatus <- status

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
 ### Partitioned estimator , same as Lin, Lawless & Cook estimator
 cumhazP <- c(cumsum(1/Gc[jump1])/nid)
 cumhazP <- cbind(timeJ,cumhazP)
# }}}

  if (is.null(times)) stop("time for recurrent events regression must be given\n")

  ### setting up regression setting with Y(t) =\int_0^t 1/G(s) dN_i(s)
  Ydata <- Y <- cumhazPiid <- sumstrata((xx$status!=0)*(xx$time<times)/Gc,xx$id,nid)
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
        ploglik <- sum(weights * (Y - p)^2)
        if (model == "exp") {
            Dlogl <- weights *  X * c(Y - p)
            D2logl <- c(weights * p) * X2
        }
        else {
            Dlogl <- weights * X * c(Y - p)
            D2logl <- c(weights) * X2
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
###	     formula = formula, formC = formC,
###        exit = exit, cens.weights = cens.weights, cens.strata = cens.strata,
###        cens.nstrata = cens.nstrata, model.frame = m, n = length(exit),
###        nevent = nevent, ncluster = nid, Y = Y))

    if (se) {# {{{

       Gcdata <- suppressWarnings(predict(cr,data,times=dexit,individual.time=TRUE,se=FALSE,km=km,tminus=TRUE)$surv)
       Gcdata[Gcdata<0.000001] <- 0.00001
       data$Hst <- revcumsumstrata((dexit<times)*(dstatus %in% cause)/Gcdata,data$id__,nid)
       HstX <- Xorig*c(data$Hst)
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
    }
    else {
        MGCiid <- 0
    }# }}}


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

simRecurrentCox <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,X=NULL,r1=NULL,r2=NULL,rd=NULL,rc=NULL, 
    model=c("not-random","random"),fraily=TRUE,var.z=0.5,death.code=3,alpha=1,...)
{# {{{
  if (is.null(r1)) r1 <- rep(1,n)
  if (is.null(r2)) r2 <- rep(1,n)
  if (is.null(rd)) rd <- rep(1,n)
  if (is.null(rc)) rc <- rep(1,n)

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


## Fast recurrent marginal mean when death is possible
##'
##' Fast Marginal means of recurrent events. Using the Lin and Ghosh (2000) standard errors.  
##' Fitting two models for death and recurent events these are
##' combined to prducte the estimator 
##' \deqn{ \int_0^t  S(u|x=0) dR(u|x=0) } the mean number of recurrent events, here
##' \deqn{ S(u|x=0) }  is the probability of survival for the baseline group, and 
##' \deqn{ dR(u|x=0) }  is the hazard rate of an event among survivors for the baseline. 
##' Here \deqn{ S(u|x=0) }  is estimated by \deqn{ exp(-\Lambda_d(u|x=0) }  with 
##'  \deqn{\Lambda_d(u|x=0) } being the cumulative baseline for death.
##' 
##' Assumes no ties in the sense that jump times needs to be unique, this is particularly so for the stratified version.
##' 
##' @param recurrent phreg object with recurrent events
##' @param death     phreg object with deaths
##' @param fixbeta   to force the estimation of standard errors to think of regression coefficients as known/fixed
##' @param km  if true then uses Kaplan-Meier for death, otherwise exp(- Nelson-Aalen ) 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' 
##' @references 
##'             Ghosh and Lin (2002) Nonparametric Analysis of Recurrent events and death, 
##'             Biometrics, 554--562.
##' @examples
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##' rr <- simRecurrent(1000,base1,death.cumhaz=dr)
##' rr$x <- rnorm(nrow(rr)) 
##' rr$strata <- floor((rr$id-0.01)/500)
##' 
##' ##  to fit non-parametric models with just a baseline 
##' xr <- phreg(Surv(entry,time,status)~cluster(id),data=rr)
##' dr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)
##' par(mfrow=c(1,3))
##' bplot(dr,se=TRUE)
##' title(main="death")
##' bplot(xr,se=TRUE)
##' ### robust standard errors 
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=4)
##' 
##' ## marginal mean of expected number of recurrent events 
##' out <- recurrentMarginal(xr,dr)
##' bplot(out,se=TRUE,ylab="marginal mean",col=2)
##' 
##' ########################################################################
##' ###   with strata     ##################################################
##' ########################################################################
##' xr <- phreg(Surv(entry,time,status)~strata(strata)+cluster(id),data=rr)
##' dr <- phreg(Surv(entry,time,death)~strata(strata)+cluster(id),data=rr)
##' par(mfrow=c(1,3))
##' bplot(dr,se=TRUE)
##' title(main="death")
##' bplot(xr,se=TRUE)
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)
##'
##' out <- recurrentMarginal(xr,dr)
##' bplot(out,se=TRUE,ylab="marginal mean",col=1:2)
##'
##' ########################################################################
##' ###   cox case        ##################################################
##' ########################################################################
##' xr <- phreg(Surv(entry,time,status)~x+cluster(id),data=rr)
##' dr <- phreg(Surv(entry,time,death)~x+cluster(id),data=rr)
##' par(mfrow=c(1,3))
##' bplot(dr,se=TRUE)
##' title(main="death")
##' bplot(xr,se=TRUE)
##' rxr <-   robust.phreg(xr)
##' bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)
##'
##' out <- recurrentMarginal(xr,dr)
##' bplot(out,se=TRUE,ylab="marginal mean",col=1:2)
##' 
##' ########################################################################
##' ###   CIF  #############################################################
##' ########################################################################
##' ### use of function to compute cumulative incidence (cif) with robust standard errors
##'  data(bmt)
##'  bmt$id <- 1:nrow(bmt)
##'  xr  <- phreg(Surv(time,cause==1)~cluster(id),data=bmt)
##'  dr  <- phreg(Surv(time,cause!=0)~cluster(id),data=bmt)
##' 
##'  out <- recurrentMarginal(xr,dr,km=TRUE)
##'  bplot(out,se=TRUE,ylab="cumulative incidence")
##' 
##' @aliases tie.breaker recmarg recurrentMarginalIPCW recurrentMarginalAIPCW recurrentMarginalAIPCWdata 
##' @export
recurrentMarginal <- function(recurrent,death,fixbeta=NULL,km=TRUE,...)
{# {{{
  xr <- recurrent
  dr <- death 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if ((xr$no.opt) | is.null(xr$coef)) fixbeta<- 1 else fixbeta <- 0

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## survival at t- to also work in competing risks situation
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
    St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazR <-  cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
  mu <- cumhazDR[,2]
# }}}

  ### robust standard errors 
  ### 1. sum_k ( int_0^t S(s)/S_0^r(s) dM_k.^r(s) )^2
 resIM1 <-  squareintHdM(xr,ft=St,fixbeta=fixbeta)
 ### 2. mu(t)^2 * sum_k ( int_0^t 1/S_0^d(s) dM_k.^d(s) )^2
 resIM2 <-  squareintHdM(dr,ft=NULL,fixbeta=fixbeta)
  ### 3. sum_k( int_0^t mu(s) /S_0^d(s) dM_k.^d(s))^2
 resIM3 <-  squareintHdM(dr,ft=mu,fixbeta=fixbeta)

 varA <-  resIM1$varInt+mu^2*resIM2$varInt+resIM3$varInt 

## covariances between different terms  13 23  12 12
## to allow different strata for xr and dr, but still nested strata
 if ((xr$nstrata>1 & dr$nstrata==1)) {
    cM1M3 <- covIntH1dM1IntH2dM2(resIM1,resIM3,fixbeta=fixbeta,mu=NULL)
    cM1M2 <- covIntH1dM1IntH2dM2(resIM1,resIM2,fixbeta=fixbeta,mu=mu)
 } else  {
    cM1M3 <- covIntH1dM1IntH2dM2(resIM3,resIM1,fixbeta=fixbeta,mu=NULL)
    cM1M2 <- covIntH1dM1IntH2dM2(resIM2,resIM1,fixbeta=fixbeta,mu=mu)
 }
 cM2M3 <- covIntH1dM1IntH2dM2(resIM2,resIM3,fixbeta=fixbeta,mu=mu)

 varA <- varA+2*cM1M3$cov12A-2*cM1M2$cov12A-2*cM2M3$cov12A 
### varA <- varA-2*cM1M3$cov12A+2*cM1M2$cov12A+2*cM2M3$cov12A 

 cov12aa <- cov13aa <- cov23aa <- 0

 if (fixbeta==0) {
    varA <-varA + cM2M3$covbeta - cM1M3$covbeta + cM1M2$covbeta 
 }

 varrs <- data.frame(mu=mu,cumhaz=mu,se.mu=varA^.5,time=xr$time,
		     se.cumhaz=varA^.5,strata=xr$strata,St=St)
 varrs <- varrs[c(xr$cox.prep$jumps)+1,]

 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,se.mu=varrs$se.mu,times=varrs$time,
     St=varrs$St,
     cumhaz=cbind(varrs$time,varrs$mu),se.cumhaz=cbind(varrs$time,varrs$se.mu),
     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
     strata.name=xr$strata.name,strata.level=recurrent$strata.level)
 class(out) <- "recurrent"
 return(out)
}# }}}

##' @export
summary.recurrent <- function(object,times=NULL,strata=NULL,estimates=FALSE,...) {# {{{
 if (is.null(times)) times <- object$times

if (object$nstrata==1) {
   where <- fast.approx(c(0,object$times),times,type="left")
   mu <- c(0,object$mu)[where]
   se.mu <- c(0,object$se.mu)[where]
   stratao <- 0
} else {
 nstrata <- object$nstrata
 if (is.null(strata))  {
    where <- indexstratarightR(object$times,object$strata,
    rep(times,each=nstrata),rep((nstrata-1):0,length(times)),nstrata,type="left")

 times <- rep(times,each=nstrata)
 strata <- rep((nstrata-1):0,length(times))

 } else where <- indexstratarightR(object$times,object$strata,times,strata,nstrata,type="left")

   mu <- object$mu[where]
   se.mu <- object$se.mu[where]
   stratao <- object$strata[where] 
}

 se.logmu=se.mu/mu
 lower <- exp(log(mu) - 1.96*se.logmu)
 upper <- exp(log(mu) + 1.96*se.logmu)

 out <- data.frame(times=times,mu=mu,se.mu=se.mu,lower=lower,upper=upper,
 strata=stratao)
 names(out) <- c("times","mean","se-mean","CI-2.5%","CI-97.5%","strata")
 if (estimates) {
         attr(out,"where") <- where
	 attr(out,"estimates") <- cbind(object$cumhaz,object$strata)[where,]
 }
 return(out)
}# }}}

##' @export
summaryTimeobject <-function(mutimes,mu,se.mu=NULL,times=NULL,type="log",...) {# {{{
 if (is.null(times)) times <- mutimes

 where <- fast.approx(c(0,mutimes),times,type="left")

 ##  see if object is vector or matrix
 if (is.matrix(mu)) mu <- rbind(0,mu)[where,] else mu <- c(0,mu)[where]
 if (!is.null(se.mu)) {
     if (is.matrix(se.mu)) se.mu <- rbind(0,se.mu)[where,] else se.mu <- c(0,se.mu)[where]
 se.logmu=se.mu/mu
 if (type=="log") {
 lower <- exp(log(mu) - 1.96*se.logmu)
 upper <- exp(log(mu) + 1.96*se.logmu)
 } else {
 lower <- mu - 1.96*se.mu
 upper <- mu + 1.96*se.mu
 }
 } else {se.mu <- se.logmu <- lower <- upper <- NA}


 out <- data.frame(times=times,mu=mu,se.mu=se.mu,lower=lower,upper=upper)
 names(out) <- c("times","mean","se-mean","CI-2.5%","CI-97.5%")
 return(out)
}# }}}

##' @export
recurrentMarginalAIPCW <- function(formula,data=data,cause=1,cens.code=0,death.code=2,
	  cens.model=~1,km=TRUE,times=NULL,augment.model=~Nt,...)
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
    id.orig <- id+1;

   ### setting up with artificial names
 data$status__ <-  status 
 data$id__ <-  id
 ## lave Countcause
 data <- count.history(data,status="status__",id="id__",types=cause,multitype=TRUE)
 data$Count1__ <- data[,paste("Count",cause[1],sep="")]
 data$death__ <- (status %in% death.code)*1
 data$entry__ <- entry 
 data$exit__ <- exit 
 data$statusC__ <- (status %in% cens.code)*1
 data$status__cause <- (status %in% cause)*1

  xr <- phreg(Surv(entry__,exit__,status__cause)~Count1__+death__+cluster(id__),data=data,no.opt=TRUE,no.var=1)

  formC <- update.formula(cens.model,Surv(entry__,exit__,statusC__)~ . +cluster(id__))
  cr <- phreg(formC,data=data)
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
 dr <- phreg(formD,data=data,no.opt=TRUE,no.var=1)

 ### augmenting partioned estimator computing \hat H_i(s,t) for fixed t
 data$Gctrr <- exp(-cpred(cr$cumhaz,exit)[,2])

 ### cook-lawless ghosh-lin

 xr0 <- phreg(form1,data=data,no.opt=TRUE)
 clgl  <- recurrentMarginal(xr0,dr)
### bplot(clgl,se=1); print(cpred(clgl$cumhaz,times)); print(cpred(clgl$se.cumhaz,times)); 

  ####  First \mu_ipcw(t) \sum_i I(T_i /\ t \leq C_i)/G_c(T_i /\ t ) N_(T_i /\ t) {{{
  x <- xr
  xx <- xr$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## stay with N(D_i) when t is large so no -1 when death
  ## xx$X[,1] er Count1 dvs N(t-)
  Nt <- revcumsumstrata(xx$X[,1]*xx$sign,xx$strata,xx$nstrata)
  Nt <- Nt/Gc
  ## counting N(D) forward in time skal ikke checke ud når man dør N_(D_i) er i spil efter D_i
  NtD <- cumsumstrata(xx$X[,1]*(xx$X[,2]==1)*(xx$sign==1)/Gc,xx$strata,xx$nstrata)
  jump1 <- xx$jumps+1
  timeJ <- xx$time[jump1]
  avNtD <- (NtD+Nt)[jump1]/nid
  strataN1J <- xx$strata[jump1]

  varIPCW1 <- NULL
  seIPCW1 <- NULL

  ### IPCW estimator 
  cumhaz <- cbind(timeJ,avNtD)
# }}}

  ### Partitioned estimator , same as Lin, Lawless & Cook estimator {{{
  cumhazP <- c(cumsumstrata(1/Gc[jump1],strataN1J,xx$nstrata)/nid)
  cumhazP <- cbind(timeJ,cumhazP)

  ### variance of partitioned estimator 

  ### calculate E(H,s,t) = E(H,t) - E(H,s) 
  ### E(H,t) = 1/(S(t)*n) \sum_ \int Y_i(s)/G_c(s) dN_{1i} = 1/(S(t)) (\mu_p(t) - \mu_p(s)) 
  ### \hat H_i for all subjects, and look at together with Hst, eHst 
  ### for censoring martingale 
  data$Nt <- data$Count1__
  data$Nt2 <- data$Nt^2
  data$expNt <- exp(-data$Nt)
  data$NtexpNt <- data$Nt*exp(-data$Nt)

  Gcdata <- exp(-cpred(cr$cumhaz,exit)[,2])

  form <- as.formula(paste("Surv(entry__,exit__,statusC__)~+1"))
  desform <- update.formula(augment.model,~Hst + . + cluster(id__))
  form[[3]] <- desform[[2]]
  nterms <- length(all.vars(form[[3]]))-1

  if (!is.null(times)) {
       semuPA <-  muPA <- semuPA.times <-  muPA.times <- rep(0,length(times))
       ww <- fast.approx(timeJ,times,type="left")
       muP.times <- cumhazP[ww,2]
       semuP.times <- clgl$se.mu[ww]

  for (i in seq_along(times)) {
     timel <- times[i]
     data$Hst <- revcumsumstrata((exit<timel)*(status %in% cause)/Gcdata,id,nid)
     cr2 <- phreg(form,data=data,no.opt=TRUE,no.var=1)
     nterms <- cr2$p-1

     dhessian <- cr2$hessianttime
     ###  matrix(apply(dhessian,2,sum),3,3)
     timeb <- which(cr$cumhaz[,1]<timel)
     ### take relevant \sum H_i(s,t) (e_i - \bar e)
     covts <- dhessian[timeb,1+1:nterms,drop=FALSE]
     ### construct relevant \sum (e_i - \bar e)^2
     Pt <- dhessian[timeb,-c((1:(nterms+1)),(1:(nterms))*(nterms+1)+1),drop=FALSE]
     ###  matrix(apply(dhessian[,c(5,6,8,9)],2,sum),2,2)
     gammahat <- .Call("CubeVec",Pt,covts,1,PACKAGE="mets")$XXbeta
     S0 <- cr$S0[timeb]
     gammahat[is.na(gammahat)] <- 0
     gammahat[gammahat==Inf] <- 0
     Gctb <- Gc[cr$jumps][timeb]
     augment.times <- sum(apply(gammahat*cr2$U[timeb,1+1:nterms,drop=FALSE],1,sum))/nid
     mterms <- length(terms)
     mterms <- nterms
     ###
     varZ <- matrix(apply(Pt/Gctb^2,2,sum),mterms,mterms)
     gamma <- .Call("CubeVec",matrix(c(varZ),nrow=1),matrix(apply(covts/Gctb,2,sum),nrow=1),1,PACKAGE="mets")$XXbeta
     gamma <- c(gamma)
     gamma[is.na(gamma)] <- 0
     gamma[gamma=Inf] <- 0
     augment <- sum(apply(gamma*t(cr2$U[timeb,1+1:nterms,drop=FALSE])/Gctb,2,sum))/nid
     ###
     muPA[i] <- muP.times[i]+augment
     semuPA[i] <- (semuP.times[i]^2 +(gamma %*% varZ %*% gamma)/nid^2)^.5
     muPA.times[i] <- muP.times[i]+augment.times
     semuPA.times[i] <- (semuP.times[i]^2+sum(gammahat*.Call("CubeVec",Pt,gammahat,0,PACKAGE="mets")$XXbeta)/(nid^2))^.5
  }
  }
# }}}

  return(list(censoring.weights=Gctb,muP.all=cumhazP,Gcjump=Gc[jump1],gamma=gamma,gamma.time=gammahat,times=times,
  muP=muP.times,semuP=semuP.times, muPAt=muPA.times,semuPAt=semuPA.times, muPA=muPA,semuPA=semuPA))
}# }}}

##' @export
recurrentMarginalAIPCWdata <- function(rr,times,km=TRUE,terms=1,idt=1,x.design=NULL,
   id="id",start="start",stop="stop",status="status",death="death",cause=1,...)
{# {{{

 if (missing(times)) stop("times of estimation must be given")

 # to avoid R check warning 
 revnr <- NULL

 formsort <- paste("~",id,"+",start)
 dsort(rr) <- as.formula(formsort) 
 rr$revnr <- NULL
 rr$cens <- 0
 rr <- count.history(rr,status=status,id=id)

 nid <- max(rr[,id])
 rr$revnr2 <-  c(revcumsumstrata(rep(1,nrow(rr)),rr[,id]-1,nid))
 rr$cens[rr$revnr2==1 & rr$death==0] <- 1

###
formC <- as.formula(paste("Surv(",start,",",stop,",cens)~cluster(",id,")",sep=""))
formD <- as.formula(paste("Surv(",start,",",stop,",",death,")~cluster(",id,")",sep=""))
form1L <- as.formula(paste("Surv(",start,",",stop,",",status,"==",cause,")~Count",cause,"+death+cens+cluster(",id,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",cause,")~cluster(",id,")",sep=""))

 xr <- phreg(form1L,data=rr,no.opt=TRUE,no.var=1)
 cr <- phreg(formC,data=rr,no.opt=TRUE,no.var=1)
 dr <- phreg(formD,data=rr,no.opt=TRUE,no.var=1)

 ### augmenting partioned estimator computing \hat H_i(s,t) for fixed t
 rr$Gctrr <- exp(-cpred(cr$cumhaz,rr[,stop])[,2])

 ### cook-lawless ghosh-lin
 xr0 <- phreg(form1,data=rr,no.opt=TRUE)
 clgl  <- recurrentMarginal(xr0,dr)

  ### censoring weights
  strat <- cr$strata[cr$jumps+1]
  x <- cr
  xx <- x$cox.prep
  S0iD <- S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## survival at t- to also work in competing risks situation
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
    Gc      <- exp(-cumhazD)
  } else Gc <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  S0it <- revcumsumstrata(xx$sign,xx$strata,xx$nstrata)

  ####  First \mu_ipcw(t) \sum_i I(T_i /\ t \leq C_i)/G_c(T_i /\ t ) N_(T_i /\ t) {{{
  x <- xr
  xx <- xr$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## stay with N(D_i) when t is large so no -1 when death
  ## xx$X[,1] er Count1 dvs N(t-)
  Nt <- revcumsumstrata(xx$X[,1]*xx$sign,xx$strata,xx$nstrata)
  Nt <- Nt/Gc
  ## counting N(D) forward in time skal ikke checke ud når man dør N_(D_i) er i spil efter D_i
  NtD <- cumsumstrata(xx$X[,1]*(xx$X[,2]==1)*(xx$sign==1)/Gc,xx$strata,xx$nstrata)
  jump1 <- xx$jumps+1
  timeJ <- xx$time[jump1]
  avNtD <- (NtD+Nt)[jump1]/nid
  strataN1J <- xx$strata[jump1]

  varIPCW1 <- NULL
  seIPCW1 <- NULL

  ### IPCW estimator 
  cumhaz <- cbind(timeJ,avNtD)
# }}}

  ### Partitioned estimator , same as Lin, Lawless & Cook estimator {{{
  cumhazP <- c(cumsumstrata(1/Gc[jump1],strataN1J,xx$nstrata)/nid)
  cumhazP <- cbind(timeJ,cumhazP)

  ### variance of partitioned estimator 

  ### calculate E(H,s,t) = E(H,t) - E(H,s) 
  ### E(H,t) = 1/(S(t)*n) \sum_ \int Y_i(s)/G_c(s) dN_{1i} = 1/(S(t)) (\mu_p(t) - \mu_p(s)) 
  ### \hat H_i for all subjects, and look at together with Hst, eHst 
  ### for censoring martingale 
  rr$Count1s <- rr$Count1^2
  rr$eCount1 <- exp(-rr$Count1)
  rr$CeCount1 <- rr$Count1*exp(-rr$Count1)

  if (!is.null(times)) {
       semuPA <-  muPA <- semuPA.times <-  muPA.times <- rep(0,length(times))

       ww <- fast.approx(timeJ,times,type="left")
       muP.times <- cumhazP[ww,2]
       semuP.times <- clgl$se.mu[ww]
       Aterms <- c("Count1","Count1s","eCount1","CeCount1")[terms]
       modP <- NULL
       if (length(Aterms)>=1) modP <- Aterms
       if (!is.null(x.design)) modP <- c(modP,x.design)
       nterms <- length(modP)
       modP <- paste(modP,collapse="+") 
       form <- as.formula(paste("Surv(",start,",",stop,",cens)~Hst+",modP,"+cluster(id)"))

  for (i in seq_along(times)) {
     timel <- times[i]
     rr$Hst <- revcumsumstrata((rr[,stop]<timel)*(rr[,status]==1)/rr$Gctrr,rr$id-1,nid)
     cr2 <- phreg(form,data=rr,no.opt=TRUE,no.var=1)
     nterms <- cr2$p-1

     dhessian <- cr2$hessianttime
     ###  matrix(apply(dhessian,2,sum),3,3)
     timeb <- which(cr$cumhaz[,1]<timel)
     ### take relevant \sum H_i(s,t) (e_i - \bar e)
     covts <- dhessian[timeb,1+1:nterms,drop=FALSE]
     ### construct relevant \sum (e_i - \bar e)^2
     Pt <- dhessian[timeb,-c((1:(nterms+1)),(1:(nterms))*(nterms+1)+1),drop=FALSE]
     ###  matrix(apply(dhessian[,c(5,6,8,9)],2,sum),2,2)
     gammahat <- .Call("CubeVec",Pt,covts,1,PACKAGE="mets")$XXbeta
     S0 <- cr$S0[timeb]
     gammahat[is.na(gammahat)] <- 0
     gammahat[gammahat==Inf] <- 0
     Gctb <- Gc[cr$jumps][timeb]
     augment.times <- sum(apply(gammahat*cr2$U[timeb,1+1:nterms,drop=FALSE],1,sum))/nid
     ###
     varZ <- matrix(apply(Pt/Gctb^2,2,sum),nterms,nterms)
     gamma <- .Call("CubeVec",matrix(c(varZ),nrow=1),matrix(apply(covts/Gctb,2,sum),nrow=1),1,PACKAGE="mets")$XXbeta
     gamma <- c(gamma)
     gamma[is.na(gamma)] <- 0
     gamma[gamma=Inf] <- 0
     augment <- sum(apply(gamma*t(cr2$U[timeb,1+1:nterms,drop=FALSE])/Gctb,2,sum))/nid
     ###
     muPA[i] <- muP.times[i]+augment
     semuPA[i] <- (semuP.times[i]^2 +(gamma %*% varZ %*% gamma)/nid^2)^.5
     muPA.times[i] <- muP.times[i]+augment.times
     semuPA.times[i] <- (semuP.times[i]^2+sum(gammahat*.Call("CubeVec",Pt,gammahat,0,PACKAGE="mets")$XXbeta)/(nid^2))^.5
  }
  }
# }}}

  return(list(censoring.weights=Gctb,muP.all=cumhazP,Gcjump=Gc[jump1],gamma=gamma,gamma.time=gammahat,times=times,
  muP=muP.times,semuP=semuP.times, muPAt=muPA.times,semuPAt=semuPA.times, muPA=muPA,semuPA=semuPA))
}# }}}

#####' @export
###recurrentMarginalAIPCW <- function(rr,times,km=TRUE,terms=1,idt=1,x.design=NULL,
###   id="id",start="start",stop="stop",status="status",death="death",cause=1,...)
###{# {{{
###
### if (missing(times)) stop("times of estimation must be given")
###
### # to avoid R check warning 
### revnr <- NULL
###
### formsort <- paste("~",id,"+",start)
### dsort(rr) <- as.formula(formsort) 
### rr$revnr <- NULL
### rr$cens <- 0
### rr <- count.history(rr,status=status,id=id)
###
### nid <- max(rr[,id])
### rr$revnr2 <-  c(revcumsumstrata(rep(1,nrow(rr)),rr[,id]-1,nid))
### rr$cens[rr$revnr2==1 & rr$death==0] <- 1
###
######
###formC <- as.formula(paste("Surv(",start,",",stop,",cens)~cluster(",id,")",sep=""))
###formD <- as.formula(paste("Surv(",start,",",stop,",",death,")~cluster(",id,")",sep=""))
###form1L <- as.formula(paste("Surv(",start,",",stop,",",status,"==",cause,")~Count",cause,"+death+cens+cluster(",id,")",sep=""))
###form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",cause,")~cluster(",id,")",sep=""))
###
### xr <- phreg(form1L,data=rr,no.opt=TRUE,no.var=1)
### cr <- phreg(formC,data=rr,no.opt=TRUE,no.var=1)
### dr <- phreg(formD,data=rr,no.opt=TRUE,no.var=1)
###
### ### augmenting partioned estimator computing \hat H_i(s,t) for fixed t
### rr$Gctrr <- exp(-cpred(cr$cumhaz,rr[,stop])[,2])
###
### ### cook-lawless ghosh-lin
### xr0 <- phreg(form1,data=rr,no.opt=TRUE)
### clgl  <- recurrentMarginal(xr0,dr)
###
###  ### censoring weights
###  strat <- cr$strata[cr$jumps+1]
###  x <- cr
###  xx <- x$cox.prep
###  S0iD <- S0i2 <- S0i <- rep(0,length(xx$strata))
###  S0i[xx$jumps+1] <-  1/x$S0
###  S0i2[xx$jumps+1] <- 1/x$S0^2
###  ## survival at t- to also work in competing risks situation
###  if (!km) { 
###    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
###    Gc      <- exp(-cumhazD)
###  } else Gc <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
###  S0it <- revcumsumstrata(xx$sign,xx$strata,xx$nstrata)
###
###  ####  First \mu_ipcw(t) \sum_i I(T_i /\ t \leq C_i)/G_c(T_i /\ t ) N_(T_i /\ t) {{{
###  x <- xr
###  xx <- xr$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
###  S0i[xx$jumps+1] <-  1/x$S0
###  S0i2[xx$jumps+1] <- 1/x$S0^2
###  ## stay with N(D_i) when t is large so no -1 when death
###  ## xx$X[,1] er Count1 dvs N(t-)
###  Nt <- revcumsumstrata(xx$X[,1]*xx$sign,xx$strata,xx$nstrata)
###  Nt <- Nt/Gc
###  ## counting N(D) forward in time skal ikke checke ud når man dør N_(D_i) er i spil efter D_i
###  NtD <- cumsumstrata(xx$X[,1]*(xx$X[,2]==1)*(xx$sign==1)/Gc,xx$strata,xx$nstrata)
###  jump1 <- xx$jumps+1
###  timeJ <- xx$time[jump1]
###  avNtD <- (NtD+Nt)[jump1]/nid
###  strataN1J <- xx$strata[jump1]
###
###  varIPCW1 <- NULL
###  seIPCW1 <- NULL
###
###  ### IPCW estimator 
###  cumhaz <- cbind(timeJ,avNtD)
#### }}}
###
###  ### Partitioned estimator , same as Lin, Lawless & Cook estimator {{{
###  cumhazP <- c(cumsumstrata(1/Gc[jump1],strataN1J,xx$nstrata)/nid)
###  cumhazP <- cbind(timeJ,cumhazP)
###
###  ### variance of partitioned estimator 
###
###  ### calculate E(H,s,t) = E(H,t) - E(H,s) 
###  ### E(H,t) = 1/(S(t)*n) \sum_ \int Y_i(s)/G_c(s) dN_{1i} = 1/(S(t)) (\mu_p(t) - \mu_p(s)) 
###  ### \hat H_i for all subjects, and look at together with Hst, eHst 
###  ### for censoring martingale 
###  rr$Count1s <- rr$Count1^2
###  rr$eCount1 <- exp(-rr$Count1)
###  rr$CeCount1 <- rr$Count1*exp(-rr$Count1)
###
###  if (!is.null(times)) {
###       semuPA <-  muPA <- semuPA.times <-  muPA.times <- rep(0,length(times))
###
###       ww <- fast.approx(timeJ,times,type="left")
###       muP.times <- cumhazP[ww,2]
###       semuP.times <- clgl$se.mu[ww]
###       Aterms <- c("Count1","Count1s","eCount1","CeCount1")[terms]
###       modP <- NULL
###       if (length(Aterms)>=1) modP <- Aterms
###       if (!is.null(x.design)) modP <- c(modP,x.design)
###       nterms <- length(modP)
###       modP <- paste(modP,collapse="+") 
###       form <- as.formula(paste("Surv(",start,",",stop,",cens)~Hst+",modP,"+cluster(id)"))
###
###  for (i in seq_along(times)) {
###     timel <- times[i]
###     rr$Hst <- revcumsumstrata((rr[,stop]<timel)*(rr[,status]==1)/rr$Gctrr,rr$id-1,nid)
###     cr2 <- phreg(form,data=rr,no.opt=TRUE,no.var=1)
###     nterms <- cr2$p-1
###
###     dhessian <- cr2$hessianttime
###     ###  matrix(apply(dhessian,2,sum),3,3)
###     timeb <- which(cr$cumhaz[,1]<timel)
###     ### take relevant \sum H_i(s,t) (e_i - \bar e)
###     covts <- dhessian[timeb,1+1:nterms,drop=FALSE]
###     ### construct relevant \sum (e_i - \bar e)^2
###     Pt <- dhessian[timeb,-c((1:(nterms+1)),(1:(nterms))*(nterms+1)+1),drop=FALSE]
###     ###  matrix(apply(dhessian[,c(5,6,8,9)],2,sum),2,2)
###     gammahat <- .Call("CubeVec",Pt,covts,1,PACKAGE="mets")$XXbeta
###     S0 <- cr$S0[timeb]
###     gammahat[is.na(gammahat)] <- 0
###     gammahat[gammahat==Inf] <- 0
###     Gctb <- Gc[cr$jumps][timeb]
###     augment.times <- sum(apply(gammahat*cr2$U[timeb,1+1:nterms,drop=FALSE],1,sum))/nid
###     mterms <- length(terms)
###     mterms <- nterms
###     ###
###     varZ <- matrix(apply(Pt/Gctb^2,2,sum),mterms,mterms)
###     gamma <- .Call("CubeVec",matrix(c(varZ),nrow=1),matrix(apply(covts/Gctb,2,sum),nrow=1),1,PACKAGE="mets")$XXbeta
###     gamma <- c(gamma)
###     gamma[is.na(gamma)] <- 0
###     gamma[gamma=Inf] <- 0
###     augment <- sum(apply(gamma*t(cr2$U[timeb,1+1:nterms,drop=FALSE])/Gctb,2,sum))/nid
###     ###
###     muPA[i] <- muP.times[i]+augment
###     semuPA[i] <- (semuP.times[i]^2 +(gamma %*% varZ %*% gamma)/nid^2)^.5
###     muPA.times[i] <- muP.times[i]+augment.times
###     semuPA.times[i] <- (semuP.times[i]^2+sum(gammahat*.Call("CubeVec",Pt,gammahat,0,PACKAGE="mets")$XXbeta)/(nid^2))^.5
###  }
###  }
#### }}}
###
###  return(list(censoring.weights=Gctb,muP.all=cumhazP,Gcjump=Gc[jump1],gamma=gamma,gamma.time=gammahat,times=times,
###  muP=muP.times,semuP=semuP.times, muPAt=muPA.times,semuPAt=semuPA.times, muPA=muPA,semuPA=semuPA))
###}# }}}

##' @export
recurrentMarginalIPCW <- function(rr,km=TRUE,times=NULL,...)
{# {{{

 # to ovoid R check warning 
 death <- revnr <- NULL

 rr$revnr <- NULL
 rr$cens <- 0
 rr <- count.history(rr)
 dsort(rr) <- ~id-start
 nid <- max(rr$id)
 rr$revnr <- cumsumstrata(rep(1,nrow(rr)),rr$id-1,nid)
 dsort(rr) <- ~id+start
 rr <- dtransform(rr,cens=1,revnr==1 & death==0)

  xr <- phreg(Surv(entry,time,status==1)~Count1+death+cluster(id),data=rr,no.opt=TRUE,no.var=1)
  cr <- phreg(Surv(entry,time,cens)~cluster(id),data=rr)
  dr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)

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
  x <- xr
  xx <- xr$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## stay with N(D_i) when t is large so no -1 when death
  signm <- xx$sign
  signm[xx$X[,2]==1 & xx$sign==-1] <- 0
  Nt <- revcumsumstrata(xx$X[,1]*xx$sign,xx$strata,xx$nstrata)
  Nt <- Nt/St
  ## counting N(D) forward in time skal ikke checke ud når man dør N_(D_i) er i spil efter D_i
  NtD <- cumsumstrata(xx$X[,1]*(xx$X[,2]==1)*(xx$sign==1)/St,xx$strata,xx$nstrata)
  risk <- revcumsumstrata(xx$sign,xx$strata,xx$nstrata)
  timeJ <- xx$time[xx$jumps+1]
  avNtD <- (NtD+Nt)[xx$jumps+1]/nid
  xxJ <- xx$jumps+1

  cumhaz <- cbind(timeJ,avNtD)

  return(list(cumhaz=cumhaz))
}# }}}

##' @export
recmarg <- function(recurrent,death,Xr=NULL,Xd=NULL,km=TRUE,...)
{# {{{
  xr <- recurrent
  dr <- death 

  if (!is.null(Xr)) rr <- exp(sum(xr$coef * Xr)) else rr <- 1
  if (!is.null(Xd)) rrd <- exp(sum(dr$coef * Xd)) else rrd <- 1

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD*rrd)
  } else St <- exp(rrd*c(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  ###
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
  mu <- rr*cumhazDR[,2]
# }}}

 varrs <- data.frame(mu=mu,time=xr$time,strata=xr$strata,St=St)
 varrs  <-  varrs[c(xr$cox.prep$jumps)+1,]

 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,time=varrs$time,
	     St=varrs$St,cumhaz=cbind(varrs$time,varrs$mu),
             strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
	     strata.name=xr$strata.name,strata.level=recurrent$strata.level)
 return(out)
}# }}}

##' @export
squareintHdM <- function(phreg,ft=NULL,fixbeta=NULL,...)
{# {{{
###  sum_k ( int_0^t f(s)/S_0^r(s) dM_k.^r(s) )^2
###  strata "r" from object and "k" id from cluster 
  if (!inherits(phreg,"phreg")) stop("Must be phreg object\n"); 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if ((phreg$no.opt) | is.null(phreg$coef)) fixbeta<- 1 else fixbeta <- 0

  x <- phreg
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
  U[xx$jumps+1,] <- x$U
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  if (is.null(ft))  ft <- rep(1,length(xx$time))
  cumS0i2 <- c(cumsumstrata(ft*S0i2,xx$strata,xx$nstrata))
  if (fixbeta==0) {
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  Ht <- apply(ft*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  } else Ht <- NULL
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
  else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  ### also weights 
  w <- c(xx$weights)
  xxx <- (ft*S0i-rr*cumS0i2)
  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
  ss <-  revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare*cumS0i2^2
  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*cumS0i2
  varA1 <- c(ssf+ss-2*covv)

  vbeta <- betaiidR <- NULL
  if (fixbeta==0) {# {{{
     invhess <- -solve(x$hessian)
     MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-EdLam0)*rr*c(xx$weights)
     UU <- apply(MGt,2,sumstrata,id,mid)
     betaiidR <- UU %*% invhess
     vbeta <- crossprod(betaiidR)
     varbetat <-   rowSums((Ht %*% vbeta)*Ht)
     ### writing each beta for all individuals 
     betakt <- betaiidR[id+1,,drop=FALSE]
     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- c(covk2)*cumS0i2
     covv <- covk1-covk2
     varA1 <- varA1+varbetat-2*apply(covv*Ht,1,sum)
  }# }}}

  return(list(xx=xx,Ht=Ht,varInt=varA1,xxx=xxx,rr=rr,
	      cumhaz=cumhaz,cumS0i2=cumS0i2,mid=mid,id=id,
	      betaiid=betaiidR,vbeta=vbeta,covv=covv))
} # }}}

##' @export
covIntH1dM1IntH2dM2 <- function(square1,square2,fixbeta=1,mu=NULL)
{# {{{

 ### strata and id same for two objects 
 xx <- square1$xx; xx2 <- square2$xx
 xxxR <- square1$xxx;     xxxD1 <- square2$xxx
 rrR  <- square1$rr;       rrD1 <- square2$rr
 id   <- id1 <- square1$id; id2 <- square2$id
 mid  <- square1$mid; w <- c(xx$weights)

 if (is.null(mu)) mu <- rep(1,length(xx$strata))

 cov12 <-  c(cumsumidstratasumCov(xxxR,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)
 cov122 <- c(revcumsumidstratasumCov(w*rrR,w*rrD1,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(square1$cumS0i2*square2$cumS0i2)
 cov123 <- covfridstrataCov(xxxR,w*rrR,xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(square2$cumS0i2)
 cov124 <- covfridstrataCov(xxxD1,w*rrD1,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(square1$cumS0i2)
 cov12A <- c(cov12+cov122+cov123+cov124)
 
 test <- 0
 if (test==1) {# {{{
	 print("________cov cov ___________________________"); 
	 print(summary(c(xxxR))); print(summary(c(xxxD1))); 
	 print(summary(c(rrD1))); print(summary(c(rrR)))
	 print(summary(c(square1$cumS0i2))); print(summary(c(square2$cumS0i2)));
	 print("-----------"); 
	 print(summary(cov12)); 
	 print(summary(cov122)); 
	 print(summary(c(cov123))); 
	 print(summary(c(cov124))); 
	 print("______________________________________"); 
	 jumps <- c(square1$xx$jumps,square2$xx$jumps)+1
	 print(summary(jumps))
         print(summary(cov12[jumps])); 
	 print(summary(cov122[jumps])); 
	 print(summary(c(cov123[jumps]))); 
	 print(summary(c(cov124[jumps]))); 
 }# }}}

 cov12aa <- 0
 if (fixbeta==0) {
 ### covariances between different terms and  beta's 
 # {{{
	 betaiidR <- square1$betaiid; betaiidD <- square2$betaiid
	 HtR <- square1$Ht; HtD <- square2$Ht
	 covbetaRD <- t(betaiidR) %*% betaiidD
	 covbeta <-   -1*rowSums((HtR %*% covbetaRD)*HtD)

	 ### cov12 wrt betaD and betaR
	 betakt <- betaiidD[id1+1,,drop=FALSE]
	 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id1,mid,xx$strata,xx$nstrata,type="sum")
	 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id1,mid,xx$strata,xx$nstrata,type="lagsum")
	 covk2 <- c(covk2)*c(square1$cumS0i2)
	 covRD12 <- apply((covk1-covk2)*HtD,1,sum)
	 betakt <- betaiidR[id2+1,,drop=FALSE]
	 covk1 <-apply(xxxD1*betakt,2,cumsumidstratasum,id2,mid,xx2$strata,xx$nstrata,type="sum")
	 covk2 <-apply(w*rrD1*betakt,2,revcumsumidstratasum,id2,mid,xx2$strata,xx$nstrata,type="lagsum")
	 covk2 <- c(covk2)*c(square2$cumS0i2)
	 covRD21 <- apply((covk1-covk2)*HtR,1,sum)
	 cov12aa <- 2*(covbeta + covRD12+covRD21)
	 test <- 0
	 if (test==1) {
		 print("--------------")
		 print(summary(2*mu*covbeta))
		 print(summary(2*mu*covRD12))
		 print(summary(2*mu*covRD21))
		 print("--------------")
	 }
 } # }}}

 cov12 <- (cov12A-cov12aa)*mu

 return(list(cov=cov12,cov12A=cov12A*mu,covbeta=cov12aa*mu))
} # }}}

##' @export
tie.breaker <- function(data,stop="time",start="entry",status="status",id=NULL,ddt=NULL,exit.unique=TRUE)
{# {{{

   if (!is.null(id)) id <- data[,id]
   ord <- 1:nrow(data)
   stat <- data[,status]
   time <- data[,stop]
   dupexit <- duplicated(time)
   time1 <- data[stat==1,stop]
   time0 <- data[stat!=1,stop]
   lt0 <- length(time0)
   ddp <- duplicated(c(time0,time1))
   if (exit.unique) ties <-ddp[(lt0+1):nrow(data)] else ties <- duplicated(c(time1))
   nties <- sum(ties)
   ordties <- ord[stat==1][ties]
   if (is.null(ddt)) {
	   abd <- abs(diff(data[,stop]))
	   abd <- min(abd[abd>0])
	   ddt <- abd*0.5
   }
   time[ordties] <- time[ordties]+runif(nties)*ddt

   data[ordties,stop] <- time[ordties]
   ties <- (ord %in% ordties)
   if (!is.null(id)) {
   lagties <- dlag(ties)
   ### also move next start time if id the same 
   change.start <- lagties==TRUE & id==dlag(id)
   change.start[is.na(change.start)] <- FALSE
   ocs <- ord[change.start]
   data[ocs,start] <- data[ocs-1,stop]
   data[,"tiebreaker"] <- FALSE
   data[ocs,"tiebreaker"] <- TRUE
   }
   
   return(data)
 } # }}}


##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Must give hazard of death and recurrent events.  Possible with two
##' event types and their dependence can be specified but the two recurrent events need
##' to have the same random effect,  simRecurrentII more flexible !  
##'
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param death.cumhaz cumulative hazard of death 
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param cens rate of exponential on total time i.e. on death time-scale 
##' @param max.recurrent limits number recurrent events to 100
##' @param dhaz rate for death hazard if it is extended to time-range of first event 
##' @param dependence  =0 independence, =1 all share same random effect with variance var.z
##'                    =2 random effect exp(normal) with correlation structure from cor.mat,
##'                    first random effect is z1 and shared for a possible second cause,  second random effect is for death 
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##'
##'  ######################################################################
##'  ### simulating simple model that mimicks data 
##'  ######################################################################
##'  rr <- simRecurrent(5,base1,death.cumhaz=dr)
##'  dlist(rr,.~id,n=0)
##'
##'  rr <- simRecurrent(1000,base1,death.cumhaz=dr)
##'  par(mfrow=c(1,3))
##'  showfitsim(causes=1,rr,dr,base1,base1)
##'
##' ######################################################################
##' ### simulating simple model 
##' ### random effect for all causes (Z shared for death and recurrent) 
##' ######################################################################
##'
##'  rr <- simRecurrent(1000,base1,death.cumhaz=dr,dependence=1,var.gamma=0.4)
##'  ### marginals do fit after input after integrating out
##'  par(mfrow=c(2,2))
##'  showfitsim(causes=1,rr,dr,base1,base1)
##'
##' @aliases showfitsim  simRecurrentGamma covIntH1dM1IntH2dM2 squareintHdM 
##' @export
simRecurrent <- function(n,cumhaz,death.cumhaz=NULL,gap.time=FALSE,cens=NULL,
	 max.recurrent=100,dhaz=NULL,dependence=0,var.z=2,cor.mat=NULL,...) 
{# {{{
  status <- fdeath <-  dtime <- NULL ## to avoid R-check 

  ### drawing relative risk frailty terms to generate dependence
  if (dependence==0) { z1 <- z2 <- zd <- rep(1,n) # {{{
     } else if (dependence==1) {
###	      zz <- rgamma(n,1/var.gamma[1])*var.gamma[1]
	      zz <- exp(rnorm(n,1)*var.z[1]^.5)
	      z1 <- zz; z2 <- zz; zd <- zz
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(n*2),n,2)
	      z <- (z%*% chol(covv))
	      z1 <- exp(z[,1]); zd <- exp(z[,2])
	      apply(exp(z),2,mean); cov(exp(z))
      } else if (dependence==3) {
	      zz <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- zz; z2 <- zz; zd <- rep(1,n) 
      }      
  # }}}

  cumhaz <- rbind(c(0,0),cumhaz)

  ## extend cumulative for death to full range  of cause 1
  if (!is.null(death.cumhaz)) {
    out <- extendCums(cumhaz,death.cumhaz)
    cumhaz <- out$cum1
    cumhazd <- out$cum2
  }

  ll <- nrow(cumhaz)
  max.time <- tail(cumhaz[,1],1)
  rc <- 1

### recurrent first time
  tall <- rchaz(cumhaz,z1)
  tall$id <- 1:n
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- rchaz(cumhazd,zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  nrr <- n
  i <- 1; 
  while (any((tt$time<tt$dtime) & (tt$status!=0) & (i < max.recurrent))) {
	  i <- i+1
	  still <- subset(tt,time<dtime & status!=0)
	  ## start at where we are or "0" for gaptime
          tt <- rchaz(cumhaz,z1[still$id],entry=(1-gap.time)*still$time)
	  if (gap.time) { 
		  tt$entry <- still$time
		  tt$time <- tt$time+still$time
	  }
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt,row.names=NULL)
	  nrr <- nrr+nt
  }
  dsort(tall) <- ~id+entry+time

  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz

  return(tall)
  }# }}}


##' @export
simRecurrentGamma <- function(n,haz=0.5,death.haz=0.1,haz2=0.1,max.recurrent=100,var.z=2,times=5000) 
{# {{{

  status <- dtime <- NULL ## to avoid R-check 

  max.time <- times
  cumhaz1 <- rbind(c(0,0),c(times,times*haz))
  cumhaz2 <- rbind(c(0,0),c(times,times*haz2))
  death.cumhaz <- rbind(c(0,0),c(times,death.haz))
  z <- rgamma(1/var.z)*var.z

  cumhaz <- cbind(times,cumhaz1+cumhaz2)

### recurrent first time
  tall <- rchaz(cumhaz,z)
  tall$id <- 1:n
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- rchaz(cumhazd,n)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
  } else { tall$dtime <- max.time; tall$fdeath <- 0; cumhazd <- NULL }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=1,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  i <- 1; 
  while (any((tt$time<tt$dtime) & (tt$status!=0) & (i < max.recurrent))) {
	  i <- i+1
	  still <- subset(tt,time<dtime & status!=0)
          tt <- rchaz(cumhaz,z[still$id],entry=still$time)
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=1,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt,row.names=NULL)
  }
  dsort(tall) <- ~id+entry+time

  ### cause 2 is there then decide if jump is 1 or 2
  if (!is.null(haz2)) {# {{{
      p2t <- haz2/(haz+haz2)
      tall$p2t <- p2t
      tall$status <- (1+rbinom(nrow(tall),1,p2t))*(tall$status>=1)
  }# }}}

  tall$start <- tall$entry
  tall$stop  <- tall$time
  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  ### haz*haz2*(var.z+1)

  return(tall)
}# }}}

##' Simulation of recurrent events data based on cumulative hazards II 
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Must give hazard of death and two recurrent events.  Possible with two
##' event types and their dependence can be specified but the two recurrent events need
##' to share random effect. Based on drawing the from cumhaz and cumhaz2 and 
##' taking the first event rather
##' the cumulative and then distributing it out. Key advantage of this is that 
##' there is  more flexibility wrt random effects 
##'
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param cumhaz2  cumulative hazard of recurrent events  of type 2
##' @param death.cumhaz cumulative hazard of death 
##' @param r1 potential relative risk adjustment of rate 
##' @param r2 potential relative risk adjustment of rate
##' @param rd potential relative risk adjustment of rate
##' @param rc potential relative risk adjustment of rate
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param dhaz rate for death hazard if it is extended to time-range of first event 
##' @param haz2 rate of second cause  if it is extended to time-range of first event 
##' @param dependence 0:independence; 1:all share same random effect with variance var.z; 2:random effect exp(normal) with correlation structure from cor.mat; 3:additive gamma distributed random effects, z1= (z11+ z12)/2 such that mean is 1 , z2= (z11^cor.mat(1,2)+ z13)/2, z3= (z12^(cor.mat(2,3)+z13^cor.mat(1,3))/2, with z11 z12 z13 are gamma with mean and variance 1 , first random effect is z1 and for N1 second random effect is z2 and for N2 third random effect is for death  
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param cens rate of censoring exponential distribution
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##'
##'  cor.mat <- corM <- rbind(c(1.0, 0.6, 0.9), c(0.6, 1.0, 0.5), c(0.9, 0.5, 1.0))
##' 
##' ######################################################################
##' ### simulating simple model that mimicks data 
##' ### now with two event types and second type has same rate as death rate
##' ######################################################################
##' set.seed(100)
##' rr <- simRecurrentII(1000,base1,base4,death.cumhaz=dr)
##' dtable(rr,~death+status)
##' par(mfrow=c(2,2))
##' showfitsim(causes=2,rr,dr,base1,base4)
##'
##' @export
simRecurrentII <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,r1=NULL,r2=NULL,rd=NULL,rc=NULL,
    gap.time=FALSE,max.recurrent=100,dhaz=NULL,haz2=NULL,dependence=0,var.z=0.22,cor.mat=NULL,cens=NULL,...) 
  {# {{{

  status <- fdeath <- dtime <- NULL # to avoid R-check 

  if (dependence==0) { z <- z1 <- z2 <- zd <- rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- z; z2 <- z; zd <- z
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
      } else if (dependence==4) {
	      zz <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- zz; z2 <- zz; zd <- rep(1,n) 
	      z <- z1
      }      else stop("dependence 0-4"); # }}}

   if (is.null(r1)) r1 <- rep(1,n)
   if (is.null(r2)) r2 <- rep(1,n)
   if (is.null(rd)) rd <- rep(1,n)
   if (is.null(rc)) rc <- rep(1,n)

  cumhaz <- rbind(c(0,0),cumhaz)

  ## extend cumulative for death to full range  of cause 1
  if (!is.null(death.cumhaz)) {
     out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz),NULL)
     cumhaz <- out$cum1
     cumhaz2 <- out$cum2
     cumhazd <- out$cum3
  } else {
     out <- extendCums(list(cumhaz,cumhaz2),NULL)
     cumhaz <- out$cum1
     cumhaz2 <- out$cum2
  }

  ll <- nrow(cumhaz)
  max.time <- tail(cumhaz[,1],1)

### recurrent first time
  tall1 <- rchaz(cumhaz,z1*r1)
  tall2 <- rchaz(cumhaz2,z2*r2)
  tall <- tall1 
  tall$status <- ifelse(tall1$time<tall2$time,tall1$status,2*tall2$status)
  tall$time <- ifelse(tall1$time<tall2$time,tall1$time,tall2$time)
  tall$id <- 1:n
  tall$rr2 <- tall2$rr
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- rchaz(cumhazd,zd*rd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  ### setting aside memory 
  tt1 <- tt2 <- tt
  i <- 1; 
  while (any((tt$time<tt$dtime) & (tt$status!=0) & (i < max.recurrent))) {
	  i <- i+1
	  still <- subset(tt,time<dtime & status!=0)
	  nn <- nrow(still)
          tt1 <- rchaz(cumhaz,r1[still$id]*z1[still$id],entry=(1-gap.time)*still$time)
          tt2 <- rchaz(cumhaz2,r2[still$id]*z2[still$id],entry=(1-gap.time)*still$time)
	  tt <- tt1
          tt$status <- ifelse(tt1$time<=tt2$time,tt1$status,2*tt2$status)
          tt$time <-   ifelse(tt1$time<=tt2$time,tt1$time,tt2$time)
	  tt$rr2 <- tt2$rr
	  if (gap.time) {
		  tt$entry <- still$time
		  tt$time  <- tt$time+still$time
	  }
          ###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt[1:nn,],row.names=NULL)
  } 
  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"z") <- z

  return(tall)
  }# }}}


##' @export
showfitsim <- function(causes=2,rr,dr,base1,base4,which=1:3) 
{# {{{
if (1 %in% which) {
  drr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)
  basehazplot.phreg(drr,ylim=c(0,8))
  lines(dr,col=2)
}
###
if (2 %in% which) {
  xrr <- phreg(Surv(entry,time,status==1)~cluster(id),data=rr)
  basehazplot.phreg(xrr,add=TRUE)
###  basehazplot.phreg(xrr)
  lines(base1,col=2)
  if (causes>=2) {
	  xrr2 <- phreg(Surv(entry,time,status==2)~cluster(id),data=rr)
	  basehazplot.phreg(xrr2,add=TRUE)
	  lines(base4,col=2)
  }
  }
if (3 %in% which) {
  meanr1 <-   recurrentMarginal(xrr,drr)
  basehazplot.phreg(meanr1,se=TRUE)
  if (causes>=2) {
	  meanr2 <-   recurrentMarginal(xrr2,drr)
	  basehazplot.phreg(meanr2,se=TRUE,add=TRUE,col=2)
  }
}
}# }}}

##' Simulation of illness-death model 
##'
##' Simulation of illness-death model 
##'
##' simMultistate with different death intensities from states 1 and 2 
##'
##' Must give cumulative hazards on some time-range 
##'
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of going from state 1 to 2.
##' @param cumhaz2  cumulative hazard of going from state 2 to 1. 
##' @param death.cumhaz cumulative hazard of death from state 1. 
##' @param death.cumhaz2 cumulative hazard of death from state 2.
##' @param rr  relative risk adjustment for cumhaz
##' @param rr2  relative risk adjustment for cumhaz2
##' @param rd  relative risk adjustment for death.cumhaz
##' @param rd2  relative risk adjustment for death.cumhaz2
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param dependence 0:independence; 1:all share same random effect with variance var.z; 2:random effect exp(normal) with correlation structure from cor.mat; 3:additive gamma distributed random effects, z1= (z11+ z12)/2 such that mean is 1 , z2= (z11^cor.mat(1,2)+ z13)/2, z3= (z12^(cor.mat(2,3)+z13^cor.mat(1,3))/2, with z11 z12 z13 are gamma with mean and variance 1 , first random effect is z1 and for N1 second random effect is z2 and for N2 third random effect is for death  
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param cens rate of censoring exponential distribution
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' dr2 <- drcumhaz
##' dr2[,2] <- 1.5*drcumhaz[,2]
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##' cens <- rbind(c(0,0),c(2000,0.5),c(5110,3))
##'
##' iddata <- simMultistate(100,base1,base1,dr,dr2,cens=cens)
##' dlist(iddata,.~id|id<3,n=0)
##'  
##' ### estimating rates from simulated data  
##' c0 <- phreg(Surv(start,stop,status==0)~+1,iddata)
##' c3 <- phreg(Surv(start,stop,status==3)~+strata(from),iddata)
##' c1 <- phreg(Surv(start,stop,status==1)~+1,subset(iddata,from==2))
##' c2 <- phreg(Surv(start,stop,status==2)~+1,subset(iddata,from==1))
##' ###
##' par(mfrow=c(2,3))
##' bplot(c0)
##' lines(cens,col=2) 
##' bplot(c3,main="rates 1-> 3 , 2->3")
##' lines(dr,col=1,lwd=2)
##' lines(dr2,col=2,lwd=2)
##' ###
##' bplot(c1,main="rate 1->2")
##' lines(base1,lwd=2)
##' ###
##' bplot(c2,main="rate 2->1")
##' lines(base1,lwd=2)
##'  
##' @aliases extendCums 
##' @export
simMultistate <- function(n,cumhaz,cumhaz2,death.cumhaz,death.cumhaz2,
		    rr=NULL,rr2=NULL,rd=NULL,rd2=NULL,
		    gap.time=FALSE,max.recurrent=100,
		    dependence=0,var.z=0.22,cor.mat=NULL,cens=NULL,...) 
{# {{{

  fdeath <- dtime <- NULL # to avoid R-check 
  status <- dhaz <- NULL; dhaz2 <- NULL

  if (dependence==0) { z <- z1 <- z2 <- zd  <- zd2 <-  rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
###	      z <- exp(rnorm(n,1)*var.z[1]^.5)
	      z1 <- z; z2 <- z; zd <- z
	      if (!is.null(cor.mat)) { zd <- rep(1,n); }
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
###	      print(summary(z))
###	      print(cor(z))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
###	      print(summary(z))
###	      print(cor(z))
      } else stop("dependence 0-3"); # }}}

  ## covariate adjustment 
  if (is.null(rr))  rr <- z1; 
  if (is.null(rr2)) rr2 <- z2; 
  if (is.null(rd))  rd  <- zd; 
  if (is.null(rd2)) rd2 <- zd2; 

  ll <- nrow(cumhaz)
  ### extend of cumulatives
  cumhaz <- rbind(c(0,0),cumhaz)
  cumhaz2 <- rbind(c(0,0),cumhaz2)
  death.cumhaz <- rbind(c(0,0),death.cumhaz)
  death.cumhaz2 <- rbind(c(0,0),death.cumhaz2)

  haz <- haz2 <- NULL
  ## range max of cumhaz and cumhaz2 

  if (!is.null(cens)) {
	  if (is.matrix(cens))  {
             out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz,death.cumhaz2,cens),NULL)
   	     cens <- out$cum5
	  }
  } else {
     out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz,death.cumhaz2),NULL)
  }
  cumhaz <- out$cum1
  cumhaz2 <- out$cum2
  cumhazd <- out$cum3
  cumhazd2 <- out$cum4
  max.time <- tail(cumhaz[,1],1)

  tall <- rcrisk(cumhaz,cumhazd,rr,rd,cens=cens)
  tall$id <- 1:n
  ### fixing the first time to event
  tall$death <- 0
  ### cause 2 is death state 3, cause 1 is state 2
  tall <- dtransform(tall,status=3,status==2)
  tall <- dtransform(tall,death=1,status==3)
  tall <- dtransform(tall,status=2,status==1)
  ### dead or censored
  deadid <- (tall$status==3 | tall$status==0)
  tall$from <- 1
  tall$to <- tall$status
  ## id's that are dead: tall[deadid,]
  ## go furhter with those that are not yet dead  or censored
  tt <- tall[!deadid,,drop=FALSE]
  ## also check that we are before max.time
  tt <- subset(tt,tt$time<max.time)

  i <- 1; 
  while ( (nrow(tt)>0) & (i < max.recurrent)) {# {{{
	  i <- i+1
	  nn <- nrow(tt)

	  z1r <- rr[tt$id]
	  zdr <- rd[tt$id]
	  z2r <- rr2[tt$id]
	  zd2r <- rd2[tt$id]

	  if (i%%2==0) { ## in state 2
	  ## out of 2 for those in 2
          tt1 <- rcrisk(cumhaz2,cumhazd2,z2r,zd2r,entry=tt$time,cens=cens)
          tt1$death <- 0
	  ### status 2 is death state 3, status 1 is state 1
	  tt1 <- dtransform(tt1,status=3,status==2)
	  tt1 <- dtransform(tt1,death=1,status==3)
	  tt1$from <- 2
	  tt1$to <- tt1$status
	  ## take id from tt
	  tt1$id <-  tt$id
	  ### add to data 
	  tall <- rbind(tall,tt1,row.names=NULL)

          deadid <- (tt1$status==3 | tt1$status==0)
	  ### those that are still under risk 
	  tt <- tt1[!deadid,,drop=FALSE]
	  ## also keep only those before max.time
          tt <- subset(tt,tt$time<max.time)
		  
	  } else { ## in state 1
	  ## out of 1 for those in 1
          tt1 <- rcrisk(cumhaz,cumhazd,z1r,zdr,entry=tt$time,cens=cens)

          tt1$death <- 0
	  ### status 2 is death state 3, status 1 is state 2
	  tt1 <- dtransform(tt1,status=3,status==2)
	  tt1 <- dtransform(tt1,death=1,status==3)
	  tt1 <- dtransform(tt1,status=2,status==1)
	  tt1$from <- 1
	  tt1$to <- tt1$status
	  tt1$id <-  tt$id

	  ### add to data 
	  tall <- rbind(tall,tt1,row.names=NULL)

	  ## take id from tt
          deadid <- (tt1$status==3 | tt1$status==0)
	  ### those that are still under risk 
	  tt <- tt1[!deadid,,drop=FALSE]
	  ### also only keep those before max.time
          tt <- subset(tt,tt$time<max.time)
	  }

  }  # }}}

  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"death.cumhaz2") <- cumhazd2
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"cens.cumhaz") <- cens
  attr(tall,"z") <- z

  return(tall)
  }# }}}

##' @export
extendCums <- function(cumA,cumB,haza=NULL)
{# {{{
 ## setup as list to run within loop
 if (!is.null(cumB)) {cumA <- list(cumA,cumB); } else cumA <- c(cumA,cumB)

 maxx <- unlist(lapply(cumA,function(x) tail(x,1)[1]))
 mm <- which.max(maxx)
 nn <- length(cumA)

for (i in seq(nn)[-mm]) {
  cumB <- as.matrix(cumA[[i]]); 
  cumB <- rbind(c(0,0),cumB); 

  ### linear extrapolation of mortality using given dhaz or 
  if (tail(cumB[,1],1)<maxx[mm]) {
      tailB <- tail(cumB,1)
      cumlast <- tailB[2]
      timelast <- tailB[1]
      if (is.null(haza)) hazb <- cumlast/timelast else hazb <- haza[i]
      cumB <- rbind(cumB,c(maxx[mm],cumlast+hazb*(maxx[mm]-timelast))) 
  }
  cumA[[i]] <- cumB
}
 cumA[[mm]] <- rbind(c(0,0),cumA[[mm]])

  return( setNames(cumA,paste("cum",seq(nn),sep="")))
}# }}}

##' Simulation of recurrent events data based on cumulative hazards: Two-stage model  
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Model is constructed such that marginals are on specified form by linear approximations
##' of cumulative hazards that are on a specific form to make them equivalent to marginals
##' after integrating out over survivors. Therefore E(dN_1 | D>t) = cumhaz, 
##' E(dN_2 | D>t) = cumhaz2,  and hazard of death is death.cumhazard 
##'
##' Must give hazard of death and two recurrent events.  Hazard of death is death.cumhazard  two
##' event types and their dependence can be specified but the two recurrent events need
##' to share random effect. 
##' 
##' Random effect for  death Z.death=(Zd1+Zd2), Z1=(Zd1^nu1) Z12,  Z2=(Zd2^nu2) Z12^nu3
##' \deqn{Z.death=Zd1+Zd2}  gamma distributions 
##' \deqn{Zdj}  gamma distribution  with mean parameters (sharej), vargamD,  share2=1-share1
##' \deqn{Z12}  gamma distribution with mean 1 and variance vargam12
##' 
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param cumhaz2  cumulative hazard of recurrent events  of type 2
##' @param death.cumhaz cumulative hazard of death 
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param nu powers of random effects where nu > -1/shape 
##' @param share1 how random effect for death splits into two parts 
##' @param vargamD variance of random effect  for death 
##' @param vargam12 shared random effect for N1 and N2 
##' @param cens rate of censoring exponential distribution
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##'
##' rr <- simRecurrentTS(1000,base1,base4,death.cumhaz=dr)
##' dtable(rr,~death+status)
##' showfitsim(causes=2,rr,dr,base1,base4)
##'
##' @export
simRecurrentTS <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,
		    nu=rep(1,3),share1=0.3,vargamD=2,vargam12=0.5,
		    gap.time=FALSE,max.recurrent=100,cens=NULL,...) 
{# {{{

k <- 1
nu1 <- nu[1]; nu2 <- nu[2]; nu3 <- nu[3]
###nu1 <- 1; nu2 <- 1; nu3 <- 0.4
share2 <- (1-share1)
vargam <- vargamD
vargam12 <- 0.5
agam1 <- share1/vargam
agam2 <- share2/vargam
betagam=1/vargam
gamma1  <- rep(rgamma(n,agam1)*vargam,each=k)
gamma2  <- rep(rgamma(n,agam2)*vargam,each=k)
agam12 <- 1/vargam12
betagam12 <- 1/vargam12
gamma12 <- rep(rgamma(n,agam12)*vargam12,each=k)
agamD <- agam1+agam2
z1 <- (gamma1^nu1)*gamma12
z2 <- (gamma2^nu2)*gamma12^nu3
gamD <- gamma1+gamma2
zd <- gamD
egamma12nu3 <- (gamma(agam12+nu3)/gamma(agam12))*1/(betagam12)^nu3
zs <- cbind(z1,z2,zd)

  status <- fdeath <- dtime <- NULL # to avoid R-check 
  dhaz <- haz2 <- dhaz <- NULL

 ll <- nrow(cumhaz)
 max.time <- tail(cumhaz[,1],1)

 ################################################################
 ### approximate hazards to make marginals fit (approximately)
 ################################################################
 orig.death <- death.cumhaz
 base1 <- death.cumhaz
 gt <- exp(vargam*base1[,2]) 
 dtt <- diff(c(0,base1[,1]))
 lams <- (diff(c(0,base1[,2]))/dtt)*gt
 death.cumhaz <- cbind(base1[,1],cumsum(dtt*lams))

 base1 <- cumhaz
 dbase1 <- cpred(rbind(c(0,0),death.cumhaz),base1[,1])[,2]
 dtt <- diff(c(0,base1[,1]))
 gt <- (gamma(agam1+nu1)/gamma(agam1))*(1/(betagam+dbase1))^nu1
 lams <- (diff(c(0,base1[,2]))/dtt)*(1/gt)
 cumhaz <- cbind(base1[,1],cumsum(dtt*lams))

 base1 <- cumhaz2
 dbase1 <- cpred(rbind(c(0,0),death.cumhaz),base1[,1])[,2]
 dtt <- diff(c(0,base1[,1]))
 gt <- (gamma(agam2+nu2)/gamma(agam2))*(1/(betagam+dbase1))^nu2
 lams <-(1/egamma12nu3)*(diff(c(0,base1[,2]))/dtt)*(1/gt)
 cumhaz2 <- cbind(base1[,1],cumsum(dtt*lams))

 cumhaz <- rbind(c(0,0),cumhaz)
 cumhaz2 <- rbind(c(0,0),cumhaz2)
 death.cumhaz <- rbind(c(0,0),death.cumhaz)

## range max of cumhaz and cumhaz2 
  out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz),NULL)
  cumhaz <- out$cum1
  cumhaz2 <- out$cum2
  cumhazd <- out$cum3
  max.time <- tail(cumhaz[,1],1)

### recurrent first time
  tall1 <- rchaz(cumhaz,rr=z1)
  tall2 <- rchaz(cumhaz2,rr=z2)
  tall <- tall1 
  tall$status <- ifelse(tall1$time<tall2$time,tall1$status,2*tall2$status)
  tall$time <- ifelse(tall1$time<tall2$time,tall1$time,tall2$time)
  tall$id <- 1:n
  tall$rr2 <- tall2$rr
### death time simulated
  if (!is.null(death.cumhaz)) {# {{{
	  timed   <- rchaz(cumhazd,rr=zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }# }}}

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  ### setting aside memory 
  tt1 <- tt2 <- tt
  i <- 1; 
  while (any((tt$time<tt$dtime) & (tt$status!=0) & (i < max.recurrent))) {
	  i <- i+1
	  still <- subset(tt,time<dtime & status!=0)
	  nn <- nrow(still)
          tt1 <- rchaz(cumhaz,rr=z1[still$id],entry=still$time)
          tt2 <- rchaz(cumhaz2,rr=z2[still$id],entry=still$time)
	  tt <- tt1
          tt$status <- ifelse(tt1$time<=tt2$time,tt1$status,2*tt2$status)
          tt$time <-   ifelse(tt1$time<=tt2$time,tt1$time,tt2$time)
	  tt$rr2 <- tt2$rr
          ###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt[1:nn,],row.names=NULL)
  }
  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"zs") <- zs

  attr(tall,"gamma.death") <- c(agam1,agam2,betagam,vargamD)
  attr(tall,"gamma.N12") <-   c(agam12,betagam12,vargam12)

  return(tall)
  }# }}}


##' Counts the number of previous events of two types for recurrent events processes
##'
##' Counts the number of previous events of two types for recurrent events processes
##'
##' @param data data-frame
##' @param status name of status 
##' @param id  id 
##' @param types types of the events (code) related to status
##' @param names.count name of Counts, for example Count1 Count2 when types=c(1,2)
##' @param lag if true counts previously observed, and if lag=FALSE counts up to know
##' @param multitype if multitype then count number of types also when types=c(1,2) for example
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##'
##' ######################################################################
##' ### simulating simple model that mimicks data 
##' ### now with two event types and second type has same rate as death rate
##' ######################################################################
##'
##' rr <- simRecurrentII(1000,base1,base4,death.cumhaz=dr)
##' rr <-  count.history(rr)
##' dtable(rr,~"Count*"+status,level=1)
##'
##' @aliases count.historyVar 
##' @export
count.history <- function(data,status="status",id="id",types=1:2,names.count="Count",lag=TRUE,multitype=FALSE)
{# {{{
stat <- data[,status]

clusters <- data[,id]
if (is.numeric(clusters)) {
      clusters <- fast.approx(unique(clusters), clusters) - 1
      max.clust <- max(clusters)
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

data[,"lbnr__id"] <- cumsumstrata(rep(1,nrow(data)),clusters,max.clust+1) 
if (!multitype) {
for (i in types)  {
if (lag==TRUE)
data[,paste(names.count,i,sep="")] <- 
   cumsumidstratasum((stat==i),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
   else 
data[,paste(names.count,i,sep="")] <- 
   cumsumidstratasum((stat==i),rep(0,nrow(data)),1,clusters,max.clust+1)$sum 
}
} else {
if (lag==TRUE)
data[,paste(names.count,types[1],sep="")] <- 
   cumsumidstratasum((stat %in% types),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
   else 
data[,paste(names.count,types[1],sep="")] <- 
   cumsumidstratasum((stat %in% types),rep(0,nrow(data)),1,clusters,max.clust+1)$sum 
}


return(data)
}# }}}

##' @export
count.historyVar <- function(data,var="status",id="id",names.count="Count",lag=TRUE)
{# {{{
vvar <- data[,var]

clusters <- data[,id]
if (is.numeric(clusters)) {
      clusters <- fast.approx(unique(clusters), clusters) - 1
      max.clust <- max(clusters)
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

data[,"lbnr__id"] <- cumsumstrata(rep(1,nrow(data)),clusters,max.clust+1) 
if (lag==TRUE)
data[,names.count] <- 
   cumsumidstratasum(vvar,rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
   else 
data[,names.count] <- 
   cumsumidstratasum(vvar,rep(0,nrow(data)),1,clusters,max.clust+1)$sum 

return(data)
}# }}}


##' Estimation of probability of more that k events for recurrent events process
##'
##' Estimation of probability of more that k events for recurrent events process
##' where there is terminal event, based on this also estimate of variance of recurrent events. The estimator is based on cumulative incidence of exceeding "k" events.
##' In contrast the probability of exceeding k events can also be computed as a 
##' counting process integral, and this is implemented in prob.exceedRecurrent
##'
##' @param data data-frame
##' @param type type of evnent (code) related to status
##' @param status name of status 
##' @param death  name of death indicator 
##' @param start start stop call of Hist() of prodlim 
##' @param stop start stop call of Hist() of prodlim 
##' @param id  id 
##' @param times time at which to get probabilites P(N1(t) >= n)
##' @param exceed n's for which which to compute probabilites P(N1(t) >= n)
##' @param cifmets if true uses cif of mets package rather than prodlim 
##' @param strata to stratify according to variable, only for cifmets=TRUE, when strata is given then only consider the output in the all.cifs
##' @param all.cifs if true then returns list of all fitted objects in cif.exceed 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @references 
##'             Scheike, Eriksson, Tribler (2019) 
##'             The mean, variance and correlation for bivariate recurrent events
##'             with a terminal event,  JRSS-C
##'
##' @examples
##'
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##'
##' cor.mat <- corM <- rbind(c(1.0, 0.6, 0.9), c(0.6, 1.0, 0.5), c(0.9, 0.5, 1.0))
##' rr <- simRecurrentII(1000,base4,cumhaz2=base4,death.cumhaz=dr,cens=2/5000)
##' rr <-  count.history(rr)
##' dtable(rr,~death+status)
##' 
##' oo <- prob.exceedRecurrent(rr,1)
##' bplot(oo)
##' 
##' par(mfrow=c(1,2))
##' with(oo,plot(time,mu,col=2,type="l"))
##' ###
##' with(oo,plot(time,varN,type="l"))
##' 
##' 
##' ### Bivariate probability of exceeding 
##' oo <- prob.exceedBiRecurrent(rr,1,2,exceed1=c(1,5),exceed2=c(1,2))
##' with(oo, matplot(time,pe1e2,type="s"))
##' nc <- ncol(oo$pe1e2)
##' legend("topleft",legend=colnames(oo$pe1e2),lty=1:nc,col=1:nc)
##' 
##' 
##' \donttest{
##' ### do not test to avoid dependence on prodlim 
##' ### now estimation based on cumualative incidence, but do not test to avoid dependence on prodlim 
##' ### library(prodlim)
##' pp <- prob.exceed.recurrent(rr,1,status="status",death="death",start="entry",stop="time",id="id")
##' with(pp, matplot(times,prob,type="s"))
##' ###
##' with(pp, matlines(times,se.lower,type="s"))
##' with(pp, matlines(times,se.upper,type="s"))
##' }
##' @export
##' @aliases prob.exceedRecurrent prob.exceedBiRecurrent prob.exceedRecurrentStrata prob.exceedBiRecurrentStrata summaryTimeobject
prob.exceed.recurrent <- function(data,type,status="status",death="death",
 start="start",stop="stop",id="id",times=NULL,exceed=NULL,cifmets=TRUE,
 strata=NULL,all.cifs=FALSE,...)
{# {{{
### setting up data 
stat <-     data[,status]
dd   <-     data[,death]
tstop <-    data[,stop]
tstart <-   data[,start]
clusters <- data[,id]

if (sum(stat==type)==0) stop("none of type events")
if (!is.null(strata) & !cifmets) stop("strata only for cifmets=TRUE\n")

if (is.numeric(clusters)) {
      clusters <- fast.approx(unique(clusters), clusters) - 1
      max.clust <- max(clusters)
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

 count <- cumsumstrata((stat==type),clusters,max.clust+1)
### count  <- cumsumidstratasum((stat==type),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
 mc <- max(count)+1
 idcount <- clusters*mc + count
 idcount <- cumsumstrata(rep(1,length(idcount)),idcount,mc*(max.clust+1))

if (is.null(times)) times <- sort(unique(tstop[stat==type]))
if (is.null(exceed)) exceed <- sort(unique(count))

if (!cifmets) {
   if (is.null(strata)) form <- as.formula(paste("Hist(entry=",start,",",stop,",statN)~+1",sep=""))
   else form <- as.formula(paste("Hist(entry=",start,",",stop,",statN)~+",strata,sep="")) 
}
else {
   if (is.null(strata)) form <- as.formula(paste("Event(",start,",",stop,",statN)~+1",sep=""))
   else form <- as.formula(paste("Event(",start,",",stop,",statN)~strata(",strata,")",sep=""))
}

cif.exceed <- NULL
if (all.cifs) cif.exceed <- list() 
probs.orig <- se.probs <- probs <- matrix(0,length(times),length(exceed))
se.lower <-  matrix(0,length(times),length(exceed))
se.upper <-  matrix(0,length(times),length(exceed))
i <- 1
for (n1 in exceed[-1]) {# {{{
	i <- i+1
	### first time that get to n1
	keep <- (count<n1 ) | (count==n1 & idcount==1)
	### status, censoring, get to n1, or die
        statN <- rep(0,nrow(data))
	statN[count==n1] <- 1
	statN[dd==1] <- 2
	statN <- statN[keep]
        if (!cifmets) 
	pN1 <-  suppressWarnings(prodlim::prodlim(form,data=data[keep,]))
        else pN1 <-  suppressWarnings(cif(form,data=data[keep,]))
	if (all.cifs) cif.exceed[[i-1]] <- pN1

	if (sum(statN)==0) {
		se.lower[,i] <- se.upper[,i] <- se.probs[,i] <- probs[,i] <- rep(0,length(times)) } else  {

                if (!cifmets) {
			mps  <- summary(pN1,times=times,cause=1)
			mps  <- suppressWarnings(summary(pN1,times=times,cause=1)$table)
			if (is.list(mps)) mps <- mps$"1"
			probs.orig[,i] <- ps <- mps[,5]
			mm <- which.max(ps)
			probs[,i] <- ps
			probs[is.na(ps),i] <- ps[mm]
			se.probs[,i] <- mps[,6]
			se.probs[is.na(ps),i] <- se.probs[mm,i]
			se.lower[,i] <- mps[,7] 
			se.lower[is.na(ps),i] <- se.lower[mm,i]
			se.upper[,i] <- mps[,8]
			se.upper[is.na(ps),i] <- se.upper[mm,i]
	        } else {
			where <- fast.approx(c(0,pN1$times),times,type="left")
	   	        probs[,i] <- c(0,pN1$mu)[where]
			se.probs[,i] <- c(0,pN1$se.mu)[where]
			se.lower[,i] <- probs[,i]-1.96*se.probs[,i] 
			se.upper[,i] <- probs[,i]+1.96*se.probs[,i] 
		}

	}
	if (i==2) { probs[,1]    <- 1-probs[,2]; 
                    se.probs[,1] <- se.probs[,2]; 
                    se.lower[,1] <- 1-se.lower[,2]; 
                    se.upper[,1] <- 1-se.upper[,2]; 
	}
}# }}}

dp <- -t(apply(cbind(probs[,-1],0),1,diff))
meanN <- apply(probs[,-1,drop=FALSE],1,sum)
meanN2 <- apply(t(exceed[-1]^2 * t(dp)),1,sum)
 
colnames(probs) <- c(paste("N=",exceed[1],sep=""),paste("exceed>=",exceed[-1],sep=""))
colnames(se.probs) <- c(paste("N=",exceed[1],sep=""),paste("exceed>=",exceed[-1],sep=""))

return(list(time=times,times=times,prob=probs,se.prob=se.probs,meanN=meanN,probs.orig=probs.orig[,-1],
	    se.lower=se.lower,se.upper=se.upper,meanN2=meanN2,varN=meanN2-meanN^2,exceed=exceed[-1],formula=form,
	    cif.exceed=cif.exceed))
}# }}}

##' @export
prob.exceedRecurrent <- function(data,type,km=TRUE,status="status",death="death",
                start="start",stop="stop",id="id",names.count="Count",...)
{# {{{

formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type,")~cluster(",id,")",sep=""))
###
form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type,")~strata(",names.count,type,")+cluster(",id,")",sep=""))

dr <- phreg(formdr,data=data)
base1   <- phreg(form1,data=data)
base1.2 <- phreg(form1C,data=data)

###cc <- base1$cox.prep
###risk <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
######### risk stratified after count 1
###cc <- base1.2$cox.prep
###risk1 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###pstrata <- risk1/risk
###pstrata[risk1==0] <- 0

### marginal int_0^t G(s) P(N1(t-)==k|D>t) \lambda_{1,N1=k}(s) ds 
### strata og count skal passe sammen
  # {{{
  strat <- dr$strata[dr$jumps]
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- base1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  risktot <- x$S0
  mu <- c(cumsumstrata(St*S0i,xx$strata,xx$nstrata))
###
  x <- base1.2
  xx <- x$cox.prep
  lss <- length(xx$strata)
###  S0i2 <- S0i <- rep(0,lss)
###  S0i[xx$jumps+1] <-  1/x$S0
  riskstrata <- x$S0
  xstrata <- xx$strata
  vals1 <- sort(unique(data[,paste("Count",type,sep="")]))
  valjumps <- vals1[xx$strata+1]
  fk <- (valjumps+1)^2-valjumps^2
###  EN2 <- c(cumsumstrata(fk*St*pstrata*S0i,rep(0,lss),1))
  EN2 <- c(cumsumstrata(fk*St*S0i,rep(0,lss),1))
  pcumhaz <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
# }}}
  EN2     <- EN2[xx$jumps+1]
  cumhaz <- pcumhaz[xx$jumps+1,]
  mu     <- mu[xx$jumps+1]
  pstrata <- riskstrata/risktot

  exceed.name <- paste("Exceed>=",vals1+1,sep="")

  out=list(cumhaz=cumhaz,time=cumhaz[,1],varN=EN2-mu^2,mu=mu,
	nstrata=base1.2$nstrata,strata=base1.2$strata[xx$jumps+1],
	jumps=1:nrow(cumhaz),riskstrata=pstrata,risktot=risktot,
	strat.cox.name=base1.2$strata.name,
	strat.cox.level=base1.2$strata.level,exceed=vals1+1,
        strata.name=exceed.name,strata.level=exceed.name)

### use recurrentMarginal estimator til dette via strata i base1 
### strata og count skal passe sammen
### see beregning via recurrent marginal function
###  base1$cox.prep$strata <- base1.2$cox.prep$strata
###  base1$cox.prep$nstrata <- base1.2$cox.prep$nstrata
###  base1$nstrata <- base1.2$cox.prep$nstrata
###  base1$strata <- base1.2$strata
###  base1$strata.name <- base1.2$strata.name
###  base1$strata.level <- base1.2$strata.level

###  mm <- recurrentMarginal(base1,dr,km=km,...)
###  out=c(mm,list(varN=EN2-mu^2))

  return(out)
}# }}}

##' @export
prob.exceedRecurrentStrata <- function(data,type,km=TRUE,status="status",death="death",
                start="start",stop="stop",id="id",names.count="Count",strata=NULL,...)
{# {{{

if (is.null(strata)) {
formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))
## bring count as covariate to use later and get sorted as data 
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type,")~",names.count,type,"+cluster(",id,")",sep=""))
} else {
formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~strata(",strata,")+cluster(",id,")",sep=""))
## bring count as covariate to use later and get sorted as data 
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type,")~",names.count,type,"+strata(",strata,")+cluster(",id,")",sep=""))
}

dr      <- phreg(formdr,data=data,no.opt=TRUE,no.var=1)
base1   <- phreg(form1,data=data,no.opt=TRUE,no.var=1)

### marginal int_0^t G(s) P(N1(t-)==k|D>t) \lambda_{1,N1=k}(s) ds 
### strata og count skal passe sammen
  # {{{
  strat <- dr$strata[dr$jumps]
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- base1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  risktot <- x$S0
  mu <- c(cumsumstrata(St*S0i,xx$strata,xx$nstrata))
  ###
  vals1 <- xx$X[,1]
  fk <- (vals1+1)^2-vals1^2
  EN2 <- c(cumsumstrata(fk*St*S0i,xx$strata,xx$nstrata))
  inc <- St[xx$jumps+1]*S0i[xx$jumps+1]

  ## new-strata names  
  valjump <- vals1[xx$jumps+1]
  xxs <- xx$strata[xx$jumps+1]
  newstrata <- mystrata(list(id=xxs,exceed=valjump))
  nnn <- !duplicated(newstrata$sindex)
  nnstrata <- attr(newstrata,"nlevel")
  newstrata <- newstrata$sindex
  exceed <- valjump[nnn]+1
  if (!is.null(strata)) 
	  exceed.levels <- paste(base1$strata.level[xxs[nnn]+1],
				 paste("Exceed",exceed,sep=">="),sep="-") 
  else exceed.levels <- paste("Exceed",exceed,sep=">=")

  newstrata <- as.numeric(strata(xx$strata[xx$jumps+1],valjump))-1
  nnstrata <- length(unique(newstrata))
  pcumhaz <- cbind(x$jumptimes,cumsumstrata(inc,newstrata,nnstrata))
# }}}

  EN2     <- EN2[xx$jumps+1]
  mu     <- mu[xx$jumps+1]
  cumhaz <- pcumhaz
  pstrata <- NULL 

  out=list(cumhaz=cumhaz,time=cumhaz[,1],varN=EN2-mu^2,mu=mu,
	nstrata=nnstrata,strata=newstrata,
	strat.cox.name=base1$strata.name,
	strat.cox.level=base1$strata.level,exceed=exceed,
	jumps=1:nrow(cumhaz),riskstrata=pstrata,risktot=risktot,
        strata.name="",strata.level=exceed.levels)


  return(out)
}# }}}

##' @export
prob.exceedBiRecurrent <- function(data,type1,type2,km=TRUE,status="status",death="death",
      start="start",stop="stop",id="id",names.count="Count",exceed1=NULL,exceed2=NULL)
{# {{{

formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~cluster(",id,")",sep=""))
form2 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~cluster(",id,")",sep=""))
###
###form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))
###form2C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~
###  strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))

form2Ccc <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~
  ",names.count,type1,"+",names.count,type2,"+","
  strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))
form1Ccc <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~
  ",names.count,type1,"+",names.count,type2,"+","
  strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))

### stratified and with counts in covariate matrix 
bb2.12 <- phreg(form2Ccc,data=data,no.opt=TRUE,no.var=1)
bb1.12 <- phreg(form1Ccc,data=data,no.opt=TRUE,no.var=1)

dr <- phreg(formdr,data=data)
base1   <- phreg(form1,data=data,no.var=1)
base2   <- phreg(form2,data=data,no.var=1)

cc <- base1$cox.prep
risk1 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###### risk stratified after count 1 og count2
cc <- bb1.12$cox.prep
risk1.12 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
pstrata1 <- risk1.12/risk1
pstrata1[1] <- 0

cc <- base2$cox.prep
risk2 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###### risk stratified after count 1 og count2
cc <- bb2.12$cox.prep
risk2.12 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
pstrata2 <- risk2.12/risk2
pstrata2[1] <- 0

### marginal int_0^t G(s) P(N1(t-)==k|D>t) \lambda_{1,N1=k}(s) ds 
### strata og count skal passe sammen

  # {{{
  strat <- dr$strata[dr$jumps]
  Gt <- exp(-dr$cumhaz[,2])
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- base1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  mu <- c(cumsumstrata(St*S0i,rep(0,lss),1))
###

  x <- bb1.12
  xx1 <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx1$jumps+1] <-  1/x$S0
  dcumhaz1 <- cbind(xx1$time,pstrata1*St*S0i)
###              cumsumstrata(pstrata1*St*S0i,xx1$strata,xx1$nstrata))
  x <- bb2.12
  xx2 <- x$cox.prep
  lss <- length(xx2$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx2$jumps+1] <-  1/x$S0
  dcumhaz2 <- cbind(xx2$time,pstrata2*St*S0i)
###   cumsumstrata(pstrata2*St*S0i,xx2$strata,xx2$nstrata))

  n1 <- length(xx1$jumps)
  n2 <- length(xx2$jumps)
  ojumps <- order(c(xx1$jumps,xx2$jumps))
  jumps <- sort(c(xx1$jumps,xx2$jumps))
  times <- xx1$time[jumps+1]

  dcumhaz1 <-  dcumhaz1[jumps+1,] 
  dcumhaz2 <-  dcumhaz2[jumps+1,] 
  x1 <- xx1$X[jumps+1,]
  x2 <- xx2$X[jumps+1,]

###  dcumhaz1 <- dcumhaz1[xx1$jumps+1,]
###  dcumhaz2 <- dcumhaz2[xx2$jumps+1,]
###  x1 <- xx1$X[xx1$jumps+1,]
###  x2 <- xx2$X[xx2$jumps+1,]
# }}}

  if (is.null(exceed1)) exceed1 <- 1:max(x1[,1])
  if (is.null(exceed2)) exceed2 <- 1:max(x1[,2])

  pe1e2 <- matrix(0,n1+n2,length(exceed1)*length(exceed2))
  m <- 0; nn <- c()
  for (i in exceed1) 
  for (j in exceed2)  {
	  m <- m+1
	  strat1 <- (x1[,2]>=j)*(x1[,1]==(i-1))
	  strat2 <- (x2[,1]>=i)*(x2[,2]==(j-1))
	  pe1e2[,m] <- cumsum(strat1*dcumhaz1[,2]) + cumsum(strat2*dcumhaz2[,2])
	  nn <- c(nn,paste("N_1(t)>=",i,",N_2(t)>=",j,sep="")) 
  }

  colnames(pe1e2) <- nn

  out=list(time=times,pe1e2=pe1e2,x1=x1,x2=x2,
	   nstrata=base1$nstrata,
	   strata.name=base1$strata.name,strata.level=base1$strata.levels)
  class(out) <- "BiRecurrent"

  return(out)
}# }}}

##' @export
plot.BiRecurrent <- function(x,stratas=NULL,add=FALSE,...)
{# {{{

   strat <- x$strata
   ## all strata
   if (is.null(stratas)) stratas <- 0:(x$nstrata-1) 

   for (s in stratas)  {
   if (add==FALSE)
    with(x, matplot(time[strata==s],pe1e2[strata==s,],type="s",...))
    else 
    with(x, matlines(time[strata==s],pe1e2[strata==s,],type="s",...))
   nc <- ncol(x$pe1e2)
   legend("topleft",colnames(x$pe1e2),lty=1:nc,col=1:nc)
   }
}# }}}


##' @export
prob.exceedBiRecurrentStrata <- function(data,type1,type2,km=TRUE,status="status",
  death="death",start="start",stop="stop",id="id",names.count="Count",
  strata=NULL,twinstrata=FALSE,exceed1=NULL,exceed2=NULL)
{# {{{

if (is.null(strata)) {
formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))
## use count and status as covariates to use later and get sorted as data 
form <- as.formula(paste("Surv(",start,",",stop,",",status,"!=0)~",status,"+",names.count,type1,"+",names.count,type2,"+cluster(",id,")",sep=""))
} else {
formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~strata(",strata,")+cluster(",id,")",sep=""))
## use count and status as covariates to use later and get sorted as data 
form <- as.formula(paste("Surv(",start,",",stop,",",status,"!=0)~",status,"+",
     names.count,type1,"+",names.count,type2,"+","strata(",strata,")+cluster(",id,")",sep=""))
	if (twinstrata) { ## to allow different strata for the two twins
	form <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~",status,"+",
	  names.count,type1,"+",names.count,type2,"+","strata(",strata,type2,")+cluster(",id,")",sep=""))
	}
}

dr     <- phreg(formdr,data=data)
base   <- phreg(form,data=data,no.opt=TRUE)

  # {{{
  strat <- dr$strata[dr$jumps]
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))

  xx <- base$cox.prep
  lss <- length(xx$strata)
  xxjump <- xx$jumps+1
  S0i <-  1/base$S0
  times <- base$jumptimes

  ## jumps for N1 and N2, sorted 
  xxstrata <- xx$strata[xxjump]
  St <- St[xxjump]
  statusj <- xx$X[xxjump,1]
  count1 <- xx$X[xxjump,2]
  count2 <- xx$X[xxjump,3]

  if (is.null(exceed1)) { exceed1 <- sort(unique(count1))+1 }
  if (is.null(exceed2)) { exceed2 <- sort(unique(count2))+1 }
  n <- length(xxjump)

  m <- 1; nn <- c(); 
  pe1e2 <- matrix(0,n,length(exceed1)*length(exceed2))
  for (i in exceed1) 
  for (j in exceed2)  {
	  escape <- (count2>=j)*(count1==(i-1))*(statusj==1)+(count1>=i)*(count2==(j-1))*(statusj==2)
	  pe1e2[,m] <- cumsumstrata(escape*St*S0i,xxstrata,xx$nstrata)  
	  nn <- c(nn,paste("N_1(t)>=",i,",N_2(t)>=",j,sep="")) 
	  m <- m+1
  }

  colnames(pe1e2) <- nn
# }}}


  out=list(time=times,pe1e2=pe1e2,strata=xxstrata,nstrata=xx$nstrata,
	   cumhazard=cbind(times,pe1e2), jumps=1:length(times), nstrata=xx$nstrata,
	   strata.name=base$strata.name,strata.level=base$strata.levels)
  class(out) <- "BiRecurrent"

  return(out)
}# }}}


##' Estimation of covariance for bivariate recurrent events with terminal event
##'
##' Estimation of probability of more that k events for recurrent events process
##' where there is terminal event 
##'
##' @param data data-frame
##' @param type1 type of first event (code) related to status
##' @param type2 type of second event (code) related to status
##' @param status name of status 
##' @param death  name of death indicator 
##' @param start start stop call of Hist() of prodlim 
##' @param stop start stop call of Hist() of prodlim 
##' @param id  id 
##' @param names.count name of count for number of previous event of different types, here generated by count.history()
##' @author Thomas Scheike
##' @references 
##'             Scheike, Eriksson, Tribler (2019) 
##'             The mean, variance and correlation for bivariate recurrent events
##'             with a terminal event,  JRSS-C
##'
##' @examples
##'
##' ########################################
##' ## getting some data to work on 
##' ########################################
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##' rr <- simRecurrentII(1000,base1,cumhaz2=base4,death.cumhaz=dr)
##' rr <- count.history(rr)
##' rr$strata <- 1
##' dtable(rr,~death+status)
##' 
##' covrp <- covarianceRecurrent(rr,1,2,status="status",death="death",
##'                         start="entry",stop="time",id="id",names.count="Count")
##' par(mfrow=c(1,3)) 
##' plot(covrp)
##' 
##' ### with strata, each strata in matrix column, provides basis for fast Bootstrap
##' covrpS <- covarianceRecurrentS(rr,1,2,status="status",death="death",
##'         start="entry",stop="time",strata="strata",id="id",names.count="Count")
##' 
##' @aliases plot.covariace.recurrent  covarianceRecurrentS Bootcovariancerecurrence BootcovariancerecurrenceS 
##' @export
covarianceRecurrent <- function(data,type1,type2,status="status",death="death",
                   	 start="start",stop="stop",id="id",names.count="Count")
{# {{{

formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~cluster(",id,")",sep=""))
form2 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~cluster(",id,")",sep=""))
form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~strata(",names.count,type2,")+cluster(",id,")",sep=""))
form2C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~strata(",names.count,type1,")+cluster(",id,")",sep=""))

dr <- phreg(formdr,data=data)
base1   <- phreg(form1,data=data)
base1.2 <- phreg(form1C,data=data)
base2   <- phreg(form2,data=data)
base2.1 <- phreg(form2C,data=data)

marginal.mean1 <- recmarg(base1,dr)
marginal.mean2 <- recmarg(base2,dr)

cc <- base2$cox.prep
risk <- c(revcumsumstrata(cc$sign,cc$strata,cc$nstrata))
###### risk stratified after count 1
cc <- base2.1$cox.prep
risk1 <- c(revcumsumstrata(cc$sign,cc$strata,cc$nstrata))
ssshed1 <- risk1/risk
ssshed1[is.na(ssshed1)] <- 1
sshed1   <- list(cumhaz=cbind(cc$time,ssshed1),
	      strata=cc$strata,nstrata=cc$nstrata,
	      jumps=1:length(cc$time),
	      strata.name=paste("prob",type1,sep=""),
	      strata.level=base2.1$strata.level)
riskstrata <- .Call("riskstrataR",cc$sign,cc$strata,cc$nstrata)$risk
nrisk <- apply(riskstrata,2,revcumsumstrata,rep(0,nrow(riskstrata)),1)
ntot <- apply(nrisk,1,sum)
vals1 <- sort(unique(data[,paste("Count",type1,sep="")]))
mean1risk <- apply(t(nrisk)*vals1,2,sum)/ntot
mean1risk[is.na(mean1risk)] <- 0


cc <- base1$cox.prep
S0 <- rep(0,length(cc$strata))
risk <- c(revcumsumstrata(cc$sign,cc$strata,cc$nstrata))
###
cc <- base1.2$cox.prep
S0 <- rep(0,length(cc$strata))
risk2 <- c(revcumsumstrata(cc$sign,cc$strata,cc$nstrata))
ssshed2 <- risk2/risk
ssshed2[is.na(ssshed2)] <- 1
###
sshed2  <- list(cumhaz=cbind(cc$time,ssshed2),
	      strata=cc$strata,nstrata=cc$nstrata,
              jumps=1:length(cc$time),
	      strata.name=paste("prob",type2,sep=""),
	      strata.level=base1.2$strata.level)
###
riskstrata <- .Call("riskstrataR",cc$sign,cc$strata,cc$nstrata)$risk
nrisk <- apply(riskstrata,2,revcumsumstrata,rep(0,nrow(riskstrata)),1)
ntot <- apply(nrisk,1,sum)
vals2 <- sort(unique(data[,paste("Count",type2,sep="")]))
mean2risk <- apply(t(nrisk)*vals2,2,sum)/ntot
mean2risk[is.na(mean2risk)] <- 0

mu1 <- cpred(rbind(c(0,0),marginal.mean1$cumhaz),cc$time)[,2]
mu2 <- cpred(rbind(c(0,0),marginal.mean2$cumhaz),cc$time)[,2]


out <- list(based=dr,base1=base1,base2=base2,
	    base1.2=base1.2,base2.1=base2.1,
	    marginal.mean1=marginal.mean1,marginal.mean2=marginal.mean2,
	    prob1=sshed1,prob2=sshed2,
	    mean1risk=mean1risk,mean2risk=mean2risk)

### marginal sum_k int_0^t G(s) k P(N1(t-)==k|D>t) \lambda_{2,N1=k}(s) ds 
### strata og count skal passe sammen
  # {{{
  strat <- dr$strata[dr$jumps]
  Gt <- exp(-dr$cumhaz[,2])
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
###
  x <- base1.2
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  xstrata <- xx$strata
  cumhazDR <- cbind(xx$time,cumsumstrata(vals2[xstrata+1]*St*ssshed2*S0i,rep(0,lss),1))
  x <- base1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazIDR <- cbind(xx$time,cumsumstrata(St*mean2risk*S0i,rep(0,lss),1))
  mu1.i <- cumhazIDR[,2]
  mu1.2 <- cumhazDR[,2]
###
  x <- base2.1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazDR <- cbind(xx$time,cumsumstrata(vals1[xx$strata+1]*St*ssshed1*S0i,rep(0,lss),1))
###
  x <- base2
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazIDR <- cbind(xx$time,cumsumstrata(St*mean1risk*S0i,rep(0,lss),1))
  mu2.i <- cumhazIDR[,2]
  mu2.1 <- cumhazDR[,2]
# }}}

  out=c(out,list(EN1N2= mu1.2+mu2.1,mu2.1=mu2.1,mu1.2=mu1.2, 
		 mu2.i=mu2.i,mu1.i=mu1.i,
		 EIN1N2=mu2.i+mu1.i,EN1EN2=mu1*mu2,time=cc$time))

  class(out) <- "covariance.recurrent"

  return(out)
}# }}}

##' @export
plot.covariance.recurrent <- function(x,main="Covariance",these=1:3,...) 
{# {{{

legend <- NULL # to avoid R-check 

if ( 1 %in% these) {
	nna <- (!is.na(x$mu1.2)) & (!is.na(x$mu1.i))
	mu1.2n <- x$mu1.2[nna]
	mu1.in <- x$mu1.i[nna]
	time <-  x$time[nna]
	plot(time,mu1.2n,type="l",ylim=range(c(mu1.2n,mu1.in)),...) 
	lines(time,mu1.in,col=2) 
	legend("topleft",c(expression(integral(N[2](s)*dN[1](s),0,t)),"independence"),lty=1,col=1:2) 
	title(main=main)
}
###
if (2 %in% these) {
	nna <- (!is.na(x$mu2.1)) & (!is.na(x$mu2.i))
	mu2.1n <- x$mu2.1[nna]
	mu2.in <- x$mu1.i[nna]
	time <- x$time[nna]
	plot(time,mu2.1n,type="l",ylim=range(c(mu2.1n,mu2.in)),...) 
	lines(time,mu2.in,col=2) 
	legend("topleft",c(expression(integral(N[1](s)*dN[2](s),0,t)),"independence"),lty=1,col=1:2) 
	title(main=main)
}
###
if (3 %in% these) {
	nna <- (!is.na(x$EN1N2)) & (!is.na(x$EIN1N2)) & (!is.na(x$EN1EN2))
	EN1N2n <- x$EN1N2[nna]
	EIN1N2n <- x$EIN1N2[nna]
	EN1EN2n <- x$EN1EN2[nna]
	time <- x$time[nna]
	plot(time,EN1N2n,type="l",lwd=2,ylim=range(c(EN1N2n,EN1EN2n,EIN1N2n)),...) 
	lines(time,EN1EN2n,col=2,lwd=2) 
	lines(time,EIN1N2n,col=3,lwd=2) 
	legend("topleft",c("E(N1N2)", "E(N1) E(N2) ", "E_I(N1 N2)-independence"),lty=1,col=1:3)
	title(main=main)
}

} # }}}

meanRisk <- function(base1,base1.2)
{# {{{
cc <- base1.2$cox.prep
S0 <- rep(0,length(cc$strata))
mid <- max(cc$id)+1
risk2 <- revcumsumidstratasum(cc$sign,cc$id,mid,cc$strata,cc$nstrata)$sumidstrata

means <- .Call("meanriskR",cc$sign,cc$id,mid,cc$strata,cc$nstrata)
mean2risk <- means$meanrisk
mean2risk[is.na(mean2risk)] <- 0
risk <- means$risk
ssshed2 <- risk2/risk
ssshed2[is.na(ssshed2)] <- 0
vals2 <- unique(cc$id)

means2  <- list(cumhaz=cbind(cc$time,mean2risk),
	      strata=cc$strata,nstrata=cc$nstrata,
              jumps=1:length(cc$time),
	      strata.name="meansrisk",
	      strata.level=base1.2$strata.level)
sshed2  <- list(cumhaz=cbind(cc$time,ssshed2),
	      strata=cc$id,nstrata=mid,
              jumps=1:length(cc$time),
	      strata.name="prob",
	      strata.level=paste(vals2),real.strata=cc$strata)

return(list(meanrisk=means2,vals=vals2,probs=sshed2,jumps=cc$jumps+1))
}# }}}

intN2dN1 <- function(dr,base1,base1.2,pm)
{# {{{
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
###
  x <- base1.2
  xx <- x$cox.prep
  xstrata <- xx$id
  jumps <- xx$jumps+1
  mid <- max(xx$id)+1
  ## risk after both id~Count and strata 
  risk2 <- revcumsumidstratasum(xx$sign,xx$id,mid,xx$strata,xx$nstrata)$sumidstrata
  S0i <-  1/risk2[jumps]
  vals <- xx$id[jumps]
  St <- St[jumps]
  probs <- pm$probs$cumhaz[jumps,2]
  cumhazDR <- cbind(xx$time[jumps],cumsumstrata(St*vals*probs*S0i,xx$strata[jumps],xx$nstrata))
  x <- base1
  xx <- x$cox.prep
  S0i <-  c(1/x$S0)
  meanrisk <- pm$meanrisk$cumhaz[jumps,2]
  cumhazIDR <- cbind(xx$time[jumps],cumsumstrata(St*meanrisk*S0i,xx$strata[jumps],xx$nstrata))
  mu1.i <- cumhazIDR[,2]
  mu1.2 <- cumhazDR[,2]
  return(list(cumhaz=cumhazDR,cumhazI=cumhazIDR,mu1.i=mu1.i,mu1.2=mu1.2,
	      time=cumhazDR[,1],
	      strata=xx$strata[jumps],nstrata=xx$nstrata,jumps=1:length(mu1.i),
	      strata.name="intN2dN1",strata.level=x$strata.level))
}# }}}

recmarg2 <- function(recurrent,death,...)
{# {{{
  xr <- recurrent
  dr <- death 

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <- 1/x$S0
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
  ###
  x <- xr
  xx <- x$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i <-  1/x$S0
  jumps <- xx$jumps+1
  cumhazDR <- cbind(xx$time[jumps],cumsumstrata(St[jumps]*S0i,xx$strata[jumps],xx$nstrata))
  mu <- cumhazDR[,2]
# }}}

 varrs <- data.frame(mu=mu,time=cumhazDR[,1],strata=xr$strata[jumps],St=St[jumps])
 out <- list(mu=varrs$mu,times=varrs$time,St=varrs$St,cumhaz=cumhazDR,
     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
     strata.name=xr$strata.name)
 return(out)
}# }}}

##' @export
covarianceRecurrentS <- function(data,type1,type2,times=NULL,status="status",death="death",
                   	 start="start",stop="stop",id="id",names.count="Count",
			 strata="NULL",plot=0,output="matrix")
{# {{{


if (is.null(times)) times <- seq(0,max(data[,stop]),length=100)

### passing strata as id to be able to use for stratified calculations 
if (is.null(strata)) stop("must give strata, for example one strata\n"); 
## uses Counts1 as cluster to pass to risk set calculations 


formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~strata(",strata,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~strata(",strata,")",sep=""))
form2 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~strata(",strata,")",sep=""))
###
form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~cluster(",names.count,type2,")+strata(",strata,")",sep=""))
form2C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~cluster(",names.count,type1,")+strata(",strata,")",sep=""))


dr <- phreg(formdr,data=data)
###basehazplot.phreg(dr)
###
base1   <- phreg(form1,data=data)
base1.2 <- phreg(form1C,data=data)
###
base2   <- phreg(form2,data=data)
base2.1 <- phreg(form2C,data=data)

rm1 <- recmarg2(base1,dr)
rm2 <- recmarg2(base2,dr)
if (plot==1) {
basehazplot.phreg(rm1)
basehazplot.phreg(rm2)
}


pm1 <- meanRisk(base2,base2.1)
pm2 <- meanRisk(base1,base1.2)
if (plot==1) {
basehazplot.phreg(pm1$meanrisk)
basehazplot.phreg(pm2$meanrisk)
}


### marginal sum_k int_0^t G(s) k P(N1(t-)==k|D>t) \lambda_{2,N1=k}(s) ds 
###          sum_k int_0^t G(s) E(N1(t-)==k|D>t)   \lambda_{2}(s) ds 
iN2dN1 <-  intN2dN1(dr,base1,base1.2,pm2)
iN1dN2 <-  intN2dN1(dr,base2,base2.1,pm1)

###print("hej")
if (plot==1) {
par(mfrow=c(2,2))
basehazplot.phreg(iN2dN1)
plot(iN2dN1$time,iN2dN1$cumhazI[,2],type="l")
basehazplot.phreg(iN1dN2)
}

#### writing output in matrix form for each strata for the times
mu1g <- matrix(0,length(times),rm1$nstrata)
mu2g <- matrix(0,length(times),rm1$nstrata)
mu1.2 <- matrix(0,length(times),rm1$nstrata)
mu2.1 <- matrix(0,length(times),rm1$nstrata)
mu1.i <- matrix(0,length(times),rm1$nstrata)
mu2.i <- matrix(0,length(times),rm1$nstrata)
mu1 <- matrix(0,length(times),rm1$nstrata)
mu2 <- matrix(0,length(times),rm1$nstrata)

if (output=="matrix") {
all <- c()
i <- 1
### going through strata 
for (i in 1:rm1$nstrata) {
	j <- i-1
	mu1[,i]   <- cpred(rm1$cumhaz[rm1$strata==j,],times)[,2]
	mu2[,i]   <- cpred(rm2$cumhaz[rm2$strata==j,],times)[,2]
	mu1.2[,i] <- cpred( iN2dN1$cumhaz[rm1$strata==j,],times)[,2]
	mu1.i[,i] <- cpred(iN2dN1$cumhazI[rm1$strata==j,],times)[,2]
	mu2.1[,i] <- cpred( iN1dN2$cumhaz[rm2$strata==j,],times)[,2]
	mu2.i[,i] <- cpred(iN1dN2$cumhazI[rm2$strata==j,],times)[,2]
}
mu1mu2 <- mu1*mu2

out=c(list(EN1N2=mu1.2+mu2.1, EIN1N2=mu1.i+mu2.i, EN1EN2=mu1mu2,
	   mu1.2=mu1.2, mu1.i=mu1.i, mu2.1=mu2.1, mu2.i=mu2.i,
	   mu1=mu1,mu2=mu2, nstrata=rm1$nstrata, time=times))

} else out <- list(iN2dN1=iN2dN1,iN1dN2=iN1dN2,rm1=rm1,rm2=rm2,pm1=pm1,pm2=pm2)


  return(out)
}# }}} 

##' @export
BootcovariancerecurrenceS <- function(data,type1,type2,status="status",death="death",
	 start="start",stop="stop",id="id",names.count="Count",times=NULL,K=100)
{# {{{


  if (is.null(times)) times <- seq(0,max(data[,stop]),length=100)
  mu1.2 <- matrix(0,length(times),K)
  mu2.1 <- matrix(0,length(times),K)
  mu1.i <- matrix(0,length(times),K)
  mu2.i <- matrix(0,length(times),K)
  mupi <- matrix(0,length(times),K)
  mupg <- matrix(0,length(times),K)
  mu1mu2 <- matrix(0,length(times),K)
  n <- length(unique(data[,id]))

  formid <- as.formula(paste("~",id))
  rrb <- blocksample(data, size = n*K, formid)
  rrb$strata <- floor((rrb[,id]-0.01)/n)
  rrb$jump <- (rrb[,status] %in% c(type1,type2)) | (rrb[,death]==1)
  rrb <- tie.breaker(rrb,status="jump",start=start,stop=stop,id=id)

  mm <- covarianceRecurrentS(rrb,type1,type2,status=status,death=death,
 	               start=start,stop=stop,id=id,names.count=names.count,
		       strata="strata",times=times)

  mm <- c(mm,list(se.mui=apply(mm$EIN1N2,1,sd),se.mug=apply(mm$EN1N2,1,sd)))
  return(mm)

}# }}}

##' @export
Bootcovariancerecurrence <- function(data,type1,type2,status="status",death="death",
	 start="start",stop="stop",id="id",names.count="Count",times=NULL,K=100)
{# {{{

  strata <- NULL # to avoid R-check 

  if (is.null(times)) times <- seq(0,max(data[,stop]),length=100)
  mu1.2 <- matrix(0,length(times),K)
  mu2.1 <- matrix(0,length(times),K)
  mu1.i <- matrix(0,length(times),K)
  mu2.i <- matrix(0,length(times),K)
  mupi <- matrix(0,length(times),K)
  mupg <- matrix(0,length(times),K)
  mu1mu2 <- matrix(0,length(times),K)
  n <- length(unique(data[,id]))

  formid <- as.formula(paste("~",id))
  rrb <- blocksample(data, size = n*K, formid)
  rrb$strata <- floor((rrb[,id]-0.01)/n)
## rrb$jump <- (rrb[,status]!=0) | (rrb[,death]==1)
  rrb$jump <- (rrb[,status] %in% c(type1,type2)) | (rrb[,death]==1)
  rrb <- tie.breaker(rrb,status="jump",start=start,stop=stop,id=id)

  for (i in 1:K)
  {
     rrbs <- subset(rrb,strata==i-1)
     errb <- covarianceRecurrent(rrbs,type1,type2,status=status,death=death,
 	                          start=start,stop=stop,id=id,names.count=names.count)
     all <- cpred(cbind(errb$time,errb$EIN1N2,errb$EN1N2,errb$EN1EN2,
			errb$mu1.2,errb$mu2.1,errb$mu1.i,errb$mu2.i),times)
     mupi[,i]   <- all[,2]
     mupg[,i]   <- all[,3]
     mu1mu2[,i] <- all[,4]
     mu1.2[,i]  <- all[,5]; 
     mu2.1[,i]  <- all[,6]; 
     mu1.i[,i]  <- all[,7]; 
     mu2.i[,i]  <- all[,8]
  }

return(list(mupi=mupi,mupg=mupg,mu1mu2=mu1mu2,time=times,
	    EN1N2=mupg,EIN1N2=mupi,EN1EN=mu1mu2,
	    mu1.2=mu1.2,mu1.i=mu1.i,mu2.1=mu2.1,mu2.i=mu1.i,
	    mup=apply(mupi,1,mean),mug=apply(mupg,1,mean),
	    dmupg=apply(mupg-mupi,1,mean),mmu1mu2=apply(mu1mu2,1,mean), 
	    se.mui=apply(mupi,1,sd),se.mug=apply(mupg,1,sd)
	    ))
}# }}}

iidCovarianceRecurrent <-  function (rec1,death,xrS,xr,means)
{# {{{
    ### makes iid decompition for covariance under independence between events
    axr <- rec1
    adr <- death
    St <- exp(-adr$cum[, 2])
    timesr <- axr$cum[, 1]
    timesd <- adr$cum[, 1]
    times <- c(timesr[-1], timesd[-1])
    or <- order(times)
    times <- times[or]
    meano <- cbind(means$time,means$mean2risk)
###    imeano <- sindex.prodlim(means$time, times, strict = FALSE)
    imeano <- fast.approx(means$time, times, type="left")
    meano <- meano[imeano,2]
    keepr <- order(or)[1:length(timesr[-1])]
###    rid <- sindex.prodlim(timesd, times, strict = FALSE)
    rid <- fast.approx(timesd, times, type="left")
###    rir <- sindex.prodlim(timesr, times, strict = FALSE)
    rir <- fast.approx(timesr, times, type="left")
    Stt <- St[rid]
    ariid <- axr$cum[rir, 2]
    mu <- cumsum(meano * Stt * diff(c(0, ariid)))
    muS <- cumsum( Stt * diff(c(0, ariid)))
    nc <- length(axr$B.iid)
    muiid <- matrix(0, length(times), nc)

   cc <- xrS$cox.prep
   rrs <- .Call("riskstrataR",cc$sign*cc$strata,cc$id,max(cc$id)+1)$risk
   rr <- .Call("riskstrataR",cc$sign,cc$id,max(cc$id)+1)$risk
   ntot <- revcumsumstrata(cc$sign,rep(0,nrow(rr)),1)
   rr  <- apply(rr,2,revcumsumstrata,rep(0,nrow(rr)),1)
   rrs <- apply(rrs,2,revcumsumstrata,rep(0,nrow(rr)),1)

###   xrid <- sindex.prodlim(cc$time, times, strict = FALSE)
   xrid <- fast.approx(cc$time, times, type="left")
   rr <- rr[xrid,]
   rrs <- rrs[xrid,]
   ntot <- ntot[xrid]
   rrs <- apply(Stt* diff(c(0,ariid))*rrs,2,cumsum)
   rrcum <- apply(Stt*meano*diff(c(0,ariid))*rr,2,cumsum)
   miid <- (rrs-rrcum)/ntot

    for (i in 1:nc) {
        mriid <- axr$B.iid[[i]]
        mdiid <- adr$B.iid[[i]]
        mriid <- mriid[rir]
        mdiid <- mdiid[rid]
        dmridd <- diff(c(0, mriid))
        dmdidd <- diff(c(0, mdiid))
        muiid[, i] <- cumsum(Stt * meano* dmridd) - mu * cumsum(dmdidd) + cumsum(mu * dmdidd) +  miid[,i]
    }
    var1 <- apply(muiid^2, 1, sum)
    se.mu <- var1[keepr]^0.5
    mu = mu[keepr]
    timeso <- times
    times <- times[keepr]

    out = list(iidtimes=timeso,muiid=muiid,times=times, 
        mu = mu, var.mu = var1[keepr], se.mu = se.mu, St = St, Stt = Stt[keepr], 
	cumhaz=cbind(times,mu),se.cumhaz=cbind(times,se.mu), 
	nstrata=1,strata=rep(0,length(mu)),jumps=1:length(mu)) 

}# }}} 

simMarginalMeanCox <- function(n,cens=3/5000,k1=0.1,k2=0,bin=1,Lam1=NULL,Lam2=NULL,LamD=NULL,beta1=rep(0,2),betad=rep(0,2),betac=rep(0,2),X=NULL,...)
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

 rr <- simRecurrentCox(n,scalecumhaz(Lam1,k1),cumhaz2=scalecumhaz(Lam1,k2),
		       death.cumhaz=LamD,X=X,cens=cens,r1=r1,rd=rd,rc=rc,...)

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

