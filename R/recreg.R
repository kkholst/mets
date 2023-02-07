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
          case.weights=NULL,cause=1,death.code=2,cens.code=0,Gc=NULL,cens.model=~+1,augmentation=NULL,
	  cox.prep=FALSE,wcomp=NULL,augment.model=NULL,ftime.augment=NULL,...) { # {{{
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
    if ((!is.null(augment.model)) & (length(other)>1)  & (length(whereC)>0)) {## {{{

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
###   if (!is.null(augmentation)) varmc <- var.augment.times

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
###	     formula = formula, formC = formC,
###        exit = exit, cens.weights = cens.weights, cens.strata = cens.strata,
###        cens.nstrata = cens.nstrata, model.frame = m, n = length(exit),
###        nevent = nevent, ncluster = nid, Y = Y))

    if (se) {# {{{

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

simGLcox <- function(n,base1,drcumhaz,var.z=0,r1,rd,rc,fz,model=c("frailty","twostage"),cens=NULL)
{# {{{
## setting up baselines for simulations 
base1 <- predictCumhaz(rbind(0,as.matrix(base1)),1:round(tail(base1[,1],1)) )
cumD <- predictCumhaz(rbind(0,as.matrix(drcumhaz)),base1[,1])
St <- exp(-cumD[,2])
Stm <- cbind(base1[,1],St)
###
dbase1 <- diff(c(0,base1[,2]))
dcum <- cbind(base1[,1],dbase1)
maxtime <- tail(base1[,1],1)

if (is.null(fz)) fz <- function(x) x

 if (var.z>0) {
	 z <- rgamma(n,1/var.z)*var.z 
	 fzz <- fz(z)
	 mza <- mean(fzz)
	 if (n<10000) {
	    zl <- rgamma(100000,1/var.z)*var.z 
	    fzl <- fz(z)
	    mza <- mean(fzl)
	 } 
	 fzz <- fzz/mza
 }  else fzz <- z <- rep(1,n)

if (var.z==0) model <- "frailty"
if (model[1]=="twostage") type <- 2 else type <- 1

 ## survival censoring given X, Z, either twostage or frailty-model 
 dd <- .Call("_mets_simSurvZ",as.matrix(rbind(c(0,1),Stm)),rd,z,var.z,type)
 dd <- data.frame(time=dd[,1],status=(dd[,1]<maxtime))
 if (!is.null(cens)) cens <- rexp(n)/(rc*cens)
 dd$status <- ifelse(dd$time<cens,dd$status,0)
 dd$time <- pmin(dd$time,cens)

 ## draw recurrent process given X,Z with rate:
 ##  1/S(t|X,Z) exp(X^t beta_1) d \Lambda_1(t)
 dcum <- cbind(base1[,1],dbase1)
 ll <- .Call("_mets_simGL",as.matrix(rbind(0,dcum)),c(1,St),r1,rd,z,fzz,dd$time,type,var.z,100)
 colnames(ll) <- c("id","start","stop","death")
 ll <- data.frame(ll)
 ll$status <- 1; 
 ll$death <- dd$status[ll$id+1]
 ids <- countID(ll)
 ll <- cbind(ll,ids[,c(2,4,5)]); 
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

