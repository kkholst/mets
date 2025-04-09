##' CIF regression
##'
##' CIF logistic-link for propodds=1 default and CIF Fine-Gray (cloglog) regression for propodds=NULL. The
##' FG model can also be called using the cifregFG function that has propodds=NULL.
##'
##' For FG model:
##' \deqn{
##' \int (X - E ) Y_1(t) w(t) dM_1
##' }
##' is computed and summed over clusters and returned multiplied with inverse
##' of second derivative as iid.naive. Here \deqn{w(t) = G(t) (I(T_i \wedge t < C_i)/G_c(T_i \wedge t))} and
##' \deqn{E(t) = S_1(t)/S_0(t)} and \deqn{S_j(t) = \sum X_i^j Y_{i1}(t) w_i(t) \exp(X_i^T \beta)}.
##'
##' The iid decomposition of the beta's, however, also have a censoring term that is also
##' is computed and added (still scaled with inverse second derivative)
##' \deqn{
##' \int (X - E ) Y_1(t) w(t) dM_1 + \int q(s)/p(s) dM_c
##' }
##' and returned as the iid 
##'
##' For logistic link standard errors are slightly to small since uncertainty from recursive baseline is not considered, so for smaller
##' data-sets it is recommended to use the prop.odds.subdist of timereg that is also more efficient due to use of different weights for
##' the estimating equations. Alternatively, one can also bootstrap the standard errors.
##'
##' @param formula formula with 'Event' outcome
##' @param data data frame
##' @param propodds to fit logit link model, and propodds=NULL to fit Fine-Gray model
##' @param cause of interest
##' @param cens.code code of censoring
##' @param no.codes certain event codes to be ignored when finding competing causes
##' @param ... Additional arguments to recreg 
##' @author Thomas Scheike
##' @examples
##' ## data with no ties
##' library(mets)
##' data(bmt,package="timereg")
##' bmt$time <- bmt$time+runif(nrow(bmt))*0.01
##' bmt$id <- 1:nrow(bmt)
##'
##' ## logistic link  OR interpretation
##' or=cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
##' summary(or)
##' par(mfrow=c(1,2))
##' plot(or)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' por <- predict(or,nd)
##' plot(por)
##'
##' ## Fine-Gray model
##' fg=cifregFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
##' summary(fg)
##' plot(fg)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pfg <- predict(fg,nd,se=1)
##' plot(pfg,se=1)
##' 
##' ## not run to avoid timing issues
##' ## gofFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
##' 
##' sfg <- cifregFG(Event(time,cause)~strata(tcell)+platelet+age,data=bmt,cause=1)
##' summary(sfg)
##' plot(sfg)
##'
##' ### predictions with CI based on iid decomposition of baseline and beta
##' ### these are used in the predict function above
##' fg <- cifregFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
##' Biid <- IIDbaseline.cifreg(fg,time=20)
##' pfg1 <- FGprediid(Biid,nd)
##' pfg1
##' @aliases vecAllStrata diffstrata IIDbaseline.cifreg FGprediid indexstratarightR gofFG cifregO cifregFG IIDbaseline.cifregO
##' @export
cifreg  <- function(formula,data,propodds=1,cause=1,cens.code=0,no.codes=NULL,...)
{# {{{
    cl <- match.call()
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster","offset","strataAugment")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Event")) stop("Expected a 'Event'-object, with codes for terminal events (death.code if any), censoring (cens.code), and event of interest (cause)")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- rep(0,nrow(Y))
        status <- Y[,2]
    } else {
        entry <- Y[,1]
        exit <- Y[,2]
        status <- Y[,3]
    }
    ## default version of codes, otherwise call recregN directly
    all.codes <-  unique(status)
    codes <- c(cause,cens.code) 
    if (!is.null(no.codes)) codes <- c(codes,no.codes) 
    mcodes <- match(codes,all.codes,nomatch=0)
    death.code <- all.codes[-mcodes]

    cif <- recreg(formula,data,propodds=propodds,cause=cause,cens.code=cens.code,death.code=death.code,...)
    class(cif) <- c("cifreg","phreg")
    cif$call <- cl
return(cif)
} # }}}

##' @export
cifregFG  <- function(formula,data,propodds=NULL,...)
{# {{{
cif <- cifreg(formula,data,propodds=propodds,...)
return(cif)
} # }}}

##' @export
IIDbaseline.cifreg <- function(x,...)
{# {{{
   out <- IIDbaseline.recreg(x,...)
   return(out)
} # }}}

##' @export
cifregO <- function(formula,data=data,cause=1,cens.code=0,cens.model=~1,weights=NULL,offset=NULL,Gc=NULL,propodds=1,...)
{ # {{{
    cl <- match.call() # {{{
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster", "offset")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Event")) stop("Expected a 'Event'-object")
    if (ncol(Y) == 2) {
        exit <- Y[, 1]
        entry <- NULL ## rep(0,nrow(Y))
        status <- Y[, 2]
    } else {
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
    } else {
        pos.cluster <- NULL
    }
    if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
        ts <- survival::untangle.specials(Terms, "strata")
        pos.strata <- ts$terms
        Terms <- Terms[-ts$terms]
        strata <- m[[ts$vars]]
        strata.name <- ts$vars
    } else {
        strata.name <- NULL
        pos.strata <- NULL
    }
    if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
        ts <- survival::untangle.specials(Terms, "offset")
        Terms <- Terms[-ts$terms]
        offset <- m[[ts$vars]]
    }
    X <- model.matrix(Terms, m)
    if (!is.null(intpos <- attributes(Terms)$intercept)) {
        X <- X[, -intpos, drop = FALSE]
    }
    if (ncol(X) == 0) X <- matrix(nrow = 0, ncol = 0)

    ## }}}

    res <- c(
        cifreg01(data, X, exit, status, id, strata, offset, weights, strata.name,
            cens.model = cens.model,
            cause = cause, cens.code = cens.code, Gc = Gc, propodds = propodds, ...
        ),
        list(
            call = cl, model.frame = m, formula = formula, strata.pos = pos.strata,
            cluster.pos = pos.cluster, n = nrow(X), nevent = sum(status == cause)
        )
    )

    class(res) <- c("cifreg", "phreg")
    return(res)
} # }}}

cifreg01 <- function(data,X,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
              strata.name=NULL,beta,stderr=TRUE,method="NR",no.opt=FALSE,propodds=1,profile=0,
              case.weights=NULL,cause=1,cens.code=0,Gc=NULL,cens.model=~+1,augmentation=0,cox.prep=FALSE,
	      adm.cens.code=NULL,adm.cens.time=NULL,...) {# {{{
    ##  setting up weights, strata, beta and so forth before the action starts# {{{
    p <- ncol(X)
    if (missing(beta)) beta <- rep(0,p)
    if (p==0) X <- cbind(rep(0,length(exit)))

    cause.jumps <- which(status %in% cause)
    max.jump <- max(exit[cause.jumps])
    other <- which((!(status %in% c(cens.code,cause,adm.cens.code)) ) & (exit< max.jump))

    entry <- NULL
    n <- length(exit)
    trunc <- (!is.null(entry))
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

    if (!trunc) entry <- rep(0,length(exit))
    if (is.null(offset)) offset <- rep(0,length(exit))
    if (is.null(weights)) weights <- rep(1,length(exit))
    if (is.null(case.weights)) case.weights <- rep(1,length(exit))
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
    whereC <- which(status %in% cens.code)
    time <- exit
    statusC <- (status %in% cens.code)
    data$id <- id
    data$exit <- exit
    data$statusC <- statusC
    cens.strata <- cens.nstrata <- NULL

    if (length(whereC)>0) {# {{{
    if (is.null(Gc)) {
        kmt <- TRUE
        if (inherits(cens.model,"formula")) {
            formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id))
            cens.model <- phreg(formC,data)
        }
        if (cens.model$p>0) kmt <- FALSE
        Pcens.model <- suppressWarnings(predict(cens.model,data,times=exit,tminus=TRUE,individual.time=TRUE,se=FALSE,km=kmt))
        Stime <- Pcens.model$surv <- c(Pcens.model$surv)
        ## strata from original data
        nCstrata <- cens.model$nstrata
        cens.strata <- cens.model$strata[order(cens.model$ord)]
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

    trunc <- FALSE
    Zcall <- cbind(status,cens.strata,Stime,adm.cens.time) ## to keep track of status and Censoring strata
    ## setting up all jumps of type "cause", need S0, S1, S2 at jumps of "cause"
    stat1 <- 1*(status %in% cause)
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
            if (is.null(dim(Gts))) Gts <- matrix(Gts,1,1)
            ### back to km
            Gts <- suppressWarnings(apply(Gts,2,function(x) exp(cumsum(log(1-x)))))
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
        entryo <- exit[other]
        ido <- id[other]
        stratao <- strata[other]
        ###
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
            DUt <-  .Call("XXMatFULL",DUt,p,PACKAGE="mets")$XXf
	    if (ncol(DUt)>0) DUt <- DUt+DUadj
            if (profile==1) {
		Ut <- Ut+c(ploglik)*Dwbeta
		## not implemented
		DUt <- DUt
            }
            ploglik <- pow*ploglik
        }

        U  <- apply(Ut,2,sum)
        DUt <- caseweightsJ*weightsJ*DUt
	if (ncol(DUt)!=p*p) {
        DU <- matrix(0,p,p);
        DU[lower.tri(DU,diag=TRUE)] <- -apply(DUt,2,sum)
        DU<- DU+t(DU)
        diag(DU) <- diag(DU)/2
	} else  DU <- -matrix(apply(DUt,2,sum),p,p)
        ploglik <- sum(ploglik)
        U <- U+augmentation

        out <- list(ploglik=ploglik,gradient=U,hessian=-DU,cox.prep=xx2,
                    hessiantime=DUt,weightsJ=weightsJ,
		    caseweightsJ=caseweightsJ,
                    jumptimes=jumptimes,
		    strata.jumps=strataJ,
		    strata=strataJ,
		    nstrata=nstrata,
                    time=jumptimes,
		    S0=S0/(caseweightsJ*weightsJ),
		    S2S0=S2S0,E=E,U=Ut,X=Xj,Gjumps=Gjumps)


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
            if (!stderr) return(cc)
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

    if (is.null(adm.cens.time)) {
    if (length(other)>=1) { ## martingale part for type-2 after T
    ### xx2 data all data
        otherxx2 <- which(!(xx2$Z[,1] %in% c(cause,cens.code,adm.cens.code)))
        rr0 <- xx2$sign
        jumpsC <-  which(xx2$Z[,1] %in% cens.code)
        strataCxx2 <- xx2$Z[,2]
        S0iC2  <-  S0iC <- rep(0,length(xx2$status))
        S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
        S0iC[jumpsC] <- 1/S0rrr[jumpsC]
        S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
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
    }

    if (!is.null(adm.cens.time)) {
    if (length(other)>=1) { ## martingale part for type-2 after T
        otherxx2 <- which(!(xx2$Z[,1] %in% c(cause,cens.code,adm.cens.code)))
        rr0 <- xx2$sign
        jumpsC <-  which(xx2$Z[,1] %in% cens.code)
        strataCxx2 <- xx2$Z[,2]
        S0iC2  <-  S0iC <- rep(0,length(xx2$status))
        S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
        S0iC[jumpsC] <- 1/S0rrr[jumpsC]
        S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
        Gcxx2 <- exp(cumsumstrata(log(1-S0iC),strataCxx2,nCstrata))
        Gstart <- rep(1,nCstrata)

        dstrata <- mystrata(data.frame(strataCxx2,xx2$strata))
        ndstrata <- attr(dstrata,"nlevel")
        lastt <- tailstrata(dstrata-1,ndstrata)
	if (!is.null(adm.cens.time)) {
            act <- xx2$Z[otherxx2,4] 
	    dact <- dstrata[otherxx2]
            wherea <- indexstratarightR(xx2$time,dstrata-1,act,dact-1,ndstrata,type="left")
	} else wherea <- lastt[dstrata][otherxx2]
        ll <-  cumsum2strata(Gcxx2,S0i,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)
        Htsj <- ll$res[wherea]-ll$res[otherxx2]
        fff <- function(x) {
            cx  <- cumsum2strata(Gcxx2,x*S0i,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)
            cx <- cx$res[wherea]-cx$res[otherxx2]
            return(cx)
        }
        EHtsj  <- apply(E,2,fff)
        rrx2 <- rr[otherxx2]*xx2$weights[otherxx2]/xx2$Z[otherxx2,3]
        MGt2  <- -(Z[otherxx2,,drop=FALSE]*Htsj-EHtsj)*rrx2
        UU2 <- apply(MGt2,2,sumstrata,xx2$id[otherxx2],mid+1)
        UU  <-  UU+UU2
    }
    }
    ## }}}


    if ((length(other)>=1) & (length(whereC)>0) & (is.null(adm.cens.time))) { ### Censoring adjustment for jumps of other type but only for KM-case {{{

        EHtsja <- Xos <- matrix(0,nrow(Z),ncol(Z));
        Xos[otherxx2,] <- Z[otherxx2,]*rrx2
        rrx <- rep(0,nrow(Z))
        rrx[otherxx2] <- rrx2
        rrsx <- cumsumstrata(rrx,strataCxx2,nCstrata)
        Xos <- apply(Xos,2,cumsumstrata,strataCxx2,nCstrata)
	Htsja <- rep(0,nrow(Z))
###	Htsja[otherxx2] <- Htsj
###	EHtsja[otherxx2,] <- EHtsj
        q2 <- (Xos*c(Htsj)-EHtsj*c(rrsx))

        sss <- headstrata(dstrata-1,ndstrata)
        fff <- function(x) {
            gtstart <- x[sss]
            cx  <- cumsum2strata(x,S0iC2,dstrata-1,ndstrata,strataCxx2,nCstrata,gtstart)$res
            return(cx)
        }
        EdLam0q2 <- apply(q2,2,fff)

        ### Martingale  as a function of time and for all subjects to handle strata
        MGc <- q2*S0iC-EdLam0q2
        MGc <- apply(MGc,2,sumstrata,xx2$id,mid+1)
        ## }}}
    } else MGc <- 0

    iH <- - tryCatch(solve(opt$hessian),error=function(e) matrix(0,nrow(opt$hessian),ncol(opt$hessian)) )
    Uiid <-  (UU+MGc) %*% iH
    UUiid <- UU %*% iH
    var1 <-  crossprod(UUiid)
    varmc <-  crossprod(Uiid)
    ### end if (p>0)
    } else {varmc <- var1 <- 0; MGc <- iH <- UUiid <- Uiid <- NULL}
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


  if (!is.null(Uiid)) {
     colnames(Uiid) <- names(beta.s)
     if (nrow(Uiid) == nrow(data)) rownames(Uiid) <- rownames(data)
  }

    out <- list(coef=beta.s,var=varmc,se.coef=diag(varmc)^.5,iid.naive=UUiid,
                iid=Uiid,ncluster=nid,
                ihessian=iH,hessian=opt$hessian,
		hessianttime=opt$hessiantime,
		strata.jumps=opt$strata.jumps,
		var1=var1,se1.coef=diag(var1)^.5,
                ploglik=opt$ploglik,gradient=opt$gradient,
                cumhaz=cumhaz, se.cumhaz=se.cumhaz,MGciid=MGc,
		strata.call=strata.call,
                strata=xx2$strata, nstrata=nstrata,strata.name=strata.name,strata.level=strata.level,
		propodds=propodds,
                S0=opt$S0,E=opt$E,S2S0=opt$S2S0,time=opt$time,Ut=opt$U,
                jumps=jumps,exit=exit,p=p,
		id=id.orig,call.id=call.id,
		no.opt=no.opt,##n=nrow(X),nevent=length(jumps),
                Pcens.model=Pcens.model,Gjumps=Gjumps,cens.code=cens.code,cause=cause
                )

    if (cox.prep) out <- c(out,list(cox.prep=xx2))

    return(out)
}# }}}

##' @export
IC.cifreg <- function(x, ...) {# {{{
  res <- with(x, iid * NROW(iid))
  return(res)
}
# }}}

##' @export
plot.cifreg <- function(x,se=FALSE,ylab=NULL,...) { ## {{{
   if (inherits(x,"cifreg") & is.null(ylab)) ylab <- "Cumulative Hazard"
if (!se) baseplot(x,se=se,ylab=ylab,...)
else {
    warning("Standard errors approximative (but too small), use predict and type='cumhaz' \n")
    baseplot(x,se=se,ylab=ylab,...)
}
} ## }}}

##' @export
predict.cifreg <- function(object,newdata,se=FALSE,times=NULL,np=50,...) { ## {{{
if (!is.null(object$propodds) & se) {se <- FALSE; warning("standard errors not computed for logit link\n"); }
if (!se) out <- predict.phreg(object,newdata,se=se,times=times,...)
else {
  out <- predictrecreg(object,newdata,times=times,np=np,...)
}
class(out) <- c("predictcifreg",class(object)[1])
return(out)
} ## }}}

##' @export
summary.predictcifreg <- function(object,times=NULL,type=c("cif","cumhaz","surv")[1],...) {# {{{
ret <- summary.predictrecreg(object,type=type[1],times=times,...)
return(ret)
}# }}}

##' @export
plot.predictcifreg <- function(x,se=FALSE,ylab=NULL,type="cif",...)
{ ## {{{
  if (inherits(x,"predictcifreg") & is.null(ylab)) ylab <- "Probability"
  plotpredictphreg(x,se=se,ylab=ylab,type=type[1],...)
} ## }}}

##' @export
gofFG <- function(formula,data,cause=1,cens.code=0,cens.model=NULL,...)
{# {{{
fgform <- update(formula, paste("Surv(fgstart, fgstop, fgstatus) ~ .+cluster(id)"))
## assumes simple Surv(time,status) is given 
vars <- all.vars(formula)
data$id <- 1:nrow(data)
formid <- update.formula(formula,paste(".~.+id"))

## process types to get type of interst and other types 
status <- data[,vars[2]]
types <- unique(status)
mm <- match(c(cens.code,cause),types)
mm <- mm[!is.na(mm)]
statusS <- status
statusS[!(status %in% c(cens.code,cause))] <- types[-mm][1]
typesF <- c(cens.code,cause,types[-mm][1])
###
###data[,vars[2]] <- factor(statusS,typesF,c("censoring","cause","ocause"))
data[,vars[2]] <- factor(statusS,typesF,typesF)

Xs <- vars[-(1:2)]
modP <- paste(Xs,collapse="+") 
formid <- as.formula(paste("Surv(",vars[1],",",vars[2],")~",modP,"+id"))
if (!is.null(cens.model)) {
Cstrata <- as.character(cens.model)
formid <- as.formula(paste("Surv(",vars[1],",",vars[2],")~",modP,"+",Cstrata[2],"+id"))
}

fgdata <- finegray(formid,data=data)
fgcph <- phreg(fgform,data=fgdata,weights=fgdata$fgwt,...)
ggmg <- gof(fgcph)
return(ggmg)
}# }}}

##' @export
indexstratarightR <- function(timeo,stratao,jump,js,nstrata,type="right")# {{{
{
###    if (any(stratao < 0 | stratao >= nstrata))
###        stop("time-strata strata not between 0-(nstrata-1)\n")
###    if (any(js < 0 | js >= nstrata))
###        stop("jump-strata strata not between 0-(nstrata-1)\n")
    mm <- cbind(timeo, stratao, 1:length(timeo), 0)
    mm <- rbind(mm, cbind(jump, js, 1:length(jump), 1))
    ord <- order(mm[, 1], mm[, 4])
    mm <- mm[ord, ]
    if (type == "right") right <- 1 else right <- 0
    ires <- .Call("indexstrataR", mm[, 2], mm[, 3], mm[, 4], nstrata,right)
    res <- ires$res
    or2 <- order(res[,2])
    res <- res[or2,1]
    reso <- which(res==0) 
    if (length(reso)>0) {
       jso <- js[reso]
       for (s in 1:nstrata) if (length(jso)>0) res[reso[jso==s-1]] <- ires$maxmin[s] 
    }
    return(res)
} ## }}}

##' @export IIDbaseline.cifreg 
IIDbaseline.cifregO <- function(x,time=NULL,fixbeta=NULL,...)
{# {{{
###  sum_i int_0^t 1/S_0(s) dM_{ki}(s) - P(t) \beta_k
###  with possible strata and cluster "k", and i in clusters 
  if (length(class(x))!=2) stop("Must be cifreg object\n"); 
  if (!inherits(x,c("cifreg","recreg"))) stop("Must be cifreg/recreg object\n"); 
  if (is.null(time)) stop("Must give time for iid of baseline")

  if (!is.null(x$propodds))  stop("Only for Fine-Gray-model") 
  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if ((x$no.opt) | is.null(x$coef)) fixbeta<- 1 else fixbeta <- 0

  if (is.null(x$cox.prep)) stop("call cifreg with cox.prep=TRUE\n"); 

###  read relevant data from cox.prep argument of cifreg
# {{{
  xx2 <- x$cox.prep
  btimexx <- c(1*(xx2$time < time))

  status <- xx2$Z[,1]
  cause.jumps <- xx2$jumps+1 
  exit <- xx2$time
  max.jump <- max(exit[cause.jumps])+1
  if (inherits(x,c("cifreg"))) 
  other <- which((!(status %in% c(x$cens.code,x$cause)) ) )
  else 
  other <- which((status %in% x$death.code) & (xx2$sign==1) )
  whereC <- which( (status %in% x$cens.code) & xx2$sign==1)

  jumps <- xx2$jumps+1
  jumptimes <- xx2$time[jumps]
  strata1jumptimes <- xx2$strata[jumps]
  Xj <- xx2$X[jumps,,drop=FALSE]
  p <- ncol(Xj)
  strataCxx2 <- xx2$Z[,2]
  nCstrata <- max(strataCxx2)+1
# }}}

### iid version given G_c {{{
    ##iid robust phreg
    S0i <- rep(0,length(xx2$strata))
    S0i[jumps] <- 1/x$S0
    Z <- xx2$X
    U <- E <- matrix(0,nrow(Z),p)
    E[jumps,] <- x$E
    U[jumps,] <- x$U
    cumhazA <- cumsumstratasum(S0i,xx2$strata,xx2$nstrata,type="all")
    cumhaz <- c(cumhazA$sum)

    if (fixbeta==0) {
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx2$strata,xx2$nstrata)
	  Ht <- apply(E*S0i*btimexx,2,cumsumstrata,xx2$strata,xx2$nstrata)
    } else Ht <- NULL
  if (fixbeta==0) rr <- c(xx2$sign*exp(Z %*% coef(x) + xx2$offset)) else rr <- c(xx2$sign*exp(xx2$offset))

### Martingale  as a function of time and for all subjects to handle strata
   mid <- max(xx2$id)
  if (fixbeta==0)  {
    MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx2$weights)
    UU <- apply(MGt,2,sumstrata,xx2$id,mid+1)
  } else  { MGt <- 0 ; UU <- 0}

    MGAiid <- NULL
    S0i2 <- rep(0,length(xx2$strata))
    ww <- xx2$caseweights*xx2$weights
    S0i2[jumps] <- 1/(x$S0^2*ww[jumps])
    MGAiid <- matrix(0,length(S0i2),1)
    MGAiid2 <- matrix(0,length(S0i2),1)
    cumhazAA <- cumsumstrata(S0i2*btimexx,xx2$strata,xx2$nstrata)
    MGAiid <- S0i*btimexx-cumhazAA*rr*c(xx2$weights)

  if (length(other)>=1) { ## martingale part for type-2 after T
    ### xx2 data all data
        otherxx2 <- other  
        rr0 <- xx2$sign
        jumpsC <- whereC 
        S0iC2  <-  S0iC <- rep(0,length(xx2$status))
        S0rrr <- revcumsumstrata(rr0,strataCxx2,nCstrata)
        S0iC[jumpsC] <- 1/S0rrr[jumpsC]
        S0iC2[jumpsC] <- 1/S0rrr[jumpsC]^2
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

        ll <-  cumsum2strata(Gcxx2,S0i2*btimexx,strataCxx2,nCstrata,xx2$strata,xx2$nstrata,Gstart)
        HBtsj <- ll$res[lastt][dstrata]-ll$res
        MGAiid2[otherxx2,] <- -HBtsj[otherxx2,,drop=FALSE]*rrx2
	MGAiid <- MGAiid+MGAiid2
    }
    ## }}}

 
    if ((length(other)>=1) & (length(whereC)>0)) { ### Censoring adjustment for jumps of other type but only for KM-case {{{

        Xos <- matrix(0,nrow(Z),ncol(Z));
        Xos[otherxx2,] <- Z[otherxx2,]*rrx2
        rrx <- rep(0,nrow(Z))
        rrx[otherxx2] <- rrx2
        rrsx <- cumsumstrata(rrx,strataCxx2,nCstrata)
        Xos <- apply(Xos,2,cumsumstrata,strataCxx2,nCstrata)
        q2 <- (Xos*c(Htsj)-EHtsj*c(rrsx))

	qB2 <- rrsx*c(HBtsj) 

        sss <- headstrata(dstrata-1,ndstrata)
        fff <- function(x) {
            gtstart <- x[sss]
            cx  <- cumsum2strata(x,S0iC2,dstrata-1,ndstrata,strataCxx2,nCstrata,gtstart)$res
            return(cx)
        }
        EdLam0q2 <- apply(q2,2,fff)
        EBdLam0q2 <- apply(qB2,2,fff)

        ### Martingale  as a function of time and for all subjects to handle strata
        MGc <- q2*S0iC-EdLam0q2*c(xx2$sign)
        MGc <- apply(MGc,2,sumstrata,xx2$id,mid+1)

        MGBc <- qB2*S0iC-EBdLam0q2*c(xx2$sign)
        ## }}}
    } else { MGc <- 0; MGBc <- 0}

  if (fixbeta==0) { betaiid <-  (UU+MGc) %*% x$ihessian} else betaiid <- NULL
   MGAiid <- MGAiid+MGBc

 ### \hat beta - \beta = \sum_i \beta_i  (iid) 
 ### iid after baseline:
 ### \hat A_s-A_s=\sum_{i clusters} \sum_{j \in i(j)=i, s(j)=s} \int_0^t 1/S_0 dM^s_j - P^s(t) \sum_i \beta_i
 ### = \sum_{i clusters} ( \sum_{j \in i(j)=i, s(j)=s} \int_0^t 1/S_0 dM^s_j - P^s(t) \beta_i ) 

 if (fixbeta==0) {# {{{
    Htlast <- tailstrata(xx2$strata,xx2$nstrata)
    HtS <- Ht[Htlast,,drop=FALSE]
 } ## }}}


## sum after id's within strata and order 
 MGAiids <- c()
 cumhaz.time <- c()
 sus <- sort(unique(xx2$strata))
 for (i in sus)  { 
	 wi <- which(xx2$strata==i)
         MGAiidl <- sumstrata(MGAiid[wi],xx2$id[wi],mid+1)
	 cumhazs <- x$cumhaz[x$strata[x$jumps]==i,,drop=FALSE]
cumhaz.time <- c(cumhaz.time,cpred(rbind(0,cumhazs),time)[,2])

        if (fixbeta==0) {
           UU <-  apply(HtS[i+1,]*t(betaiid),2,sum)
           MGAiidl <- MGAiidl - UU
         }
         MGAiids <- cbind(MGAiids,MGAiidl)
 }
 colnames(MGAiids) <- paste("strata",sus,sep="")
 names(cumhaz.time) <- paste("strata",sus,sep="")

 return(list(time=time,base.iid=MGAiids, nstrata=xx2$nstrata, beta.iid=betaiid,
	     strata.call=x$strata.call,id=x$id,call.id=x$call.id,
	     coef=coef(x),cumhaz=x$cumhaz,cumhaz.strata=x$strata[x$jumps],
	     cumhaz.time=cumhaz.time,strata.time=sus,
             nstrata=x$nstrata,strata.name=x$strata.name,strata.level=x$strata.level,
	     model.frame=x$model.frame,formula=x$formula))
} # }}}


##' @export
FGprediid <- function(iidBase,newdata,conf.type=c("log","cloglog","plain"),model="FG")
{# {{{
  des <- readPhreg(iidBase,newdata)
  strata <- des$strata
  if (!is.null(iidBase$beta.iid))  { 
	  fixbeta <- 0; beta.iid <- iidBase$beta.iid; X <- des$X; p <- ncol(beta.iid); 
  } else { fixbeta <- 1; beta.iid <- 0; X <- matrix(0,1,1); p <- 1 }

  sus <- sort(unique(strata))
  At <- iidBase$cumhaz.time[match(sus,iidBase$strata.time)]

   if (missing(X)) X <- matrix(0,1,p)
   if (ncol(X)!=p) stop("X and coef does not match \n"); 

   Ft <- function(coef,Xi=rep(0,length(p)),type="log") {
	   base <- coef[1]
	   p <- coef[-1]
   if (model=="FG") {
       if (type=="log")     y <- log(1-exp(-base*exp(Xi %*% p)))
       if (type=="plain")   y <- 1-exp(-base*exp(Xi %*% p))
       if (type=="cloglog") y <- log(-log(1-exp(-base*exp(Xi %*% p))))
   } else { ## Ghosh-Lin model
       if (type=="log")     y <- log(base*exp(Xi %*% p))
       if (type=="plain")   y <- base*exp(Xi %*% p)
   }
       return(y)
   }
   Ftback <- function(p,type="log")
   {
       if (type=="log")     y <- exp(p)
       if (type=="plain")   y <- p
       if (type=="cloglog") y <- exp(-exp(p)) 
       return(y)
   }

   preds <- matrix(0,length(strata),4)

   k <- 0
   for (i in sus) {
      wheres <- which(strata==i) 
      wi <- match(i,iidBase$strata.time)
      At <- iidBase$cumhaz.time[wi]
	if (fixbeta==0)  {
	   iidAB <- cbind(iidBase$base.iid[,wi],iidBase$beta.iid)
	   covv <- crossprod(iidAB)
	   coeff <- c(At,iidBase$coef)
	   for (j in wheres)  {
	      Xj <- X[j,]
	      eud <- estimate(coef=coeff,vcov=covv,f=function(p) Ft(p,Xi=Xj,type=conf.type[1]))
	      cmat <- eud$coefmat
	      cmat <- c(cmat[,-5])
	      cicmat <- Ftback(cmat[c(1,3:4)],type=conf.type[1])
	      cmat[c(1,3:4)] <- cicmat
	      preds[j,] <- c(cmat)
	   } 
	} else {
	      iidAB <- cbind(iidBase$base.iid[,wi])
	      covv <- crossprod(iidAB)
	      coeff <- c(At)
	      eud <- estimate(coef=coeff,vcov=covv,f=function(p) Ft(p,Xi=1,type=conf.type[1]))
	      cmat <- eud$coefmat
	      cmat <- c(cmat[,-5])
	      cicmat <- Ftback(cmat[c(1,3:4)],type=conf.type[1])
	      cmat[c(1,3:4)] <- cicmat
	      preds[wheres,] <- matrix(c(cmat),length(wheres),4,byrow=TRUE)
	}
}

colnames(preds) <- c("pred",paste("se",conf.type[1],sep="-"),"lower","upper")

return(preds)
}# }}}

##' @export
strataC <- survival::strata

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
##' \deqn{E(\beta^p,t)} is given. Assumes that no strata for baseline of ine-Gay model that is augmented. 
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
##' set.seed(100)
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
##' summary(fgaugS2)
##'
##' @aliases strataC  simul.cifs setup.cif drop.strata
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
        if (is.numeric(id))
            id <-  fast.approx(ids,id)-1
        else  {
            id <- as.integer(factor(id,labels=seq(nid)))-1
        }
    } else id <- as.integer(seq_along(exit))-1;

    p <- ncol(X)
    beta <- NULL
    if (is.null(beta)) beta <- rep(0,p)
    if (p==0) X <- cbind(rep(0,length(exit)))
    if (is.null(strata)) {
        strata <- rep(0,length(exit));
        nstrata <- 1;
        strata.level <- NULL; }
    else {
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

    if (is.null(strataC)) {
        strataC <- rep(0,length(exit))
        nstrataC <- 1
        strataC.level <- NULL
    } else {
        strataC.level <- levels(strataC)
        ustrataC <- sort(unique(strataC))
        nstrataC <- length(ustrataC)
        strataC.values <- ustrataC
        if (is.numeric(strataC))
            strataC <-  fast.approx(ustrataC,strataC)-1
        else  {
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
        if (is.numeric(id))
            id <-  fast.approx(ids,id)-1
        else  {
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
    jumps1 <- which(xxstatus %in% cause)
    jumpsD <- which(!(xxstatus %in% cens.code))
    rr <- c(dd$sign*exp(dd$offset))
    ## S0 after strata
    S0 = c(revcumsumstrata(rr,xxstrata,nstrata))
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
        cumhazD <- cumsumstratasum(S00i,xxstrataC,nstrataC)
        Gc      <- exp(-c(cumhazD$sum))
        Gcm      <- exp(-c(cumhazD$lagsum))
    } else { 
        logGc <- cumsumstratasum(log(1-S00i),xxstrataC,nstrataC)
	Gcm <- c(exp(logGc$lagsum))
	Gc  <- c(exp(logGc$sum))
    }

    cif1 <- cumsumstrata(Stm*S01i,xxstrata,nstrata)
    cif2 <- cumsumstrata(Stm*S02i,xxstrata,nstrata)
    ## }}}

    Et <- matrix(0,nrow(Z),ncol(Z))
    Et[jumps1,] <- E

    Lam1fg <- -log(1-cif1)
    ## to deal with cif1=1 in which case cif2=0
    Lam1fg[is.na(Lam1fg)] <- 0
    laststrata <- tailstrata(xxstrata,nstrata)
    gtstart <- Lam1fg[laststrata]
    dLam1fg <- c(diffstrata(Lam1fg,xxstrata,nstrata))

    ## Gc~strataC, dLam1fg~strata
    tailcstrata <- tailstrata(xxstrataC,nstrataC)
    Gcstart <- Gc[tailcstrata]

    dstrata <- mystrata2index(cbind(xxstrataC,xxstrata))
    ndstrata <- attr(dstrata,"nlevel")
    lastt <- tailstrata(dstrata-1,ndstrata)

    ### ## \int_t^\infty G_c^j(t) d\Lambda_1^k(t)
    G0start <- rep(1,nstrataC)
    cLam1fg  <- cumsum2strata(Gc,dLam1fg,xxstrataC,nstrataC,xxstrata,nstrata,G0start)
    RLam1fg <- cLam1fg$res[lastt][dstrata]-cLam1fg$res

    ## E(s) from FG without strata
    ## \int_0^t  G_c^j(s) E(s) d\Lambda_1^k(s)
    fff <- function(x) {
        cx  <- cumsum2strata(Gc,x*dLam1fg,xxstrataC,nstrataC,xxstrata,nstrata,G0start)
        return(cx$res[lastt][dstrata]-cx$res)
    }
    ERLam1fg0  <- apply(Et,2,fff)

    cif2GS <-  c(cif2/(Gc*St))
    cif2GS[St<0.00001] <- 0
    cif2GS[Gc<0.00001] <- 0
    gt <-  RLam1fg*cif2GS
    ERLam1fg  <- ERLam1fg0*cif2GS

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
                    jumps1=jumps1,jumps=jumps,other=other,
                    nstrata=nstrata,nstrataC=nstrataC,dstrata=dstrata,ndstrata=ndstrata,
                    cif1=cif1,cif2=cif2,St=St,Gc=Gc,strata=xxstrata,strataC=xxstrataC,time=dd$time)

    ## drop strata's from formula and run wiht augmention term
    formulans <- drop.strata(formula)

    if (nstrataC==1) cens.model <- ~+1 else cens.model <- ~strata(strataCC)
    data$strataCC <- cens.strata

    fga <- tryCatch(cifreg(formulans,data=data,cause=cause,
                  propodds=NULL,augmentation=augment$augment,cens.model=cens.model,...),error=function(x) NULL) 

    if (!is.null(fga)) {
    ## adjust SE and var based on augmentation term
    fga$var.orig <- fga$var
    fga$augment <- augment$augment
    fga$iid <- fga$iid.naive + MGiid %*% fga$ihessian
    fga$var <- crossprod(fga$iid)
    fga$se.coef <-  diag(fga$var)^.5
    fga$MGciid <- MGiid
    } else  {
    fga$augment <- augment$augment
    fga$MGciid <- MGiid
    }

    return(fga)
}# }}})

##' @export
simul.cifs <- function(n,rho1,rho2,beta,rc=0.5,depcens=0,rcZ=0.5,bin=1,type=c("cloglog","logistic"),rate=1,Z=NULL) {# {{{
    p=length(beta)/2
    tt <- seq(0,6,by=0.1)
    if (length(rate)==1) rate <- rep(rate,2)
    Lam1 <- rho1*(1-exp(-tt/rate[1]))
    Lam2 <- rho2*(1-exp(-tt/rate[2]))

    if (length(bin)==1) bin <- rep(bin,2)
    if (length(rcZ)==1) rcZ <- c(rcZ,0)

    if (is.null(Z)) 
    Z=cbind((bin[1]==1)*(2*rbinom(n,1,1/2)-1)+(bin[1]==0)*rnorm(n),(bin[2]==1)*(rbinom(n,1,1/2))+(bin[2]==0)*rnorm(n))
    colnames(Z) <- paste("Z",1:2,sep="")
    p <- ncol(Z)

    cif1 <- setup.cif(cbind(tt,Lam1),beta[1:p],Znames=colnames(Z),type=type[1])
    cif2 <- setup.cif(cbind(tt,Lam2),beta[(p+1):(2*p)],Znames=colnames(Z),type=type[1])
    data <- sim.cifsRestrict(list(cif1,cif2),n,Z=Z)

    if (depcens==0) censor=pmin(rexp(n,1)*(1/rc),6) else censor=pmin(rexp(n,1)*(1/(rc*exp(Z %*% rcZ))),6)

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
