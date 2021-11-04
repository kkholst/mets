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
##' and cause (code), censoring-codes (cens.code) and death-codes (death.code) indentifies these. See example. 
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
##' ll <- recreg(Event(start,stop,statusG)~x+cluster(id),data=rr,cause=1)
##' summary(ll)
##' 
##' ## censoring stratified after quartiles of x
##' lls <- recreg(Event(start, stop, statusG) ~ x+cluster(id), data=rr,cause=1, cens.model=~strata(gx))
##' summary(lls)
##' 
##' @aliases EventCens  strataAugment
##' @export
recreg <- function(formula,data=data,cause=1,death.code=c(2),cens.code=0,cens.model=~1,weights=NULL,offset=NULL,Gc=NULL,...)
{# {{{
    cl <- match.call()# {{{
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster","offset","strataAugment")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (class(Y)=="EventCens") stop("Change to Event call, see example, EventCens disabled")
    if (class(Y)!="Event") stop("Expected a 'Event'-object")
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

   if (!is.null(stratapos <- attributes(Terms)$specials$strataAugment)) {
    ts <- survival::untangle.specials(Terms, "strataAugment")
    Terms  <- Terms[-ts$terms]
    strataAugment <- as.numeric(m[[ts$vars]])-1
    strataAugment.name <- ts$vars
  }  else { strataAugment <- NULL; strataAugment.name <- NULL}

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
		      cens.model=cens.model, cause=cause, strata.name=strata.name, strataA=strataAugment,
		      death.code=death.code,cens.code=cens.code,Gc=Gc,...),
             list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster,n=nrow(X),nevent=sum(status==cause))
             )

    class(res) <- c("phreg","recreg")
    return(res)
}# }}}

recreg01 <- function(data,X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,strataA=NULL,
          strata.name=NULL,beta,stderr=1,method="NR",no.opt=FALSE, propodds=NULL,profile=0,
          case.weights=NULL,cause=1,death.code=2,cens.code=0,Gc=NULL,cens.model=~+1,augmentation=0,cox.prep=FALSE,...) {# {{{ setting up weights, strata, beta and so forth before the action starts# {{{ p <- ncol(X)
    p <- ncol(X)
    if (missing(beta)) beta <- rep(0,p)
    if (p==0) X <- cbind(rep(0,length(exit)))

    cause.jumps <- which(status==cause)
    max.jump <- max(exit[cause.jumps])
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
        if (class(cens.model)[1]=="formula") {
            formC <- update.formula(cens.model,Surv(entry__,exit__,statusC)~ . +cluster(id))
            cens.model <- phreg(formC,data)
        }
        if (cens.model$p>0) kmt <- FALSE
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
    }# }}}


    Zcall <- cbind(status,cens.strata,Stime,cens,strata,strataA) ## to keep track of status and Censoring strata
    ## setting up all jumps of type "cause", need S0, S1, S2 at jumps of "cause"
    stat1 <- 1*(status==cause)
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
            GtsAl<- Gts <- apply(Gts,2,function(x) exp(cumsum(log(1-x))))
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
        entryo <- exit[other]
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
        jumpsC <- which((xx2$Z[,4]==cens.code) & xx2$sign==1)
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
    } else {varmc <- var1 <- 0; MGc <- iH <- UUiid <- Uiid <- NULL}
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
	cause=cause,
	strataA=strataA,nstrataA=nstrataA,augment=augment.new,MGAc=MGAc
	)

if (cox.prep) out <- c(out,list(cox.prep=xx2))

    return(out)
}# }}}

####' @export
###EventCens <- function(time,time2=TRUE,cause=NULL,cens=NULL,cens.code=0,...) {# {{{
###    out <- cbind(time,time2,cause,cens)
###    colnames(out) <- c("entry","exit","cause","cens")
###    class(out) <- "EventCens"
###    attr(out,"cens.code") <- cens.code
###    return(out)
###}# }}}

##' @export
strataAugment <- survival:::strata

simRecurrentCox <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,X=NULL,r1=NULL,r2=NULL,rd=NULL,rc=NULL,...)
{# {{{
  if (is.null(r1)) r1 <- rep(1,n)
  if (is.null(r2)) r2 <- rep(1,n)
  if (is.null(rd)) rd <- rep(1,n)
  if (is.null(rc)) rc <- rep(1,n)

 ################################################################
 ### approximate hazards to make marginals fit (approximately)
 ################################################################

 ## addapt to make recurrent mean on cox form with this baseline
 base1 <- cumhaz
 if (is.null(death.cumhaz)) stop("Modification for death, otherwise use simRecurrentII\n")
 if (is.null(X)) stop("X must be given to link with simulated data\n"); 

 St <- exp(-Cpred(rbind(c(0,0),death.cumhaz),base1[,1])[,2])
 rds <- unique(rd)
 dtt <- diff(c(0,base1[,1]))
 dbase1 <- diff(c(0,base1[,2]))
 data <- c()
 XX <- c()

 nks <- 0
 for (rdss in rds) {
    lam1ms <- (dbase1)/St^rdss
    where <- which(rd==rdss)
    nk <- length(where)
    cumhaz1 <- cbind(base1[,1],cumsum(lam1ms))

    datss <- simRecurrentII(nk,cumhaz1,cumhaz2,death.cumhaz=t(t(death.cumhaz)*c(1,rdss)),
		   r1=r1[where],r2=NULL,rd=NULL,rc=rc[where],...)
    Xs <- X[where,,drop=FALSE][datss$id,,drop=FALSE]
    XX <- rbind(XX,Xs)
    datss$id <- datss$id+nks
    nks <- nks+nk
    data <- rbind(data,as.matrix(datss))
 }

 rownames(XX) <- NULL
 rownames(data) <- NULL

 return(list(data=data,X=XX))
}# }}}

simMarginalMeanCox <- function(n,cens=3/5000,k1=0.1,k2=0,bin=1,Lam1=NULL,Lam2=NULL,LamD=NULL,beta1=rep(0,2),betad=rep(0,2),betac=rep(0,2),...)
{# {{{
###

### to avoid R-check error
 revnr <- death <- status <- NULL

 if (bin==1) X <- matrix(rbinom(n*2,1,0.5),n,2) else  X <- matrix(rnorm(n*2),n,2)
 colnames(X) <- paste("X",1:2,sep="")
 r1 <- exp( X %*% beta1)
 rd <- exp( X %*% betad)
 rc <- exp( X %*% betac)

 if (is.null(Lam2)) Lam2 <- Lam1; 

 rr <- simRecurrentCox(n,t(t(Lam1)*c(1,k1)),cumhaz2=t(t(Lam1)*c(1,k2)),
		       death.cumhaz=LamD,X=X,cens=cens,r1=r1,rd=rd,rc=rc,...)
 rr <- as.data.frame(cbind(rr$data,rr$X))
 dsort(rr) <- ~id+start
 nid <- max(rr$id)
 rr$revnr <- revcumsumstrata(rep(1,nrow(rr)),rr$id-1,nid)
 rr$statusG <- rr$status 
 rr <- dtransform(rr,statusG=3,death==1)

 if (bin==0) dcut(rr,breaks=4) <- X1g~X1 else rr$X1g <- rr$X1
 if (bin==0) dcut(rr,breaks=4) <- X2g~X2 else rr$X2g <- rr$X2
 return(rr)
}# }}}


