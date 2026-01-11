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
##' @param no.codes certain event codes to be ignored when finding competing causes, can be used with administrative censoring.
##' @param death.code can also specify death.code (in addition to cause) to overrule default which takes all remaining codes (minus cause,cens.code,no.codes)
##' @param ... Additional arguments to recreg 
##' @author Thomas Scheike
##' @examples
##' ## data with no ties
##' library(mets)
##' data(bmt,package="mets")
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
##' ## approximate standard errors 
##' por <-mets:::predict.phreg(or,nd)
##' plot(por,se=1)
##'
##' ## Fine-Gray model
##' fg=cifregFG(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
##' summary(fg)
##' ##fg=recreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,death.code=2)
##' ##summary(fg)
##' plot(fg)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pfg <- predict(fg,nd,se=1)
##' plot(pfg,se=1)
##' 
##' ## bt <- iidBaseline(fg,time=30)
##' ## bt <- IIDrecreg(fg$cox.prep,fg,time=30)
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
##' Biid <- iidBaseline(fg,time=20)
##' pfg1 <- FGprediid(Biid,nd)
##' pfg1
##' @aliases vecAllStrata diffstrata FGprediid indexstratarightR gofFG cifregFG
##' @export
cifreg  <- function(formula,data,propodds=1,cause=1,cens.code=0,no.codes=NULL,death.code=NULL,...)
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
    ## default death.code are all other codes 
    if (is.null(death.code)) death.code <- all.codes[-mcodes]

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
iidBaseline.cifreg <- function(object,time=NULL,...)
{# {{{
  out <- iidBaseline.recreg(object,time=time,...)
  return(out)
} # }}}

iidBaseline.recreg <- function(object,time=NULL,ft=NULL,fixbeta=NULL,beta.iid=object$iid,tminus=FALSE,...)
{ ## {{{
	if (is.null(object$cox.prep)) stop("must call cifreg/recreg with cox.prep=TRUE\n")
	out <- IIDrecreg(object$cox.prep,object,time=time,fixbeta=fixbeta,beta.iid=beta.iid,adm.cens=object$adm.cens,tminus=tminus,...)
	out$design <- object$design
	return(out)
}  ## }}}

##' @export
IC.cifreg <- function(x, time=NULL, ...) {# {{{
  if (!is.null(time)) {
    iid <- iidBaseline.recreg(x, time=time,...)
    res <- iid$base.iid
    attr(res, "strata.level") <- iid$strata.level
    attr(res, "coef") <- iid$cumhaz.time
    attr(res, "time") <- time
    return(res * NROW(res))
  }
  res <- with(x, iid * NROW(iid))
  return(res)
}
# }}}

##' @export
estimate.cifreg <- function(x, ..., time = NULL, baseline.args = list()) {
  if (NCOL(model.matrix(x))==0L & is.null(time)) stop("Non-parametric model; need `time` argument")
  ic <- do.call(IC, c(list(x, time = time), baseline.args))
  cc <- attr(ic, "coef")
  if (is.null(cc)) cc <- coef(x)
  lab <- names(cc)
  if (!is.null(time)) {
    lab <- x$strata.level
    if (is.null(lab)) lab <- "cif"
  }
  b <- lava::estimate(
      coef = cc, IC = ic, labels = lab,
      function(x) 1 - exp(-x)
  )
  return(lava::estimate(b, ...))
}

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
##' library(mets)
##' set.seed(100)
##' rho1 <- 0.2; rho2 <- 10
##' n <- 100
##' beta=c(0.0,-0.1,-0.5,0.3)
##' dats <- simul.cifs(n,rho1,rho2,beta,rc=0.2)
##' dtable(dats,~status)
##' dsort(dats) <- ~time
##' fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
##' summary(fg)
##' plot(fg);  lines(attr(dats,"Lam1"),col=2)
##'
##' fgaugS <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fg$E)
##' summary(fgaugS)
##' fgaugS2 <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fgaugS$E)
##' summary(fgaugS2)
##'
##' @aliases strataC drop.strata
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

  call.id <- id;
  conid <- construct_id(id,nrow(X),as.data=TRUE)
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid

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


