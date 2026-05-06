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


##' @export
IC.cifreg <- function(x, time=NULL,all=FALSE,...) {# {{{
  if (!is.null(time)) {
    iid <- iidBaseline.recreg(x, time=time,...)
    if (all)  {
       res <- with(iid, cbind(beta.iid,base.iid) * NROW(base.iid))
       coefs <- c(iid$coef,iid$cumhaz.time) 
       colnames(res) <- names(coefs)
       attr(res, "coef") <- coefs
    } else { 
        res <- with(iid, base.iid * NROW(base.iid))
        attr(res, "coef") <- c(iid$cumhaz.time)
    }

    attr(res, "strata.level") <- iid$strata.level
    attr(res, "time") <- time
    return(res)
  }
  res <- with(x, iid * NROW(iid))
  return(res)
}
# }}}

##' @export
estimate.cifreg <- function(x, ..., time = NULL, all=FALSE, baseline.args = list()) {
  if (NCOL(model.matrix(x))==0L & is.null(time)) stop("Non-parametric model; need `time` argument")
  ic <- do.call(IC, c(list(x, time = time,all=all), baseline.args))
  cc <- attr(ic, "coef")
  if (is.null(cc)) cc <- coef(x)
###  lab <- names(cc)
  if (!is.null(time)) {
    lab <- x$strata.level
    if (is.null(lab)) lab <- "cif"
  }
  b <- lava::estimate( coef = cc, IC = ic)
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
predict.cifreg <- function(object,newdata=NULL,se=FALSE,times=NULL,np=50,...) { ## {{{
	call.times <- times
if (!is.null(object$propodds) & se) {se <- FALSE; warning("standard errors not computed for logit link\n"); }
if (!se) out <- predict.phreg(object,newdata,se=se,times=times,...)
else {
  out <- predictrecreg(object,newdata,times=times,np=np,...)
}
class(out) <- c("predictcifreg",class(object)[1])
out$call.times <- call.times
return(out)
} ## }}}

##' @export
summary.predictcifreg <- function(object,times=NULL,type=c("cif","cumhaz","surv")[1],...) {# {{{
ret <- summary.predictrecreg(object,type=type[1],times=times,...)
return(ret)
}# }}}

##' @export
print.predictcifreg <- function(x,...) {# {{{
ret <- summary.predictrecreg(x,...)
cat("Predictions displayed, for rows:\n")
print(ret$rows)
if (!is.null(call.times))  {
  cat("t- Predictions based on predict object, for times:\n")
  print(out$times)
}
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

#####' @export
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
  if (!inherits(iidBase,"iidBaseline")) stop("must be a iidBaseline object with iid decomposition \n"); 
	if (!is.null(newdata)) xx <- update_design(iidBase$design,data = newdata,response=FALSE) else xx <- iidBase$design
	X <- xx$x
	des <- list()
	des$X <- X
	if (!is.null(xx$strata)) strataNew <- as.numeric(xx$strata)-1 else strataNew <- rep(0,nrow(X))
###	des$strata <- strataNew

###  des <- readPhreg(iidBase,newdata)
  strata <- strataNew


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

###strataC <- survival::strata

