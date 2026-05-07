##' Percentage of Years Lost Due to a Cause Regression
##'
##' Estimates the percentage of the restricted mean time lost (RMTL) that is attributable 
##' to a specific cause and models how covariates affect this percentage using IPCW regression.
##' 
##' Let the total years lost be \eqn{Y = t - \min(T, t)} and the years lost due to cause 1 be 
##' \eqn{Y_1 = I(\epsilon=1) (t - \min(T, t))}. The function models the ratio:
##' \deqn{ \text{logit}\left( \frac{E(Y_1 | X)}{E(Y | X)} \right) = X^T \beta }
##' 
##' Estimation is based on a binomial regression IPCW response estimating equation:
##' \deqn{ X \left( \Delta^{\text{ipcw}}(t) \left( Y \cdot \text{expit}(X^T \beta) - Y_1 \right) \right) = 0 }
##' where \eqn{\Delta^{\text{ipcw}}(t) = I(\min(t,T) < C) / G_c(\min(t,T))} is the IPCW adjustment.
##'
##' The function supports three types of estimators:
##' \itemize{
##'   \item \code{"I"}: Classical outcome IPCW regression (no augmentation).
##'   \item \code{"II"}: Adds a censoring augmentation term \eqn{X \int E(Z(t)| T>s)/G_c(s) d \hat M_c} 
##'     to improve efficiency (requires an initial estimate of \eqn{\beta}).
##'   \item \code{"III"}: Adds a more complex augmentation term separating the expectations of 
##'     \eqn{Y} and \eqn{Y_1} for further efficiency gains.
##' }
##'
##' The variance is based on the squared influence functions (IID). A "naive" variance 
##' (assuming known censoring) is also provided for comparison.
##'
##' @param formula Formula with an outcome (see \code{coxph}). The first covariate on the RHS 
##'   is typically the treatment or group indicator. Can include \code{cluster(id)}.
##' @param data Data frame containing the variables.
##' @param cause Numeric code of the cause of interest.
##' @param time Time point \eqn{t} for the analysis. Required.
##' @param beta Starting values for optimization (default NULL, uses zeros).
##' @param type Type of estimator: \code{"I"} (IPCW only), \code{"II"} (IPCW + augmentation), 
##'   or \code{"III"} (IPCW + complex augmentation). Default is \code{"III"}.
##' @param offset Offsets for the partial likelihood.
##' @param weights Weights for the score equations.
##' @param cens.weights External censoring weights (if provided, \code{cens.model} is ignored).
##' @param cens.model Formula for the censoring model (default \code{~+1}, stratified KM). 
##'   Can include \code{strata()} for stratified censoring.
##' @param se Logical; if TRUE, computes standard errors based on IPCW (default TRUE).
##' @param relative.to.causes If not NULL, compares the RMTL of the specified \code{cause} 
##'   to the RMTL of the causes in this vector (the denominator becomes the sum of these causes).
##' @param kaplan.meier Logical; if TRUE, uses Kaplan-Meier for IPCW weights; if FALSE, 
##'   uses \eqn{\exp(-\text{cumulative hazard})}.
##' @param cens.code Censoring code (default 0).
##' @param no.opt Logical; if TRUE, skips optimization and uses \code{beta} directly.
##' @param method Optimization method: \code{"nr"} (Newton-Raphson) or \code{"nlm"}.
##' @param augmentation Initial augmentation term (used for type "II" and "III").
##' @param outcome Outcome type: \code{"rmtl"} (years lost) or \code{"cif"} (cumulative incidence).
##' @param model Link function: \code{"logit"} (default), \code{"exp"}, or \code{"lin"}.
##' @param Ydirect Matrix with two columns (numerator, denominator) to use directly as the response.
##' @param ... Additional arguments passed to lower-level functions.
##' @return An object of class \code{"binreg"} and \code{"ratio"} containing:
##'   \item{coef}{Coefficient estimates.}
##'   \item{se.coef}{Standard errors.}
##'   \item{var}{Variance-covariance matrix.}
##'   \item{iid}{Influence function decomposition (with censoring adjustment).}
##'   \item{iidI}{Influence function without censoring adjustment.}
##'   \item{naive.var}{Variance assuming known censoring.}
##'   \item{time}{Time point used.}
##'   \item{cause}{Cause of interest.}
##'   \item{Causes}{Set of causes considered in the denominator.}
##'   \item{Yipcw}{IPCW-adjusted response matrix.}
##'   \item{coefI, varI}{Results from the initial (type "I") fit.}
##'   \item{augmentation}{Final augmentation term used.}
##' @author Thomas Scheike
##' @references 
##' Scheike, T. & Tanaka, S. (2025). Restricted mean time lost ratio regression: Percentage of restricted mean time lost due to specific cause. WIP.
##' @seealso \code{\link{resmeanIPCW}}, \code{\link{binreg}}
##' @examples
##' data(bmt); bmt$time <- bmt$time+runif(408)*0.001
##' 
##' rmtl30 <- rmstIPCW(Event(time,cause!=0)~platelet+tcell+age, bmt, time=30, cause=1, outcome="rmtl")
##' rmtl301 <- rmstIPCW(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=1)
##' rmtl302 <- rmstIPCW(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=2)
##' 
##' estimate(rmtl30)
##' estimate(rmtl301)
##' estimate(rmtl302)
##' 
##' ## Percentage of total RMTL due to cause 1
##' rmtlratioI <- rmtlRatio(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=1)
##' summary(rmtlratioI)
##' 
##' newdata <- data.frame(platelet=1, tcell=1, age=1)
##' pp <- predict(rmtlratioI, newdata)
##' pp
##' 
##' ## Percentage of total cumulative incidence due to cause 1
##' cifratio <- binregRatio(Event(time,cause)~platelet+tcell+age, bmt, time=30, cause=1, model="cif")
##' summary(cifratio)
##' pp <- predict(cifratio, newdata)
##' pp
##' @aliases rmtlRatio
##' @export
binregRatio <- function(formula,data,cause=1,time=NULL,beta=NULL,type=c("III","II","I"),
	   offset=NULL,weights=NULL,cens.weights=NULL,cens.model=~+1,se=TRUE,relative.to.causes=NULL,
	   kaplan.meier=TRUE,cens.code=0,no.opt=FALSE,method="nr",augmentation=NULL,
	   outcome=c("rmtl","cif"),model=c("logit","exp","lin"),Ydirect=NULL,...)
{# {{{
  cl <- match.call()# {{{
    m <- match.call(expand.dots = TRUE)[1:3]
    des <- proc_design(
        formula, data = data, specials = c("offset", "weights", "cluster"),
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
    id      <- des$cluster

 call.id <- id;
 conid <- construct_id(id,nrow(X))
 name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
 orig.id <- id

  ## take offset and weight first from formula, but then from arguments
  if (is.null(des.offset)) {
	  if (is.null(offset)) offset <- rep(0,length(exit)) 
  } else offset <- des.offset
  if (is.null(des.weights)) {
	  if (is.null(weights)) weights <- rep(1,length(exit)) 
  } else weights <- des.weights
# }}}

  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
  statusC <- 1*(status %in% cens.code) 
  statusE <- (status %in% cause) & (exit<= time) 
  if (sum(statusE)==0) warning("No events of type 1 before time \n"); 
  kmt <- kaplan.meier

  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
  competing  <-  (length(Causes)>1) 
  if (!is.null(relative.to.causes)) Causes <- relative.to.causes
  data$id__ <- id
  data$exit <- exit
  data$statusC <- statusC 
  cens.strata <- cens.nstrata <- NULL 

 nevent <- sum((status %in% cause)*(exit<=time))
 ## if event before time or alive, then uncensored, equality for both censored and events  
 obs <- (exit<=time & (!statusC)) | (exit>=time)

  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else { se <- FALSE; resC <- formC <- NULL}
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  p <- ncol(X)
  if (is.null(beta)) beta <- rep(0,p)
  if (is.null(augmentation))  augmentation=rep(0,p)
  X <-  as.matrix(X)
  X2  <- .Call("vecCPMat",X)$XX

 if (!is.null(Ydirect)) Y <-  Ydirect*obs/cens.weights else {
     if (outcome[1]=="cif") Y <- cbind( c((status %in% Causes)*(exit<=time)/cens.weights) , c((status %in% cause)*(exit<=time)/cens.weights) )
     else { 
            Y <- cbind(c((status %in% Causes)*(time-pmin(exit,time))*obs)/cens.weights, c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights)
     }
  }
  Yipcw <- Y
  cens.weights.origsort <- obs/cens.weights 

obj <- function(pp,all=FALSE)
{ # {{{
lp <- c(X %*% pp+offset)

if (model[1]=="logit") {
	p <- expit(lp)
        D2logl <- c(Y[,1]*weights*p/(1+exp(lp)))*X2
} else if (model[1]=="exp") {
	p <- exp(lp)
        D2logl <- c(Y[,1]*p*weights)*X2
} else {
	p <- lp
        D2logl <- c(Y[,1]*weights)*X2
}
Dlogl <- weights*X*c(Y[,1]*p-Y[,2])
D2log <- apply(D2logl,2,sum)
ploglik <- 0   ###sum(weights*(Y[,1]*p-Y[,2])^2)

gradient <- apply(Dlogl,2,sum)+augmentation
np <- length(pp)
hessian <- matrix(.Call("XXMatFULL",matrix(D2log,nrow=1),np,PACKAGE="mets")$XXf,np,np)

  if (all) {
      ploglik <- sum(weights*(Y[,1]*p-Y[,2])^2)
      ihess <- pinv(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      beta.iid <-  apply(beta.iid,2,sumstrata,id,max(id)+1)
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
	 id=id,Dlogl=Dlogl,iid=beta.iid,robvar=robvar,var=robvar,se.robust=diag(robvar)^.5)
      return(val)
  }  
 structure(0,gradient=-gradient/nid,hessian=-hessian/nid)
}# }}}

 if (model[1]=="exp") control <- list(stepsize=0.5)  else control <- NULL
  

  ## first run without pseudo-value augmentation
  p <- ncol(X)
  opt <- NULL
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
	  tim <- system.time(opt <- lava::NR(beta,obj,control=control))
	  opt$timing <- tim
	  opt$estimate <- opt$par
      } else {
	  opt <- nlm(obj,beta,...)
	  opt$method <- "nlm"
      }
      cc <- opt$estimate; 
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
  } else val <- c(list(coef=beta),obj(beta,all=TRUE))

 coefI <- val$coef
 gradientI <- val$gradient
 hessianI <- val$hessian
 ihessianI <- val$ihessian
 iidI <- val$iid

 if (se) {## {{{ censoring adjustment of variance 
    ### order of sorted times
    ord <- resC$ord
    X <-  X[ord,,drop=FALSE]
    X2 <-  X2[ord,,drop=FALSE]
    status <- status[ord]
    exit <- exit[ord]
    weights <- weights[ord]
    offset <- offset[ord]
    cens.weights <- cens.weights[ord]
    Y <- Y[ord,]
    lp <- c(X %*% val$coef+offset)
    if (model[1]=="logit")  p <- expit(lp) else if (model[1]=="exp") p <- exp(lp)  else p <- lp
    Yo <- Y[,1]*p-Y[,2]
    id <- id[ord]

    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ## compute function h(s) = \sum_i X_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  apply(X*Yo,2,revcumsumstrata,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    btime <- 1*(exit<time)
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    IhdLam0 <- apply(h*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*h[xx$jumps+1,]/c(resC$S0)
    MGt <- (U[,drop=FALSE]-IhdLam0)*c(xx$weights)

    ### Censoring Variance Adjustment  \int h^2(s) / y.(s) d Lam_c(s) estimated by \int h^2(s) / y.(s)^2  d N.^C(s) 
    mid <- max(id)
    MGCiidI <- MGCiid <- apply(MGt,2,sumstrata,xx$id,mid+1)

   if (type[1]=="II") { ##  pseudo-value type augmentation
    hYt  <-  revcumsumstrata(Yo,xx$strata,xx$nstrata)
    IhdLam0 <- cumsumstrata(hYt*S0i2*btime,xx$strata,xx$nstrata)
    U <- rep(0,length(xx$strata))
    U[xx$jumps+1] <- (resC$jumptimes<time)*hYt[xx$jumps+1]/c(resC$S0)
    MGt <- X*c(U-IhdLam0)*c(xx$weights)
    MGtiid <- apply(MGt,2,sumstrata,xx$id,mid+1)
    augmentation  <-  apply(MGtiid,2,sum) + augmentation
    ###
    EXt  <-  apply(X,2,revcumsumstrata,xx$strata,xx$nstrata)
    IEXhYtdLam0 <- apply(EXt*c(hYt)*S0i*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*hYt[xx$jumps+1]*EXt[xx$jumps+1,]/c(resC$S0)^2
    MGt2 <- (U[,drop=FALSE]-IEXhYtdLam0)*c(xx$weights)
    ###
    MGCiid2 <- apply(MGt2,2,sumstrata,xx$id,mid+1)
    ### Censoring Variance Adjustment 
    MGCiid <- MGCiid+(MGtiid-MGCiid2)
   }
   if (type[1]=="III") { ##  pseudo-value type augmentation
    hYt  <-  apply(Y,2,revcumsumstrata,xx$strata,xx$nstrata)
    IhdLam0 <- apply(hYt*c(S0i2)*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- rep(0,length(xx$strata))
    U[xx$jumps+1] <- (resC$jumptimes<time)*(hYt[xx$jumps+1,1]*p[xx$jumps+1]-hYt[xx$jumps+1,2])/c(resC$S0)
    MGt <- X*c(U-IhdLam0[,1]*p+IhdLam0[,2])*c(xx$weights)

    MGtiid <- apply(MGt,2,sumstrata,xx$id,mid+1)
    augmentation  <-  apply(MGtiid,2,sum) + augmentation
    ###
    EXt  <-  apply(X,2,revcumsumstrata,xx$strata,xx$nstrata)
    EXFt  <-  apply(X*p,2,revcumsumstrata,xx$strata,xx$nstrata)

    IEXhY1tdLam0 <- apply(EXt*hYt[,2]*S0i*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    IEXhYtdLam0 <- apply(EXFt*hYt[,1]*S0i*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)

    U <- matrix(0,nrow(xx$X),ncol(X))
    U[xx$jumps+1,] <- (resC$jumptimes<time)*(hYt[xx$jumps+1,1]*EXFt[xx$jumps+1,]-hYt[xx$jumps+1,2]*EXt[xx$jumps+1,]) /c(resC$S0)^2
    MGt2 <- (U[,drop=FALSE]-IEXhYtdLam0+IEXhY1tdLam0)*c(xx$weights)
    ###
    MGCiid2 <- apply(MGt2,2,sumstrata,xx$id,mid+1)
    ### Censoring Variance Adjustment 
    MGCiid <- MGCiid+(MGtiid-MGCiid2)
   }
   ## use data ordered by time (keeping track of id also)
   ## since X; Y and so forth are ordered in time
   ###   id <- xx$id
   }  else {
	 MGCiidI <-  MGCiid <- 0
  }## }}}

  ## first run without pseudo-value augmentation, then run with augmentation, if type="II"
 if (type[1]!="I") {
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
	  tim <- system.time(opt <- lava::NR(coefI,obj,...))
	  opt$timing <- tim
	  opt$estimate <- opt$par
      } else {
	  opt <- nlm(obj,beta,...)
	  opt$method <- "nlm"
      }
      cc <- opt$estimate; 
###	      if (!se) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
 }

 if (length(val$coef)==length(colnames(X))) names(val$coef) <- colnames(X)
  val <- c(val,list(time=time,formula=formula,formC=formC,
    exit=exit, cens.weights=cens.weights, cens.strata=cens.strata, cens.nstrata=cens.nstrata, 
    n=length(exit),nevent=nevent,ncluster=nid))


  val$call <- cl
  val$MGciid <- MGCiid
  val$id <- orig.id
  val$call.id <- call.id
  val$name.id <- name.id
  val$nid <- nid
  val$iid.naive <- val$iid
  val$naive.var <- NULL 
  val$iidI <- iidI
  if (se) val$iid  <- val$iid+(MGCiid %*% val$ihessian)
  if (se) val$iidI  <- iidI+(MGCiidI %*% ihessianI)
  if (!is.null(call.id)) {
  val$iid <- nameme(val$iid,name.id)
  val$iidI <- nameme(val$iidI,name.id)
  }
  robvar <- crossprod(val$iid)
  val$var <-  val$robvar <- robvar
  val$se.robust <- diag(robvar)^.5
  val$se.coef <- diag(val$var)^.5
  val$cause <- cause
  val$cens.code <- cens.code 
  val$augmentation <- augmentation
  val$model <- model[1]
  val$outcome <- outcome[1]
  val$Yipcw <- Yipcw
  val$cenw.weights <- cens.weights.origsort 
  val$Causes <- Causes
  val$nevent <- nevent

  val$coefI <- coefI
  val$varI <- crossprod(val$iidI)
  val$se.coefI <- diag(val$varI)^.5
  val$gradientI <-  gradientI
  val$hessianI <- hessianI
  val$ihessianI <- ihessianI
  val$design <- des

  class(val) <- c("binreg","ratio")
  return(val)
}# }}}

##' @export
rmtlRatio <- function(formula,data,outcome=c("rmtl"),...)
{# {{{
   out <- binregRatio(formula,data,outcome=outcome[1],...)
   return(out)
}# }}}

