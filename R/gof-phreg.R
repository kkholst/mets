##' Goodness-of-Fit for Cox PH Regression (Proportionality)
##'
##' Performs cumulative score process residual tests for the proportional hazards (PH) 
##' assumption in Cox regression. The test statistics are based on the cumulative 
##' score process:
##' \deqn{ U(t) = \int_0^t (X_i - E(t) ) d \hat M_i(s) }
##' where \eqn{\hat M_i(s)} are the martingale residuals.
##' 
##' P-values are computed using the Lin, Wei, and Ying (1993) resampling method, which 
##' simulates the asymptotic distribution of the supremum of the score process under 
##' the null hypothesis of proportional hazards.
##'
##' The function supports two types of simulation:
##' \itemize{
##'   \item \strong{Standard}: Uses \eqn{dN_i} (counting process increments) for simulation.
##'   \item \strong{Robust}: Uses \eqn{\hat M_i(t)} (martingale residuals) adjusted for clustering 
##'     if a \code{cluster()} term is present in the model.
##' }
##'
##' @param object A fitted \code{phreg} object (from \code{mets} or \code{survival}).
##' @param n.sim Number of simulations for the resampling procedure (default 1000).
##' @param silent Logical; if TRUE, suppresses timing estimates for long jobs.
##' @param robust Logical; if TRUE, uses robust martingale-based simulations. 
##'   If NULL, defaults to TRUE if a cluster term is detected in the model call.
##' @param ... Additional arguments passed to lower-level functions.
##' @return An object of class \code{"gof.phreg"} containing:
##'   \item{jumptimes}{Event times used in the process.}
##'   \item{supUsim}{Matrix of simulated supremum values for each covariate.}
##'   \item{res}{Matrix with observed supremum (\code{Sup|U(t)|}) and p-values.}
##'   \item{supU}{Observed supremum values.}
##'   \item{pvals}{Vector of p-values for each covariate.}
##'   \item{score}{Cumulative score process values over time.}
##'   \item{simUt}{Simulated score processes.}
##'   \item{type}{Type of test performed ("prop").}
##'   \item{robust}{Logical flag indicating if robust simulation was used.}
##' @author Thomas Scheike and Klaus K. Holst
##' @references 
##' Lin, D. Y., Wei, L. J., & Ying, Z. (1993). Checking the Cox model with cumulative sums of martingale-based residuals. Biometrika, 80(3), 557-572.
##' @seealso \code{\link{gofM_phreg}}, \code{\link{gofZ_phreg}}
##' @examples
##' data(sTRACE)
##' 
##' m1 <- phreg(Surv(time,status==9)~vf+chf+diabetes, data=sTRACE) 
##' gg <- gof(m1)
##' gg
##' par(mfrow=c(1,3))
##' plot(gg)
##' 
##' m1 <- phreg(Surv(time,status==9)~strata(vf)+chf+diabetes, data=sTRACE) 
##' gg <- gof(m1)
##' 
##' ## Robust simulations with cluster
##' sTRACE$id <- 1:500
##' m1 <- phreg(Surv(time,status==9)~vf+chf+diabetes+cluster(id), data=sTRACE) 
##' gg <- gof(m1)
##' gg
##' @export
gof.phreg  <- function(object,n.sim=1000,silent=1,robust=NULL,...)
{# {{{

### test for proportionality 
p <- length(object$coef)
if (p==0) stop("Goodness of fit for proportionality, must contain regression effects\n");
nnames <- names(object$coef)
ii <- pinv(object$hessian)
jumptimes <- object$jumptimes
Pt <- object$hessianttime
U <-  object$U
Pt <- apply(Pt,2,cumsum)
Ut <- apply(U,2,cumsum)

nd <- nrow(object$U)
Pt <-  .Call("XXMatFULL",Pt,p,PACKAGE="mets")$XXf
Pt <-  .Call("CubeMat",Pt,ii,PACKAGE="mets")$XXX
sup <- matrix(0,n.sim,nrow(ii))
hatti <- matrix(0,nd,nrow(ii))
obs <- apply(abs(Ut),2,max)

if (is.null(robust)) {
   if (!is.null(object$call.id)) robust <- TRUE else robust <- FALSE
   if (inherits(object,c("cifreg"))) {
	   robust <- FALSE  ## approximative test using MG framework
###	   if (is.null(object$cox.prep)) stop("Must be called with cox.prep=TRUE\n")
   }
}


   ### cluster call or robust \hat M_i(t) based  
if (robust) {# {{{
   xx <- object$cox.prep
   Z <- xx$X
   rrw <- c(xx$sign * exp(Z %*% coef(object) + xx$offset)*xx$weights)
   UdN <- xx$weights[xx$jumps+1]*object$U 
   ### also weights 
   nn <- nrow(Z)

   tt <- system.time(simcox1<- .Call("PropTestCoxClust",UdN,Pt,rrw,xx$X,
     object$S0,object$E,
     10,obs,nn,xx$id,xx$strata,xx$nstrata,object$strata.jumps,xx$jumps))

   prt <- n.sim*tt[3]/(10*60)
   if (prt>1 & silent==0) cat(paste("Predicted time minutes",signif(prt,2),"\n"))
   simcox <- .Call("PropTestCoxClust",UdN,Pt,rrw,xx$X,
	     object$S0,object$E,
	     n.sim,obs,nn,xx$id,xx$strata,xx$nstrata,object$strata.jumps,xx$jumps)
} else {# }}}
###  or dN_i based  # {{{
   tt <- system.time(simcox1<-.Call("PropTestCox",U,Pt,10,obs,PACKAGE="mets"))
   prt <- n.sim*tt[3]/(10*60)
   if (prt>1 & silent==0) cat(paste("Predicted time minutes",signif(prt,2),"\n"))
   simcox <-  .Call("PropTestCox",U,Pt,n.sim,obs,PACKAGE="mets")
}# }}}

sup <-  simcox$supUsim
res <- cbind(obs,simcox$pval)
colnames(res) <- c("Sup|U(t)|","pval")
rownames(res) <- nnames 

if (silent==0) {
cat("Cumulative score process test for Proportionality:\n")
prmatrix(round(res,digits=2))
}

out <- list(jumptimes=object$jumptimes,supUsim=sup,res=res,supU=obs,
    pvals=simcox$pval,score=Ut,simUt=simcox$simUt,type="prop",robust=robust)
class(out) <- "gof.phreg"
return(out)
}# }}}

##' Goodness-of-Fit for Cox Covariates (Model Matrix)
##'
##' Tests the functional form of covariates in a Cox PH model by computing cumulative 
##' residuals against a user-specified model matrix. This helps detect non-linear effects 
##' or time-varying coefficients (interaction with time).
##'
##' The test statistic is:
##' \deqn{ U(t) = \int_0^t K^T d \hat M }
##' where \eqn{K} is the model matrix (e.g., a set of basis functions for a continuous covariate) 
##' and \eqn{\hat M} are the martingale residuals.
##'
##' P-values are based on the Lin, Wei, and Ying (1993) resampling method. The plot shows 
##' whether the residuals are consistent with the model across the range of the covariate.
##'
##' @param formula Formula for the Cox regression model.
##' @param data Data frame.
##' @param offset Offset vector.
##' @param weights Weights vector.
##' @param modelmatrix Matrix for cumulating residuals. Typically constructed using 
##'   \code{cumContr()} or manually (e.g., quartiles of a continuous covariate).
##' @param n.sim Number of simulations (default 1000).
##' @param silent Logical; suppresses timing estimates if TRUE.
##' @param ... Additional arguments passed to \code{phreg}.
##' @return An object of class \code{"gof.phreg"} containing:
##'   \item{jumptimes}{Event times.}
##'   \item{supUsim}{Simulated supremum values.}
##'   \item{res}{Matrix with observed supremum and p-values for each column of \code{modelmatrix}.}
##'   \item{score}{Cumulative score process.}
##'   \item{simUt}{Simulated processes.}
##'   \item{Utlast, pval.last}{Supremum and p-value for the final time point (covariate direction).}
##'   \item{type}{Type of test ("modelmatrix").}
##' @author Thomas Scheike and Klaus K. Holst
##' @references 
##' Lin, D. Y., Wei, L. J., & Ying, Z. (1993). Checking the Cox model with cumulative sums of martingale-based residuals. Biometrika, 80(3), 557-572.
##' @seealso \code{\link{gof.phreg}}, \code{\link{gofZ_phreg}}, \code{\link{cumContr}}
##' @examples
##' data(TRACE)
##' set.seed(1)
##' TRACEsam <- blocksample(TRACE, idvar="id", replace=FALSE, 100)
##' dcut(TRACEsam) <- ~. 
##' mm <- model.matrix(~-1+factor(wmicat.4), data=TRACEsam)
##' m1 <- gofM_phreg(Surv(time,status==9)~vf+chf+wmi, data=TRACEsam, modelmatrix=mm)
##' summary(m1)
##' if (interactive()) {
##' par(mfrow=c(2,2))
##' plot(m1)
##' }
##' 
##' ## Cumulative sums in covariates via design matrix
##' mm <- mets:::cumContr(TRACEsam$wmi, breaks=10, equi=TRUE)
##' m1 <- gofM_phreg(Surv(time,status==9)~strata(vf)+chf+wmi, data=TRACEsam,
##'          modelmatrix=mm, silent=0)
##' summary(m1)
##' @export
gofM_phreg  <- function(formula,data,offset=NULL,weights=NULL,modelmatrix=NULL,
			n.sim=1000,silent=1,...)
{# {{{

if (is.null(modelmatrix)) stop(" must give matrix for cumulating residuals\n"); 

cox1 <- phreg(formula,data,offset=NULL,weights=NULL,Z=modelmatrix,cumhaz=FALSE,...) 
p <- length(cox1$coef)
if (p==0) stop("Goodness of fit, must contain regression effects\n");

## put modelmatrix on data, and take y, strata from design
datl <- as.data.frame(modelmatrix)
namesmm <- paste("names",1:ncol(modelmatrix),sep="")
names(datl) <- namesmm
datl[,"y__"] <- cox1$design$y
if (!is.null(cox1$design$strata)) {
datl[,"strata__"] <- cox1$design$strata
formM <- y__~strata(strata__)
formM <- update(formM, paste(". ~ . +",paste(namesmm, collapse = " + ")))
} else {
formM <- as.formula(paste("y__~",paste(namesmm, collapse = " + ")))
}

offsets <- cox1$X %*% cox1$coef
if (!is.null(offset)) offsets <- offsets*offset
coxM <- phreg(formM,datl,offset=offsets,weights=weights,no.opt=TRUE,cumhaz=FALSE,no.var=1,...)
nnames <- colnames(modelmatrix)

Ut <- apply(coxM$U,2,cumsum)
jumptimes <- coxM$jumptimes
U <- coxM$U
Ubeta <- cox1$U
ii <- -pinv(cox1$hessian)
EE <- .Call("vecMatMat",coxM$E,cox1$E,PACKAGE="mets")$vXZ; 
Pt <- cox1$ZX - EE
Pt <- apply(Pt,2,cumsum)
betaiid <- t(ii %*% t(Ubeta))
obs <- apply(abs(Ut),2,max)
simcox <-  .Call("ModelMatrixTestCox",U,Pt,betaiid,n.sim,obs,PACKAGE="mets")

sup <-  simcox$supUsim
res <- cbind(obs,simcox$pval)
colnames(res) <- c("Sup_t |U(t)|","pval")
rownames(res) <- nnames 

if (silent==0) {
   cat("Cumulative score process test for modelmatrix:\n")
   prmatrix(round(res,digits=2))
}

 ## pvals efter z i model.matrix sup_z | M(z,tau) | 
 Utlast <- max(abs(tail(Ut,1)))
 maxlast <- apply(abs(simcox$last),1,max)
 pval.last <- mean(maxlast>=Utlast)
 res.last <- matrix(c(Utlast,pval.last),1,2)
 colnames(res.last) <- c("Sup_z |U(tau,z)|","pval")
 rownames(res.last) <- "matrixZ"

out <- list(jumptimes=jumptimes,supUsim=simcox$supUsim,res=res,supU=obs,
	    pvals=simcox$pval,score=Ut,simUt=simcox$simUt,
	    simUtlast=simcox$last,Utlast=Utlast,pval.last=pval.last,
	    res.last=res.last, type="modelmatrix")
class(out) <- "gof.phreg"

return(out)
}# }}}

##' Goodness-of-Fit for Cox Covariates (Linearity)
##'
##' Tests the functional form of continuous covariates in a Cox PH model to check for 
##' linearity. It computes cumulative residuals evaluated at a grid of covariate values \eqn{z}.
##' 
##' The test statistic is:
##' \deqn{ U(z, \tau) = \int_0^\tau K(z)^T d \hat M }
##' where \eqn{K(z)} is a design matrix based on indicator functions \eqn{I(Z_i \leq z_l)} 
##' for a grid of points \eqn{z_l}.
##'
##' The p-value is valid but depends on the chosen grid. As the number of break points increases, 
##' this test converges to the original Lin, Wei, and Ying test for linearity.
##'
##' @param formula Formula for the Cox regression.
##' @param data Data frame.
##' @param vars Vector of variable names to test. If NULL, automatically detects continuous 
##'   covariates with more than 2 levels.
##' @param offset Offset vector.
##' @param weights Weights vector.
##' @param breaks Number of break points for the grid (default 50).
##' @param equi Logical; if TRUE, uses equidistant breaks; if FALSE, uses quantiles.
##' @param n.sim Number of simulations (default 1000).
##' @param silent Logical; suppresses timing estimates.
##' @param ... Additional arguments passed to \code{gofM_phreg}.
##' @return An object of class \code{"gof.phreg"} with type "Zmodelmatrix" containing:
##'   \item{res}{Matrix of p-values for each tested variable.}
##'   \item{Zres}{List of \code{gof.phreg} objects, one for each variable.}
##'   \item{type}{Type of test ("Zmodelmatrix").}
##' @author Thomas Scheike and Klaus K. Holst
##' @seealso \code{\link{gofM_phreg}}, \code{\link{cumContr}}
##' @examples
##' data(TRACE)
##' set.seed(1)
##' TRACEsam <- blocksample(TRACE, idvar="id", replace=FALSE, 100)
##' 
##' ## Test linearity of continuous covariates
##' \donttest{ ## Reduce Ex.Timings
##' m1 <- gofZ_phreg(Surv(time,status==9)~strata(vf)+chf+wmi+age, data=TRACEsam)
##' summary(m1) 
##' plot(m1, type="z")
##' }
##' @aliases cumContr 
##' @export
gofZ_phreg  <- function(formula,data,vars=NULL,offset=NULL,weights=NULL,breaks=50,equi=FALSE,
			n.sim=1000,silent=1,...)
{# {{{
if (is.null(vars)) {
     vars <- NULL
      ## find strata var's 
      yxzf <- procform(formula,x=NULL,z=NULL,data=data,do.filter=FALSE)
      avars <- all.vars(formula[-2])
      svar <- grep("strata",yxzf$predictor)
    if (length(svar)>=1) {
	avars <- avars[-svar]
     }
     ## check that it is not a factor and that there are more than 2 levels
     for (vv in avars) {
	  if (length(unique(data[,vv]))>2 & !is.factor(data[,vv]))
	  vars  <-  c(vars,vv)
     }
 }

 if (length(vars)==0) stop("Checking functional form of covariates, must contain regression effects of continous covariates\n")
 gres <- list()
 res <- matrix(0,length(vars),2)
 colnames(res) <- c("Sup_z |U(tau,z)|","pval")
 rownames(res) <- vars

i <- 1
for (vv in vars) {
 modelmatrix <- cumContr(data[,vv],breaks=breaks,equi=equi)
 lres <- gofM_phreg(formula,data,modelmatrix=modelmatrix) 
 lres$xaxs <- attr(modelmatrix,"breaks")

 res[i,] <- c(lres$Utlast,lres$pval.last)
 i <- i+1

 lres <- list(lres)
 names(lres) <- vv
 gres <- c(gres,lres)
}

out <- list(res=res,Zres=gres,type="Zmodelmatrix")
class(out) <- c("gof.phreg")

return(out)
}# }}}

cumContr <- function(data,breaks=4,probs=NULL,equi=TRUE,na.rm=TRUE,unique.breaks=TRUE,...)
 {# {{{
 if (is.vector(data)) {
        if (is.list(breaks))
            breaks <- unlist(breaks)
        if (length(breaks) == 1) {
            if (!is.null(probs)) {
                breaks <- quantile(data, probs, na.rm = na.rm,...)
	        breaks <- breaks[-1]
            }
            else {
                if (!equi) {
                  probs <- seq(0, 1, length.out = breaks + 1)
                  breaks <- quantile(data, probs, na.rm = na.rm, ...)
	          if (unique.breaks) breaks <- unique(breaks)
	          breaks <- breaks[-1]
                }
                if (equi) {
                  rr <- range(data, na.rm = na.rm)
                  breaks <- seq(rr[1], rr[2], length.out = breaks + 1)
	          breaks <- breaks[-1]
                }
            }
        }
        if (sum(duplicated(breaks)) == 0) {
		i <- 0; 
		gm <- matrix(0,length(data),length(breaks))
            for (bb in breaks)  {
		    i <- i+1 
		    gm[,i] <- (data <= bb)*1
	    }
	} else {
            wd <- which(duplicated(breaks))
            mb <- min(diff(breaks[-wd]))
            breaks[wd] <- breaks[wd] + (mb/2) * seq(length(wd))/length(wd)
            i <- 0; gm <- matrix(0,length(data),length(breaks))
            for (bb in breaks)  {
		    i <- i+1 
		    gm[,i] <- (data <= bb)*1
	    }
            warning(paste("breaks duplicated"))
        }
        colnames(gm) <- paste("<=",breaks,sep="")
	attr(gm,"breaks") <- breaks
        return(gm)
    }
 }# }}}

##' @export
plot.gof.phreg <-  function(x,col=3,type=NULL,...)
{# {{{

if (is.null(type)) {
	if (x$type=="prop") type <- "time"
	if (x$type=="modelmatrix" ) type <- "modelmatrix"
	if (x$type=="Zmodelmatrix") type <- "z"
}

if (type=="time" || type=="modelmatrix") {
p <- ncol(x$score)
for (i in 1:p)
{
simU <- x$simUt[,(0:49)*p+i]
rsU <- max(abs(simU))
rsU <- max(rsU,abs(x$score[,i]))
plot(x$jumptimes,x$score[,i],type="s",ylim=c(-rsU,rsU),xlab="",ylab="")
title(main=rownames(x$res)[i])
matlines(x$jumptimes,simU,type="s",lwd=0.3,col=col)
lines(x$jumptimes,x$score[,i],type="s",lwd=1.5)
}
} else {
	if (type=="modelmatrix") {
		obsz <- c(tail(x$score,1))
		times <- 1:length(obsz)
		rsU <- max(max(abs(obsz)),max(abs(x$simUtlast[1:50,])))
		plot(times,obsz,type="l",ylim=c(-rsU,rsU),xlab="",ylab="")
		matlines(times,t(x$simUtlast[1:50,]),type="l",lwd=0.3,col=col)
	        ## redraw with thick to make observed clear 
	        lines(times,obsz,lwd=2,col=1)
	} else {
	   for (i in 1:length(x$Zres))
	   {
	    xr <- x$Zres[[i]]
	    obsz <- c(tail(xr$score,1))
	    times <- xr$xaxs
	    rsU <- max(max(abs(obsz)),max(abs(xr$simUtlast[1:50,])))
	    plot(times,obsz,type="l",ylim=c(-rsU,rsU),xlab="",ylab="")
	    title(main=rownames(x$res)[i])
	    matlines(times,t(xr$simUtlast[1:50,]),type="l",lwd=0.3,col=col)
	    ## redraw with thick to make observed clear 
	    lines(times,obsz,lwd=2,col=1)
           }
	}
}

}# }}}

##' @export
summary.gof.phreg <-  function(object,...)
{# {{{
if (object$type=="prop")
     cat("Cumulative score process test for Proportionality:\n")
else cat("Cumulative residuals versus modelmatrix :\n")
print(object$res)

if (!is.null(object$res.last)) {
   cat("\n")
   cat("Cumulative score process versus covariates (discrete z via model.matrix):\n")
   print(object$res.last)
}

} # }}}

##' @export
print.gof.phreg <-  function(x,...)
{# {{{
if (x$type=="prop")
     cat("Cumulative score process test for Proportionality:\n")
else cat("Cumulative residuals versus modelmatrix :\n")
print(x$res)

if (!is.null(x$res.last)) {
   cat("\n")
   cat("Cumulative score process versus covariates (discrete z via model.matrix):\n")
   print(x$res.last)
}

} # }}}

