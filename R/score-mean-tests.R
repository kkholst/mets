##' Pepe-Mori Test for Marginal Mean Comparison
##'
##' Performs score-test type tests for proportionality of marginal means in competing 
##' risks data and recurrent events, as presented in Ghosh and Lin (2000). The test 
##' is based on an IPCW (Inverse Probability of Censoring Weighting) formulation.
##' 
##' The function computes several tests:
##' \enumerate{
##'   \item \strong{Pepe-Mori Test}: Tests for equality of marginal mean functions between groups.
##'   \item \strong{Ratio of AUC}: Compares the area under the curve of marginal means.
##'   \item \strong{Difference of AUC}: Tests for difference in areas under the curve.
##'   \item \strong{Score Test}: Tests for proportionality (equivalent to Gray's test for CIF).
##'   \item \strong{Proportionality Test}: Tests the proportional hazards assumption.
##' }
##' 
##' The Pepe-Mori test uses weights based on the number at risk in each group to 
##' construct a weighted integral of the difference in marginal means.
##'
##' @param formula Formula with an \code{Event} object on the left-hand side and 
##'   covariates (typically with \code{strata()} for group comparison) on the right. 
##'   Can include \code{cluster(id)} for correlated data.
##' @param data Data frame containing all variables referenced in the formula.
##' @param cause Cause of interest (default 1).
##' @param cens.code Censoring code (default 0).
##' @param death.code Code for death (terminating event, default 2).
##' @param death.code.prop Code for other causes of death for Fine-Gray regression model.
##' @param time Upper limit for Pepe-Mori and AUC integrals. If NULL, defaults to 
##'   the maximum event time for the cause of interest.
##' @param beta Starting values for the score test (default NULL, uses zeros).
##' @param ... Additional arguments passed to lower-level functions.
##' @return An object of class \code{"marginalTest"} containing:
##'   \item{pepe.mori}{Pepe-Mori test results with compare p-value.}
##'   \item{RatioAUC}{Ratio of AUC test results with compare p-value.}
##'   \item{difAUC}{Difference of AUC test results with compare p-value.}
##'   \item{prop.test}{Proportionality test results.}
##'   \item{score.test}{Score test results (equivalent to Gray's test).}
##'   \item{score.iid}{Influence function for the score test.}
##'   \item{time}{Upper time limit used.}
##'   \item{RAUCl, RAUCe}{Raw and transformed AUC estimates.}
##' @author Thomas Scheike
##' @references 
##' Ghosh, D. and Lin, D. Y. (2000). Nonparametric Analysis of Recurrent Events and Death. Biometrics, 56, 554--562.
##' @seealso \code{\link{test_logrankRecurrent}}, \code{\link{test_conc}}
##' @examples
##' data(bmt,package="mets")
##' bmt$time <- bmt$time+runif(nrow(bmt))*0.01
##' bmt$id <- 1:nrow(bmt)
##' dcut(bmt) <- age.f~age
##'      
##' fg=cifregFG(Event(time,cause)~tcell,data=bmt,cause=1)
##'  
##' ## computing tests for difference  for CIF
##' pmt <- test_marginalMean(Event(time,cause)~strata(tcell)+cluster(id),data=bmt,cause=1,
##' 			 death.code=1:2,death.code.prop=2,cens.code=0,time=40)
##' summary(pmt) 
##'  
##' pmt$pepe.mori
##' pmt$RatioAUC
##' pmt$prop.test
##' ## score test equialent to Gray's test but variance estimated differently 
##' pmt$score.test
##'  
##' ### age-groups  
##' pmt <- test_marginalMean(Event(time,cause)~strata(age.f)+cluster(id),data=bmt,cause=1,
##' 			 death.code=1:2,death.code.prop=2,cens.code=0)
##' summary(pmt) 
##'  
##' ## having a look at the cumulative incidences 
##' cifs <- cif(Event(time,cause)~strata(age.f)+cluster(id),data=bmt,cause=1)
##' plot(cifs) 
##'  
##' ## recurrent events   
##' data(hfactioncpx12)
##' hf <- hfactioncpx12
##' pmt <- test_marginalMean(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,
##' 			 cause=1,death.code=2,cens.code=0)
##' summary(pmt) 
##'  
##' @export
test_marginalMean <- function(formula,data,cause=1,cens.code=0,...,death.code=2,death.code.prop=NULL,time=NULL,beta=NULL) { ## {{{ 
cl <- match.call()
m <- match.call(expand.dots = TRUE)[1:3]
des <- proc_design(formula, data = data, specials = c("offset","weights","cluster","strata","marks"), intercept = FALSE)
Y <- des$y
if (!inherits(Y, c("Event", "Surv"))) {
stop("Expected a 'Surv' or 'Event'-object")
}
if (ncol(Y) == 2) {
exit <- Y[, 1]
entry <- rep(0, nrow(Y))
status <- Y[, 2]
}
else {
entry <- Y[, 1]
exit <- Y[, 2]
status <- Y[, 3]
}
X <- des$x
strata <- des$strata
specials = c("offset", "weights", "cluster", "strata")
Terms <- terms(formula, specials, data = data)
ts <- survival::untangle.specials(Terms, "strata")
if (!is.null(strata))
strata.name <- ts$vars
else strata.name <- NULL
pos.strata <- NULL
des.weights <- des$weights
des.offset <- des$offset
id <- des$cluster
pos.cluster <- NULL

if (is.null(death.code.prop)) death.code.prop <- death.code

data$strata__ <- strata
data$entry__ <- entry
data$exit__ <- exit
nstrata <- as.numeric(strata)

###drop.specials(formula,"strata")
formR <- as.formula(paste("Event(entry__,exit__,status__)~factor(strata__)+cluster(id__)"))
data$status__  <- status
data$id__  <- id

## proportionality test  
proptest <- recreg(formR,data,cause=cause,death.code=death.code.prop,cens.model=~strata(strata__))

nlev <- nlevels(strata)
if (is.null(beta)) beta <- rep(0,nlev-1)
## score test  for prop-factor 
proptest0 <- recreg(formR,data,cause=cause,death.code=death.code.prop,
		   cens.model=~strata(strata__),beta=beta,no.opt=TRUE)
gradient <- matrix(proptest0$gradient,length(proptest0$gradient),1)
iidbeta <-  IC(proptest0)
n <- nrow(iidbeta)
iidU0 <- iidbeta %*% proptest0$hessian
score.test <- estimate(coef=gradient,IC=iidU0,null=0)

### varU0 <- crossprod(iidU0)
### logrank.robust <- t(gradient) %*% solve(varU0) %*% gradient
### p.logrank.robust <- 1-pchisq(logrank.robust,nrow(gradient))
### score.test <- list(logrank.robust=logrank.robust,p.logrank.robust=p.logrank.robust,
###		    test.statistic=gradient,iid=iidU0)

if (is.null(time)) time <- max(exit[status %in% cause])

## {{{ pepe-mori test 

ddataid <- data.frame(id=id)
cid <- countID(ddataid)
###ns <- c(table(strata[cid$reverse==1]))

data$statusC__ <- 1*(status %in% cens.code) 
data$statusD__ <- 1*(status %in% death.code)
formD <- as.formula(paste("Event(entry__,exit__,statusD__)~+1"))

if (nlev==2) {
	formC <- as.formula(paste("Event(entry__,exit__,statusC__)~strata(strata__)"))
	kms <- km(formC,data)
	kms0 <- phreg(formC,data,no.opt=TRUE)
	kms0$jumptimes
	kms0$strata.jumps
        Yr <- vecAllStrata(kms0$S0,kms0$strata.jumps,kms0$nstrata)
	minRisk <- apply(Yr,1,min)

        if (is.null(time)) {
               kms0 <- phreg(formC,data,no.opt=TRUE)
	       kms0$jumptimes
	       kms0$strata.jumps
               Yr <- vecAllStrata(kms0$S0,kms0$strata.jumps,kms0$nstrata)
	       minRisk <- apply(Yr,1,min)
               time <- max(kms0$jumptimes[minRisk>=5])
        }

        ns <- c(table(strata[cid$reverse==1]))
	n1 <- ns[1]
	n2 <- ns[2]
	## using number at risk for weighting 
	kmss <- cbind(kms$time,t(kms$surv))
	wt <- (n1+n2)*kmss[,2]*kmss[,3]/(n1*kmss[,2]+n2*kmss[,3])
	wt[is.na(wt)] <- 0
	Wt <- cumsum(diff(c(0,kms$time))*wt)

	Wtmark <- lin_approx(pmin(exit,time),rbind(0,cbind(kms$time,Wt)))
	Wfinal <- tail(Wt,1)
	###
	data$pmmark__ <- Wfinal-Wtmark

	pepe.mori <- recregIPCW(formR,data,cause=cause,death.code=death.code,
	    times=time,marks=data$pmmark__,cens.model=~strata(strata__),model="lin")
        pepe.mori <- estimate(pepe.mori,null=0)
	pepe.mori <- estimate(pepe.mori,lava::contr(2))
} else {

###  formRC <- as.formula(paste("Event(entry__,exit__,status__)~factor(cstrata__)+cluster(id__)"))
###	contr <- contr.iid <- c()
###        for (i in 1:(nlev-1)) {
###
###   	   data$cstrata__  <-  1*(nstrata==i)
###           formC <- as.formula(paste("Event(entry__,exit__,statusC__)~strata(cstrata__)"))
###
###           kms <- km(formC,data)
###           kmss <- cbind(kms$time,t(kms$surv))
###
###           ns <- c(table(data$cstrata__[cid$reverse==1]))
###	   n1 <- ns[1]
###	   n2 <- ns[2]
###	   ## using number at risk for weighting 
###	   kmss <- cbind(kms$time,t(kms$surv))
###	   wt <- (n1+n2)*kmss[,2]*kmss[,3]/(n1*kmss[,2]+n2*kmss[,3])
###	   wt[is.na(wt)] <- 0
###	   Wt <- cumsum(diff(c(0,kms$time))*wt)
###	   Wtmark <- lin_approx(pmin(exit,time),rbind(0,cbind(kms$time,Wt)))
###	   Wfinal <- tail(Wt,1)
###	   ###
###	   data$pmmark__ <- Wfinal-Wtmark
###
###	   pmOut <- recregIPCW(formRC,data,cause=cause,death.code=death.code,
###	       times=time,marks=data$pmmark__,cens.model=~strata(cstrata__),model="lin")
###
###        contr <- c(contr,coef(pmOut)[2])
###	contr.iid <- cbind(contr.iid,pmOut$iid[,2])
###   }
###   var.contr <- crossprod(contr.iid)
###   ###  pepe.mori <- estimate(coef=contr,vcov=var.contr,null=0)
###   pepe.mori <- estimate(coef=contr,IC=contr.iid*nrow(contr.iid),null=0)
###   p <- length(contr)
###   pepe.mori <- estimate(pepe.mori,lava::contr(1:p))
	pepe.mori <- list()
	pepe.mori$error <- "only for two levels"
        pepe.mori$compare$p.value  <-  NA
}

## }}}

## Ratio of AUC
data$pmmarkAUC__ <- time-exit
## linear model for stability when fitting
RAUCl<- recregIPCW(formR,data,cause=cause,death.code=death.code,
	  times=time,cens.model=~strata(strata__),marks=data$pmmarkAUC__,
	  model="lin")
f <- function(p) p[-1]/p[1]
## reparametrize as baseline on log-scale, and log-ratio contrasts 
RAUCe <- estimate(RAUCl,function(p) c(log(p[1]),log((p[1]+p[-1])/p[1])))
p <- length(coef(RAUCl))
RAUCet <- estimate(RAUCe,as.list(2:p))
RAUClt <- estimate(RAUCl,as.list(2:p))

proptest <- estimate(proptest,null=0)

out <- list(pepe.mori=pepe.mori,RatioAUC=RAUCet,difAUC=RAUClt,
	    RAUCl=RAUCl,RAUCe=RAUCe,
	    score.test=score.test,score.iid=iidU0,prop.test=proptest,time=time)
class(out) <- "marginalTest"
return(out)
} ## }}} 

##' @export
print.marginalTest  <- function(x,...) {# {{{
  print(summary(x),...)
}# }}}

##' @export
summary.marginalTest <- function(object,...) {# {{{

p.values <- c(object$time,
        object$pepe.mori$compare$p.value,
        object$RatioAUC$compare$p.value,
        object$prop.test$compare$p.value,
        object$score.test$compare$p.value)
p.values <- matrix(p.values,ncol=1)
colnames(p.values) <- "p-value"
rownames(p.values) <- c("time","Pepe-Mori","Ratio-AUC","Proportionality","Proportionality-score-test")

res <- p.values
class(res) <- "summary.maginalTest"
return(res)
}# }}}

##' @export
print.summary.marginalTest <- function(x,...) { ## {{{
  cat("coeffients:\n")
  printCoefmat(x$coef,...)
  cat("\n")
} ## }}}


## proc_design <- mets:::proc_design
