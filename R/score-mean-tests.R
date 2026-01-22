##' Pepe-Mori test, score-test type test, proportionality test 
##'
##' for competing risks data and recurrent events the pepe-mori test uses the weights suggested in Ghosh and Lin (2000)
##' @param formula with Event object
##' @param data data frame for computation
##' @param cause of interest (1 default)
##' @param cens.code censoring code 
##' @param ... Additional arguments to lower level funtions
##' @param death.code  code for death for marginal mean of recurrent events 
##' @param death.code.prop code for other causes of death for Fine-Gray regression model 
##' @param death.code codes for death (terminating event, 2 default)
##' @param time upper limit in pepe-mori and AUC integrals  otherwise max event time
##' @author Thomas Scheike
##' 
##' @examples
##' library(mets)
##' data(bmt,package="mets")
##' bmt$time <- bmt$time+runif(nrow(bmt))*0.01
##' bmt$id <- 1:nrow(bmt)
##' dcut(bmt) <- age.f~age
##' str(bmt)      
##'      
##' fg=cifregFG(Event(time,cause)~tcell,data=bmt,cause=1)
##'  
##' ## computing tests for difference  for CIF
##' pmt <- test_marginalMean(Event(time,cause)~strata(tcell)+cluster(id),data=bmt,cause=1,
##' 			 death.code=1:2,death.code.prop=2,cens.code=0)
##' pmt$pepe.mori
##' pmt$RatioAUC
##' pmt$prop.test
##' ## score test equialent to Gray's test but variance estimated differently 
##' pmt$score.test
##'  
##' ### age-groups  
##' pmt <- test_marginalMean(Event(time,cause)~strata(age.f)+cluster(id),data=bmt,cause=1,
##' 			 death.code=1:2,death.code.prop=2,cens.code=0)
##' pmt$pepe.mori
##' pmt$RatioAUC
##' pmt$prop.test
##' ## score test equialent to Gray's test but variance estimated differently 
##' pmt$score.test
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
##' pmt$pepe.mori
##' pmt$RatioAUC
##' pmt$prop.test
##' pmt$score.test
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

if (is.null(time)) time <- max(exit[status %in% cause])
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
 varU0 <- crossprod(iidU0)/n^2
 logrank.robust <- t(gradient) %*% solve(varU0) %*% gradient
 p.logrank.robust <- 1-pchisq(logrank.robust,nrow(gradient))

 score.test <- list(logrank.robust=logrank.robust,p.logrank.robust=p.logrank.robust,
		    test.statistic=gradient)

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
	kmss <- cbind(kms$time,t(kms$surv))

        ns <- c(table(strata[cid$reverse==1]))
	n1 <- ns[1]
	n2 <- ns[2]
	## using number at risk for weighting 
	kmss <- cbind(kms$time,t(kms$surv))
	wt <- (n1+n2)*kmss[,2]*kmss[,3]/(n1*kmss[,2]+n2*kmss[,3])
	Wt <- cumsum(diff(c(0,kms$time))*wt)
	Wtmark <- lin.approx(exit,rbind(0,cbind(kms$time,Wt)))
	Wfinal <- tail(Wt,1)
	###
	data$pmmark__ <- Wfinal-Wtmark

	pepe.mori <- recregIPCW(formR,data,cause=cause,death.code=death.code,
	    time=time,marks=data$pmmark__,cens.model=~strata(strata__),model="lin")
        pepe.mori <- estimate(pepe.mori,null=0)
} else {

  formRC <- as.formula(paste("Event(entry__,exit__,status__)~factor(cstrata__)+cluster(id__)"))
	contr <- contr.iid <- c()
        for (i in 1:(nlev-1)) {

   	   data$cstrata__  <-  1*(nstrata==i)
           formC <- as.formula(paste("Event(entry__,exit__,statusC__)~strata(cstrata__)"))

           kms <- km(formC,data)
           kmss <- cbind(kms$time,t(kms$surv))

           ns <- c(table(data$cstrata__[cid$reverse==1]))
	   n1 <- ns[1]
	   n2 <- ns[2]
	   ## using number at risk for weighting 
	   kmss <- cbind(kms$time,t(kms$surv))
	   wt <- (n1+n2)*kmss[,2]*kmss[,3]/(n1*kmss[,2]+n2*kmss[,3])
	   Wt <- cumsum(diff(c(0,kms$time))*wt)
	   Wtmark <- lin.approx(exit,rbind(0,cbind(kms$time,Wt)))
	   Wfinal <- tail(Wt,1)
	   ###
	   data$pmmark__ <- Wfinal-Wtmark

	   pmOut <- recregIPCW(formRC,data,cause=cause,death.code=death.code,
	       time=time,marks=data$pmmark__,cens.model=~strata(cstrata__),model="lin")

        contr <- c(contr,coef(pmOut)[2])
	contr.iid <- cbind(contr.iid,pmOut$iid[,2])
   }
   var.contr <- crossprod(contr.iid)
   pepe.mori <- estimate(coef=contr,vcov=var.contr,null=0)
}

##### computing S0 for the two strata 
###pt0S <- recreg(formula,data,cause=cause,death.code=death.code.prop,
###		   cens.model=~strata(strata__),beta=0,no.opt=TRUE)
###Yr <- vecAllStrata(pt0S$S0,pt0S$strata.jumps,pt0S$nstrata)
###nss <- pt0S$S0[headstrata(pt0S$strata.jumps,pt0S$nstrata)[,1]]
###risk02 <- Yr[,2]==0
###risk01 <- Yr[,1]==0
###Yr[risk01,1] <- ns[1]
###Yr[risk02,2] <- ns[2]
###head(Yr)
###wt <- (sum(ns)/(n1*n2))*Yr[,1]*Yr[,2]/(Yr[,1]+Yr[,2])
###Wt2 <- cumsum(diff(c(0,pt0S$cumhaz[,1]))*wt)
###Wt2 <- lin.approx(exit,rbind(0,cbind(pt0S$cumhaz[,1],Wt2)))
###Wfinal2 <- tail(Wt2,1)
###data$pmmark2__ <- Wfinal2-Wt2
###
###pmOut2 <- recregIPCW(formR,data,cause=cause,death.code=death.code,
###    time=time,marks=data$pmmark2__,cens.model=~strata(strata__),model="lin")

## }}}

## Ratio of AUC
data$pmmarkAUC__ <- time-exit
RAUC<- recregIPCW(formR,data,cause=cause,death.code=death.code,
		  time=time,cens.model=~strata(strata__),marks=data$pmmarkAUC__)
RAUC <- estimate(RAUC,null=0)

proptest <- estimate(proptest,null=0)

out <- list(pepe.mori=pepe.mori,RatioAUC=RAUC,score.test=score.test,prop.test=proptest,time=time)
return(out)
} ## }}} 

## proc_design <- mets:::proc_design
