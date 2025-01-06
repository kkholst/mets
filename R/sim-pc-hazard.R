#' Simulation of Piecewise constant hazard model (Cox).
#' 
#' Simulates data from piecwise constant baseline hazard that can also be of
#' Cox type. Censor data at highest value of the break points.
#' 
#' For a piecewise linear cumulative hazard the inverse is easy to compute with 
#' and delayed entry x we compute 
#' \deqn{\Lambda^{-1}(\Lambda(x) + E/RR)}, 
#' where RR are the relative risks and E is exponential with mean 1.
#' This quantity has survival function 
#' \deqn{P(T > t | T>x) = exp(-RR (\Lambda(t) - \Lambda(x)))}. 
#' 
#' @param cumhazard cumulative hazard, or piece-constant rates for periods defined by first column of input.
#' @param rr relative risk for simulations, alternatively when rr=1 specify n 
#' @param n number of simulation if rr not given 
#' @param entry delayed entry time for simuations.
#' @param cum.hazard specifies wheter input is cumulative hazard or rates.
#' @param cause name of cause 
#' @param extend  to extend piecewise constant with constant rate. Default is average rate over time from cumulative (when TRUE), if numeric then uses given rate.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' chaz <-  c(0,1,1.5,2,2.1)
#' breaks <- c(0,10,   20,  30,   40)
#' cumhaz <- cbind(breaks,chaz)
#' n <- 10
#' X <- rbinom(n,1,0.5)
#' beta <- 0.2
#' rrcox <- exp(X * beta)
#' 
#' pctime <- rchaz(cumhaz,n=10)
#' pctimecox <- rchaz(cumhaz,rrcox,entry=runif(n))
#' 
#' @export 
#' @aliases simrchaz addCums lin.approx
rchaz <- function(cumhazard,rr,n=NULL,entry=NULL,cum.hazard=TRUE,cause=1,extend=FALSE)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE
###  if (is.null(rr) & is.null(n)) stop("must give rr or n\n"); 
  if (!is.null(n)) rr <- rep(1,n)
  n <- length(rr)

  breaks <- cumhazard[,1]
  rates <- cumhazard[,2][-1]
  mm <- tail(breaks,1)
  if (cum.hazard==FALSE) {
        cumh <- cumsum(c(0,diff(breaks)*rates))
        cumhazard <- cbind(breaks,cumh)
  } else cumh <- cumhazard[,2] 
   ttt <- rexp(n)/rr
   if (cumhazard[1,2]>0)  { ## start cumulative hazard with a 0
###	   warning("Safest to start with cumulative hazard 0 to avoid problems\n"); 
	   cumhazard <- rbind(c(0,0),cumhazard)
	   cumh <- c(0,cumh)
   }
   ###
   if (!is.null(entry)) {
	   if (length(entry)==1) entry <- rep(entry,n) else entry <- entry
	   cumentry <- lin.approx(entry,cumhazard,x=1)
   } else { entry <- cumentry <- rep(0,n) }
   ###
   ttte <- ttt+cumentry
   rrx <- lin.approx(ttte,cumhazard,x=-1)
   rrx <- ifelse(rrx>mm,mm,rrx)
   status <- rep(0,n)
   status <- ifelse(rrx<mm,cause,status)
   tcum <- tail(cumhazard,1)
   extend.rate <- NULL
   if (is.logical(extend))
	   if (extend) extend.rate <- tcum[2]/tcum[1]
   if (is.numeric(extend)) { extend.rate <- extend; extend <- TRUE;}
   if (extend) {
	   whoe <- which(status==0)
	   nn <- length(whoe)
	   if (nn>0) {
	     textend <- rexp(nn)/(rr[whoe]*extend.rate)
	     rrx[whoe] <- mm+textend
	     status[whoe] <- 1
	   }
   }
   dt <- data.frame(entry=entry,time=rrx,status=status,rr=rr)
   colnames(dt) <- c("entry","time","status","rr"); 
   attr(dt,"cumhaz") <- cumhazard
   attr(dt,"extend.rate") <- extend.rate
   return(dt)
}# }}}

#' @export
lin.approx <- function(x2,xfx,x=1)
{# {{{
   ### x=1   gives  f(x2) 
   ### x=-1  gives  f^-1(x2) 
   breaks <- xfx[,x]
   fx     <- xfx[,-x]
   ri <- fast.approx(breaks,x2,type="left")
   maxindex <- which(ri==length(breaks))
   rip1 <- ri+1
   rip1[maxindex] <- length(breaks)
   rrr <- (x2-breaks[ri])/(breaks[rip1]-breaks[ri])
   rrr[maxindex] <- 0
   res <- rrr*(fx[rip1]-fx[ri])+fx[ri]
   res[is.na(res)] <- tail(fx,1)
   return(res)
}# }}}

#' @export
addCums <- function(cumB,cumA,max=NULL)
{# {{{
 ## take times
 times <- sort(unique(c(cumB[,1],cumA[,1])))
 if (!is.null(max)) times <- times[times<max]
 cumBjx <- lin.approx(times,cumB,x=1)
 cumAjx <- lin.approx(times,cumA,x=1)
 cumBA <- cumBjx+cumAjx
 return(cbind(times,cumBA))
}# }}}

#' @export
simrchaz <- function(cumhazard,rr,n=NULL,cens=NULL,rrc=NULL,...)
{# {{{
###   adds censoring to to rchaz call
   if (!is.null(n)) rr <- rep(1,n)
   n <- length(rr)  

   ptt <- rchaz(cumhazard,rr,...)

   if (is.null(rrc)) {
	   rrc <- rep(1,n)
   }
   if (is.matrix(cens)) {
	   pct <- rchaz(cens,rrc,...)
	   pct <- pct$time
   } else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
	   }
   }
   dt <- data.frame(time=pmin(ptt$time,pct), status=ifelse(ptt$time<pct,ptt$status,0))

   return(dt)
}# }}}

#' Simulation of Piecewise constant hazard models with two causes (Cox).
#' 
#' Simulates data from piecwise constant baseline hazard that can also be of
#' Cox type. Censor data at highest value of the break points for either of the
#' cumulatives.
#' 
#' @param cumhaz1 cumulative hazard of cause 1
#' @param cumhaz2 cumulative hazard of cause 1
#' @param rr1 number of simulations or vector of relative risk for simuations.
#' @param rr2 number of simulations or vector of relative risk for simuations.
#' @param n number of simulation if rr not given 
#' @param cens to censor further , rate or cumumlative hazard
#' @param rrc retlativ risk for censoring.
#' @param extend to extend the cumulative hazards to largest end-point 
#' @param ... arguments for rchaz 
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' data(bmt); 
#' cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
#' cox2 <- phreg(Surv(time,cause==2)~tcell+platelet,data=bmt)
#'
#' X1 <- bmt[,c("tcell","platelet")]
#' n <- 100
#' xid <- sample(1:nrow(X1),n,replace=TRUE)
#' Z1 <- X1[xid,]
#' Z2 <- X1[xid,]
#' rr1 <- exp(as.matrix(Z1) %*% cox1$coef)
#' rr2 <- exp(as.matrix(Z2) %*% cox2$coef)
#'
#' d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2)
#' dd <- cbind(d,Z1)
#'
#' scox1 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
#' scox2 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
#' par(mfrow=c(1,2))
#' plot(cox1); plot(scox1,add=TRUE)
#' plot(cox2); plot(scox2,add=TRUE)
#' cbind(cox1$coef,scox1$coef,cox2$coef,scox2$coef)
#' 
#' @export 
#' @aliases cause.pchazard.sim  rcrisks rchazl
rcrisk <-function(cumhaz1,cumhaz2,rr1,rr2,n=NULL,cens=NULL,rrc=NULL,extend=TRUE,...)
{#'# {{{
 
 if (extend)  {
    out <- extendCums(list(cumhaz1,cumhaz2),NULL)
    cumhaz1 <- out$cum1
    cumhaz2 <- out$cum2
 }

if (!is.null(n)) { rr1 <- rep(1,n); rr2 <- rep(1,n) }
n <- length(rr1); 
if (missing(rr2)) rr2 <-rep(1,n)

ptt1 <- rchaz(cumhaz1,rr1,cum.hazard=TRUE,...)
ptt2 <- rchaz(cumhaz2,rr2,cum.hazard=TRUE,...)
ptt <- data.frame(time=pmin(ptt1$time,ptt2$time),status=ifelse(ptt1$time<=ptt2$time,ptt1$status,ptt2$status*2),entry=ptt1$entry)

if (!is.null(cens)) {
	if (is.null(rrc)) rrc <- rep(1,n)
	if (is.matrix(cens)) {
		pct <- rchaz(cens,rrc,...)
		pct <- pct$time
	} else {
		if (is.numeric(cens)) pct<- rexp(n)/cens  else {
			chaz <-sum(ptt1$status+ptt2$status)/sum(ptt1$time+ptt2$time)  ## hazard averate T haz 
			pct<- rexp(n)/chaz 
		}
	}
	ptt <- data.frame(time=pmin(ptt$time,pct),status=ifelse(ptt$time<pct,ptt$status,0),entry=ptt$entry)
}

return(ptt)
}# }}}

#' @export 
rcrisks <-function(cumhazs,rrs,n=NULL,cens=NULL,rrc=NULL,entry=NULL,causes=NULL,extend=TRUE,...)
{#'# {{{
  if (extend)  {
    cumhazs <- extendCums(cumhazs,NULL)
 }

  status <- NULL
  if (!is.null(n)) rrs <- matrix(1,n,length(cumhazs)) 
  n <- nrow(rrs); 

  if (!is.null(cens)) cens <- list(cens)
  if (!is.null(cens)) { 
	  if (is.null(rrc)) rrc <- rep(1,n); 
	  rrs <- cbind(rrs,rrc)
  }
  cumss  <-  c(cumhazs,cens)
  for (i in seq_along(cumss))  cumss[[i]] <- rbind(0,cumss[[i]])
  cumss <- extendCums(cumss,NULL)

  ## first time
  ptt <- rchaz(cumss[[1]],rrs[,1],entry=entry)

  ## other times and always min 
  if (length(cumss)>=2)
  for (i in 2:length(cumss)) {
  cum <- cumss[[i]]
  ptt1 <- rchaz(cumss[[i]],rrs[,i],entry=entry)
  ptt <- data.frame(time=pmin(ptt$time,ptt1$time),entry=ptt1$entry,
                    status=ifelse(ptt$time<=ptt1$time,ptt$status,ptt1$status*i))
  }
  if (!is.null(cens)) ptt <- dtransform(ptt,status=0,status==length(cumss))

  if (!is.null(causes))  {
	  ptt$status <- c(0,causes)[ptt$status+1]
  }

return(ptt)
}# }}}

#' @export 
rchazl <- function(cumhaz,rr,...)
{# {{{
     l <- length(cumhaz)
     ## simulate first 
     tall <- rchaz(cumhaz[[1]],rr[,1],...)
     if (l>=2) 
     for (i in 2:l) {
	  tall2 <- rchaz(cumhaz[[i]],rr[,i],...)
          tall$status <- ifelse(tall$time<tall2$time,tall$status,i*tall2$status)
          tall$time <- pmin(tall$time,tall2$time)
	  }
 return(tall)
}
# }}}

#' @export 
cause.pchazard.sim<-function(cumhaz1,cumhaz2,rr1,rr2,cens=NULL,rrc=NULL,...)
{
	rcrisk(cumhaz1,cumhaz2,rr1,rr2,cens=NULL,rrc=NULL,...)
}

#' Simulation of output from Cox model.
#' 
#' Simulates data that looks like fit from Cox model. Censor data automatically
#' for highest value of the event times by using cumulative hazard. 
#' 
#' @param cox output form coxph or cox.aalen model fitting cox model.
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed
#' covariates).
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' given then takes average rate of in simulated data from cox model.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param entry delayed entry variable for simulation.
#' @param rr possible vector of relative risk for cox model.
#' @param Z possible covariates to use instead of sampling from data.
#' @param extend to extend possible stratified baselines to largest end-point 
#' @param ... arguments for rchaz, for example entry-time.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' data(sTRACE)
#' nsim <- 100
#' coxs <-  phreg(Surv(time,status==9)~strata(chf)+vf+wmi,data=sTRACE)
#' sim3 <- sim.phreg(coxs,nsim,data=sTRACE)
#' cc <-   phreg(Surv(time, status)~strata(chf)+vf+wmi,data=sim3)
#' cbind(coxs$coef,cc$coef)
#' plot(coxs,col=1); plot(cc,add=TRUE,col=2)
#' 
#' @aliases sim.phreg read.phreg read.fit sim.base simulate.cox sim.phregs
#' @export sim.cox
#' @usage sim.cox(cox,n,data=NULL,cens=NULL,rrc=NULL,entry=NULL,rr=NULL,Z=NULL,extend=TRUE,...)
sim.cox <- function(cox,n,data=NULL,cens=NULL,rrc=NULL,entry=NULL,rr=NULL,Z=NULL,extend=TRUE,...)
{# {{{
   des <-  read.fit(cox,n,data=data,...)
   Z <- des$Z
   cumhazard <- des$cum
   if (is.null(rr)) rr <- des$rr
   if (!is.null(des$strataname)) {
       stratacox <- Z[,des$strataname]
   }
   ids <- 1:n
   lentry <- NULL

if (!inherits(cox,"phreg")) {
    if (cumhazard[1,2]>0) cumhazard <- rbind(c(0,0),cumhazard)
    ptt <- rchaz(cumhazard,rr,entry=entry) 
    ptt <- cbind(ptt,Z)
} else {
   ptt <- data.frame()
   stratj <- cox$strata.jumps
   if (cox$nstrata>1) {

	cumhazs <- list()
 	for (j in 0:(cox$nstrata-1)) {
	    cumhazardj <- rbind(c(0,0),cox$cumhaz[stratj==j,])
	    cumhazs[[j+1]] <- cumhazardj
	}
        if (extend) cumhazs <- extendCums(cumhazs,NULL)

 	for (j in 0:(cox$nstrata-1)) {
###	cumhazardj <- rbind(c(0,0),cox$cumhaz[stratj==j,])
    	    cumhazardj <- cumhazs[[j+1]]
	if (!is.null(entry)) lentry <- entry[stratacox==j]
		pttj <- rchaz(cumhazardj,rr[stratacox==j],entry=lentry) 
		Zj <- Z[des$strataid==j,,drop=FALSE]
		pttj$id <- ids[stratacox==j]
		ptt  <-  rbind(ptt,pttj)
	}
	dsort(ptt) <- ~id
	drm(ptt) <- ~id
	} else {
		if (cumhazard[1,2]>0) cumhazard <- rbind(c(0,0),cumhazard)
		ptt <- rchaz(cumhazard,rr,entry=entry) 
	}
	ptt <- cbind(ptt,Z)
}

if (!is.null(cens))  {
if (is.null(rrc)) rrc <- rep(1,n)
if (is.matrix(cens)) {
	pct <- rchaz(cens,rrc,entry=entry)
	pct <- pct$time
}
else {
	if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	pct<- rexp(n)/chaz 
	if (!is.null(entry)) pct  <- entry + pct
    }
}
ptt$time <- pmin(ptt$time,pct)
ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
} 

	attr(ptt,"id") <- des$id
	return(ptt)
}# }}}

#' @export sim.phreg
#' @usage sim.phreg(cox,n,data,rr=NULL,entry=NULL,extend=NULL,cens=NULL,...)
sim.phreg <- function(cox,n,data=NULL,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,...)
{# {{{

   scox1 <- read.phreg(cox,n,data=data,...)
   dat <- scox1$data
   dat$orig.id <- scox1$id

   if (is.null(strata))  strata <- scox1$strata 
   if (is.null(rr)) rr <- scox1$rr  
   
   cumhaz <- cox$cum
   cumhaz <- basecumhaz(cox,only=1)
   if (!is.null(extend))  {
      if (is.numeric(extend)) 
      if (length(extend)!=length(cumhaz)) extend <- rep(extend[1],length(cumhaz))
      cumhaz <- extendCums(cumhaz,NULL,haza=extend)
   }
   ids <- 1:n
   lentry <- NULL

   ptt <- c()
   for (i in seq(length(cumhaz))) {
      whichi <- which(strata==i-1)
      cumhazj <- rbind(0,cumhaz[[i]])
      if (!is.null(entry)) lentry <- entry[whichi]
      simj <- rchaz(cumhazj,rr[whichi],entry=lentry) 
      simj$id <- ids[whichi]
      ptt  <-  rbind(ptt,simj)
    }
    dsort(ptt) <- ~id
    ptt <- cbind(ptt,dat)

if (!is.null(cens))  {
if (is.null(rrc)) rrc <- rep(1,n)
if (is.matrix(cens)) {
	pct <- rchaz(cens,rrc,entry=entry)
	pct <- pct$time
} else {
	if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	pct<- rexp(n)/chaz 
    }
    if (!is.null(entry)) pct  <- entry + pct
}
ptt$time <- pmin(ptt$time,pct)
ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
} 

return(ptt)
}# }}}

#' @export read.phreg
#' @usage read.phreg(cox,n,data=NULL,Z=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
read.phreg <- function(cox,n,data=NULL,Z=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
{# {{{
if (!inherits(cox,"phreg")) stop("must be phreg object\n"); 

## give data so that design can be constructed based on model-formula
if (is.null(Z)) {
cid <- countID(data.frame(id=cox$id))
whereid <- which(cid$Countid==1)
if (drawZ==TRUE) xid <- sample(whereid,n,replace=TRUE) else xid <- id
vars <- all.vars(update(cox$formula,-1~.))
dataid <- data[xid,vars,drop=FALSE] 
###    ms <- match(cox$strata.name,names(cox$model.frame))
###    stratname <-  substring(cox$strata.name,8,nchar(cox$strata.name)-1)
} else { xid <- 1:nrow(Z); n <- nrow(Z); dataid <- Z}

desX <- readPhreg(cox,dataid)
Z <- desX$X
strata <- desX$strata
nz <- ncol(Z)
if (nz>0) rr <- exp(as.matrix(Z) %*% coef(cox)) else rr <- rep(1,nrow(Z))
cumhaz <- rbind(c(0,0),cox$cumhaz)
###   if (drawZ==TRUE) {
###           xid <- sample(1:nrow(Z),n,replace=TRUE)
###   } else xid <- 1:nrow(Z)
###   if (drawZ==FALSE) xid <- 1:nrow(Z) 
   if (cox$nstrata>1) stratid <- strata else stratid <- NULL
###   rr <- rr[xid]
###   Z <- Z[xid,,drop=FALSE]
   if (cox$nstrata>1) {
       ms <- match(cox$strata.name,names(cox$model.frame))
       stratname <-  substring(cox$strata.name,8,nchar(cox$strata.name)-1)
###       Z[,stratname] <- stratid
   } else stratname <- NULL
   model <-c(class(cox),is.null(cox$propodds))

out <- list(Z=Z,cumhaz=cumhaz,rr=rr,id=xid,model=model,strata=strata,data=dataid, stratname=stratname)

return(out)
} ## }}}

#' @export sim.base
#' @usage sim.base(cumhaz,n,stratajump=NULL,cens=NULL,entry=NULL,strata=NULL,rr=NULL,rc=NULL,extend=TRUE,...)
sim.base <- function(cumhaz,n,stratajump=NULL,cens=NULL,entry=NULL,strata=NULL,rr=NULL,rc=NULL,extend=TRUE,...)
{# {{{

## stratajump that indentifies baselines of each strata
if (!is.null(stratajump)) {
us <- unique(stratajump)
nus <- length(us)
}
if (is.null(rr)) rr <- rep(1,n)
if (is.null(rc)) rc <- rep(1,n)
id <- 1:n
lentry <- NULL

if (is.null(stratajump))  {
  if (cumhaz[1,2]>0) cumhaz <- rbind(c(0,0),cumhaz)
  ptt <- rchaz(cumhaz,rr,entry=entry) 
} else {
ptt <- c()
cumhazs <- list()
us <- unique(stratajump)
ss  <- seq_along(us)
for (j in ss) {
    i <- us[j]
    jjs <- which(stratajump==i)
    cumhazardj <- cumhaz[jjs,]
    cumhazs[[j]] <- cumhazardj
}
if (extend) cumhazs <- extendCums(cumhazs,NULL)

## strata among n simulations
for (i in seq_along(unique(strata))) {
	cumhazardj <- cumhazs[[i]]
	j <- unique(strata)[i]
	js <- which(strata==j)
	idj <- id[js]
	if (!is.null(entry)) lentry <- entry[js]
	ns <- length(js)
	rrj <- rr[js]
	pttj <- cbind(rchaz(cumhazardj,rrj,entry=lentry),j)
	colnames(pttj)[5] <- "strata"
	pttj$id <- idj
	ptt  <-  rbind(ptt,pttj)
}
dsort(ptt) <- ~id
drm(ptt) <- ~id
} 

if (!is.null(cens))  {
   if (is.matrix(cens))  pct <- rchaz(cens,rc,entry=entry)$time  else  pct<- rexp(n)/(rc*cens) 
   ptt$time <- pmin(ptt$time,pct)
   ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
} 

return(ptt)
}# }}}

extendit <- function(cox,extend=NULL)
{# {{{
    cumhaz <- basecumhaz(cox,only=1)
    if (!is.null(extend))  {
       if (length(extend)!=length(cumhaz)) extend <- rep(extend[1],length(cumhaz))
       cumhaz <- extendCums(cumhaz,NULL,haza=extend)
   }
    return(cumhaz)
}# }}}

#' @export sim.phregs
sim.phregs <- function(coxs,n,data=NULL,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,...)
{# {{{

   scox1 <- read.phreg(coxs[[1]],n,data=data)
   datas <- scox1$data
   stratam <-  scox1$strata
   rrm <- scox1$rr 

   cumhazl <- list()
   cumhazl[[1]] <- extendit(coxs[[1]],extend=extend)

   if (length(coxs)>1) 
   for (i in 2:length(coxs)) {
      coxn <- read.phreg(coxs[[i]],n,data=data,drawZ=FALSE,id=scox1$id)
      ind <-  match(names(datas),names(coxn$data),nomatch=0)
      ind <- ind[ind!=0]
      if (length(ind)>1)  datas <- cbind(datas,coxn$data[,-ind,drop=FALSE]) else datas <- cbind(datas,coxn$data)
      rrm <- cbind(rrm,coxn$rr)
      stratam <- cbind(stratam,coxn$strata)
      cumhazl[[i]] <- extendit(coxs[[i]],extend=extend)
   }
   datas$orig.id <- scox1$id
   if (is.null(rr)) rr <- rrm
   if (is.null(strata)) strata <- stratam
   lentry <- NULL

   ## simulate first  time
   simdata <- sim.phreg.base(cumhazl[[1]],n,rr=rr[,1],strata=strata[,1])
   l <- length(coxs)
   if (l>=2) 
   for (i in 2:l) {
      tall2 <- sim.phreg.base(cumhazl[[i]],n,rr=rr[,i],strata=strata[,i])
      simdata$status <- ifelse(simdata$time<tall2$time,simdata$status,i*tall2$status)
      simdata$time <- pmin(simdata$time,tall2$time)
   }
   ptt <- cbind(simdata,datas)

if (!is.null(cens))  {
if (is.null(rrc)) rrc <- rep(1,n)
if (is.matrix(cens)) {
	pct <- rchaz(cens,rrc,entry=entry)
	pct <- pct$time
} else {
	if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	pct<- rexp(n)/chaz 
    }
    if (!is.null(entry)) pct  <- entry + pct
}
ptt$time <- pmin(ptt$time,pct)
ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
} 

return(ptt)
}# }}}

sim.phreg.base <- function(cox.baseline,n,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,...)
{# {{{
   if (is.null(strata))  strata <- scox1$strata 
   if (is.null(rr)) rr <- scox1$rr  
   
   cumhaz <- cox.baseline
   if (!is.null(extend))  {
      if (is.numeric(extend)) 
      if (length(extend)!=length(cumhaz)) extend <- rep(extend[1],length(cumhaz))
      cumhaz <- extendCums(cumhaz,NULL,haza=extend)
   }
   ids <- 1:n
   lentry <- NULL

   ptt <- c()
   for (i in seq(length(cumhaz))) {
      whichi <- which(strata==i-1)
      cumhazj <- rbind(0,cumhaz[[i]])
      if (!is.null(entry)) lentry <- entry[whichi]
      simj <- rchaz(cumhazj,rr[whichi],entry=lentry) 
      simj$id <- ids[whichi]
      ptt  <-  rbind(ptt,simj)
    }
    dsort(ptt) <- ~id

if (!is.null(cens))  {
if (is.null(rrc)) rrc <- rep(1,n)
if (is.matrix(cens)) {
	pct <- rchaz(cens,rrc,entry=entry)
	pct <- pct$time
} else {
	if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	pct<- rexp(n)/chaz 
    }
    if (!is.null(entry)) pct  <- entry + pct
}
ptt$time <- pmin(ptt$time,pct)
ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
} 

return(ptt)
}# }}}

#' Simulation of cause specific from Cox models.
#' 
#' Simulates data that looks like fit from cause specific Cox models. 
#' Censor data automatically. When censoring is given in the  list of causes this
#' will give censoring that looks like the data.  Covariates are drawn from data-set
#' with replacement. This gives covariates like the data. 
#' 
#' 
#' @param coxs list of cox models.
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed
#' covariates).
#' @param cens specifies censoring model, if NULL then only censoring for 
#'   	       each cause at end of last event of this type. 
#' 	       if "is.matrix" then uses cumulative. 
#'             hazard given, if "is.scalar" then uses rate for exponential, and if not
#'             given then takes average rate of in simulated data from cox model.
#'             But censoring can also be given as a cause.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param ... arguments for rchaz, for example entry-time
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' nsim <- 100; data(bmt)
#' 
#' cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=bmt)
#' cox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=bmt)
#' coxs <- list(cox1,cox2)
#' dd <- sim.cause.cox(coxs,nsim,data=bmt)
#' scox1 <- phreg(Surv(time,status==1)~strata(tcell)+platelet,data=dd)
#' scox2 <- phreg(Surv(time,status==2)~tcell+strata(platelet),data=dd)
#' cbind(cox1$coef,scox1$coef)
#' cbind(cox2$coef,scox2$coef)
#' par(mfrow=c(1,2))
#' plot(cox1); plot(scox1,add=TRUE); 
#' plot(cox2); plot(scox2,add=TRUE); 
#' 
#' @export sim.cause.cox
#' @usage sim.cause.cox(coxs,n,data=NULL,cens=NULL,rrc=NULL,...)
sim.cause.cox <- function(coxs,n,data=NULL,cens=NULL,rrc=NULL,...)
{# {{{
if (!is.list(coxs)) stop("Cox models in list form\n"); 

 ptt <- sim.cox(coxs[[1]],n,data=data)
 simcovs <- ptt[,(5:ncol(ptt))]
 dt <- ptt

 if (length(coxs)>=2)
 for (i in 2:length(coxs)) {
    ptt1 <- sim.cox(coxs[[i]],n,data=data,id=attr(ptt,"id"),Z=simcovs)
    dt <- data.frame(time=pmin(dt$time,ptt1$time),
    status=ifelse(dt$time<=ptt1$time,dt$status,ptt1$status*i))
 }
 dt <- cbind(dt,simcovs)

 if (!is.null(cens)) {
 if (is.matrix(cens)) {
     pct <- rchaz(cens,rrc)
     pct <- pct$time
 }
 else {
   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
      pct<- rexp(n)/chaz 
   }
}

dt <- cbind(data.frame(time=pmin(ptt$time,pct),status=ifelse(ptt$time<pct,ptt$status,0)),simcovs)
}

return(dt)
}# }}}

#' @export
simsubdist <- function(cumhazard,rr,n=NULL,entry=NULL,type="cloglog",startcum=c(0,0),...)
{# {{{
  ## Fine-Gray model cloglog F1= 1-exp(-cum(t)*rr)
  ## logistic                F1= cum(t)*rr/(1+cum(t)*rr)
  ## rr                      F1= cum(t)*rr,  rr=exp(X^t beta) 
  if (!is.null(n)) rr <- rep(1,n)
  entry=NULL

  logit <- function(p) log(p/(1-p))

  if (cumhazard[1,2]>0)  cumhazard <- rbind(startcum,cumhazard)
  breaks <- cumhazard[,1]
  cumh <- cumhazard[,2] 
  mm <- tail(breaks,1)
  n <- length(rr)

  if (type=="cloglog") {
      F1tau <- 1-exp(-tail(cumh,1)*rr)
      ttt <- -log(1-runif(n)*F1tau)/rr
  } else if (type=="logistic") {
     F1tau <- tail(cumh,1)*rr/(1+tail(cumh,1)*rr)
     v <- runif(n)*F1tau
     ttt <- exp(logit(v))/rr; 
  }  else if (type=="rr" | type=="cif") {
     F1tau <- tail(cumh,1)
     ttt <- runif(n)*F1tau
     ## rr only affects binomial draw 
     F1tau <- F1tau*rr 
  } else stop(" cloglog or logistic or give function (fun=) \n"); 
  ###
   entry <- cumentry <- rep(0,n)
   ttte <- ttt+cumentry
   rrx <- lin.approx(ttt,cumhazard,x=-1)
   timecause <- rrx
   ###
   rrx <- ifelse(rrx>mm,mm,rrx)
   status <- rbinom(n,1,F1tau) 
   rrx[status==0] <- mm
   dt <- data.frame(entry=entry,time=rrx,status=status,rr=rr,F1tau=F1tau,timecause=timecause)
   attr(dt,"cumhaz") <- cumhazard
   return(dt)
}# }}}

#' @export
invsubdist <- function(F1,u,entry=NULL,cond=1,ptrunc=NULL)
{# {{{
  ### computes inverse subdist F1^{-1}(u) 
  ### computes inverse subdist F1^{-1}(u)
  F1max <- tail(F1[,2],1)
  mmax <- tail(F1[,1],1)

  if (!is.null(entry))  F1entry <- subdist(F1,entry)[,2];  
  if (!is.null(entry)) if (is.null(ptrunc)) ptrunc <- 1- F1entry
  if (cond==0 & !is.null(entry))  u<-u*ptrunc+F1entry
  if (cond==1 & !is.null(entry))  u<-u*(F1max-F1entry)+F1entry 
  if (cond==1 & is.null(entry))   u<-u*F1max
  rrx  <- lin.approx(u,F1,x=-1)
  status <- rep(1,length(u))
  status[u>=F1max] <- 0
  dt <- data.frame(time=rrx,status=status,y=u)
  if (!is.null(entry)) dt <- cbind(dt,entry)
  attr(dt,"F1") <- F1
  attr(dt,"Fmax") <- F1max
  return(dt)
}# }}}

#' @export
subdist <- function(F1,times)
{# {{{
  rrx <-   lin.approx(times,F1,x=1)
  dt <- cbind(times,rrx)
  colnames(dt) <- c("time","subdist")
  return(dt)
}# }}}

#' Simulation of output from Cumulative incidence regression model 
#' 
#' Simulates data that looks like fit from fitted cumulative incidence model
#' 
#' @param cif output form prop.odds.subdist or ccr (cmprsk), can also call invsubdist with 
#'    with cumulative and linear predictor 
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed
#' covariates).
#' @param Z to use these covariates for simulation rather than drawing new ones. 
#' @param drawZ to random sample from Z or not 
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' given then takes average rate of in simulated data from cox model.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param cumstart to start cumulatives at time 0 in 0. 
#' @param ... arguments for invsubdist
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' data(bmt)
#' 
#' scif <-  cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,prop=NULL)
#' summary(scif)  
#' plot(scif)
#' ################################################################
#' #  simulating several causes with specific cumulatives 
#' ################################################################
#'
#' cif1 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=1,prop=NULL)
#' cif2 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=2,prop=NULL)
#' # dd <- sim.cifsRestrict(list(cif1,cif2),200,data=bmt)
#
#' dd <- sim.cifs(list(cif1,cif2),200,data=bmt)
#' scif1 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=1)
#' scif2 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=2)
#'    
#' par(mfrow=c(1,2))   
#' plot(cif1); plot(scif1,add=TRUE,col=2)
#' plot(cif2); plot(scif2,add=TRUE,col=2)
#' @export sim.cif
#' @aliases sim.cif sim.cifs subdist pre.cifs sim.cifsRestrict simsubdist invsubdist
#' @usage sim.cif(cif,n,data=NULL,Z=NULL,drawZ=TRUE,cens=NULL,rrc=NULL,cumstart=c(0,0),...)
sim.cif <- function(cif,n,data=NULL,Z=NULL,drawZ=TRUE,cens=NULL,rrc=NULL,cumstart=c(0,0),...)
{# {{{
## also extracts coefficients and baseline from coxph, cox.aalen, phreg
## and uses them assuming that the model is cloglog unless otherwise
if (!inherits(cif,"defined"))  {
des <- 	read.fit(cif,n,data=data,Z=Z,drawZ=drawZ,...)
Z <- des$Z; cumhaz <- des$cum; rr <- des$rr; 
id <- des$id
cumhaz <- rbind(cumstart,cumhaz)
znames <- names(Z); 
model <- des$model
} else {
cumhaz <- cif$cumhaz
if (is.null(Z)) stop("must give Z")
rr <- exp(as.matrix(Z) %*% cif$coef)
if (is.data.frame(Z)) znames <- names(Z) else znames <- colnames(Z)
model <- cif$model
id <- 1:nrow(Z)
}

if (!inherits(cif,"phreg")) {
if (model=="logistic2" | model=="logistic") ptt <- simsubdist(cumhaz,rr,type="logistic",...) else ptt <- simsubdist(cumhaz,rr,type="cloglog",...)
    ptt <- cbind(ptt,Z)
} else { ### phreg class# {{{
	ptt <- data.frame()
	if (cif$nstrata>1) {
		stratj <- cif$strata[cif$jumps]
		for (j in 0:(cif$nstrata-1)) {
                        cumhazardj <- rbind(c(0,0),cif$cumhaz[stratj==j,])
			if (model[3]) 
		          pttj <- simsubdist(cumhazardj,rr[des$strataid==j],type="cloglog") else 
	                  pttj <- simsubdist(cumhazardj,rr[des$strataid==j],type="logistic") 
			Zj <- Z[des$strataid==j,,drop=FALSE]
			pttj <- cbind(pttj,Zj)
			ptt  <-  rbind(ptt,pttj)
			ptt  <-  rbind(ptt,pttj)
		}
	} else {
	if (model[3]) ptt <- simsubdist(cumhaz,rr,type="cloglog") else 
	              ptt <- simsubdist(cumhaz,rr,type="logistic") 
		ptt <- cbind(ptt,Z)
	}
 }# }}}

  ### adds censoring 
   if (!is.null(cens))  {# {{{
      if (is.null(rrc)) rrc <- rep(1,n)
      if (is.matrix(cens)) {
	   pct <- rchaz(cens,rrc,...)
	   pct <- pct$time
      }
      else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
           }
      }
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
   } # }}}

   attr(ptt,"model") <- model
   attr(ptt,"id") <-  id
   attr(ptt,"znames") <- znames
   return(ptt)
}# }}}

#' @export sim.cifs
sim.cifs <- function(cifs,n,data=NULL,Z=NULL,cens=NULL,rrc=NULL,max.times=NULL,causes=c(1,2),...)
{# {{{

if (!is.list(cifs)) stop("Cif models in list form\n"); 

  ## must consider all models out to last observation times or max.times
  cifs <- pre.cifs(cifs,max.times=max.times)

  tau <- tail(cifs[[1]]$cum[,1],1)
  sim1 <- sim.cif(cifs[[1]],n,data=data,Z=Z)
  Z <- sim1[,attr(sim1,"znames")]
  if (!inherits(cifs[[2]],"defined"))  {
      sim2p <- read.fit(cifs[[2]],1,data=data,Z=NULL)
      Z2 <- data[attr(sim1,"id"),names(sim2p$Z)]
      sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2,drawZ=FALSE)
  } else { 
        Z2 <- Z[,attr(cifs[[2]],"znames")]
	sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2)
  }

  ptot <- sim1$F1tau+sim2$F1tau
  ###
  rt <- rbinom(n,1,pmin(ptot,1))
  rb <- rbinom(n,1,sim1$F1tau/ptot)
  cause=ifelse(rb==1,1,2)
  time=ifelse(cause==causes[1],sim1$timecause,sim2$timecause)
  cause <- rt*cause
  time[cause==0] <- tau

  ptt <- data.frame(time=time,status=cause,cause=cause,ptot=ptot)
  Zcovs <- cbind(Z,Z2)
  samecovs <-  match(names(Z2),names(Z)) 
  Ze <- Zcovs[,-samecovs]
  ptt <- cbind(ptt,Ze)

   if (!is.null(cens))  {# {{{
      if (is.null(rrc)) rrc <- rep(1,n)
      if (is.matrix(cens)) {
	   pct <- rchaz(cens,rrc,...)
	   pct <- pct$time
      }
      else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
           }
      }
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
   } # }}}

   return(ptt)
}# }}}

#' @export sim.cifsRestrict
sim.cifsRestrict <- function(cifs,n,data=NULL,Z=NULL,cens=NULL,rrc=NULL,max.times=NULL,causes=c(1,2),...)
{# {{{

if (!is.list(cifs)) stop("Cif models in list form\n"); 

  ## must consider all models out to last observation times or max.times
  cifs <- pre.cifs(cifs,max.times=max.times)

  tau <- tail(cifs[[1]]$cum[,1],1)
  sim1 <- sim.cif(cifs[[1]],n,data=data,Z=Z)
  Z <- sim1[,attr(sim1,"znames")]
  if (!inherits(cifs[[2]],"defined"))  {
      sim2p <- read.fit(cifs[[2]],1,data=data,Z=NULL)
      Z2 <- data[attr(sim1,"id"),names(sim2p$Z)]
      sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2,drawZ=FALSE)
  } else { 
        Z2 <- Z[,attr(cifs[[2]],"znames")]
	sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2)
  }


  ptot <- sim1$F1tau+sim2$F1tau*(1-sim1$F1tau)
  ###
  rt <- rbinom(n,1,pmin(ptot,1))
  rb <- rbinom(n,1,sim1$F1tau/ptot)
  cause=ifelse(rb==1,1,2)
  time=ifelse(cause==causes[1],sim1$timecause,sim2$timecause)
  cause <- rt*cause
  time[cause==0] <- tau

  ptt <- data.frame(time=time,status=cause,cause=cause,ptot=ptot)
  ptt <- cbind(ptt,Z)
  samecovs <-  match(names(Z2),names(Z)) 
  Ze <- Z2[,-samecovs]
  ptt <- cbind(ptt,Ze)

   if (!is.null(cens))  {# {{{
      if (is.null(rrc)) rrc <- rep(1,n)
      if (is.matrix(cens)) {
           cum.hazard <- TRUE
	   pct <- rchaz(cens,rrc,...)
	   pct <- pct$time
      }
      else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
           }
      }
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
   } # }}}

   return(ptt)
}# }}}

#' @export pre.cifs
pre.cifs <- function(cifs,n,data=NULL,pprint=FALSE,max.times=NULL,...)
{# {{{

if (!is.list(cifs)) stop("Cif models in list form\n"); 

## iff needed put cumulative baseline in argument cifs$cum
## extracts maximum time of all cumulatives
maxtimes <- rep(0,length(cifs))
  for (i in 1:length(cifs)) { 
     if (inherits(cifs[[i]],"coxph")) {
	     stop("Use phreg of mets instead\n"); 
     } else  if (inherits(cifs[[i]],"crr")) {
	cifs[[i]]$cum <- rbind(c(0,0),cbind(cifs[[i]]$uftime,cumsum(cifs[[i]]$bfitj)))
        maxtimes[i] <- max(cifs[[i]]$uftime)
     } else if (inherits(cifs[[i]],"phreg")) {
	cifs[[i]]$cumhaz <- cifs[[i]]$cumhaz
        maxtimes[i] <- max(cifs[[i]]$cumhaz[,1])
     }  else if (inherits(cifs[[i]],"defined")) {
	cifs[[i]]$cumhaz <- cifs[[i]]$cumhaz
        maxtimes[i] <- max(cifs[[i]]$cumhaz[,1])
     } else {
        maxtimes[i] <- max(cifs[[i]]$cum[,1])
     }
  }

  mmtimes <- min(maxtimes)
  if (is.null(max.times)) mtimes <- mmtimes else mtimes <- max.times
  if (mtimes > mmtimes)  {
	  warning("max.time is further out than max for some cumulatives\n"); 
	  cat(" Max times for cumulatives in cifs \n"); 
	  print(maxtimes) 
  }

  for (i in 1:length(cifs)) { 
	cums <- cifs[[i]]$cum
        keep <- cums[,1]<mtimes
        Fmm <- subdist(cums,mtimes)
	cums <- cums[keep,,drop=FALSE]
	cums <- rbind(cums,Fmm)
	cifs[[i]]$cum <- cums
	cifs[[i]]$cumhaz <- cums
	if (pprint) {
	print(head(cums))
	print(tail(cums))
	}
  }

  return(cifs)
}# }}}

## reads a coxph, cox.aalen, phreg, crr, comprisk, prop.odds.subdist object
## and returns cumulative hazard, linear predictor of a resample of size
## n of data, for coxph cox.aalen, comprisk the design Z is constructed
## using the data, but Z can also be given which is needed for crr
## if drawZ is true the covariates in Z are used but multplied coefficient 
## based on column names, if fixZ=TRUE the matrix Z is multiplied coefficients 
## as is Z %*% coef(x)  to form linear predictor
read.fit <- function(cox,n,data=NULL,Z=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
{# {{{

if (inherits(cox,"coxph"))
{# {{{
   ## basehaz(cox,centered=FALSE)[,c(2,1)] 
   sfit <- survival::survfit(cox, se.fit=FALSE)   
   zcoef <- ifelse(is.na(coef(cox)), 0, coef(cox))
   offset <- sum(cox$means * zcoef)
   chaz <- sfit$cumhaz * exp(-offset)
   base <- cbind(sfit$time,chaz)  
   mt <- model.frame(cox)
   Y <- model.extract(mt, "response")
   if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
   if (attr(Y, "type") == "right") {
    time <- Y[, "time"]; 
    status <- Y[, "status"]
   } else stop("Expected right-censored data.");
  if (is.null(Z)) Z <- na.omit(model.matrix(cox))
  nn <- colnames(Z)
  if (fixZ) Z <- Z else Z <- Z[,names(coef(cox)),drop=FALSE] 
  rownames(Z) <- NULL
  jtime <- sort(time[status==1])
  cumhazard <- rbind(c(0,0),cpred(base,jtime))
  rr <- exp( Z %*% matrix(coef(cox),ncol=1))
  if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
  if (!is.null(id)) xid <- id
  rr <- rr[xid]
  Z <- Z[xid,,drop=FALSE]
  colnames(Z) <- nn
  model <- "fg"
}# }}}
if (inherits(cox,"cox.aalen"))
{# {{{
   formula <- attr(cox, "Formula")
###   Z1 <- na.omit(model.matrix(cox,data=data))
   p <- nrow(cox$gamma)
   if (is.null(Z)) {
   Z <- na.omit(get_all_vars(formula,data=data))
   nz <- ncol(Z)
   Z <- Z[,seq(nz-p+1,nz)]
   }
   lrr <- as.matrix(Z) %*% cox$gamma
   cumhazard <- cox$cum
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   model <- "fg"
   if (cox$prop.odds==TRUE) model <- "logistic"
}# }}}
if (inherits(cox,"phreg"))
{# {{{
   p <- length(cox$coef)
   if (is.null(Z)) {
      Z <- cox$model.frame[,-1,drop=FALSE]
      if (cox$nstrata>1) {
	   ms <- match(cox$strata.name,names(Z))
           stratname <-  substring(cox$strata.name,8,nchar(cox$strata.name)-1)
	   strata <- cox$strata
	   Z  <-  Z[,-ms,drop=FALSE]
      }
   }
   nz <- ncol(Z)
   if (fixZ) Z <- Z else Z <- Z[,names(cox$coef),drop=FALSE] 
   lrr <- as.matrix(Z) %*% cox$coef
   cumhazard <- rbind(c(0,0),cox$cumhaz)
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   if (cox$nstrata>1) stratid <- strata[xid] else stratid <- NULL
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   if (cox$nstrata>1) {
	   Z[,stratname] <- stratid
   }
   model <-c(class(cox),is.null(cox$propodds))
   ## if (cox$prop.odds==TRUE) cox$model <- "logistic"
}# }}}
if (inherits(cox,"crr"))
{# {{{
   p <- length(cox$coef)
   if (is.null(Z)) stop("must give covariates for crr simulations\n");
   nz <- ncol(Z)
   rownames(Z) <- NULL
   if (fixZ) Z <- Z else Z <- Z[,names(cox$coef),drop=FALSE] 
   lrr <- as.matrix(Z) %*% cox$coef
   cumhazard <- rbind(c(0,0),cbind(cox$uftime,cumsum(cox$bfitj)))
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   model <- "fg"
}# }}}
if (inherits(cox,"comprisk"))
{# {{{
   p <- length(cox$gamma)
   formula <- attr(cox, "Formula")
###   Z1 <- na.omit(model.matrix(cox,data=data))
   p <- nrow(cox$gamma)
   if (is.null(Z)) {
	   Z <- na.omit(get_all_vars(formula,data=data))
	   nz <- ncol(Z)
	   Z <- Z[,seq(nz-p+1,nz),drop=FALSE]
   }
   nz <- ncol(Z)
###   if (fixZ) Z <- Z else Z <- Z[,names(cox$gamma),drop=FALSE] 
   lrr <- as.matrix(Z) %*% cox$gamma
###   lrr <- as.matrix(Z) %*% cox$gamma
   cumhazard <- cox$cum
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   if (!(cox$model %in% c("fg","logistic2"))) stop("musg be fg or logistic2"); 
   model <- cox$model
}# }}}

out <- list(Z=Z,cum=cumhazard,rr=rr,id=xid,model=model)
if (inherits(cox,"phreg"))
   if (cox$nstrata>1) out <- c(out,list(strataid=stratid,strataname=stratname))
return(out)

}# }}}

