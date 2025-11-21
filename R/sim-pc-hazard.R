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
#' library(mets)
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
#' @aliases simrchaz addCums lin.approx simCens
rchaz <- function(cumhazard,rr,n=NULL,entry=NULL,cum.hazard=TRUE,cause=1,extend=FALSE)
{# {{{
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
lin.approx <- function (x2, xfx, x = 1)
{ ## {{{ 
    breaks <- xfx[, x]
    fx <- xfx[, -x]
    ri <- fast.approx(breaks, x2, type = "left")
    maxindex <- which(ri == length(breaks))
    rip1 <- ri + 1
    rip1[maxindex] <- length(breaks)
    rrr <- (x2 - breaks[ri])/(breaks[rip1] - breaks[ri])
    rrr[maxindex] <- 0
    res <- rrr * (fx[rip1] - fx[ri]) + fx[ri]
    res[is.na(res)] <- tail(fx, 1)
    res[is.na(rrr)] <- fx[ri][is.na(rrr)]
    return(res)
} ## }}}

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
simrchaz <- function(cumhazard,rr,n=NULL,cens=NULL,rrc=NULL,entry=NULL,...)
{# {{{
###   adds censoring to to rchaz call
   if (!is.null(n)) rr <- rep(1,n)
   n <- length(rr)  

   ptt <- rchaz(cumhazard,rr,entry=entry,...)

   if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
      dt <- data.frame(time=pmin(ptt$time,pct),status=ifelse(ptt$time<pct,ptt$status,0))
   } else dt <- ptt[,c("time","status")]

   return(dt)
}# }}}

#' @export
simCens <- function(cens,rrc=NULL,n=NULL,entry=NULL,...)
{# {{{
   if (is.null(rrc) & is.null(n)) stop("must give either rr or n\n"); 
   if (is.null(rrc)) rrc <- rep(1,n)
   n <- length(rrc)  
   if (is.matrix(cens)) {
	   pct <- rchaz(cens,rrc,entry=entry,...)
	   pct <- pct$time
   } else if (is.numeric(cens)) pct<- rexp(n)/(rrc*cens)  
   else stop("cens must be cumulative hazard or constant rate\n"); 

   if (!is.null(entry)) pct <- entry+pct
   return(pct)
} ## }}}

#' Simulation of Piecewise constant hazard models with two causes (Cox).
#' 
#' Simulates data from piecwise constant baseline hazard that can also be of
#' Cox type. Censor data at highest value of the break points for either of the
#' cumulatives, see also sim.phregs  
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
      pct <- simCens(cens,rrc=rrc,n=n,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

return(ptt)
}# }}}

#' @export 
rcrisks <-function(cumhazs,rrs,n=NULL,cens=NULL,rrc=NULL,entry=NULL,causes=NULL,extend=NULL,...)
{#'# {{{
    if (!is.null(extend))  {
      if (is.numeric(extend)) 
      if (length(extend)!=length(cumhazs)) extend <- rep(extend[1],length(cumhaz))
      cumhaz <- extendCums(cumhaz,NULL,haza=extend)
   }

###  if (extend)  {
###    cumhazs <- extendCums(cumhazs,NULL)
### }

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
#' @param rr possible vector of relative risk for cox model.
#' @param strata possible vector of strata
#' @param entry delayed entry variable for simulation.
#' @param extend to extend possible stratified baselines to largest end-point 
#' given then takes average rate of in simulated data from cox model.
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param ... arguments for rchaz, for example entry-time.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' library(mets)
#' data(sTRACE)
#' nsim <- 100
#' coxs <-  phreg(Surv(time,status==9)~strata(chf)+vf+wmi,data=sTRACE)
#' sim3 <- sim.phreg(coxs,nsim,data=sTRACE)
#' cc <-   phreg(Surv(time,status)~strata(chf)+vf+wmi,data=sim3)
#' cbind(coxs$coef,cc$coef)
#' plot(coxs,col=1); plot(cc,add=TRUE,col=2)
#' 
#' @aliases draw.phreg sim.base simulate.cox sim.phregs setup.phreg
#' @export sim.phreg 
#' @usage sim.phreg(cox,n,data,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,...)
sim.phreg <- function(cox,n,data=NULL,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,...)
{# {{{
   scox1 <- draw.phreg(cox,n,data=data,...)
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

 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }


return(ptt)
}# }}}

#' @export draw.phreg
#' @usage draw.phreg(cox,n,data=NULL,Z=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
draw.phreg <- function(cox,n,data=NULL,Z=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
{# {{{
###if (!inherits(cox,"phreg")) stop("must be phreg object\n"); 

## give data so that design can be constructed based on model-formula
if (is.null(Z)) {
cid <- countID(data.frame(id=cox$id))
whereid <- which(cid$Countid==1)
if (drawZ==TRUE) xid <- sample(whereid,n,replace=TRUE) else xid <- id
## TMP vars <- all.vars(update(cox$formula,-1~.))
vars <- all.vars(cox$formula)
dataid <- data[xid,vars,drop=FALSE] 

desX <- readPhreg(cox,dataid)
Z <- desX$X
strata <- desX$strata
} else { xid <- 1:nrow(Z); n <- nrow(Z); dataid <- Z; strata <- rep(0,n)}

nz <- ncol(Z)
if (nz>0) rr <- exp(as.matrix(Z) %*% cox$coef) else rr <- rep(1,nrow(Z))
cumhaz <- rbind(c(0,0),cox$cumhaz)
   if (cox$nstrata>1) stratid <- strata else stratid <- NULL
   if (cox$nstrata>1) {
       stratname <-  substring(cox$strata.name,8,nchar(cox$strata.name)-1)
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
sim.phregs <- function(coxs,n,data=NULL,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,...)
{# {{{
   scox1 <- draw.phreg(coxs[[1]],n,data=data)
   datas <- scox1$data
   stratam <-  scox1$strata
   rrm <- scox1$rr 

   cumhazl <- list()
   cumhazl[[1]] <- basecumhaz(coxs[[1]],only=1)

   if (length(coxs)>1) 
   for (i in 2:length(coxs)) {
      coxn <- draw.phreg(coxs[[i]],n,data=data,drawZ=FALSE,id=scox1$id)
      ind <-  match(names(datas),names(coxn$data),nomatch=0)
      ind <- ind[ind!=0]
      if (length(ind)>1)  datas <- cbind(datas,coxn$data[,-ind,drop=FALSE]) else datas <- cbind(datas,coxn$data)
      rrm <- cbind(rrm,coxn$rr)
      stratam <- cbind(stratam,coxn$strata)
      cumhazl[[i]] <- basecumhaz(coxs[[i]],only=1)
   }
   datas$orig.id <- scox1$id
   if (is.null(rr)) rr <- rrm
   if (is.null(strata)) strata <- stratam
   lentry <- NULL

   if (!is.null(extend)) {
         flat <- unlist(cumhazl, recursive = FALSE)
         lengths <- lengths(cumhazl)

	restore <- function(flat, lengths) {
	  ends <- cumsum(lengths)
	  starts <- c(1, head(ends + 1, -1))
	  mapply(function(s, e) flat[s:e], starts, ends, SIMPLIFY = FALSE)
	}
	cumhazl <- extendCums(flat,NULL,haza=extend)
        cumhazl <- restore(cumhazl,lengths)
   }

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

 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

return(ptt)
}# }}}

#' @export
setup.phreg  <- function(cumhazard,coef,Znames=NULL,strata=NULL)
{# {{{
    cox <- list()
    cox$cumhaz <- cumhazard
    cox$coef <- coef
    if (is.null(strata)) { strata <- rep(0,nrow(cumhazard)); nstrata <- 1} else nstrata <- max(strata)+1
    cox$strata <- strata
    cox$nstrata <- nstrata
    class(cox) <- c("setup","phreg")
    attr(cox,"znames") <- Znames
    return(cox)
}# }}}

sim.phreg.base <- function(cox.baseline,n,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,...)
{# {{{
   if (is.null(strata))  strata <- rep(0,n) 
   if (is.null(rr)) rr <- rep(1,n)
   
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

 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
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
#' with replacement. This gives covariates like the data.  Calls sim.phregs
#' 
#' 
#' @param coxs list of cox models.
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed covariates).
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
#' library(mets)
#' nsim <- 100; data(bmt)
#' 
#' cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet+age,data=bmt)
#' cox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=bmt)
#' coxs <- list(cox1,cox2)
#' ## just calls sim.phregs !
#' dd <- sim.phregs(coxs,nsim,data=bmt,extend=c(0.001))
#' scox1 <- phreg(Surv(time,status==1)~strata(tcell)+platelet+age,data=dd)
#' scox2 <- phreg(Surv(time,status==2)~tcell+strata(platelet),data=dd)
#'
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
dt <- sim.phregs(coxs,n,data=data,rrc=rrc,cens=cens,)
return(dt)
}# }}}

#' @export
simsubdist <- function(cumhazard,rr,n=NULL,entry=NULL,type="cloglog",startcum=c(0,0),U=NULL,...)
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

  if (is.null(U)) U <- runif(n)

  if (type=="cloglog") {
      F1tau <- 1-exp(-tail(cumh,1)*rr)
      ttt <- -log(1-U*F1tau)/rr
  } else if (type=="logistic") {
     F1tau <- tail(cumh,1)*rr/(1+tail(cumh,1)*rr)
     v <- U*F1tau
     ttt <- exp(logit(v))/rr; 
  }  else if (type=="rr" | type=="cif") {
     F1tau <- tail(cumh,1)
     ttt <- U*F1tau
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
#' @param strata possible strata 
#' @param drawZ to random sample from Z or not 
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' given then takes average rate of in simulated data from cox model.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param cumstart to start cumulatives at time 0 in 0. 
#' @param U uniforms to use for drawing of timing for cumulative incidence. 
#' @param pU uniforms to use for drawing event type (F1,F2,1-F1-F2). 
#' @param type of model logistic,cloglog,rr 
#' @param ... arguments for simsubdist (for example Uniform variable for realizations)
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' library(mets)
#' data(bmt)
#' nsim <- 10000
#' 
#' ## logit cumulative incidence regression model 
#' cif <- cifreg(Event(time,cause)~platelet+age,data=bmt,cause=1)
#' estimate(cif)  
#' plot(cif,col=1)
#' simbmt <- sim.cif(cif,nsim,data=bmt)
#' dtable(simbmt,~status)
#' scif <- cifreg(Event(time,status)~platelet+age,data=simbmt,cause=1)
#' estimate(scif)
#' plot(scif,add=TRUE,col=2)
#' 
#' ## Fine-Gray cloglog cumulative incidence regression model 
#' cif <- cifregFG(Event(time,cause)~strata(tcell)+age,data=bmt,cause=1)
#' estimate(cif)  
#' plot(cif,col=1)
#' simbmt <- sim.cif(cif,nsim,data=bmt)
#' scif <- cifregFG(Event(time,status)~strata(tcell)+age,data=simbmt,cause=1)
#' estimate(scif)
#' plot(scif,add=TRUE,col=2)
#' 
#' ################################################################
#' #  simulating several causes with specific cumulatives 
#' ################################################################
#' cif1 <-  cifreg(Event(time,cause)~strata(tcell)+age,data=bmt,cause=1)
#' cif2 <-  cifreg(Event(time,cause)~strata(platelet)+tcell+age,data=bmt,cause=2)
#' cifss <-  list(cif1,cif2)
#' simbmt <- sim.cifs(list(cif1,cif2),nsim,data=bmt,extend=0.005)
#' dtable(simbmt,~status)
#' scif1 <-  cifreg(Event(time,status)~strata(tcell)+age,data=simbmt,cause=1)
#' scif2 <-  cifreg(Event(time,status)~strata(platelet)+tcell+age,data=simbmt,cause=2)
#' cbind(cif1$coef,scif1$coef)   
#' cbind(cif2$coef,scif2$coef)   
#'     
#' par(mfrow=c(1,2))   
#' ## Cause 1 follows the model    
#' plot(cif1); plot(scif1,add=TRUE,col=1:2,lwd=2)
#' ## Cause 2 : second cause is modified with restriction to satisfy F1+F2<= 1, so scaled down     
#' plot(cif2); plot(scif2,add=TRUE,col=1:2,lwd=2)
#'    
#' @aliases sim.cif sim.cifs sim.cif.base simul.cifs setup.cif subdist pre.cifs sim.cifsRestrict simsubdist invsubdist
#' @export sim.cif
sim.cif <- function(cif,n,data=NULL,Z=NULL,strata=NULL,drawZ=TRUE,cens=NULL,rrc=NULL,cumstart=c(0,0),U=NULL,pU=NULL,type=NULL,...)
{# {{{
## also extracts coefficients and baseline from coxph, cox.aalen, phreg
## and uses them assuming that the model is cloglog unless otherwise
if (!inherits(cif,"defined"))  {
des <- 	draw.phreg(cif,n,data=data,Z=Z,drawZ=drawZ,...)
if (is.null(strata)) strata <- des$strata
dats <- des$data[,-(1:2)]
Z <- des$Z; cumhaz <- des$cum; rr <- des$rr; 
id <- des$id
cumhaz <- rbind(cumstart,cumhaz)
znames <- names(dats); 
###model <- des$model
} else {
cumhaz <- cif$cumhaz
if (is.null(Z)) stop("must give Z")
rr <- exp(as.matrix(Z) %*% cif$coef)
if (is.data.frame(Z)) znames <- names(Z) else znames <- colnames(Z)
###model <- cif$model
id <- 1:nrow(Z)
if (is.null(strata)) strata <- rep(0,n)
}

if (!is.null(type)) model <- type else {
if (is.null(cif$propodds)) model <- "cloglog" else model <- "logistic" 
}
if (!is.list(cif$cumhaz)) cumhazl <- basecumhaz(cif,only=1) else cumhazl <- cif$cumhaz

if (!inherits(cif,"phreg")) {
if (model=="logistic2" | model=="logistic") ptt <- simsubdist(cumhaz,rr,type="logistic",U=U,...) else ptt <- simsubdist(cumhaz,rr,type="cloglog",U=U,...)
    ptt <- cbind(ptt,Z)
} else { ### phreg class# {{{
	ptt <- sim.cif.base(cumhazl,n,rr=rr,strata=strata,U=U,type=model)
        ptt <- cbind(ptt,dats)
 }# }}}

 ### adds censoring 
 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

   attr(ptt,"model") <- model
   attr(ptt,"id") <-  id
   attr(ptt,"znames") <- znames
   attr(ptt,"cumhaz") <-  cumhazl
   return(ptt)
}# }}}

#' @export sim.cifs
sim.cifs <- function(cifs,n,data=NULL,rr=NULL,strata=NULL,Z=NULL,cens=NULL,rrc=NULL,max.times=NULL,causes=c(1,2),U=NULL,pU=NULL,extend=TRUE,type=NULL,restrict=TRUE,...)
{# {{{
if (!is.list(cifs)) stop("Cif models in list form\n"); 
if (length(cifs)!=2) stop("Only two models\n"); 

if (!is.null(type)) { 
   if (length(type)==1) type <- rep(type,2)
   model1 <- type[1]; model2 <- type[2];
} else {
   if (is.null(cifs[[1]]$propodds)) model1 <- "cloglog" else model1 <- "logistic" 
   if (is.null(cifs[[2]]$propodds)) model2 <- "cloglog" else model2 <- "logistic" 
}
###if (!is.list(cifs[[1]]$cumhaz)) cumhazl <- basecumhaz(cif,only=1) else cumhazl <- cif$cumhaz

   scox1 <- draw.phreg(cifs[[1]],n,data=data,Z=Z)
   datas <- as.data.frame(scox1$data)
   stratam <-  scox1$strata
   rrm <- scox1$rr 

   cumhazl <- list()
   cumhazl[[1]] <- basecumhaz(cifs[[1]],only=1)

   if (length(cifs)>1) 
   for (i in 2:length(cifs)) {
      coxn <- draw.phreg(cifs[[i]],n,data=data,drawZ=FALSE,id=scox1$id,Z=Z)
      ind <-  match(names(datas),names(coxn$data),nomatch=0)
      ind <- ind[ind!=0]
      if (length(ind)>1)  datas <- cbind(datas,coxn$data[,-ind,drop=FALSE]) else datas <- cbind(datas,coxn$data)
      rrm <- cbind(rrm,coxn$rr)
      stratam <- cbind(stratam,coxn$strata)
      cumhazl[[i]] <- basecumhaz(cifs[[i]],only=1)
   }
   datas$orig.id <- scox1$id
   if (is.null(rr)) rr <- rrm
   if (is.null(strata)) strata <- stratam
   lentry <- NULL

   if (!is.null(extend)) {
         flat <- unlist(cumhazl, recursive = FALSE)
         lengths <- lengths(cumhazl)

	restore <- function(flat, lengths) {
	  ends <- cumsum(lengths)
	  starts <- c(1, head(ends + 1, -1))
	  mapply(function(s, e) flat[s:e], starts, ends, SIMPLIFY = FALSE)
	}
	cumhazl <- extendCums(flat,NULL,haza=extend)
        cumhazl <- restore(cumhazl,lengths)
   }

  cifs <- pre.cifs(cifs,max.times=max.times)
  cifs[[1]]$cumhaz <- cumhazl[[1]]
  cifs[[2]]$cumhaz <- cumhazl[[2]]

  tau <- tail(cumhazl[[1]][[1]],1)[1]
  ## simulate first  time
  sim1 <- sim.cif.base(cumhazl[[1]],n,rr=rr[,1],strata=strata[,1],U=U,type=model1)
  sim2 <- sim.cif.base(cumhazl[[2]],n,rr=rr[,2],strata=strata[,2],U=U,type=model2)

  ## drawing which cause  that is seen if any 
  F1tau <- sim1$F1tau
  if (restrict) F2tau <- sim2$F1tau * (1-F1tau) else F2tau <- sim2$F1tau
  ptot <- F1tau+F2tau
  if (!is.null(pU)) {
      rt <- 1*(pU< pmin(ptot,1))
      rb <- 1*(pU< sim1$F1tau) 
  } else {
      rt <- rbinom(n,1,pmin(ptot,1))
      rb <- rbinom(n,1,sim1$F1tau/ptot)
  }
  cause=ifelse(rb==1,1,2)
  time=ifelse(cause==causes[1],sim1$timecause,sim2$timecause)
  cause <- rt*cause
  time[cause==0] <- tau

  ptt <- data.frame(time=time,status=cause,cause=cause,ptot=ptot)
  ptt <- cbind(ptt,datas)

### adds censoring 
 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

   return(ptt)
}# }}}

#' @export 
sim.cif.base <- function(baseline,n,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,U=NULL,type=c("cloglog","logistic","rr")[1],...)
{# {{{
   if (is.null(strata))  strata <- rep(0,n) 
   if (is.null(rr)) rr <- rep(1,n)
   
   cumhaz <- baseline
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
      simj <- simsubdist(cumhazj,rr[whichi],type=type,U=U)
      simj$id <- ids[whichi]
      ptt  <-  rbind(ptt,simj)
    }
    dsort(ptt) <- ~id

 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }


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

##' @export
simul.cifs <- function(n,rho1,rho2,beta,rc=0.5,depcens=0,rcZ=0.5,bin=1,type=c("cloglog","logistic"),rate=1,Z=NULL,U=NULL,pU=NULL) {# {{{
    p=length(beta)/2
    tt <- seq(0,6,by=0.1)
    if (length(rate)==1) rate <- rep(rate,2)
    Lam1 <- rho1*(1-exp(-tt/rate[1]))
    Lam2 <- rho2*(1-exp(-tt/rate[2]))

    if (length(bin)==1) bin <- rep(bin,2)
    if (length(rcZ)==1) rcZ <- c(rcZ,0)

    if (is.null(Z)) 
    Z=cbind((bin[1]==1)*(2*rbinom(n,1,1/2)-1)+(bin[1]==0)*rnorm(n),(bin[2]==1)*(rbinom(n,1,1/2))+(bin[2]==0)*rnorm(n))
    p <- ncol(Z)
    colnames(Z) <- paste("Z",1:p,sep="")

    cif1 <- setup.cif(cbind(tt,Lam1),beta[1:p],Znames=colnames(Z),type=type[1])
    cif2 <- setup.cif(cbind(tt,Lam2),beta[(p+1):(2*p)],Znames=colnames(Z),type=type[1])

    data <- sim.cifs(list(cif1,cif2),n,Z=Z,U=U,pU=pU,type=type[1])

    if (!is.null(rc)) {
    if (depcens==0) censor=pmin(rexp(n,1)*(1/rc),6) else censor=pmin(rexp(n,1)*(1/(rc*exp(Z %*% rcZ))),6)
    } else censor <- 6 

    status=data$status*(data$time<=censor)
    time=pmin(data$time,censor)
    data <- data.frame(time=time,status=status)
    data <- cbind(data,Z)
    attr(data,"Lam1") <- cbind(tt,Lam1)
    attr(data,"Lam2") <- cbind(tt,Lam2)
    return(data)

}# }}}

#' @export
setup.cif  <- function(cumhazard,coef,Znames=NULL,type="logistic")
{# {{{
    cif <- list()
    cif$cumhaz <- cumhazard
    cif$se.cumhaz <- cumhazard
    cif$coef <- coef
    cif$model <- type
    cif$strata <- rep(0,nrow(cumhazard))
    cif$jumps <- 1:nrow(cumhazard)

    cif$nstrata <- 1
    class(cif) <- "defined"
    attr(cif,"znames") <- Znames
    return(cif)
}# }}}

##' @export
kumarsim <- function (n,rho1=0.71,rho2=0.40,rate = c(6.11,24.2),
	      beta=c(-0.67,0.59,-0.55,0.25,0.68,0.18,0.45,0.31),
	      labels= c("gp","dnr","preauto","ttt24(24,300]"),
	      depcens=0,type = c("logistic", "cloglog"),restrict=0 )
{# {{{
    p = length(beta)/2
    tt <- seq(0, 150, by = 1)
    if (length(rate) == 1)
    rate <- rep(rate, 2)
    Lam1 <- rho1 * (1 - exp(-tt/rate[1]))
    Lam2 <- rho2 * (1 - exp(-tt/rate[2]))
    
    ## fully saturated model for kumar covariates 
    Zdist <- c(0.21064815,0.02083333,0.05555556,0.01504630,
	       0.13888889,0.15393519, 0.02662037, 0.04398148, 
	       0.04745370, 0.02430556, 0.02199074, 0.03935185,
               0.02430556, 0.05555556, 0.01273148, 0.10879630)
    Zs <- expand.grid(gp=c(0,1),dnr=c(0,1),preauto=c(0,1),ttt24=c(0,1))
    Zs <- dsort(Zs,~ttt24+gp+dnr+preauto)

    samn <- sample(1:16,n,replace=TRUE,prob=c(Zdist))
    Z <- Zs[samn,]

    colnames(Z) <- labels
    cif1 <- setup.cif(cbind(tt, Lam1), beta[1:4], Znames = colnames(Z), type = type[1])
    cif2 <- setup.cif(cbind(tt, Lam2), beta[5:8], Znames = colnames(Z), type = type[1])
    if (restrict==0) 
    data <- sim.cifs(list(cif1, cif2), n, Z = Z)
    else {
    ## keep model 2 on logistic form
    data <- sim.cifs(list(cif2, cif1), n, Z = Z)
    data$status21 <- data$status
    data$status21[data$status==1] <- 2
    data$status21[data$status==2] <- 1
    data$status <- data$status21
    }

    ## kumar censoring, cox model 
    c0 <- list()     
    c0$cumhaz <- cbind(c(0,20,60,90,160),
		       c(0, 0.07797931, 0.28512764, 0.76116180, 1.95720759))
    c0$coef <- c(1.8503113,-0.6976226,0.5828763,-0.2000003)

    if (depcens == 0)
        censor  <- rchaz(c0$cumhaz,n=n)
    else {
        rrc <- exp(as.matrix(Z) %*% c0$coef)
        censor  <- rchaz(c0$cumhaz,rrc)
    }
    status = data$status * (data$time <= censor$time)
    time = pmin(data$time, censor$time)
    data <- data.frame(time = time, status = status)
    return(cbind(data, Z))
}# }}}

##' @export
kumarsimRCT <- function (n,rho1=0.71,rho2=0.40,rate = c(6.11,24.2),
	      beta=c(-0.67,0.59,-0.55,0.25,0.68,0.18,0.45,0.31),
	      labels= c("gp","dnr","preauto","ttt24(24,300]"),
	      nocens=0,addcens=1,rct=1,
	      type = c("logistic", "cloglog"),restrict=1,
	      censpar=c(1,1,1,1), F1par=c(1,1,1,1), F2par=c(1,1,1,1),
	      treatmodel=c(-0.18,-0.16,0.06,0.24) 
)
{# {{{
    p = length(beta)/2
    tt <- seq(0, 150, by = 1)
    if (length(rate) == 1)
    rate <- rep(rate, 2)
    Lam1 <- rho1 * (1 - exp(-tt/rate[1]))
    Lam2 <- rho2 * (1 - exp(-tt/rate[2]))
    
    ## fully saturated model for kumar covariates 
    Zdist <- c(0.21064815,0.02083333,0.05555556,0.01504630,
	       0.13888889,0.15393519, 0.02662037, 0.04398148, 
	       0.04745370, 0.02430556, 0.02199074, 0.03935185,
               0.02430556, 0.05555556, 0.01273148, 0.10879630)
    Zs <- expand.grid(gp=c(0,1),dnr=c(0,1),preauto=c(0,1),ttt24=c(0,1))
    Zs <- dsort(Zs,~ttt24+gp+dnr+preauto)

    samn <- sample(1:16,n,replace=TRUE,prob=c(Zdist))
    Z <- Zs[samn,]
    ## randomized gp instead
    if (rct==1) Z[,1] <- rbinom(n,1,0.5)
    ## randomized gp given other covariates 
    if (rct==2) {
	    lp <- as.matrix(cbind(1,Z[,-1])) %*% treatmodel 
	    Z[,1] <- rbinom(n,1, expit(lp))
    }
    colnames(Z) <- labels
    cif1 <- setup.cif(cbind(tt, Lam1), F1par*beta[1:4], Znames = colnames(Z), type = type[1])
    cif2 <- setup.cif(cbind(tt, Lam2), F2par*beta[5:8], Znames = colnames(Z), type = type[1])
    if (restrict==0) 
    data <- sim.cifs(list(cif1, cif2), n, Z = Z)
    else {
    ## keep model 2 on logistic form
    data <- sim.cifs(list(cif2, cif1), n, Z = Z)
    data$status21 <- data$status
    data$status21[data$status==1] <- 2
    data$status21[data$status==2] <- 1
    data$status <- data$status21
    }

    if (nocens==0) {
    ## kumar censoring, cox model 
    c0 <- list()     
    c0$cumhaz <- cbind(c(0,20,60,90,160),
		       c(0, 0.07797931, 0.28512764, 0.76116180, 1.95720759))
    c0$coef <-censpar* c(1.8503113,-0.6976226,0.5828763,-0.2000003)

	if (addcens)  {
	  ## draw from cox gp model  
         rrc <- exp(as.matrix(Z) %*% c(c(1,0,0,0)*c0$coef))
         censorgp  <- rchaz(c0$cumhaz,rrc)
	 ## draw from cox other covariates   
         rrc <- exp(as.matrix(Z) %*% c(c(0,1,1,1) * c0$coef))
         censoroc  <- rchaz(c0$cumhaz,rrc)
	 ## miniumum of the two censoring times 
	 censor <- censorgp
	 censor$time <- pmin(censoroc$time,censorgp$time)
	 censor$status <- pmax(censoroc$status,censorgp$status)
	} else {
        rrc <- exp(as.matrix(Z) %*% c0$coef)
        censor  <- rchaz(c0$cumhaz,rrc)
	}
    status = data$status * (data$time <= censor$time)
    time = pmin(data$time, censor$time)
    data <- data.frame(time = time, status = status)
    }
    return(cbind(data, Z))
}# }}}


