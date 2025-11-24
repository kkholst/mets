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
#' @aliases simrchaz addCums lin.approx simCens extendCums 
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
rchazl <- function (cumhaz, rr, causes=NULL,...)
{# {{{
    l <- length(cumhaz)
    tall <- rchaz(cumhaz[[1]], rr[, 1], ...)
    if (l >= 2)
        for (i in 2:l) {
            tall2 <- rchaz(cumhaz[[i]], rr[, i], ...)
            tall$status <- ifelse(tall$time < tall2$time, tall$status, i * tall2$status)
            tall$time <- pmin(tall$time, tall2$time)
        }
 if (!is.null(causes)) {
      where <- which(tall$status!=0) 
      tall$status[where] <- causes[tall$status[where]]
 }
    return(tall)
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

##' @export
extendCums <- function(cumA,cumB,extend=NULL)
{# {{{
## setup as list to run within loop
if (!is.null(cumB)) {cumA <- list(cumA,cumB); } 

###print(str(cumA))
###print("____________________")
###print(lapply(cumA, function(x) is.data.frame(x) | is.matrix(x)))

## to also work with strata version where each list contains a list of cumHaz for strata
matrixlist <- any(unlist(lapply(cumA, function(x) is.data.frame(x) | is.matrix(x))))
nn <- length(cumA)
nl <- lengths(cumA)
if (!matrixlist) cumA <- unlist(cumA, recursive = FALSE)

restore <- function(flat, lengths) {
  ends <- cumsum(lengths)
  starts <- c(1, head(ends + 1, -1))
  mapply(function(s, e) flat[s:e], starts, ends, SIMPLIFY = FALSE)
}

 maxx <- unlist(lapply(cumA,function(x) tail(x,1)[1]))
 mm <- which.max(maxx)
 haza <- NULL
 if (is.numeric(extend)) 
 if (length(extend)!=length(cumA)) haza <- rep(extend,length(cumA)) 

 ## extend all that are not at maxtime
for (i in seq(nn)[-mm]) {
  cumB <- as.matrix(cumA[[i]]); 
  cumB <- rbind(c(0,0),cumB); 

  ### linear extrapolation of mortality using given dHaz/dt or haza given rate
  if (tail(cumB[,1],1)<maxx[mm]) {
      tailB <- tail(cumB,1)
      cumlast <- tailB[2]
      timelast <- tailB[1]
      if (is.null(haza)) hazb <- cumlast/timelast else hazb <- haza[i]
      cumB <- rbind(cumB,c(maxx[mm],cumlast+hazb*(maxx[mm]-timelast))) 
  }
  cumA[[i]] <- cumB
}
 cumA[[mm]] <- rbind(c(0,0),cumA[[mm]])

 if (!matrixlist) cumA <- restore(cumA,nl)
 return( setNames(cumA,paste("cum",seq(nn),sep="")))
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
   } else if (is.numeric(cens)) {
	   pct<- rexp(n)/(rrc*cens)  
           if (!is.null(entry)) pct <- entry+pct
   } else stop("cens must be cumulative hazard or constant rate\n"); 

   return(pct)
} ## }}}

#' Simulation of Piecewise constant hazard models with two causes (Cox).
#' 
#' Simulates data from piecwise constant baseline hazard that can also be of
#' Cox type. Censor data at highest value of the break points for either of the
#' cumulatives, see also sim.phregs  
#' 
#' @param cumA cumulative hazard of cause 1, or list of multiple cumulative hazards
#' @param cumB cumulative hazard of cause 2 or NULL when cumA is a list 
#' @param rr1 number of simulations or vector of relative risk for simuations, or matrix with columns equal to number of hazards in list
#' @param rr2 number of simulations or vector of relative risk for simuations.
#' @param n number of simulation if rr not given, must be given when rr is not given  
#' @param cens to censor further , rate or cumumlative hazard
#' @param rrc retlativ risk for censoring.
#' @param extend to extend the cumulative hazards to largest end-point 
#' @param causes to assign status values for each of the causes, vector of integers 
#' @param ... arguments for rchaz 
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' library(mets)
#' data(bmt); 
#' n <- 100
#' cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
#' cox2 <- phreg(Surv(time,cause==2)~tcell+platelet,data=bmt)
#'
#' X1 <- bmt[,c("tcell","platelet")]
#' xid <- sample(1:nrow(X1),n,replace=TRUE)
#' Z1 <- X1[xid,]
#' Z2 <- X1[xid,]
#' rr1 <- exp(as.matrix(Z1) %*% cox1$coef)
#' rr2 <- exp(as.matrix(Z2) %*% cox2$coef)
#'
#' d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2,cens=2/70)
#' dd <- cbind(d,Z1)
#'
#' d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2,cens=cbind(c(1,30,68),c(.01,1,3)))
#' dd <- cbind(d,Z1)
#'
#' scox0 <- phreg(Surv(time,status==0)~tcell+platelet,data=dd)
#' plot(scox0); lines(cbind(c(1,30,68),c(.01,1,3)))
#'
#' scox1 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
#' scox2 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
#' par(mfrow=c(1,2))
#' plot(cox1); plot(scox1,add=TRUE)
#' plot(cox2); plot(scox2,add=TRUE)
#' cbind(cox1$coef,scox1$coef,cox2$coef,scox2$coef)
#' 
#' @export 
#' @aliases  rchazl 
rcrisk <-function(cumA,cumB,rr1=NULL,rr2=NULL,n=NULL,cens=NULL,rrc=NULL,extend=TRUE,causes=NULL,...)
{#'# {{{
 if (!is.null(cumB)) {
	cumA <- list(cumA,cumB); 
        rr <- cbind(rr1,rr2);
 } else cumA <- c(cumA,cumB)

if (!is.null(n)) { 
   rr <- matrix(1,n,length(cumA)); 
} else {
n <- length(rr1); 
rr <- cbind(rr1,rr2)
}

if (!is.null(extend))  cumA <- extendCums(cumA,NULL,extend=extend)

 l <- length(cumA)
 ## simulate first 
 ptt <- rchaz(cumA[[1]],rr[,1],...)
 if (l>=2) for (i in 2:l) {
 ptt2 <- rchaz(cumA[[i]],rr[,i],...)
 ptt$status <- ifelse(ptt$time<ptt2$time,ptt$status,i*ptt2$status)
 ptt$time <- pmin(ptt$time,ptt2$time)
 }
 if (!is.null(causes)) {
      where <- which(ptt$status!=0) 
      ptt$status[where] <- causes[ptt$status[where]]
 }

 ## add censoring 
if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
}

return(ptt)
}# }}}

#' Simulation of output from Cox model.
#' 
#' Simulates data that looks like fit from Cox model. Censor data automatically
#' for highest value of the event times by using cumulative hazard. 
#' 
#' @param cox output form coxph or cox.aalen model fitting cox model.
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed covariates).
#' @param Z give design matrix instead of data
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
#' set.seed(100)
#' sim3 <- sim.phreg(coxs,nsim,data=sTRACE)
#' head(sim3)
#' cc <-   phreg(Surv(time,status)~strata(chf)+vf+wmi,data=sim3)
#' cbind(coxs$coef,cc$coef)
#' plot(coxs,col=1); plot(cc,add=TRUE,col=2)
#' 
#' Z <- sim3[,c("vf","chf","wmi")]
#' strata <- sim3[,c("chf")]
#' rr <- exp(as.matrix(Z[,-2]) %*% coef(coxs))
#' sim4 <- sim.phreg(coxs,nsim,data=NULL,rr=rr,strata=strata)
#' sim4 <- cbind(sim4,Z)
#' cc <-   phreg(Surv(time,status)~strata(chf)+vf+wmi,data=sim4)
#' cbind(coxs$coef,cc$coef)
#' plot(coxs,col=1); plot(cc,add=TRUE,col=2)
#' 
#' @aliases draw.phreg setup.phreg
#' @export sim.phreg 
#' @usage sim.phreg(cox,n,data=NULL,Z=NULL,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,...)
sim.phreg <- function(cox,n,data=NULL,Z=NULL,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,...)
{# {{{

if  (!is.null(data)) {
   scox1 <- draw.phreg(cox,n,data=data,...)
   dat <- scox1$data
   dat$orig.id <- scox1$id

   if (is.null(strata))  strata <- scox1$strata 
   if (is.null(rr)) rr <- scox1$rr  
   n <- length(rr)
   id <- 1:length(rr)
} else {
	if (!is.null(Z)) {
	  if (is.data.frame(Z)) znames <- names(Z) else znames <- colnames(Z)
	} else znames <- NULL
	if (is.null(rr) & (!is.null(Z))) rr <- exp(as.matrix(Z) %*% cif$coef)
	if (is.null(rr) & is.null(Z)) rr <- rep(1,n)
	id <- 1:length(rr)
	n <- length(rr)
	dat <- NULL
}
if (is.null(strata)) strata <- rep(0,n)
   
   if (inherits(cox,c("phreg","cifreg"))) cumhaz <- basecumhaz(cox,only=1)
   else {
        if (!is.list(cox)) stop("must be phreg or list of hazards\n") else cumhaz <- cox
   }
   if (!is.null(extend))  cumhaz <- extendCums(cumhaz,NULL,extend=extend)
   ids <- 1:n
   lentry <- NULL

   ptt <- c()
   for (i in unique(strata)) {
      whichi <- which(strata==i)
      cumhazj <- rbind(0,cumhaz[[i+1]])
      if (!is.null(entry)) lentry <- entry[whichi]
      simj <- rchaz(cumhazj,rr[whichi],entry=lentry) 
      simj$id <- ids[whichi]
      ptt  <-  rbind(ptt,simj)
    }
    dsort(ptt) <- ~id
    if (!is.null(dat)) ptt <- cbind(ptt,dat)

 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

return(ptt)
}# }}}

#' @export draw.phreg
#' @usage draw.phreg(cox,n,data=NULL,Z=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
draw.phreg <- function(cox,n,data=NULL,Z=NULL,strata=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
{# {{{
###if (!inherits(cox,"phreg")) stop("must be phreg object\n"); 

## give data so that design can be constructed based on model-formula
if (!is.null(data)) {
cid <- countID(data.frame(id=cox$id))
whereid <- which(cid$Countid==1)
if (drawZ==TRUE) xid <- sample(whereid,n,replace=TRUE) else xid <- id
## TMP vars <- all.vars(update(cox$formula,-1~.))
vars <- all.vars(cox$formula)
dataid <- data[xid,vars,drop=FALSE] 

desX <- readPhreg(cox,dataid)
Z <- desX$X
strata <- desX$strata
} else {  ## Z and strata
   xid <- 1:nrow(Z); 
   n <- nrow(Z); 
   dataid <- Z; 
   if (is.null(strata)) strata <- rep(0,n)
}

nz <- ncol(Z)
if (nz>0) rr <- exp(as.matrix(Z) %*% cox$coef) else rr <- rep(1,nrow(Z))
cumhaz <- rbind(c(0,0),cox$cumhaz)
   if (cox$nstrata>1) {
      stratname <-  substring(cox$strata.name,8,nchar(cox$strata.name)-1)
   } else stratname <- NULL
   model <-c(class(cox),is.null(cox$propodds))

out <- list(Z=Z,cumhaz=cumhaz,rr=rr,id=xid,model=model,
	    strata=strata,data=dataid,stratname=stratname)

return(out)
} ## }}}

#' @export
setup.phreg  <- function(cumhazard,coef,Znames=NULL,strata=NULL)
{# {{{
    cox <- list()
    cox$cumhaz <- cumhazard
    cox$coef <- coef
    if (is.null(strata)) { strata <- rep(0,nrow(cumhazard)); nstrata <- 1} else nstrata <- max(strata)+1
    cox$strata <- strata
    cox$nstrata <- nstrata
    cox$strata.name <- ""
    cox$jumps <- 1:nrow(cumhazard)
    class(cox) <- c("setup","phreg")
    attr(cox,"znames") <- Znames
    return(cox)
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
#' @param rr possible vector of relative risk for cox model.
#' @param strata possible vector of strata 
#' @param entry delayed entry variable for simulation.
#' @param extend to extend possible stratified baselines to largest end-point 
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
#' data(bmt)
#' nsim <- 100; 
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
#' @export sim.phregs
sim.phregs <- function(coxs,n,data=NULL,rr=NULL,strata=NULL,entry=NULL,extend=NULL,cens=NULL,rrc=NULL,...)
{# {{{
   scox1 <- draw.phreg(coxs[[1]],n,data=data)
   datas <- scox1$data
   stratam <-  scox1$strata
   rrm <- scox1$rr 

   cumhazl <- list()
   cumhazl[[1]] <- basecumhaz(coxs[[1]],only=1)
i=2
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

   if (!is.null(extend)) cumhazl <- extendCums(cumhazl,NULL,extend=extend)

   ## simulate first  time
   simdata <- sim.phreg(cumhazl[[1]],n,rr=rr[,1],strata=strata[,1])
   l <- length(coxs)
   if (l>=2) 
   for (i in 2:l) {
      tall2 <- sim.phreg(cumhazl[[i]],n,rr=rr[,i],strata=strata[,i])
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
#' @param rr possible vector of relative risk for cox model.
#' @param strata possible vector of strata 
#' @param drawZ to random sample from Z or not 
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' given then takes average rate of in simulated data from cox model.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param cumstart to start cumulatives at time 0 in 0. 
#' @param U uniforms to use for drawing of timing for cumulative incidence. 
#' @param pU uniforms to use for drawing event type (F1,F2,1-F1-F2). 
#' @param type of model logistic,cloglog,rr 
#' @param extend  to extend piecewise constant with constant rate. Default is average rate over time from cumulative (when TRUE), if numeric then uses given rate.
#' @param ... arguments for simsubdist (for example Uniform variable for realizations)
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' library(mets)
#' data(bmt)
#' nsim <- 100
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
#' ## can be off due to restriction F1+F2<= 1    
#' cbind(cif2$coef,scif2$coef)   
#'     
#' par(mfrow=c(1,2))   
#' ## Cause 1 follows the model    
#' plot(cif1); plot(scif1,add=TRUE,col=1:2,lwd=2)
#' # Cause 2:second cause is modified with restriction to satisfy F1+F2<= 1, so scaled down     
#' plot(cif2); plot(scif2,add=TRUE,col=1:2,lwd=2)
#'    
#' @aliases sim.cif sim.cifs simul.cifs setup.cif subdist simsubdist invsubdist
#' @export sim.cif
sim.cif <- function(cif,n,data=NULL,Z=NULL,rr=NULL,strata=NULL,drawZ=TRUE,cens=NULL,rrc=NULL,
		    cumstart=c(0,0),U=NULL,pU=NULL,type=NULL,extend=NULL,...)
{# {{{
## also extracts coefficients and baseline from cifreg
if (!is.null(data))  {
   des <- 	draw.phreg(cif,n,data=data,Z=Z,drawZ=drawZ,...)
   dats <- des$data
   dats$orig.id <- des$id

   if (is.null(strata))  strata <- des$strata 
   if (is.null(rr)) rr <- des$rr  
   znames <- names(dats); 
    n <- nrow(dats)
} else {
	if (is.null(Z) & is.null(rr)) stop("must give Z or rr ")
	if (!is.null(Z)) {
	  if (is.data.frame(Z)) znames <- names(Z) else znames <- colnames(Z)
	} else znames <- NULL
	if (is.null(rr)) rr <- exp(as.matrix(Z) %*% cif$coef)
	n <- length(rr)
	orig.id <- 1:n
	dats <- NULL
}
if (is.null(strata)) strata <- rep(0,n)

if (inherits(cif,c("phreg","cifreg"))) cumhaz <- basecumhaz(cif,only=1,extend=extend)
else {
if (!is.list(cif)) stop("must be cifreg or list of hazards\n") else cumhaz <- cif
}
if (!is.null(extend))   cumhaz <- extendCums(cumhaz,NULL,extend=extend)

if (is.null(type)) {
if (is.null(cif$propodds)) type <- "cloglog" else type <- "logistic" 
}
ids <- 1:n

## {{{ simulation of cif 
    ptt <- c()
    for (i in unique(strata)) {
      whichi <- which(strata==i)
      cumhazj <- rbind(0,cumhaz[[i+1]])
      if (!is.null(U)) Ui <- U[whichi] else Ui <- NULL
      simj <- simsubdist(cumhazj,rr[whichi],type=type,U=Ui)
      simj$id <- ids[whichi]
      ptt  <-  rbind(ptt,simj)
    }
    dsort(ptt) <- ~id
### }# }}}
 if (!is.null(dats)) ptt <- cbind(ptt,dats)

 ### adds censoring 
 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

   attr(ptt,"type") <- type
   attr(ptt,"znames") <- znames
   attr(ptt,"cumhaz") <-  cumhaz
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

   if (!is.null(extend)) cumhazl <- extendCums(cumhazl,NULL,extend=extend)

  cifs[[1]]$cumhaz <- cumhazl[[1]]
  cifs[[2]]$cumhaz <- cumhazl[[2]]

  tau <- tail(cumhazl[[1]][[1]],1)[1]
  ## simulate first  time
  sim1 <- sim.cif(cumhazl[[1]],n,rr=rr[,1],strata=strata[,1],U=U,type=model1)
  sim2 <- sim.cif(cumhazl[[2]],n,rr=rr[,2],strata=strata[,2],U=U,type=model2)

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

    data <- sim.cifs(list(cif1,cif2),n,Z=Z,U=U,pU=pU,type=type[1],extend=NULL)

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
    if (type=="logistic") cif$propodds <- 1 else cif$propodds <- NULL
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
	      depcens=0,type = c("logistic", "cloglog"),restrict=TRUE,...)
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
    data <- sim.cifs(list(cif2, cif1), n, Z = Z,extend=NULL,restrict=restrict,...)
    stat12 <- which(data$status %in% c(1,2))
    data$status[stat12] <- c(2,1)[data$status[stat12]]

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
	      type = c("logistic", "cloglog"),restrict=TRUE,
	      censpar=c(1,1,1,1), F1par=c(1,1,1,1), F2par=c(1,1,1,1),
	      treatmodel=c(-0.18,-0.16,0.06,0.24),... )
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
    data <- sim.cifs(list(cif2, cif1), n, Z = Z,extend=NULL,restrict=restrict,...)
    stat12 <- which(data$status %in% c(1,2))
    data$status[stat12] <- c(2,1)[data$status[stat12]]

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

######################## ##################### #####################
## illness death competing risks with two causes of death
## First simple model, then with covariates addedd 
## effect of transition into 2 being early for cause 3, and 4
######################## ##################### #####################
simMultistateII <- function(cumhaz,death.cumhaz,death.cumhaz2,n=NULL,
		    rr=NULL,rd=NULL,rd2=NULL,gamma23=0,gamma24=0,early2=10000,
		    gap.time=FALSE,max.recurrent=100,cens=NULL,rrc=NULL,extend=TRUE,...) 
{# {{{

  status <- NULL

  ## cumhaz is cumulative out of state 1 to state 2 
  if (!is.null(n)) rr <- matrix(1,n,1) 
  n <- nrow(rr); 
  ## covariate adjustment 
  if (is.null(rd))   rd  <- matrix(1,n,length(death.cumhaz))
  if (is.null(rd2))  rd2  <- matrix(1,n,length(death.cumhaz2))

  ####
  if (!is.null(cens)) cens <- list(cens)
  cumss  <-  c(list(cumhaz),death.cumhaz,death.cumhaz2,cens)
  ### extend of cumulatives
  for (i in seq_along(cumss))  cumss[[i]] <- rbind(0,cumss[[i]])
  ## extend range of all cumulatives to max range, by linear extension
  cumss <- extendCums(cumss,NULL,extend=extend)
  maxtime <- tail(cumss[[1]],1)[1]

 if (!is.null(cens)) { 
	  cens <- cumss[[1+length(death.cumhaz)+length(death.cumhaz2)+1]]
	  if (is.null(rrc)) rrc <- rep(1,n); 
	  ctime <- rchaz(cens,rrc)$time
  } else ctime <- rep(maxtime,n)

  ## hazards out of 1 and out of 2
  chaz1 <- cumss[1:3]; rr1 <- cbind(rr,rd)
  chaz2 <- cumss[4:5]; rr2 <- cbind(rd2)

  ## time out of state 1 
  tall <- rchazl(chaz1,rr1,causes=c(2,3,4))
  tall$id <- 1:n
  tall$status[tall$time>ctime] <- 0; 
  tall$time[tall$time>ctime] <- ctime[tall$time>ctime] 
  tall$from <- 1
  tall$to <- tall$status

  ## simulating out of 2 
  tall2 <- subset(tall,status==2)
  rr23t <- exp(gamma23*(tall2$time<early2))
  rr24t <- exp(gamma24*(tall2$time<early2))
  sim2 <- rchazl(chaz2,rr2[tall2$id,]*cbind(rr23t,rr24t),causes=c(3,4),entry=tall2$time)
  sim2$id <- tall2$id
  ctime2 <- ctime[tall2$id]
  sim2$status[sim2$time>ctime2] <- 0; 
  sim2$time[sim2$time>ctime2] <- ctime2[sim2$time>ctime2] 
  sim2$from <- 2
  sim2$to <- sim2$status
  tall <- rbind(tall,sim2)

  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  return(tall)
  }# }}}

##' Simulation of illness-death model 
##'
##' Simulation of illness-death model 
##'
##' simMultistate with different death intensities from states 1 and 2 
##'
##' Must give cumulative hazards on some time-range 
##'
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of going from state 1 to 2.
##' @param cumhaz2  cumulative hazard of going from state 2 to 1. 
##' @param death.cumhaz cumulative hazard of death from state 1. 
##' @param death.cumhaz2 cumulative hazard of death from state 2.
##' @param rr  relative risk adjustment for cumhaz
##' @param rr2  relative risk adjustment for cumhaz2
##' @param rd  relative risk adjustment for death.cumhaz
##' @param rd2  relative risk adjustment for death.cumhaz2
##' @param rrc  relative risk adjustment for censoring 
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param dependence 0:independence; 1:all share same random effect with variance var.z; 2:random effect exp(normal) with correlation structure from cor.mat; 3:additive gamma distributed random effects, z1= (z11+ z12)/2 such that mean is 1 , z2= (z11^cor.mat(1,2)+ z13)/2, z3= (z12^(cor.mat(2,3)+z13^cor.mat(1,3))/2, with z11 z12 z13 are gamma with mean and variance 1 , first random effect is z1 and for N1 second random effect is z2 and for N2 third random effect is for death  
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param cens rate of censoring exponential distribution
##' @param extend to extend hazards to max-time 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##' library(mets)
##' data(CPH_HPN_CRBSI)
##' dr <- CPH_HPN_CRBSI$terminal
##' base1 <- CPH_HPN_CRBSI$crbsi 
##' base4 <- CPH_HPN_CRBSI$mechanical
##' dr2 <- scalecumhaz(dr,1.5)
##' cens <- rbind(c(0,0),c(2000,0.5),c(5110,3))
##'
##' iddata <- simMultistate(100,base1,base1,dr,dr2,cens=cens)
##' dlist(iddata,.~id|id<3,n=0)
##'  
##' ### estimating rates from simulated data  
##' c0 <- phreg(Surv(start,stop,status==0)~+1,iddata)
##' c3 <- phreg(Surv(start,stop,status==3)~+strata(from),iddata)
##' c1 <- phreg(Surv(start,stop,status==1)~+1,subset(iddata,from==2))
##' c2 <- phreg(Surv(start,stop,status==2)~+1,subset(iddata,from==1))
##' ###
##' par(mfrow=c(2,3))
##' plot(c0)
##' lines(cens,col=2) 
##' plot(c3,main="rates 1-> 3 , 2->3")
##' lines(dr,col=1,lwd=2)
##' lines(dr2,col=2,lwd=2)
##' ###
##' plot(c1,main="rate 1->2")
##' lines(base1,lwd=2)
##' ###
##' plot(c2,main="rate 2->1")
##' lines(base1,lwd=2)
##'  
##' @aliases simMultistateII
##' @export
simMultistate <- function(n,cumhaz,cumhaz2,death.cumhaz,death.cumhaz2,
		    rr=NULL,rr2=NULL,rd=NULL,rd2=NULL,rrc=NULL,
		    gap.time=FALSE,max.recurrent=100,
		    dependence=0,var.z=0.22,cor.mat=NULL,cens=NULL,extend=TRUE,...) 
{# {{{

  fdeath <- dtime <- NULL # to avoid R-check 
  status <- dhaz <- NULL; dhaz2 <- NULL

  if (dependence==0) { z <- z1 <- z2 <- zd  <- zd2 <-  rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
###	      z <- exp(rnorm(n,1)*var.z[1]^.5)
	      z1 <- z; z2 <- z; zd <- z
	      if (!is.null(cor.mat)) { zd <- rep(1,n); }
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
###	      print(summary(z))
###	      print(cor(z))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
###	      print(summary(z))
###	      print(cor(z))
      } else stop("dependence 0-3"); # }}}

  ## covariate adjustment 
  if (is.null(rr))  rr <- z1; 
  if (is.null(rr2)) rr2 <- z2; 
  if (is.null(rd))  rd  <- zd; 
  if (is.null(rd2)) rd2 <- zd2; 

  ll <- nrow(cumhaz)
  ### extend of cumulatives
  cumhaz <- rbind(c(0,0),cumhaz)
  cumhaz2 <- rbind(c(0,0),cumhaz2)
  death.cumhaz <- rbind(c(0,0),death.cumhaz)
  death.cumhaz2 <- rbind(c(0,0),death.cumhaz2)

  haz <- haz2 <- NULL
  ## range max of cumhaz and cumhaz2 

  if (!is.null(cens)) {
      if (is.matrix(cens))  {
      out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz,death.cumhaz2,cens),NULL,extend=extend)
      cens <- out$cum5
      }
  } else {
     out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz,death.cumhaz2),NULL,extend=extend)
  }
  cumhaz <- out$cum1
  cumhaz2 <- out$cum2
  cumhazd <- out$cum3
  cumhazd2 <- out$cum4
  max.time <- tail(cumhaz[,1],1)

  tall <- rcrisk(cumhaz,cumhazd,rr,rd,cens=cens,extend=NULL)

  tall$id <- 1:n
  ### fixing the first time to event
  tall$death <- 0
  ### cause 2 is death state 3, cause 1 is state 2
  tall <- dtransform(tall,status=3,status==2)
  tall <- dtransform(tall,death=1,status==3)
  tall <- dtransform(tall,status=2,status==1)
  ### dead or censored
  deadid <- (tall$status==3 | tall$status==0)
  tall$from <- 1
  tall$to <- tall$status
  ## id's that are dead: tall[deadid,]
  ## go furhter with those that are not yet dead  or censored
  tt <- tall[!deadid,,drop=FALSE]
  ## also check that we are before max.time
  tt <- subset(tt,tt$time<max.time)

  i <- 1; 
  while ( (nrow(tt)>0) & (i < max.recurrent)) {# {{{
	  i <- i+1
	  nn <- nrow(tt)

	  z1r <- rr[tt$id]
	  zdr <- rd[tt$id]
	  z2r <- rr2[tt$id]
	  zd2r <- rd2[tt$id]

	  if (i%%2==0) { ## in state 2
	  ## out of 2 for those in 2
          tt1 <- rcrisk(cumhaz2,cumhazd2,z2r,zd2r,entry=tt$time,cens=cens,extend=NULL)
          tt1$death <- 0
	  ### status 2 is death state 3, status 1 is state 1
	  tt1 <- dtransform(tt1,status=3,status==2)
	  tt1 <- dtransform(tt1,death=1,status==3)
	  tt1$from <- 2
	  tt1$to <- tt1$status
	  ## take id from tt
	  tt1$id <-  tt$id
	  ### add to data 
	  tall <- rbind(tall,tt1,row.names=NULL)

          deadid <- (tt1$status==3 | tt1$status==0)
	  ### those that are still under risk 
	  tt <- tt1[!deadid,,drop=FALSE]
	  ## also keep only those before max.time
          tt <- subset(tt,tt$time<max.time)
		  
	  } else { ## in state 1
	  ## out of 1 for those in 1
          tt1 <- rcrisk(cumhaz,cumhazd,z1r,zdr,entry=tt$time,cens=cens,extend=NULL)

          tt1$death <- 0
	  ### status 2 is death state 3, status 1 is state 2
	  tt1 <- dtransform(tt1,status=3,status==2)
	  tt1 <- dtransform(tt1,death=1,status==3)
	  tt1 <- dtransform(tt1,status=2,status==1)
	  tt1$from <- 1
	  tt1$to <- tt1$status
	  tt1$id <-  tt$id

	  ### add to data 
	  tall <- rbind(tall,tt1,row.names=NULL)

	  ## take id from tt
          deadid <- (tt1$status==3 | tt1$status==0)
	  ### those that are still under risk 
	  tt <- tt1[!deadid,,drop=FALSE]
	  ### also only keep those before max.time
          tt <- subset(tt,tt$time<max.time)
	  }

  }  # }}}

  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"death.cumhaz2") <- cumhazd2
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"cens.cumhaz") <- cens
  attr(tall,"z") <- z

  return(tall)
  }# }}}


