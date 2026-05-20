
## {{{ hazard simulation

##' Simulation of Piecewise Constant Hazard Model (Cox)
##'
##' Simulates data from a piecewise constant baseline hazard that can also be of Cox type. 
##' Censors data at the highest value of the break points.
##' 
##' For a piecewise linear cumulative hazard, the inverse is easy to compute. With delayed 
##' entry \eqn{x}, we compute:
##' \deqn{\Lambda^{-1}(\Lambda(x) + E/RR)}
##' where \eqn{RR} are the relative risks and \eqn{E} is exponential with mean 1.
##' This quantity has survival function:
##' \deqn{P(T > t | T>x) = \exp(-RR (\Lambda(t) - \Lambda(x)))}
##'
##' @param cumhazard Cumulative hazard matrix (columns: time, cumulative hazard), or 
##'   piece-constant rates for periods defined by the first column of input.
##' @param rr Relative risk vector for simulations, or alternatively when \code{rr=1} 
##'   specify \code{n}.
##' @param n Number of simulations if \code{rr} not given.
##' @param entry Delayed entry time for simulations (optional).
##' @param cause Name/code of the event cause (default 1).
##' @param extend To extend piecewise constant with constant rate beyond the last break point. 
##'   Default is \code{FALSE}. If \code{TRUE}, extends with average rate over time from 
##'   cumulative. If numeric, uses the given rate.
##' @return A data frame containing:
##'   \item{entry}{Entry times.}
##'   \item{time}{Event/censoring times.}
##'   \item{status}{Event status (1=event, 0=censored).}
##'   \item{rr}{Relative risks used.}
##'   
##'   Attributes include:
##'   \item{cumhaz}{The cumulative hazard used.}
##'   \item{extend.rate}{The extension rate if used.}
##' @author Thomas Scheike
##' @keywords survival simulation
##' @examples
##' chaz <-  c(0,1,1.5,2,2.1)
##' breaks <- c(0,10,   20,  30,   40)
##' cumhaz <- cbind(breaks,chaz)
##' n <- 10
##' X <- rbinom(n,1,0.5)
##' beta <- 0.2
##' rrcox <- exp(X * beta)
##' 
##' pctime <- rchaz(cumhaz,n=10)
##' pctimecox <- rchaz(cumhaz,rrcox,entry=runif(n))
#' @export 
#' @aliases sim_rchaz lin_approx 
rchaz <- function(cumhazard,rr,n=NULL,entry=NULL,cause=1,extend=FALSE)
{# {{{
  if (!is.null(n)) rr <- rep(1,n)
  n <- length(rr)

  breaks <- cumhazard[,1]
  mm <- tail(breaks,1)
  cumh <- cumhazard[,2] 
   ttt <- rexp(n)/rr
   if (cumhazard[1,2]>0)  { ## start cumulative hazard with a 0
###   warning("Safest to start with cumulative hazard 0 to avoid problems\n"); 
      cumhazard <- rbind(c(0,0),cumhazard)
      cumh <- c(0,cumh)
   }
   ###
   if (!is.null(entry)) {
	   if (length(entry)==1) entry <- rep(entry,n) 
	   cumentry <- lin_approx(entry,cumhazard,x=1)
	   if (any(entry>tail(breaks,1)))  stop("Some entry times further out than last cumulative hazard time\n"); 
   } else { entry <- cumentry <- rep(0,n) }
   ###
   ttte <- ttt+cumentry
   rrx <- lin_approx(ttte,cumhazard,x=-1)
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

##' Multiple Cause Piecewise Constant Hazard Simulation
##'
##' Simulates data from multiple piecewise constant baseline hazards for competing risks. 
##' Takes the minimum of all cause-specific event times and assigns the corresponding cause.
##'
##' @param cumhaz List of cumulative hazard matrices, one for each cause.
##' @param rr Matrix of relative risks (rows = subjects, columns = causes).
##' @param causes Vector of cause codes to assign (default NULL, uses 1,2,...).
##' @param ... Additional arguments passed to \code{rchaz}.
##' @return A data frame with event times and status indicating the cause of the first event.
##' @author Thomas Scheike
##' @seealso \code{\link{rchaz}}
##' @export 
rchazl <- function (cumhaz, rr, causes=NULL,...)
{# {{{
    l <- length(cumhaz)
    tall <- rchaz(cumhaz[[1]], rr[, 1], ...)
    if (l >= 2)
        for (i in 2:l) {
            tall2 <- rchaz(cumhaz[[i]], rr[, i], ...)
###            tall$status <- ifelse(tall$time < tall2$time, tall$status, i * tall2$status)
	    tall$status <- ifelse(tall$time < tall2$time, tall$status,
                      ifelse(tall2$status == 0, 0, i))
            tall$time <- pmin(tall$time, tall2$time)
        }
 if (!is.null(causes)) {
      where <- which(tall$status!=0) 
      tall$status[where] <- causes[tall$status[where]]
 }
    return(tall)
}# }}}

#' @export
lin_approx <- function (x2, xfx, x = 1)
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

addCums <- function(cumB,cumA,max=NULL)
{# {{{
 ## take times
 times <- sort(unique(c(cumB[,1],cumA[,1])))
 if (!is.null(max)) times <- times[times<max]
 cumBjx <- lin_approx(times,cumB,x=1)
 cumAjx <- lin_approx(times,cumA,x=1)
 cumBA <- cumBjx+cumAjx
 return(cbind(times,cumBA))
}# }}}

extendCums <- function(cumA,cumB,extend=NULL)
{# {{{
## setup as list to run within loop
if (!is.null(cumB)) { cumA <- list(cumA,cumB); } 

## to also work with strata version where each list contains a list of cumHaz for strata
###matrixlist <- any(unlist(lapply(cumA, function(x) is.data.frame(x[[1]]) | is.matrix(x[[1]]))))
## any stratified components 
basecumhaz <- any(unlist(lapply(cumA, function(x) inherits(x,"basecumhaz"))))
nn <- length(cumA)
nl <- lengths(cumA)
if (basecumhaz) 
	cumA <- unlist(cumA, recursive = FALSE)

restore <- function(flat, lengths) {
  ends <- cumsum(lengths)
  starts <- c(1, head(ends + 1, -1))
  mapply(function(s, e) flat[s:e], starts, ends, SIMPLIFY = FALSE)
}

 maxx <- unlist(lapply(cumA,function(x) tail(x,1)[1]))
 mm <- which.max(maxx)
 haza <- NULL
 if (is.numeric(extend)) 
    haza <- rep(extend, length.out = length(cumA))

 ## extend all that are not at maxtime
for (i in seq(length(cumA))[-mm]) {
  cumB <- cumA[[i]]; 
  if (cumB[1,1]!=0 & cumB[1,2]!=0) cumB <- rbind(0,cumB); 

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

 if (basecumhaz) cumA <- restore(cumA,nl)
 return( setNames(cumA,paste("cum",seq(nn),sep="")))
}# }}}

#' @export
simrchaz <- function(cumhazard,rr,n=NULL,cens=NULL,rrc=NULL,entry=NULL,...)
{# {{{
###   adds censoring to to rchaz call
   if (!is.null(n)) rr <- rep(1,n)
   n <- length(rr)  

   if (is.list(cumhazard)) ptt <- rchazl(cumhazard,rr,entry=entry,...) else 
	   ptt <- rchaz(cumhazard,rr,entry=entry,...)

   if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
      dt <- data.frame(time=pmin(ptt$time,pct),status=ifelse(ptt$time<pct,ptt$status,0))
   } else dt <- ptt[,c("time","status")]

   return(dt)
}# }}}

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
#' cumulatives, see also sim_phregs  
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
#' par(mfrow=c(1,3))
#' scox0 <- phreg(Surv(time,status==0)~tcell+platelet,data=dd)
#' plot(scox0); lines(cbind(c(1,30,68),c(.01,1,3)),col=2)
#' ##
#' scox1 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
#' scox2 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
#' plot(cox1); plot(scox1,add=TRUE,col=2)
#' plot(cox2); plot(scox2,add=TRUE,col=2)
#' cbind(cox1$coef,scox1$coef,cox2$coef,scox2$coef)
#' 
#' # 3 causes and censoring 
#' d3 <-  rcrisk(list(cox1$cum,cox2$cum,cox1$cum),NULL,n=100,cens=cbind(c(1,30,68),c(.01,1,3)))
#' dtable(d3,~status)
#' 
#' @export 
rcrisk <-function(cumA,cumB,rr1=NULL,rr2=NULL,n=NULL,
		  cens=NULL,rrc=NULL,extend=TRUE,causes=NULL,...)
{# {{{
	if (!is.null(cumB)) {
		cumA <- list(cumA,cumB);
	} else if (!is.list(cumA)) cumA <- list(cumA)
	l <- length(cumA)

    ## --- resolve n and build rr matrix -------------------------------
    ## collect any supplied rr vectors into a list, dropping NULLs
    rr_supplied <- Filter(Negate(is.null), list(rr1, rr2))

    if (length(rr_supplied) == 0) {
        ## no rr given: need n explicitly
        if (is.null(n)) stop("supply 'n' or at least one of 'rr1', 'rr2'\n")
        rr <- matrix(1, n, l)
    } else {
        ## derive n from the first supplied rr
        n <- length(rr_supplied[[1]])
        ## build matrix: supplied columns first, remainder filled with 1
        rr_cols <- vector("list", l)
        rr_cols[[1]] <- if (!is.null(rr1)) rr1 else rep(1, n)
        rr_cols[[2]] <- if (!is.null(rr2)) rr2 else rep(1, n)
        if (l > 2)
            for (i in seq(3, l)) rr_cols[[i]] <- rep(1, n)
        rr <- do.call(cbind, rr_cols)
    }

	if (!is.null(extend))  cumA <- extendCums(cumA,NULL,extend=extend)

	l <- length(cumA)
	## simulate all hazards  for all causes
        ptt <- rchazl(cumA,rr,causes=causes,...)

	## add censoring 
	if (!is.null(cens)) {
		pct <- simCens(cens,rrc=rrc,n=n,...)
		ptt$time <- pmin(ptt$time,pct)
		ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
	}

	return(ptt)
}# }}}

##' Simulation of Output from Cox Model
##'
##' Simulates data that looks like fit from a Cox model. Automatically censors data 
##' for the highest value of the event times by using cumulative hazard.
##'
##' @param cox Output from \code{coxph} or \code{phreg} model fitting.
##' @param n Number of simulations.
##' @param data Data frame to extract covariates for simulations (draws from observed covariates).
##' @param Z Design matrix instead of data.
##' @param rr Vector of relative risks for Cox model.
##' @param strata Vector of strata.
##' @param entry Delayed entry variable for simulation.
##' @param extend Extend possible stratified baselines to largest endpoint.
##' @param cens Censoring specification (matrix = cumulative hazard, scalar = rate).
##' @param rrc Relative risks for Cox-type censoring.
##' @param ... Arguments for \code{rchaz} (e.g., entry-time).
##' @return Data frame with simulated event times, status, and covariates.
##' @author Thomas Scheike
##' @keywords survival simulation
#' @examples
#' data(sTRACE)
#' nsim <- 100
#' coxs <-  phreg(Surv(time,status==9)~strata(chf)+vf+wmi,data=sTRACE)
#' set.seed(100)
#' sim3 <- sim_phreg(coxs,nsim,data=sTRACE)
#' head(sim3)
#' cc <-   phreg(Surv(time,status)~strata(chf)+vf+wmi,data=sim3)
#' cbind(coxs$coef,cc$coef)
#' plot(coxs,col=1); plot(cc,add=TRUE,col=2)
#' 
#' Z <- sim3[,c("vf","chf","wmi")]
#' strata <- sim3[,c("chf")]
#' rr <- exp(as.matrix(Z[,-2]) %*% coef(coxs))
#' sim4 <- sim_phreg(coxs,nsim,data=NULL,rr=rr,strata=strata)
#' sim4 <- cbind(sim4,Z)
#' cc <-   phreg(Surv(time,status)~strata(chf)+vf+wmi,data=sim4)
#' cbind(coxs$coef,cc$coef)
#' plot(coxs,col=1); plot(cc,add=TRUE,col=2)
#' 
#' @export sim_phreg 
sim_phreg <- function(cox,n,data=NULL,Z=NULL,rr=NULL,strata=NULL,
		      entry=NULL,extend=TRUE,cens=NULL,rrc=NULL,...)
{# {{{
	if  (!is.null(data)) {
		scox1 <- draw_phreg(cox,n,data=data,onlyX=TRUE,...)
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
		if (is.null(rr) & (!is.null(Z))) rr <- exp(as.matrix(Z) %*% cox$coef)
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
		if (!is.null(entry)) 
			lentry <- entry[whichi]
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

##### add correct names to entry,time,status
if (inherits(cox,"phreg"))  {
varsY <- all.vars(update(drop.specials(cox$formula,"cluster"),.~1)) 
if (length(varsY)==3) ptt[,varsY] <- cbind(ptt$entry,ptt$time,ptt$status)
if (length(varsY)==2)  {
	if (!is.null(entry)) { 
		varsY <- c("entry",varsY);  
	        ptt[,varsY] <- cbind(ptt$entry,ptt$time,ptt$status) 
	} else ptt[,varsY] <- cbind(ptt$time,ptt$status)
}
ptt <- dkeep(ptt,varsY)
}

if (!is.null(dat)) ptt <- cbind(ptt,dat)

return(ptt)
}# }}}


draw_phreg <- function(cox,n,data=NULL,Z=NULL,strata=NULL,
                       drawZ=TRUE,fixZ=FALSE,id=NULL,onlyX=TRUE)
{# {{{
###if (!inherits(cox,"phreg")) stop("must be phreg object\n"); 

## give data so that design can be constructed based on model-formula
if (!is.null(data)) {
cid <- countID(data.frame(id=cox$id))
whereid <- which(cid$Countid==1)
if (drawZ==TRUE) xid <- sample(whereid,n,replace=TRUE) else xid <- id
if (onlyX) vars <- all.vars(update(drop.specials(cox$formula,"cluster"),-1~.)) else vars <- all.vars(cox$formula)
dataid <- data[xid,vars,drop=FALSE] 

desX <- readPhreg(cox,dataid,data=FALSE)
Z <- desX$X
strata <- desX$strata
###   xx <- update_design(cox,data = dataid,response=FALSE) 
###   Z <- xx$x
###   if (!is.null(xx$strata)) strata <- as.numeric(xx$strata)-1 else strata <- rep(0,nrow(Z))
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

out <- list(Z=Z,cumhaz=cumhaz,rr=rr,id=xid,model=model,strata=strata,data=dataid,stratname=stratname)

return(out)
} ## }}}

draw_phregs <- function(coxs,n,data,onlyX=TRUE,...) { ## {{{ 
   scox1 <- draw_phreg(coxs[[1]],n,data=data,onlyX=onlyX,...)
   datas <- scox1$data
   stratam <-  scox1$strata
   rrm <- scox1$rr

   if (length(coxs)>1) 
   for (i in 2:length(coxs)) {
      coxn <- draw_phreg(coxs[[i]],n,data=data,drawZ=FALSE,id=scox1$id,onlyX=onlyX,...)
      coxndata <- coxn$data
      ind <-  match(colnames(datas),colnames(coxndata),nomatch=0)
      ind <- ind[ind!=0]
      if (length(ind)>0)  datas <- cbind(datas,coxndata[,-ind,drop=FALSE]) else datas <- cbind(datas,coxndata)
      rrm <- cbind(rrm,coxn$rr)
      stratam <- cbind(stratam,coxn$strata)
   }
   datas <- data.frame(datas)
   datas$orig.id <- scox1$id
   out <- list(data=datas,rr=rrm,strata=stratam)

   return(out)
} ## }}} 

setup_phreg  <- function(cumhazard,coef,Znames=NULL,strata=NULL)
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


##' Simulation of Cause-Specific Cox Models
##'
##' Simulates data that looks like fit from cause-specific Cox models. Censors data 
##' automatically. When censoring is given in the list of causes, this provides censoring 
##' that looks like the data.
##'
##' @param coxs List of Cox models.
##' @param n Number of simulations.
##' @param data Data frame to extract covariates.
##' @param rr Relative risks.
##' @param strata Strata vector.
##' @param entry Delayed entry.
##' @param extend Extend baselines to largest endpoint.
##' @param cens Censoring specification.
##' @param rrc Relative risks for censoring.
##' @param ... Arguments for \code{rchaz}.
##' @return Data frame with simulated event times, status, and covariates.
##' @author Thomas Scheike
##' @keywords survival simulation
#' @examples
#' data(bmt)
#' nsim <- 100; 
#' 
#' cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet+age,data=bmt)
#' cox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=bmt)
#' coxs <- list(cox1,cox2)
#' ## just calls sim_phregs !
#' dd <- sim_phregs(coxs,nsim,data=bmt,extend=c(0.001))
#' scox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet+age,data=dd)
#' scox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=dd)
#'
#' cbind(cox1$coef,scox1$coef)
#' cbind(cox2$coef,scox2$coef)
#' par(mfrow=c(1,2))
#' plot(cox1); plot(scox1,add=TRUE); 
#' plot(cox2); plot(scox2,add=TRUE); 
#' 
#' @export sim_phregs
sim_phregs <- function(coxs,n,data=NULL,rr=NULL,strata=NULL,
                       entry=NULL,extend=TRUE,cens=NULL,rrc=NULL,...)
{# {{{
   out <- draw_phregs(coxs,n,data)
   if (is.null(rr)) rr <- out$rr
   if (is.null(strata)) strata <- out$strata
   lentry <- NULL

   cumhazl <- list()
   cumhazl[[1]] <- basecumhaz(coxs[[1]],only=1)
   if (length(coxs)>1) 
   for (i in 2:length(coxs)) 
      cumhazl[[i]] <- basecumhaz(coxs[[i]],only=1)
   if (!is.null(extend)) cumhazl <- extendCums(cumhazl,NULL,extend=extend)

   ## simulate first  time
   simdata <- sim_phreg(cumhazl[[1]],n,data=NULL,rr=rr[,1],strata=strata[,1],entry=entry)
   l <- length(coxs)
   if (l>=2) 
   for (i in 2:l) {
      tall2 <- sim_phreg(cumhazl[[i]],n,rr=rr[,i],strata=strata[,i],entry=entry)
      simdata$status <- ifelse(simdata$time<tall2$time,simdata$status,i*tall2$status)
      simdata$time <- pmin(simdata$time,tall2$time)
   }
   ptt <- simdata

 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,entry=entry,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

##### add correct names to entry,time,status
if (inherits(coxs[[1]],"phreg"))  {
varsY <- all.vars(update(drop.specials(coxs[[1]]$formula,"cluster"),.~1)) 
if (length(varsY)==3) ptt[,varsY] <- cbind(ptt$entry,ptt$time,ptt$status)
if (length(varsY)==2)  {
	if (!is.null(entry)) { 
		varsY <- c("entry",varsY);  
	        ptt[,varsY] <- cbind(ptt$entry,ptt$time,ptt$status) 
	} else ptt[,varsY] <- cbind(ptt$time,ptt$status)
}
ptt <- dkeep(ptt,varsY)
}

if (!is.null(data)) ptt <- cbind(ptt,out$data)

return(ptt)
}# }}}


## }}} 

## {{{ cumulative incidence 

##' @export sim_cifs
sim_cifs <- function(cifs,n,data=NULL,rr=NULL,strata=NULL,Z=NULL,
                     cens=NULL,rrc=NULL,max.times=NULL,causes=c(1,2),
                     U=NULL,pU=NULL,extend=TRUE,type=NULL,restrict=TRUE,entry=NULL,...)
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

   out <- draw_phregs(cifs,n,data,Z=Z)
   if (is.null(rr)) rr <- out$rr
   if (is.null(strata)) strata <- out$strata

   cumhazl <- list()
   cumhazl[[1]] <- basecumhaz(cifs[[1]],only=1)
   if (length(cifs)>1) 
   for (i in 2:length(cifs)) 
      cumhazl[[i]] <- basecumhaz(cifs[[i]],only=1)
  if (!is.null(extend)) cumhazl <- extendCums(cumhazl,NULL,extend=extend)

  tau <- tail(cumhazl[[1]][[1]],1)[1]

  ## simulate first  time
  sim1 <- sim_cif(cumhazl[[1]],n,rr=rr[,1],strata=strata[,1],U=U,type=model1,entry=entry)
  sim2 <- sim_cif(cumhazl[[2]],n,rr=rr[,2],strata=strata[,2],U=U,type=model2,entry=entry)

   ## drawing which cause  is observed, 0,1,2
  if (!is.null(entry)) {
     if (restrict) {
	   F1tau <- sim1$F1tau- sim1$F1entry
	   F2tau <- sim2$F1tau-sim2$F1entry
	   F2tau <- F2tau*(1-sim1$F1tau) 
	   F2entry <- sim2$F1entry*(1-sim1$F1tau) 
           St <- (1-sim1$F1entry-F2entry)
     } else { 
	   F1tau <- sim1$F1tau; 
	   F2tau <- sim2$F1tau; 
           St <- (1-sim1$F1entry-sim2$F1entry)
     }
           F1tau <- F1tau/St 
           F2tau <- F2tau/St
  } else {
     F1tau <- sim1$F1tau
     if (restrict) F2tau <- sim2$F1tau * (1-F1tau) else F2tau <- sim2$F1tau
  }

  ptot <- F1tau+F2tau
  if (!is.null(pU)) {
      rt <- 1*(pU< pmin(ptot,1))
      rb <- 1*(pU< F1tau) 
  } else {
      rt <- rbinom(n,1,pmin(ptot,1))
      rb <- rbinom(n,1,F1tau/ptot)
  }
  cause=ifelse(rb==1,1,2)
  time=ifelse(cause==causes[1],sim1$timing,sim2$timing)
  cause <- rt*cause
  time[cause==0] <- tau

  if (is.null(entry)) entry <- 0
  ptt <- data.frame(entry=entry,time=time,status=cause,cause=cause,ptot=ptot)

### adds censoring 
 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }


##### add correct names to entry,time,status
if (inherits(cifs[[1]],"phreg"))  {
varsY <- all.vars(update(drop.specials(cifs[[1]]$formula,"cluster"),.~1)) 
if (length(varsY)==2) ptt[,varsY] <- cbind(ptt$time,ptt$status)
if (length(varsY)==3) ptt[,varsY] <- cbind(ptt$entry,ptt$time,ptt$status)
if (!is.null(entry) & length(varsY)==2) varsY <- c("entry",varsY)
###varsY <- c(varsY,"F1tau","F1entry","timing")
ptt <- dkeep(ptt,varsY)
}
if (!is.null(data)) ptt <- cbind(ptt,out$data)

return(ptt)
}# }}}

sim_subdist <- function(cumhazard,rr,n=NULL,entry=NULL,Sentry=NULL,type="cloglog",startcum=c(0,0),U=NULL,...)
{# {{{
  ## Fine-Gray model cloglog F1= 1-exp(-cum(t)*rr)
  ## logistic                F1= cum(t)*rr/(1+cum(t)*rr)
  ## rr                      F1= cum(t)*rr,  rr=exp(X^t beta) 
  if (!is.null(n)) rr <- rep(1,n)

  logit <- function(p) log(p/(1-p))

  if (cumhazard[1,2]>0)  cumhazard <- rbind(startcum,cumhazard)
  breaks <- cumhazard[,1]
  cumh <- cumhazard[,2] 
  mm <- tail(breaks,1)
  n <- length(rr)
  if (!is.null(entry)) {
       cumhentry <- lin_approx(entry,cumhazard,x=1)
  } else cumhentry <- 0
  F1entry <- 0

  if (is.null(U)) U <- runif(n)

  if (type=="cloglog") {
      if (!is.null(entry)) F1entry <- 1-exp(-cumhentry*rr) 
      F1tau <- 1-exp(-tail(cumh,1)*rr)-F1entry
      ttt <- -log(1-U*F1tau-F1entry)/rr
  } else if (type=="logistic") {
     if (!is.null(entry)) F1entry <- cumhentry*rr/(1+cumhentry*rr) 
     F1tau <- tail(cumh,1)*rr/(1+tail(cumh,1)*rr)-F1entry
     v <- U*F1tau+F1entry 
     ttt <- exp(logit(v))/rr; 
  }  else if (type=="rr" | type=="cif") {
     if (!is.null(entry)) F1entry <- cumhentry
     F1tau <- tail(cumh,1) - F1entry
     ttt <- U*F1tau-F1entry
     ## rr only affects binomial draw 
     F1tau <- (F1tau)*rr 
     F1entry <- (F1entry)*rr 
  } else stop(" cloglog or logistic or give function (fun=) \n"); 
  ###
   rrx <- lin_approx(ttt,cumhazard,x=-1)
  if (is.null(entry)) entry <- 0
  timecause <- rrx
  ###
  rrx <- ifelse(timecause>mm,mm,rrx)

   if (!is.null(entry))  {
      if (is.null(Sentry)) Sentry <- 1-F1entry 
   } else Sentry <- 1

   ## drawing event type using conditional distribution
   status <- rbinom(n,1,F1tau/Sentry) 
   rrx[status==0] <- mm
   dt <- data.frame(entry=entry,time=rrx,status=status,rr=rr,
		    F1tau=F1tau+F1entry,F1entry=F1entry,
		    timing=timecause)
   attr(dt,"cumhaz") <- cumhazard
   return(dt)
}# }}}

calcCIF <- function(cumhaz,tau,rr,entry=NULL,type=c("cloglog","logistic","rr")) { ## {{{

  if (!is.null(entry)) {
       cumhentry <- lin_approx(entry,cumhaz,x=1)
  } else F1entry <- NULL

  if (type[1]=="cloglog") {
      if (!is.null(entry)) F1entry <- 1-exp(-cumhentry*rr) 
      F1tau <- 1-exp(-tail(cumh,1)*rr) 
  } else if (type[1]=="logistic") {
     if (!is.null(entry)) F1entry <- cumhentry*rr/(1+cumhentry*rr) 
     F1tau <- tail(cumh,1)*rr/(1+tail(cumh,1)*rr) 
  }  else if (type[1]=="rr" | type=="cif") {
     if (!is.null(entry)) F1entry <- cumhentry
     F1tau <- tail(cumh,1)
  } else stop(" cloglog or logistic or give function (fun=) \n"); 

  return(list(F1tau=F1tau,rr=rr,F1entry=F1entry,entry=entry))
} ## }}} 

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
  rrx  <- lin_approx(u,F1,x=-1)
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
  rrx <-   lin_approx(times,F1,x=1)
  dt <- cbind(times,rrx)
  colnames(dt) <- c("time","subdist")
  return(dt)
}# }}}

##' Simulation of Output from Cumulative Incidence Regression Model
##'
##' Simulates data that looks like fit from a fitted cumulative incidence model 
##' (Fine-Gray or logistic).
##'
##' @param cif Output from \code{prop.odds.subdist} or \code{ccr} (\code{cmprsk}), 
##'   or call \code{invsubdist} with cumulative and linear predictor.
##' @param n Number of simulations.
##' @param data Data frame to extract covariates.
##' @param Z Design matrix instead of data.
##' @param rr Relative risks.
##' @param strata Strata vector.
##' @param drawZ Logical; randomly sample from Z.
##' @param cens Censoring specification.
##' @param rrc Relative risks for censoring.
##' @param entry Delayed entry time.
##' @param Sentry Survival related to delayed entry.
##' @param cumstart Start cumulatives at time 0.
##' @param U Uniforms for drawing timing.
##' @param pU Uniforms for drawing event type.
##' @param type Model type: \code{"logistic"}, \code{"cloglog"}, or \code{"rr"}.
##' @param extend Extend piecewise constant with constant rate.
##' @param ... Arguments for \code{sim_subdist}.
##' @return Data frame with simulated event times, status, and covariates.
##' @author Thomas Scheike
##' @keywords survival simulation
##' @examples
#' data(bmt)
#' nsim <- 100
#' 
#' ## logit cumulative incidence regression model 
#' cif <- cifreg(Event(time,cause)~platelet+age,data=bmt,cause=1)
#' estimate(cif)  
#' plot(cif,col=1)
#' simbmt <- sim_cif(cif,nsim,data=bmt)
#' dtable(simbmt,~cause)
#' scif <- cifreg(Event(time,cause)~platelet+age,data=simbmt,cause=1)
#' estimate(scif)
#' plot(scif,add=TRUE,col=2)
#' 
#' ## Fine-Gray cloglog cumulative incidence regression model 
#' cif <- cifregFG(Event(time,cause)~strata(tcell)+age,data=bmt,cause=1)
#' estimate(cif)  
#' plot(cif,col=1)
#' simbmt <- sim_cif(cif,nsim,data=bmt)
#' scif <- cifregFG(Event(time,cause)~strata(tcell)+age,data=simbmt,cause=1)
#' estimate(scif)
#' plot(scif,add=TRUE,col=2)
#' 
#' ################################################################
#' #  simulating several causes with specific cumulatives 
#' ################################################################
#' cif1 <-  cifreg(Event(time,cause)~strata(tcell)+age,data=bmt,cause=1)
#' cif2 <-  cifreg(Event(time,cause)~strata(platelet)+tcell+age,data=bmt,cause=2)
#' cifss <-  list(cif1,cif2)
#' simbmt <- sim_cifs(list(cif1,cif2),nsim,data=bmt,extend=0.005)
#' scif1 <-  cifreg(Event(time,cause)~strata(tcell)+age,data=simbmt,cause=1)
#' scif2 <-  cifreg(Event(time,cause)~strata(platelet)+tcell+age,data=simbmt,cause=2)
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
#' @seealso [simul_cifs()]
#' @aliases sim_cif sim_cifs subdist sim_subdist invsubdist
#' @export sim_cif
sim_cif <- function(cif,n,data=NULL,Z=NULL,rr=NULL,strata=NULL,
                    drawZ=TRUE,cens=NULL,rrc=NULL,entry=NULL,Sentry=NULL,
		    cumstart=c(0,0),U=NULL,pU=NULL,type=NULL,extend=NULL,...)
{# {{{
## also extracts coefficients and baseline from cifreg
if (!is.null(data))  {
   des <- 	draw_phreg(cif,n,data=data,Z=Z,drawZ=drawZ,...)
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

lentry <- Slentry <- NULL

## {{{ simulation of cif 
    ptt <- c()
    for (i in unique(strata)) {
      whichi <- which(strata==i)
      cumhazj <- rbind(0,cumhaz[[i+1]])
      if (!is.null(entry)) lentry <- entry[whichi]
      if (!is.null(Sentry)) Slentry <- Sentry[whichi]
      if (!is.null(U)) Ui <- U[whichi] else Ui <- NULL
      simj <- sim_subdist(cumhazj,rr[whichi],type=type,U=Ui,entry=lentry,Sentry=Slentry)
      simj$id <- ids[whichi]
      ptt  <-  rbind(ptt,simj)
    }
    dsort(ptt) <- ~id
### }# }}}

 ### adds censoring 
 if (!is.null(cens)) {
      pct <- simCens(cens,rrc=rrc,n=n,...)
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
 }

##### add correct names to entry,time,status if "phreg/cifreg/recreg" object
if (inherits(cif,"phreg"))  {
varsY <- all.vars(update(drop.specials(cif$formula,"cluster"),.~1)) 
if (length(varsY)==2) ptt[,varsY] <- cbind(ptt$time,ptt$status)
if (length(varsY)==3) ptt[,varsY] <- cbind(ptt$entry,ptt$time,ptt$status)
if (!is.null(entry) & length(varsY)==2) varsY <- c("entry",varsY)
varsY <- c(varsY,"F1tau","F1entry","timing")
ptt <- dkeep(ptt,varsY)
} 
if (!is.null(dats)) ptt <- cbind(ptt,dats)

attr(ptt,"type") <- type
attr(ptt,"znames") <- znames
attr(ptt,"cumhaz") <-  cumhaz
return(ptt)
}# }}}


##' @export
simul_cifs <- function(n,rho1,rho2,beta,rc=0.5,depcens=0,rcZ=0.5,
               bin=1,type=c("cloglog","logistic"),rate=1,entry=NULL,
               Z=NULL,U=NULL,pU=NULL,...) {# {{{
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

    cif1 <- setup_cif(cbind(tt,Lam1),beta[1:p],Znames=colnames(Z),type=type[1])
    cif2 <- setup_cif(cbind(tt,Lam2),beta[(p+1):(2*p)],Znames=colnames(Z),type=type[1])

    data <- sim_cifs(list(cif1,cif2),n,Z=Z,U=U,pU=pU,type=type[1],extend=NULL,entry=entry,...)

    if (is.null(entry)) entry <- 0
    if (!is.null(rc)) {
    if (depcens==0) censor=pmin(entry+rexp(n,1)*(1/rc),6) else censor=pmin(entry+rexp(n,1)*(1/(rc*exp(Z %*% rcZ))),6)
    } else censor <- 6 

    data$status=data$status*(data$time<=censor)
    data$time=pmin(data$time,censor)
    data <- cbind(data,Z)
    attr(data,"Lam1") <- cbind(tt,Lam1)
    attr(data,"Lam2") <- cbind(tt,Lam2)
    return(data)
}# }}}

simul_cifsRA <- function (n, rho1, rho2, beta, rc = 0.5, depcens.R = 0, rcZ = 0.5,pCA=0.5,pCR=0.5,
    bin = 1, type = c("cloglog", "logistic"), rate = 1, Z = NULL,rc2=0.3,depcens.Adm=0)
{# {{{
    p = length(beta)/2
    tt <- seq(0, 6, by = 0.1)
    if (length(rate) == 1)
        rate <- rep(rate, 2)
    Lam1 <- rho1 * (1 - exp(-tt/rate[1]))
    Lam2 <- rho2 * (1 - exp(-tt/rate[2]))
    if (length(bin) == 1)
        bin <- rep(bin, 2)
    if (length(rcZ) == 1)
        rcZ <- c(rcZ, 0)
    if (is.null(Z))
        Z = cbind((bin[1] == 1) * (2 * rbinom(n, 1, 1/2) - 1) +
            (bin[1] == 0) * rnorm(n), (bin[2] == 1) * (rbinom(n,
            1, 1/2)) + (bin[2] == 0) * rnorm(n))
    colnames(Z) <- paste("Z", 1:2, sep = "")
    p <- ncol(Z)
    cif1 <- setup_cif(cbind(tt, Lam1), beta[1:p], Znames = colnames(Z),
        type = type[1])
    cif2 <- setup_cif(cbind(tt, Lam2), beta[(p + 1):(2 * p)],
        Znames = colnames(Z), type = type[1])
    data <- sim_cifs(list(cif1, cif2), n, Z = Z,extend=NULL)

    admcens <- rbinom(n,1,pCA)
    if (depcens.Adm==1) rrA <- exp( Z %*% rcZ) else rrA <- 0
    censorA = admcens*(rrA+runif(n)*(6-rrA))+ 6*(admcens==0)

    statusA = data$status * (data$time <= censorA)
    timeA = pmin(data$time, censorA)
    statusA[statusA==0] <- 7

    Rcens <- rbinom(n,1,pCR)
    if (depcens.R == 0)  
        censorR = Rcens*runif(n)*6+(Rcens==0)*6
    else censorR = Rcens*pmin(rexp(n, 1) * (1/(rc * exp(Z %*% rcZ))),6) + (Rcens==0)*6
    censor <- ifelse( censorR < censorA, censorR,censorA)

    status = data$status * (data$time <= censor)
    time = pmin(data$time, censor)
    status[status==0 & (censorR > censorA)] <- 7
    ## extra censoring
    data <- data.frame(time=time,status=status,cens.time=censor,censorA=censorA,censorR=censorR,
		       timeA=timeA,statusA=statusA)
    return(cbind(data, Z))
}# }}}


setup_cif  <- function(cumhazard,coef,Znames=NULL,type="logistic")
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


## }}}

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
    cif1 <- setup_cif(cbind(tt, Lam1), beta[1:4], Znames = colnames(Z), type = type[1])
    cif2 <- setup_cif(cbind(tt, Lam2), beta[5:8], Znames = colnames(Z), type = type[1])
    data <- sim_cifs(list(cif2, cif1), n, Z = Z,extend=NULL,restrict=restrict,...)
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
    cif1 <- setup_cif(cbind(tt, Lam1), F1par*beta[1:4], Znames = colnames(Z), type = type[1])
    cif2 <- setup_cif(cbind(tt, Lam2), F2par*beta[5:8], Znames = colnames(Z), type = type[1])
    data <- sim_cifs(list(cif2, cif1), n, Z = Z,extend=NULL,restrict=restrict,...)
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
        rrc <- exp(as.matrix(Z) %*% c(censpar*c0$coef))
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

##' Illness-Death Competing Risks with Two Causes of Death
##'
##' Simulates data from an illness-death model with two causes of death from the 
##' illness state. Covariate effects can be introduced via relative risk terms.
##'
##' @param cumhaz Cumulative hazard from state 1 to 2.
##' @param death.cumhaz Cumulative hazard of death from state 1.
##' @param death.cumhaz2 Cumulative hazard of death from state 2.
##' @param n Number of simulations.
##' @param rr Relative risks.
##' @param rd Relative risks for death from state 1.
##' @param rd2 Relative risks for death from state 2.
##' @param gamma23 Early effect parameters for death causes.
##' @param gamma24 Early effect parameters for death causes.
##' @param early2 Time threshold for early effect.
##' @param gap.time Gap time indicator.
##' @param max.recurrent Maximum recurrent events.
##' @param cens Censoring specification.
##' @param rrc Censoring relative risks.
##' @param extend Extend hazards.
##' @param ... Additional arguments.
##' @return Data frame with multi-state event history.
##' @author Thomas Scheike
##' @keywords survival simulation
##' @export
sim_multistateII <- function(cumhaz,death.cumhaz,death.cumhaz2,n=NULL,
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
###  chaz1 <- cumss[1:3]; 
###  chaz2 <- cumss[4:5]; 
  idx1 <- 1:(1 + length(death.cumhaz))
  idx2 <- (max(idx1) + 1):(max(idx1) + length(death.cumhaz2))
  chaz1 <- cumss[idx1]; chaz2 <- cumss[idx2]
  rr1 <- cbind(rr,rd)
  rr2 <- cbind(rd2)

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

##' Simulation of Illness-Death Model
##'
##' Simulates data from a full illness-death model with reversible transitions 
##' and multiple causes of death. Supports various dependence structures via 
##' shared frailties.
##'
##' @param n Number of IDs.
##' @param cumhaz Cumulative hazard from state 1 to 2.
##' @param cumhaz2 Cumulative hazard from state 2 to 1.
##' @param death.cumhaz Cumulative hazard of death from state 1.
##' @param death.cumhaz2 Cumulative hazard of death from state 2.
##' @param rr12 Relative risk for 1->2.
##' @param rr21 Relative risk for 2->1.
##' @param rd13 Relative risk for death 1->3.
##' @param rd23 Relative risk for death 2->3.
##' @param rrc Relative risk for censoring.
##' @param gap.time Gap time indicator. If true simulates gap-times with specified cumulative hazard.
##' @param max.recurrent Maximum recurrent events.
##' @param dependence Dependence structure (0-3).
##' @param var.z Variance of random effects.
##' @param cor.mat Correlation matrix.
##' @param cens Censoring rate.
##' @param extend Extend hazards.
##' @param ... Additional arguments.
##' @return Data frame with multi-state event history.
##' @author Thomas Scheike
##' @keywords survival simulation
##' @examples
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##' data(CPH_HPN_CRBSI)
##' dr <- CPH_HPN_CRBSI$terminal
##' base1 <- CPH_HPN_CRBSI$crbsi 
##' base4 <- CPH_HPN_CRBSI$mechanical
##' dr2 <- scalecumhaz(dr,1.5)
##' cens <- rbind(c(0,0),c(2000,0.5),c(5110,3))
##'
##' iddata <- sim_multistate(100,base1,base1,dr,dr2,cens=cens)
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
##' @export
sim_multistate <- function(n,cumhaz,cumhaz2,death.cumhaz,death.cumhaz2,
		    rr12=NULL,rr21=NULL,rd13=NULL,rd23=NULL,rrc=NULL,
		    gap.time=FALSE,max.recurrent=100,
		    dependence=0,var.z=0.5,cor.mat=NULL,cens=NULL,extend=TRUE,...) 
{# {{{

  fdeath <- dtime <- NULL # to avoid R-check 
  status <- dhaz <- NULL; dhaz2 <- NULL

  if (dependence==0) { z <- z1 <- z2 <- zd  <- zd2 <-  rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
###	      z <- exp(rnorm(n,1)*var.z[1]^.5)
	      zd2 <- z; z1 <- z; z2 <- z; zd <- z
	      if (!is.null(cor.mat)) { zd <- rep(1,n); }
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; zd2 <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(4*n,1),n,4)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3]+z[,4]^cor.mat[1,4])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3]+z[,4]^cor.mat[2,4])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3]+z[,4]^cor.mat[3,4])
              zd2<- (z[,4]^cor.mat[4,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3]+z[,4]^cor.mat[4,4])
	      z <- cbind(z1,z2,zd,zd2)
      } else stop("dependence 0-3"); # }}}

  ## covariate adjustment 
  if (is.null(rr12))  rr12 <- z1  else rr12 <- rr12*z1
  if (is.null(rr21))  rr21 <- z2  else rr21 <- rr21*z2
  if (is.null(rd13))  rd13 <- zd  else rd13 <- rd13*zd
  if (is.null(rd23))  rd23 <- zd2 else rd23 <- rd23*zd2

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

  tall <- rcrisk(cumhaz,cumhazd,rr12,rd13,cens=cens,extend=NULL)

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

	  z1r <- rr12[tt$id]
	  zdr <- rd13[tt$id]
	  z2r <- rr21[tt$id]
	  zd2r <- rd23[tt$id]

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


