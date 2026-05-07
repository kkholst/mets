##' Estimation of Concordance in Bivariate Competing Risks Data
##'
##' Estimates the bivariate cumulative incidence function (concordance) for paired 
##' data (e.g., twins, family members) in the presence of competing risks. The function 
##' handles both the IPCW (Inverse Probability of Censoring Weighting) estimator and 
##' the Aalen-Johansen estimator (via \code{prodlim}).
##' 
##' The concordance function \eqn{C(t)} is defined as the probability that both members 
##' of a pair experience the event of interest by time \eqn{t}:
##' \deqn{ C(t) = P(T_1 \leq t, T_2 \leq t, \epsilon_1 = k, \epsilon_2 = k) }
##' where \eqn{T_i} is the event time and \eqn{\epsilon_i} is the cause of failure for individual \eqn{i}.
##' 
##' The function supports:
##' \itemize{
##'   \item Stratified analysis (e.g., by zygosity in twin studies).
##'   \item Clustering for robust standard errors.
##'   \item Both IPCW and Aalen-Johansen estimation methods.
##'   \item Resampling for uniform standard errors.
##'   \item IID decomposition for further inference (e.g., casewise concordance tests).
##' }
##'
##' @param formula Formula with an \code{Event} object on the left-hand side. The right-hand 
##'   side specifies the covariate structure, including \code{strata()} for grouping (e.g., MZ/DZ) 
##'   and \code{id()} for pairing.
##' @param data Data frame containing the variables.
##' @param cause Vector of cause codes for which to estimate the bivariate cumulative incidence 
##'   (default \code{c(1,1)}).
##' @param cens Censoring code (default 0).
##' @param causes Vector of all possible causes (optional, inferred from data if missing).
##' @param indiv Variable indicating individual within a pair (optional, inferred from \code{id}).
##' @param strata Variable for stratification (optional, can be specified in formula).
##' @param id Clustering variable (pair ID). Required.
##' @param num Variable for numbering individuals within pairs (optional, auto-generated if missing).
##' @param max.clust Maximum number of clusters to use for IID decomposition in \code{timereg::comp.risk}. 
##'   If NULL, uses all clusters. Useful for large datasets to speed up computation.
##' @param marg Optional marginal cumulative incidence object (from \code{comp.risk}) to compute 
##'   standard errors for same-cluster comparisons in subsequent \code{casewise.test()}.
##' @param se.clusters Vector of cluster indices or column name in \code{data} for standard error calculation. 
##'   Defaults to the \code{id} variable.
##' @param wname Name of an additional weight variable for paired competing risks data.
##' @param prodlim Logical; if TRUE, uses the \code{prodlim} (Aalen-Johansen) estimator instead of 
##'   the IPCW estimator based on \code{comp.risk}. These are equivalent in the absence of covariates.
##' @param messages Control amount of output (0 = silent, 1 = messages).
##' @param model Type of competing risk model for \code{comp.risk} (default "fg" for Fine-Gray).
##' @param return.data If 1, returns the reshaped data; if 2, returns only the data; otherwise returns the model.
##' @param uniform Logical; if TRUE, computes uniform standard errors based on resampling.
##' @param conservative Logical; if TRUE, uses conservative standard errors (recommended for large datasets).
##' @param resample.iid Logical; if TRUE, returns IID residual processes for further computations.
##' @param ... Additional arguments passed to \code{timereg::comp.risk}.
##' @return An object of class \code{"bicomprisk"} (or \code{"bicomprisk.strata"} if stratified) containing:
##'   \item{model}{List of fitted models for each stratum (or a single model).}
##'   \item{strata}{Names of strata (if applicable).}
##'   \item{N}{Number of strata (if applicable).}
##'   \item{time}{Event times.}
##'   \item{P1}{Bivariate cumulative incidence estimates.}
##'   \item{se.P1}{Standard errors of the estimates.}
##'   \item{P1.iid}{IID decomposition (if \code{resample.iid=TRUE}).}
##'   \item{clusters}{Cluster assignments (if \code{marg} or \code{se.clusters} provided).}
##' @author Thomas Scheike, Klaus K. Holst
##' @references 
##' Scheike, T. H.; Holst, K. K. & Hjelmborg, J. B. (2014). Estimating twin concordance for bivariate competing risks twin data. Statistics in Medicine, 33, 1193-1204.
##' @seealso \code{\link{bicompriskData}}, \code{\link{test_casewise}}, \code{\link{casewise}}
##' @examples
##' library("timereg")
##' 
##' ## Simulated data example
##' prt <- sim_nordic_random(2000,delayed=TRUE,ptrunc=0.7,
##'		      cordz=0.5,cormz=2,lam0=0.3)
##' ## Bivariate competing risk, concordance estimates
##' p11 <- bicomprisk(Event(time,cause)~strata(zyg)+id(id),data=prt,cause=c(1,1))
##' 
##' p11mz <- p11$model$"MZ"
##' p11dz <- p11$model$"DZ"
##' par(mfrow=c(1,2))
##' ## Concordance
##' plot(p11mz,ylim=c(0,0.1));
##' plot(p11dz,ylim=c(0,0.1));
##' 
##' ## Entry time, truncation weighting
##' ### Other weighting procedure
##' prtl <-  prt[!prt$truncated,]
##' prt2 <- ipw2(prtl,cluster="id",same.cens=TRUE,
##'      time="time",cause="cause",entrytime="entry",
##'      pairs=TRUE,strata="zyg",obs.only=TRUE)
##' 
##' prt22 <- fast.reshape(prt2,id="id")
##' 
##' prt22$event <- (prt22$cause1==1)*(prt22$cause2==1)*1
##' prt22$timel <- pmax(prt22$time1,prt22$time2)
##' ipwc <- timereg::comp.risk(Event(timel,event)~-1+factor(zyg1),
##'   data=prt22,cause=1,n.sim=0,model="rcif2",times=50:90,
##'   weights=prt22$weights1,cens.weights=rep(1,nrow(prt22)))
##' 
##' p11wmz <- ipwc$cum[,2]
##' p11wdz <- ipwc$cum[,3]
##' lines(ipwc$cum[,1],p11wmz,col=3)
##' lines(ipwc$cum[,1],p11wdz,col=3)
##' @export
bicomprisk <- function(formula, data, cause=c(1,1), cens=0, causes, indiv,
 strata=NULL, id,num, max.clust=1000, marg=NULL,se.clusters=NULL,wname=NULL,
 prodlim=FALSE,messages=TRUE,model,return.data=0,uniform=0,conservative=1,resample.iid=1,...) {
## {{{

  mycall <- match.call()
  formulaId <- unlist(Specials(formula,"id"))
  formulaIndiv <- Specials(formula,"indiv")
  formulaStrata <- unlist(Specials(formula,"strata"))
  formulaSt <- paste("~.-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")",
                     "-indiv(",paste(formulaIndiv,collapse="+"),")")
  formula <- update(formula,formulaSt)

  ### ts  11/10
  ## if (substr(as.character(formula)[2],1,4)=="Hist") {
  ##     stop("Since version : The left hand side of the formula must be specified as
  ##     Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
  ## }


  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  indiv <- formulaIndiv
  if (!is.null(formulaIndiv)) {
    mycall$indiv <- indiv
  }
  if (missing(id)) stop("Missing 'id' variable")

  ### setting up cluster level for iid decomposition
  if (is.null(se.clusters) & is.null(marg))  lse.clusters <- data[,c(id)]
  else {
    if (is.null(se.clusters)) {
      lse.clusters <- marg$clusters
###      if (!is.null(max.clust)) {   }
    } else {
      if (is.character(se.clusters)) {
        lse.clusters <- data[,se.clusters]
      } else {
        lse.clusters <- se.clusters
      }
    }
  }
  if (length(lse.clusters)!=nrow(data)) stop("'se.clusters' and 'data' does not match!")
  se.clusters.call <- se.clusters

  data <- data.frame(cbind(data,lse.clusters))

  timevar <- terms(formula)[[2]]
  ##  hh <- with(data,eval(timevar))
  if (is.call(timevar)) {
    causes <- timevar[[3]]
    timevar <- timevar[[2]]
  }
  timevar <- as.character(timevar)
  causes <- as.character(causes)

  if (!is.null(strata)) {
    fac <- interaction(data[,strata])
    dd <- split(data,fac)
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq_len(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        idx <- which(fac==names(dd)[i])
        mycall$se.clusters <- lse.clusters[idx]
        mycall$formula <- formula
        mycall$data <- dd[[i]]
        eval(mycall)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("bicomprisk.strata","twinlm.strata")
      res$N <- length(dd)
      return(res)
    }
  }

  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  ### adds weights,
  if (!is.null(wname)) covars <- c(covars,wname)

  indiv2 <- covars2 <- NULL

  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  ##which(data[,id]%in%names(idtab==2))
  data <- data[which(data[,id]%in%names(idtab==2)),]
  if (missing(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,num] <- unlist(lapply(idtab,seq_len))
  }

  oldreshape <- 0
  if (oldreshape==1) sep="." else sep=""
  timevar2 <- paste(timevar,1:2,sep=sep)
  causes2 <- paste(causes,1:2,sep=sep)
  if (length(covars)>0)
    covars2 <- paste(covars,1,sep=sep)
  for (i in seq_len(length(indiv)))
  indiv2 <- c(indiv2, paste(indiv[i],1:2,sep=sep))

  if (oldreshape==1)
  ww0 <- reshape(data[,c(timevar,causes,covars,indiv,id,num,"lse.clusters")],
         direction="wide",idvar=id,timevar=num)[,c(id,"lse.clusters.1",timevar2,causes2,indiv2,covars2)]
  else
  ww0 <- fast.reshape(data[,c(timevar,causes,covars,indiv,id,num,"lse.clusters")],id=id,num=data$num,labelnum=TRUE)[,c(id,"lse.clusters1",timevar2,causes2,indiv2,covars2)]
  ww0 <- na.omit(ww0)

  status <- rep(0,nrow(ww0))
  time <- ww0[,timevar2[1]]

  ## {{{ (i,j) causes
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cause[2])
  if (length(idx2)>0) {
      status[idx2] <- 1
      time[idx2] <- apply(ww0[idx2,timevar2[1:2],drop=FALSE],1,max)
  }

  ##(0,0), (0,j)
  idx2 <- which(ww0[,causes2[1]]==cens & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 0
      time[idx2] <- ww0[idx2,timevar2[1]]
  }

  ##(ic,0), (ic,j)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- ww0[idx2,timevar2[1]]
  }

  ##(i,0)
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cens)
  if (length(idx2)>0) {
      status[idx2] <- 0
      time[idx2] <- ww0[idx2,timevar2[2]]
  }

  ##(ic,jc)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- apply(ww0[idx2,timevar2[1:2],drop=FALSE],1,min)
  }

  ##(0,jc),(i,jc)
  idx2 <- which((ww0[,causes2[1]]==cens | ww0[,causes2[1]]==cause[1]) & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- ww0[idx2,timevar2[2]]
  }

  mydata0 <- mydata <- data.frame(time,status,ww0[,covars2],ww0[,indiv2])
  names(mydata) <- c(timevar,causes,covars,indiv2)
  ## }}}

  if (return.data==2) return(list(data=mydata)) else {
  if (!prodlim) {
    ff <- paste("Event(",timevar,",",causes,",cens.code=",cens,") ~ 1",sep="")
    if (!is.null(wname)) covars <- covars[-which(covars %in% c(wname))]
    if (length(c(covars,indiv))>0) {
     xx <- c(covars,indiv2)
     for (i in seq_len(length(xx)))
	     xx[i] <- paste("const(",xx[i],")",sep="")
      ff <- paste(c(ff,xx),collapse="+")
      if (missing(model)) model <- "fg"
    }
    if (missing(model)) model <- "fg"
    ### clusters for iid construction
    lse.clusters <- NULL
    if (!is.null(se.clusters.call)) {
        lse.clusters <- ww0[,"lse.clusters1"]
    }

###    if (!(is.null(wname))) mydata <- ipw2(mydata,time=timevar,cause=causes) #,cens.code=cens)

    if (is.null(wname)) {
	    add<- timereg::comp.risk(as.formula(ff),data=mydata,
	    cause=1,n.sim=0,resample.iid=resample.iid,model=model,conservative=conservative,
	    clusters=lse.clusters, max.clust=max.clust,...)
    } else {
	    add<-timereg::comp.risk(as.formula(ff),data=mydata,
	    cause=1,n.sim=0,resample.iid=resample.iid,model=model,conservative=conservative,
	    clusters=lse.clusters, max.clust=max.clust,
	    weights=mydata[,wname]*mydata$indi.weights,cens.weights=rep(1,nrow(mydata)),...)
    }

    padd <- predict(add,X=1,se=1,uniform=uniform,resample.iid=resample.iid)
    padd$cluster.names <- lse.clusters
  } else {
      if (!requireNamespace("prodlim",quietly=TRUE)) stop("prodlim requested but not installed")
    ff <- as.formula(paste("prodlim::Hist(",timevar,",",causes,")~",paste(c("1",covars,indiv2),collapse="+")))
    padd <- prodlim::prodlim(ff,data=mydata)
  }
###  class(padd) <- c("bicomprisk",class(padd))
 if (return.data==1) return(list(comp.risk=padd,data=mydata)) else return(padd)
  }
## }}}
}

##' @export
bicompriskData <- function(formula, data, cause=c(1,1), cens=0, causes, indiv,
 strata=NULL, id,num, se.clusters=NULL,wname=NULL, ...) {
## {{{

  mycall <- match.call()
  formulaId <- unlist(Specials(formula,"id"))
  formulaIndiv <- Specials(formula,"indiv")
  formulaStrata <- unlist(Specials(formula,"strata"))
  formulaSt <- paste("~.-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")",
                     "-indiv(",paste(formulaIndiv,collapse="+"),")")
  formula <- update(formula,formulaSt)

  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  indiv <- formulaIndiv
  if (!is.null(formulaIndiv)) {
    mycall$indiv <- indiv
  }
  if (missing(id)) stop("Missing 'id' variable")

  ### setting up cluster level for iid decomposition
  if (is.null(se.clusters))  lse.clusters <- data[,c(id)]
  if (length(lse.clusters)!=nrow(data)) stop("'se.clusters' and 'data' does not match!")
  se.clusters.call <- se.clusters


  data <- data.frame(cbind(data,lse.clusters))

  timevar <- terms(formula)[[2]]
  ##  hh <- with(data,eval(timevar))
  if (is.call(timevar)) {
    causes <- timevar[[3]]
    timevar <- timevar[[2]]
  }
  timevar <- as.character(timevar)
  causes <- as.character(causes)

  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  covars <- c(covars,strata)
  ### adds weights,
  if (!is.null(wname)) covars <- c(covars,wname,strata)

  indiv2 <- covars2 <- NULL

  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  ##which(data[,id]%in%names(idtab==2))
  data <- data[which(data[,id]%in%names(idtab==2)),]
  if (missing(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,num] <- unlist(lapply(idtab,seq_len))
  }

  sep <- ""
  timevar2 <- paste(timevar,1:2,sep=sep)
  causes2 <- paste(causes,1:2,sep=sep)
  if (length(covars)>0)
    covars2 <- paste(covars,1,sep=sep)
  for (i in seq_len(length(indiv)))
  indiv2 <- c(indiv2, paste(indiv[i],1:2,sep=sep))

  ww0 <- fast.reshape(data[,c(timevar,causes,covars,indiv,id,num,"lse.clusters")],id=id,num=data$num,labelnum=TRUE)[,c(id,"lse.clusters1",timevar2,causes2,indiv2,covars2)]
  ww0 <- na.omit(ww0)

  status <- rep(0,nrow(ww0))
  time <- ww0[,timevar2[1]]

  ## {{{ (i,j) causes
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cause[2])
  if (length(idx2)>0) {
      status[idx2] <- 1
      time[idx2] <- apply(ww0[idx2,timevar2[1:2],drop=FALSE],1,max)
  }

  ##(0,0), (0,j)
  idx2 <- which(ww0[,causes2[1]]==cens & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 0
      time[idx2] <- ww0[idx2,timevar2[1]]
  }

  ##(ic,0), (ic,j)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- ww0[idx2,timevar2[1]]
  }

  ##(i,0)
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cens)
  if (length(idx2)>0) {
      status[idx2] <- 0
      time[idx2] <- ww0[idx2,timevar2[2]]
  }

  ##(ic,jc)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- apply(ww0[idx2,timevar2[1:2],drop=FALSE],1,min)
  }

  ##(0,jc),(i,jc)
  idx2 <- which((ww0[,causes2[1]]==cens | ww0[,causes2[1]]==cause[1]) & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- ww0[idx2,timevar2[2]]
  }


  mydata0 <- mydata <- data.frame(time,status,ww0[,covars2],ww0[,indiv2],ww0[,id])

  names(mydata) <- c(timevar,causes,covars,indiv2,id)
  ## }}}

  return(mydata)
} ## }}}


## plot.bicomprisk <- function(x,add=FALSE,...) {
##   if ("predict.timereg"%in%class(a)) {
##     if (!add) { plot.predict.timereg(x,...) }
##     else {
##       with(x,lines(time,P1,...))
##     }
##   } else {
##     plot(x,...)
##   }
##   return(invisible(x))
##}
