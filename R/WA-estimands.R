##' While-Alive estimands for recurrent events 
##'
##' Considers the ratio of means \deqn{E(N(min(D,t)))/E(min(D,t))} and the
##' the mean of the events per time unit \deqn{E(N(min(D,t))/min(D,t))} both based on
##' IPCW etimation. RMST estimator equivalent to Kaplan-Meier based estimator.
##'
##' @param formula Event formula first covariate on rhs must be a factor giving the treatment
##' @param data data frame 
##' @param time for estimation 
##' @param cens.code of censorings 
##' @param cause of events 
##' @param death.code of terminal events 
##' @param trans possible power for mean of events per time-unit
##' @param cens.formula censoring model, default is to use strata(treatment)
##' @param augmentR covariates for model of mean ratio
##' @param augmentC covariates for censoring augmentation
##' @param type augmentation for call of binreg, when augmentC is given default is "I" and otherwise "II"
##' @param marks possible marks for composite outcome situation for model for counts with marks
##' @param ...  arguments for binregATE 
##' @author Thomas Scheike
##' @export
WA_recurrent <- function(formula,data,time=NULL,cens.code=0,cause=1,death.code=2,
	 trans=NULL,cens.formula=NULL,augmentR=NULL,augmentC=NULL,type=NULL,marks=NULL,...)
{ ## {{{
  cl <- match.call() ## {{{
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (!inherits(Y,"Event")) stop("Expected a 'Event'-object")
  if (ncol(Y)==2) {
    stop("must give start stop formulation \n"); 
    exit <- Y[,1]
    entry <- NULL ## rep(0,nrow(Y))
    status <- Y[,2]
  } else {
###    stop("only right censored data, will not work for delayed entry\n"); 
    entry <- Y[,1]
    exit <- Y[,2]
    status <- Y[,3]
  }
  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    pos.cluster <- ts$terms
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  } else pos.cluster <- NULL
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms  <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  ### possible handling of id to code from 0:(antid-1)
  if (!is.null(id)) {
      ids <- sort(unique(id))
      nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
     orig.id <- id
   } else { orig.id <- NULL; nid <- nrow(X); 
             id <- as.integer(seq_along(exit))-1; ids <- NULL
  }
  ### id from call coded as numeric 1 -> 

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
  data$id__ <- id ## }}}

  ## use sorted id for all things 
  cid <- countID(data,"id__",sorted=TRUE)
  data$id__ <- cid$indexid
  ### take last record for everybody to use for RMST
  rrR <- subset(data,cid$reverseCountid==1)

## first var on rhs of formula
vars <- all.vars(formula)
treat.name <- vars[4]
treatvar <- data[,treat.name]
if (!is.factor(treatvar)) stop(paste("treatment=",treat.name," must be coded as factor \n",sep="")); 
nlev <- nlevels(treatvar)
nlevs <- levels(treatvar)
ntreatvar <- as.numeric(treatvar)-1
treat.formula <- treat.model <- as.formula(paste(treat.name,"~+1",sep=""))

if (is.null(cens.formula)) cens.formula <- as.formula( paste("~strata(",treat.name,")",collapse=""))
formC <- as.formula( paste("Event(",vars[1],",",vars[2],",",vars[3],"%in% cens.code)~+cluster(id__)",collapse=""))
formD <- as.formula( paste("Event(",vars[2],",",vars[3]," )~-1+",treat.name,"+cluster(id__)",collapse=""))
form1 <- as.formula( paste("Event(",vars[2],",",vars[3],")~-1+",treat.name,"+cluster(id__)",collapse=""))

## take out intercept, to get mean in treated/non-treated
formrec <- update(formula, reformulate(c(".", "-1")))

if (!is.null(augmentR)) {
   varsR <- c(attr(terms(augmentR), "term.labels"))
   form1X <- update(form1, reformulate(c(".", varsR)))
} else form1X <- form1

###print(formD); print(formrec); 
###print(form1); print(form1X); print(cens.formula);

## ratio of means ## {{{
dd <- resmeanIPCW(formD,data=rrR,cause=death.code,cens.code=cens.code,cens.model=cens.formula,time=time,model="l")
ddN <- recregIPCW(formrec,data=data,cause=cause,death.code=death.code,cens.code=cens.code,
		  cens.model=cens.formula,times=time,model="l",marks=marks)
cc <- c(ddN$coef,dd$coef)
cciid <- cbind(ddN$iid,dd$iid)
ratio.means  <- estimate(coef=cc,vcov=crossprod(cciid),f=function(p) (p[1:2]/p[3:4]))
ratio.means.log <- estimate(coef=cc,vcov=crossprod(cciid),f=function(p) log(p[1:2]/p[3:4]))

RAW <- list(iid=cciid,coef=cc,time=time,rmst=dd,meanN=ddN,ratio.means=ratio.means,ratio.means.log=ratio.means.log)
## }}}

data <- count.history(data,id="id__",lag=TRUE,types=cause,status=vars[3],marks=marks,multitype=TRUE)

nameCount <- paste("Count",cause[1],sep="")
formulaCount <- update.formula(formula,.~+1)
cform <- as.formula(paste("~",nameCount,"+cluster(id__)",sep=""))
formulaCount <- update.formula(formulaCount,cform)

## While-Alive mean of events per time-unit 
## with possible marks for death.codes 
if (any(cause %in% death.code)) {
	wd <- match(cause,death.code,nomatch=0)
	mark.codes <- death.code[wd]
} else mark.codes <- NULL
dataDmin <- evalTerminal(formulaCount,data=data,time=time,death.code=death.code, mark.codes=mark.codes,marks=marks)

### setting new response , Ratio of composite outcome
rrR[,"ratio__"] <- dataDmin[cid$reverseCountid==1,"ratio"]
if (!is.null(trans)) {
     rrR[,"ratio__"] <- rrR[,"ratio__"]^trans
}
Yr <- rrR[,"ratio__"]

if (is.null(type)) {
    if (is.null(augmentC)) type <- "II" else type <- "I"
}

outae <- binregATE(form1X,rrR,cause=death.code,time=time,treat.model=treat.formula,
       cens.code=cens.code,Ydirect=Yr,outcome="rmst",model="lin",cens.model=cens.formula,type=type[1],...) 
ET <- list(riskDR=outae)

#ids <- countID(data,"id__",sorted=TRUE)
### assume id is ordered 
data[,"ratio__"] <- outae$Yipcw[cid$indexid+1]

if (!is.null(augmentC)) { ## {{{
dc0 <- dynCensAug(formC,subset(data,ntreatvar==0),augmentC=augmentC,response="ratio__",time=time)
dc1 <- dynCensAug(formC,subset(data,ntreatvar==1),augmentC=augmentC,response="ratio__",time=time)
###nn <- table(rr$treat)
nid <- nrow(rrR)

treats <- function(treatvar) { treatvar <- droplevels(treatvar)# {{{
        nlev <- nlevels(treatvar)
        nlevs <- levels(treatvar)
        ntreatvar <- as.numeric(treatvar)
        return(list(nlev = nlev, nlevs = nlevs, ntreatvar = ntreatvar))
}  

fittreat <- function(treat.model, data, id, ntreatvar, nlev) {
if (nlev == 2) {
    treat.model <- drop.specials(treat.model, "cluster")
    treat <- glm(treat.model, data, family = "binomial")
    iidalpha <- lava::iid(treat, id =id)
    lpa <- treat$linear.predictors
    pal <- expit(lpa)
    pal <- cbind(1 - pal, pal)
    ppp <- (pal/pal[, 1])
    spp <- 1/pal[, 1]
}
else {
    treat.modelid <- update.formula(treat.model, . ~ . + cluster(id__))
    treat <- mlogit(treat.modelid, data)
    iidalpha <- lava::iid(treat)
    pal <- predict(treat, data, se = 0, response = FALSE)
    ppp <- (pal/pal[, 1])
    spp <- 1/pal[, 1]
}
Xtreat <- model.matrix(treat.model, data)
tvg2 <- 1 * (ntreatvar >= 2)
pA <- c(mdi(pal, 1:length(ntreatvar), ntreatvar))
pppy <- c(mdi(ppp, 1:length(ntreatvar), ntreatvar))
Dppy <- (spp * tvg2 - pppy)
Dp <- c()
for (i in seq(nlev - 1)) Dp <- cbind(Dp, Xtreat * ppp[,
    i + 1] * Dppy/spp^2)
DPai <- -1 * Dp/pA^2
out <- list(iidalpha = iidalpha, pA = pA, Dp = Dp, pal = pal,
    ppp = ppp, spp = spp, id = id, DPai = DPai)
return(out)
} # }}}

expit <- function(x) 1/(1 + exp(-x))
idW <- rrR[,"id__"]
treatsvar <- rrR[,treat.name]
treats <- treats(treatsvar)

fitt <- fittreat(treat.model, rrR, idW, treats$ntreatvar, treats$nlev)
iidalpha0 <- fitt$iidalpha
wPA <- c(fitt$pA)

riskDRC <- outae$riskDR+c(dc0$augment,dc1$augment)/nid

MGC0 <- sumstrata(dc0$MGCiid/c(wPA[dc0$id+1]),dc0$id,nid)*dc0$n
MGC1 <- sumstrata(dc1$MGCiid/c(wPA[dc1$id+1]),dc1$id,nid)*dc1$n
MGC <- cbind(MGC0,MGC1)
ccaugment <- apply(MGC,2,sum)

riskDRC <- outae$riskDR+ccaugment/nid
iidDRC <- outae$riskDR.iid+MGC/nid
varDRC <- crossprod(iidDRC)
se.riskDRC <- diag(varDRC)^.5

riskDRC <- list(riskDRC=riskDRC,iid=iidDRC,var=varDRC,coef=riskDRC,se=se.riskDRC)
ET <- list(riskDRC=riskDRC,riskDR=outae)
} ## }}}

out <- list(time=time,id=cid$indexid,orig.id=orig.id,trans=trans,cause=cause,cens.code=cens.code,death.code,
	    RAW=RAW,ET=ET,augmentR=augmentR,augmentC=augmentC)
class(out) <- "WA"
return(out)
} ## }}}

##' @export
print.WA  <- function(x,type="log",...) {# {{{
  print(summary(x),type=type,...)
}# }}}

##' @export
summary.WA <- function(object,type="p",...) {# {{{

rmst <- estimate(object$RAW$rmst)
rmst.test <- estimate(rmst,contrast=rbind(c(1,-1)))
rmst.log <- estimate(rmst,function(p) log(p))
rmst.test.log <- estimate(rmst.log,contrast=rbind(c(1,-1)))

meanNtD <- estimate(object$RAW$meanN)
meanNtD.test <- estimate(meanNtD,contrast=rbind(c(1,-1)))
meanNtD.log <- estimate(meanNtD,function(p) log(p))
meanNtD.test.log <- estimate(meanNtD.log,contrast=rbind(c(1,-1)))

eer <- estimate(object$RAW$ratio.means)
eedr <- estimate(eer,contrast=rbind(c(1,-1)))
ee <- estimate(coef=object$ET$riskDR$riskDR,vcov=object$ET$riskDR$var.riskDR)
eed <- estimate(ee,contrast=rbind(c(1,-1)))
eer.log <- object$RAW$ratio.means.log
eedr.log <- estimate(eer.log,contrast=rbind(c(1,-1)))
eelog <-  estimate(ee,function(p) log(p))
eedlog <- estimate(eelog,contrast=rbind(c(1,-1)))

res <- list(rmst=rmst,rmst.test=rmst.test,meanNtD=meanNtD,meanNtD.test=meanNtD.test,
	    ratio=eer,test.ratio=eedr,meanpt=ee,test.meanpt=eed)
reslog <- list(rmst=rmst.log,rmst.test=rmst.test.log,
	       meanNtD=meanNtD.log,meanNtD.test=meanNtD.test.log,
    ratio=eer.log,test.ratio=eedr.log,meanpt=eelog,test.meanpt=eedlog)
class(res) <- "summary.WA"
class(reslog) <- "summary.WA"
attr(res,"log") <- (type!="p")
attr(reslog,"log") <- (type!="p")

if (type=="p") return(res) else return(reslog)
}# }}}

##' @export
print.summary.WA <- function(x,...)
{# {{{

if (attr(x,"log")) 
cat(paste("While-Alive summaries, log-scale:",x$time,"\n\n"))
else
cat(paste("While-Alive summaries:",x$time,"\n\n"))

cat("RMST,  E(min(D,t)) \n")
print(x$rmst)
cat(" \n")
print(x$rmst.test)

cat("mean events, E(N(min(D,t))): \n")
print(x$meanNtD)
cat(" \n")
print(x$meanNtD.test)
cat("_______________________________________________________ \n")

cat("Ratio of means E(N(min(D,t)))/E(min(D,t)) \n")
print(x$ratio)
cat(" \n")
print(x$test.ratio)
cat("_______________________________________________________ \n")

cat("Mean of Events per time-unit E(N(min(D,t))/min(D,t)) \n")
print(x$meanpt)
cat(" \n")
print(x$test.meanpt)
} # }}}

dynCensAug <- function(formC,data,augmentC=~+1,response="Yipcw",time=NULL,Z=NULL) { ## {{{

if (is.null(time)) stop("must give time of response \n")
   data$Y__ <- data[,response]
   varsC <- c("Y__",attr(terms(augmentC), "term.labels"))
   formCC <- update(formC, reformulate(c(".", varsC)))
###	www <- rep(1,nrow(data))
###	if (treat.specific.cens==1)  www <-data$WW1__;  
   cr2 <- phreg(formCC, data = data, no.opt = TRUE, no.var = 1,Z=Z)
   xx <- cr2$cox.prep
   icoxS0 <- rep(0,length(cr2$S0))
   icoxS0[cr2$S0>0] <- 1/cr2$S0[cr2$S0>0]
   S0i <- rep(0, length(xx$strata))
   S0i[xx$jumps + 1] <- icoxS0
   km <- TRUE
   if (!km) {
        cumhazD <- c(cumsumstratasum(S0i, xx$strata, xx$nstrata)$lagsum)
        St <- exp(-cumhazD)
   } else St <- c(exp(cumsumstratasum(log(1 - S0i), xx$strata, xx$nstrata)$lagsum))

   nterms <- cr2$p-1
   dhessian <- cr2$hessianttime
   dhessian <-  .Call("XXMatFULL",dhessian,cr2$p,PACKAGE="mets")$XXf
   ###  matrix(apply(dhessian,2,sum),3,3)
   timeb <- which(cr2$cumhaz[,1]<time)
   ### take relevant \sum H_i(s,t) (e_i - \bar e)
   covts <- dhessian[timeb,1+1:nterms,drop=FALSE]
   ### construct relevant \sum (e_i - \bar e)^2
   Pt <- dhessian[timeb,-c((1:(nterms+1)),(1:(nterms))*(nterms+1)+1),drop=FALSE]
   ###  matrix(apply(dhessian[,c(5,6,8,9)],2,sum),2,2)
   gammatt <- .Call("CubeVec",Pt,covts,1,PACKAGE="mets")$XXbeta
   S0 <- cr2$S0[timeb]
   gammatt[is.na(gammatt)] <- 0
   gammatt[gammatt==Inf] <- 0
   Gctb <- St[cr2$cox.prep$jumps+1][timeb]
   augmentt.times <- apply(gammatt*cr2$U[timeb,1+1:nterms,drop=FALSE],1,sum)
   augment.times <- sum(augmentt.times)
   if (!is.null(Z)) {
             Zj <- cr2$cox.prep$Z[cr2$cox.prep$jumps+1,][timeb]
             Xaugment.times <- apply( augmentt.times*Zj,2,sum)
   }

   p <- 1
   #### iid magic  for censoring augmentation martingale #
   ### int_0^infty gamma(t) (e_i - ebar(s)) 1/G_c(s) dM_i^c
   xx <- cr2$cox.prep
   nid <- max(xx$id)+1
   jumpsC <- xx$jumps+1
   rr0 <- xx$sign
   S0i <- rep(0,length(xx$strata))
   S0i[jumpsC] <- c(1/(icoxS0*St[jumpsC]))
   S0i[jumpsC] <- icoxS0
   xxtime <- 1*c(xx$time<time)

   pXXA <- ncol(cr2$E)-1
   EA <- cr2$E[timeb,-1,drop=FALSE]
   gammasEs <- .Call("CubeMattime",gammatt,EA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
   gammasE <- matrix(0,length(xx$strata),1)
   gammattt  <-    matrix(0,length(xx$strata),pXXA*1)
   jumpsCt <- jumpsC[timeb]
   gammasE[jumpsCt,] <- gammasEs
   gammattt[jumpsCt,] <- gammatt
   gammaEsdLam0 <- apply(gammasE*S0i*xxtime,2,cumsumstrata,xx$strata,xx$nstrata)
   gammadLam0 <-   apply(gammattt*S0i*xxtime,2,cumsumstrata,xx$strata,xx$nstrata)
   XgammadLam0 <- .Call("CubeMattime",gammadLam0,xx$X[,-1],pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
   Ut <- Et <- matrix(0,length(xx$strata),1)
   Ut[jumpsCt,] <- augmentt.times
   MGCt <- Ut[,drop=FALSE]-(XgammadLam0-gammaEsdLam0)*c(rr0)
   MGCiid <- apply(MGCt,2,sumstrata,xx$id,nid)
   MGCiid <- MGCiid/nid
   #

   nid <- max(cr2$id)
   ids <- headstrata(cr2$id-1,nid)
   ids <- cr2$call.id[ids]

   res <- list(MGCiid=MGCiid,gammat=gammatt,augment=augment.times, id=ids,n=nid)
} ## }}}

##' Evaluates piece constant covariates at min(D,t) where D is a terminal event
##'
##' returns X(min(D,t)) and min(D,t) and their ratio. for censored observation 0. 
##' to use with the IPCW models implemented. 
##'
##' @param formula formula with 'Event' outcome and X to evaluate at min(D,t)
##' @param data data frame
##' @param death.code codes for death (terminating event, 2 default)
##' @param time for evaluation 
##' @param marks for terminal events to add marks*I(D \leq t,\epsilon \in mark.codes)  to X(min(D,t))
##' @param mark.codes gives death codes for which to add mark value
##' @author Thomas Scheike
##' @export
evalTerminal <- function(formula,data=data,death.code=2,time=NULL,marks=NULL,mark.codes=NULL)
{# {{{
    cl <- match.call()# {{{
    m <- match.call(expand.dots = TRUE)[1:3]
    special <- c("strata", "cluster","offset")
    Terms <- terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y,"Event")) stop("Expected a 'Event'-object")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- NULL ## rep(0,nrow(Y))
        status <- Y[,2]
    } else {
        entry <- Y[,1]
        exit <- Y[,2]
        status <- Y[,3]
    }
    id <- strata <- NULL
    if (!is.null(attributes(Terms)$specials$cluster)) {
        ts <- survival::untangle.specials(Terms, "cluster")
        pos.cluster <- ts$terms
        Terms  <- Terms[-ts$terms]
        id <- m[[ts$vars]]
    } else pos.cluster <- NULL
###    if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
###        ts <- survival::untangle.specials(Terms, "strata")
###        pos.strata <- ts$terms
###        Terms  <- Terms[-ts$terms]
###        strata <- m[[ts$vars]]
###        strata.name <- ts$vars
###    }  else { strata.name <- NULL; pos.strata <- NULL}
###    if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
###        ts <- survival::untangle.specials(Terms, "offset")
###        Terms  <- Terms[-ts$terms]
###        offset <- m[[ts$vars]]
###    }
    X <- model.matrix(Terms, m)
    if (!is.null(intpos  <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
    ## }}}

   if (is.null(time)) time <- max(exit)+1

   if (!is.null(id)) {
	   call.id <- id
        ids <- unique(id)
        nid <- length(ids)
        if (is.numeric(id))
            id <-  fast.approx(ids,id)-1
        else  {
            id <- as.integer(factor(id,labels=seq(nid)))-1
        }
    } else { call.id <- id <- as.integer(seq_along(entry))-1;  nid <- nrow(X); }
    ## orginal id coding into integers 1:...
   id <- id+1

   call.marks <- marks
   if (is.null(marks)) marks <- rep(1,length(id))

   dd <- data.frame(id=id)
   dd <- countID(dd,sorted=TRUE)
   ## new id 1,2,.... and so on, referring to rows of data
   id <- dd$indexid+1

 ###
 indexD  <- which(exit <= time & (status %in% death.code))
 indexDM  <- which(exit <= time & (status %in% mark.codes))
 idD <- id[indexD]
 idDM <- id[indexDM]
 obsid <- MarkDid <- Dmintid <- rep(0,nid)
 Dmintid[idD] <- exit[indexD]
 MarkDid[idDM] <- marks[idDM]
 ### 
 indexA<- which(entry <= time & time <= exit)
 idA <- id[indexA]
 Dmintid[idA] <- time
 obsid[idA] <- 1
 Dmint <- Dmintid[id]
 obsid[idA] <- 1
 obsid[idD] <- 1
 ###
 XminDtid <- matrix(0,nid,ncol(X))
 colnames(XminDtid) <- colnames(X)
 XminDtid[idD,] <- X[indexD,]
 XminDtid[idA,]  <- X[indexA,]
 if (length(idDM)>0) XminDtid  <- XminDtid+MarkDid
 XminDt <- XminDtid[id,,drop=FALSE]
 ratio <- rep(0,length(exit))
 obs <- obsid[id]
 ratio[obs==1] <- XminDt[obs==1]/Dmint[obs==1]

 nX <- colnames(XminDt) <- paste(colnames(X),"minDt",sep="") 
 if (ncol(X)==1) nR <- "ratio" else nR <- paste(colnames(X),"ratio",sep="") 
 dd <- data.frame(cbind(XminDt,Dmint,id,call.id,obs,ratio))
 colnames(dd) <- c(nX,"minDt","nid","call.id","uncensored",nR)

 return(dd)
} # }}}

