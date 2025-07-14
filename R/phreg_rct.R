##' Lu-Tsiatis More Efficient Log-Rank for Randomized studies with baseline covariates
##'
##' Efficient implementation of the Lu-Tsiatis improvement using baseline covariates, extended to competing risks and recurrent events. Results
##' almost equivalent with the speffSurv function of the speff2trial function in the survival case. A dynamic 
##' censoring augmentation regression is also computed to gain even more from the censoring augmentation. Furhter, we also deal with twostage
##' randomizations. The function was implemented to deal with recurrent events (start,stop) + cluster, and  more examples in vignette. 
##'
##' @param formula formula with 'Surv' or 'Event' outcome (see \code{coxph}) and
##'   treatment (randomization 0/1)
##' @param data data frame
##' @param cause to use for competing risks, recurrent events data
##' @param cens.code to use for competing risks, recurrent events data
##' @param weights to be used for phreg
##' @param typesR augmentations used for randomization
##' @param typesC augmentations used for censoring
##' @param weights weights for score equation
##' @param augmentR0 formula for the randomization augmentation (~age+sex)
##' @param augmentR1 formula for the randomization augmentation (~age+sex)
##' @param augmentC formula for the censoring augmentation (~age+sex)
##' @param treat.model propensity score model, default is ~+1, assuming an RCT study
##' @param RCT if false will use propensity score adjustment for marginal model
##' @param treat.var in case of twostage randomization, this variable is 1 for
##'   the treatment times, if start,stop then default assumes that only one
##'   treatment at first record
##' @param km use Kaplan-Meier for the censoring weights (stratified on
##'   treatment)
##' @param level of confidence intervals
##' @param cens.model default is censoring model ~strata(treatment) but any
##'   model can be used to make censoring martingales
##' @param estpr estimates propensity scores
##' @param pi0 possible fixed propensity scores for randomizations
##' @param base.augment TRUE to covariate augment baselines (only for R0
##'   augmentation)
##' @param return.augmentR0 to return augmentation data
##' @param mlogit if TRUE then forces use of this function for propensity
##'   scores, default for binary treatment is glm
##' @param ... Additional arguments to phreg function
##' @author Thomas Scheike
##' @references Lu, Tsiatis (2008), Improving the efficiency of the log-rank
##'   test using auxiliary covariates, Biometrika, 679--694
##' 
##' Scheike et al. (2024), WIP, Two-stage randomization for recurrent events, 
##' @examples
##' ## Lu, Tsiatis simulation
##' data <- mets:::simLT(0.7,100)
##' dfactor(data) <- Z.f~Z
##' 
##' out <- phreg_rct(Surv(time,status)~Z.f,data=data,augmentR0=~X,augmentC=~factor(Z):X)
##' summary(out)
##' @export
phreg_rct <- function(formula,data,cause=1,cens.code=0,
     typesR=c("R0","R1","R01"),typesC=c("C","dynC"),weights=NULL,
     augmentR0=NULL,augmentR1=NULL,augmentC=NULL,treat.model=~+1,RCT=TRUE,
     treat.var=NULL,km=TRUE,level=0.95,cens.model=NULL,estpr=1,pi0=0.5,
     base.augment=FALSE,return.augmentR0=FALSE,mlogit=FALSE,...) {# {{{
  Z <- typeII <- NULL
  cl <- match.call()# {{{
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (!inherits(Y,c("Event","Surv"))) stop("Expected a 'Surv' or 'Event'-object")
  if (ncol(Y)==2) {
    exit <- Y[,1]
    entry <- NULL ## rep(0,nrow(Y))
    status <- Y[,2]
    mstatus <- matrix(status,ncol=1)
  } else {
    entry <- Y[,1]
    exit <- Y[,2]
    status <- Y[,3]
    mstatus <- matrix(status,ncol=1)
  }
  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    pos.cluster <- ts$terms
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  } else pos.cluster <- NULL

 call.id <- id;
 conid <- construct_id(id,length(exit))
 name.id <- conid$name.id; id <- conid$id+1; nid <- conid$nid

 data$id__  <-  id
 data$cid__ <- cumsumstrata(rep(1,length(id)),id-1,nid)
 expit <- lava::expit

sides <- function(formula,vars) {# {{{
lhs <- update(formula,.~+1)
rhs <- update(formula,-1~.)
if (all.vars(lhs)[1]==".") 
   formula <- update.formula(formula,as.formula(paste(vars,"~.")))
res <- list(formula=formula,lhs=lhs,rhs=rhs)
}
# }}}

ssform <- sides(formula,"")
varss <- all.vars(ssform$lhs)

## first varaible on rhs of formula
## also candidate for treat variable
streat.name <-  all.vars(ssform$rhs)[1]

treatform <- sides(treat.model,streat.name)
treat.formula <- treatform$formula
treat.name <-  all.vars(treat.formula)[1]
## }}}

treats <- function(treatvar) {# {{{
treatvar <- droplevels(treatvar)
nlev <- nlevels(treatvar)
nlevs <- levels(treatvar)
###treatvar <- as.numeric(treatvar)
ntreatvar <- as.numeric(treatvar)
return(list(nlev=nlev,nlevs=nlevs,ntreatvar=ntreatvar))
}
# }}}

fittreat <- function(treat.model,data,id,ntreatvar,nlev,mlogit=FALSE) {# {{{
if (nlev==2 & !mlogit) {
   treat.model <- drop.specials(treat.model,"cluster")
   treat <- glm(treat.model,data,family="binomial")
   iidalpha <- lava::iid(treat,id=id)
   lpa <- treat$linear.predictors 
   pal <- expit(lpa)
   pal <-cbind(1-pal,pal)
   ppp <- (pal/pal[,1])
   spp <- 1/pal[,1]
} else {  
   treat.modelid <- update.formula(treat.model,.~.+cluster(id__))
   treat <- mlogit(treat.modelid,data)
   iidalpha <- lava::iid(treat)
   pal <- predict(treat,data,se=0,response=FALSE)
   ppp <- (pal/pal[,1])
   spp <- 1/pal[,1]
}

   ###########################################################
   ### computes derivative of D (1/Pa) propensity score 
   ###########################################################
   Xtreat <- model.matrix(treat.model,data)
   tvg2 <- 1*(ntreatvar>=2)
   pA <- c(mdi(pal, 1:length(ntreatvar), ntreatvar))
   pppy <- c(mdi(ppp,1:length(ntreatvar), ntreatvar))
   Dppy <-  (spp*tvg2-pppy) 
   DpA <- c()
   for (i in seq(nlev-1)) DpA <- cbind(DpA,Xtreat*ppp[,i+1]*Dppy/spp^2);  
   DPai <- -1*DpA/pA^2

   ## Dp binomial
   Dp <- Xtreat*(pal[,1]*(1-pal[,1]))

out <- list(iidalpha=iidalpha,pA=pA,Dp=Dp,DpA=Dp,pal=pal,ppp=ppp,spp=spp,id=id,DPai=DPai)
return(out)
} # }}}

if (!is.null(treat.var)) { # {{{
	## time-changing weights
	weightWT <- data[,treat.var]
	whereW <- which(weightWT==1)
	CountW <- cumsumstrata(weightWT,id-1,nid)
        dataW <- data[whereW,]; 
	CountWW <- CountW[whereW]
        idW <- id[whereW]; } 
else { ## only first record is associated with treatment  
	weightWT <- rep(1,nrow(data))
	cid <- cumsumstrata(weightWT,id-1,nid)
	## first record of each subject is treatment 
	whereW <- which(cid==1)
	CountW <- cumsumstrata((cid==1)*1,id-1,nid)
	## constant weights 
	CountWW <- CountW[whereW]
	dataW <- data[whereW,]
        idW <- id[whereW]
	if (!is.null(augmentR1)) warning("must specify treat.var to indicate where treatment is given\nassumes now that treatment is given only at start\n")
} 
# }}}

treatvar <- dataW[,treat.name]
if (!is.factor(treatvar)) stop(paste("treatment=",treat.name," must be factor \n",sep="")); 
treats <- treats(treatvar)

if (estpr[1]==1 ) {
   fitt <- fittreat(treat.formula,dataW,idW,treats$ntreatvar,treats$nlev,mlogit=mlogit)
   pi0 <- fitt$pal[,-1]
   ## p(A) 
   wPA <- c(fitt$pA)
   DPai <- fitt$DPai
} else {
   ## assumes constant fixed prob over groups
   wPA <- ifelse((treats$ntreatvar)[idW]==2,pi0[1],1-pi0[1])         
   pi0 <- rep(pi0[1],length(idW))
}

## construct multiplicative weights, with possible start stop structure
## put propensity score weights at time of weight change, only when RCT=FALSE
ww <-  rep(1,nrow(data))
ww[whereW] <- wPA
wwt <- exp(cumsumstrata(log(ww),id-1,nid))

## set propensity score weights for Cox model's below 
if (!RCT) ww <- 1/wwt  else ww <- rep(1,nrow(data))

if (!is.null(weights)) ww <- weights*ww

rsss <- all.vars(formula)
if (ncol(Y)==2) 
rformulaS <-as.formula( paste("Surv(",rsss[1],",",rsss[2],"==",cause,")~."))
else 
rformulaS <-as.formula( paste("Surv(",rsss[1],",",rsss[2],",",rsss[3],"==",cause,")~."))
formula <- update(formula,rformulaS)

## change id from call to id__
formula <- drop.specials(formula,"cluster")
formula <- update(formula, .~.+cluster(id__))

if (RCT) {
### ... for phreg
fit0 <- phreg(formula,data=data,...)
if (fit0$p>0) eaM <- ea <- ea.iid <- (lava::iid(fit0) %*% fit0$hessian)
else ea <- eaM  <- ea.iid <- matrix(0,max(fit0$id),1)
} else {
fit0 <- phreg_IPTW(formula,data=data,treat.model=treat.formula,treat.var=treat.var,estpr=estpr,pi0=pi0,...)
if (fit0$p>0) {
ea <- ea.iid <- fit0$IID %*% fit0$hessian
## iid without Taylor expansion in weights, to use for censoring augmentation
eaM <- fit0$beta.iid %*% fit0$hessian
} else  ea <- eaM  <- ea.iid <- matrix(0,max(fit0$id),1)
}


data.augR0 <- NULL
AugR0 <- AugR1 <- AugR01 <- rep(0,ncol(ea))
AugR0.iid <-  AugR1.iid <-  AugR01.iid <- matrix(0,nrow(ea),ncol(ea))
if (!is.null(augmentR0)) {# {{{
   ## design 
   ff0 <- sides(augmentR0,all.vars(ssform$rhs)[1])
   dataW0 <- subset(dataW,CountWW==1)
   idW0 <- idW[CountWW==1]
   XR <- model.matrix(augmentR0,dataW0) 
   Z0 <- dataW0[,all.vars(ff0$formula)[1]]
   if (is.factor(Z0)) Z0 <- as.numeric(Z0)-1
   ## order after idW0
   XR[idW0, ] <-  XR 
   Z0[idW0] <- Z0
   strata0 <- piW0 <- rep(0,length(idW0))
   piW0[idW0] <- pi0[CountWW==1]
   ## strata0 0 ??  
   strata0 <- fit0$strata.call[idW0]

   XRpi <- (Z0-piW0)*XR
   XR0pi <- XRpi 
   if (estpr[1]==1) {
       Dp0 <- matrix(0,nid,ncol(fitt$Dp))
       Dp0[idW0,]  <- fitt$Dp[CountWW==1,]
    }
   if (return.augmentR0)  data.augR0 <- list(ea=ea,XRpi=XRpi,XR=XR,A=Z0,p0=piW0)

   augR0fit <- lm(ea~-1+XRpi)
   XRgamma <- augR0fit$fitted.values/(Z0-piW0)
   AugR0.iid <- matrix(augR0fit$fitted.values,ncol=ncol(ea))
   AugR0 <- apply(AugR0.iid,2,sum)


   ## iid term for predicted P(treat=1)
      if (estpr[1]==1) {
	 Dp0f <- crossprod(Dp0,XRgamma) 
         iid.treat <- fitt$iidalpha %*% Dp0f
	 ### same for RCT \hat \pi(,X)
         AugR0.iid <- AugR0.iid - iid.treat
      } 
      ## oucome model iid 
      if (!RCT) {
         DpO <- apply(XRpi,2,sum)
         iid.outcome <- iid(augR0fit)
###         iid.outcome <- iid.outcome %*% DpO
###         AugR0.iid <- AugR0.iid-iid.outcome
      }

} # }}}

if (!is.null(augmentR1)) {# {{{
   ff1 <- sides(augmentR1,all.vars(ssform$rhs)[2])
   dataW1 <- subset(dataW,CountWW==2)
   idW1 <- idW[CountWW==2]
   XR11 <- model.matrix(augmentR1,dataW1) #[,-1,drop=FALSE]
   XR1 <- matrix(0,nid,ncol(XR11))
   XR1[idW1,] <- XR11
   Z1 <- dataW1[,all.vars(ff1$formula)[1]]
   if (is.factor(Z1)) Z1 <- as.numeric(Z1)-1
   piW1 <- pi0[CountWW==2]
   Z1p <- (Z1-piW1)
   Z1p1  <-  rep(0,nid)
   Z1p1[idW1] <- Z1p
   XR1pi <- Z1p1*XR1

   if (estpr[1]==1) {
	 Dp1 <- matrix(0,nid,ncol(fitt$Dp))
	 Dp1[idW1,]  <- fitt$Dp[CountWW==2,]
   }

   augR1fit <- lm(ea~-1+XR1pi)
   XRgamma <- XR1 %*%  coef(augR1fit)
   AugR1.iid <- augR1fit$fitted.values 
   AugR1 <- apply(AugR1.iid,2,sum)

   ## iid term for predicted P(treat=1)
      if (estpr[1]==1) {
	 Dp1f <- crossprod(Dp1,XRgamma) 
         iid.treat <- fitt$iidalpha %*% Dp1f
	 ### same for RCT \hat \pi(,X)
         AugR1.iid <- AugR1.iid - iid.treat
      } 

      ## oucome model iid 
      if (!RCT) {
###         DpO <- apply(XR1pi,2,sum)
###         iid.outcome <- iid(augR1fit)
###         iid.outcome <- iid.outcome %*% rep(DpO,ncol(ea))
###         AugR1.iid <- AugR1.iid-iid.outcome
      }

} # }}}

if (!is.null(augmentR0) & !is.null(augmentR1)) {# {{{
   XRbpi <- cbind(XRpi,XR1pi)
   XRb <- cbind(XR,XR1)
   xxi <- solve(crossprod(XRbpi)) 

   augR01fit <- lm(ea~-1+XRbpi)
   XRgamma <- XRb %*%  coef(augR01fit)
   AugR01.iid <- augR01fit$fitted.values 
   AugR01 <- apply(AugR01.iid,2,sum)

   ## iid term for predicted P(treat=1)
      if (estpr[1]==1) {
         XRgamma0 <- XRpi %*% coef(augR01fit)[1:ncol(XRpi),]
         XRgamma1 <- XR1pi %*% coef(augR01fit)[-(1:ncol(XRpi)),]
         Dp01f <- crossprod(Dp0,XRgamma0)+crossprod(Dp1,XRgamma1)
         iid.treat <- fitt$iidalpha %*% Dp01f
         AugR01.iid <- AugR01.iid - iid.treat
      } 
      ## oucome model iid 
      if (!RCT) {
         DpO <- apply(XRbpi,2,sum)
###         iid.outcome <- iid(augR01fit)
###         iid.outcome <- iid.outcome %*% rep(DpO,ncol(ea))
###         AugR01.iid <- AugR1.iid-iid.outcome
      }

} else {AugR01.iid <- AugR0.iid+AugR1.iid; AugR01 <- AugR0+AugR1;  }
# }}}

varC.improve <- 0; formulaC <- NULL
AugC <- AugC.times <- AugClt <- rep(0,ncol(ea)) 
AugC.iid <- AugClt.iid <- matrix(0,nrow(ea),ncol(ea))

## compute regression augmentation for censoring martingale 
if ((!is.null(augmentC))) {## {{{

xxx <- fit0$cox.prep
### X fra GL + tail-death 
rr <- c(exp(xxx$X %*% coef(fit0)+ xxx$offset)*xxx$weights)
Zrr <- xxx$X*rr
S0i <- rep(0,nrow(xxx$X))
jumps <- xxx$jumps+1
S0i[jumps] <- 1/fit0$S0
E <- U <- matrix(0,nrow(Zrr),ncol(Zrr))
U[jumps,] <- fit0$U
E[jumps,] <- fit0$E
ZEdN <- apply(U,2,revcumsumstrata,xxx$id,nid)
cumhaz <- cumsumstrata(S0i,xxx$strata,xxx$nstrata)
EdLam0 <- apply(E*S0i,2,cumsumstrata,xxx$strata,xxx$nstrata)

mmA <- cbind(status,weightWT,CountW)

## formulac with or without start,stop formulation
if (is.null(cens.model)) 
  cens.model <- as.formula(paste("~strata(",treat.name,")+cluster(id__)"))
else cens.model <- update.formula(cens.model,.~.+cluster(id__))
rrss <- all.vars(ssform$lhs)
if (length(rrss)==2) 
formulaC <-as.formula( paste("Surv(",rrss[1],",",rrss[2],"==",cens.code,")~."))
else 
formulaC <-as.formula( paste("Surv(",rrss[1],",",rrss[2],",",rrss[3],"==",cens.code,")~."))
formulaC <- update.formula(formulaC,cens.model)
###
varsC <- attr(terms(augmentC),"term.labels")
formCC <- update(formulaC, reformulate(c(".", varsC)))

Cfit0 <- phreg(formCC,data=data,no.opt=TRUE,no.var=1,Z=mmA,...)

  ### computing weights for censoring terms
  x <- Cfit0
  xx <- x$cox.prep
  S0i <- rep(0,length(xx$strata))
  jumpsC <- xx$jumps+1
  S0i[jumpsC] <-  1/x$S0
  ## G_c survival at t- 
  if (!km) {
    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
    St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  ###
  rr0 <- c(xx$sign)
  XXA <- xx$X
  EA <- Cfit0$E
  UA <- Cfit0$U
  dhessian <- Cfit0$hessianttime
  hesst <-  .Call("XXMatFULL",dhessian,Cfit0$p,PACKAGE="mets")$XXf

CovZXstrata <- function(X,Ej,Z,Sign,strata,nstrata,jumps) 
{# {{{
		strata  <- c(strata); Sign <- c(Sign)
         ###	Ej <- Ej[jumps,,drop=FALSE]; Ej <- Ej
		ZE <- apply(Z*Sign,2,revcumsumstrata,strata,nstrata)[jumps,,drop=FALSE]; 
		XZ  <- .Call("vecMatMat",X,Z)$vXZ;  
		XZ <- apply(XZ*Sign,2,revcumsumstrata,strata,nstrata)[jumps,,drop=FALSE]; 
		EXZ  <- .Call("vecMatMat",Ej,ZE)$vXZ;  
		out <- XZ-EXZ
		return(out)
}# }}}

fid <- headstrata(xxx$id,nid)

  ## {{{
  if (any(CountWW==2)) {
  rrd <- rr-rr[fid][xxx$id+1]
  Zrrd <- Zrr-Zrr[fid,,drop=FALSE][xxx$id+1,]

  ## place to start with Lam_o(T_R) 
  jumpsW <- which((xx$Z[,2]==1)*(xx$Z[,3]==2)*(xx$sign==-1)==1)
  EdLamTR <- XcumTR <- matrix(0,length(xxx$strata),ncol(xxx$X))
  XcumTR[jumpsW,] <-  Zrrd[jumpsW,]*cumhaz[jumpsW]
  XcumTR <- apply(XcumTR,2,cumsumstrata,xx$id,nid)

  ###  dd <-   cbind(XcumTR,xx$id,xxx$id,xx$status,xx$sign,xxx$sign,xx$Z,xx$time,rr)
  ###  jumpsW; dd[dd[,4]==61,]; dd[dd[,4]==20,]; dd[dd[,4]==95,]; cumhaz[jumpsW]; 

  EdLamTR[jumpsW,] <- EdLam0[jumpsW,]*rrd[jumpsW]
  EdLamTR <- apply(EdLamTR,2,cumsumstrata,xx$id,nid)

  ###  XcumTR <- xx$X*c(cumTR)
  covXTRsZ <-   CovZXstrata(XXA,EA,XcumTR,rr0,xx$strata,xx$nstrata,jumpsC) 
  covELTRsZ <-  CovZXstrata(XXA,EA,EdLamTR,rr0,xx$strata,xx$nstrata,jumpsC) 
  covttt <- covXTRsZ-covELTRsZ
  } else  covttt <- 0
  ## }}}

###  dd <-   cbind(cumTR,xxx$id,xx$sign,xx$Z,xx$time) 

covXsZ <-   CovZXstrata(XXA,EA,Zrr,rr0,xx$strata,xx$nstrata,jumpsC) 
covXsrr <-  CovZXstrata(XXA,EA,as.matrix(rr,ncol=1),rr0,xx$strata,xx$nstrata,jumpsC) 
covXsUs3 <- .Call("vecMatMat",covXsrr,EdLam0[jumpsC,,drop=FALSE])$vXZ;  
covXsUs2 <- covXsZ*cumhaz[jumpsC]-covXsUs3 
### U(infty)= UU
Uinfiid <- -eaM[xx$id+1,,drop=FALSE]
cZEdN <- ZEdN[fid,,drop=FALSE][xx$id+1,,drop=FALSE]-ZEdN
Us1 <- Uinfiid-cZEdN
covXsUs1 <- CovZXstrata(XXA,EA,Us1,rr0,xx$strata,xx$nstrata,jumpsC) 
## scale with Y_(s) because hessiantime is also scaled with this 
covXsYs <- (covXsUs1+covXsUs2-1*covttt)/c(Cfit0$S0); 

p <- ncol(ea)
pXXA <- ncol(XXA)
gammat <-  -.Call("CubeMattime",hesst,covXsYs,pXXA,pXXA,pXXA,p,1,0,0,PACKAGE="mets")$XXX
### solve(matrix(hesst[1,],3,3)) %*% matrix(covXsYs[1,],3,2)
gammat[is.na(gammat)] <- 0
gammat[gammat==Inf] <- 0
augmentt <- .Call("CubeMattime",gammat,UA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
AugCdyn <-  apply(augmentt,2,sum)
gain.times <- .Call("CubeMattime",covXsYs,gammat,pXXA,p,pXXA,p,0,1,0,PACKAGE="mets")$XXX
vardynC.improve  <- matrix(apply(gain.times,2,sum),p,p)

## regress U(s)=\int_s^\infty (Z-E) w(s) dM(s) on agument-model among survivors 
## U(s) = U(\infty) - \int_0^s (Z-E) w(s)  dM(s)
## sum (e_i - \bar e) U(s) Y_i(s)
  ## Lu-Tsiatis augmentation 
  out1 <- iidBaseline(Cfit0,ft=1/St,time=0,fixbeta=0)
  Hiid <- (out1$beta.iid %*% Cfit0$hessian)
  xxi <- solve(crossprod(Hiid))
  ###
  for (i in 1:ncol(ea)) {
     gamma <- xxi %*% crossprod(Hiid, eaM[,i])
     AugClt.iid[,i] <- Hiid %*%  gamma
     AugClt[i] <- sum(AugClt.iid[,i])
  }

Gcj <- St[jumpsC]
varZdN <- matrix(apply(hesst/c(Gcj^2),2,sum),pXXA,pXXA)
covXYdN <- matrix(apply(covXsYs/c(Gcj),2,sum),p,pXXA,byrow=TRUE) 
gamma <- -1*.Call("CubeMattime",matrix(varZdN,nrow=1),matrix(covXYdN,nrow=1),pXXA,pXXA,p,pXXA,1,0,1,PACKAGE="mets")$XXX
gamma <- matrix(gamma,p,pXXA,byrow=TRUE)
gamma[is.na(gamma)] <- 0; gamma[gamma==Inf] <- 0
augment <- c(gamma %*% apply(UA/c(Gcj),2,sum))
var.Clt.improve <-  gamma %*% t(covXYdN) ###  /(nid^2)

#### iid magic  for censoring augmentation martingale{{{
### int_0^infty gamma (e_i - ebar(s)) 1/G_c(s) dM_i^c
S0iG <- S0i <- rep(0,length(xx$strata))
S0iG[jumpsC] <- 1/c(x$S0*Gcj)
S0i[jumpsC] <-  1/x$S0
U <- E <- matrix(0,nrow(xx$X),pXXA)
E[jumpsC,] <- EA; 
U[jumpsC,] <- UA/c(Gcj)
cumhaz <- cumsumstrata(S0iG,xx$strata,xx$nstrata)
EdLam0 <- apply(E*S0iG,2,cumsumstrata,xx$strata,xx$nstrata)

gammasEs <- .Call("CubeMattime",gammat,EA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
gammasE <- matrix(0,nrow(XXA),p)
gammatt  <-    matrix(0,nrow(XXA),pXXA*p)
gammasE[jumpsC,] <- gammasEs
gammatt[jumpsC,] <- gammat
gammaEsdLam0 <- apply(gammasE*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
gammadLam0 <-   apply(gammatt*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
XgammadLam0 <- .Call("CubeMattime",gammadLam0,xx$X,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
Ut <- Et <- matrix(0,nrow(XXA),p)
Ut[jumpsC,] <- augmentt
MGCtt <- Ut[,drop=FALSE]-(XgammadLam0-gammaEsdLam0)*c(rr0)
MGCttiid <- apply(MGCtt,2,sumstrata,xx$id,nid)
AugCdyn.iid <-  MGCttiid ## %*% iH
# }}}

} else {
iid.augment <- iid.augment.times <- augment <- augment.times <- NULL 
var.augment.times <- var.augment <- NULL
var.augment.times.iid <- var.augment.iid <- NULL
Uiid.augment.times <- Uiid.augment <- NULL
time.gammat <- gamma <- gammat <- NULL
ftime.gamma <- NULL
Gcj <- NULL
} ## }}}

## Fitting all models with augmentation terms 
if (fit0$p>0) {
coefMarg <- estimate(fit0,vcov=fit0$var,level=level)$coefmat
var.names <- rownames(coefMarg)
rownames(coefMarg) <- paste("Marginal",rownames(coefMarg),sep="-")
coefs <- coefMarg
} else coefs <- NULL

if (length(typesR)!=3) typesRR <- typesR else 
typesRR <- typesR[c(!is.null(augmentR0),!is.null(augmentR1),!is.null(augmentR0) & !is.null(augmentR1)) ]

if (length(typesC)!=2) typesCC <- typesC else 
typesCC <- typesC[rep(!is.null(augmentC),length(typesC))]
if (is.null(augmentC)) {
	AugCdyn <- AugClt <- AugCdyn.iid <- AugClt.iid <-  0
	typesCC <- "none"
}
if (length(typesRR)==0) { typesRR <- "none" }


baselinecox <- cumhazR0 <- cumhaz <-  se.cumhaz <- se.R0cumhaz <- NULL
iidn <- c()
varbeta <- iid <- fitt <- list(); j <- 0
for (typeR in typesRR) 
for (typeC in typesCC) {# {{{
if (typeR!=typeC) {
   j <- j+1
   AugR <- (typeR=="R0")*AugR0+ (typeR=="R1")*AugR1+(typeR=="R01")*AugR01 +0 
   AugR.iid <- 0+(typeR=="R0")*AugR0.iid+(typeR=="R1")*AugR1.iid + (typeR=="R01")*AugR01.iid
   AugC <- (typeC=="C")*AugClt+(typeC=="dynC")*AugCdyn+0
   Aug <- AugR+AugC

   if (fit0$p>0) iid[[j]] <- -1*(ea.iid-AugR.iid ) %*% fit0$ihessian else iid[[j]] <- 0
   var.beta <- crossprod(iid[[j]])
   AugC.iid <- 0+(typeC=="C")*AugClt.iid+ (typeC=="dynC")*AugCdyn.iid
   if (typeC=="dynC") {
       var.beta <- var.beta - fit0$ihessian %*% vardynC.improve%*% fit0$ihessian
   }
   if (typeC=="C") {
       var.beta <- var.beta - fit0$ihessian %*% var.Clt.improve%*% fit0$ihessian
   }
   if (fit0$p>0) varbeta[[j]] <- var.beta

   if (fit0$p>0) 
   fitts <- phreg(formula,data=data,augmentation=Aug,no.var=1,weights=ww,...)
   else fitts <- fit0

   ## only baseline augment with R0 augmentation
   if (base.augment & typeR=="R0" & (typeC!="dync" | typeC!="C")) { ## {{{
      xxx <- fitts$cox.prep
      xx <- crossprod(XR0pi)
      xxi <- solve(xx)
      XRpit <- XR0pi[xxx$id+1,]
      if (fit0$p>0) rr <- c(exp(xxx$X %*% coef(fitts)+ xxx$offset)*xxx$weights)
      else rr <- c(exp(xxx$offset)*xxx$weights)
      rrs <- rr*c(xxx$sign)
      XRE <- apply(XRpit*rrs,2,revcumsumstrata,xxx$strata,xxx$nstrata)
      S0i2 <- S0i <- rep(0,nrow(xxx$X))
      jumps <- xxx$jumps+1
      S0i[jumps] <- 1/fitts$S0
      ww <- xxx$caseweights*xxx$weights
      S0i2[jumps] <- 1/(fitts$S0^2*ww[xxx$jumps+1])
      ###
      U <- matrix(0,nrow(XRE),ncol(XRE))
      U[jumps,] <- XRpit[jumps,]*S0i[jumps]
      NXRE <- apply(U,2,cumsumstrata,xxx$strata,xxx$nstrata)[jumps,]
      XREdLam0 <- apply(XRE*S0i2,2,cumsumstrata,xxx$strata,xxx$nstrata)[jumps,]
      ns <- c(sumstrata(rep(1,nid),xxx$strata[headstrata(xxx$id,nid)],xxx$nstrata))[xxx$strata[jumps,]+1]
      if (RCT) covBase <- (NXRE-XREdLam0) else covBase <- (NXRE-XREdLam0)

      gamR0Base <- (covBase) %*% xxi
      XRs <- apply(XR0pi,2,sum)
      R0baseline.augment <-  c(XRs %*% t(gamR0Base))
      R0baseline.reduction<- apply(covBase*gamR0Base,1,sum)
      ## using augmented estimator of beta 
      if (fit0$p>0) fitr0 <- robust.phreg(fitts,beta.iid=iid[[j]])
      else fitr0 <- robust.phreg(fit0)
      cumhazR0 <- cumhaz <- fitts$cumhaz
      cumhazR0[,2] <- cumhaz[,2]-R0baseline.augment
      se.R0cumhaz <- se.cumhaz <- fitr0$robse.cumhaz
      se.R0cumhaz[,2] <- (se.cumhaz[,2]^2-R0baseline.reduction)^.5
      baselinecox <- list(phreg=fitts,beta.iid=iid[[j]],XR0pi=XR0pi,gamR0Base=gamR0Base)
   } else baselinecox <- cumhazR0 <- cumhaz <-  se.cumhaz <- se.R0cumhaz <- NULL  ## }}}

   if (fit0$p>0)  {
     coeffitt <- estimate(coef=coef(fitts),vcov=var.beta,level=level)$coefmat
     nnn <- paste(typeR,typeC,sep="_")
     iidn <- c(iidn,nnn)
     rownames(coeffitt) <- paste(nnn,var.names,sep=":")
     coefs <- rbind(coefs,coeffitt)
     if (!is.null(augmentC)) iid[[j]] <- iid[[j]] - AugC.iid%*% fit0$ihessian
   }
   } 
}
names(iid) <- iidn
# }}}

namesortme <- function(iid,name.id) { ## {{{
if (is.matrix(iid))  
	if (nrow(iid)==length(name.id)) {
		rownames(iid) <- name.id
		oid <- order(name.id)
		iid <- iid[oid,]
}
return(iid)
} ## }}}

### sort iid after name.id and put as rownames
ea.iid <- namesortme(ea.iid,name.id)
for (l in 1:length(iid)) iid[[l]] <- namesortme(iid[[l]],name.id)
AugC.iid <- namesortme(AugC.iid,name.id)
AugClt.iid <- namesortme(AugClt.iid,name.id)

###AugR0.iid <- namesortme(AugR0.iid,name.id)

out <- list(marginal=fit0,AugR0=AugR0,AugR1=AugR1,AugR01=AugR01,AugCdyn=AugCdyn,
    AugClt=AugClt,
    coefs=coefs,iid=iid,AugC.iid=AugC.iid,AugClt.iid=AugClt.iid,Cox.iid=ea.iid,
    formula=formula,formulaC=formulaC,treat.model=treat.model,
    id=id,call.id=call.id,name.id=name.id,
    cumhaz=cumhazR0,se.cumhaz=se.R0cumhaz,
    cumhaz.noAug=cumhaz,se.cumhaz.noAug=se.cumhaz, 
    strata=fit0$strata.jumps,nstrata=fit0$nstrata,jumps=seq_along(fit0$strata.jumps),
    strata.name=fit0$strata.name,strata.level=fit0$strata.level,
    baselinecox=baselinecox,data.augR0=data.augR0,var.beta=varbeta)
class(out) <- "phreg_rct"
return(out)
} ## }}} 

##' @export
plot.phreg_rct  <- function(x,...)  { ## {{{
	if (!is.null(x$cumhaz)) baseplot(x,...) 
} ## }}}

##' @export
print.phreg_rct <- function(x,...) {# {{{
  print(summary(x),...)
}# }}}

##' @export
summary.phreg_rct <- function(object,...) {# {{{
if (!is.null(object$coefs)) res  <- object$coefs else res <- "only baseline"
class(res) <- "summary.phreg_rct"
return(res)
}# }}}

simLT <- function(rho,n,beta=0,betac=0,ce=1,betao=0)
{# {{{
Sigma <- matrix(c(1,rho,rho,1),2,2)
M <- t(chol(Sigma))
# M %*% t(M)
Z <- matrix(rnorm(2*n),2,n) # 2 rows, N/2 columns
XY <- t(M %*% Z) 
###
X <- XY[,1]
Y <- XY[,2]
if (betao!=0) px <- expit(X*betao) else px <- 0.5
Z <- rbinom(n,1,px)
TT <- -exp(Z*beta)*log(1-pnorm(Y))
C <- exp(Z*betac)*rexp(n)*ce
status <- (TT<C)
mean(status)
time <- pmin(TT,C)
data <- data.frame(time=time,status=status,X=X,Z=Z)
return(data)
}# }}}

simLTTS <- function(rho,n,beta=c(0,0),betac=0,ce=1,cr=0.5,betao=0.4)
{# {{{

## to avoid R-check
TR <- Z1.f <- count2 <- X1 <- NULL

sigma <- matrix(rho,4,4)
diag(sigma) <- 1
m <- t(chol(sigma))
# m %*% t(m)
z <- matrix(rnorm(4*n),4,n) # 2 rows, n/2 columns
xy <- t(m %*% z) 
zz <- matrix(rnorm(4*n),4,n) # 2 rows, n/2 columns
xyz <- t(m %*% zz) 
###
x <- xy[,1]
y <- xy[,2]
tr <- xy[,3]
x1 <- xyz[,1]
y1 <- xyz[,2]
if (betao!=0) px <- lava::expit(x*betao) else px <- 0.5
if (betao!=0) px1 <- lava::expit(x1*betao+tr) else px1 <- 0.5
z0 <- rbinom(n,1,px)
z1 <- rbinom(n,1,px1)
tt0 <- -exp(z0*beta[1])*log(1-pnorm(y))
tt1 <- -exp(z0*beta[1]+z1*beta[2])*log(1-pnorm(y1))  ## log(1-pnorm(y1))
tr <- exp(z0*beta[1])*rexp(n)*cr 
tt <- ifelse(tt0<tr,tt0,tr)+(tt0>tr)*(tt1)
c <- exp(z0*betac)*rexp(n)*ce
status <- (tt<c)
time <- pmin(tt,c)
data <- data.frame(time=time,status=status,X=x,X1=x1,Z0=z0,Z1=z1,TR=tr)
data <- event.split(data,cuts="TR")
###data <- EventSplit(data,cuts="TR")
data <- dtransform(data,status=2,TR==time)
data <- dtransform(data,Z1=0,start<TR)
data$count2 <- 1*(data$start==data$TR)
data$Count2 <- 1*(data$start==data$TR)
data$cw <- 1
data$Zt.f <- factor(data$Z0)
data$Z0.f <- factor(data$Z0)
data$Z1.f <- factor(data$Z1)
data$Xt <- data$X
data <- dtransform(data,Zt.f=Z1.f,count2==1)
data <- dtransform(data,Xt=X1,count2==1)
data$X0 <- data$X
data$X1t <- 0
data <- dtransform(data,X1t=X1,count2==1)

return(data)
}# }}}

simLUCox <- function(n,cumhaz,death.cumhaz=NULL,cumhaz2=NULL,r1=NULL,r2=NULL,rd=NULL,
		     rc=NULL,rtr=NULL,ztr=0,betatr=0.3,A0=NULL,Z0=NULL,efTR=0,
   gap.time=FALSE,max.recurrent=100,dhaz=NULL,haz2=NULL,dependence=4,var.z=1,cor.mat=NULL,cens=NULL,...) 
{# {{{
  status <- death <- fdeath <- dtime <- NULL # to avoid R-check 

  if (dependence==0) { z <- z1 <- z2 <- zd <- rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- z; z2 <- z; zd <- z 
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
      } else if (dependence==4) {
	      zz <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- zz; z2 <- zz; 
	      zd <- rep(1,n) 
	      z <- z1
      }      else stop("dependence 0-4"); # }}}

      if (!is.null(Z0)) {
	      if (dependence==1) zd <- Z0
	      z <- z1 <- z2 <- Z0
      }

   if (is.null(r1)) r1 <- rep(1,n)
   if (is.null(r2)) r2 <- rep(1,n)
   if (is.null(rd)) rd <- rep(1,n)
   if (is.null(rtr)) rtr <- rep(1,n)
   if (is.null(rc)) rc <- rep(1,n)

  ## extend cumulative for death to full range  of cause 1
  if (!is.null(death.cumhaz)) {
	  if (!is.null(cumhaz2))  {
	     out <- extendCums(list(cumhaz,cumhaz2,death.cumhaz),NULL)
	     cumhaz <- out$cum1
	     cumhaz2 <- out$cum2
	     cumhazd <- out$cum3
	  } else  {
	     out <- extendCums(list(cumhaz,death.cumhaz),NULL)
	     cumhaz <- out$cum1
	     cumhazd <- out$cum2
	     cumhaz2 <- NULL
     }
  }  else cumhaz <- rbind(c(0,0),cumhaz)

  ll <- nrow(cumhaz)
  max.time <- tail(cumhaz[,1],1)

### recurrent first time
  tall1 <- rchaz(cumhaz,r1*z1)
  if (ztr==1)  {
  ## should deal with selection, we know that TR has not happened, when A1(t)=0
  cum1A00 <- cbind(cumhaz[,1],cumsum(c(0,diff(cumhaz[,2]))*(1+var.z*cumhaz2[,2])))
  rtr1 <- exp(betatr)
  cum1A01 <- cbind(cumhaz[,1],cumsum(c(0,diff(cumhaz[,2]))*(1+var.z*cumhaz2[,2]*rtr1)))
  tallA0 <- rchaz(cum1A00,z1[A0==0]*r1[A0==0])
  tallA1 <- rchaz(cum1A01,z1[A0==1]*r1[A0==1])
  tall1[A0==0,] <- tallA0
  tall1[A0==1,] <- tallA1
  } else { cum1A00 <-cum1A01 <-  cumhaz; }
  tall <- tall1
  tall$id <- 1:n
  tall$tr <- 0
  tall$sel <- 1
  tall$timetr <- 0


  ## if tr occurs put tr=1
  if (!is.null(cumhaz2))  {
if (ztr==1) tall2 <- rchaz(cumhaz2,z1*rtr)
if (ztr==0) tall2 <- rchaz(cumhaz2,rtr)
       tall$status <- ifelse(tall1$time<tall2$time,tall1$status,2*tall2$status)
       tall$time <- ifelse(tall1$time<tall2$time,tall1$time,tall2$time)
       tall$tr <- ifelse(tall2$time<tall1$time,1,0)
       tall$timetr <- ifelse(tall2$time<tall1$time,tall2$time,0)
if (ztr==1) {
       ww <- (tall2$time<=tall1$time)
       baseTR <- cpred(cumhaz2,tall2$time[ww])
       basesel <- rtr[ww]*baseTR[,2]
       tall$sel[ww] <- (1+var.z*basesel)/(1+var.z)
  }
  }

### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- rchaz(cumhazd,zd*rd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/(rc*cens)
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }

### fixing the first time to event, if death or censoring occurs before
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  i <- 1
  while (any((tt$time<tt$dtime) & (tt$status!=0) & (i < max.recurrent))) {
	  i <- i+1
	  still <- subset(tt,time<dtime & status!=0)
	  nn <- nrow(still)
	  tr <- still$tr
 	  sel <- still$sel
 	  timetr <- still$timetr
	  r1tr <- rep(1,nn)
###	  if (efTR!=0) r1tr <- exp(efTR*timetr); 
	  rrt <- r1tr*r1[still$id]*(1-tr)+r1tr*r1[still$id]*r2[still$id]*tr
        if (ztr==0) tt1 <- rchaz(cumhaz,rrt*z1[still$id],entry=(1-gap.time)*still$time)
	if (ztr==1)  {

             ## sets up dataframe, overwrites below 
             tt1 <- rchaz(cumhaz,rrt*z1[still$id],entry=(1-gap.time)*still$time)
		if (any(tr==0)) {
			rz <- r1[still$id]*z1[still$id]
			A0s <- A0[still$id]
			if (any(A0s[tr==0]==0)) {
                          ww0 <-  A0s[tr==0]==0
			  tallA0 <- rchaz(cum1A00,rz[ww0],entry=still$time[ww0])
			  tt1[A0s[tr==0]==0,] <- tallA0
			}
			if (any(A0s[tr==0]==1)) {
                          ww1 <-  A0s[tr==0]==1
			  tallA0 <- rchaz(cum1A00,rz[ww1],entry=still$time[ww1])
			  tt1[A0s[tr==0]==1,] <- tallA0
			}
		}
		if (any(tr==1)) {
		    rrt <- r1[still$id]*r2[still$id]*sel[still$id]*z1[still$id]
		    rrt <- rrt[tr==1]
		    tt11 <- rchaz(cumhaz,rrt,entry=(1-gap.time)*still$time[tr==1])
		    tt1[tr==1,] <- tt11
		}
	}
	tt <- tt1
	tt$tr <- tr; tt$sel <- sel; tt$timetr <- timetr; 
	  if (!is.null(cumhaz2))  {
              if (ztr==0)  tt2 <- rchaz(cumhaz2,rtr[still$id],entry=(1-gap.time)*still$time)
              if (ztr==1)  tt2 <- rchaz(cumhaz2,rtr[still$id]*z1[still$id],entry=(1-gap.time)*still$time)
	  ## modify to t2=tr  only when tr==0
          tt$status <- ifelse(tt2$time<=tt1$time & tr==0,2*tt2$status,tt1$status)
          tt$time <-   ifelse(tt2$time<=tt1$time & tr==0,tt2$time,tt1$time)
          tt$timetr <-   ifelse(tt2$time<=tt1$time & tr==0,tt2$time,timetr)
	  if (ztr==1) {
             ww <- (tt2$time<=tt1$time & tr==0)
             Base2TR <- cpred(cumhaz2,tt$time[ww])
             basesel <- r1[still$id][ww]*Base2TR[,2]
             tt$sel[ww] <- (1+var.z*basesel)/(1+var.z)
	  }
          tt$tr[tt2$time<=tt1$time & tr==0] <- 1
	  }
	  if (gap.time) {
		  tt$entry <- still$time
		  tt$time  <- tt$time+still$time
	  }
          ###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tt <- tt[,colnames(tall)]
	  tall <- rbind(tall,tt[1:nn,],row.names=NULL)
  } 
  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  tall$statusD <- tall$status
  tall <- dtransform(tall,statusD=3,death==1)

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"z") <- z1

  return(tall)
  }# }}}

simFrail <- function(n,cumhaz,r1=NULL,rc=NULL,dependence=1,var.z=1,cens=NULL,...) 
{# {{{

 cor.mat <- death <- NULL

  if (dependence==0) { z <- z1 <- z2 <- zd <- rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- z; z2 <- z; zd <- z 
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
      } else if (dependence==4) {
	      zz <- rgamma(n,1/var.z[1])*var.z[1]
	      z1 <- zz; z2 <- zz; zd <- rep(1,n) 
	      z <- z1
      }      else stop("dependence 0-4"); # }}}

   if (is.null(r1)) r1 <- rep(1,n)
   if (is.null(rc)) rc <- rep(1,n)

  ll <- nrow(cumhaz)

 tall   <- rchaz(cumhaz,z1*r1)
  if (!is.null(cens)) { 
     ctime <- rexp(n)/(rc*cens)
     tall$status[tall$time>ctime] <- 0; 
     tall$time[tall$time>ctime] <- ctime[tall$time>ctime] 
  }
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"z") <- z1

  return(tall)
  }# }}}

