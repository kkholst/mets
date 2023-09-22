##' 2 Stage Randomization for Survival Data or competing Risks Data 
##'
##' Under two-stage randomization we can estimate the average treatment effect E(Y(i,j)) of treatment regime (i,j). 
##' The estimator can be agumented in different ways: using the two randomizations and the dynamic censoring augmetatation.
##' The treatment's must be given as factors. 
##'
##' The solved estimating eqution is 
##' \deqn{  (  I(min(T_i,t) < G_i)/G_c(min(T_i ,t)) I(T \leq t, \epsilon=1 ) - AUG_0 - AUG_1 + AUG_C  -  p(i,j)) = 0 }
##' where  using the covariates from augmentR0
##' \deqn{ AUG_0 = \frac{A_0(i) - \pi_0(i)}{ \pi_0(i)} X_0 \gamma_0}
##' and  using the covariates from augmentR1
##' \deqn{ AUG_1 = \frac{A_0(i)}{\pi_0(i)} \frac{A_1(j) - \pi_1(j)}{ \pi_1(j)} X_1 \gamma_1}
##' and   the censoring augmentation is 
##' \deqn{  AUG_C =  \int_0^t \gamma_c(s)^T (e(s) - \bar e(s))  \frac{1}{G_c(s) } dM_c(s) }
##' where 
##' \deqn{ \gamma_c(s)} is chosen to minimize the variance given the dynamic  covariates specified by augmentC.
##'
##' In the observational case, we can use propensity score modelling and outcome modelling (using linear regression).
##' 
##' Standard errors are estimated using the influence function  of all estimators and tests of differences can therefore be computed
##' subsequently.
##'
##' @param formula formula with outcome (see \code{coxph})
##' @param data data frame
##' @param cause cause of interest
##' @param time  time of interest 
##' @param cens.code gives censoring code
##' @param beta starting values 
##' @param treat.model0 logistic treatment model for 1st randomization
##' @param treat.model1 logistic treatment model for 2ndrandomization
##' @param augmentR0 augmentation model for  1st randomization
##' @param augmentR1 augmentation model for  2nd randomization
##' @param augmentC augmentation model for censoring 
##' @param cens.model stratification for censoring model based on observed covariates
##' @param estpr estimate randomization probabilities using model
##' @param response.name can give name of response variable, otherwise reads this as first variable of treat.model1
##' @param response.code code of status of survival data that indicates a response at which 2nd randomization is performed
##' @param offset not implemented 
##' @param weights not implemented 
##' @param cens.weights can be given 
##' @param kaplan.meier  for censoring weights, rather than exp cumulative hazard
##' @param no.opt not implemented 
##' @param method not implemented 
##' @param augmentation not implemented 
##' @param outcome can be c("cif","rmst","rmst-cause")
##' @param model  not implemented, uses linear regression for augmentation
##' @param Ydirect use this Y instead of outcome constructed inside the program (e.g. I(T< t, epsilon=1)), see binreg for more on this
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ddf <- gsim(n,covs=covs,null=null,cens=cens,ce=ce)
##' 
##' bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),ddf$datat,time=2,cause=c(1),
##' 		cens.code=0,treat.model0=A0.f~+1, treat.model1=A1.f~A0.f,
##' 		augmentR1=~X11+X12+TR, augmentR0=~X01+X02,augmentC=~A0+X01+X02+A1+X11+X12+TR,
##' 		response.code=2,...)
##' summary(bb) 
##' @export
binregTSR <- function(formula,data,cause=1,time=NULL,
      cens.code=0,response.code=NULL,
      augmentR0=NULL,treat.model0=~+1, augmentR1=NULL,treat.model1=~+1, 
      augmentC=NULL, cens.model=~+1, estpr=c(1,1),response.name=NULL,
      offset=NULL,weights=NULL,cens.weights=NULL,beta=NULL,
      kaplan.meier=TRUE,no.opt=FALSE,method="nr",augmentation=NULL,
      outcome=c("cif","rmst","rmst-cause"),model="exp",Ydirect=NULL,...)
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
          orig.id <- id
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
  } else { orig.id <- NULL; nid <- nrow(X); id <- as.integer(seq_along(exit))-1; ids <- NULL}
  ### id from call coded as numeric 1 -> 
  id <- id+1
  orig.id <- id.orig <- id;
  nid <- length(unique(id))
###  print(id); print("=TSR===============")

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
# }}}

  if (is.null(time)) stop("Must give time for construction of survival outome modelling \n"); 
  statusC <- (status %in%cens.code) 
  statusE <- (status %in% cause) & (exit<= time) 
  if (is.null(response.code)) stop("Give codes for 2nd treatment start (response), 2nd randomization time \n"); 
  statusR <- (status %in% response.code) & (exit<= time) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 
  if (sum(statusR)==0) stop("No Responses before time \n"); 
  kmt <- kaplan.meier

  response__  <-  rid__  <- Count1__ <- NULL

  data$status__ <-  status 
  data$id__ <-  id
  data$Count1__ <- cumsumstrata(rep(1,length(entry)),id-1,nid)
  data$entry__ <- entry 
  data$exit__ <- exit 
  data$statusC__ <- statusC
  data$status__cause <- statusE
  data$rid__ <- revcumsumstrata(rep(1,length(entry)),id-1,nid)
  data$response__  <-  cumsumstrata(statusR,id-1,nid) 

  cens.strata <- cens.nstrata <- NULL 
  if (is.null(cens.weights))  {
      formC <- update.formula(cens.model,Surv(entry__,exit__,statusC__)~ . +cluster(id__))
      resC <- phreg(formC,data)
      if (resC$p>0) kmt <- FALSE
      exittime <- pmin(exit,time)
      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt)$surv)
      ## strata from original data 
      cens.strata <- resC$strata[order(resC$ord)]
      cens.nstrata <- resC$nstrata
  } else formC <- NULL

  MG.se <- (sum(data$statusC__ & (data$exit__ < time))>0)*1
  expit  <- function(z) 1/(1+exp(-z)) ## expit


  ################################################################
  ### iid data consisting of last record, ordered after id
  ################################################################
  id <- id[data$rid__==1]
  oid <- order(id)
  dataiid <- data[data$rid__==1,][oid,]
###  Xorig <- X <- as.matrix(X)
###  Xdata <- X <- X[data$rid__==1,,drop=FALSE][oid,]
  offset <- offset[data$rid__==1][oid]
  weights <- weights[data$rid__==1][oid]
  status <- status[data$rid__==1][oid]
  exit <- exit[data$rid__==1][oid]
###  X2 <- .Call("vecMatMat", X, X)$vXZ
  ph <- 1
  if (is.null(beta)) beta <- rep(0,ncol(X))
  ## take iid vession of data 
  cens.weights <- cens.weights[data$rid__==1][oid]
  response <- (data$response__[data$rid__==1] >0)[oid]


  if (is.null(beta)) beta <- rep(0,ncol(X))
###  p <- ncol(X)
###  X <-  as.matrix(X)
###  X2  <- .Call("vecCPMat",X)$XX
  ucauses  <-  sort(unique(status))
  ccc <- which(ucauses %in% cens.code)
  if (length(ccc)>=1) Causes <- ucauses[-ccc] else Causes <- ucauses
  obs <- (exit<=time & (status %in% Causes)) | (exit>=time)
  p <- 1

  if (!is.null(Ydirect)) Y <-  Ydirect*obs/cens.weights else {
    if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time)/cens.weights)
    else if (outcome[1]=="rmst") Y <-  c(pmin(exit,time)*obs)/cens.weights 
    else if (outcome[1]=="rmst-cause") Y <- c((status %in% cause)*(time-pmin(exit,time))*obs)/cens.weights
    else stop("outcome not defined") 
  }
  dataiid$Y__ <- Y

 if (is.null(augmentation))  augmentation=rep(0,p)
 nevent <- sum((status %in% cause)*(exit<=time))

treats <- function(treatvar) {# {{{
nlev <- nlevels(treatvar)
nlevs <- levels(treatvar)
###treatvar <- as.numeric(treatvar)
ntreatvar <- as.numeric(treatvar)
return(list(nlev=nlev,nlevs=nlevs,ntreatvar=ntreatvar))
}
# }}}

## fitting treatment models for 1st and 2nd randomization 
fittreat <- function(treat.model,data,id,ntreatvar,nlev) {# {{{
if (nlev==2) {
   treat.model <- drop.specials(treat.model,"cluster")
   treat <- glm(treat.model,data,family="binomial")
   iidalpha <- lava::iid(treat,id=id)
   lpa <- treat$linear.predictors 
   pal <- expit(lpa)
   pal <-cbind(1-pal,pal)
   ppp <- (pal/pal[,1])
   spp <- 1/pal[,1]
} else {  
   treat.modelid <- update.formula(treat.model,.~.+cluster(id))
   treat <- mlogit(treat.modelid,data)
   iidalpha <- lava::iid(treat)
   pal <- predictmlogit(treat,data,se=0,response=FALSE)
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
   Dp <- c()
   for (i in seq(nlev-1)) Dp <- cbind(Dp,Xtreat*ppp[,i+1]*Dppy/spp^2);  
   DPai <- -1*Dp/pA^2

out <- list(iidalpha=iidalpha,pA=pA,lpa=lpa,pal=pal,ppp=ppp,spp=spp,id=id,DPai=DPai)
return(out)
} # }}}

dataR0 <- dataiid
## response==1 data
dataR1 <- subset(dataiid,response__==1)

### treatment is rhs of treat.model 
treat.name0 <-  all.vars(treat.model0)[1]
treatvar0 <- dataR0[,treat.name0]
## make into factor
if (!is.factor(treatvar0))  treatvar0 <- as.factor(treatvar0) 
## treatvar, 1,2,...,nlev or 1,2
treats0 <- treats(treatvar0)

## fit treatment model for first randomization
idR0 <-  id.orig[data$Count1__==1] ## first record with first randomization info
idR0 <- dataiid[,"id__"]
fitt0 <- fittreat(treat.model0,dataR0,idR0,treats0$ntreatvar,treats0$nlev)
iidalpha0 <- fitt0$iidalpha

###
treat.name1 <-  all.vars(treat.model1)[1]
## only  look among those randomized 
treatvar1 <- dataR1[,treat.name1]
###treatvar1 <- dataiid[,treat.name1]
if (!is.factor(treatvar1))  treatvar1 <- factor(treatvar1) 

## fit treatment model for second randomization   using combined model
treats1 <- treats(treatvar1)
###idR1 <- id.orig[data$response__==1 & data$rid__==1] ## first record with first randomization info
idR1 <- dataR1[,"id__"]
dataR1[,treat.name1] <- treatvar1
fit1 <- fittreat(treat.model1,dataR1,idR1,treats1$ntreatvar,treats1$nlev)
## put iid in matrix after nid, and scale to nid
iidalpha1 <- apply(fit1$iidalpha,2,sumstrata,sort(idR1)-1,nid)*(nrow(dataR1)/nid)

## read if there are multiple responses recorded, to make 2nd randomization depend on this 

if (length(all.vars(treat.model1)) >=3 & is.null(response.name)) response.name <- all.vars(treat.model1)[2]
if (!is.null(response.name)) {
respvar <- as.factor(dataR1[,response.name])
if (!is.factor(respvar)) stop(paste("Response=",respvar," must be coded as factor \n",sep="")); 
## treatvar, 1,2,...,nlev or 1,2
## drop possible "0" empty level among dataR1 with responses
respvar <- droplevels(respvar)
nlevresp <- nlevels(respvar)
nlevsr <- levels(respvar)
nrespvar <- as.numeric(respvar)
} else {
nlevresp <- 1 
nlevsr <- "Response"
nrespvar <- rep(1 ,nrow(dataR1))
response.name <- "response"
dataiid[,response.name] <- 1
}

## to fit stratified models 
treat.model1s <- update.formula(treat.model1,~+1)
### stratify on A0.f 
treat.model1s <- update.formula(treat.model1s,as.formula(paste("~",treat.name0)) )
treat.model1s <- treat.model1
###print(treat.model1s)

## fit stratified models for response
tv1 <- list(); fitt1 <- list()
combof <- combo <- NULL
for (nnresp in seq(nlevresp)) {
   k <- nnresp
   dataresp <- subset(dataR1,nrespvar==nnresp)
   tv1[[k]] <- treats(treatvar1[nrespvar==nnresp])
   if (k==1)  {combof <-  tv1[[k]]$nlevs; 
               combo <- seq(tv1[[k]]$nlev) }
   else {
	   combof <- expand.grid(combof,tv1[[k]]$nlevs);
	   combo <- expand.grid(combo,seq(tv1[[k]]$nlev))
   }
###   idresp <- idR1[nrespvar==nnresp]
###   fitt1[[k]] <- fittreat(treat.model1s,dataresp,idresp,tv1[[k]]$ntreatvar,tv1[[k]]$nlev)
}
if (nlevresp==1) { combof <- expand.grid(combof); combo <- expand.grid(combo) }

A0 <- as.numeric(dataiid[,treat.name0])
A1 <- as.numeric(dataiid[,treat.name1])
llR0 <- llR <- llR1 <- llR01 <- c()
riskG0 <- riskG <- riskG1 <- riskG01 <- c()
riskG0C <- riskGC <- riskG1C <- riskG01C <- c()
riskG0.iid <- riskG.iid <- riskG1.iid <- riskG01.iid <- c()
riskG0C.iid <- riskGC.iid <- riskG1C.iid <- riskG01C.iid <- c()
nc <- nrow(combo)
Augment <- Augment.times <-  c()
dynCgammat <- list()

DA1 <- c()

rnames <- c()
cown <- 0
### looping over A0.f response*A1.f
for (i in seq(treats0$nlev)) { ## {{{
for (v in seq(nc))  {
	cown <- cown+1
        namesiv <- paste( paste(treat.name0,"=",treats0$nlevs[i],", ",sep=""),
	      paste(response.name,"*",treat.name1,"=",
		    paste(combof[v,],collapse=","),sep="",collapse=""),sep="")
        rnames <- c(rnames,namesiv)
	## first only one response
	k <- 1
	j <- combo[v,]
	dataij <- dataiid
	pA1j <- rep(1,nrow(dataij))
	## possibly extend 
	if (estpr[2]==1) pA1j[fit1$id] <- fit1$pA else pA1j[fit1$id] <- 0.5
	pA0i <- fitt0$pA
	if (estpr[1]==0) pA0i <- 0.5
	responsetype <- as.numeric(as.factor(dataiid[,response.name]))

	A1j <- (apply(outer(A1,j,FUN="==")*outer(responsetype,seq(nlevresp),FUN="==") ,1,sum) >=1)
	W <- ((A0==i)/pA0i)*((response==0)+(response!=0)*(A1j)/pA1j)
	DW0 <- W*pA0i
	dataij$W__ <- W
	dataij$YW__ <- W*Y
	dataij$W0 <- ((A0==i)-pA0i)/pA0i
	dataij$W1 <- ((A0==i)/pA0i)*(response!=0)*((A1j)-pA1j)/pA1j
	dataij$W11 <- ((A0==i)/pA0i)*((A1j))/pA1j

       ### possibly estimate censoring weights depending on A0.f and A1.f, with weights

	## simple {{{
	ff <- as.formula(YW__ ~ +1)
        ll0 <- lm(ff,data=dataij)
        sll0 <- estimate(ll0)
	## iid ordered after id i data
	iid0 <- ll0$residuals
	h0 <- DW0*Y
        DaPsia <-  apply(fitt0$DPai*h0,2,sum,drop=FALSE)
        iidpala0 <- c(DaPsia %*% t(iidalpha0))

        ddR1 <- response*(A1j)*Y*(A0==i)/pA0i
        ddR1 <- ddR1[response>=1]
        DaPsia <-  apply(fit1$DPai*c(ddR1),2,sum,drop=FALSE)
        iidpala1 <- c(DaPsia %*% t(iidalpha1))
        ###
        risk.iid <-  (iid0+estpr[1]*iidpala0+estpr[2]*iidpala1)/nid
        llR <- c(coef(ll0)[1],crossprod(risk.iid)^.5)
        riskG <- rbind(riskG,llR)
        riskG.iid <- cbind(riskG.iid,risk.iid)
	DA1 <- cbind(DA1,DaPsia)
        # }}}

	### using 1st stage randomization
	if (!is.null(augmentR0)) {# {{{
              varR0 <- all.vars(augmentR0)
	      newnR0 <- paste(c("int",varR0),"__W0",sep="")
	      newaugR0 <- as.formula(paste("~",paste(c("int",varR0),"__W0",sep="",collapse="+")))
	      ff0 <- update.formula(newaugR0,YW__~.)
	      dataij[,newnR0] <- cbind(1,dataij[,varR0]) * dataij$W0
              ###	      
              llR0 <- lm(ff0,data=dataij)
	      iid0 <- llR0$residuals
              ###
	      h0 <- DW0*Y - as.matrix(cbind(1,dataij[,varR0])) %*% coef(llR0)[-1]
              DaPsia <-  apply(fitt0$DPai*(A0==i)*c(h0),2,sum,drop=FALSE)
	      iidpala0 <- c(DaPsia %*% t(fitt0$iidalpha))
              ###
	      ddR1 <- (A1j)*Y*(A0==i)/pA0i
	      ddR1 <- ddR1[response>=1]
	      DaPsia <-  apply(fit1$DPai*c(ddR1),2,sum,drop=FALSE)
	      iidpala1 <- c(DaPsia %*% t(iidalpha1))

	      miid <- estimate(llR0) 
	      miid <-  miid$IC[,-1]/nid
              Dma <-  apply(dataij[,newnR0,drop=FALSE],2,sum,drop=FALSE)
	      miid <-  -c(Dma %*% t(miid))

              risk.iid <- (iid0+estpr[1]*iidpala0+estpr[2]*iidpala1+0*miid)/nid
	      llR0 <- c(coef(llR0)[1],crossprod(risk.iid)^.5)
	      riskG0 <- rbind(riskG0,llR0)
	      riskG0.iid <- cbind(riskG0.iid,risk.iid)
	} # }}}

	### using 2nd stage randomization
	if (!is.null(augmentR1)) {# {{{
              varR1 <- all.vars(augmentR1)
	      newnR1 <- paste(c("int",varR1),"__W1",sep="")
	      newaugR1 <- as.formula(paste("~",paste(c("int",varR1),"__W1",sep="",collapse="+")))
	      ff1 <- update.formula(newaugR1,YW__~.)
	      dataij[,newnR1] <- cbind(1,dataij[,varR1]) * dataij$W1

              llR1 <- lm(ff1,data=dataij)
              sllR1 <- estimate(llR1)
	      iid0 <- llR1$residuals
	      ###   ( W Y -  W1 * h1) 
	      h0 <- (llR1$residuals+coef(llR1)[1])*pA0i
              DaPsia <-  apply(fitt0$DPai*(A0==i)*h0,2,sum,drop=FALSE)
	      iidpala0 <- c(DaPsia %*% t(fitt0$iidalpha))
              ###
	      pred1 <- as.matrix(cbind(1,dataij[,varR1])) %*% coef(llR1)[-1]
	      h1 <- Y - pred1 
	      ddR1 <- (A1j)*h1*(A0==i)/pA0i
	      ddR1 <- ddR1[response>=1]
              DaPsia <-  apply(fit1$DPai*c(ddR1),2,sum,drop=FALSE)
	      iidpala1 <- c(DaPsia %*% t(iidalpha1))

	      miid <-  sllR1$IC[,-1]/nid
              Dma <-  apply(dataij[,newnR1,drop=FALSE],2,sum,drop=FALSE)
	      miid <-  -c(Dma %*% t(miid))

              risk.iid <- (iid0+estpr[1]*iidpala0+estpr[2]*iidpala1+0*miid)/nid
	      llR1 <- c(coef(llR1)[1],crossprod(risk.iid)^.5)
	      riskG1 <- rbind(riskG1,llR1)
	      riskG1.iid <- cbind(riskG1.iid,risk.iid)
	} # }}}

	### using both randomizations
	if ((!is.null(augmentR0)) & (!is.null(augmentR1)) ) {# {{{
	      ff01 <- update(ff0, reformulate(c(".",newnR1)))
              llR01 <- lm(ff01,data=dataij)
	      iid0 <- llR01$residuals
              ###
              varR01 <- all.vars(ff01)[-1]
	      h00 <- as.matrix(cbind(1,dataij[,varR0])) %*% coef(llR01)[newnR0]
	      h01 <- as.matrix(dataij[,newnR1]*pA0i) %*% coef(llR01)[newnR1]
	      h0 <- DW0*Y - h00-h01
	      h11 <- as.matrix(cbind(1,dataij[,varR1])) %*% coef(llR01)[newnR1]
	      h1 <-  Y - h11

              ###
              DaPsia0 <-  apply(fitt0$DPai*(A0==i)*c(h0),2,sum,drop=FALSE)
	      iidpala0 <- c(DaPsia0 %*% t(fitt0$iidalpha))
              ###
	      ddR1 <- (A1j)*h1*(A0==i)/pA0i
	      ddR1 <- ddR1[response>=1]
              DaPsia <-  apply(fit1$DPai*c(ddR1),2,sum,drop=FALSE)
	      iidpala1 <- c(DaPsia %*% t(iidalpha1))

              miid <- estimate(llR01) 
	      miid <-  miid$IC[,-1]/nid
              Dma <-  apply(dataij[,c(newnR0,newnR1),drop=FALSE],2,sum,drop=FALSE)
	      miid <-  -c(Dma %*% t(miid))

              risk.iid <- (iid0+estpr[1]*iidpala0+estpr[2]*iidpala1+0*miid)/nid
	      llR01 <- c(coef(llR01)[1],crossprod(risk.iid)^.5)
	      riskG01 <- rbind(riskG01,llR01)
	      riskG01.iid <- cbind(riskG01.iid,risk.iid)
	} # }}}

 if (MG.se) {## {{{ censoring adjustment of variance one Common KM (could be stratified)
    xx <- resC$cox.prep
    S0i2 <- S0i <- rep(0,length(xx$strata))
    S0i[xx$jumps+1]  <- 1/resC$S0
    S0i2[xx$jumps+1] <- 1/resC$S0^2
    ### set up response assuming that Y from basic data is ordered after id
    Yt <- Y[xx$id+1]
    Wt <- W[xx$id+1]
    Yt <- Yt*Wt*xx$sign
    ## compute function h(s) = \sum_i Y_i(t) I(s \leq T_i \leq t) 
    ## to make \int h(s)/Ys  dM_i^C(s) 
    h  <-  revcumsumstrata(Yt,xx$strata,xx$nstrata)
    ### Cens-Martingale as a function of time and for all subjects to handle strata 
    ## to make \int h(s)/Ys  dM_i^C(s)  = \int h(s)/Ys  dN_i^C(s) - dLambda_i^C(s)
    btime <- 1*(xx$time<time)
    IhdLam0 <- apply(h*S0i2*btime,2,cumsumstrata,xx$strata,xx$nstrata)
    U <- matrix(0,nrow(xx$X),1)
    U[xx$jumps+1,] <- (resC$jumptimes<=time)*h[xx$jumps+1,]/c(resC$S0)
    MGt <- U[,drop=FALSE]*c(xx$weights)-IhdLam0*xx$sign*c(xx$weights)
    ### Censoring Variance Adjustment  
    MGCiid <- apply(MGt,2,sumstrata,xx$id,nid)

    riskG.iid[,cown] <-     riskG.iid[,cown]+MGCiid/nid 
    riskG0.iid[,cown] <-   riskG0.iid[,cown]+MGCiid/nid 
    riskG1.iid[,cown] <-   riskG1.iid[,cown]+MGCiid/nid 
    riskG01.iid[,cown] <- riskG01.iid[,cown]+MGCiid/nid 

  }  else { MGCiid <- 0 }
 ## }}}

      ### computing censoring augmentation
      if ((!is.null(augmentC)) & MG.se) {# {{{

	data$YW__ <- Y[data$id__]*W[data$id__]

        varsC <- c("YW__",attr(terms(augmentC), "term.labels"))
        formCC <- update(formC, reformulate(c(".", varsC)))
        cr2 <- phreg(formCC, data = data, no.opt = TRUE, no.var = 1, ...)
        xx <- cr2$cox.prep
        S0i <- rep(0, length(xx$strata))
        S0i[xx$jumps + 1] <- 1/cr2$S0
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
     augment.times <- sum(augmentt.times)/nid

     Augment.times <- c(Augment.times,augment.times)
     augmentt <- augment.times 

   #### iid magic  for censoring augmentation martingale ## {{{
   ### int_0^infty gamma(t) (e_i - ebar(s)) 1/G_c(s) dM_i^c
   xx <- cr2$cox.prep
   jumpsC <- xx$jumps+1
   rr0 <- xx$sign
   S0i <- rep(0,length(xx$strata))
   S0i[jumpsC] <- c(1/(cr2$S0*St[jumpsC]))
   S0i[jumpsC] <- c(1/(cr2$S0))

   pXXA <- ncol(cr2$E)-1
   EA <- cr2$E[timeb,-1]
   gammasEs <- .Call("CubeMattime",gammatt,EA,pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
   gammasE <- matrix(0,length(xx$strata),1)
   gammattt  <-    matrix(0,length(xx$strata),pXXA*1)
   jumpsCt <- jumpsC[timeb]
   gammasE[jumpsCt,] <- gammasEs
   gammattt[jumpsCt,] <- gammatt
   gammaEsdLam0 <- apply(gammasE*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
   gammadLam0 <-   apply(gammattt*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
   XgammadLam0 <- .Call("CubeMattime",gammadLam0,xx$X[,-1],pXXA,p,pXXA,1,0,1,0,PACKAGE="mets")$XXX
   Ut <- Et <- matrix(0,length(xx$strata),1)
   Ut[jumpsCt,] <- augmentt.times
   MGCt <- Ut[,drop=FALSE]-(XgammadLam0-gammaEsdLam0)*c(rr0)
   MGCiid <- apply(MGCt,2,sumstrata,xx$id,nid)

   riskG[cown,]      <- riskG[cown,]+augmentt
   riskG.iid[,cown]  <- riskG.iid[,cown]+MGCiid/nid 
   riskG0[cown,]     <- riskG0[cown,]+augmentt
   riskG0.iid[,cown] <- riskG0.iid[,cown]+MGCiid/nid 
   riskG1[cown,]     <- riskG1[cown,]+augmentt
   riskG1.iid[,cown] <- riskG1.iid[,cown]+MGCiid/nid 
   riskG01[cown,]    <- riskG01[cown,]+augmentt
   riskG01.iid[,cown]<- riskG01.iid[,cown]+MGCiid/nid 

   dynCgammat  <-  c(dynCgammat,list(gammatt)) 
   ## }}}

    } ## }}}

} } # }}} loop over A0, response*A1.f

varG0 <- varG1 <- varG01 <- NULL
varG <- crossprod(riskG.iid)
if (MG.se) {
   riskG <- cbind(riskG,diag(varG)^.5)[,c(1,3)]
   rownames(riskG) <- rnames
   colnames(riskG) <- c("coef","se","se-fixed-Gc")[1:2]
if ((!is.null(augmentR0))) {
   varG0 <- crossprod(riskG0.iid)
   riskG0 <- cbind(riskG0,diag(varG0)^.5)[,c(1,3)]
   colnames(riskG0) <- c("coef","se","se-fixed-Gc")[1:2]
   rownames(riskG0) <- rnames
} 
if ((!is.null(augmentR1))) {
   varG1 <- crossprod(riskG1.iid)
   riskG1 <- cbind(riskG1,diag(varG1)^.5)[,c(1,3)]
   colnames(riskG1) <- c("coef","se","se-fixed-Gc")[1:2]
   rownames(riskG1) <- rnames
} 
if ((!is.null(augmentR0)) & (!is.null(augmentR0)) ) {
   varG01 <- crossprod(riskG01.iid)
   riskG01 <- cbind(riskG01,diag(varG01)^.5)[,c(1,3)]
   colnames(riskG01) <- c("coef","se","se-fixed-Gc")[1:2]
   rownames(riskG01) <- rnames
}
}

if (!is.null(augmentC) & MG.se) names(Augment.times) <- rnames

###pdiff <- function(x) lava::contr(lapply(seq(x-1), function(z) seq(z,x)))
###contrast <- -pdiff(length(nlevs))
###nncont <- c()
###for (k in seq_along(nlevs[-length(nlevs)])) nncont <-c(nncont, paste("treat:",nlevs[-seq(k)],"-",nlevs[k],sep="")) 
###rownames(contrast) <- nncont
###
###mm <- estimate(coef=val$riskDR,vcov=val$var.riskDR,contrast=contrast)
###val$difriskDR <- mm$coef 
###names(val$difriskDR) <-  rownames(contrast) 
###val$var.difriskDR <- mm$vcov 
###val$se.difriskDR <- diag(val$var.difriskDR)^.5
###

riskG <- list( riskG=riskG, riskG0=riskG0,riskG1=riskG1,riskG01=riskG01)
riskG.iid <- list(riskG0.iid=riskG0.iid,riskG1.iid=riskG1.iid,riskG01.iid=riskG01.iid,
		  riskG.iid=riskG.iid)
varG <- list(varG=varG, varG0=varG0, varG1=varG1, varG01=varG01)
val <- list( riskG.iid=riskG.iid,CensAugment.times=Augment.times,
             dynCens.coef=dynCgammat, riskG=riskG,varG=varG)

  class(val) <- "binregTSR"
  return(val)
}# }}}


##' @export
print.binregTSR  <- function(x,...) {# {{{
	print(summary(x),...)
}# }}}

##' @export
summary.binregTSR <- function(object,...) {# {{{
res <- object$riskG
class(res) <- "summary.binregTSR"
return(res)
}# }}}

##' @export
print.summary.binregTSR  <- function(x,...) {# {{{
    cat("Simple estimator :\n")
    print(x$riskG)
    cat("\n")

    cat("First Randomization Augmentation :\n")
    print(x$riskG0)
    cat("\n")

   cat("Second Randomization Augmentation :\n")
   print(x$riskG1)
   cat("\n")

   cat("1st and 2nd Randomization Augmentation :\n")
   print(x$riskG01)
   cat("\n")
} # }}}

   simdataTSR <- function(n,base12,base1d,base2d,...) { ## {{{
   A0 <- rbinom(n,1,0.5)
   A1 <- rbinom(n,1,0.5)
   X0 <- matrix(rbinom(2*n,1,0.5),n,2)
   expit <- function(x)  return(1/(1+exp(-x)))
   betaX1 <- c(2.7,-2.7)
   X1 <- matrix(rbinom(2*n,1,expit(X0 %*% betaX1)),n,2)
   cov(cbind(X0,X1))
   ###
   beta13 <- c(-0.3,0.4,-0.5)
   beta14 <- c(-0.3,0.2,-0.5)
   beta23 <- c(-0.3,0.4,-0.5)
   beta24 <- c(-0.3,0.2,-0.5)
   betaR <- c(-0.3,0.4,-0.5)
   gamma23 <- 0.3
   gamma24 <- -0.3
   ###
   rr <- exp( cbind(A0,X0) %*% betaR)
   rr13 <- exp( cbind(A0,X0) %*% beta13)
   rr14 <- exp( cbind(A0,X0) %*% beta14)
   rr23 <- exp( cbind(A1,X1) %*% beta23)
   rr24 <- exp( cbind(A1,X1) %*% beta24)
   rd <- cbind(rr13,rr14)
   rd2 <- cbind(rr23,rr24)
   ###
   iddata <- simMultistateII(base12,base1d,base2d,rr=rr,rd=rd,rd2=rd2,
			     early2=1800,gamma23=gamma23,gamma24=gamma24,...)
   covX0 <- cbind(A0[iddata$id],X0[iddata$id,])
   colnames(covX0) <- c("A0","X01","X02")
   covX1 <- cbind(A1[iddata$id],X1[iddata$id,])

   iddata  <-  count.history(iddata,types=2)
   ## only covariates after having gone to stage 2
   covX1 <- covX1*c(iddata$Count2)
   colnames(covX1) <- c("A1","X11","X12")
   iddata <- cbind(iddata,covX0,covX1)
   iddata$R3 <- (iddata$entry<1800)*c(iddata$Count2)
   iddata$R4 <- (iddata$entry<1800)*c(iddata$Count2)
   ###

   cid  <-  countID(iddata)
   iddata$cid <- cid$Countid
   iddata$tid <- cid$Totalid
   return(iddata) 
   } ## }}}

######################## ##################### #####################
## illness death competing risks with two causes of death
## First simple model, then with covariates addedd 
## effect of transition into 2 being early for cause 3, and 4
######################## ##################### #####################
simMultistateII <- function(cumhaz,death.cumhaz,death.cumhaz2,n=NULL,
		    rr=NULL,rd=NULL,rd2=NULL,gamma23=0,gamma24=0,early2=10000,
		    gap.time=FALSE,max.recurrent=100,cens=NULL,rrc=NULL,...) 
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
  cumss <- extendCums(cumss,NULL)
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
  tall <- rcrisks(chaz1,rr1,causes=c(2,3,4))
  tall$id <- 1:n
  tall$status[tall$time>ctime] <- 0; 
  tall$time[tall$time>ctime] <- ctime[tall$time>ctime] 
  tall$from <- 1
  tall$to <- tall$status

  ## simulating out of 2 
  tall2 <- subset(tall,status==2)
  rr23t <- exp(gamma23*(tall2$time<early2))
  rr24t <- exp(gamma24*(tall2$time<early2))
  sim2 <- rcrisks(chaz2,rr2[tall2$id,]*cbind(rr23t,rr24t),causes=c(3,4),entry=tall2$time)
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

   gsim <- function(n,null=1,cens=NULL,ce=9,covs=1,
	    beta0=c(0.1,0.5,-0.5),beta1=c(0.4,0.3,0.5,-0.5),betaR=c(-0.3,-0.5,0.5)) 
   {# {{{
 
   X0 <- matrix(rbinom(2*n,1,0.5),n,2)
   expit <- function(x)  return(1/(1+exp(-x)))
   betaX1 <- 0*c(2.7,-2.7)
   X1 <- matrix(rbinom(2*n,1,expit(X0 %*% betaX1)),n,2)

   O10 <- 0.5; O20 <- 4; OR1 <- 1; OR2 <- 2; 
   if (null==0) Os11 <- 8 
   if (null==1) Os11 <- 3
   Os12 <- 3
   if (null==0) Os21 <- 6
   if (null==1) Os21 <- 3
   Os22 <- 3
   ###
   pr1 <- 0.5
   pr2 <- 0.5
   ###
   p0 <- 0.5
   A0 <- rbinom(n,1,p0)+1
   A1 <- rbinom(n,1,p0)+1
   ###
   if (covs==1) {
   rr01 <- exp( cbind(0,X0) %*% beta0)
   rr02 <- exp( cbind(1,X0) %*% beta0)
   rr11 <- exp( cbind(0,0,X1) %*% beta1 )
   rr12 <- exp( cbind(0,1,X1) %*% beta1 )
   rr21 <- exp( cbind(1,0,X1) %*% beta1 )
   rr22 <- exp( cbind(1,1,X1) %*% beta1 )
   rrr1 <- exp( cbind(0,X0) %*% betaR)
   rrr2 <- exp( cbind(1,X0) %*% betaR)
   } else {
	   rr21 <- rr22  <- rr01 <- rr02 <- rr21 <- rr22 <- rr11 <- rr12 <- rrr1 <- rrr2 <- 1
   }
   ###
   TR1 <- runif(n)*OR1*rrr1
   TR2 <- runif(n)*OR2*rrr2
   T10 <- rexp(n)*O10*rr01
   T20 <- rexp(n)*O20*rr02
   ### response if (TR1 < T10) or (TR2 < T20)
   R1 <- (TR1 < T10)*1 
   R2 <- (TR2 < T20)*1 
   ###
   Ts11 <- rexp(n)*Os11*exp(-3*TR1)*rr11
   Ts12 <- rexp(n)*Os12*exp(-3*TR1)*rr12
   Ts21 <- rexp(n)*Os21*exp(-1*TR2)*rr21
   Ts22 <- rexp(n)*Os22*exp(-1*TR2)*rr22
   ###
   Tt11 <- TR1 + Ts11
   Tt12 <- TR1 + Ts12
   Tt21 <- TR2 + Ts21
   Tt22 <- TR2 + Ts22
   ###
   TT <- (A0==1)*( (1-R1)*T10 + R1*((A1==1)*Tt11 + (A1==2)*Tt12))+ 
         (A0==2)*( (1-R2)*T20 + R2*((A1==1)*Tt21 + (A1==2)*Tt22))  
###	     
 TT11  <-  (1-R1)*T10 + R1*Tt11
 TT12  <-  (1-R1)*T10 + R1*Tt12
 TT21  <-  (1-R2)*T20 + R2*Tt21
 TT22  <-  (1-R2)*T20 + R2*Tt22
 TTt <- cbind(TT11,TT12,TT21,TT22)
 if (is.null(cens)) cc <- max(TT)+1 else cc <- rexp(n)*ce
  ###
  time <- pmin(TT,cc)
  status <- (TT<cc)*1
  ###  
  data <- data.frame(time=time,T10=T10,Tt11=Tt11,Tt12=Tt12,TT=TT,status=status,A1=A1,A0=A0,R1=R1,R2=R2,R=(A0==1)*R1+(A0==2)*R2,id=1:n,TR=(A0==1)*TR1+(A0==2)*TR2)
  dfactor(data) <- A0.f~A0
  dfactor(data) <- A1.f~A1
  datat <- event.split(data,cuts="TR",name.start="entry")

  datat <- dtransform(datat,status=2,time==TR)
  datat  <-  count.history(datat,types=2)
  datat  <-  dtransform(datat,A1=0,Count2==0)
  datat$response <- (datat$Count2==1)*1

  covX0 <- X0[datat$id,]
  colnames(covX0) <- c("X01","X02")
  covX1 <- X1[datat$id,]
  ###
  covX1 <- covX1*c(datat$Count2)
  colnames(covX1) <- c("X11","X12")
  datat <- cbind(datat,covX0,covX1)
  datat <- dtransform(datat,TR=0,Count2==0) 

  res <- list(data=data,datat=datat,TTt=TTt)
  return(res)
  } # }}}


