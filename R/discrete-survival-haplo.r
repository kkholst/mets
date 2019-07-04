##' Discrete time to event haplo type analysis 
##'
##' Can be used for logistic regression when time variable is "1" for all id. 
##'
##' Cycle-specific logistic regression of haplo-type effects with known 
##' haplo-type probabilities. Given observed genotype G and unobserved haplotypes H
##' we here mix out over the possible haplotypes using that P(H|G) is provided. 
##' 
##' \deqn{
##' S(t|x,G)) = E( S(t|x,H) | G)  = \sum_{h \in G} P(h|G) S(t|z,h) 
##' }
##' so survival can be computed by mixing out over possible h given g.
##'
##' Survival is based on logistic regression for the discrete hazard function of the
##' form 
##' \deqn{
##' logit(P(T=t| T \geq t, x,h)) = \alpha_t + x(h) \beta
##' }
##' where x(h) is a regression design of x and haplotypes \eqn{h=(h_1,h_2)}
##' 
##' Likelihood is maximized and standard errors assumes that P(H|G) is known. 
##' 
##' The design over the possible haplotypes is constructed by merging X with Haplos and  
##' can be viewed by design.only=TRUE
##' 
##' @param X design matrix data-frame (sorted after id and time variable) with id time response  and desnames 
##' @param y name of response (binary response with logistic link) from X
##' @param time.name to sort after time  for X
##' @param Haplos (data.frame with id, haplo1, haplo2 (haplotypes (h)) and  p=P(h|G)) haplotypes given as factor.  
##' @param id name of id variale from X
##' @param desnames names for design matrix
##' @param designfunc function that computes design given haplotypes h=(h1,h2) x(h) 
##' @param beta starting values 
##' @param no.opt optimization TRUE/FALSE 
##' @param method NR, nlm 
##' @param stderr to return only estimate 
##' @param designMatrix  designMatrix not implemented 
##' @param design.only to return only design matrices for haplo-type analyses. 
##' @param covnames names of covariates to extract from object for regression
##' @param fam family of models, now binomial default and only option 
##' @param ... Additional arguments to lower level funtions lava:::NR  optimizer or nlm
##' @author Thomas Scheike
##' @examples
##' ## some haplotypes of interest
##' types <- c("DCGCGCTCACG","DTCCGCTGACG","ITCAGTTGACG","ITCCGCTGAGG")
##' 
##' ## some haplotypes frequencies for simulations 
##' data(hapfreqs)
##'
##' www <-which(hapfreqs$haplotype %in% types)
##' hapfreqs$freq[www]
##'
##' baseline=hapfreqs$haplotype[9]
##' baseline
##'
##' designftypes <- function(x,sm=0) {# {{{
##' hap1=x[1]
##' hap2=x[2]
##' if (sm==0) y <- 1*( (hap1==types) | (hap2==types))
##' if (sm==1) y <- 1*(hap1==types) + 1*(hap2==types)
##' return(y)
##' }# }}}
##'
##' tcoef=c(-1.93110204,-0.47531630,-0.04118204,-1.57872602,-0.22176426,-0.13836416,
##' 0.88830288,0.60756224,0.39802821,0.32706859)
##' 
##' data(hHaplos)
##' data(haploX)
##' 
##' haploX$time <- haploX$times
##' Xdes <- model.matrix(~factor(time),haploX)
##' colnames(Xdes) <- paste("X",1:ncol(Xdes),sep="")
##' X <- dkeep(haploX,~id+y+time)
##' X <- cbind(X,Xdes)
##' Haplos <- dkeep(ghaplos,~id+"haplo*"+p)
##' desnames=paste("X",1:6,sep="")   # six X's related to 6 cycles 
##' out <- haplo.surv.discrete(X=X,y="y",time.name="time",
##'          Haplos=Haplos,desnames=desnames,designfunc=designftypes) 
##' names(out$coef) <- c(desnames,types)
##' out$coef
##' 
##' @aliases simGlm  predictSurvd
##' @export
haplo.surv.discrete <- function (
###	formula = formula(data),data=sys.parent(),
    X=NULL,y="y",time.name="time",Haplos=NULL,id="id", desnames=NULL, designfunc=NULL,beta=NULL, 
    no.opt=FALSE,method="NR",stderr=TRUE, designMatrix=NULL,design.only=FALSE,
###    geno.type = NULL, geno.setup=NULL,haplo.freq = NULL,
###    haplo.design=NULL,haplo.baseline=NULL,alpha=NULL, 
    covnames=NULL,fam=binomial,...)
{ ## {{{ 

###	browser()

  ## X: y, Xdes, id
  ## Haplos, haplo1, haplo2, id, p (HgivenG)
  bothid <- intersect(X$id,Haplos$id)
  X <- subset(X,id %in% bothid)
  Haplos <- subset(Haplos,id %in% bothid)
  ## new iid starts at 1 
  Haplos$id <- fast.approx(bothid,Haplos$id) 
  iid <- Haplos$id-1
  X$id <- fast.approx(bothid,X$id) 

  Xhap <- merge(X,Haplos,by.x="id",by.y="id")
  Xhap <- dsort2(Xhap,~id+"haplo*"+time)
  Haplos <- dsort2(Haplos,~id+"haplo*")

  time <- Xhap[,time.name]
  tn <- match(time.name,names(Xhap))
  Xhap <- Xhap[,-tn]

  response <- Xhap[,y]
  yn <- match(y,names(Xhap))
  Xhap <- Xhap[,-yn]

  mm <-  grep("haplo*",names(Xhap))
  Xhaps <- Xhap[,mm]
  X <- as.matrix(apply(Xhap[,mm],1,designfunc))
  ###
  if (ncol(X)==nrow(Xhap))  X <- t(X)
  colnames(X) <- paste("haplo",1:ncol(X),sep="")
  X <- as.matrix(cbind(Xhap[,desnames],X))

  ## X(i)^T %*% X(i) for each row 
  X2  <- .Call("vecMatMat",X,X)$vXZ

## creates sub-index for haplo-types within each id
   nmm <- names(Xhap)[mm]
   lll <- lapply(Xhap[,c("id",nmm)],as.numeric)
   stratidhap <-    as.numeric(survival::strata(lll))
###   stratidhap <- fast.approx(unique(stratidhap),stratidhap) 
###   print(cbind(Xhap$id,stratidhap))

   nidhap <- length(unique(stratidhap))
   nid <- length(unique(Haplos$id))

   hgiveng <- Haplos$p

if (!design.only) {

if (is.null(beta)) beta <- rep(0,ncol(X))

expit  <- function(z) 1/(1+exp(-z)) ## expit


obj <- function(pp,all=FALSE)
{ # {{{

lp <- X %*% pp
## plp <- family$linkinv(lp)
plp <- expit(lp)
nplp <- 1-plp
###lognp <- log(nplp)

### logpht <-  (response - plp)/family$variance(plp)
logpht <- log(plp)*response+log(nplp)*(1-response)
pht <-c(exp(logpht))
Dlogpht <-  X* c(response-plp)
D2logpht <- c(plp/(1+exp(lp)))*X2

ph <- c(exp(sumstrata(logpht,stratidhap-1,nidhap)))
pg <- c(sumstrata(ph*hgiveng,iid,nid))
logl <- log(pg)

## Derivative 
Dlogph  <- apply(Dlogpht,2,sumstrata,stratidhap-1,nidhap)
Dph  <- c(ph)*Dlogph
Dpg <- apply(Dph*hgiveng,2,sumstrata,iid,nid)
Dlogl    <- Dpg/pg
DpgDpg  <- .Call("vecMatMat",Dpg,Dpg)$vXZ

## 2nd Derivative 
D2logph <- apply(D2logpht,2,sumstrata,stratidhap-1,nidhap)
DphDlogph <- .Call("vecMatMat",Dph,Dlogph)$vXZ
D2ph <- ph*D2logph+DphDlogph
D2pg <-apply(D2ph*hgiveng,2,sumstrata,iid,nid)
D2logi <- (pg*D2pg-DpgDpg)/pg^2
D2log <- apply(D2logi,2,sum)
D2log <- matrix(D2log,ncol(X),ncol(X))


ploglik <- sum(logl)
gradient <- apply(Dlogl,2,sum)
hessian <- D2log

  if (all) {
      ihess <- solve(hessian)
      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
      robvar <- crossprod(beta.iid)
      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
		  iid=beta.iid,robvar=robvar,var=ihess,
		  se=diag(ihess)^.5,se.robust=diag(robvar)^.5)
      return(val)
  }  
 structure(-ploglik,gradient=-gradient,hessian=hessian)
}# }}}

  p <- ncol(X)
  opt <- NULL
  if (p>0) {
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
          tim <- system.time(opt <- lava::NR(beta,obj,...))
          opt$timing <- tim
          opt$estimate <- opt$par
      } else {
          opt <- nlm(obj,beta,...)
          opt$method <- "nlm"
      }
      cc <- opt$estimate; 
      if (!stderr) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }


###
###  ## {{{ defining names with haplostuff
###  ## best guess  for standard designs  
###  haploX.names<-c()
###  if (px!=dimxih) {
###    haploef.x<-dimxih-px; 
###    if (haploef.x>0) haploX.names<-rep("Haplo effect",haploef.x)
###  }
###  haploZ.names<-c()
###  if (pz!=dimzih) {
###    haploef.z<-dimzih-pz; 
###    if (haploef.z>0) haploZ.names<-rep("Haplo effect",haploef.z)
###  }
###  ## }}}
###
###
###  ## {{{ adding names 
###  if (is.null(covnames)==TRUE)  covnamesXuse<-c(covnamesXin,haploX.names)  
###  else covnamesXuse<-c(covnamesX)
###
###  if (is.null(covnamesZ)==TRUE)  covnamesZuse<-c(covnamesZin,haploZ.names)  
###  else covnamesZuse<-c(covnamesZ)
####print(covnamesZuse); print(covnamesXuse)
###
###  if (length(covnamesXuse)==dimxih) {
###     colnames(ud$cum)<-colnames(ud$var.cum) <-
###     colnames(ud$robvar.cum) <-c("time",covnamesXuse); 
###  }
###
###  if ((length(covnamesZuse)==dimzih))  {
###      names(ud$score)<-covnamesZuse;
###      if (fix.beta==0) colnames(ud$test.procProp)<-c("time",covnamesZuse)
###      if (resample.iid==1) colnames(gamiid)<-covnamesZuse; 
###      if ((sim>=1) &  (fix.beta==0)) names(ud$pval.Prop)<-covnamesZuse
###  }
###  if ((sim>=1) & (length(covnamesXuse)==dimxih)) {
###    names(ud$conf.band)<- names(ud$pval.testBeq0)<-
###    names(ud$pval.testBeqC)<- names(ud$obs.testBeq0)<- 
###    names(ud$obs.testBeqC)<- colnames(ud$sim.testBeq0)<-covnamesXuse 
###  } 
###  if (sim==0) {
###    ud$pval.Prop<- ud$conf.band<- ud$pval.testBeq0<- ud$pval.testBeqC<-
###      ud$obs.testBeq0<- ud$obs.testBeqC<- ud$sim.testBeq0<- NULL 
###  }
###
###  if ((length(covnamesZuse)==dimzih)) {
###     rownames(ud$gamma)<-covnamesZuse;
###     colnames(ud$gamma)<-"estimate"; 
###
###     if (fix.beta==1) {
###        namematrix(ud$var.gamma,covnamesZuse); 
###        namematrix(ud$robvar.gamma,covnamesZuse); 
###        namematrix(ud$D2linv,covnamesZuse); 
###     }
###  } 
###  if (fix.beta==1) {  ud$var.gamma<-matrix(0,pz,pz); 
###                      ud$robvar.gamma<-matrix(0,pz,pz);
###                    }
###  ## }}}
###
###  attr(ud, "Call") <- sys.call()
###  attr(ud, "Formula") <- formula
###  attr(ud, "id") <- id.call
}  else val <- NULL

  val <- c(list(Xhap=Xhap,X=X,Haplos=Haplos),val)
  class(val) <- "survd"
  return(val)
} ## }}} 

##' @export
simGlm <- function(coef=NULL,n=100,Xglm=NULL,times=NULL)
{# {{{
	  
  Z <- Xglm  
  if (!is.null(Z)) n <- nrow(Z) 

  if (!is.null(Z)) data <- Z else data <- data.frame(id=1:n)

  if (!is.null(times)) {
     timesf <- data.frame(times=rep(times,n),id=rep(1:n,each=length(times)))
     data <- merge(data,timesf,by.x="id",by.y="id")
     mt <- model.matrix(~factor(times),data)
     nm <- match(c("id","times"),names(data))
     Z <- cbind(mt,data[,-nm])
  }

  expit  <- function(z) 1/(1+exp(-z)) ## expit

  p <- c(expit(as.matrix(Z) %*% coef))
  y <- rbinom(length(p),1,p)

  data <- cbind(y,data)
  data <- count.history(data,status="y",id="id",types=1)
  data <- subset(data,data$Count1<=0)

  attr(data,"coef") <- beta
  return(data)
 }# }}}

##' @export
predictSurvd <- function(hsd,Z,times=1:6,se=FALSE,type="prob")
{# {{{
  if (!is.null(Z)) n <- nrow(Z) 
  if (!is.null(Z)) data <- Z else data <- data.frame(id=1:n)
  Z <- data.frame(Z)
  Z$id <- 1:n

  if (!se) {# {{{{{{
	  data <- Z
	  if (!is.null(times)) {
	     timesf <- data.frame(times=rep(times,n),id=rep(1:n,each=length(times)))
	     data <- merge(data,timesf,by.x="id",by.y="id")
	     mt <- model.matrix(~factor(times),data)
	     nm <- match(c("id","times"),names(data))
	     Z <- cbind(mt,data[,-nm])
	  }
	  ccc <- hsd$coef
	  if (ncol(Z)!=length(c(ccc))) {
		  print(head(Z))
		  print(ccc)
		  stop("dimension of Z not consistent with length of coefficients"); 
	  }

	  p <- c(expit(as.matrix(Z) %*% ccc))

	  pred <- data.frame(p=p,id=data$id,times=data$times)
	  survt <- 1-exp(cumsumstrata(log(1-pred$p),data$id-1,6))
	  if (type=="prob") survt <- 1-survt
	  if (type=="hazard") survt <- p
	  pred <- cbind(pred,survt)
# }}}
  } else {# {{{
    ccc <- hsd$coef

    expit <- function(p) exp(p)/(1+exp(p))
    Ft <- function(p)
    {
	   xp <- as.matrix(Zi) %*% p
	   lam <- expit(xp)
	   st <- cumprod(1-lam)
	   if (type=="prob") st <- 1-st 
	   if (type=="hazard") st <- lam
	   return(st)
    }

  preds <- c()
  for (i in 1:nrow(Z)) {
     Zi <- data.frame(Z[i,,drop=FALSE])
     data <- Zi
     if (!is.null(times)) {
        timesf <- data.frame(times=rep(times,n),id=rep(1:n,each=length(times)))
        data <- merge(data,timesf,by.x="id",by.y="id")
        mt <- model.matrix(~factor(times),data)
        nm <- match(c("id","times"),names(data))
        Zi <- cbind(mt,data[,-nm])
     }
     eud <- estimate(coef=hsd$coef,vcov=hsd$var,f=function(p) Ft(p))
     cmat <- data.frame(eud$coefmat)
     cmat$id <- i
     cmat$times <- times
     names(cmat)[1:4] <- c("pred","se","lower","upper")
     preds <- rbind(preds,cmat)
  } 

  }# }}}

return(preds)
}# }}}


##' ## uses HaploSurvival package of github install via devtools
##' ## devtools::install_github("scheike/HaploSurvival")
##' ## this is only used for simulations 
##' ## out <- simHaplo(1,100,tcoef,hapfreqs)

###simHaplo <- function(i,n,tcoef,hapfreqs)
###{ ## {{{  
###
###   haplos <- sample(19,2*n,replace=TRUE,prob=hapfreqs$freq)
###   haplos <- matrix(haplos,n,2) 
###   hap1 <- hapfreqs$haplotype[haplos[,1]]
###   hap2 <- hapfreqs$haplotype[haplos[,2]]
###   ###
###   X <- t(apply(cbind(hap1,hap2),1,designftypes))
###   X <- data.frame(X,id=1:n)
###   sud <- simGlm(coef=tcoef,Xglm=X,n=n,times=1:6)
###   ## known haplotypes
###   ssud <- glm(y~factor(times)+X1+X2+X3+X4+X5,data=sud,family=binomial())
###
###   genotype <- c()
###   for (i in 1:11)
###     genotype <- cbind(genotype,substr(hap1,i,i), substr(hap2,i,i) )
###
###   setup <- geno.setup(genotype,haplo.baseline=baseline,sep="") 
###   wh <- match(setup$uniqueHaploNames,hapfreqs$haplotype)
###   wwf <- wh[!is.na(wh)]
###   ghaplos <- matrix(unlist(setup$HPIordered),byrow=TRUE,ncol=2)
###   ghaplos <- cbind(rep(1:n,setup$nPossHaps),ghaplos)
###   ghaplos <- data.frame(ghaplos)
###   names(ghaplos) <- c("id","haplo1","haplo2")
###   haploff <- rep(0,length(setup$uniqueHaploNames))
###   haploff[!is.na(wh)] <- hapfreqs[wwf,"freq"]
######
###   hap1f <- haploff[ghaplos[,2]]
###   hap2f <- haploff[ghaplos[,3]]
###   hap12f <- hap1f*hap2f
###   hapsshed  <- hap12f
###   ghaplos$p <- hapsshed
###   ghaplos <- subset(ghaplos,p>0)
###   ## back to characters for indentification of design
###   ghaplos$haplo1 <- as.factor(setup$uniqueHaploNames[ghaplos$haplo1])
###   ghaplos$haplo2 <- as.factor(setup$uniqueHaploNames[ghaplos$haplo2])
###   ptot <- sumstrata(ghaplos$p,ghaplos$id-1,n)
###   ghaplos$p <- ghaplos$p/ptot[ghaplos$id]
###
###sud$time <- sud$times
###Xdes <- model.matrix(~factor(time),sud)
###colnames(Xdes) <- paste("X",1:ncol(Xdes),sep="")
###X <- dkeep(sud,~id+y+time)
###X <- cbind(X,Xdes)
###Haplos <- dkeep(ghaplos,~id+"haplo*"+p)
###dtable(Haplos,~"haplo*",level=1)
######Haplos
######
###y <- "y"
###time.name="time"
###desnames=paste("X",1:6,sep="")
######
######
###mm <- system.time(
###mud <- haplo.surv.discrete(X=X,y="y",time.name="time",##design.only=TRUE,
###		  Haplos=Haplos,designfunc=designftypes,desnames=desnames)
###)
###
###   ### max haplo-type
###   Haplos$nhaplo1 <- as.numeric(Haplos$haplo1)
###   Haplos$nhaplo2 <- as.numeric(Haplos$haplo2)
###   dsort(Haplos) <- ~id+nhaplo1+nhaplo2-p
###   Haplos <- count.history(Haplos,status="p",types="1")
###   mHaplos <- subset(Haplos,lbnr__id==1)
###   bothid <- intersect(X$id, mHaplos$id)
###   X <- subset(X, id %in% bothid)
###   mHaplos <- subset(mHaplos, id %in% bothid)
###   Xhap <- merge(X, mHaplos, by.x = "id", by.y = "id")
###   mm <- grep("haplo*", names(Xhap))
###   X <- t(as.matrix(apply(Xhap[, mm], 1, designftypes)))
###   ###
###   mmud <- glm(y~factor(time)+X,data=Xhap,family=binomial())
###
###ud <- list(coef=mud$coef,se=mud$se,se.robust=mud$se.robust,mcoef=mmud$coef,kcoef=ssud$coef)
###return(ud)
###} ## }}} 
###
