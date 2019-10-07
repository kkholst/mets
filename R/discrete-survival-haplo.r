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
##' @param designMatrix  gives response and designMatrix directly not implemented (mush contain: p, id, idhap)
##' @param response gives response and design directly designMatrix not implemented 
##' @param idhap name of id-hap variable to specify different haplotypes for different id 
##' @param design.only to return only design matrices for haplo-type analyses. 
##' @param covnames names of covariates to extract from object for regression
##' @param fam family of models, now binomial default and only option 
##' @param weights weights following id for GLM 
##' @param offsets following id  for GLM
##' @param weightsidhap weights following id-hap for GLM (WIP)
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
##' @aliases simTTP  predictSurvd plotSurvd
##' @export
haplo.surv.discrete <- function (X=NULL,y="y",time.name="time",Haplos=NULL,id="id",desnames=NULL,designfunc=NULL,
    beta=NULL,no.opt=FALSE,method="NR",stderr=TRUE,designMatrix=NULL,response=NULL,idhap=NULL,design.only=FALSE,
### geno.type = NULL, geno.setup=NULL,haplo.freq = NULL, haplo.design=NULL,haplo.baseline=NULL,alpha=NULL, 
    covnames=NULL,fam=binomial,weights=NULL,offsets=NULL,idhapweights=NULL,...)
{ ## {{{ 
  cond=NULL

if (is.null(designMatrix)) { 
if (!is.null(Haplos)) { ## with haplo-types {{{
  ## X: y, Xdes, id
  ## Haplos, haplo1, haplo2, id, p (HgivenG)
  bothid <- intersect(X$id,Haplos$id)
  X <- subset(X,id %in% bothid)
  Haplos <- subset(Haplos,id %in% bothid)
  ## new iid starts at 1 
  Haplos$id <- fast.approx(bothid,Haplos$id) 
  iiid <- Haplos$id-1
  X$id <- fast.approx(bothid,X$id) 
  Xo <- X

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
  if (!is.null(designfunc)) {
	  Xo <- X <- as.matrix(apply(Xhap[,mm],1,designfunc))
	  if (ncol(X)==nrow(Xhap))  X <- t(X)
	  colnames(X) <- paste("haplo",1:ncol(X),sep="")
	  X <- as.matrix(cbind(Xhap[,desnames],X))
  } else Xo <- X <- as.matrix(Xhap[,desnames])

  ## X(i)^T %*% X(i) for each row 
  X2  <- .Call("vecMatMat",X,X)$vXZ

## creates sub-index for haplo-types within each id
   nmm <- names(Xhap)[mm]
###   lll <- lapply(Xhap[,c("id",nmm)],as.numeric)
###   stratidhap <-    as.numeric(survival::strata(lll))
   ms <- mystrata(Xhap[,c("id",nmm)])
   stratidhap <- ms$sindex
###nidhap <- length(unique(stratidhap))
   nidhap <- attr(ms,"nlevel")
   nid <- length(unique(Haplos$id))

  ## weights will follow id 
  if (is.null(weights))  wiid <- rep(1,nid)  else wiid <- weights
  ## idhap weights will follow id-hap
  if (is.null(idhapweights))  whapiid <- rep(1,nidhap)  else whapiid <- idhapweights
  ## offsets can follow both Haplo or X design and will then appear in mixed design
  if (is.null(offsets))  offiid <- rep(0,nrow(X)) else offiid <- offsets

  ### to make the optimizer more flexible and use for interval censored data 
  if (is.null(cond)) cond <- rep(0,nidhap)

   hgiveng <- Haplos$p
# }}}
} else {  ## standard glm {{{ 
  ## X: y, Xdes, id

  id.name <- id
  ## id going from 1 to #id's
  id <- X$id <- fast.approx(unique(X[,id]),X[,id]) 
  Xhap <- Xo <- X
  iiid <- unique(X$id)-1
  nid <- length(iiid)
  stratidhap <- id 
  nidhap <- length(unique(stratidhap))

  ###
  response <- Xhap[,y]
  yn <- match(y,names(Xhap))
  Xhap <- Xhap[,-yn]

  X <- as.matrix(Xhap[,desnames])
  ## X(i)^T %*% X(i) for each row 
  X2  <- .Call("vecMatMat",X,X)$vXZ

  ## weights will follow id 
  if (is.null(weights))  wiid <- rep(1,nid)  else wiid <- weights
  ## idhap weights will follow id-hap
  if (is.null(idhapweights))  wph <- rep(1,nidhap)  else wph <- idhapweights
  ## offsets can follow both Haplo or X design and will then appear in mixed design
  if (is.null(offsets))  offiid <- rep(0,nrow(X)) else offiid <- offsets

  if (is.null(cond)) cond <- rep(0,nidhap)
  hgiveng <- rep(1,nid)

}# }}}
} else {# {{{
	hgiveng <- designMatrix$p
	X <- designMatrix[,desnames]
	id <- designMatrix[,id]
        iiid <- unique(id)-1

        stratidhap <- designMatrix[,idhap]
        nidhap <- length(unique(stratidhap))
        nid <- length(iiid)
	design.only <- FALSE

        ## weights will follow id 
        if (is.null(weights))  wiid <- rep(1,nid)  else wiid <- weights
        ## idhap weights will follow id-hap
        if (is.null(idhapweights))  wph <- rep(1,nidhap)  else wph <- idhapweights
        ## offsets can follow both Haplo or X design and will then appear in mixed design
        if (is.null(offsets))  offiid <- rep(0,nrow(X)) else offiid <- offsets

        if (is.null(cond)) cond <- rep(0,nidhap)

}# }}}

if (!design.only) {

	if (is.null(beta)) beta <- rep(0,ncol(X))
	expit  <- function(z) 1/(1+exp(-z)) ## expit

	obj <- function(pp,all=FALSE)
	{ # {{{

	lp <- X %*% pp 
	## plp <- family$linkinv(lp)
	plp <- expit(lp+ offiid)
	nplp <- 1-plp
	###lognp <- log(nplp)

	### logpht <-  (response - plp)/family$variance(plp)
	logpht <- log(plp)*response+log(nplp)*(1-response)
	pht <-c(exp(logpht))
	Dlogpht <-  X* c(response-plp)
	D2logpht <- c(plp/(1+exp(lp)))*X2

	ph <- c(exp(sumstrata(logpht,stratidhap-1,nidhap)))
	pg <- c(sumstrata(ph*hgiveng,iiid,nid))
	logl <- wiid*log(pg)

	## Derivative 
	Dlogph  <- apply(Dlogpht,2,sumstrata,stratidhap-1,nidhap)
	Dph  <- c(ph)*Dlogph
	Dpg <- apply(Dph*hgiveng,2,sumstrata,iiid,nid)# {{{}}}
	Dlogl    <- wiid*Dpg/pg
	DpgDpg  <- .Call("vecMatMat",Dpg,Dpg)$vXZ

	## 2nd Derivative 
	D2logph <- apply(D2logpht,2,sumstrata,stratidhap-1,nidhap)
	DphDlogph <- .Call("vecMatMat",Dph,Dlogph)$vXZ
	D2ph <- ph*D2logph+DphDlogph
	D2pg <-apply(D2ph*hgiveng,2,sumstrata,iiid,nid)
	D2logi <- wiid*(pg*D2pg-DpgDpg)/pg^2
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
			  id=iiid,
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


}  else val <- NULL

  val <- c(list(Xhap=Xhap,X=X,Haplos=Haplos),val)
  class(val) <- "survd"
  return(val)
} ## }}} 

##' @export
mystrata <- function(ll) {# {{{
	   id <- ll[[1]]
           ss <- id
	   for (j in seq(2,length(ll))) {
		   if (!is.factor(ll[[j]])) ll[[j]] <- factor(ll[[j]])
                   mm <- 1/nlevels(ll[[j]])
		   ## two decimals for each level
		   dec <- as.numeric(ll[[j]])/(nlevels(ll[[j]]))-mm/2
		   ss <- ss+dec/100^{j-2}
	   }
	   uss <- unique(ss)
	   nindex <- length(uss)
           sindex <- fast.approx(uss,ss)
	   dd <- data.frame(id=id,sindex=sindex)
           attr(dd,"nlevel") <- nindex
	   return(dd)
} # }}}

##' @export
summary.survd <- function(x,...) return(lava:::estimate(x,...)) 

##' @export
print.survd <- function(x,...) return(lava:::estimate(x,...)) 

##' @export
vcov.survd <- function(x) return(x$var) 

##' @export
coef.survd <- function(x) return(x$coef)

##' @export
simTTP <- function(coef=NULL,n=100,Xglm=NULL,times=NULL)
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
predictSurvd <- function(ds,Z,times=1:6,se=FALSE,type="prob")
{# {{{
  if (!is.null(Z)) n <- nrow(Z) 
  if (!is.null(Z)) data <- Z else data <- data.frame(id=1:n)
  Z <- data.frame(Z)
  Z$id <- 1:n
  ccc <- ds$coef

  if (!se) {# {{{{{{
	  data <- Z
	  if (!is.null(times)) {
	     timesf <- data.frame(times=rep(times,n),id=rep(1:n,each=length(times)))
	     data <- merge(data,timesf,by.x="id",by.y="id")
	     mt <- model.matrix(~factor(times),data)
	     nm <- match(c("id","times"),names(data))
	     Z <- cbind(mt,data[,-nm])
	  }
	  if (ncol(Z)!=length(c(ccc))) {
		  print(head(Z))
		  print(ccc)
		  stop("dimension of Z not consistent with length of coefficients"); 
	  }

	  p <- c(expit(as.matrix(Z) %*% ccc))

	  preds <- data.frame(p=p,id=data$id,times=data$times)
	  survt <- exp(cumsumstrata(log(1-preds$p),data$id-1,6))
	  if (type=="prob") pred <- 1-survt
	  if (type=="surv") pred <- survt
	  if (type=="hazard") pred <- p
	  if (type=="rrm") { ## restricted residual mean 
		  ll <- length(survt)
	        pred <- cumsum(c(1,survt[-ll]))
	  }
	  preds <- cbind(preds,pred)
# }}}
  } else {# {{{

    expit <- function(p) exp(p)/(1+exp(p))
    Ft <- function(p)
    {
	   xp <- as.matrix(Zi) %*% p
	   lam <- expit(xp)
	   st <- cumprod(1-lam)
	   if (type=="prob") st <- 1-st 
	   if (type=="surv") st <- st 
	   if (type=="hazard") st <- lam
           if (type=="rrm") { ## restricted residual mean 
		ll <- length(st)
	        st <- cumsum(c(1,st[-ll]))
	   }
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
     if (is.null(ds$var)) covv <- vcov(ds)  else covv <- ds$var
     eud <- estimate(coef=ds$coef,vcov=covv,f=function(p) Ft(p))
     cmat <- data.frame(eud$coefmat)
     cmat$id <- i
     cmat$times <- times
     names(cmat)[1:4] <- c("pred","se","lower","upper")
     preds <- rbind(preds,cmat)
  } 

  }# }}}

return(preds)
}# }}}

## }}} 

##' @export
plotSurvd <- function(ds,ids=NULL,add=FALSE,se=FALSE,cols=NULL,ltys=NULL,...)
{# {{{

 if (is.null(ids)) ids <- unique(ds$id)
 if (is.null(cols)) cols <- 1:length(ids)
 if (is.null(ltys)) ltys <- 1:length(ids)

  k <- 1
  fplot <- 0
  for (i in ids) {
	  timei <- ds$time[ds$id==i]
	  predi <- ds$pred[ds$id==i]

 	  if (fplot==0) {
	  if (!add) plot(timei,predi,type="s",col=cols[k],lty=ltys[k],...)
	  if (add) lines(timei,predi,type="s",col=cols[k],lty=ltys[k],...)
	  fplot <- 1
	  } else lines(timei,predi,type="s",col=cols[k],lty=ltys[k],...)

          if (se) {
	  loweri <- ds$lower[ds$id==i]
	  upperi <- ds$upper[ds$id==i]
	  plotConfRegion(timei,cbind(loweri,upperi),col=cols[k])
	  }
	  k <- k+1
  }

} ## }}} 


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
###   sud <- simTTP(coef=tcoef,Xglm=X,n=n,times=1:6)
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

##' @export
Int_discrete <- function (time, time2 , id = NULL, max.time=6, cens.code = 0, ...)
{# {{{
 ## T in interval  [time,time2] so both time and time2 possible, 
 ## T>time same as T in [time+1,Inf]
 ## same as P(T > time-1)  - P(T > time2) = P(T > time-1) ( 1 - P(T > time2 | T> time-1))
 ##  P(T > time2+1 | T> time-1) = prod_{j=time}^time2 (1- lambda(j))
 otime2 <- time2
 time2[time2==Inf] <- max.time 
 if (is.null(id)) id <- seq_along(time)

 dtp <- c()
 for (i in seq_along(time))
 {
    if (otime2[i]!=Inf & time2[i]<=max.time) 
    loc <- cbind(seq(time[i],time2[i]),id[i],1)
    else loc <- c()
 if (time[i]>1) loc<- rbind(loc,cbind(seq(1,min(time[i]-1,max.time)),id[i],0))
 dtp <- rbind(dtp,loc)
 }
 colnames(dtp) <- c("time","id","cond")
 dtp <- data.frame(dtp)
 dtp$status <- 0

 dsort(dtp) <- ~id+time
 return(dtp)

}# }}}

##' @export
interval.surv.discrete <- function (formula,data,beta=NULL,no.opt=FALSE,method="NR",max.time=6,
				    stderr=TRUE,fam=binomial,weights=NULL,offsets=NULL,...)
{ ## {{{ 

  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")

 if (ncol(Y)==2) {
	time2 <- eventtime <- Y[,2]
###	status <- delta  <- Y[,2]
	entrytime <- Y[,1]
	left <- 0
    } else {
	time2 <- eventtime <- Y[,2]
	status <- delta  <- Y[,3]
	entrytime <- Y[,1]
	left <- 1
	if (max(entrytime)==0) left <- 0
    }

  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    pos.cluster <- ts$terms
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  } else pos.cluster <- NULL
  if (!is.null(attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}

  datint <- Int_discrete(entrytime,time2,max.time=max.time)
  dfactor(datint) <- cond~cond

  Xcov <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
  Xcov <- Xcov[,-intpos,drop=FALSE]
  if (ncol(Xcov)>0) Xcov.names <- colnames(Xcov) else Xcov.names <- NULL
  if (ncol(Xcov)==0) Xcov <- NULL 
  ###

  if (!is.null(Xcov)) {
	  Xcov <- Xcov[datint$id,]
	  datint <- cbind(datint,Xcov)
  } 

  ## also take covariates from X and put on data 
  Xtime <- model.matrix(~factor(time),datint)
  Xtime.names <- colnames(Xtime)
  X <- cbind(Xtime,Xcov)
  Xnames <- c(Xtime.names,Xcov.names)

  id <- datint$id
  ## covariates 
###  X <- X[datint$id,]
  X2  <- .Call("vecMatMat",X,X)$vXZ

 ## id with cond 
  ms  <-  mystrata(datint[,c("id","cond")])

## id going from 1 to #id's
  iiid <- ms$sindex -1
  nid <- attr(ms,"nlevel")

  ## weights/offets will follow id 
  if (is.null(weights))  wiid <- rep(1,nid)  else wiid <- weights
  if (is.null(offsets))  offiid <- rep(0,nrow(X)) else offiid <- offsets
  if (is.null(beta)) beta <- rep(0,ncol(X))
  expit  <- function(z) 1/(1+exp(-z)) ## expit

  ### set up cond after iiid 
  cond <- as.numeric(datint$cond)-1
  condg <- c(sumstrata(cond,iiid,nid))
  condg[condg>0] <- 1
  response <- rep(0,nrow(X))

obj <- function(pp,all=FALSE)
	{ # {{{

	lp <- X %*% pp
	plp <- expit(lp+ offiid)
	nplp <- 1-plp

	### logpht <-  (response - plp)/family$variance(plp)
	logpgt <- log(plp)*response+log(nplp)*(1-response)
	pgt <-c(exp(logpgt))
	Dlogpgt <-  X* c(response-plp)
	D2logpgt <- c(plp/(1+exp(lp)))*X2
	pg <- c(exp(sumstrata(logpgt,iiid,nid)))
	pgg <- log(pg)*(condg==0)+ log(1-pg)*(condg==1)
	logl <- wiid*pgg

	## Derivative 
	Dlogpg  <- apply(Dlogpgt,2,sumstrata,iiid,nid)
	Dpg  <- c(pg)*Dlogpg
	Dlogl    <- wiid*((condg==0)*Dpg/pg-(condg==1)*Dpg/(1-pg))

	DpgDpg  <- .Call("vecMatMat",Dpg,Dpg)$vXZ

	## 2nd Derivative 
	D2logpg <- apply(D2logpgt,2,sumstrata,iiid,nid)
	DpgDlogp <- .Call("vecMatMat",Dpg,Dlogpg)$vXZ
	D2pg <- pg*D2logpg+DpgDlogp
	D2logi <-(condg==0)*(pg*D2pg-DpgDpg)/pg^2+(condg==1)*((1-pg)*D2pg-DpgDpg)/(1-pg)^2
	D2log <- apply(wiid*D2logi,2,sum)
	D2log <- matrix(D2log,ncol(X),ncol(X))
	D2log

	ploglik <- sum(logl)
	gradient <- apply(Dlogl,2,sum)
	hessian <- D2log

	ploglik
	gradient
	hessian


	  if (all) {
	      ihess <- solve(hessian)
	      beta.iid <- Dlogl %*% ihess ## %*% t(Dlogl) 
	      robvar <- crossprod(beta.iid)
	      val <- list(par=pp,ploglik=ploglik,gradient=gradient,hessian=hessian,ihessian=ihess,
			  iid=beta.iid,robvar=robvar,var=ihess,
			  id=iiid,
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


  if (length(val$coef)==length(Xnames)) names(val$coef) <- Xnames
  val <- c(list(X=X,condg=condg,datint=datint),val)
  class(val) <- "survd"
  return(val)
} ## }}} 



