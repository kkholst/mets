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
##' @param idhapweights weights following id-hap for GLM (WIP)
##' @param ... Additional arguments to lower level funtions lava::NR  optimizer or nlm
##' @author Thomas Scheike
##' @examples
##' ## some haplotypes of interest
##' types <- c("DCGCGCTCACG","DTCCGCTGACG","ITCAGTTGACG","ITCCGCTGAGG")
##' 
##' ## some haplotypes frequencies for simulations 
##' data(haplo)
##' hapfreqs <- haplo$hapfreqs 
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
##' ghaplos <- haplo$ghaplos
##' haploX  <- haplo$haploX
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
##' summary(out)
##' @export
haplo.surv.discrete <- function (X=NULL,y="y",time.name="time",Haplos=NULL,id="id",desnames=NULL,designfunc=NULL,
    beta=NULL,no.opt=FALSE,method="NR",stderr=TRUE,designMatrix=NULL,response=NULL,idhap=NULL,design.only=FALSE,
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
   ms <- mystrata(Xhap[,c("id",nmm)])
	   stratidhap <- ms
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

	obj <- function(pp,all=FALSE)
	{ # {{{

	lp <- X %*% pp 
	## plp <- family$linkinv(lp)
	plp <- expit(lp+ offiid)
	nplp <- 1-plp

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
  class(val) <- "haplosurvd"
  return(val)
} ## }}} 

## uses HaploSurvival package of github install via devtools
## devtools::install_github("scheike/HaploSurvival")
## this is only used for simulations 
## out <- simHaplo(1,100,tcoef,hapfreqs)

##' @export
summary.haplosurvd <- function(object,...) { ## {{{ 
out <- lava::estimate(object,...)
return(out)
} ## }}} 

##' @export
print.haplosurvd <- function(x,...) summary(x,...)

##' @export
vcov.haplosurvd <- function(object,...) return(object$var) 

##' @export
coef.haplosurvd <- function(object,...) return(object$coef)


