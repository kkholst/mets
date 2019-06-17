##' Discrete time to event haplo type analysis 
##'
##' Can be used for logistic regression when time variable is "1" for all. 
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
##' 
##' Take all names that begins with "haplo" and uses for design and sub-grouping (the h's)
##' 
##' @param X 
##' @param y  name of response (binary response with logistic link)
##' @param time.name (to construct design ~factor(time)) 
##' @param Haplos (data.frame with id, haplo1, haplo2 (haplotypes (h)) and  P(h|G))
##' @param id name of id variale 
##' @param des.names names for design matrix
##' @param no.opt optimization TRUE/FALSE 
##' @param method NR, nlm 
##' @param designfunc function that computes x(h) 
##' @param designMatrix 
##' @param design.only to return only design matrices for haplo-type analyses. 
##' @param beta 
##' @param covnames names of 
##' @param ... Additional arguments to lower level funtions lava:::NR  optimizer or nlm
##' @author Thomas Scheike
##' @export
haplo.surv.discrete <- function (
###	formula = formula(data),data=sys.parent(),
    X=NULL,y="y",time.name="time",Haplos=NULL,id="id",
    desnames=NULL, designfunc=NULL,beta=NULL, 
    no.opt=FALSE,method="NR",stderr=TRUE,
    designMatrix=NULL,design.only=FALSE,
###    geno.type = NULL, geno.setup=NULL, haplo.freq = NULL,stderr=TRUE,
###    haplo.design=NULL,haplo.baseline=NULL,alpha=NULL, 
    covnames=NULL,...)
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
   stratidhap <-    as.numeric(strata(lll))
   stratidhap <- fast.approx(unique(stratidhap),stratidhap) 
###   print(cbind(Xhap$id,stratidhap))

   nidhap <- length(unique(stratidhap))
   nid <- length(unique(Haplos$id))

   hgiveng <- Haplos$p

if (!design.only) {

if (is.null(beta)) beta <- rep(0,ncol(X))

obj <- function(pp,all=FALSE)
{ # {{{
lp <- X %*% pp
plp <- expit(lp)
nplp <- 1-plp
lognp <- log(nplp)

###
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
### structure(-ploglik,gradient=-gradient,hessian=-hessian)
### hessian <- NULL 
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
  class(val) <- "haplo.surv.discrete"
  return(val)
} ## }}} 


