
##' IPTW GLM, Inverse Probaibilty of Treatment Weighted GLM  
##'
##' Fits GLM model with treatment weights \deqn{ w(A)= \sum_a I(A=a)/P(A=a|X) }, computes
##' standard errors via influence functions that are returned as the IID argument. 
##' Propensity scores are fitted using either logistic regression (glm) or the multinomial model (mlogit) when more
##' than two categories for treatment. The treatment needs to be a factor and is identified on the rhs
##' of the "treat.model". 
##'
##' Also works with cluster argument. 
##'
##' @param formula for glm 
##' @param data data frame for risk averaging
##' @param treat.model propensity score model (binary or multinomial) 
##' @param family of glm (logistic regression)
##' @param id cluster id for standard errors 
##' @param weights may be given, and then uses weights*w(A) as the weights
##' @param estpr to estimate propensity scores and get infuence function contribution to uncertainty
##' @param pi0 fixed simple weights 
##' @param ...  arguments for glm call
##' @author Thomas Scheike
##' @export
glm_IPTW <- function(formula,data,treat.model=NULL,family=binomial(),id=NULL,weights=NULL,estpr=1,pi0=0.5,...) {# {{{
  if (!is.null(id)) {
          orig.id <- id
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
  } else { orig.id <- NULL; nid <- nrow(data); id <- 0:(nid-1); ids <- NULL}
  ### id from call coded as numeric 1 -> 
  id <- id+1
  nid <- length(unique(id))
  data$id__  <-  id

treats <- function(treatvar) {# {{{
treatvar <- droplevels(treatvar)
nlev <- nlevels(treatvar)
nlevs <- levels(treatvar)
###treatvar <- as.numeric(treatvar)
ntreatvar <- as.numeric(treatvar)
return(list(nlev=nlev,nlevs=nlevs,ntreatvar=ntreatvar))
}
# }}}

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
   treat.modelid <- update.formula(treat.model,.~.+cluster(id__))
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

out <- list(iidalpha=iidalpha,pA=pA,pal=pal,ppp=ppp,spp=spp,id=id,DPai=DPai)
return(out)
} # }}}

expit <- function(x) 1/(1+exp(-x))

treat.name <-  all.vars(treat.model)[1]
treatvar <- data[,treat.name]
if (!is.factor(treatvar)) stop(paste("treatment=",treat.name," must be coded as factor \n",sep="")); 

treats <- treats(treatvar)
if (estpr[1]==1) {
   fitt <- fittreat(treat.model,data,id,treats$ntreatvar,treats$nlev)
   iidalpha0 <- fitt$iidalpha
} else {
   ## assumes constant fixed prob over groups
   pi0 <- rep(pi0,treats$nlev)
}

wPA <- c(1/fitt$pA)
if (is.null(weights)) ww <- wPA else ww <- weights*wPA
data$ww <- ww

## fit the weighted model and bring the derivatives to bring them along
glmw <- suppressWarnings(glm(formula,data,weights=ww,family=family,...))
glm.iid <- lava::iid(glmw,id=id)
ihess <- attr(glm.iid,"bread")

###browser()
check.derivative <- 0
### for checking derivative 
if (check.derivative==1) {
###fpar <- glm(treat.model,data,family=binomial)
###mm <- as.matrix(model.matrix(treat.model,data))
###cpar <- coef(fpar)
###library(numDeriv)
###X <- model.matrix(formula,data)
###formulaS <- Event(tt,status)~tcell+platelet+age.f
###ff <- function(par) {
###pa <- 	lava::expit(mm %*% par)
###www <- 1/ifelse(data$tcell== 1, pa, 1 - pa)
######print(summary(www))
###
###pp <- binreg(Event(tt,cause)~tcell+platelet+age.f,bmt,cause=1,time=5,cens.code=5,no.opt=TRUE,
###	     weights=www, beta=coef(glmw))
###gradient  <- pp$gradient
###return(gradient)
###}
###print(ff(cpar))
###gf <- jacobian(ff,cpar)
###print(gf)
}

### iid after propensity model 
### computing  derivatives 
if (estpr[1]==1) {
X <- model.matrix(formula,data)
res <- glmw$y -predict(glmw,data,type="response") 
XY <- X*res 
###XWY <- X*ww*res
###print( apply(XWY,2,sum))
DUa  <-  t(fitt$DPai) %*% XY  
###print(DUa)
iidpal <- iidalpha0 %*% DUa ## /nid
glmw$naive.var  <- crossprod(glm.iid)
glm.iid <-   glm.iid + iidpal %*% ihess
glmw$DUa <- DUa
} 

glmw$iid <- glm.iid
glmw$var  <-  crossprod(glm.iid)
## add some arguments to comply with summary.binreg
glmw$n <- nrow(data)
glmw$ncluster <- nid
glmw$nevent <- sum(glmw$y)
class(glmw) <- c("binreg","glm")

return(glmw)
}# }}}

##' @export
predictGLM <- function(object,newdata,id=NULL,fun=NULL,link.conf=TRUE,...) {# {{{
    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, xlev = object$xlevels)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    }
    offsets <- rep(0,nrow(X))
    if (!is.null(object$offset)) offsets <- m[, grep("offset", names(m))] 

if (!is.null(object$family$linkinv)) linkinv <- object$family$linkinv  else linkinv <- function(x) x
f <- function(p) { pp <- X %*% p+offsets; return(linkinv(pp)); }
fl <- function(p) { pp <- X %*% p+offsets; return(pp); }
if (!is.null(fun))  f <- fun

if (!is.null(id)) coef <- estimate(object,id=id,...) else coef <- estimate(object,...)

if (link.conf) { 
     if (!is.null(id)) resl <- estimate(object,f=fl,id=id,...) else resl <- estimate(object,f=fl,...)
     res <- linkinv(resl$coefmat[,c(1,3,4)]) 
} else {
if (!is.null(id)) res <- estimate(object,f=f,id=id,...) else res <- estimate(object,f=f,...)
res <- res$coefmat
}

return(list(coef=coef,pred=res))
}
# }}}

##' Reporting OR (exp(coef)) from glm with binomial link and glm predictions
##' 
##' Reporting OR from glm with binomial link  and glm predictions
##' 
##' @param object glm output 
##' @param id possible id for cluster corrected standard errors
##' @param fun possible function for non-standard predictions based on object
##' @param ... arguments of estimate of lava for example level=0.95 
##' @author Thomas Scheike
##' @export
##' @examples
##' data(sTRACE)
##' sTRACE$id <- sample(1:100,nrow(sTRACE),replace=TRUE)
##'
##' model <- glm(I(status==9)~sex+factor(diabetes)+age,data=sTRACE,family=binomial)
##' summaryGLM(model)
##' summaryGLM(model,id=sTRACE$id)
##'
##' nd <- data.frame(sex=c(0,1),age=67,diabetes=1)
##' predictGLM(model,nd)
##' @aliases predictGLM 
summaryGLM <- function(object,id=NULL,fun=NULL,...) {# {{{

f <- function(p) { pp <- exp(pp); return(pp); }

if (!is.null(id)) coef <- estimate(object,id=id,...) else coef <- estimate(object,...)

if (!is.null(id)) resl <- estimate(object,id=id,...) else resl <- estimate(object,...)
res <- exp(resl$coefmat[,c(1,3,4)]) 

if (!is.null(fun))  {
if (!is.null(id)) resl <- estimate(object,id=id,f=fun,...) else resl <- estimate(object,f=fun,...)
fout <- resl$coefmat
} else fout <- NULL

return(list(coef=coef,or=res,fout=NULL))
}
# }}}

waldTest <- function(object,...)
{# {{{
Vars <-attr(terms(formula(object)),"term.labels")
## remove possible cluster term
clusterTerm<- grep("^cluster[(][A-z0-9._:]*",Vars,perl=TRUE)
if (length(clusterTerm)==1) Vars <- Vars[-clusterTerm]
Varf <- names(attr(attr(object$model,"terms"),"dataClasses"))[grep("factor",attr(attr(object$model,"terms"),"dataClasses"))]
txtVar<-names(coef(object))
###
test.matrix <- matrix(0,length(Vars),3)
colnames(test.matrix) <- c("Chisq","df","p.value")
rownames(test.matrix) <- Vars
k <- 1
for (ff in Vars) {
	pos <- grep(ff,txtVar,fixed=TRUE)
        wt <- estimate(object,as.list(pos),...)$compare
###	wt <- estimate(object,as.list(3:5))
	test.matrix[k,] <- unlist(wt[1:3])
	k <- k+1
}
return(test.matrix)
}# }}}


