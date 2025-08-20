
##' Computes mediation weights 
##'
##' Computes mediation weights for either binary or multinomial mediators.
##' The important part is that the influence functions can be obtained to compute standard errors. 
##'
##' @param fit either glm-binomial or mlogit (mets package) 
##' @param data data frame with data
##' @param var is NULL reads mediator and exposure from formulae in the fit.
##' @param name.weight name of weights
##' @param id.name name of id variable, important for SE computations
##' @param ... Additional arguments to 
##' @author Thomas Scheike
##' @export
medweight <- function(fit,data=data,var=NULL,name.weight="weights",id.name="id",...)
{# {{{
vvars <- all.vars(fit$formula)
fexp <- data[,vvars[2]]
fmed <- data[,vvars[1]]

if (is.null(var))  { 
if (!is.factor(fexp)) stop("exposure variable must be a factor\n")
 var <- vvars[2] } else {
if (!is.factor(data[,var])) stop("first predictor of weightmodel is exposure variable and must be a factor\n")
}
if (!is.factor(fmed)) stop("mediator variable must be a factor\n")

nlev <- nlevels(fexp)
levexp <- levels(fexp)
ss <- as.numeric(fexp)
ind <- rep(1:nrow(data),each=nlev)
varn <- paste(var,0:1,sep="")
data[,varn[1]] <- data[,var]

wdata <- data[ind,]
mlev <- nlevels(fmed)
ovar <- c(1+t(outer(ss,(nlev-1):0,FUN="+") %%nlev))
wdata[,varn[2]]  <- factor(ovar,labels=levels(fexp))
## set original variable to a going through exposure levels to use predict
wdata[,var] <- wdata[,varn[2]]  

## predicting using predict function
if (inherits(fit,"glm")) { 
	pp1 <- predict(fit,wdata,type="response"); 
	y <- as.numeric(wdata[,vvars[1]])-1
	pp <- pp1*y + (1-pp1)*(1-y)
	wdata[,"basep1"] <- pp
} else if (inherits(fit,"mlogit")) pp <- predict(fit,wdata,se=0) 
else stop("Considers only binary or multinomial regression modelss glm mlogit \n");

countid <- rep(1:nlev,nrow(data))
weight0 <- rep(pp[countid==1],each=nlev)
weight1 <- pp/weight0
wdata[,name.weight] <- weight1
wdata[,id.name] <- rep(data[,id.name],each=nlev)

attr(wdata,"exposure") <- vvars[[2]]
attr(wdata,"mediator") <- vvars[[1]]
attr(wdata,"nlev-exposure") <- nlev 
attr(wdata,"countid") <- countid
return(wdata)
}# }}}

##' Mediation analysis in survival context 
##'
##' Mediation analysis in survival context  with robust standard errors taking the weights into account
##' via influence function computations. Mediator and exposure must be factors.  This is based on numerical
##' derivative wrt parameters for weighting.  See vignette for more examples. 
##'
##' @param survmodel with mediation model (binreg, aalenMets, phreg) 
##' @param weightmodel mediation model
##' @param data for computations
##' @param wdata weighted data expansion for computations
##' @param id name of id variable, important for SE computations
##' @param silent to be silent
##' @param ... Additional arguments to survival model 
##' @author Thomas Scheike
##' @examples
##' library(mets)
##' n <- 400
##' dat <- kumarsimRCT(n,rho1=0.5,rho2=0.5,rct=2,censpar=c(0,0,0,0),
##'           beta = c(-0.67, 0.59, 0.55, 0.25, 0.98, 0.18, 0.45, 0.31),
##'     treatmodel = c(-0.18, 0.56, 0.56, 0.54),restrict=1)
##' dfactor(dat) <- dnr.f~dnr
##' dfactor(dat) <- gp.f~gp
##' drename(dat) <- ttt24~"ttt24*"
##' dat$id <- 1:n
##' dat$ftime <- 1
##' 
##' weightmodel <- fit <- glm(gp.f~dnr.f+preauto+ttt24,data=dat,family=binomial)
##' wdata <- medweight(fit,data=dat)
##'
##' ### fitting models with and without mediator
##' aaMss2 <- binreg(Event(time,status)~gp+dnr+preauto+ttt24+cluster(id),data=dat,time=50,cause=2)
##' aaMss22 <- binreg(Event(time,status)~dnr+preauto+ttt24+cluster(id),data=dat,time=50,cause=2)
##'
##' ### estimating direct and indirect effects (under strong strong assumptions) 
##' aaMss <- binreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),
##'                 data=wdata,time=50,weights=wdata$weights,cause=2)
##' ## to compute standard errors , requires numDeriv
##' library(numDeriv)
##' ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
##' summary(ll)
##' ## not run bootstrap (to save time)
##' ## bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=500)
##'
##' @aliases  BootmediatorSurv 
##' @export
mediatorSurv <- function(survmodel,weightmodel,data=data,wdata=wdata,id="id",silent=TRUE,...)
{# {{{
if (is.null(weightmodel$iid)) {
idvar <- data[,id]
iidp <- iid(weightmodel,id=idvar)
ihessian <- attr(iidp,"bread")
} else { iidp <- weightmodel$iid; ihessian <- weightmodel$ihessian }

if (is.null(survmodel$iid)) iids <- iid(survmodel) else iids <- survmodel$iid

evar <- all.vars(weightmodel$formula)[2]
mvar <- all.vars(weightmodel$formula)[1]
fexp <- data[,evar]
nlev <- nlevels(fexp)
countid <- rep(1:nlev,nrow(data))

DPbetaScoreA <- function(p) {# {{{
if (inherits(weightmodel,"glm")) {
weightmodel$coefficients <- p
pp1 <- predict(weightmodel,wdata,type="response")
y <- as.numeric(wdata[,mvar])-1
pp <- pp1*y + (1-pp1)*(1-y)
} else if (inherits(weightmodel,"mlogit"))  {
weightmodel$coef <- p
pp <- predict(weightmodel,wdata,se=0)
} else stop("must be glm or mlogit weights\n") 

countid <- attr(wdata,"countid")
weight0 <- rep(pp[countid==1],each=nlev)
weight1 <- pp/weight0

model <- tail(class(survmodel),1)
## starting values coef(survmodel) and ... arguments + new weights 
ll <- c(list(formula=survmodel$formula),list(data=wdata),list(weights=weight1),list(...))
## add time variable for binreg
if (model=="binreg")  ll <- c(ll,list(time=survmodel$time),list(cause=survmodel$cause),list(cens.code=survmodel$cens.code))
## fine-gray or logit link
if (model=="cifreg")  ll <- c(ll,list(propodds=survmodel$propodds),list(cause=survmodel$cause),list(cens.code=survmodel$cens.code))
## do not optimize, just compute score for these weights at fixed survival parameter 
if (model!="aalenMets")  ll <- c(ll,list(no.opt=TRUE))
if (model!="aalenMets")  ll <- c(ll,list(beta=coef(survmodel)))

## evaluate score for the new weights 
mscore <- do.call(model,ll)
if (model=="aalenMets")  {
	## compute score for weight change  with fixed parameter
	score <- mscore$intZHdN - mscore$intZHZ %*% coef(survmodel)
} else score <- mscore$gradient

return(score)
}# }}}

###iiDD <- t(iids) %*% ( as.matrix(iidp)  %*% solve(ihessian))
###print(iiDD)
if (!silent) cat("start numerical derivative wrt weight-parameters\n")
DD <- numDeriv::jacobian(DPbetaScoreA,coef(weightmodel)) # /nrow(data)
iiDD <- survmodel$ihessian %*% DD
iid.w <-    iidp  %*% t(iiDD)
iid.tot <- iids +  iid.w
iid.tot2 <- iids -  iid.w
robvar <- crossprod(iid.tot)
robvar2 <- crossprod(iid.tot2)

if (inherits(survmodel,"binreg")) model <- survmodel$model
if (inherits(survmodel,"phreg")) model <- "exp"
if (inherits(survmodel,"aalenMets")) model <- "lin"

out <- list(coef=survmodel$coef,iid=iid.tot,
	    var=robvar,se=diag(robvar)^.5,
	    var2=robvar2,se2=diag(robvar2)^.5,
       n=survmodel$n,ncluster=survmodel$ncluster,nevent=survmodel$nevent,iDPbeta=iiDD,
       iid.w=iid.w,iid.surv=iids,model=model)
   class(out) <- c("binreg","medSurv")
   return(out)
}# }}}

##' @export
BootmediatorSurv <- function(survmodel,weightmodel,data=data,id="id",silent=TRUE,k.boot=100,...)
{# {{{

p <- length(coef(survmodel))

mediatorAn <- function(bdata,...) 
{# {{{

if (inherits(weightmodel,"glm")) {
bfit <- suppressWarnings(glm(weightmodel$formula,data=bdata,family=binomial))
} else if (inherits(weightmodel,"mlogit"))  {
bfit <- mlogit(weightmodel$formula,data=bdata)
} else stop("must be glm or mlogit weights\n") 

wdatab <- medweight(bfit,bdata)

model <- tail(class(survmodel),1)
ll <- c(list(formula=survmodel$formula),list(data=wdatab),list(weights=wdatab$weights),list(...))
if (model=="binreg")  ll <- c(ll,list(time=survmodel$time),list(cause=survmodel$cause),list(cens.code=survmodel$cens.code))
if (model=="cifreg")  ll <- c(ll,list(propodds=survmodel$propodds),list(cause=survmodel$cause),list(cens.code=survmodel$cens.code))

coefb <- do.call(model,ll)$coef
return(coefb)
}# }}}

parb <- c()
if (k.boot>=1)  {
for (i in 1:k.boot) {
bdata <- blocksample(data,nrow(data),idvar=id)
parb <- rbind(parb,mediatorAn(bdata,...))
}
bvar <- cov(parb)
} else bvar <- matrix(0,p,p)

out <- list(coef=survmodel$coef,var=bvar,se.coef=diag(bvar)^.5,parb=parb,
       n=survmodel$n,ncluster=survmodel$ncluster,nevent=survmodel$nevent)
   class(out) <- c("binreg","medSurv")
   return(out)
}# }}}

