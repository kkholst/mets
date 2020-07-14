##' Double CIF Fine-Gray model with two causes
##'
##' Estimation based on derived hazards and recursive estimating equations.
##' fits two parametrizations
##' 1)
##'  \deqn{
##' F_1(t,X) = 1 - \exp( \exp( X^T \beta ) \Lambda_1(t))
##' }
##' and
##'  \deqn{
##' F_2(t,X_2) = 1 - \exp( \exp( X_2^T \beta_2 ) \Lambda_2(t))
##' }
##' or restricted version
##' 2)
##'  \deqn{
##' F_1(t,X) = 1 - \exp( \exp( X^T \beta ) \Lambda_1(t))
##' }
##' and
##'  \deqn{
##' F_2(t,X_2,X) = ( 1 - \exp( \exp( X_2^T \beta_2 ) \Lambda_2(t)) ) (1 - F_1(\infty,X))
##' }
##'
##' @param formula formula with 'Event'
##' @param data data frame
##' @param offset offsets for cox model
##' @param weights weights for Cox score equations
##' @param X2 specifies the regression design for second CIF model
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' res <- 0
##' data(bmt)
##' bmt$age2 <- bmt$age
##' newdata <- bmt[1:19,]
##' if (interactive()) par(mfrow=c(5,3))
##'
##' ## same X1 and X2
##' pr2 <- doubleFGR(Event(time,cause)~age+platelet,data=bmt,restrict=res)
##' if (interactive()) {
##'   bplotdFG(pr2,cause=1)
##'   bplotdFG(pr2,cause=2,add=TRUE)
##' }
##' pp21 <- predictdFG(pr2,newdata=newdata)
##' pp22 <- predictdFG(pr2,newdata=newdata,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##' pp21 <- predictdFG(pr2)
##' pp22 <- predictdFG(pr2,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##'
##' pr2 <- doubleFGR(Event(time,cause)~strata(platelet),data=bmt,restrict=res)
##' if (interactive()) {
##'   bplotdFG(pr2,cause=1)
##'   bplotdFG(pr2,cause=2,add=TRUE)
##' }
##' pp21 <- predictdFG(pr2,newdata=newdata)
##' pp22 <- predictdFG(pr2,,newdata=newdata,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##' pp21 <- predictdFG(pr2)
##' pp22 <- predictdFG(pr2,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##'
##' ## different X1 and X2
##' pr2 <- doubleFGR(Event(time,cause)~age+platelet+age2,data=bmt,X2=3,restrict=res)
##' if (interactive()) {
##'   bplotdFG(pr2,cause=1)
##'   bplotdFG(pr2,cause=2,add=TRUE)
##' }
##' pp21 <- predictdFG(pr2,newdata=newdata)
##' pp22 <- predictdFG(pr2,newdata=newdata,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##' pp21 <- predictdFG(pr2)
##' pp22 <- predictdFG(pr2,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##'
##' ### uden X1
##' pr2 <- doubleFGR(Event(time,cause)~age+platelet,data=bmt,X2=1:2,restrict=res)
##' if (interactive()) {
##'   bplotdFG(pr2,cause=1)
##'   bplotdFG(pr2,cause=2,add=TRUE)
##' }
##' pp21 <- predictdFG(pr2,newdata=newdata)
##' pp22 <- predictdFG(pr2,newdata=newdata,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##' pp21 <- predictdFG(pr2)
##' p22 <- predictdFG(pr2,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##'
##' ### without X2
##' pr2 <- doubleFGR(Event(time,cause)~age+platelet,data=bmt,X2=0,restrict=res)
##' if (interactive()) {
##'   bplotdFG(pr2,cause=1)
##'   bplotdFG(pr2,cause=2,add=TRUE)
##' }
##' pp21 <- predictdFG(pr2,newdata=newdata)
##' pp22 <- predictdFG(pr2,newdata=newdata,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##' pp21 <- predictdFG(pr2)
##' pp22 <- predictdFG(pr2,cause=2)
##' if (interactive()) {
##'   plot(pp21)
##'   plot(pp22,add=TRUE,col=2)
##' }
##'
##' @aliases bplotdFG predictdFG
##' @export
doubleFGR <- function(formula,data,offset=NULL,weights=NULL,X2=NULL,...) {# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
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
  if (!is.null(attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}
###  if (!is.null(attributes(Terms)$specials$offset)) {
###    ts <- survival::untangle.specials(Terms, "offset")
###    pos.offset <- ts$terms
###    Terms  <- Terms[-ts$terms]
###    offset <- m[[ts$vars]]
###  }  else pos.offset <- NULL
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  X2call <- X2;
  if (!is.null(X2)) {
	 if (X2[1]==0)  X2 <- matrix(nrow=0,ncol=0)
	 else { X2 <-  X[,X2call,drop=FALSE];
                X <- X[,-X2call,drop=FALSE];
	 }
  } else X2 <- X
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  if (ncol(X2)==0) X2 <- matrix(nrow=0,ncol=0)
  res <- c(doubleFG01R(X,X2,entry,exit,status,id,strata,offset,weights,strata.name,X2call=X2call,...),
   list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster))
  class(res) <- c("doubleFG")

  res
}# }}}

doubleFG01R <- function(X,X2, entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
             strata.name=NULL,cumhaz=TRUE,
             beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
	     restrict=0,case.weights=NULL,X2call=NULL,...) {# {{{
  p1 <- ncol(X);
  p2 <- ncol(X2)
  p <- p1+p2
  if (missing(beta))  beta <- rep(0,p)
  if (p1==0) X <- cbind(rep(0,length(exit)))
  if (p2==0) X2 <- cbind(rep(0,length(exit)))
  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }
  if (is.null(offset)) offset <- rep(0,length(exit))
  if (is.null(weights)) weights <- rep(1,length(exit))
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given
  if (!is.null(Z)) Zcall <- Z

  ## possible casewights to use for bootstrapping and other things
  if (is.null(case.weights)) case.weights <- rep(1,length(exit))

  trunc <- (!is.null(entry))
  if (!trunc) entry <- rep(0,length(exit))

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1;
   ## orginal id coding into integers
   id.orig <- id+1;

   dd <- .Call("FastCoxPrepStrata",entry,exit,status,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")
   if (!is.null(X2call))  {
   dd2 <- .Call("FastCoxPrepStrata",entry,exit,status,X2,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")
   dd$X2 <- dd2$X
   dd$XX2 <- dd2$XX
   } else {
      dd$X2 <- dd$X
      dd$XX2 <- dd$XX
   }
   if (p1==0)  dd$X <- matrix(nrow=0,ncol=0)
   if (p2==0)  dd$X2 <- matrix(nrow=0,ncol=0)

   cause <- status[dd$ord+1]
   dd$cause <- cause
   dd$nstrata <- nstrata

	obj <- function(pp,U=FALSE,all=FALSE) {# {{{
	  val <- with(dd, doubleFGstrataR(pp,X,XX,X2,XX2,sign,cause,jumps,strata,nstrata,weights,offset,ZX,caseweights,restrict))

	  if (all) {
	      val$time <- dd$time
	      val$cause <- dd$cause
	      val$ord <- dd$ord+1
	      val$jumps <- dd$jumps+1
	      val$jumptimes <- val$time[val$jumps]
	      val$weightsJ <- dd$weights[val$jumps]
	      val$case.weights <- dd$case.weights[val$jumps]
	      val$strata.jumps <- val$strata[val$jumps]
	      val$nevent <- length(val$S0)
	      val$nstrata <- dd$nstrata
	      val$strata <- dd$strata
	      return(val)
	  }
	 with(val,structure(-ploglik,gradient=-gradient,hessian=-hessian))
	}# }}}

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
      names(cc) <- c(colnames(X),colnames(X2))
      if (!stderr) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(rep(0,2),all=TRUE)
  }

  se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL
  II <- NULL
  ### computes Breslow estimator
  if (cumhaz==TRUE) { # {{{
	 if (no.opt==FALSE & p!=0) {
               II <- - tryCatch(solve(val$hessian),error=
	              function(e) matrix(0,nrow(val$hessian),ncol(val$hessian)) )
	 } else II <- matrix(0,p,p)
	 strata <- val$strata[val$jumps]
	 nstrata <- val$nstrata
	 jumptimes <- val$jumptimes

	 ## Brewslow estimator
###	 cumhaz <- cbind(jumptimes,cumsumstrata(1/val$S0,strata,nstrata))
	 cumhaz <- cbind(jumptimes,val$base12$base)

###	 if ((no.opt==FALSE & p!=0)) {
###	     DLambeta.t <- apply(val$E/c(val$S0),2,cumsumstrata,strata,nstrata)
###	     varbetat <-   rowSums((DLambeta.t %*% II)*DLambeta.t)
###	 ### covv <-  apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
###	 } else varbetat <- 0
###	 var.cumhaz <- cumsumstrata(1/val$S0^2,strata,nstrata)+varbetat
###	 se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)
	 se.cumhaz <- cumhaz

###	 colnames(cumhaz1)    <- c("time","cumhaz")
###	 colnames(cumhaz2)    <- c("time","cumhaz")
###	 colnames(se.cumhaz) <- c("time","se.cumhaz")
 } # }}}
 else {cumhaz <- se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL}

 II <- - tryCatch(solve(val$hessian),error=
        function(e) matrix(0,nrow(val$hessian),ncol(val$hessian)) )

  res <- c(val,
           list(cox.prep=dd,
		strata.call=strata.call, strata.level=strata.level,
                entry=entry,
                exit=exit,
                status=status,
                p=p,
                X=X,
                X2=X2,
		X2call=X2call,
		offsets=offset,
		weights=weights,
                id=id.orig,
		opt=opt,
		var=II,
		cumhaz=cumhaz, se.cumhaz=se.cumhaz,
		lcumhaz=lcumhaz, lse.cumhaz=lse.cumhaz,
		II=II,
		strata.name=strata.name,propodds=propodds,restrict=restrict))
  class(res) <- "doubleFG"
  res
}# }}}

## double FG in in R version
doubleFGstrataR <- function(beta, X, XX, X2, XX2, Sign, cause, Jumps, strata, nstrata, weights, offsets, ZX, caseweights,restrict)
{# {{{
	p=length(beta)
	p1  <-  ncol(X)
	p2 <- ncol(X2)
	if (nrow(X)==0) X <- matrix(1,length(strata),1)
	if (nrow(X2)==0) X2 <- matrix(1,length(strata),1)
	strata=c(strata)
###	X1 <- X[,1:p,drop=FALSE]
###	X2 <- X[,(p+1):(2*p),drop=FALSE]
###	print(c(p1,p2))
	if (p1==0) beta1 <- 0 else beta1 <- beta[1:p1]
	if (p2==0) beta2 <- 0 else beta2 <- beta[(p1+1):p]
	Xb1 = c(X %*% beta1+offsets)
	Xb2 = c(X2 %*% beta2+offsets)
	eXb1 = c(exp(Xb1)*weights);
	eXb2 = c(exp(Xb2)*weights);
	if (nrow((Sign))==length(eXb1)) { ## Truncation
		eXb1 = c(Sign)*eXb1;
		eXb2 = c(Sign)*eXb2;
	}

	S01 = c(revcumsumstrata(eXb1,strata,nstrata))
	S02 = c(revcumsumstrata(eXb2,strata,nstrata))
	E1=apply(eXb1*as.matrix(X),2,revcumsumstrata,strata,nstrata)/S01;
	E2=apply(eXb2*as.matrix(X2),2,revcumsumstrata,strata,nstrata)/S02;
	Jumps=Jumps+1
	causeJ <- cause[Jumps]
	Jumps1 <- Jumps[causeJ==1]
	Jumps2 <- Jumps[causeJ==2]

	## both E's, SO's, X same
	S0J <- cbind(S01,S02)[Jumps,]
	eXbJ <- cbind(eXb1,eXb2)[Jumps,]
	EJ <- cbind(E1,E2)[Jumps,]
	XJ <- X[Jumps,]
	X2J <- X2[Jumps,]

	E1 = E1[Jumps1,,drop=FALSE];
	E2 = E2[Jumps2,,drop=FALSE];
        E21=.Call("vecMatMat",E1,E1)$vXZ;
        E22=.Call("vecMatMat",E2,E2)$vXZ;

	XX21=apply(XX*eXb1,2,revcumsumstrata,strata,nstrata)/S01;
	XX22=apply(XX2*eXb2,2,revcumsumstrata,strata,nstrata)/S02;
	XX21 = XX21[Jumps1,,drop=FALSE];
	XX22 = XX22[Jumps2,,drop=FALSE];

	weightsJ=weights[Jumps];
	caseweightsJ=caseweights[Jumps];

	## compute recursive weights given beta1, beta2
        base12 <- .Call("cumsumstrataDFGR",weightsJ,S0J,causeJ,strata[Jumps],nstrata,eXbJ)
	if (restrict>0) {
		for (i in 1:restrict) {
		   Lam1inf <- tailstrata(strata[Jumps],nstrata);
		   Lam1inf  <- base12$base[Lam1inf,1]
		   base12 <- .Call("cumsumstrataDFGRestrictR",weightsJ,S0J,causeJ,strata[Jumps],nstrata,eXbJ,Lam1inf);
	        }
	}

###	DLam  <-  .Call("DLambetaDFGR",weightsJ,S0J,causeJ,EJ,XJ,strata[Jumps],nstrata,eXbJ);

	grad1 = (X[Jumps1,,drop=FALSE]-E1);         ## Score
	grad2 = (X2[Jumps2,,drop=FALSE]-E2);        ## Score

	### weights for Fg1 an Fg2 models
	pow1 <- base12$pow1[causeJ==1]
	pow2 <- base12$pow2[causeJ==2]
	cw1 <- (caseweightsJ*weightsJ)[causeJ==1]
	cw2 <- (caseweightsJ*weightsJ)[causeJ==2]
	S012 = S01[Jumps1]/(pow1*cw1);        ## S0 with weights
	S022 = S02[Jumps2]/(pow2*cw2);        ## S0 with weights

	val =  sum(pow1*cw1*(Xb1[Jumps1]-log(S012)))+
	       sum(pow2*cw2*(Xb2[Jumps2]-log(S022))); ## Partial log-likelihood

	grad21= grad1*(pow1*cw1);                   ## score  with weights
	grad22= grad2*(pow2*cw2);                   ## score  with weights
	grad <- c(apply(grad21,2,sum),apply(grad22,2,sum))
	gradient <- grad
	if (p1==0) gradient <- grad[-1]
	if (p2==0) gradient  <-  grad[1:p1]


	## no weights
###	val2 = caseweightsJ*weightsJ*val;    ## Partial log-likelihood with weights
	val2 = val;                          ## Partial log-likelihood with weights

 	hesst1 = -(XX21-E21);               ## hessian contributions in jump times
	hesst2 = -(XX22-E22);               ## hessian contributions in jump times
###	hess  = matrix(apply(hesst,2,sum),p,p);
	hesst12 = hesst1*(cw1)*pow1;        ## hessian over time with weights
	hesst22 = hesst2*(cw2)*pow2;        ## hessian over time with weights

	## missing some derivative terms for hessian (due to pow=w(beta,\Lam(beta,t-))

	## setup hessian matrix
	hess12 = matrix(apply(hesst12,2,sum),p1,p1);         ## hessian with weights
	hess22 = matrix(apply(hesst22,2,sum),p2,p2);         ## hessian with weights
	hess2 <- matrix(0,p1+p2,p1+p2)
	pd2 <- p1
	p <- length(beta)
	if (p1>0) hess2[1:p1,1:p1] <- hess12
	if (p2>0) hess2[(p1+1):p,(p1+1):p] <- hess22

	out=list(jumps=Jumps, ploglik=sum(val2),U=grad2,base12=base12,
		 gradient=matrix(gradient,1,p), hessian=hess2,
		 ##hessianttime=hesst2, S2S0=XX2,
		 E=EJ, S0=S0J
		 )
	return(out)
}# }}}

##' @export
coef.doubleFG <- function(object,...) object$coef

##' @export
vcov.doubleFG <- function(object,...) object$var

##' @export
summary.doubleFG <- function(object,...) estimate(object)

##' @export
print.doubleFG <- function(x,...) estimate(x)

##' @export
bplotdFG <- function(x,cause=1,...)  {# {{{
	x$cumhaz <- x$cumhaz[,c(1,cause+1)]
	basehazplot.phreg(x,...)
}# }}}

##' @export
predictdFG <- function(x,cause=1,se=FALSE,times=NULL,...)  {# {{{
	cumhaz <- x$cumhaz
	coef <- x$coef
	p <- x$p
        X2call <- x$X2call;

	### finding the right coefficients
	if (p>0) {
              if (!is.null(X2call))  {
		 if (X2call[1]==0)  {
			 p2 <- 0; p1 <- p;
		         cc1 <- 1:p; cc2 <- 0
		 } else {
	             p2 <- length(X2call);
	             p1 <- p-p2;
		     cc1 <- (1:p)[-X2call]
		     cc2 <- (1:p)[X2call]
		     if (p2==0) x$X <- x$X
		     if (p1==0)  x$X <- x$X2
		     if (p1!=0 & p2!=0) x$X <- cbind(x$X,x$X2)
		 }
         } else {
	         model.frame2 <- x$model.frame
		 p1 <- p/2
		 cc1 <- 1:p1
		 cc2 <- (p1+1):p
		 p <- p1
	 }
	if (cause==1) {
	      if (is.null(X2call)) x$coef <- x$coef[cc1]
	      else x$coef[cc2] <- 0
	}
	if (cause==2) {
	      if (is.null(X2call)) x$coef <- x$coef[cc2]
	      else x$coef[cc1] <- 0
	}
	}

	x$cumhaz <- x$cumhaz[,c(1,cause+1)]

###	if (p>0) {
###           ## sets coefficients to 0 for other cause
###	   print("====================")
###	   print(x$coef)
###	   print(x$p)
###	   print(head(x$model.frame))
###	   print(head(x$X))
###	   print(x$cumhaz)
###	}
	class(x) <- c("phreg","cifreg")
	pll <- predict(x,se=se,times=times,...)

	if (x$restrict>0 & cause==2) {
		x$cumhaz <- cumhaz[,1:2]
                if (is.null(X2call)) x$coef <- coef[cc1] else { x$coef <- coef;  x$coef[cc2] <- 0; }
		times <- max(cumhaz[,1])
		mm <- "individual.time" %in% names(list(...))
		if (mm) times <- rep(times,nrow(x$model.frame))
		p1ll <- predict(x,se=se,times=times,...)
		cif1 <- c(1-p1ll$cif)
		pll$cif <- pll$cif*cif1
	}
	return(pll)
}# }}}


##' @export
simul.cifs <- function(n,rho1,rho2,beta,rc=0.5,depcens=0,rcZ=0.5,bin=1,type=
		       c("cloglog","logistic"),rate=1) {# {{{
p=length(beta)/2
tt <- seq(0,6,by=0.1)
if (length(rate)==1) rate <- rep(rate,2)
Lam1 <- rho1*(1-exp(-tt/rate[1]))
Lam2 <- rho2*(1-exp(-tt/rate[2]))

if (length(bin)==1) bin <- rep(bin,2)

Z=cbind((bin[1]==1)*(2*rbinom(n,1,1/2)-1)+(bin[1]==0)*rnorm(n),(bin[2]==1)*(rbinom(n,1,1/2))+(bin[2]==0)*rnorm(n))
colnames(Z) <- paste("Z",1:2,sep="")

cif1 <- setup.cif(cbind(tt,Lam1),beta[1:2],Znames=colnames(Z),type=type[1])
cif2 <- setup.cif(cbind(tt,Lam2),beta[3:4],Znames=colnames(Z),type=type[1])
###
data <- sim.cifsRestrict(list(cif1,cif2),n,Z=Z)

if (depcens==0) censor=pmin(rexp(n,1)*(1/rc),6) else censor=pmin(rexp(n,1)*(1/(rc*exp(rcZ*Z[,1]))),6)

status=data$status*(data$time<=censor)
time=pmin(data$time,censor)
data <- data.frame(time=time,status=status)
return(cbind(data,Z))

}# }}}

simul.mod <- function(n,rho1,rho2,beta,rc=0.5,k=1,depcens=0) {# {{{
p=length(beta)/2
tt <- seq(0,6,by=0.1)
Lam1 <- rho1*(1-exp(-tt))
Lam2 <- rho2*(1-exp(-tt))

Z=cbind(2*rbinom(n,1,1/2)-1,rnorm(n))
colnames(Z) <- paste("Z",1:2,sep="")
cif1 <- setup.cif(cbind(tt,Lam1),beta[1:2],Znames=colnames(Z),type="cloglog")
cif2 <- setup.cif(cbind(tt,Lam2),beta[3:4],Znames=colnames(Z),type="cloglog")
###
data <- sim.cifsRestrict(list(cif1,cif2),n,Z=Z)

censhaz  <-  cbind(tt,k*tt)
if (depcens==1) {
datc <- rchaz(censhaz,exp(Z[,1]*rc))
} else datc <- rchaz(censhaz,n=n)

data$time <- pmin(data$time,datc$time)
data$status <- ifelse(data$time<datc$time,data$status,0)
dsort(data) <- ~time

return(data)
}# }}}

##' @export
strataC <- survival:::strata

##' Augmentation for Fine-Gray model based on stratified NPMLE Cif (Aalen-Johansen) 
##'
##' Computes  the augmentation term for each individual as well as the sum
##' \deqn{
##' A(\beta) = \int H(t,X,\beta) \frac{F_2^*(t,s)}{S^*(t,s)} \frac{1}{G_c(t)} dM_c
##' }
##' with 
##' \deqn{
##' H(t,X,\beta) = \int_t^\infty (X - E(\beta,t) ) G_c(t) d\Lambda_1^*i(t,s)
##' }
##' using a KM for \deqn{G_c(t)} and a working model for cumulative baseline
##' related to \deqn{F_1^*(t,s)} and \deqn{s} is strata, 
##' \deqn{S^*(t,s) = 1 - F_1^*(t,s) - F_2^*(t,s)}, and
##' \deqn{E(\beta^p,t)} is given.
##'
##' After a couple of iterations we end up with a solution of 
##' \deqn{
##' \int (X - E(\beta) ) Y_1(t) w(t) dM_1 + A(\beta)
##' }
##' the augmented FG-score. 
##'
##' Standard errors computed under assumption of correct \deqn{G_c} model.
##'
##' @param formula formula with 'Event', strata model for CIF given by strata, and strataC specifies censoring strata
##' @param data data frame
##' @param offset offsets for cox model
##' @param data data frame
##' @param E from FG-model 
##' @param cause of interest 
##' @param cens.code code of censoring 
##' @param km to use Kaplan-Meier
##' @param case.weights weights for FG score equations (that follow dN_1) 
##' @param weights weights for FG score equations
##' @param offset offsets for FG   model
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' rho1 <- 0.2; rho2 <- 0.9
##' n <- 200
##' beta=c(0.0,-0.1,-0.5,0.3)
##' dats <- simul.cifs(n,rho1,rho2,beta)
##' dsort(dats) <- ~time
##' fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
##' 
##' fgaugS <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fg$E)
##' fgaugS2 <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fgaugS$E)
##' 
##' @aliases strataC  simul.cifs
##' @export
FG_AugmentCifstrata <- function(formula,data=data,E=NULL,cause=NULL,cens.code=0,km=TRUE,case.weights=NULL,weights=NULL,offset=NULL,...)
{# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset","strataC")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (class(Y)!="Event") stop("Expected a 'Event'-object")
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
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else strata.name <- NULL
  if (!is.null(stratapos <- attributes(Terms)$specials$strataC)) {
    ts <- survival::untangle.specials(Terms, "strataC")
    Terms  <- Terms[-ts$terms]
    strataC <- as.numeric(m[[ts$vars]])-1
    strataC.name <- ts$vars
  }  else { strataC <- NULL; strataC.name <- NULL}

  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms  <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
  X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  id.orig <- id; 
  if (!is.null(id)) {
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(exit))-1; 

  p <- ncol(X)
  beta <- NULL
  if (is.null(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(exit)))

  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }

if (is.null(strataC)) { strataC <- rep(0,length(exit)); nstrataC <- 1; strataC.level <- NULL; } else {
	  strataC.level <- levels(strataC)
	  ustrataC <- sort(unique(strataC))
	  nstrataC <- length(ustrataC)
	  strataC.values <- ustrataC
      if (is.numeric(strataC)) strataC <-  fast.approx(ustrataC,strataC)-1 else  {
      strataC <- as.integer(factor(strataC,labels=seq(nstrataC)))-1
    }
  }

  cens.strata <- strataC
  cens.nstrata <- nstrataC 

  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
  strata.call <- strata
  Z <- NULL
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z

  ## possible casewights to use for bootstrapping and other things
  if (is.null(case.weights)) case.weights <- rep(1,length(exit)) 

  trunc <- (!is.null(entry))
  if (!trunc) entry <- rep(0,length(exit))
  call.id <- id

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1; 
   ## orginal id coding into integers 
   id.orig <- id+1; 


  statusC <- (status==cens.code) 
  statusE <- (status==cause) 
  if (sum(statusE)==0) stop("No events of type 1 before time \n"); 

  ## sorting after time and statusC, but event times unique after order
  Zcall <- cbind(status,strata)
  dd <- .Call("FastCoxPrepStrata",entry,exit,statusC,X,id, 
	     trunc,strataC,weights,offset,Zcall,case.weights,PACKAGE="mets")

  Z <- dd$X
  jumps <- dd$jumps+1
  xxstrataC <- c(dd$strata)
  xxstatus  <- dd$Z[,1]
  xxstrata  <- dd$Z[,2]
  other <- which((!(xxstatus %in% c(cens.code,cause)) ) )               
  jumps1 <- which(xxstatus==cause)
  jumpsD <- which(xxstatus!=cens.code)
  rr <- c(dd$sign*exp(dd$offset))
  ## S0 after strata
  S0 = c(revcumsumstrata(rr,strata,nstrata))
  ## S0 after strataC
  S00C = c(revcumsumstrata(rr,xxstrataC,nstrataC))

  ## censoring MG, strataC
  stratJumps <- dd$strata[jumps]
  S00i <- rep(0,length(dd$strata))
  S00i[jumps] <-  1/S00C[jumps]

  ## cif calculation, uses strata {{{
  S0Di <- S02i <- S01i <- rep(0,length(dd$strata))
  S01i[jumps1] <-  1/S0[jumps1]
  S02i[other] <-  1/S0[other]
  S0Di[jumpsD] <-  1/S0[jumpsD]

  ## strata-Cif G_T(t)
  if (!km) { 
     cumhazD <- cumsumstratasum(S0Di,xxstrata,nstrata)
     Stm <- exp(-cumhazD$lagsum)
     St <- exp(-cumhazD$sum)
  } else { 
     StA <- cumsumstratasum(log(1-S0Di),xxstrata,nstrata)
     Stm <- exp(StA$lagsum)
     St <- exp(StA$sum)
  }

  ## G_c(t-) 
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S00i,xxstrataC,nstrataC)$lagsum)
    Gc      <- exp(-cumhazD)
  } else Gc <- c(exp(cumsumstratasum(log(1-S00i),xxstrataC,nstrataC)$lagsum))
 cif1 <- cumsumstrata(Stm*S01i,xxstrata,nstrata)
 cif2 <- cumsumstrata(Stm*S02i,xxstrata,nstrata)
# }}}

 Et <- matrix(0,nrow(Z),ncol(Z))
 Et[jumps1,] <- E

 Lam1fg <- -log(1-cif1)
 laststrata <- tailstrata(xxstrata,nstrata)
 gtstart <- Lam1fg[laststrata]
 dLam1fg <- c(diffstrata(Lam1fg,xxstrata,nstrata))

 ## Gc~strataC, dLam1fg~strata
 tailcstrata <- tailstrata(xxstrataC,nstrataC)
 Gcstart <- Gc[tailcstrata]

 dstrata <- mystrata(data.frame(cbind(xxstrataC,xxstrata)))
 ndstrata <- attr(dstrata,"nlevel")
 lastt <- tailstrata(dstrata-1,ndstrata)

 ### ## \int_t^\infty G_c^j(t) d\Lambda_1^k(t)
 G0start <- rep(1,nstrataC)
 cLam1fg  <- cumsum2strata(Gc,dLam1fg,xxstrataC,nstrataC,xxstrata,nstrata,G0start)$res
 lastt <- tailstrata(dstrata-1,ndstrata)
 RLam1fg <- cLam1fg[lastt][dstrata]-cLam1fg

 ## E(s) from FG without strata
 ## \int_0^t  G_c^j(s) E(s) d\Lambda_1^k(s)
 fff <- function(x) {
   cx  <- cumsum2strata(Gc,x*dLam1fg,xxstrataC,nstrataC,xxstrata,nstrata,G0start)$res
    return(cx[lastt][dstrata]-cx)
 }
 ERLam1fg0  <- apply(Et,2,fff)

 gt <-  RLam1fg*c(cif2/(Gc*St))
 gt[gt==Inf] <- 0
 gt[is.na(gt)] <- 0
 ERLam1fg  <- ERLam1fg0*c(cif2/(Gc*St))
 ERLam1fg[ERLam1fg==Inf] <- 0
 ERLam1fg[is.na(ERLam1fg)] <- 0

 sss <- headstrata(dstrata-1,ndstrata)
 gtstart <- gt[sss]
 E1dLam0 <- cumsum2strata(gt,S00i,dstrata-1,ndstrata,xxstrataC,nstrataC,gtstart)$res

 fff <- function(x) {
    gtstart <- x[sss]
    cx  <- cumsum2strata(x,S00i,dstrata-1,ndstrata,xxstrataC,nstrataC,gtstart)$res
    return(cx)
 }
 E2dLam0 <- apply(ERLam1fg,2,fff)

 U1 <- matrix(0,nrow(Z),1)
 U2 <- matrix(0,nrow(Z),ncol(Z))
 U1[jumps,] <- gt[jumps]
 U2[jumps,] <- ERLam1fg[jumps,]

 ### Martingale  as a function of time and for all subjects to handle strata 
 MG1t <- Z*c(U1[,,drop=FALSE]-E1dLam0)*rr*c(dd$weights)
 MG2t <- (U2[,,drop=FALSE]-E2dLam0)*rr*c(dd$weights)
 MGt <- MG1t-MG2t
 MGiid <- apply(MGt,2,sumstrata,dd$id,max(dd$id)+1)
 augment <- apply(MGt,2,sum)

 augment <- list(MGiid=MGiid,augment=augment,id=id,id.orig=id.orig,
###  dd=dd,div=div,Et=Et,E=E, gt=gt, ERLam1fg=ERLam1fg,RLam1fg=RLam1fg,ERLam1fg0=ERLam1fg0,
     jumps1=jumps1,jumps=jumps,other=other,
     nstrata=nstrata,nstrataC=nstrataC,dstrata=dstrata,ndstrata=ndstrata,
     cif1=cif1,cif2=cif2,St=St,Gc=Gc,strata=xxstrata,strataC=xxstrataC,time=dd$time)

 ## drop strata's from formula and run wiht augmention term
 # {{{
 drop.strata <- function(x) {
   mm <- unlist(Specials(x,"strata"))
   for (i in mm) x <- update(x, paste(".~.-strata(",i,")"))
   mm <- unlist(Specials(x,"strataC"))
   for (i in mm) x <- update(x, paste(".~.-strataC(",i,")"))
   return(x)
 }

 formulans <- drop.strata(formula)
# }}}

  if (nstrataC==1) cens.model <- ~+1 else cens.model <- ~strata(strataCC)
  data$strataCC <- cens.strata

 fga <- cifreg(formulans,data=data,cause=cause,
     propodds=NULL,augmentation=augment$augment,cens.model=cens.model,...)

 ## adjust SE and var based on augmentation term
 fga$var.orig <- fga$var
 fga$augment <- augment$augment
 fga$iid <- fga$iid + MGiid %*% fga$ihessian 
 fga$var <- crossprod(fga$iid)
 fga$se.coef <-  diag(fga$var)^.5
 fga$MGciid <- MGiid

 return(fga)
}# }}})

