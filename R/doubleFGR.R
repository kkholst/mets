##' Double CIF Fine-Gray model with two causes
##'
##' Estimation based on derived hazards and recursive esitmating euqations.
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

##' Augmentation for Fine-Gray model with two causes
##'
##' @examples
##' rho1 <- 0.1; rho2 <- 0.9
##' n <- 200
##' beta=c(0.3,-0.3,0.1,0.1)
##' dats <- simul.cifs(n,rho1,rho2,beta)
##' dsort(dats) <- ~time
##' fgcm <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL,cens.model=~strata(Z1))
##' 
##' fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
##' cr <- phreg(Surv(time,status==0)~+1,data=dats)
##' dtable(dats,~status)
##' dfg <- doubleFGR(Event(time,status)~Z1+Z2,data=dats,restrict=3)
##' fgaug <- augmentationFG(dats,fg,fgcm,dfg,cr)
##' 
##' fgaug2 <- augmentationFG(dats,fgaug,fgcm,dfg,cr)
##' 
##' @aliases simul.cifs 
##' @export
augmentationFG <- function(dats,fg,fgcm,dfg,cr,cens.model=FALSE,status="status") 
{# {{{
## dats sorted after time

if (!is.null(dfg)) {

timec <- cr$cumhaz[,1]
times1 <- fgcm$cumhaz[,1]
###
dup <- duplicated(timec)
censrows <- which(dats[,status]==0)[!dup]
ll <- length(censrows)
censrows <- censrows[-ll]
Gtc <- exp(-cr$cumhaz[,2])[!dup]
Gtc <- Gtc[-ll]
###
timec <- timec[!dup]
timec <- timec[-ll]

### take E from fgcm to get correct E
pcif1 <- predictdFG(dfg,newdata=dats,cause=1,times=c(0,times1),se=FALSE)
pcif2 <- predictdFG(dfg,newdata=dats,cause=2,times=c(0,times1),se=FALSE)
pGct1  <-  predict(cr,newdata=dats,times=c(times1),se=FALSE)$surv
cum1 <- pcif1$cumhaz
cif1 <- pcif1$cif
dcum <- t(apply(pcif1$cumhaz,1,diff))
nc <- ncol(dcum)
###########################################################
### H(t) = \int_t^infty e(t) G_c(t) d \Lambda_1(t,X)    ###
### H(t) = \int_t^infty  G_c(t) d \Lambda_1(t,X)        ###
### evaluated first at times1 and then at timec 
### we can compute this  as
###     exp(X_i^T beta)  \int_t^infty  G_c(t) d \Lambda_1(t) 
###     exp(X_i^T beta)  \int_t^infty  E(t) G_c(t) d \Lambda_1(t) 
### so only need the two integrals once ! 
###########################################################
whereH <- sindex.prodlim(c(0,times1),timec)

intAug <- function(dcum,pGct1,E,whereH) {# {{{
	HGt<-    t(apply(dcum*pGct1,1,revcumsum))
	He1Gt <- t(apply(dcum*pGct1*E[,1],1,revcumsum))
	He2Gt <- t(apply(dcum*pGct1*E[,2],1,revcumsum))
	nc <- ncol(HGt)
        ##### \int_s^infty      G_c(t)  d \Lambda_1(t,X)
        ###### \int_s^infty E(t) G_c(t)  d \Lambda_1(t,X)
	### do not start at 0 but  and evaluated at timec
	dHGt <-         cbind(HGt[,1],HGt)[,whereH]     
	dHe1Gt <-   cbind(He1Gt[,1],He1Gt)[,whereH]
	dHe2Gt <-   cbind(He2Gt[,1],He2Gt)[,whereH] 

return(list(dHGt=dHGt,dHe1Gt=dHe1Gt,dHe2Gt=dHe2Gt))
}# }}}

Hfg   <- intAug(dcum,pGct1,fg$E,whereH)
Hfg$timec <- timec
Hfg$dcum <- dcum
Hfg$Gct1 <- pGct1 
Hfg$whereH <- whereH

F1t <- pcif1$cif; F2t <- pcif2$cif; GTt <- 1-F1t-F2t
F1t <- F1t[,whereH]; F2t <- F2t[,whereH]; GTt <- GTt[,whereH]; 
pGct1  <-  predict(cr,newdata=dats,times=timec,tminus=TRUE,se=FALSE)$surv
Xalls <- as.matrix(fg$model.frame[,-1])
###RRall <- c(exp( Xalls %*% ccr$coef) )
RR1   <- rep(1,nrow(Xalls))

intdMc<-function(dHGt,dHe1,dHe2,Gtt,RRall,Xalls,timec,GTt,F1t,F2t,censrows)
{# {{{
aug <- matrix(0,length(censrows),ncol(Xalls))
ww <- 1
for (ww in seq_along(timec)) {
i <- censrows[ww]
xi <- Xalls[i,]; GTtx <- GTt[i,ww]; F2tx <-  F2t[i,ww]; 
Gct <- Gtt[i,ww]; 
dHe <-    c(dHe1[i,ww],dHe2[i,ww]) 
aug1 <-  (F2tx/GTtx)*(xi*dHGt[i,ww]-dHe)*(1/Gct)
#### compensator part 
takei <- i:nrow(Xalls)
F2dGt <- F2t[takei,ww]/GTt[takei,ww]
Xs <-  Xalls[takei,,drop=FALSE]
RR <-  RRall[takei]; 
Gctt  <- Gtt[takei,ww]
avau21  <-  apply(F2dGt*dHGt[takei,ww]*Xs*RR/Gctt,2,sum)/sum(RR)
HtXsGt <- RR*(F2dGt/Gctt)*cbind(dHe1[takei,ww],dHe2[takei,ww]) 
avau22 <- apply(HtXsGt,2,sum)/sum(RR)
avau1 <- (avau21-avau22)
aug[ww,] <- aug1-avau1
}
augt <- aug
aug <- apply(augt,2,sum,na.rm=TRUE)
return(list(augt=augt,aug=aug,maug=aug/nrow(Xalls)))
}# }}}

aug <- intdMc(Hfg$dHGt,Hfg$dHe1Gt,Hfg$dHe2Gt,pGct1,RR1,Xalls,timec,GTt,F1t,F2t,censrows)

### agumentation using Gc for augmentation term, using E
fgas <- cifreg(fg$formula,data=dats,cause=1,propodds=NULL,beta=fg$coef,augmentation=aug$aug)

} else {
   fgas <- fgasX <- fgascm <-  data.frame(coef=rep(NA,2))
   aug <- augX  <-  augXX <- data.frame(aug=rep(NA,2))
   S0aug <- 0; S1aug <- 0
}

fgas$aug <- aug$aug
fgas$maug <- aug$maug
res= fgas
res$Hfg <- Hfg

###	  fgdrE=fgdrE$coef, ## augX=augX, fgascm=fgascm$coef, 

return(res)
}# }}}


##' @export
simul.cifs <- function(n,rho1,rho2,beta,rc=0.5,depcens=0,rcZ=0.5,bin=0) {# {{{
p=length(beta)/2
tt <- seq(0,6,by=0.1)
Lam1 <- rho1*(1-exp(-tt))
Lam2 <- rho2*(1-exp(-tt))

if (bin==0) Z=cbind(2*rbinom(n,1,1/2)-1,rnorm(n))
else Z=cbind(2*rbinom(n,1,1/2)-1,rbinom(n,1,1/2))
colnames(Z) <- paste("Z",1:2,sep="")
cif1 <- setup.cif(cbind(tt,Lam1),beta[1:2],Znames=colnames(Z),type="cloglog")
cif2 <- setup.cif(cbind(tt,Lam2),beta[3:4],Znames=colnames(Z),type="cloglog")
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


