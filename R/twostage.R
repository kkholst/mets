##' @title Twostage survival model for multivariate survival data
##'
##' @description
##' Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data
##' using marginals that are on Cox form. The dependence can be modelled via
##' \enumerate{
##' \item  Regression design on dependence parameter.
##' \item  Random effects, additive gamma model.
##' }
##'
##' If clusters contain more than two subjects, we use a composite likelihood
##' based on the pairwise bivariate models, for full MLE see twostageMLE.
##'
##' The two-stage model is constructed such that
##' given the gamma distributed random effects it is assumed that the survival functions
##' are indpendent, and that the marginal survival functions are on Cox form (or additive form)
##' \deqn{
##' P(T > t| x) = S(t|x)= exp( -exp(x^T \beta) A_0(t) )
##' }
##'
##' One possibility is to model the variance within clusters via a regression design, and
##' then one can specify a regression structure for the independent gamma distributed
##' random effect for each cluster, such that the variance is given by
##' \deqn{
##'  \theta = h( z_j^T \alpha)
##' }
##' where \eqn{z} is specified by theta.des, and a possible link function var.link=1 will
##' will use the exponential link \eqn{h(x)=exp(x)}, and var.link=0 the identity link \eqn{h(x)=x}.
##' The reported standard errors are based on the estimated information from the
##' likelihood assuming that the marginals are known (unlike the twostageMLE and for the
##' additive gamma model below).
##'
##' Can also fit a structured additive gamma random effects model, such
##' as the ACE, ADE model for survival data.  In this case the
##' random.design specificies the random effects for each subject within a cluster. This is
##' a matrix of 1's and 0's with dimension n x d.  With d random effects.
##' For a cluster with two subjects, we let the random.design rows be
##'  \eqn{v_1} and \eqn{v_2}.
##' Such that the random effects for subject
##' 1 is \deqn{v_1^T (Z_1,...,Z_d)}, for d random effects. Each random effect
##' has an associated parameter \eqn{(\lambda_1,...,\lambda_d)}.
##' By construction subjects 1's random effect are Gamma distributed with
##' mean \eqn{\lambda_j/v_1^T \lambda}
##' and variance \eqn{\lambda_j/(v_1^T \lambda)^2}. Note that the random effect
##' \eqn{v_1^T (Z_1,...,Z_d)} has mean 1 and variance \eqn{1/(v_1^T \lambda)}.
##' It is here asssumed that  \eqn{lamtot=v_1^T \lambda} is fixed within clusters
##' as it would be for the ACE model below.
##'
##' Based on these parameters the relative contribution (the heritability, h) is
##' equivalent to  the expected values of the random effects: \eqn{\lambda_j/v_1^T \lambda}
##'
##' The DEFAULT parametrization (var.par=1) uses the variances of the random effecs
##' \deqn{
##' \theta_j  = \lambda_j/(v_1^T \lambda)^2
##' }
##' For alternative parametrizations one can specify how the parameters relate to \eqn{\lambda_j}
##' with the argument var.par=0.
##'
##' For both types of models the basic model assumptions are that
##' given the random effects of the clusters the survival distributions within a cluster
##' are independent and ' on the form
##' \deqn{
##' P(T > t| x,z) = exp( -Z \cdot Laplace^{-1}(lamtot,lamtot,S(t|x)) )
##' }
##' with the inverse laplace of the gamma distribution with mean 1 and variance 1/lamtot.
##'
##' The parameters \eqn{(\lambda_1,...,\lambda_d)} are related to the parameters of the model
##' by a regression construction \eqn{pard} (d x k), that links the \eqn{d}
##' \eqn{\lambda} parameters
##' with the (k) underlying \eqn{\theta} parameters
##' \deqn{
##' \lambda = theta.des  \theta
##' }
##' here using theta.des to specify these low-dimension association. Default is a diagonal matrix.
##' This can be used to make structural assumptions about the variances of the random-effects
##' as is needed for the ACE model for example.
##'
##' The case.control option that can be used with the pair specification of the pairwise parts
##' of the estimating equations. Here it is assumed that the second subject of each pair is the proband.
##'
##' @references
##'
##' Twostage estimation of additive gamma frailty models for survival data.
##' Scheike (2019), work in progress
##'
##' Shih and Louis (1995) Inference on the association parameter in copula models for bivariate
##' survival data, Biometrics, (1995).
##'
##' Glidden (2000), A Two-Stage estimator of the dependence
##' parameter for the Clayton Oakes model, LIDA, (2000).
##'
##' Measuring early or late dependence for bivariate twin data
##' Scheike, Holst, Hjelmborg (2015), LIDA
##'
##' Estimating heritability for cause specific mortality based on twins studies
##' Scheike, Holst, Hjelmborg (2014), LIDA
##'
##' Additive Gamma frailty models for competing risks data, Biometrics (2015)
##' Eriksson and Scheike (2015),
##'
##' @examples
##' data(diabetes)
##'
##' # Marginal Cox model  with treat as covariate
##' margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
##' ### Clayton-Oakes, MLE
##' fitco1<-twostageMLE(margph,data=diabetes,theta=1.0)
##' summary(fitco1)
##'
##' ### Plackett model
##' mph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
##' fitp <- survival.twostage(mph,data=diabetes,theta=3.0,Nit=40,
##'                clusters=diabetes$id,var.link=1,model="plackett")
##' summary(fitp)
##'
##' ### Clayton-Oakes
##' fitco2 <- survival.twostage(mph,data=diabetes,theta=0.0,detail=0,
##'                  clusters=diabetes$id,var.link=1,model="clayton.oakes")
##' summary(fitco2)
##' fitco3 <- survival.twostage(margph,data=diabetes,theta=1.0,detail=0,
##'                  clusters=diabetes$id,var.link=0,model="clayton.oakes")
##' summary(fitco3)
##'
##' ### without covariates but with stratafied
##' marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
##' fitpa <- survival.twostage(marg,data=diabetes,theta=1.0,
##'                 clusters=diabetes$id,model="clayton.oakes")
##' summary(fitpa)
##'
##' fitcoa <- survival.twostage(marg,data=diabetes,theta=1.0,clusters=diabetes$id,
##'                  model="clayton.oakes")
##' summary(fitcoa)
##'
##' ### Piecewise constant cross hazards ratio modelling
##' ########################################################
##'
##' d <- subset(simClaytonOakes(2000,2,0.5,0,stoptime=2,left=0),!truncated)
##' udp <- piecewise.twostage(c(0,0.5,2),data=d,method="optimize",
##'                           id="cluster",timevar="time",
##'                           status="status",model="clayton.oakes",silent=0)
##' summary(udp)
##'
##' \donttest{ ## Reduce Ex.Timings
##' ### Same model using the strata option, a bit slower
##' ########################################################
##' ## makes the survival pieces for different areas in the plane
##' ##ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
##' ##ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
##' ##ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
##' ##ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")
##'
##' ## everything done in one call
##' ud <- piecewise.data(c(0,0.5,2),data=d,timevar="time",status="status",id="cluster")
##' ud$strata <- factor(ud$strata);
##' ud$intstrata <- factor(ud$intstrata)
##'
##' ## makes strata specific id variable to identify pairs within strata
##' ## se's computed based on the id variable across strata "cluster"
##' ud$idstrata <- ud$id+(as.numeric(ud$strata)-1)*2000
##'
##' marg2 <- timereg::aalen(Surv(boxtime,status)~-1+factor(num):factor(intstrata),
##'                data=ud,n.sim=0,robust=0)
##' tdes <- model.matrix(~-1+factor(strata),data=ud)
##' fitp2 <- survival.twostage(marg2,data=ud,se.clusters=ud$cluster,clusters=ud$idstrata,
##'                 model="clayton.oakes",theta.des=tdes,step=0.5)
##' summary(fitp2)
##'
##' ### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
##' ud$stratas <- ud$strata;
##' ud$stratas[ud$strata=="0.5-2,0-0.5"] <- "0-0.5,0.5-2"
##' tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
##' fitp3 <- survival.twostage(marg2,data=ud,clusters=ud$idstrata,se.cluster=ud$cluster,
##'                 model="clayton.oakes",theta.des=tdes2,step=0.5)
##' summary(fitp3)
##'
##' ### same model using strata option, a bit slower
##' fitp4 <- survival.twostage(marg2,data=ud,clusters=ud$cluster,se.cluster=ud$cluster,
##'                 model="clayton.oakes",theta.des=tdes2,step=0.5,strata=ud$strata)
##' summary(fitp4)
##' }
##'
##' \donttest{ ## Reduce Ex.Timings
##' ### structured random effects model additive gamma ACE
##' ### simulate structured two-stage additive gamma ACE model
##' data <- simClaytonOakes.twin.ace(4000,2,1,0,3)
##' out <- twin.polygen.design(data,id="cluster")
##' pardes <- out$pardes
##' pardes
##' des.rv <- out$des.rv
##' head(des.rv)
##' aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data,robust=0)
##' ts <- survival.twostage(aa,data=data,clusters=data$cluster,detail=0,
##' 	       theta=c(2,1),var.link=0,step=0.5,
##' 	       random.design=des.rv,theta.des=pardes)
##' summary(ts)
##' }
##'
##' @keywords survival
##' @author Thomas Scheike
##' @param margsurv Marginal model
##' @param data data frame
##' @param method Scoring method "nr", for lava NR optimizer
##' @param detail Detail
##' @param clusters Cluster variable
##' @param silent Debug information
##' @param weights Weights
##' @param theta Starting values for variance components
##' @param theta.des design for dependence parameters, when pairs are given the indeces of the
##' theta-design for this pair, is given in pairs as column 5
##' @param var.link Link function for variance:  exp-link.
##' @param baseline.iid to adjust for baseline estimation, using phreg function on same data.
##' @param model model
##' @param marginal.trunc marginal left truncation probabilities
##' @param marginal.survival optional vector of marginal survival probabilities
##' @param strata strata for fitting, see example
##' @param se.clusters for clusters for se calculation with iid
##' @param numDeriv to get numDeriv version of second derivative, otherwise uses sum of squared scores for each pair
##' @param random.design random effect design for additive gamma model, when pairs are given the
##' indeces of the pairs random.design rows are given as columns 3:4
##' @param pairs matrix with rows of indeces (two-columns) for the pairs considered in the pairwise
##' composite score, useful for case-control sampling when marginal is known.
##' @param dim.theta dimension of the theta parameter for pairs situation.
##' @param numDeriv.method uses simple to speed up things and second derivative not so important.
##' @param additive.gamma.sum for two.stage=0, this is specification of the lamtot in the models via
##' a matrix that is multiplied onto the parameters theta (dimensions=(number random effects x number
##' of theta parameters), when null then sums all parameters.
##' @param var.par is 1 for the default parametrization with the variances of the random effects,
##' var.par=0 specifies that the \eqn{\lambda_j}'s are used as parameters.
##' @param no.opt for not optimizng
##' @param ... Additional arguments to maximizer NR of lava.
##' and ascertained sampling
##' @aliases survival.twostage twostage.aalen twostage.cox.aalen twostage.coxph twostage.phreg randomDes readmargsurv 
##' @export survival.twostage
survival.twostage <- function(margsurv,data=parent.frame(),
    method="nr",detail=0,clusters=NULL,
    silent=1,weights=NULL,theta=NULL,theta.des=NULL,
    var.link=1,baseline.iid=1,model="clayton.oakes",
    marginal.trunc=NULL,marginal.survival=NULL,strata=NULL,
    se.clusters=NULL,numDeriv=1,random.design=NULL,pairs=NULL,dim.theta=NULL,
    numDeriv.method="simple",additive.gamma.sum=NULL,var.par=1,no.opt=FALSE,...)
{#{{{
  requireNamespace("timereg")
# seting up design and variables
iid <- 1
two.stage <- 1; rate.sim <- 1; sym=1; var.func <- NULL
if (model=="clayton.oakes" || model=="gamma") dep.model <- 1
else if (model=="plackett" || model=="or") dep.model <- 2
else stop("Model must by either clayton.oakes or plackett \n");
start.time <- NULL; ptrunc <- NULL; psurvmarg <- NULL; status <- NULL; score.iid <- NULL
fix.baseline <- 0; convergence.bp <- 1;  ### to control if baseline profiler converges

  if ((!is.null(margsurv)) | (!is.null(marginal.survival))) fix.baseline <- 1
  antpers <- nrow(data); RR <-  rep(1,antpers);


  if (!is.null(margsurv)) {
     rrr <-  readmargsurv(margsurv,data,clusters)
     psurvmarg <- rrr$psurvmarg; ptrunc <- rrr$ptrunc; start.time <- rrr$entry;
     time2 <- rrr$exit; status <- rrr$status; clusters <- rrr$clusters
  }

  if (is.null(psurvmarg)) psurvmarg <- rep(1,antpers);
  if (!is.null(marginal.survival)) psurvmarg <- marginal.survival
  if (!is.null(marginal.trunc)) ptrunc <- marginal.trunc
  if (is.null(ptrunc)) ptrunc <- rep(1,length(psurvmarg))
  if (is.null(weights)==TRUE) weights <- rep(1,antpers);
  if (is.null(strata)==TRUE) strata<- rep(1,antpers);
  if (length(strata)!=antpers) stop("Strata must have length equal to number of data points \n");

  # cluster set up
  cluster.call <- clusters
  out.clust <- cluster.index(clusters);
  clusters <- out.clust$clusters
  maxclust <- out.clust$maxclust
  antclust <- out.clust$cluster.size
  clusterindex <- out.clust$idclust
  clustsize <- out.clust$cluster.size
  call.secluster <- se.clusters

  if (is.null(se.clusters)) { se.clusters <- clusters; antiid <- nrow(clusterindex);} else  {
      iids <-  unique(se.clusters);
      antiid <- length(iids);
      if (is.numeric(se.clusters)) se.clusters <-  fast.approx(iids,se.clusters)-1
       else se.clusters <- as.integer(factor(se.clusters, labels = seq(antiid)))-1
  }
  if (length(se.clusters)!=length(clusters)) stop("Length of seclusters and clusters must be same\n");

  #

  ### setting design for random variables, in particular with pairs are given
  ddd <- randomDes(dep.model,random.design,theta.des,theta,antpers,additive.gamma.sum,pairs,var.link,clusterindex,dim.theta)
  random.design=ddd$random.design;clusterindex=ddd$clusterindex;
  antpairs=ddd$antpairs;
  theta=ddd$theta;ptheta=ddd$ptheta;theta.des=ddd$theta.des
  pair.structure=ddd$pair.structure; dep.model=ddd$dep.model
  dim.rv <- ddd$dim.rv; additive.gamma.sum=ddd$additive.gamma.sum

  theta.score<-rep(0,ptheta);
  Stheta<-var.theta<-matrix(0,ptheta,ptheta);
  #
  ascertained <- 0

  obj <- function(par)
  { #

     if (pair.structure==0 | dep.model!=3) Xtheta <- as.matrix(theta.des) %*% matrix(c(par),nrow=ptheta,ncol=1);
     if (pair.structure==1 & dep.model==3) Xtheta <- matrix(0,antpers,1); ## not needed

      if (var.link==1 & dep.model==3) epar <- c(exp(par)) else epar <- c(par)
      partheta <- epar

      if (var.par==1 & dep.model==3) {
	 ## from variances to
	 if (is.null(var.func)) {
	    sp <- sum(epar)
	    partheta <- epar/sp^2
         } else partheta <- epar ## par.func(epar)
      }

      if (pair.structure==0) {
	      outl<-.Call("twostageloglikeRV", # only two stage model for this option
	      icause=status,ipmargsurv=psurvmarg,
	      itheta=c(partheta), ithetades=theta.des,
	      icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
	      ivarlink=var.link,iid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
	      itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
	      irvdes=random.design,iags=additive.gamma.sum,iascertained=ascertained,
              PACKAGE="mets") #
      }
      else { # pair-structure
          outl<-.Call("twostageloglikeRVpairs", #
          icause=status,ipmargsurv=psurvmarg, itheta=c(partheta), ithetades=theta.des,
          icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex, ivarlink=var.link,
	  iiid=iid,iweights=weights,isilent=silent, idepmodel=dep.model,
	  itrunkp=ptrunc,istrata=as.numeric(strata), iseclusters=se.clusters,iantiid=antiid,
	  irvdes=random.design, iags=additive.gamma.sum, iascertained=ascertained,PACKAGE="mets")
	  #
      } #

    if (detail==3) print(c(partheta,outl$loglike))

    ## variance parametrization, and inverse.link
    if (dep.model==3) {
    if (var.par==1) {
         ## from variances to and with sum for all random effects
         if (is.null(var.func)) {
	 if (var.link==0)  {
             mm <- matrix(-epar*2*sp,length(epar),length(epar))
	     diag(mm) <- sp^2-epar*2*sp
	 } else {
            mm <- -c(epar) %o% c(epar)*2*sp
            diag(mm) <- epar*sp^2-epar^2*2*sp
	 }
	    mm <- mm/sp^4
	 } else mm  <- numDeriv::hessian(var.func,par)
      } else {
	   if (var.link==0) mm <- diag(length(epar)) else
			  mm <- diag(length(c(epar)))*c(epar)
      }
      }


    if (dep.model==3) {
       outl$score <-  t(mm) %*% outl$score
       outl$Dscore <- t(mm) %*% outl$Dscore %*% mm
       if (iid==1) { outl$score.iid <- t(t(mm) %*% t(outl$score.iid))
               if (inherits(margsurv,"phreg")) {
	          outl$D1thetal  <- t(t(mm) %*% t(outl$D1thetal))
	          outl$D2thetal  <- t(t(mm) %*% t(outl$D2thetal))
	       }
       }
    }


    outl$gradient <- outl$score 
    outl$hessian <- outl$Dscore

    if (oout==0) ret <- with(outl,structure(outl$loglike,gradient=-gradient,hessian=-hessian))
    else if (oout==1) ret <- outl$gradient 
    else if (oout==2) ret <- outl

    return(ret)
  } #

  score1 <- NULL
  theta.iid <- NULL
  p <- theta

  oout <- 0
  opt <- NULL
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
          tim <- system.time(opt <- lava::NR(p,obj,...))
          opt$timing <- tim
          opt$estimate <- opt$par
      } else {
          opt <- nlm(obj,beta,...)
          opt$method <- "nlm"
      }
      cc <- opt$estimate;  ## names(cc) <- colnames(X)
      oout <- 2
      val <- c(list(coef=cc),obj(opt$estimate))
 } else {oout <- 2; val <- c(list(coef=p),obj(p))}

     if (numDeriv>=1) {
	 p <- val$coef
         oout <- 1
         if (detail==1 ) cat("starting numDeriv for second derivative \n");
         val$hessian <- numDeriv::jacobian(obj,p,method=numDeriv.method)
         if (detail==1 ) cat("finished numDeriv for second derivative \n");
        }


    ### update p, but note that score and derivative in fact related to previous p
    if (!is.nan(sum(p))) {
      if (detail==1 && iid==1) cat("iid decomposition\n");
         if (iid==1) {  #
            score.iid <- val$score.iid
	    val$theta.iid <- score.iid
	    theta.iid <- score.iid
	    ucc <-  unique(cluster.call)
	    if (length(ucc)== nrow(theta.iid))
		    rownames(theta.iid) <- unique(cluster.call)
                ### lets iid for theta be just score to start, correction for marginal for phreg call

                if (inherits(margsurv,"phreg") & baseline.iid==1) {  # adjust for baseline when phreg is used

                  if ((margsurv$no.opt) | is.null(margsurv$coef)) fixbeta<- 1 else fixbeta <- 0
		  xx <- margsurv$cox.prep
                  id  <-  xx$id+1
		  cumhazt <- c(rrr$cum)
		  rr <- c(rrr$RR)

		  ## these refers to id given in cluster.call
		  xx <- margsurv$cox.prep
		  D1thetal<- val$D1thetal
		  D2thetal<- val$D2thetal
		  Dlamthetal <- -(D1thetal+D2thetal)
		  ## ordered as original data
                  Fbeta <- t(Dlamthetal) %*% (margsurv$X *c(cumhazt*psurvmarg))

                  ### timeordered
		  ## t- not needed because using S(T-) for survival already
                  Gbasedt <- Dlamthetal[xx$ord+1,,drop=FALSE]*psurvmarg[xx$ord+1]*rr[xx$ord+1]
                  Gbase <- apply(Gbasedt,2,revcumsumstrata,xx$strata,xx$nstrata)

                  if (fixbeta==0) {  ### iid after beta of marginal model
		     Z <- xx$X
		     U <- E <- matrix(0,nrow(xx$X),margsurv$p)
		     E[xx$jumps+1,] <- margsurv$E
		     U[xx$jumps+1,] <- margsurv$U
		     invhess <- -solve(margsurv$hessian)
		     S0i <- rep(0,length(xx$strata))
		     S0i[xx$jumps+1] <- 1/margsurv$S0
		     cumhaz <- c(cumsumstrata(S0i,xx$strata,xx$nstrata))
		     EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	             rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset))
		     ### Martingale for all subjects
		     MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx$weights)
		     UUbeta <- apply(MGt,2,sumstrata,id-1,max(id))
		     UUbeta <- UUbeta %*% invhess
		     GbaseEdLam0 <- t(Gbase) %*% (E*S0i)
		     Fbeta   <-  Fbeta -  GbaseEdLam0
		     Fbetaiid <- UUbeta %*% t(Fbeta)
	          }

	  ###  \int_0^\tau (GBase(s) / S_0(s)) dN_i(s) - dLamba_i(s)
	  xx <- margsurv$cox.prep
	  S0i  <- S0i2 <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/margsurv$S0
	  S0i2[xx$jumps+1] <- 1/margsurv$S0^2

	  GbasedN <- Gbase*S0i
	  GbasedLam0 <- apply(Gbase*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)

	  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset)) else rr <- c(xx$sign*exp( xx$offset))
	  MGt <- (GbasedN[,drop=FALSE]-GbasedLam0*rr)*c(xx$weights)
	  MGti <- apply(MGt,2,sumstrata,id-1,max(id))

	  if (fixbeta==0) MGti <-  MGti + Fbetaiid
	  theta.iid <-  -theta.iid + MGti
	  val$theta.iid <- theta.iid
	  } #
   } #
   }

  hess <- val$hessian
  if (!is.na(sum(hess))) hessi <- lava::Inverse(val$hessian) else hessi <- diag(nrow(val$hessian))

# handling output
  loglikeiid <- NULL; robvar.theta <- NULL; likepairs <- NULL; marginal.surv <- psurvmarg; 
  marginal.trunc <- ptrunc; 

  if (iid==1) {
  if (dep.model==3 & pair.structure==1) likepairs <- val$likepairs
  if (dep.model==3 & two.stage==0) {
	  hessi <- -1*hessi
	  all.likepairs <- val$all.likepairs
	  colnames(all.likepairs) <- c("surv","dt","ds","dtds","cause1","cause2")
  }
     theta.iid <- val$theta.iid %*% hessi
### rownames not set to make more robust
### if (is.null(call.secluster)) rownames(theta.iid) <- unique(cluster.call) else rownames(theta.iid) <- unique(se.clusters)
     robvar.theta  <- crossprod(theta.iid)
     loglikeiid <- val$loglikeiid
    } else { all.likepairs <- NULL}
  var.theta <- robvar.theta
  if (is.null(robvar.theta)) var.theta <- hessi

  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:length(theta),sep="")
  theta <- matrix(theta,length(c(theta)),1)
  if (length(thetanames)==nrow(theta)) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }


  ud <- list(theta=matrix(val$coef,ncol=1),coef=val$coef,score=val$gradient,hess=hess,hessi=hessi,var.theta=var.theta,
     model=model,robvar.theta=robvar.theta,
     loglike=val$loglike,loglikeiid=loglikeiid,likepairs=likepairs,
     theta.iid=theta.iid,thetanames=thetanames,
     score1=score1,Dscore=val$Dscore,
     marginal.surv=marginal.surv,marginal.trunc=marginal.trunc,
     se=diag(robvar.theta)^.5,score.iid=-score.iid,theta.des=theta.des,random.design=random.design)
  class(ud) <- "mets.twostage"
  attr(ud,"response") <- "survival"
  attr(ud,"Formula") <- formula
  attr(ud,"clusters") <- clusters
  attr(ud,"cluster.call") <- cluster.call
  attr(ud,"secluster") <- c(se.clusters)
  attr(ud,"sym")<-sym;
  attr(ud,"var.link")<-var.link;
  attr(ud,"var.par")<-var.par;
  attr(ud,"var.func")<-var.func;
  attr(ud,"ptheta")<-ptheta
  attr(ud,"antpers")<-antpers;
  attr(ud,"antclust")<-antclust;
  attr(ud,"Type") <- model
  attr(ud,"twostage") <- two.stage
  attr(ud,"additive-gamma") <- (dep.model==3)*1
  if (!is.null(marginal.trunc)) attr(ud,"trunclikeiid")<- val$trunclikeiid
  if (dep.model==3 & two.stage==0) attr(ud,"all.likepairs")<- all.likepairs
  if (dep.model==3 ) attr(ud,"additive.gamma.sum") <- additive.gamma.sum

if (dep.model==3 & pair.structure==0) {
       attr(ud, "pardes") <- theta.des
       attr(ud, "theta.des") <- theta.des
       attr(ud, "rv1") <-    random.design[1,]
    }
    if (dep.model==3 & pair.structure==1) {
       nrv <- clusterindex[1,6]
       attr(ud, "theta.des") <- matrix(theta.des[1,1:(nrv*ptheta)],nrv,ptheta)
       attr(ud, "pardes") <-    matrix(theta.des[1,1:(nrv*ptheta)],nrv,ptheta)
       attr(ud, "rv1") <-    random.design[1,]
       attr(ud, "nrv") <- clusterindex[1,6]
    }


  return(ud);
  #

} # }}}

##' @export
randomDes <- function(dep.model,random.design,theta.des,theta,antpers,ags,pairs,var.link,clusterindex,dim.theta)
{
   additive.gamma.sum <- ags

  if (!is.null(random.design)) { ### different parameters for Additive random effects 
     dep.model <- 3
     dim.rv <- ncol(random.design);
     if (is.null(theta.des)) theta.des<-diag(dim.rv);
  } else {
      random.design <- matrix(0,1,1);  dim.rv <- 1;
      additive.gamma.sum <- matrix(1,1,1);
  }

  if (is.null(theta.des)) ptheta<-1;
  if (is.null(theta.des)) theta.des<-matrix(1,antpers,ptheta); ###  else theta.des<-as.matrix(theta.des);
  ptheta<-ncol(theta.des)

  if (is.matrix(pairs)) {
       if ( ncol(pairs)==6)  ptheta<- dim.theta else ptheta<-ncol(theta.des)
  }

  theta.des <- as.matrix(theta.des)
  if (is.null(theta)==TRUE) {
         if (var.link==1) theta<- rep(-0.7,ptheta);
         if (var.link==0) theta<- rep(exp(-0.7),ptheta);
  }

  if (length(theta)!=ptheta) {
         theta<-rep(theta[1],ptheta);
  }
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta);
  antpairs <- 1; ### to define

  if (is.null(additive.gamma.sum)) additive.gamma.sum <- matrix(1,dim.rv,ptheta)
  if (!is.null(pairs)) { pair.structure <- 1;} else  pair.structure <- 0;

  if (pair.structure==1 & dep.model==3) { #
  antpairs <- nrow(pairs);

  if ( (ncol(pairs)==2))
  {
      pairs <- cbind(pairs,pairs,1,ncol(random.design))
      theta.des <- matrix(c(theta.des),1,length(c(theta.des)))
  }
   clusterindex <- cbind(pairs[,1:5]-1,ncol(random.design))
  } #

  if (pair.structure==1 & dep.model!=3) { #
       antpairs <- nrow(pairs);
       if ((ncol(pairs)==2)) { ## when only pairs are given we refer to rows of theta.des from 3 column: same as index of pair-1
       pairs <- cbind(pairs,pairs[,1])
       }
       if (ncol(pairs)!=3) stop("with pairstructure and theta.des, 3rd column of pairs must give relevant design from theta.des\n")
       ## index in c are Rindex -1
       clusterindex <- pairs-1;
  }

###  print(head(pairs)); print(head(theta.des))
###  print(dim(pairs)); print(dim(theta.des))

  return(list(random.design=random.design,clusterindex=clusterindex,
	 antpairs=antpairs,pair.structure=pair.structure,
	 dep.model=dep.model,dim.rv=dim.rv, additive.gamma.sum=additive.gamma.sum,
         theta=theta,ptheta=ptheta,theta.des=theta.des))

}



##' @export
readmargsurv <- function(margsurv,data,clusters)
{
start.time <- 0

if (inherits(margsurv,c("aalen","cox.aalen"))) {
	 formula<-attr(margsurv,"Formula");
	 beta.fixed <- attr(margsurv,"beta.fixed")
	 if (is.null(beta.fixed)) beta.fixed <- 1;
	 ldata <- timereg::aalen.des(formula,data=data,model="cox.aalen");
	 id <- attr(margsurv,"id");
	 mclusters <- attr(margsurv,"cluster.call")
	 X<-ldata$X;
	 time<-ldata$time2;
	 Z<-ldata$Z;
	 status<-ldata$status;
	 time2 <- attr(margsurv,"stop");
	 start.time <- attr(margsurv,"start")
	 antpers<-nrow(X);
         if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else {
		     Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
         if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
	 px<-ncol(X);

         if (is.null(clusters) && is.null(mclusters)) stop("No cluster variabel specified in marginal or twostage call\n");
         if (is.null(clusters)) clusters <- mclusters
	 if (nrow(X)!=length(clusters)) stop("Length of Marginal survival data not consistent with cluster length\n"); 
	 } else if (inherits(margsurv,"phreg")) { 
            ### setting up newdata with factors and strata
                antpers <- nrow(data)
		 rr <- readPhreg(margsurv,data,nr=FALSE)
		 time2 <- rr$exit
		 if (!is.null(rr$entry)) start.time <- rr$entry
		 else start.time <- rep(0,antpers);
		 status <- rr$status
	         if (is.null(clusters)) clusters <- rr$clusters
	
} else { ### coxph 
	  antpers <- margsurv$n
	  id <- 0:(antpers-1)
	  mt <- model.frame(margsurv)
	  Y <- model.extract(mt, "response")
	  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
	  if (attr(Y, "type") == "right") {
	      time2 <- Y[, "time"];
	      status <- Y[, "status"]
		start.time <- rep(0,antpers);
		} else {
		 start.time <- Y[, 1];
		 time2 <- Y[, 2];
		 status <- Y[, 3];
		}
	   Z <- matrix(1,antpers,length(coef(margsurv)));

	   if (is.null(clusters)) stop("must give clusters for coxph\n");
	   cluster.call <- clusters;
	   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
	   px <- 1; pz <- ncol(Z);
	   start <- rep(0,antpers);
	   beta.fixed <- 0
	   semi <- 1
}

  if (!is.null(start.time)) {
      if (any(abs(start.time)>0)) lefttrunk <- 1  else lefttrunk <- 0;
  } else lefttrunk <- 0

if (!is.null(margsurv))  {
  if (inherits(margsurv,c("aalen","cox.aalen")))  { 
         resi <- timereg::residualsTimereg(margsurv,data=data)
         RR  <- resi$RR
	 cum <- resi$cumhaz/RR
	 psurvmarg <- exp(-resi$cumhaz);
         ptrunc <- rep(1,length(psurvmarg));
	 if (lefttrunk==1) ptrunc <- exp(-resi$cumhazleft); 
  } else if (inherits(margsurv,"coxph")) {  
    ### some problems here when data is different from data used in margsurv
       residuals <- residuals(margsurv)
       cum <- cumhaz <- status-residuals
       psurvmarg <- exp(-cumhaz);
       cumhazleft <- rep(0,antpers)
       ptrunc <- rep(1,length(psurvmarg));
       RR<- exp(margsurv$linear.predictors+sum(margsurv$means*coef(margsurv)))
       cum <- cum/RR
        if ((lefttrunk==1)) {
         stop("Use cox.aalen function for truncation case \n");
         baseout <- survival::basehaz(margsurv,centered=FALSE);
         cumt <- cbind(baseout$time,baseout$hazard)
	 cumt <- cpred(cumt,start.time)[,2]
	 ptrunc <- exp(-cumt * RR)
	}
  } else if (inherits(margsurv,"phreg")) {

	ppsurvmarg <- predict(margsurv,data,tminus=TRUE,times=time2,individual.time=TRUE,se=FALSE)
        psurvmarg <- ppsurvmarg$surv
        cum <- ppsurvmarg$cumhaz
	RR <- ppsurvmarg$RR
        ptrunc <- rep(1,length(psurvmarg));
        if ((lefttrunk==1)) {
	  ptrunc <- predict(margsurv,data,tminus=TRUE,times=start.time,individual.time=TRUE,se=FALSE)$surv
	}
  }
}

if (is.null(clusters) &  (!inherits(margsurv,"phreg"))) stop("must give clusters")

return(list(psurvmarg=psurvmarg,ptrunc=ptrunc,entry=start.time,exit=time2,
	    status=status,clusters=clusters,cum=cum,RR=RR))

} #

##' @title Twostage survival model fitted by pseudo MLE
##'
##' @description
##' Fits Clayton-Oakes clustered  survival data
##' using marginals that are on Cox form in the likelihood for the dependence parameter
##' as in Glidden (2000). The dependence can be modelled via  a
##' \enumerate{
##' \item  Regression design on dependence parameter.
##' }
##'
##' We allow a regression structure for the indenpendent gamma distributed
##' random effects  and their variances that may depend on cluster covariates. So
##' \deqn{
##'  \theta = h( z_j^T \alpha)
##' }
##' where \eqn{z} is specified by theta.des . The link function can be the exp when var.link=1
##' @references
##'
##' Measuring early or late dependence for bivariate twin data
##' Scheike, Holst, Hjelmborg (2015), LIDA
##'
##' Twostage modelling of additive gamma frailty models for survival data.
##' Scheike and Holst, working paper
##'
##' Shih and Louis (1995) Inference on the association parameter in copula models for bivariate
##' survival data, Biometrics, (1995).
##'
##' Glidden (2000), A Two-Stage estimator of the dependence
##' parameter for the Clayton Oakes model, LIDA, (2000).
##'
##' @examples
##' data(diabetes)
##' dd <- phreg(Surv(time,status==1)~treat+cluster(id),diabetes)
##' oo <- twostageMLE(dd,data=diabetes)
##' summary(oo)
##'
##' theta.des <- model.matrix(~-1+factor(adult),diabetes)
##'
##' oo <-twostageMLE(dd,data=diabetes,theta.des=theta.des)
##' summary(oo)
##' @keywords survival
##' @author Thomas Scheike
##' @param margsurv Marginal model from phreg
##' @param data data frame
##' @param theta Starting values for variance components
##' @param theta.des design for dependence parameters, when pairs are given this is could be a
##' (pairs) x (numer of parameters)  x (max number random effects) matrix
##' @param var.link Link function for variance  if 1 then uses exp link
##' @param method type of opitmizer, default is Newton-Raphson "NR"
##' @param no.opt to not optimize, for example to get score and iid for specific theta
##' @param weights cluster specific weights, but given with length equivalent to data-set, weights for score equations
##' @param se.cluster specifies how the influence functions are summed before squared when computing the variance. Note that the id from the marginal model is used to construct MLE, and then these scores can be summed with the se.cluster argument. 
##' @param ... arguments to be passed to  optimizer
##' @export
twostageMLE <-function(margsurv,data=parent.frame(),
theta=NULL,theta.des=NULL,var.link=0,method="NR",no.opt=FALSE,weights=NULL,se.cluster=NULL,...)
{
 if (!inherits(margsurv,"phreg"))  stop("Must use phreg for this \n");

 clusters <- margsurv$cox.prep$id
 n <- nrow(margsurv$cox.prep$X)


 if (is.null(theta.des)==TRUE) ptheta<-1;
 if (is.null(theta.des)==TRUE) theta.des<-matrix(1,n,ptheta) else theta.des<-as.matrix(theta.des);
  ptheta<-ncol(theta.des);
  if (nrow(theta.des)!=n) stop("Theta design does not have correct dim");

  if (is.null(theta)==TRUE) {
      if (var.link==1) theta<- rep(0.0,ptheta);
      if (var.link==0) theta<- rep(1.0,ptheta);
  }
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta);
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta);

  max.clust <- length(unique(clusters))
  theta.iid <- matrix(0,max.clust,ptheta)

  xx <- margsurv$cox.prep
  nn <- length(xx$strata)

  if (is.null(weights)) weights <-  rep(1,nn)
  if (length(weights)!=nn) stop("Weights do not have right length")

  statusxx <- rep(0,length(xx$strata))
  statusxx[xx$jumps+1] <- 1
  xx$status <- statusxx
  mid <- max(xx$id)+1

  ## Ni.(t-), Ni.(t)
  Nsum <- cumsumstratasum(statusxx,xx$id,mid,type="all")
  ## Ni.(tau)
  Ni.tau <- sumstrata(statusxx,xx$id,mid)
  ## cumahz(T_i-)
  S0i2 <- S0i <- rep(0, length(xx$strata))
  S0i[xx$jumps + 1] <- 1/margsurv$S0
  cumhazD <- cumsumstratasum(S0i, xx$strata, xx$nstrata)$lagsum
  if (!is.null(margsurv$coef)) RR <- exp(xx$X %*% margsurv$coef) else RR  <-  rep(1,nn)
  H <- c(cumhazD * RR)

  cc <- cluster.index(xx$id)
  firstid <- cc$firstclustid+1
  if (max(cc$cluster.size)==1) stop("No clusters !, maxclust size=1\n");

  ### order after time-sorted data
  theta.des <- theta.des[xx$ord+1,,drop=FALSE]
  weightsid <- weights <- weights[xx$ord+1]
  ### clusterspecific weights
  weights <- weights[firstid]

  thetaX<- as.matrix(theta.des[firstid,,drop=FALSE])

obj <- function(par,all=FALSE)
{

if (var.link==1) epar <- c(exp(c(par))) else epar <- c(par)

thetav <- c(as.matrix(theta.des) %*%  c(epar))
thetai <-thetav[firstid]

Hs <- sumstrata(H*exp(thetav*H),xx$id,mid)
R <- sumstrata((exp(thetav*H)-1),xx$id,mid) + 1
H2 <- sumstrata(H^2*exp(thetav*H),xx$id,mid)

l1 <- sumstrata(log(1+thetav*Nsum$lagsum)*statusxx,xx$id,mid)
l2 <- sumstrata(statusxx*H,xx$id,mid)
l3 <- -((thetai)^{-1}+Ni.tau) * log(R)
logliid <- (l1 + thetai* l2 + l3)*c(weights)
logl <- sum(logliid)
ploglik <- logl

## first derivative
l1s <- sumstrata(Nsum$lagsum/(1+thetav*Nsum$lagsum)*statusxx,xx$id,mid)
l2s <- (thetai^{-2}) * log(R)
l3s <- -(thetai^{-1}+Ni.tau) * Hs / R
Dltheta <- (l1s+l2s+l3s+l2)*c(weights)
scoreiid <- thetaX* c(Dltheta)

## second derivative
D2N <- -sumstrata(Nsum$lagsum^2/(1+thetav*Nsum$lagsum)^2*statusxx,xx$id,mid)
Dhes  <- c(D2N+(2/thetai^2)*Hs/R-(2/thetai^3)*log(R)-(1/thetai+Ni.tau)*(H2*R-Hs*Hs)/R^2)
Dhes  <-  Dhes * c(weights)

if (var.link==1) {
scoreiid  <-  scoreiid * c(thetai);  
Dhes<- Dhes* thetai^2 + thetai *  Dltheta
}

gradient <- apply(scoreiid,2,sum)
hessian <-  crossprod(thetaX,thetaX*c(Dhes))
hess2 <-    crossprod(scoreiid)

val <- list(id=xx$id,score.iid=scoreiid,logl.iid=logliid,ploglik=ploglik,
            gradient=-gradient,hessian=hessian,hess2=hess2)

if (all) return(val);

with(val, structure(-ploglik,gradient=-gradient,hessian=hessian))
}

  opt <- NULL
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
          opt <- lava::NR(theta,obj,...)
          opt$estimate <- opt$par
      } else {
          opt <- nlm(obj,theta,...)
          opt$method <- "nlm"
      }
      cc <- opt$estimate;
###   names(cc) <- colnames(theta.des)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
  } else val <- c(list(coef=theta),obj(theta,all=TRUE))
  val$score <- val$gradient


  theta <- matrix(c(val$coef),length(c(val$coef)),1)
  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:length(c(theta)),sep="")
  if (length(thetanames)==length(c(theta))) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }

  hessianI <- solve(val$hessian)
  val$theta.iid.naive  <-  val$score.iid %*% hessianI

  ### iid due to Marginal model
  biid <- 1
  if (biid==1) {

  if ((margsurv$no.opt) | is.null(margsurv$coef)) fixbeta<- 1 else fixbeta <- 0

   id  <-  xx$id+1
   if (var.link==1) epar <- c(exp(theta)) else epar <- c(theta)

   thetaX<- as.matrix(theta.des[firstid,,drop=FALSE])
   thetav <- c(as.matrix(theta.des) %*%  c(epar))

   Hs <- c(sumstrata(H*exp(thetav*H),xx$id,mid))
   R <- c(sumstrata((exp(thetav*H)-1),xx$id,mid)) + 1
   H2 <- c(sumstrata(H^2*exp(thetav*H),xx$id,mid) )

   Ft <- ((1/(thetav*R[id]))*exp(thetav*H)-(1/thetav+Ni.tau[id])*(1+thetav*H)*exp(thetav*H)/R[id]+statusxx+(1+thetav*Ni.tau[id])*exp(thetav*H)*Hs[id]/R[id]^2)

  if (var.link==1) {
	Ft  <-  Ft * c(exp(theta.des %*% val$coef));
  }
  Gt <- c(RR *Ft * xx$sign * weightsid)
  Ft  <- c( Ft * H * weightsid )
  Fbeta <-  t(theta.des) %*% ( xx$X * Ft )
  Gbase <- apply(Gt* theta.des,2,revcumsumstrata,xx$strata,xx$nstrata)

  ### beta part
  if (fixbeta==0) {
	  Z <- xx$X
	  U <- E <- matrix(0,nrow(xx$X),margsurv$p)
	  E[xx$jumps+1,] <- margsurv$E
	  U[xx$jumps+1,] <- margsurv$U
          invhess <- -solve(margsurv$hessian)
	  S0i <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/margsurv$S0
	  cumhaz <- c(cumsumstrata(S0i,xx$strata,xx$nstrata))
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset))
	  ### Martingale  as a function of time and for all subjects to handle strata
	  MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx$weights)
          UUbeta <- apply(MGt,2,sumstrata,id-1,max(id))
	  UUbeta <- UUbeta %*% invhess
	  GbaseEdLam0 <- t(Gbase) %*% (E*S0i)
	  Fbeta   <-  Fbeta -  GbaseEdLam0
          Fbetaiid <- UUbeta %*% t(Fbeta)
  }

  ###  \int_0^\tau (GBase(s) / S_0(s)) dN_i(s) - dLamba_i(s)
  xx <- margsurv$cox.prep
  S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <- 1/margsurv$S0
  S0i2[xx$jumps+1] <- 1/margsurv$S0^2

  GbasedN <- Gbase*S0i
  GbasedLam0 <- apply(Gbase*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)

  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset)) else rr <- c(xx$sign*exp( xx$offset))
  MGt <- (GbasedN[,drop=FALSE]-GbasedLam0*rr)*c(xx$weights)
  MGti <- apply(MGt,2,sumstrata,id-1,max(id))

###  if (fixbeta==0) print(cbind(val$score.iid,MGti,Fbetaiid)) else print(cbind(val$score.iid,MGti))

  if (fixbeta==0) MGti <- MGti+Fbetaiid
  theta.iid <-  val$score.iid + MGti
  theta.iid  <-  theta.iid %*% hessianI
  ## take names from phreg call
  if (!is.null(margsurv$id)) rownames(theta.iid) <- unique(margsurv$id)
  val$theta.iid <- theta.iid

  }

  ### if se.cluster is there we cluster iid-decomp before squaring 
  if (!is.null(se.cluster)) 
	  if (length(se.cluster)!=length(clusters)) stop("Length of seclusters and clusters must be same\n");

  if (!is.null(se.cluster)) {
      iids <-  unique(se.cluster);
      nseclust  <- length(iids);
      if (is.numeric(se.cluster)) se.cluster <-  fast.approx(iids,se.cluster)-1
       else se.cluster <- as.integer(factor(se.cluster, labels = seq(nseclust)))-1

      val$theta.iid <-       apply(val$theta.iid,se.cluster,nseclust)
      val$theta.iid.naive <- apply(val$theta.iid.naive,se.cluster,nseclust)
  }

  var <- robvar.theta  <-  var.theta  <-  crossprod(val$theta.iid)
  naive.var  <-   crossprod(val$theta.iid.naive)

  val  <- c(val,list(theta=theta,var.theta=var.theta,robvar.theta=robvar.theta,
		     var=var,thetanames=thetanames,model="clayton.oakes",
	             se=diag(robvar.theta)^.5),
	             var.naive=naive.var)
  class(val) <- "mets.twostage"
  attr(val,"clusters") <- clusters
  attr(val,"secluster") <- c(se.cluster)
  attr(val,"var.link")<-var.link;
  attr(val,"ptheta")<-ptheta
  attr(val,"n")<-n ;
  attr(val,"response")<- "survival"
  attr(val,"additive-gamma")<-0
  attr(val,"twostage") <- "two.stage"

  return(val)
}

##' @export
summary.mets.twostage <- function(object,digits = 3,silent=0,theta.des=NULL,...)
{ #
  if (!(inherits(object,"mets.twostage"))) stop("Must be a Two-Stage object")

  var.link<-attr(object,"var.link");
  var.par  <- attr(object,"var.par");
  model <- object$model
  if ((model=="plackett" || model=="or") )         model <- "or"
  if ((model=="clayton.oakes" || model=="gamma") ) model <- "gamma"
  if ((attr(object,"additive-gamma")==1) & (silent==0)) addgam <- TRUE else addgam <- FALSE

  if ((model=="or") && (silent==0)) cat("Dependence parameter for Odds-Ratio (Plackett) model \n");
  if (attr(object,"response")=="binomial") response <- "binomial" else response <- "survival"
  if ((model=="gamma") && (silent==0)) {
	  cat("Dependence parameter for Clayton-Oakes model\n")
	  if (var.par==1 || !addgam) cat("Variance of Gamma distributed random effects \n");
	  if (var.par==0 && addgam) cat("Inverse of variance of Gamma distributed random effects \n");
  }
  if (var.link==1 && silent==0) cat("With log-link \n")
  if ((sum(abs(object$score))>0.0001) & (silent==0))  {
	  cat(" Variance parameters did not converge, allow more iterations.\n");
	  cat(paste("Score:",object$score,"  \n"));
  }
  theta <- object$theta
  if (is.null(rownames(theta)))  rownames(theta) <- paste("dependence",1:nrow(theta),sep="")

  coefs <- coef.mets.twostage(object,response=response,...);

  if (attr(object,"additive-gamma")==1 & (!is.null(object$robvar.theta))  ) {
      var.link <- attr(object,"var.link");
      var.par  <- attr(object,"var.par");
      rv1      <- attr(object,"rv1");
      if (is.matrix(rv1)) rv1 <- rv1[1,]
      theta.des <- attr(object,"pardes");
      ags <- attr(object,"additive.gamma.sum");
      ptheta <- attr(object,"ptheta")
      npar <- nrow(object$theta)
      theta <- object$theta[seq(1,ptheta),1,drop=FALSE]
      robvar.theta <- object$robvar.theta[seq(1,ptheta),seq(1,ptheta)]
      if (var.link==1) par <- theta.des %*% exp(theta) else  par <- theta.des %*% theta

      if (model=="or" || model=="plackett") var.par<-1 ## same as var.par=0 for this model

      if (attr(object,"twostage")==0) {
      }

      if (var.par==0) { #
	      if (var.link==1) {
		  fp <- function(p,d,t){ res <- exp(p*t)/(sum(rv1* c(theta.des %*% matrix(exp(p),ncol=1))))^d;
					     if (t==0) res <- res[1]; return(res); }
                  pare <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) exp(p),labels=rownames(theta))
	      } else {
		      fp <- function(p,d,t) {  res <- (p^t)/(sum(rv1* c(theta.des %*% matrix(p,ncol=1))))^d;
					     if (t==0) res <- res[1]; return(res); }
	          pare <- NULL
	      }

             e      <- lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,1),labels=rownames(theta))
             vare <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,2,1),labels=rownames(theta))
             vartot <- lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,0))
      } #
      if (var.par==1) { #

	   if (var.link==1) { #
	     fp <- function(p,d,t){  res <- exp(p*t)/(sum(rv1* c(theta.des %*% matrix(exp(p),ncol=1))))^d;
                                     if (t==0) res <- res[1]; return(res); }
             e      <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,1),labels=rownames(theta))
             vare   <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) exp(p),labels=rownames(theta))
             vartot <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) sum(exp(p)))
      } else {
              fp <- function(p,d,t) {  res <- (p^t)/(sum(rv1* c(theta.des %*% matrix(p,ncol=1))))^d;
                                     if (t==0) res <- res[1]; return(res); }
              e     <-  lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,1),labels=rownames(theta))
              pare  <-  lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,2,1),labels=rownames(theta))
              vartot <- lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) sum(p))
	      vare <- NULL
	       } #
      } #
      res <- list(estimates=coefs, type=attr(object,"Type"),h=e,vare=vare,vartot=vartot)

  } else {
	 if (var.link==1) { #
	     if (model=="or") {
		or <-  lava::estimate(coef=object$theta,vcov=object$var.theta,f=function(p) exp(p),labels=rownames(theta))
	        res <- list(estimates=coefs,or=or,type=attr(object,"Type"))
	     } else {
	        vargam <-  lava::estimate(coef=object$theta,vcov=object$var.theta,f=function(p) exp(p),labels=rownames(theta))
                res <- list(estimates=coefs,vargam=vargam, type=attr(object,"Type"))
	     }
      } else {
	     if (model=="or") {
             lor <- lava::estimate(coef=object$theta,vcov=object$var.theta,f=function(p) log(p),labels=rownames(theta))
	     res <- list(estimates=coefs,log.or=lor,type=attr(object,"Type"))
	     } else  {
             res <- list(estimates=coefs,type=attr(object,"Type"))
	     }
      } #
  }

  class(res) <- "summary.mets.twostage"
  res
} #

##' @export
coef.mets.twostage <- function(object,var.link=NULL,response="survival",...)
{ #
  pt <- attr(object,"ptheta")
  theta <- object$theta[seq(pt),1]

  if (is.null(object$robvar.theta))
	  var.theta  <-  object$var.theta[seq(1,pt),seq(1,pt),drop=FALSE] else
	  var.theta  <-   object$robvar.theta[seq(1,pt),seq(1,pt),drop=FALSE]
  se <- diag(var.theta)^.5

  if (is.null(var.link))
     if (attr(object,"var.link")==1) vlink <- 1 else vlink <- 0
     else vlink <- var.link
  res <- cbind(theta, se )
  wald <- theta/se
  waldp <- (1 - pnorm(abs(wald))) * 2
  if (response=="survival") {
       if (object$model=="plackett") {
       spearman <- alpha2spear(theta,link=vlink)
       Dspear <- numDeriv::jacobian(alpha2spear,theta,link=vlink)
       var.spearman <- Dspear %*% var.theta %*%  Dspear
       se.spearman <- diag(var.spearman)^.5
       res <- as.matrix(cbind(res, wald, waldp,spearman,se.spearman))
       if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Spearman Corr.","SE")
	  else colnames(res) <- c("Coef.", "SE","z", "P-val","Spearman Corr.","SE")
	  if ((!is.null(object$thetanames)) & (nrow(res)==length(object$thetanames))) rownames(res)<-object$thetanames
       }
       if (object$model=="clayton.oakes") {
       kendall <- alpha2kendall(theta,link=vlink)
       Dken <- numDeriv::jacobian(alpha2kendall,theta,link=vlink)
       var.kendall<- Dken %*% var.theta %*%  Dken
       se.kendall <- diag(var.kendall)^.5
       res <- as.matrix(cbind(res, wald, waldp,kendall,se.kendall))
       if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Kendall tau","SE")
       else colnames(res) <- c("Coef.", "SE","z", "P-val","Kendall tau","SE")
       if ((!is.null(object$thetanames)) & (nrow(res)==length(object$thetanames))) rownames(res)<-object$thetanames
       }
  }
  return(res)
} #

##' @export
print.mets.twostage<-function(x,digits=3,...)
{ #
  cat("\n")
  print(summary(x,silent=0));
} #

##' @export
plot.mets.twostage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
		level=0.05,
		start.time=0,stop.time=0,add.to.plot=FALSE,mains=TRUE,
                xlab="Time",ylab ="Cumulative regression function",...)
{ #
  if (!(inherits(x, 'two.stage'))) stop("Must be a Two-Stage object")
  object <- x; rm(x);

  B<-object$cum; V<-object$var.cum; p<-dim(B)[[2]];
  if (robust>=1) V<-object$robvar.cum;

  if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1
  if (stop.time==0) stop.time<-max(B[,1]);

  med<-B[,1]<=stop.time & B[,1]>=start.time
  B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
  V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs);
  Vrob<-object$robvar.cum;
  Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs);

  c.alpha<- qnorm(1-level/2)
  for (v in comp) {
    c.alpha<- qnorm(1-level/2)
    est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
    if (add.to.plot==FALSE)
      {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab)
        if (mains==TRUE) title(main=colnames(B)[v]); }
    else lines(B[,1],est,type="s");
    if (pointwise.ci>=1) {
      lines(B[,1],ul,lty=pointwise.ci,type="s");
      lines(B[,1],nl,lty=pointwise.ci,type="s"); }
    if (robust>=1) {
      lines(B[,1],ul,lty=robust,type="s");
      lines(B[,1],nl,lty=robust,type="s"); }
    abline(h=0);
  }
}  #

##' @export
matplot.mets.twostage <- function(object,...)
{ #
B <- object$baseline
matplot(B[,1],B[,-1],type="s",...)
} #

##' @export
predict.mets.twostage <- function(object,X=NULL,Z=NULL,times=NULL,times2=NULL,theta.des=NULL,diag=TRUE,...)
{ #
time.coef <- data.frame(object$cum)
if (!is.null(times)) {
cum <- cpred(object$cum,times);
cum2 <- cpred(object$cum,times);
} else { cum <- object$cum; cum2 <- object$cum }
if (!is.null(times2)) cum2 <- cpred(object$cum,times2);

if (is.null(X)) X <- 1;
if (is.null(X) & (!is.null(Z))) { Z <- as.matrix(Z);  X <- matrix(1,nrow(Z),1)}
if (is.null(Z) & (!is.null(X)))  {X <- as.matrix(X);  Z <- matrix(0,nrow(X),1); gamma <- 0}

if (diag==FALSE) {
   time.part <-  X %*% t(cum[,-1])
   time.part2 <-  X %*% t(cum2[,-1])
   if (!is.null(object$gamma)) { RR <- exp( Z %*% gamma );
       cumhaz <- t( t(time.part) * RR ); cumhaz2 <- t( t(time.part2) * RR )}
	    else { cumhaz <- time.part;  cumhaz2 <- time.part2;   }
} else {
	time.part <-  apply(as.matrix(X*cum[,-1]),1,sum)
	time.part2 <-  apply(as.matrix(X*cum2[,-1]),1,sum)
}

if (!is.null(object$gamma)) {
	RR<- exp(Z%*%gamma);
	cumhaz <- t( t(time.part) * RR );
	cumhaz2 <- t( t(time.part2) * RR )} else {
		cumhaz <- time.part;  cumhaz2 <- time.part2;
}
S1 <- exp(-cumhaz); S2 <- exp(-cumhaz2)

if (attr(object,"var.link")==1) theta  <- exp(object$theta) else theta <- object$theta
if (!is.null(theta.des)) theta <- c(theta.des %*% object$theta)

if (diag==FALSE) St1t2<- (outer(c(S1)^{-(theta)},c(S2)^{-(theta)},FUN="+") - 1)^(-(1/theta)) else
St1t2<- ((S1^{-(theta)}+S2^{-(theta)})-1)^(-(1/theta))

out=list(St1t2=St1t2,S1=S1,S2=S2,times=times,times2=times2,theta=theta)
return(out)
} #

##' @export ascertained.pairs
ascertained.pairs <-function (pairs,data,cr.models,bin=FALSE)
{
      timestatus <- all.vars(cr.models)
      ### let first event by second column and only
      ### use pairs where first is event
      apairs <- pairs
      if (bin==TRUE) fj <- ifelse(data[pairs[,1],timestatus[1]] > data[pairs[,2],timestatus[1]],1,2)
      else fj <- ifelse(data[pairs[,1],timestatus[1]] < data[pairs[,2],timestatus[1]],1,2)
      ### change order when 1st comes first
      apairs[fj==1,1] <- pairs[fj==1,2]
      apairs[fj==1,2] <- pairs[fj==1,1]
      ### only take pairs where first is a jump
      if (bin==FALSE) {
	      jmpf <- (data[apairs[,2],timestatus[2]]==1)
	      apairs <- apairs[data[apairs[,2],timestatus[2]]==1,]
              attr(pairs,"jump-first") <- jmpf
      }
      pairs <- apairs
      return(pairs)
} 

##' @export
alpha2spear <- function(theta,link=1) { #
   if (link==1) theta <- exp(theta)
if (length(theta)>1) {
   out <- c()
   for (thet in theta) {
   if (thet!=1) out <- c(out,( (thet+1)/(thet-1) -2* thet* log(thet)/ (thet-1)^2))
   else out <- c(out,0)
   }
} else { if (theta!=1) out <- ( (theta+1)/(theta-1) -2* theta* log(theta)/ (theta-1)^2) }

return(out)
} #

##' @export
alpha2kendall <- function(theta,link=0) {  #
   if (link==1) theta <- exp(theta)
   return(1/(1+2/theta))
} #

##' @export piecewise.twostage
piecewise.twostage <- function(cut1,cut2,data=parent.frame(),timevar="time",status="status",id="id",covars=NULL,covars.pairs=NULL,num=NULL,
            method="optimize",Nit=100,detail=0,silent=1,weights=NULL,
            control=list(),theta=NULL,theta.des=NULL,var.link=1,
	    step=0.5,model="plackett",data.return=0)
{ #{{{
iid <- 1
ud <- list()
if (missing(cut2)) cut2 <- cut1;
nc1 <- length(cut1); nc2 <- length(cut2)
names1 <- names2 <- c()
theta.mat <- se.theta.mat <- cor.mat <- score.mat <- se.cor.mat <- matrix(0,nc1-1,nc2-1);
clusters <- data[,id]
cluster.call <- clusters
idi <- unique(data[,id]);

if (iid==1) { theta.iid <- matrix(0,length(idi),(nc1-1)*(nc2-1));
              rownames(theta.iid) <- idi
            } else theta.iid <- NULL

thetal <- c()
k <- 0;
for (i1 in 2:nc1)
for (i2 in 2:nc2)
{
k <-(i1-2)*(nc2-1)+(i2-1)
if (silent<=0) cat(paste("Data-set ",k,"out of ",(nc1-1)*(nc2-1)),"\n");
datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data,timevar=timevar,
status=status,id=id,covars=covars,covars.pairs=covars.pairs,num=num,silent=silent)
if (silent<=-1) print("back in piecewise.twostage");
if (silent<=-1) print(summary(datalr));
if (silent<=-1) print(head(datalr));
if (silent<=-1) print(summary(datalr[,id]));
 boxlr <- list(left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]))
datalr$tstime <- datalr[,timevar]
datalr$tsstatus <- datalr[,status]
datalr$tsid <- datalr[,id]
###

if (is.null(covars))
f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(",num,")")))
else f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(",num,"):",covars)))
marg1 <- timereg::aalen(f,data=datalr,n.sim=0,robust=0)

fitlr<-  survival.twostage(marg1,data=datalr,clusters=datalr$tsid,
,model=model, Nit=Nit,detail=detail,silent=silent,weights=weights,
theta=theta,theta.des=theta.des,var.link=var.link,step=step)
###fitlr<-  survival.twostageCC(marg1,data=datalr,clusters=datalr$tsid,model=model,method=method,
###Nit=Nit,detail=detail,silent=silent,weights=weights,
###baseline.iid=0,control=control,
###theta=theta,theta.des=theta.des,var.link=var.link,step=step)
####
coef <- coef(fitlr)
theta.mat[i1-1,i2-1] <- fitlr$theta
se.theta.mat[i1-1,i2-1] <- fitlr$var.theta^.5
cor.mat[i1-1,i2-1] <- coef[1,5]
se.cor.mat[i1-1,i2-1] <- coef[1,6]
score.mat[i1-1,i2-1] <- fitlr$score
if (data.return==0)
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr)
if (data.return==1)
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr,data=datalr)
if (i2==2) names1 <- c(names1, paste(cut1[i1-1],"-",cut1[i1]))
if (i1==2) names2 <- c(names2, paste(cut2[i2-1],"-",cut2[i2]))
thetal <- c(thetal,fitlr$theta)

if ((silent<=-1) & (iid==1)) print(head(fitlr$theta.iid));
if ((silent<=-1) & (iid==1)) {
print(idi) ; print(datalr$tsid)
print(dim(fitlr$theta.iid))
print(head(fitlr$theta.iid))
print(dim(theta.iid))
print(length( idi %in% unique(datalr$tsid)))
}
if (iid==1) theta.iid[idi %in% unique(datalr$tsid),k] <-c(fitlr$theta.iid)
###if (iid==1) theta.iid[rownames(fitlr$theta.iid),k] <-  fitlr$theta.iid
}

var.thetal <- NULL
if (iid==1)  var.thetal <- t(theta.iid) %*% theta.iid

colnames(score.mat) <- colnames(cor.mat) <-  colnames(se.cor.mat)  <- colnames(se.theta.mat) <- colnames(theta.mat) <- names1;
rownames(score.mat) <- rownames(cor.mat) <-  rownames(se.cor.mat) <-  rownames(se.theta.mat) <- rownames(theta.mat) <- names2;

ud <- list(model.fits=ud,theta=theta.mat,var.theta=se.theta.mat^2,
	   se.theta=se.theta.mat,thetal=thetal,thetal.iid=theta.iid,var.thetal=var.thetal,model=model,
	   cor=cor.mat,se.cor=se.cor.mat,score=score.mat);
class(ud)<-"pc.twostage"
attr(ud,"var.link")<-var.link;
attr(ud, "Type") <- model
return(ud);
} #}}}

##' @export piecewise.data
piecewise.data <- function(cut1,cut2,data=parent.frame(),timevar="time",status="status",id="id",covars=NULL,covars.pairs=NULL,num=NULL,silent=1)
{ #
ud <- list()
if (missing(cut2)) cut2 <- cut1;
nc1 <- length(cut1); nc2 <- length(cut2)
dataud <- c()

k <- 0;
for (i1 in 2:nc1)
for (i2 in 2:nc2)
{
k <-(i1-2)*(nc2-1)+(i2-1)
if (silent<=0) cat(paste("Data-set ",k,"out of ",(nc1-1)*(nc2-1)),"\n");
 datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data,timevar=timevar,
			status=status,id=id,covars=covars,covars.pairs=covars.pairs,num=num,silent=silent)
if (silent<=-1) print(summary(datalr));
if (silent<=-1) print(head(datalr));
datalr$tstime <- datalr[,timevar]
datalr$tsstatus <- datalr[,status]
datalr$tsid <- datalr[,id]
datalr$strata <- paste( c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),collapse=",",sep="-")
datalr$intstrata <-
c(paste(c(cut1[i1-1],cut1[i1]),collapse=",",sep="-"),paste( c(cut2[i2-1],cut2[i2]),collapse=",",sep="-"))

if (silent<=-1) print(head(datalr));
dataud <- rbind(dataud,datalr)
}

return(data.frame(dataud))
} #

##' @export
summary.pc.twostage <- function(object,var.link=NULL,...)
{ #
  if (!(inherits(object,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")

  res <- list(estimates=object$theta,se=object$se.theta,cor=object$cor,se.cor=object$se.cor,
	      model=object$model,score=object$score)
  class(res) <- "summary.pc.twostage"
  attr(res,"var.link")<-attr(object,"var.link");
  attr(res, "Type") <- object$model
  res
} #

##' @export
print.pc.twostage <- function(x,var.link=NULL,...)
{ #
   if (!(inherits(x,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
   print( summary(x,var.link=var.link,...))
} #

##' @export
print.summary.pc.twostage <- function(x,var.link=NULL, digits=3,...)
{ #

  if (is.null(var.link)) { if (attr(x,"var.link")==1) vlink <- 1 else vlink <- 0; } else vlink <- var.link
  print(vlink)

  if (x$model=="plackett") cat("Dependence parameter for Plackett model \n");
  if (x$model=="clayton.oakes") cat("Dependence parameter for Clayton-Oakes model \n");

  if (max(x$score)>0.001) { cat("Score of log-likelihood for parameter estimates (too large?)\n"); print(x$score);cat("\n\n");}

  if (vlink==1) cat("log-coefficient for dependence parameter (SE) \n")  else cat("Dependence parameter (SE) \n");
  print(coefmat(x$estimate,x$se,digits=digits,...))
  cat("\n")

  if (x$model=="plackett") {cat("Spearman Correlation (SE) \n");cor.type <- "Spearman Correlation"; }
  if (x$model=="clayton.oakes") {cat("Kendall's tau (SE) \n"); cor.type <- "Kendall's tau";}

  print(coefmat(x$cor,x$se.cor,digits,...))
  cat("\n")
} #

##' @export
coefmat <- function(est,stderr,digits=3,...) { #
  myest <- round(10^digits*(est))/10^digits;
  myest <- paste(ifelse(myest<0,""," "),myest,sep="")
  mysd <- round(10^digits*(stderr))/10^digits;
  res <- matrix(paste(format(myest)," (",format(mysd),")",sep=""),ncol=ncol(est))
  dimnames(res) <- dimnames(est)
  colnames(res) <- paste("",colnames(res))
  noquote(res)
} #


##' @export
simSurvFam <- function(n,beta=0.0,theta=1,lam0=0.5,lam1=1,lam2=1,ctime=10,...) { #
xm <- rbinom(n,1,0.5); xf <- rbinom(n,1,0.5);
xb1 <- rbinom(n,1,0.5); xb2 <- rbinom(n,1,0.5);
###
zf <- rgamma(n,shape=lam1); zb <- rgamma(n,shape=lam2);
tm <- rexp(n)/(zf*exp(xm*beta)*lam0)
tf <- rexp(n)/(zf*exp(xf*beta)*lam0)
tb1 <- rexp(n)/((zf+zb)*exp(xb1*beta)*2*lam0)
tb2 <- rexp(n)/((zf+zb)*exp(xb2*beta)*2*lam0)
cm <- ifelse(tm<ctime,1,0); cf <- ifelse(tf<ctime,1,0);
cb1 <- ifelse(tb1<ctime,1,0); cb2 <- ifelse(tb2<ctime,1,0);
tm <- ifelse(tm<ctime,tm,ctime); tf <- ifelse(tf<ctime,tf,ctime)
tb1 <- ifelse(tb1<ctime,tb1,ctime); tb2 <- ifelse(tb2<ctime,tb2,ctime)
#
data.frame(xm=xm,xf=xf,xb1=xb1,xb2=xb2,timem=tm,timef=tf,timeb1=tb1,timeb2=tb2,statusm=cm,statusf=cf,
	   statusb1=cb1,statusb2=cb2,id=1:n)
} #

##' @export
object.defined <- function(object)
{
   exists(as.character(substitute(object)))
}

##' @export
twin.polygen.design <-function (data,id="id",zyg="DZ",zygname="zyg",type="ace",tv=NULL,...) { #
  ### twin case
  id <- data[,id]
  tv <- diff(c(NA,id))
  tv[tv!=0 | is.na(tv)] <- 1
  tv[tv==0] <- 2

  zygbin <- (data[,zygname]==zyg)*1
  zygdes=model.matrix(~-1+factor(zygbin),data)
  n <- length(zygbin)

  if (type=="ace") { ### ace #
  ### random effects for each cluster
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+
		(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns,1)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2","env")
  pard <- rbind(c(1,0), c(0.5,0),c(0.5,0), c(0.5,0), c(0,1))
  } #

  if (type=="ae") { ### ae #
  ###AE model
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+
		(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2")
  pard <- rbind(c(1,0), c(0.5,0),c(0.5,0), c(0.5,0))[,1,drop=FALSE]
  } #

  if (type=="dce") { ### dce #
  ### DCE
  ### random effects for each cluster
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns,1)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2","env")
  pard <- rbind(c(1,0), c(0.25,0),c(0.75,0), c(0.75,0), c(0,1))
  } #

  if (type=="ade") { ### ade #
  #ADE
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns,DZns,1)
  pard <- rbind(c(1,0),c(0.25,0),c(0.75,0),c(0.75,0),c(0,1),c(0,0.5),c(0,0.5),c(0,0.5) )
  pardes <- matrix(pard,n,16,byrow=TRUE)
  des.rv <- NULL
  } #

  if (type=="adce") { ### adce #
  } #

  if (type=="de") { ### ae #
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+
		(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2")
  pard <- rbind(c(1,0), c(0.25,0),c(0.75,0), c(0.75,0))[,1,drop=FALSE]
  } #

  if (type=="un") { ### ae #
  des.rv <- cbind(zygdes)
  colnames(des.rv) <- c("MZ","DZ")
  pard <- rbind(c(1,0),c(0,1))
  } #

res <- list(pardes=pard,des.rv=des.rv)
return(res)
} #

##' @export
ace.family.design <-function (data,id="id",member="type",mother="mother",father="father",child="child",child1="child",type="ace",...) {
#
  ### standard family case
###  nid <- table(data[,id])
  id <- data[,id]

  if (type=="ace") { ### ace #
  ### random effects for each cluster
	  mom <- 1*(data[,member]==mother)
	  mo <- cbind(mom*1,1,1,1)*(mom)
	  fa <- (data[,member]==father)
	  fad <- cbind(fa,1,1,1)*fa
	  ch1 <- (data[,member]==child)*(data[,child1]==1)
	  cc1 <- cbind(ch1,1,0,0,1,1,0,0)*ch1
	  ch2 <- (data[,member]==child)*(data[,child1]==0)
	  cc2 <- cbind(ch2,0,1,0,1,0,1,0)*ch2
	  des.rv <- cbind(cbind(mo,fad)+cc1+cc2,1)

	  colnames(des.rv) <- c("m1","m2","m3","m4","f1","f2","f3","f4","env")
	  pard <- rbind( c(0.25,0), c(0.25,0),c(0.25,0), c(0.25,0),
			 c(0.25,0), c(0.25,0),c(0.25,0), c(0.25,0), c(0,1))
  } #

  if (type=="ae") { ### ae #
  ###AE model
          mom <- 1*(data[,member]==mother)
	  mo <- cbind(mom*1,1,1,1)*(mom)
	  fa <- (data[,member]==father)
	  fad <- cbind(fa,1,1,1)*fa
	  ch1 <- (data[,member]==child)*(data[,child1]==1)
	  cc1 <- cbind(ch1,1,0,0,1,1,0,0)*ch1
	  ch2 <- (data[,member]==child)*(data[,child1]==0)
	  cc2 <- cbind(ch2,0,1,0,1,0,1,0)*ch2
	  des.rv <- cbind(cbind(mo,fad)+cc1+cc2)
	  colnames(des.rv) <- c("m1","m2","m3","m4","f1","f2","f3","f4")
	  pard <- rbind( c(0.25), c(0.25),c(0.25), c(0.25),
			 c(0.25), c(0.25),c(0.25), c(0.25))
  } #

  if (type=="dce") { ### dce #
  ### DCE
  ### random effects for each cluster
	  stop("not done yet");
          mom <- 1*(data[,member]==mother)
	  mo <- cbind(mom*1,1,1,1)*(mom)
	  fa <- (data[,member]==father)
	  fad <- cbind(fa,1,1,1)*fa
	  ch1 <- (data[,member]==child)*(data[,child1]==1)
	  cc1 <- cbind(ch1,1,0,0,1,1,0,0)*ch1
	  ch2 <- (data[,member]==child)*(data[,child1]==0)
	  cc2 <- cbind(ch2,0,1,0,1,0,1,0)*ch2
	  des.rv <- cbind(cbind(mo,fad)+cc1+cc2)
	  colnames(des.rv) <- c("m1","m2","m3","m4","f1","f2","f3","f4")
	  pard <- rbind( c(0.25), c(0.25),c(0.25), c(0.25),
			 c(0.25), c(0.25),c(0.25), c(0.25))

  } #

  if (type=="ade") { ### ade #
  #ADE
    stop("not done yet");
###    pard <- rbind(c(1,0),c(0.25,0),c(0.75,0),c(0.75,0),c(0,1),c(0,0.5),c(0,0.5),c(0,0.5) )
###    pardes <- matrix(pard,n,16,byrow=TRUE)
  } #

  if (type=="adce") { ### adce #
	  stop("not done yet");
  } #

  if (type=="de") { ### ae #
    stop("not done yet");
    pard <- rbind(c(1,0), c(0.25,0),c(0.75,0), c(0.75,0))[,1,drop=FALSE]
  } #

  if (type=="un") { ### ae #
	 stop("not done yet");
         pard <- diag(4)
  } #

res <- list(pardes=pard,des.rv=des.rv)
return(res)
} #

##' @export
make.pairwise.design  <- function(pairs,kinship,type="ace")
{ #
### makes pairwise random effects design for shared and non-shared random effects
### kinship gives shared genes for each pair

uk <- sort(unique(kinship))
iuk <- fast.approx(uk,kinship)

if (type=="ace")  {
   theta.des <-  matrix(0,length(uk),8)
   i <- 0
   for (uuk in uk)  {
     i <- i+1
     theta.des[i,] <-c( rbind(c(uuk,0), c(1-uuk,0), c(1-uuk,0), c(0,1)))
   }
   random.des <- rbind(c(1,1,0,1),c(1,0,1,1))
   rvs <- rep(4,length(kinship))
   new.pairs <- cbind(pairs,1,2,iuk,4)
}

if (type=="ae") {
  i <- 0
  theta.des <-  matrix(0,length(uk),3)
  for (uuk in uk)  {
    i <- i+1
    theta.des[i,] <- rbind(c(uuk), c(1-uuk), c(1-uuk))
  }
  random.des <- rbind(c(1,1,0),c(1,0,1))
  rvs <- rep(3,length(kinship))
  new.pairs <- cbind(pairs,1,2,iuk,3)
}

if (type=="ade") {
   stop(" not specified yet \n")
   theta.des <-  matrix(0,length(uk),8)
   i <- 0
   for (uuk in uk)  {
      i <- i+1
      theta.des[i,] <- rbind(c((uuk==1)+(uuk!=1)*uuk*0.5,0), c(1-(uuk==1)+(uuk!=1)*uuk*0.5,0),
			c(1-(uuk==1)+(uuk!=1)*uuk*0.5,0),c(0,1))
   }
   random.des <- rbind(c(1,1,0,1),c(1,0,1,1))
   rvs <- rep(4,length(kinship))
   new.pairs <- cbind(pairs,1,2,iuk,4)
}

if (type=="ad") {
   stop(" not specified yet \n")
   i <- 0
   for (uuk in uk)  {
   i <- i+1
   theta.des[i,] <- rbind(c((uuk==1)+(uuk!=1)*uuk*0.5,0), c(1-(uuk==1)+(uuk!=1)*uuk*0.5,0),
			c(1-(uuk==1)+(uuk!=1)*uuk*0.5,0),c(0,1))
   random.des <- rbind(c(1,1,0),c(1,0,1))
   new.pairs <- cbind(pairs,1,2,iuk,2)
   }
}

return(list(new.pairs=new.pairs,theta.des=theta.des,random.design=random.des))
} #


##' Relative risk for additive gamma model
##'
##' Computes the relative risk for additive gamma model at time 0
##'
##' @references
##'
##' Eriksson and Scheike (2015), Additive Gamma frailty models for competing risks data, Biometrics (2015)
##'
##' @examples
##' lam0 <- c(0.5,0.3)
##' pars <- c(1,1,1,1,0,1)
##' ## genetic random effects, cause1, cause2 and overall
##' parg <- pars[c(1,3,5)]
##' ## environmental random effects, cause1, cause2 and overall
##' parc <- pars[c(2,4,6)]
##'
##' ## simulate competing risks with two causes with hazards 0.5 and 0.3
##' ## ace for each cause, and overall ace
##' out <- simCompete.twin.ace(10000,parg,parc,0,2,lam0=lam0,overall=1,all.sum=1)
##'
##' ## setting up design for running the model
##' mm <- familycluster.index(out$cluster)
##' head(mm$familypairindex,n=10)
##' pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
##' tail(pairs,n=12)
##' #
##' kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5
##'
##' # dout <- make.pairwise.design.competing(pairs,kinship,
##' #          type="ace",compete=length(lam0),overall=1)
##' # head(dout$ant.rvs)
##' ## MZ
##' # dim(dout$theta.des)
##' # dout$random.design[,,1]
##' ## DZ
##' # dout$theta.des[,,nrow(pairs)]
##' # dout$random.design[,,nrow(pairs)]
##' #
##' # thetades <- dout$theta.des[,,1]
##' # x <- dout$random.design[,,1]
##' # x
##' ##EVaddGam(rep(1,6),x[1,],x[3,],thetades,matrix(1,18,6))
##'
##' # thetades <- dout$theta.des[,,nrow(out)/2]
##' # x <- dout$random.design[,,nrow(out)/2]
##' ##EVaddGam(rep(1,6),x[1,],x[4,],thetades,matrix(1,18,6))
##' @author Thomas Scheike
##' @export
##' @param theta theta
##' @param x1 x1
##' @param x2 x2
##' @param thetades thetades
##' @param ags ags
EVaddGam <- function(theta,x1,x2,thetades,ags)
{ #
	pars <- thetades %*% theta
	lamtot <- ags %*% theta

	mvar <- pars/lamtot
	vvar <- pars/lamtot^2

	x1mvar <- sum(x1 * mvar)
	x2mvar <- sum(x2 * mvar)
	x1x2vvar <- sum(x1*x2*vvar)

	list(x1m=x1mvar,mx2=x2mvar,
	     dN=x1x2vvar/x2mvar)
} #


##' @export
twostage.aalen <- function(object,...) survival.twostage(object,...)

##' @export
twostage.cox.aalen <- function(object,...) survival.twostage(object,...)

##' @export
twostage.coxph <- function(object,...) survival.twostage(object,...)

##' @export
twostage.phreg <- function(object,...) survival.twostage(object,...)
