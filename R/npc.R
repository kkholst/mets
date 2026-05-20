##' Non-parametric Cumulative Incidence Functions
##'
##' Functions for computing and visualizing non-parametric cumulative incidence
##' estimates, as well as dependence measures (odds ratio, relative risk) for
##' bivariate competing risks data.
##'
##' \code{npc} computes bivariate non-parametric cumulative incidence using
##' inverse-probability-of-censoring weights.
##'
##' \code{nonparcuminc} computes univariate non-parametric cumulative incidence
##' for multiple causes.
##'
##' \code{plotcr} plots cumulative incidence curves for competing risks using
##' the prodlim package.
##'
##' \code{or_cif} fits an odds-ratio model for bivariate cumulative incidence.
##'
##' \code{rr_cif} fits a relative-risk model for bivariate cumulative incidence.
##'
##' \code{random.cif} and \code{Grandom.cif} are aliases for \code{random_cif}
##' and \code{Grandom_cif} (random effects CIF models).
##'
##' \code{predictPairPlack} predicts pairwise joint probabilities under a
##' Plackett (odds-ratio) dependence model.
##'
##' \code{matplot.mets.twostage} produces matrix-plots of concordance over time
##' from a twostage object.
##'
##' @name cif-nonpar
##' @param T matrix with columns: time1, time2, status1, status2 (for \code{npc}).
##' @param cause vector of length 2 specifying causes of interest (for \code{npc}).
##' @param same.cens logical; if TRUE, uses joint censoring weights.
##' @param sep logical; if TRUE, uses separate censoring models for each subject.
##' @param cif a cumulative incidence model object (from timereg).
##' @param data a data.frame with the variables.
##' @param cif2 optional second CIF model if different from first.
##' @param times time points for evaluation.
##' @param cause1 cause for first coordinate.
##' @param cause2 cause for second coordinate.
##' @param cens.code censoring code value.
##' @param cens.model censoring model type (default \code{"KM"}).
##' @param Nit maximum number of iterations.
##' @param detail level of output detail.
##' @param clusters cluster variable name or vector.
##' @param theta dependence parameter(s).
##' @param theta.des design matrix for theta.
##' @param step step size for optimization.
##' @param sym if 1, symmetric dependence structure.
##' @param weights optional weights.
##' @param censoring.weights optional pre-computed censoring weights.
##' @param silent verbosity level.
##' @param par.func optional parameter function.
##' @param dpar.func optional derivative of parameter function.
##' @param dimpar dimension of parameter vector.
##' @param score.method optimization method (default \code{"nlminb"}).
##' @param entry optional entry time variable.
##' @param estimator estimator type.
##' @param trunkp truncation probability.
##' @param admin.cens administrative censoring time.
##' @param cif1 CIF values for subject 1 (for \code{predictPairPlack}).
##' @param status1 status for subject 1.
##' @param status2 status for subject 2.
##' @param x data matrix or competing risks object.
##' @param ... additional arguments.
##' @return For \code{npc}: matrix with columns (time, cumulative incidence).
##'   For \code{nonparcuminc}: matrix with time and cause-specific cumulative incidences.
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases npc nonparcuminc plotcr or_cif rr_cif random.cif Grandom.cif
##' @aliases predictPairPlack
NULL

##' @rdname cif-nonpar
##' @export
npc <- function(T,cause,same.cens=TRUE,sep=FALSE) {
  mtime <- apply(T[,1:2],1,max)
  ot <- order(mtime)
  mtime <- mtime[ot]
  T <- T[ot,]
  if (!sep) {
    time1 <- as.vector(T[,1:2]); status1 <- as.vector(T[,3:4])
    ud.cens1<-survival::survfit(Surv(time1,status1==0)~+1);
    Gfit1<-cbind(ud.cens1$time,ud.cens1$surv)
    Gfit2 <- Gfit1<-rbind(c(0,1),Gfit1);
  } else {
    time1 <- as.vector(T[,1]); status1 <- as.vector(T[,3])
    ud.cens1<-survival::survfit(Surv(time1,status1==0)~+1);
    time2 <- as.vector(T[,2]); status2 <- as.vector(T[,4])
    ud.cens2<-survival::survfit(Surv(time2,status2==0)~+1);
    Gfit1<-cbind(ud.cens1$time,ud.cens1$surv)
    Gfit1<-rbind(c(0,1),Gfit1);
    Gfit2<-cbind(ud.cens2$time,ud.cens2$surv)
    Gfit2<-rbind(c(0,1),Gfit2);
  }
  i1 <- fast.approx(Gfit1[,1],T[,1])
  cweights1<-fast.approx(,Gfit1[,2])[[1]]
  cweights2<-fast.approx(Gfit2[,1],T[,2],Gfit2[,2])[[1]];
  weight11 <- apply(cbind(cweights1,cweights2),1,min)

  if (same.cens) {
    conc <- (T[,3]==cause[1])*(T[,4]==cause[2])/weight11
  } else {
    conc <-(T[,3]==cause[1])*(T[,4]==cause[2])/(cweights1*cweights2);
  }
  mtime <- mtime[!is.na(conc)]
  conc <- conc[!is.na(conc)]
  cbind(mtime,cumsum(conc)/length(conc))
}

##' @rdname cif-nonpar
##' @param t vector of event/censoring times (for \code{nonparcuminc}).
##' @param status vector of status codes (for \code{nonparcuminc}).
##' @param cens censoring code (default 0).
##' @export
nonparcuminc <- function(t,status,cens=0) {
  ord <- order(t); t <- t[ord]; status <- status[ord]
  ud.cens<-survival::survfit(Surv(t,status==cens)~1)
  Gfit<-cbind(ud.cens$time,ud.cens$surv)
  Gfit<-rbind(c(0,1),Gfit);
  causes <- setdiff(unique(status),cens)
  cweight<-fast.approx(Gfit[,1],t,Gfit[,2])[[1]];
  cc <- t
  for (i in 1:length(causes)) {
    c1 <- status==causes[i]
    cc <- cbind(cc,cumsum(c1/cweight)/length(c1))
  }
  return(cc)
}
