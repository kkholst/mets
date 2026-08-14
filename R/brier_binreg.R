
##' Stratified Binomial Regression with cross-validation option 
##'
##' Fits a stratified binomial regression model for the cumulative incidence
##' function (CIF), restricted mean survival time (RMST), or restricted mean
##' time lost (RMTL), using inverse probability of censoring weighting (IPCW).
##' Each stratum may either be fit on its own observations or, when
##' \code{cv.fold = TRUE}, on the complement of its observations (i.e. a
##' leave-one-stratum-out / cross-validation style fit). Standard errors are
##' computed via an iid decomposition that accounts for the estimation of the
##' censoring distribution.
##'
##' @param formula A formula with a survival outcome, \code{\link[survival]{Surv}} or
##'   \code{Event} object, on the left-hand side optional
##'   \code{cluster()}, \code{strata()}, \code{offset()}, and \code{weights()}
##'   terms on the right-hand side. Left-hand side response can also be a numeric or a factor
##'  that is then  used as \code{as.numeric(factor)-1} for example to do logistic regression. 
##' @param data A \code{data.frame} containing the variables in \code{formula}.
##' @param cause Value(s) of the status variable that identify the event of
##'   interest (cause 1 by default).
##' @param time The fixed time-point at which the CIF, RMST, or RMTL is
##'   evaluated. Required.
##' @param beta Optional starting values for the regression coefficients.
##'   Either a vector (recycled across strata) or a matrix with one row per
##'   stratum.
##' @param type Character; either \code{"II"} (default) or \code{"I"}. When
##'   \code{"II"}, an additional augmentation term and an extra component of
##'   the censoring-related iid decomposition are computed.
##' @param offset Optional offset term; can also be specified via the formula.
##' @param weights Optional case weights; can also be specified via the
##'   formula.
##' @param cens.weights Optional pre-computed inverse probability of censoring
##'   weights. If supplied, the internal censoring model fit is skipped and
##'   \code{se} is forced to \code{FALSE}.
##' @param cens.model A formula for the censoring model (right-hand side only,
##'   fitted via \code{\link{phreg}}). Defaults to the Kaplan-Meier
##'   (\code{~+1}).
##' @param se Logical. If \code{TRUE} (default), standard errors are computed,
##'   including the contribution from estimating the censoring distribution.
##' @param kaplan.meier Logical. If \code{TRUE} (default), uses the
##'   Kaplan-Meier estimator for censoring weights when the censoring model has
##'   no covariates; otherwise a Cox-type estimator is used internally when
##'   covariates are present.
##' @param cens.code The code in the status variable denoting censoring
##'   (default 0).
##' @param no.opt Logical. If \code{TRUE}, no optimization is performed and the
##'   object is evaluated at the supplied \code{beta}.
##' @param method Optimization method. Currently only \code{"nr"} (Newton-Raphson
##'   via \code{lava::NR}) is supported.
##' @param augmentation Optional augmentation term(s) added to the score
##'   equations, either a vector (recycled across strata) or a matrix with one
##'   row per stratum.
##' @param outcome Character; the outcome scale to model: \code{"cif"}
##'   (default), \code{"rmst"}, or \code{"rmtl"}.
##' @param model Character; the link function: \code{"default"} (chooses
##'   \code{"logit"} for CIF and \code{"exp"} for RMST/RMTL), \code{"logit"},
##'   \code{"exp"}, or \code{"lin"}.
##' @param Ydirect Optional, directly supplied response vector, overriding the
##'   response otherwise constructed from \code{outcome}.
##' @param strata Optional numeric vector defining strata, used when strata are
##'   not specified via \code{strata()} in the formula.
##' @param cv.fold Logical. If \code{TRUE}, for each stratum \code{s} the
##'   regression coefficients are estimated using all observations \emph{not}
##'   in stratum \code{s} (a cross-validation-style fit), rather than the
##'   observations within stratum \code{s} (the default, \code{FALSE}).
##' @param low.memory Logical. If \code{TRUE} only saves key quanties, coef, var, iid and some design structures.
##' @param ... Additional arguments, currently passed as optimizer
##'   \code{control} settings (e.g. \code{tol}, \code{stepsize}) to the
##'   Newton-Raphson routine.
##'
##' @return An object of class \code{c("binregStrata","binreg")}, a list with
##'   components including:
##'   \item{coef}{Estimated regression coefficients (stacked across strata).}
##'   \item{iid}{Subject-level iid decomposition of the coefficients, adjusted
##'     for censoring estimation when \code{se = TRUE}.}
##'   \item{var, robvar}{Robust (sandwich) variance-covariance matrix.}
##'   \item{se.coef, se.robust}{Standard errors derived from \code{var}.}
##'   \item{hessian, ihessian}{Per-stratum Hessians and their (pseudo-)inverses.}
##'   \item{coef_list}{Per-stratum coefficient vectors.}
##'   \item{time, cause, cens.code, model, outcome}{Echoed arguments.}
##'   \item{Y, Yipcw, IPCW}{The constructed response, IPCW-weighted response,
##'     and IPCW weights.}
##'   \item{cens.weights, formC, cens.strata, cens.nstrata}{Censoring model
##'     weights and related information (when internally estimated).}
##'   \item{nstrata, strata_orig, strata_call, cv.fold}{Stratification
##'     information and whether the cross-validation fold scheme was used.}
##'   \item{n, nevent, ncluster, nid, id, name.id}{Sample size and
##'     identifier information.}
##'   \item{call, design}{The matched call and the design object.}
##'
##' @details
##' When \code{cv.fold = FALSE} (the standard case), the coefficients for
##' stratum \code{s} are estimated using only observations with
##' \code{strata == s}. When \code{cv.fold = TRUE}, the coefficients for
##' stratum \code{s} are instead estimated using observations with
##' \code{strata != s}, i.e. the complement; this also affects the
##' construction of the censoring-related iid terms (\code{MGCiid}) and the
##' per-stratum subject-level influence functions used for the robust
##' variance.
##'
##' Censoring weights are computed via a Cox model (\code{\link{phreg}}) on
##' \code{cens.model}, evaluated at \code{pmin(exit, time)}, unless
##' \code{cens.weights} is supplied directly. When \code{se = TRUE}, the
##' additional variability coming from the estimation of the censoring
##' distribution is incorporated into the standard errors via a martingale-type
##' iid decomposition (\code{MGCiid}).
##'
##' @seealso \code{\link{phreg}}, \code{\link{binreg}}
##'
##' @examples
##' \dontrun{
##' library(mets)
##' data(bmt)
##' out <- binregStrata(Event(time, cause) ~ tcell + strata(platelet),
##'                      data = bmt, cause = 1, time = 50)
##' summary(out)
##' }
##'
#' @export
binregStrata <- function(formula, data, cause=1, time=NULL, beta=NULL,
                          type=c("II","I"),
                          offset=NULL, weights=NULL, cens.weights=NULL,
                          cens.model=~+1, se=TRUE,
                          kaplan.meier=TRUE, cens.code=0, no.opt=FALSE,
                          method="nr", augmentation=NULL,
                          outcome=c("cif","rmst","rmtl"),
                          model=c("default","logit","exp","lin"),
                          Ydirect=NULL, strata=NULL,cv.fold=FALSE,low.memory=TRUE,...)
{ # {{{
  monotone <- TRUE
  cl       <- match.call()

  ## ------------------------------------------------------------------
  ## Design
  ## ------------------------------------------------------------------
  des <- proc_design(
    formula, data = data,
    specials = c("offset","weights","cluster","strata"),
    intercept = TRUE
  )

  Y <- des$y
  X <- des$x

  des.weights <- des$weights
  des.offset  <- des$offset
  id      <- des$cluster
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  ### possible handling of id to code from 0:(antid-1)
  call.id <- id 
  conid <- construct_id(id,nrow(X),namesX=rownames(X))
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
  ## id before time-sorting later 
  orig.id <- id

  ## take offset and weight first from formula, but then from arguments
  if (is.null(des.offset)) {
	  if (is.null(offset)) offset <- rep(0,nrow(X)) 
  } else offset <- des.offset
  if (is.null(des.weights)) {
	  if (is.null(weights)) weights <- rep(1,nrow(X)) 
  } else weights <- des.weights

  ## ------------------------------------------------------------------
  ## Strata  (explicit arg > formula-embedded special > default)
  ## ------------------------------------------------------------------
  des.strata <- des$strata
  if (!is.null(strata)) {
    if (!is.numeric(strata)) stop("strata must be numeric\n")
    ustrata <- sort(unique(strata)); nstrata <- length(ustrata)
    strata  <- fast.approx(ustrata, strata) - 1
  } else if (!is.null(des.strata)) {
    strata  <- des.strata
    ustrata <- sort(unique(strata)); nstrata <- length(ustrata)
    strata  <- fast.approx(ustrata, strata) - 1
  } else {
    strata  <- rep(0, nrow(X)); nstrata <- 1
  }


  ## ------------------------------------------------------------------
  ## Response
  ## ------------------------------------------------------------------
    if (inherits(Y, c("Event", "Surv"))) { ## {{{  Event time outcome 
	    if (ncol(Y) == 2) {
		exit <- Y[, 1]
		entry <- rep(0, nrow(Y))
		status <- Y[, 2]
	    } else {
		entry <- Y[, 1]
		exit <- Y[, 2]
		status <- Y[, 3]
	    }

	  if (is.null(time)) stop("Must give time for logistic modelling \n"); 
	  statusC <- 1*(status %in% cens.code) 
	  sumC  <- sum(statusC)
	  statusE <- (status %in% cause) & (exit<= time) 
	  if (sum(statusE)==0) warning("No events of type 1 before time \n"); 
	  kmt <- kaplan.meier

	  ucauses  <-  sort(unique(status))
	  ccc <- which(ucauses %in% cens.code)
	  if (length(ccc)==0) Causes <- ucauses else Causes <- ucauses[-ccc]
	  competing <- any(!(Causes %in% cause))
	  data$id__ <- id
	  data$exit <- exit
	  data$statusC <- statusC 
	  cens.strata <- cens.nstrata <- NULL 

	 nevent <- sum((status %in% cause)*(exit<=time))
	 ## if event before time or alive, then uncensored, equality for both censored and events  
	 obs <- (exit<=time & (statusC==0)) | (exit>=time)

	  if (is.null(cens.weights) & sumC > 0 )  { ## {{{
	      formC <- update.formula(cens.model,Surv(exit,statusC)~ . +cluster(id__))
	      resC <- phreg(formC,data)
	      if (resC$p>0) kmt <- FALSE
	      exittime <- pmin(exit,time)
	      cens.weights <- suppressWarnings(predict(resC,data,times=exittime,individual.time=TRUE,se=FALSE,km=kmt,tminus=TRUE)$surv)
	      ## strata from original data 
	      cens.strata <- resC$strata[order(resC$ord)]
	      cens.nstrata <- resC$nstrata
	  } else { se <- FALSE; resC <- formC <- NULL}  ## }}} 
	 if (is.null(cens.weights)) cens.weights <- rep(1,nrow(X))

	 ### setting up survival/competing risks outcomes 
	 if (!is.null(Ydirect)) Y <-  Ydirect else {
	     if (outcome[1]=="cif") Y <- c((status %in% cause)*(exit<=time))
	     else { if (!competing) {
		     if (outcome[1]=="rmst")
		     Y <-  c(pmin(exit,time))
		     else Y <-  c((time-pmin(exit,time)))
		    } else Y <- c((status %in% cause)*(time-pmin(exit,time)))
	     }
	  }

	 ## default links
	 if (model[1]=="default") {
		 if (outcome[1]=="cif") model <- "logit"
		 if (outcome[1]=="rmst") model <- "exp"
		 if (outcome[1]=="rmtl") model <- "exp"
		 if (outcome[1]=="years-lost") model <- "exp"
	 }

    }  else { ### numeric outcome or factor coded as numeric 
	    if (!(is.numeric(Y) | is.factor(Y))) stop("must be Event object, numeric, or a factor\n"); 
	    if (is.factor(Y)) Y <- as.numeric(Y)-1
	    cens.weights <- rep(1,nrow(X))
	    obs  <- 1
	    se  <- FALSE
	    cens.code <- Causes <- cause <- time  <- 0
	    nevent <- nrow(X)
	    formC  <- NULL
	    cens.strata <- cens.nstrata <- NULL

	    ## Ydirect rules over Y 
	    if (!is.null(Ydirect)) Y <- Ydirect

	 ## default links
	 if (model[1]=="default") {
		 model <- "lin"
	 }
    } ## }}} 


  expit <- function(z) 1 / (1 + exp(-z))
  p     <- ncol(X)

  if (is.null(beta))         beta         <- rep(0, p)
  if (is.null(augmentation)) augmentation <- matrix(0, nstrata, p)
  if (!is.matrix(augmentation))
    augmentation <- matrix(augmentation, nstrata, p, byrow = TRUE)

  X  <- as.matrix(X)
  X2 <- .Call("vecCPMat", X)$XX

  Yresp <- c(Y)
  Y     <- Yipcw <- c(Y * obs / cens.weights)
  IPCW  <- c(obs / cens.weights)

  ## ------------------------------------------------------------------
  ## MGCiid: censoring IID influence
  ## cv.fold: block s populated with complement rows (strata != s)
  ## standard: block s populated with stratum-s rows
  ## ------------------------------------------------------------------
  MGCiid <- 0
  if (se) { # {{{
    ord     <- resC$ord
    X       <- X[ord,  , drop = FALSE]
    X2      <- X2[ord, , drop = FALSE]
    status  <- status[ord]
    exit    <- exit[ord]
    weights <- weights[ord]
    offset  <- offset[ord]
    Y       <- Y[ord]
    strata  <- strata[ord]

    xx   <- resC$cox.prep
    S0i2 <- S0i <- rep(0, length(xx$strata))
    S0i[xx$jumps  + 1] <- 1 / resC$S0
    S0i2[xx$jumps + 1] <- 1 / resC$S0^2
    btime <- 1*(exit < time)
    mid   <- max(id)

    ## Block-expanded design: n x (nstrata*p)
    ## cv.fold: block s gets complement rows; standard: block s gets stratum-s rows
    Xs <- matrix(0, nrow(X), nstrata * p)
    for (s in 0:(nstrata - 1)) {
      idx  <- if (cv.fold) which(strata != s) else which(strata == s)
      cols <- s * p + seq_len(p)
      Xs[idx, cols] <- X[idx, ]
    }

    h           <- Mrevcumsumstrata(Xs * Y, xx$strata, xx$nstrata)
    IhdLam0     <- Mcumsumstrata(h * S0i2 * btime, xx$strata, xx$nstrata)
    U           <- matrix(0, nrow(xx$X), ncol(Xs))
    U[xx$jumps + 1,] <- (resC$jumptimes < time) *
                         h[xx$jumps + 1,] / c(resC$S0)
    MGt    <- (U - IhdLam0) * c(xx$weights)
    MGCiid <- Msumstrata(MGt, xx$id, mid + 1)

    if (type[1] == "II") { # {{{
      hYt      <- revcumsumstrata(Y,   xx$strata, xx$nstrata)
      IhdLam0  <- cumsumstrata(hYt * S0i2 * btime, xx$strata, xx$nstrata)
      U2       <- rep(0, length(xx$strata))
      U2[xx$jumps + 1] <- (resC$jumptimes < time) *
                           hYt[xx$jumps + 1] / c(resC$S0)
      MGt2    <- Xs * c(U2 - IhdLam0) * c(xx$weights)
      MGtiid  <- Msumstrata(MGt2, xx$id, mid + 1)
      augmentation <- colSums(MGtiid) + augmentation
      augmentation <- matrix(augmentation, nstrata, p, byrow = TRUE)

      EXt          <- Mrevcumsumstrata(Xs, xx$strata, xx$nstrata)
      IEXhYtdLam0  <- Mcumsumstrata(EXt * c(hYt) * S0i * S0i2 * btime,
                                     xx$strata, xx$nstrata)
      U3 <- matrix(0, nrow(xx$X), ncol(Xs))
      U3[xx$jumps + 1,] <- (resC$jumptimes < time) *
                            hYt[xx$jumps + 1] *
                            EXt[xx$jumps + 1,] / c(resC$S0)^2
      MGt3    <- (U3 - IEXhYtdLam0) * c(xx$weights)
      MGCiid2 <- Msumstrata(MGt3, xx$id, mid + 1)
      MGCiid  <- MGCiid + (MGtiid - MGCiid2)
    } # }}}

    id <- xx$id
  } else {
    MGCiid <- 0
  } # }}}

  ## ------------------------------------------------------------------
  ## Model default
  ## ------------------------------------------------------------------
  if (model[1] == "default") {
    if (outcome[1] == "cif")        model <- "logit"
    if (outcome[1] == "rmst")       model <- "exp"
    if (outcome[1] == "rmtl")       model <- "exp"
    if (outcome[1] == "years-lost") model <- "exp"
  }

  ## ------------------------------------------------------------------
  ## obj(): full objective (used for final all=TRUE call and NR)
  ## ------------------------------------------------------------------
  obj <- function(pp, all = FALSE) { # {{{
    np <- length(pp[[1]])

    ## ----------------------------------------------------------------
    ## cv.fold path: for stratum s, beta_s estimated on complement rows
    ## ----------------------------------------------------------------
    if (cv.fold) {
      ploglik_list  <- vector("list", nstrata)
      gradient_list <- vector("list", nstrata)
      D2log_list    <- vector("list", nstrata)

      for (s in 0:(nstrata - 1)) {
        idx_c <- which(strata != s)          # complement indices
        X_c   <- X[idx_c,  , drop=FALSE]
        X2_c  <- X2[idx_c, , drop=FALSE]
        Y_c   <- Y[idx_c]
        w_c   <- weights[idx_c]
        off_c <- offset[idx_c]
        aug_s <- augmentation[s + 1, ]
        lp_c  <- as.vector(X_c %*% pp[[s + 1]]) + off_c

        if (model[1] == "exp") {
          pr_c  <- exp(lp_c)
          if (monotone) {
            Dl_c  <- w_c * X_c * c(Y_c - pr_c)
            D2l_c <- c(w_c * pr_c) * X2_c
          } else {
            Dl_c  <- w_c * X_c * c((Y_c - pr_c) * pr_c)
            D2l_c <- c(w_c * pr_c^2) * X2_c
          }
          pl_c <- 0
        } else if (model[1] == "lin") {
          pr_c  <- lp_c
          Dl_c  <- w_c * X_c * c(Y_c - pr_c)
          D2l_c <- c(w_c) * X2_c
          pl_c  <- sum(w_c * (Y_c - pr_c)^2)
        } else if (model[1] == "logit") {
          pr_c <- expit(lp_c)
          if (monotone) {
            Dl_c  <- w_c * X_c * c(Y_c - pr_c)
            D2l_c <- c(w_c * pr_c / (1 + exp(lp_c))) * X2_c
          } else {
            varp  <- pr_c / (1 + exp(lp_c))
            Dl_c  <- w_c * X_c * c((Y_c - pr_c) * varp)
            D2l_c <- c(w_c * varp^2) * X2_c
          }
          pl_c <- 0
        } else stop("link must be logit, exp, or lin\n")

        ploglik_list[[s + 1]]  <- pl_c
        gradient_list[[s + 1]] <- -(colSums(Dl_c) + aug_s)
        D2log_list[[s + 1]]    <- colSums(D2l_c)
      }

      sym <- function(x)
        matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
                     PACKAGE="mets")$XXf, np, np)
      hessian <- lapply(D2log_list, sym)

      if (all) {
        ihess <- lapply(hessian, pinv)

        ## Per-stratum IID: recompute Dl_c at optimal pp[[s]] separately
        ## for each stratum so blocks are independent and correctly sized.
        ## Obs i in stratum r contributes to beta_s for all s != r;
        ## Msumstrata fills zeros for stratum-s obs (they don't contribute).
        beta.iid <- do.call(cbind, lapply(seq_len(nstrata), function(s) {
          idx_c <- which(strata != (s - 1))
          X_c   <- X[idx_c,  , drop=FALSE]
          Y_c   <- Y[idx_c]
          w_c   <- weights[idx_c]
          off_c <- offset[idx_c]
          id_c  <- id[idx_c]
          lp_c  <- as.vector(X_c %*% pp[[s]]) + off_c

          Dl_c <- switch(model[1],
            exp = {
              pr_c <- exp(lp_c)
              if (monotone) w_c * X_c * c(Y_c - pr_c)
              else          w_c * X_c * c((Y_c - pr_c) * pr_c)
            },
            lin = w_c * X_c * c(Y_c - lp_c),
            logit = {
              pr_c <- expit(lp_c)
              if (monotone) w_c * X_c * c(Y_c - pr_c)
              else {
                varp <- pr_c / (1 + exp(lp_c))
                w_c * X_c * c((Y_c - pr_c) * varp)
              }
            },
            stop("unknown model\n")
          )

          raw <- Dl_c %*% ihess[[s]]
          Msumstrata(raw, id_c, max(id) + 1)
        }))

        robvar <- crossprod(beta.iid)
        return(list(
          par      = pp,             coef      = unlist(pp),
          ploglik  = ploglik_list,   gradient  = gradient_list,
          hessian  = hessian,        ihessian  = ihess,
          id       = id,             Dlogl     = NULL,
          iid      = beta.iid,       robvar    = robvar,
          var      = robvar,         se.robust = diag(robvar)^.5
        ))
      }

      return(structure(ploglik_list,
                       gradient = gradient_list,
                       hessian  = hessian))
    }

# ----------------------------------------------------------------
# Standard non-CV
# ----------------------------------------------------------------
    loffset <- eta <- numeric(length(Y))
    for (s in 0:(nstrata - 1)) {
      idx          <- strata == s
      eta[idx]     <- as.vector(X[idx, , drop=FALSE] %*% pp[[s + 1]])
      loffset[idx] <- offset[idx]
    }
    lp <- eta + loffset

    if (model[1] == "exp") {
      pr     <- exp(lp)
      if (monotone) {
        Dlogl  <- weights * X * c(Y - pr)
        D2logl <- c(weights * pr) * X2
      } else {
        Dlogl  <- weights * X * c((Y - pr) * pr)
        D2logl <- c(weights * pr^2) * X2
      }
    } else if (model[1] == "lin") {
      pr     <- lp
      Dlogl  <- weights * X * c(Y - pr)
      D2logl <- c(weights) * X2
    } else if (model[1] == "logit") {
      pr <- expit(lp)
      if (monotone) {
        Dlogl  <- weights * X * c(Y - pr)
        D2logl <- c(weights * pr / (1 + exp(lp))) * X2
      } else {
        varp   <- pr / (1 + exp(lp))
        Dlogl  <- weights * X * c((Y - pr) * varp)
        D2logl <- c(weights * varp^2) * X2
      }
    } else stop("link functions must be logit, exp, lin\n")

    ploglik  <- sumstrata(weights * (Y - pr)^2, strata, nstrata)
    ploglik  <- matrix(ploglik, nstrata, 1)
    if (model[1] != "lin") ploglik <- 0 * ploglik
    D2log    <- Msumstrata(D2logl,  strata, nstrata)
    gradient <- Msumstrata(Dlogl,   strata, nstrata) + augmentation
    if (!is.matrix(D2log))    D2log    <- matrix(D2log,    1, length(D2log))
    if (!is.matrix(gradient)) gradient <- matrix(gradient, 1, np)

    gradient <- split(-gradient, seq_len(nrow(gradient)))
    ploglik  <- split(-ploglik,  seq_len(nrow(ploglik)))
    D2log    <- split(D2log,     seq_len(nrow(D2log)))

    sym <- function(x)
      matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
                   PACKAGE="mets")$XXf, np, np)
    hessian <- lapply(D2log, sym)

    if (all) {
      ihess    <- lapply(hessian, pinv)
      iid_list <- lapply(seq_len(nstrata), function(s) {
        idx  <- strata == (s - 1)
        id_s <- id[idx]
        Dl_s <- Dlogl[idx, , drop=FALSE]
        raw  <- Dl_s %*% ihess[[s]]
        Msumstrata(raw, id_s, max(id) + 1)
      })
      beta.iid <- do.call(cbind, iid_list)
      robvar   <- crossprod(beta.iid)
      return(list(
        par      = pp,           coef      = unlist(pp),
        ploglik  = ploglik,      gradient  = gradient,
        hessian  = hessian,      ihessian  = ihess,
        id       = id,           Dlogl     = Dlogl,
        iid      = beta.iid,     robvar    = robvar,
        var      = robvar,       se.robust = diag(robvar)^.5
      ))
    }
    structure(ploglik, gradient = gradient, hessian = hessian)
  } # }}}

  ## ------------------------------------------------------------------
  ## make_stratum_obj: per-stratum closure for NR
  ## cv.fold: use complement rows; standard: use stratum rows
  ## ------------------------------------------------------------------
  make_stratum_obj <- function(s, pp) { # {{{
    idx       <- if (cv.fold) strata != s else strata == s
    X_s       <- X[idx,  , drop=FALSE]
    X_s2      <- X2[idx, , drop=FALSE]
    Y_s       <- Y[idx]
    weights_s <- weights[idx]
    offset_s  <- offset[idx]
    aug_s     <- augmentation[s + 1, ]

    sym <- function(x, np)
      matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
                   PACKAGE="mets")$XXf, np, np)

    function(theta) {
      lp <- as.vector(X_s %*% theta) + offset_s
      if (model[1] == "exp") {
        pr <- exp(lp)
        if (monotone) {
          Dlogl  <- weights_s * X_s * c(Y_s - pr)
          D2logl <- c(weights_s * pr) * X_s2
        } else {
          Dlogl  <- weights_s * X_s * c((Y_s - pr) * pr)
          D2logl <- c(weights_s * pr^2) * X_s2
        }
        ploglik <- 0
      } else if (model[1] == "lin") {
        pr      <- lp
        Dlogl   <- weights_s * X_s * c(Y_s - pr)
        D2logl  <- c(weights_s) * X_s2
        ploglik <- sum(weights_s * (Y_s - pr)^2)
      } else if (model[1] == "logit") {
        pr <- expit(lp)
        if (monotone) {
          Dlogl  <- weights_s * X_s * c(Y_s - pr)
          D2logl <- c(weights_s * pr / (1 + exp(lp))) * X_s2
        } else {
          varp   <- pr / (1 + exp(lp))
          Dlogl  <- weights_s * X_s * c((Y_s - pr) * varp)
          D2logl <- c(weights_s * varp^2) * X_s2
        }
        ploglik <- 0
      } else stop("link functions must be logit, exp, lin\n")

      np   <- length(theta)
      grad <- -(colSums(Dlogl) + aug_s)
      D2   <- colSums(D2logl)
      hess <- sym(D2, np)
      structure(-ploglik, gradient = grad, hessian = hess)
    }
  } # }}}

  ## ------------------------------------------------------------------
  ## fast_lin_fit(): closed-form path for model=="lin". No per-stratum
  ## subsetting of X/X2/Y/weights -- one grouped pass gives every
  ## stratum's own-sums; cv.fold complements are total - own; iid
  ## residuals for all strata computed via one n x nstrata matmul.
  ## ------------------------------------------------------------------
  fast_lin_fit <- function() { ## {{{ 
    np  <- ncol(X)
    wX2 <- c(weights) * X2
    wXY <- c(weights) * X * c(Y)

    own_XtWX <- Msumstrata(wX2, strata, nstrata)
    own_XtWY <- Msumstrata(wXY, strata, nstrata)
    tot_XtWX <- colSums(wX2)
    tot_XtWY <- colSums(wXY)

    if (!is.matrix(own_XtWX)) own_XtWX <- matrix(own_XtWX, 1, length(own_XtWX))
    if (!is.matrix(own_XtWY)) own_XtWY <- matrix(own_XtWY, 1, np)

    sym <- function(v) matrix(.Call("XXMatFULL", matrix(v, nrow=1), np,
                                     PACKAGE="mets")$XXf, np, np)

    hessian <- ihess <- beta_opt <- vector("list", nstrata)
    for (s in seq_len(nstrata)) {
      vecXtWX <- if (cv.fold) tot_XtWX - own_XtWX[s, ] else own_XtWX[s, ]
      XtWY    <- if (cv.fold) tot_XtWY - own_XtWY[s, ] else own_XtWY[s, ]
      XtWY    <- XtWY + augmentation[s, ]
      H       <- sym(vecXtWX)
      iH      <- pinv(H)
      hessian[[s]]  <- H
      ihess[[s]]    <- iH
      beta_opt[[s]] <- as.vector(iH %*% XtWY)
    }

    ## residuals for ALL strata betas computed at once: n x nstrata
    BetaMat <- do.call(cbind, beta_opt)
    LP      <- X %*% BetaMat
    R       <- c(Y) - LP

    mid <- max(id)
    iid_list <- vector("list", nstrata)
    for (s in seq_len(nstrata)) {
      mask <- if (cv.fold) (strata != (s - 1)) else (strata == (s - 1))
      Dl_s <- weights * X * (R[, s] * mask)
      raw  <- Dl_s %*% ihess[[s]]
      iid_list[[s]] <- Msumstrata(raw, id, mid + 1)
    }
    beta.iid <- do.call(cbind, iid_list)
    robvar   <- crossprod(beta.iid)

    list(par = beta_opt, coef = unlist(beta_opt),
         ploglik  = lapply(seq_len(nstrata), function(s) 0),
         gradient = lapply(beta_opt, function(b) rep(0, length(b))),
         hessian = hessian, ihessian = ihess,
         id = id, Dlogl = NULL,
         iid = beta.iid, robvar = robvar, var = robvar,
         se.robust = diag(robvar)^.5)
  } ## }}} 

  ## ------------------------------------------------------------------
  ## Optimise
  ## ------------------------------------------------------------------
###  dots    <- list(...)
###  control <- if (length(dots) == 0) {
###    if (model[1] == "exp") list(tol=1e-10, stepsize=0.5) else NULL
###  } else dots[[1]]
###
###  if (!is.matrix(beta)) beta <- matrix(beta, nstrata, p, byrow = TRUE)
###  beta <- split(beta, seq_len(nrow(beta)))
###
###  NR_list <- function() {
###    Map(function(s) {
###      lava::NR(start   = beta[[s + 1]],
###                obj     = make_stratum_obj(s, beta),
###                control = control)
###    }, seq_len(nstrata) - 1L)
###  }
###
###  beta_opt <- beta   # fallback
###
###  if (p > 0) {
###    if (!no.opt) {
###      if (tolower(method) == "nr") {
###        tim      <- NR_list()
###        beta_opt <- lapply(tim, function(r) r$par)
###        cc       <- unlist(beta_opt)
###        val      <- c(list(coef = cc), obj(beta_opt, all = TRUE))
###        val$gradient <- lapply(tim, function(r) r$gradient)
###      } else {
###        stop("only NR of lava\n")
###      }
###    } else {
###      val <- c(list(coef = unlist(beta)), obj(beta, all = TRUE))
###    }
###  } else {
###    val <- obj(as.list(rep(0, nstrata)), all = TRUE)
###  }

## ------------------------------------------------------------------
  ## Optimise
  ## ------------------------------------------------------------------
  dots    <- list(...)
  control <- if (length(dots) == 0) {
    if (model == "exp") list(tol=1e-10, stepsize=0.5) else NULL
  } else dots[[1]]
  if (!is.matrix(beta)) beta <- matrix(beta, nstrata, p, byrow = TRUE)
  beta <- split(beta, seq_len(nrow(beta)))
  NR_list <- function() {
    Map(function(s) {
      lava::NR(start   = beta[[s + 1]],
                obj     = make_stratum_obj(s, beta),
                control = control)
    }, seq_len(nstrata) - 1L)
  }
  beta_opt <- beta   # fallback
  if (p > 0) {
    if (!no.opt) {
      if (model == "lin") {
        val      <- fast_lin_fit()
        beta_opt <- val$par
      } else if (tolower(method) == "nr") {
        tim      <- NR_list()
        beta_opt <- lapply(tim, function(r) r$par)
        cc       <- unlist(beta_opt)
        val      <- c(list(coef = cc), obj(beta_opt, all = TRUE))
        val$gradient <- lapply(tim, function(r) r$gradient)
      } else {
        stop("only NR of lava\n")
      }
    } else {
      if (model == "lin") {
        val <- fast_lin_fit()   # still uses supplied beta? no -- no.opt means
                                 # skip optimisation; keep supplied beta instead:
        val$par  <- beta
        val$coef <- unlist(beta)
        beta_opt <- beta
      } else {
        val <- c(list(coef = unlist(beta)), obj(beta, all = TRUE))
      }
    }
  } else {
    val <- if (model == "lin") fast_lin_fit() else obj(as.list(rep(0, nstrata)), all = TRUE)
  }




  ## ------------------------------------------------------------------
  ## Coefficient naming
  ## ------------------------------------------------------------------
  if (nstrata == 1 && length(val$coef) == length(colnames(X))) {
    names(val$coef) <- colnames(X)
  } else {
    cnames <- paste(
      rep(paste0("s", 1:nstrata), each = ncol(X)),
      rep(colnames(X), times = nstrata),
      sep = "."
    )
    names(val$coef)   <- cnames
    colnames(val$iid) <- cnames
  }

  ## ------------------------------------------------------------------
  ## Assemble output
  ## ------------------------------------------------------------------
  val <- c(val, list(
    time          = time,          formula      = formula,
    formC         = formC,         cens.weights = cens.weights,
    cens.strata   = cens.strata,   cens.nstrata = cens.nstrata,
    n             = nrow(X),  nevent       = nevent,
    ncluster      = nid
  ))
  val$call      <- cl
  val$MGciid    <- MGCiid
  val$call.id   <- call.id
  val$id        <- orig.id
  val$name.id   <- name.id
  val$nid       <- nid
  val$iid.naive <- val$iid
  val$naive.var <- NULL
  val$coef_list <- beta_opt

  ## ------------------------------------------------------------------
  ## Censoring IID adjustment, per stratum block
  ## iid[, cols_s] += MGCiid[, cols_s] %*% H_s^{-1}
  ## For cv.fold, MGCiid[, cols_s] already encodes complement rows via Xs,
  ## and ihessian[[s]] is from the complement fit
  ## ------------------------------------------------------------------
  if (se) {
    val$iid <- do.call(cbind, lapply(seq_len(nstrata), function(s) {
      cols <- (s - 1) * p + seq_len(p)
      val$iid[, cols] + MGCiid[, cols] %*% val$ihessian[[s]]
    }))
  }

  if (!is.null(name.id)) val$iid <- nameme(val$iid, name.id)
  robvar        <- crossprod(val$iid)
  val$var       <- val$robvar <- robvar
  val$se.robust <- diag(robvar)^.5
  val$se.coef   <- diag(val$var)^.5

  val$cause        <- cause
  val$cens.code    <- cens.code
  val$augmentation <- augmentation
  val$model        <- model[1]
  val$outcome      <- outcome[1]
  val$Yipcw        <- Yipcw
  val$Y            <- Yresp
  val$IPCW         <- IPCW
  val$Causes       <- Causes
  val$nevent       <- nevent
  val$design       <- des
  val$nstrata      <- nstrata
  val$strata_orig  <- strata
  val$strata_call  <- des.strata
  val$cv.fold      <- cv.fold

  if (low.memory) {
	val$design$y <- NULL
	val$design$data <- NULL
	val$design$x <- NULL
	val$design$strata <- NULL
	val$strata_orig <- NULL
	val$strata_call <- NULL
	val$MGciid <- NULL
	val$iid.naive <- NULL
	val$iid.naive <- NULL
	val$cens.strata <-  NULL
	val$cens.weights <-  NULL
	val$name.id  <-  NULL
  }

  class(val) <- c("binregStrata","binreg")
  return(val)
} # }}}


###binregStrata <- function(formula, data, cause=1, time=NULL, beta=NULL,
###                          type=c("II","I"),
###                          offset=NULL, weights=NULL, cens.weights=NULL,
###                          cens.model=~+1, se=TRUE,
###                          kaplan.meier=TRUE, cens.code=0, no.opt=FALSE,
###                          method="nr", augmentation=NULL,
###                          outcome=c("cif","rmst","rmtl"),
###                          model=c("default","logit","exp","lin"),
###                          Ydirect=NULL, strata=NULL, cv.fold=FALSE, low.memory=FALSE,...)
###{ # {{{
###  monotone <- TRUE
###  cl       <- match.call()
###
###  ## ------------------------------------------------------------------
###  ## Design
###  ## ------------------------------------------------------------------
###  des <- proc_design(
###    formula, data = data,
###    specials = c("offset","weights","cluster","strata"),
###    intercept = TRUE
###  )
###  Y <- des$y
###  if (!inherits(Y, c("Event","Surv")))
###    stop("Expected a 'Surv' or 'Event'-object")
###
###  if (ncol(Y) == 2) {
###    exit  <- Y[,1]; entry <- rep(0, nrow(Y)); status <- Y[,2]
###  } else {
###    entry <- Y[,1]; exit  <- Y[,2];           status <- Y[,3]
###  }
###
###  X           <- des$x
###  des.weights <- des$weights
###  des.offset  <- des$offset
###  id          <- des$cluster
###
###  if (ncol(X) == 0) X <- matrix(nrow=0, ncol=0)
###
###  call.id <- id
###  conid   <- construct_id(id, nrow(X), namesX = rownames(X))
###  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
###  orig.id <- id
###
###  ## ------------------------------------------------------------------
###  ## Weights / offset  (explicit arg > formula-embedded special > default)
###  ## ------------------------------------------------------------------
###  if (!is.null(offset)) {
###    offset <- offset
###  } else if (!is.null(des.offset)) {
###    offset <- des.offset
###  } else {
###    offset <- rep(0, length(exit))
###  }
###
###  if (!is.null(weights)) {
###    weights <- weights
###  } else if (!is.null(des.weights)) {
###    weights <- des.weights
###  } else {
###    weights <- rep(1, length(exit))
###  }
###
###  ## ------------------------------------------------------------------
###  ## Strata  (explicit arg > formula-embedded special > default)
###  ## ------------------------------------------------------------------
###  des.strata <- des$strata
###  if (!is.null(strata)) {
###    if (!is.numeric(strata)) stop("strata must be numeric\n")
###    ustrata <- sort(unique(strata)); nstrata <- length(ustrata)
###    strata  <- fast.approx(ustrata, strata) - 1
###  } else if (!is.null(des.strata)) {
###    strata  <- des.strata
###    ustrata <- sort(unique(strata)); nstrata <- length(ustrata)
###    strata  <- fast.approx(ustrata, strata) - 1
###  } else {
###    strata  <- rep(0, nrow(X)); nstrata <- 1
###  }
###
###
###  ## ------------------------------------------------------------------
###  ## Event / censoring indicators
###  ## ------------------------------------------------------------------
###  statusC  <- 1*(status %in% cens.code)
###  statusE  <- (status %in% cause) & (exit <= time)
###  if (sum(statusE) == 0) warning("No events of type 1 before time\n")
###
###  kmt     <- kaplan.meier
###  ucauses <- sort(unique(status))
###  ccc     <- which(ucauses %in% cens.code)
###  Causes  <- if (length(ccc) == 0) ucauses else ucauses[-ccc]
###  competing <- any(!(Causes %in% cause))
###
###  data$id__    <- id
###  data$exit    <- exit
###  data$statusC <- statusC
###
###  cens.strata <- cens.nstrata <- NULL
###  nevent      <- sum((status %in% cause) * (exit <= time))
###  obs         <- (exit <= time & statusC == 0) | (exit >= time)
###
###  ## ------------------------------------------------------------------
###  ## Censoring weights
###  ## ------------------------------------------------------------------
###  if (is.null(cens.weights)) { # {{{
###    formC    <- update.formula(cens.model, Surv(exit,statusC) ~ . + cluster(id__))
###    resC     <- phreg(formC, data)
###    if (resC$p > 0) kmt <- FALSE
###    exittime <- pmin(exit, time)
###    cens.weights <- suppressWarnings(
###      predict(resC, data, times = exittime, individual.time = TRUE,
###              se = FALSE, km = kmt, tminus = TRUE)$surv)
###    cens.strata  <- resC$strata[order(resC$ord)]
###    cens.nstrata <- resC$nstrata
###  } else {
###    se <- FALSE; resC <- formC <- NULL
###  } # }}}
###
###  ## ------------------------------------------------------------------
###  ## Response
###  ## ------------------------------------------------------------------
###  expit <- function(z) 1 / (1 + exp(-z))
###  p     <- ncol(X)
###
###  if (is.null(beta))         beta         <- rep(0, p)
###  if (is.null(augmentation)) augmentation <- matrix(0, nstrata, p)
###  if (!is.matrix(augmentation))
###    augmentation <- matrix(augmentation, nstrata, p, byrow = TRUE)
###
###  X  <- as.matrix(X)
###  X2 <- .Call("vecCPMat", X)$XX
###
###  if (!is.null(Ydirect)) {
###    Y <- Ydirect
###  } else {
###    if (outcome[1] == "cif") {
###      Y <- c((status %in% cause) * (exit <= time))
###    } else {
###      if (!competing) {
###        if (outcome[1] == "rmst") Y <- c(pmin(exit, time))
###        else                      Y <- c(time - pmin(exit, time))
###      } else {
###        Y <- c((status %in% cause) * (time - pmin(exit, time)))
###      }
###    }
###  }
###  Yresp <- c(Y)
###  Y     <- Yipcw <- c(Y * obs / cens.weights)
###  IPCW  <- c(obs / cens.weights)
###
###  ## ------------------------------------------------------------------
###  ## MGCiid: censoring IID influence
###  ## cv.fold: block s populated with complement rows (strata != s)
###  ## standard: block s populated with stratum-s rows
###  ## ------------------------------------------------------------------
###  MGCiid <- 0
###  if (se & monotone) { # {{{
###    ord     <- resC$ord
###    X       <- X[ord,  , drop = FALSE]
###    X2      <- X2[ord, , drop = FALSE]
###    status  <- status[ord]
###    exit    <- exit[ord]
###    weights <- weights[ord]
###    offset  <- offset[ord]
###    Y       <- Y[ord]
###    strata  <- strata[ord]
###
###    xx   <- resC$cox.prep
###    S0i2 <- S0i <- rep(0, length(xx$strata))
###    S0i[xx$jumps  + 1] <- 1 / resC$S0
###    S0i2[xx$jumps + 1] <- 1 / resC$S0^2
###    btime <- 1*(exit < time)
###    mid   <- max(id)
###
###    ## Block-expanded design: n x (nstrata*p)
###    ## cv.fold: block s gets complement rows; standard: block s gets stratum-s rows
###    Xs <- matrix(0, nrow(X), nstrata * p)
###    for (s in 0:(nstrata - 1)) {
###      idx  <- if (cv.fold) which(strata != s) else which(strata == s)
###      cols <- s * p + seq_len(p)
###      Xs[idx, cols] <- X[idx, ]
###    }
###
###    h           <- Mrevcumsumstrata(Xs * Y, xx$strata, xx$nstrata)
###    IhdLam0     <- Mcumsumstrata(h * S0i2 * btime, xx$strata, xx$nstrata)
###    U           <- matrix(0, nrow(xx$X), ncol(Xs))
###    U[xx$jumps + 1,] <- (resC$jumptimes < time) *
###                         h[xx$jumps + 1,] / c(resC$S0)
###    MGt    <- (U - IhdLam0) * c(xx$weights)
###    MGCiid <- Msumstrata(MGt, xx$id, mid + 1)
###
###    if (type[1] == "II") { # {{{
###      hYt      <- revcumsumstrata(Y,   xx$strata, xx$nstrata)
###      IhdLam0  <- cumsumstrata(hYt * S0i2 * btime, xx$strata, xx$nstrata)
###      U2       <- rep(0, length(xx$strata))
###      U2[xx$jumps + 1] <- (resC$jumptimes < time) *
###                           hYt[xx$jumps + 1] / c(resC$S0)
###      MGt2    <- Xs * c(U2 - IhdLam0) * c(xx$weights)
###      MGtiid  <- Msumstrata(MGt2, xx$id, mid + 1)
###      augmentation <- colSums(MGtiid) + augmentation
###      augmentation <- matrix(augmentation, nstrata, p, byrow = TRUE)
###
###      EXt          <- Mrevcumsumstrata(Xs, xx$strata, xx$nstrata)
###      IEXhYtdLam0  <- Mcumsumstrata(EXt * c(hYt) * S0i * S0i2 * btime,
###                                     xx$strata, xx$nstrata)
###      U3 <- matrix(0, nrow(xx$X), ncol(Xs))
###      U3[xx$jumps + 1,] <- (resC$jumptimes < time) *
###                            hYt[xx$jumps + 1] *
###                            EXt[xx$jumps + 1,] / c(resC$S0)^2
###      MGt3    <- (U3 - IEXhYtdLam0) * c(xx$weights)
###      MGCiid2 <- Msumstrata(MGt3, xx$id, mid + 1)
###      MGCiid  <- MGCiid + (MGtiid - MGCiid2)
###    } # }}}
###
###    id <- xx$id
###  } else {
###    MGCiid <- 0
###  } # }}}
###
###  ## ------------------------------------------------------------------
###  ## Model default
###  ## ------------------------------------------------------------------
###  if (model[1] == "default") {
###    if (outcome[1] == "cif")        model <- "logit"
###    if (outcome[1] == "rmst")       model <- "exp"
###    if (outcome[1] == "rmtl")       model <- "exp"
###    if (outcome[1] == "years-lost") model <- "exp"
###  }
###
###  ## ------------------------------------------------------------------
###  ## obj(): full objective (used for final all=TRUE call and NR)
###  ## ------------------------------------------------------------------
###  obj <- function(pp, all = FALSE) { # {{{
###    np <- length(pp[[1]])
###
###    ## ----------------------------------------------------------------
###    ## cv.fold path: for stratum s, beta_s estimated on complement rows
###    ## ----------------------------------------------------------------
###    if (cv.fold) {
###      ploglik_list  <- vector("list", nstrata)
###      gradient_list <- vector("list", nstrata)
###      D2log_list    <- vector("list", nstrata)
###
###      for (s in 0:(nstrata - 1)) {
###        idx_c <- which(strata != s)          # complement indices
###        X_c   <- X[idx_c,  , drop=FALSE]
###        X2_c  <- X2[idx_c, , drop=FALSE]
###        Y_c   <- Y[idx_c]
###        w_c   <- weights[idx_c]
###        off_c <- offset[idx_c]
###        aug_s <- augmentation[s + 1, ]
###        lp_c  <- as.vector(X_c %*% pp[[s + 1]]) + off_c
###
###        if (model[1] == "exp") {
###          pr_c  <- exp(lp_c)
###          if (monotone) {
###            Dl_c  <- w_c * X_c * c(Y_c - pr_c)
###            D2l_c <- c(w_c * pr_c) * X2_c
###          } else {
###            Dl_c  <- w_c * X_c * c((Y_c - pr_c) * pr_c)
###            D2l_c <- c(w_c * pr_c^2) * X2_c
###          }
###          pl_c <- 0
###        } else if (model[1] == "lin") {
###          pr_c  <- lp_c
###          Dl_c  <- w_c * X_c * c(Y_c - pr_c)
###          D2l_c <- c(w_c) * X2_c
###          pl_c  <- sum(w_c * (Y_c - pr_c)^2)
###        } else if (model[1] == "logit") {
###          pr_c <- expit(lp_c)
###          if (monotone) {
###            Dl_c  <- w_c * X_c * c(Y_c - pr_c)
###            D2l_c <- c(w_c * pr_c / (1 + exp(lp_c))) * X2_c
###          } else {
###            varp  <- pr_c / (1 + exp(lp_c))
###            Dl_c  <- w_c * X_c * c((Y_c - pr_c) * varp)
###            D2l_c <- c(w_c * varp^2) * X2_c
###          }
###          pl_c <- 0
###        } else stop("link must be logit, exp, or lin\n")
###
###        ploglik_list[[s + 1]]  <- pl_c
###        gradient_list[[s + 1]] <- -(colSums(Dl_c) + aug_s)
###        D2log_list[[s + 1]]    <- colSums(D2l_c)
###      }
###
###      sym <- function(x)
###        matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
###                     PACKAGE="mets")$XXf, np, np)
###      hessian <- lapply(D2log_list, sym)
###
###      if (all) {
###        ihess <- lapply(hessian, pinv)
###
###        ## Per-stratum IID: recompute Dl_c at optimal pp[[s]] separately
###        ## for each stratum so blocks are independent and correctly sized.
###        ## Obs i in stratum r contributes to beta_s for all s != r;
###        ## Msumstrata fills zeros for stratum-s obs (they don't contribute).
###        beta.iid <- do.call(cbind, lapply(seq_len(nstrata), function(s) {
###          idx_c <- which(strata != (s - 1))
###          X_c   <- X[idx_c,  , drop=FALSE]
###          Y_c   <- Y[idx_c]
###          w_c   <- weights[idx_c]
###          off_c <- offset[idx_c]
###          id_c  <- id[idx_c]
###          lp_c  <- as.vector(X_c %*% pp[[s]]) + off_c
###
###          Dl_c <- switch(model[1],
###            exp = {
###              pr_c <- exp(lp_c)
###              if (monotone) w_c * X_c * c(Y_c - pr_c)
###              else          w_c * X_c * c((Y_c - pr_c) * pr_c)
###            },
###            lin = w_c * X_c * c(Y_c - lp_c),
###            logit = {
###              pr_c <- expit(lp_c)
###              if (monotone) w_c * X_c * c(Y_c - pr_c)
###              else {
###                varp <- pr_c / (1 + exp(lp_c))
###                w_c * X_c * c((Y_c - pr_c) * varp)
###              }
###            },
###            stop("unknown model\n")
###          )
###
###          raw <- Dl_c %*% ihess[[s]]
###          Msumstrata(raw, id_c, max(id) + 1)
###        }))
###
###        robvar <- crossprod(beta.iid)
###        return(list(
###          par      = pp,             coef      = unlist(pp),
###          ploglik  = ploglik_list,   gradient  = gradient_list,
###          hessian  = hessian,        ihessian  = ihess,
###          id       = id,             Dlogl     = NULL,
###          iid      = beta.iid,       robvar    = robvar,
###          var      = robvar,         se.robust = diag(robvar)^.5
###        ))
###      }
###
###      return(structure(ploglik_list,
###                       gradient = gradient_list,
###                       hessian  = hessian))
###    }
###
#### ----------------------------------------------------------------
#### Standard non-CV
#### ----------------------------------------------------------------
###    loffset <- eta <- numeric(length(Y))
###    for (s in 0:(nstrata - 1)) {
###      idx          <- strata == s
###      eta[idx]     <- as.vector(X[idx, , drop=FALSE] %*% pp[[s + 1]])
###      loffset[idx] <- offset[idx]
###    }
###    lp <- eta + loffset
###
###    if (model[1] == "exp") {
###      pr     <- exp(lp)
###      if (monotone) {
###        Dlogl  <- weights * X * c(Y - pr)
###        D2logl <- c(weights * pr) * X2
###      } else {
###        Dlogl  <- weights * X * c((Y - pr) * pr)
###        D2logl <- c(weights * pr^2) * X2
###      }
###    } else if (model[1] == "lin") {
###      pr     <- lp
###      Dlogl  <- weights * X * c(Y - pr)
###      D2logl <- c(weights) * X2
###    } else if (model[1] == "logit") {
###      pr <- expit(lp)
###      if (monotone) {
###        Dlogl  <- weights * X * c(Y - pr)
###        D2logl <- c(weights * pr / (1 + exp(lp))) * X2
###      } else {
###        varp   <- pr / (1 + exp(lp))
###        Dlogl  <- weights * X * c((Y - pr) * varp)
###        D2logl <- c(weights * varp^2) * X2
###      }
###    } else stop("link functions must be logit, exp, lin\n")
###
###    ploglik  <- sumstrata(weights * (Y - pr)^2, strata, nstrata)
###    ploglik  <- matrix(ploglik, nstrata, 1)
###    if (model[1] != "lin") ploglik <- 0 * ploglik
###    D2log    <- Msumstrata(D2logl,  strata, nstrata)
###    gradient <- Msumstrata(Dlogl,   strata, nstrata) + augmentation
###    if (!is.matrix(D2log))    D2log    <- matrix(D2log,    1, length(D2log))
###    if (!is.matrix(gradient)) gradient <- matrix(gradient, 1, np)
###
###    gradient <- split(-gradient, seq_len(nrow(gradient)))
###    ploglik  <- split(-ploglik,  seq_len(nrow(ploglik)))
###    D2log    <- split(D2log,     seq_len(nrow(D2log)))
###
###    sym <- function(x)
###      matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
###                   PACKAGE="mets")$XXf, np, np)
###    hessian <- lapply(D2log, sym)
###
###    if (all) {
###      ihess    <- lapply(hessian, pinv)
###      iid_list <- lapply(seq_len(nstrata), function(s) {
###        idx  <- strata == (s - 1)
###        id_s <- id[idx]
###        Dl_s <- Dlogl[idx, , drop=FALSE]
###        raw  <- Dl_s %*% ihess[[s]]
###        Msumstrata(raw, id_s, max(id) + 1)
###      })
###      beta.iid <- do.call(cbind, iid_list)
###      robvar   <- crossprod(beta.iid)
###      return(list(
###        par      = pp,           coef      = unlist(pp),
###        ploglik  = ploglik,      gradient  = gradient,
###        hessian  = hessian,      ihessian  = ihess,
###        id       = id,           Dlogl     = Dlogl,
###        iid      = beta.iid,     robvar    = robvar,
###        var      = robvar,       se.robust = diag(robvar)^.5
###      ))
###    }
###    structure(ploglik, gradient = gradient, hessian = hessian)
###  } # }}}
###
###  ## ------------------------------------------------------------------
###  ## make_stratum_obj: per-stratum closure for NR
###  ## cv.fold: use complement rows; standard: use stratum rows
###  ## ------------------------------------------------------------------
###  make_stratum_obj <- function(s, pp) { # {{{
###    idx       <- if (cv.fold) strata != s else strata == s
###    X_s       <- X[idx,  , drop=FALSE]
###    X_s2      <- X2[idx, , drop=FALSE]
###    Y_s       <- Y[idx]
###    weights_s <- weights[idx]
###    offset_s  <- offset[idx]
###    aug_s     <- augmentation[s + 1, ]
###
###    sym <- function(x, np)
###      matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
###                   PACKAGE="mets")$XXf, np, np)
###
###    function(theta) {
###      lp <- as.vector(X_s %*% theta) + offset_s
###      if (model[1] == "exp") {
###        pr <- exp(lp)
###        if (monotone) {
###          Dlogl  <- weights_s * X_s * c(Y_s - pr)
###          D2logl <- c(weights_s * pr) * X_s2
###        } else {
###          Dlogl  <- weights_s * X_s * c((Y_s - pr) * pr)
###          D2logl <- c(weights_s * pr^2) * X_s2
###        }
###        ploglik <- 0
###      } else if (model[1] == "lin") {
###        pr      <- lp
###        Dlogl   <- weights_s * X_s * c(Y_s - pr)
###        D2logl  <- c(weights_s) * X_s2
###        ploglik <- sum(weights_s * (Y_s - pr)^2)
###      } else if (model[1] == "logit") {
###        pr <- expit(lp)
###        if (monotone) {
###          Dlogl  <- weights_s * X_s * c(Y_s - pr)
###          D2logl <- c(weights_s * pr / (1 + exp(lp))) * X_s2
###        } else {
###          varp   <- pr / (1 + exp(lp))
###          Dlogl  <- weights_s * X_s * c((Y_s - pr) * varp)
###          D2logl <- c(weights_s * varp^2) * X_s2
###        }
###        ploglik <- 0
###      } else stop("link functions must be logit, exp, lin\n")
###
###      np   <- length(theta)
###      grad <- -(colSums(Dlogl) + aug_s)
###      D2   <- colSums(D2logl)
###      hess <- sym(D2, np)
###      structure(-ploglik, gradient = grad, hessian = hess)
###    }
###  } # }}}
###
###  ## ------------------------------------------------------------------
###  ## Optimise
###  ## ------------------------------------------------------------------
###  dots    <- list(...)
###  control <- if (length(dots) == 0) {
###    if (model[1] == "exp") list(tol=1e-10, stepsize=0.5) else NULL
###  } else dots[[1]]
###
###  if (!is.matrix(beta)) beta <- matrix(beta, nstrata, p, byrow = TRUE)
###  beta <- split(beta, seq_len(nrow(beta)))
###
###  NR_list <- function() {
###    Map(function(s) {
###      lava::NR(start   = beta[[s + 1]],
###                obj     = make_stratum_obj(s, beta),
###                control = control)
###    }, seq_len(nstrata) - 1L)
###  }
###
###  beta_opt <- beta   # fallback
###
###  if (p > 0) {
###    if (!no.opt) {
###      if (tolower(method) == "nr") {
###        tim      <- NR_list()
###        beta_opt <- lapply(tim, function(r) r$par)
###        cc       <- unlist(beta_opt)
###        val      <- c(list(coef = cc), obj(beta_opt, all = TRUE))
###        val$gradient <- lapply(tim, function(r) r$gradient)
###      } else {
###        stop("only NR of lava\n")
###      }
###    } else {
###      val <- c(list(coef = unlist(beta)), obj(beta, all = TRUE))
###    }
###  } else {
###    val <- obj(as.list(rep(0, nstrata)), all = TRUE)
###  }
###
###  ## ------------------------------------------------------------------
###  ## Coefficient naming
###  ## ------------------------------------------------------------------
###  if (nstrata == 1 && length(val$coef) == length(colnames(X))) {
###    names(val$coef) <- colnames(X)
###  } else {
###    cnames <- paste(
###      rep(paste0("s", 1:nstrata), each = ncol(X)),
###      rep(colnames(X), times = nstrata),
###      sep = "."
###    )
###    names(val$coef)   <- cnames
###    colnames(val$iid) <- cnames
###  }
###
###  ## ------------------------------------------------------------------
###  ## Assemble output
###  ## ------------------------------------------------------------------
###  val <- c(val, list(
###    time          = time,          formula      = formula,
###    formC         = formC,         cens.weights = cens.weights,
###    cens.strata   = cens.strata,   cens.nstrata = cens.nstrata,
###    n             = length(exit),  nevent       = nevent,
###    ncluster      = nid
###  ))
###  val$call      <- cl
###  val$MGciid    <- MGCiid
###  val$call.id   <- call.id
###  val$id        <- orig.id
###  val$name.id   <- name.id
###  val$nid       <- nid
###  val$iid.naive <- val$iid
###  val$naive.var <- NULL
###  val$coef_list <- beta_opt
###
###  ## ------------------------------------------------------------------
###  ## Censoring IID adjustment, per stratum block
###  ## iid[, cols_s] += MGCiid[, cols_s] %*% H_s^{-1}
###  ## For cv.fold, MGCiid[, cols_s] already encodes complement rows via Xs,
###  ## and ihessian[[s]] is from the complement fit
###  ## ------------------------------------------------------------------
###  if (se) {
###    val$iid <- do.call(cbind, lapply(seq_len(nstrata), function(s) {
###      cols <- (s - 1) * p + seq_len(p)
###      val$iid[, cols] + MGCiid[, cols] %*% val$ihessian[[s]]
###    }))
###  }
###
###  if (!is.null(name.id)) val$iid <- nameme(val$iid, name.id)
###  robvar        <- crossprod(val$iid)
###  val$var       <- val$robvar <- robvar
###  val$se.robust <- diag(robvar)^.5
###  val$se.coef   <- diag(val$var)^.5
###
###  val$cause        <- cause
###  val$cens.code    <- cens.code
###  val$augmentation <- augmentation
###  val$model        <- model[1]
###  val$outcome      <- outcome[1]
###  val$Yipcw        <- Yipcw
###  val$Y            <- Yresp
###  val$IPCW         <- IPCW
###  val$Causes       <- Causes
###  val$nevent       <- nevent
###  val$design       <- des
###  val$nstrata      <- nstrata
###  val$strata_orig  <- strata
###  val$strata_call  <- des.strata
###  val$cv.fold      <- cv.fold
###
###  if (low.memory) {
###	val$design$y <- NULL
###	val$design$data <- NULL
###	val$design$x <- NULL
###	val$design$strata <- NULL
###	val$strata_orig <- NULL
###	val$strata_call <- NULL
###	val$MGciid <- NULL
###	val$iid.naive <- NULL
###	val$iid.naive <- NULL
###	val$cens.strata <-  NULL
###	val$cens.weights <-  NULL
###	val$name.id  <-  NULL
###  }
###
###  class(val) <- c("binregStrata","binreg")
###  return(val)
###} # }}}

## {{{  predict, 

## small helper: given flat coef vector, route predictions per stratum
strata_lp <- function(object, Z, strata_new, offset_new) {
  p         <- ncol(Z)
  beta_list <- split(object$coef, rep(seq_len(object$nstrata), each = p))
  eta       <- numeric(nrow(Z))
  for (s in seq_len(object$nstrata)) {
    idx       <- strata_new == (s - 1)
    eta[idx]  <- as.vector(Z[idx, , drop = FALSE] %*% beta_list[[s]]) + offset_new[idx]
  }
  eta
}

##' @export
predict.binregStrata <- function(object, newdata = NULL, se = TRUE,
                                  pred.iid = FALSE, level = 0.95,
                                  strata = NULL,offset=NULL, ...) {  ## {{{ 
  ## design matrix
  if (!is.null(newdata)) {
    x <- update_design(object$design, data = newdata, response = FALSE)
  } else {
    x <- object$design
  }
  Z          <- x$x

  if (!is.null(strata)) {
  strata_new <- as.numeric(strata) - 1
  } else if (!is.null(x$strata)) {
  strata_new <- as.numeric(x$strata) - 1
  } else {
  strata_new <- object$strata_orig
  }

  offset_new <- if (!is.null(offset)) {
  offset
  } else if (!is.null(x$offset)) {
  x$offset
  } else {
  rep(0, nrow(Z))
  }

  ## linear predictor routed per stratum
  lp <- strata_lp(object, Z, strata_new, offset_new)

  ## link
  p <- switch(object$model[1],
    logit = expit(lp),
    exp   = exp(lp),
    lin   = lp,
    stop("unknown model")
  )

## derivative for delta method: same shape as before, per-stratum coef block
  p_col <- ncol(Z)
  Dpv   <- matrix(0, nrow(Z), object$nstrata * p_col)
  for (s in seq_len(object$nstrata)) {
    idx  <- strata_new == (s - 1)
    cols <- (s - 1) * p_col + seq_len(p_col)
    if (object$model[1] == "logit") {
      Dpv[idx, cols] <- Z[idx, , drop = FALSE] * exp(-lp[idx]) * p[idx]^2
    } else if (object$model[1] == "exp") {
      Dpv[idx, cols] <- Z[idx, , drop = FALSE] * p[idx]
    } else {
      Dpv[idx, cols] <- Z[idx, , drop = FALSE]
    }
  }

  preds <- p
  if (se) {
    covv   <- object$var
    se_vec <- apply((Dpv %*% covv) * Dpv, 1, sum)^.5
    crit   <- qnorm(1 - (1 - level) / 2)
    preds  <- data.frame(pred  = p,
                         se    = se_vec,
                         lower = p - crit * se_vec,
                         upper = p + crit * se_vec)
  }
  if (pred.iid) {
    Piid  <- object$iid %*% t(Dpv)
    preds <- list(pred = preds, iid = Piid, id = x$cluster)
  }
  attr(preds, "id") <- x$cluster
  preds
} ## }}} 

## }}} 

## --------------------------------------------------------------------------
## summary.binregStrata
## --------------------------------------------------------------------------
##' @export
summary.binregStrata <- function(object, level = 0.95, ...) { ## {{{
  p       <- length(object$coef) / object$nstrata
  nms     <- colnames(object$design$x)
  n       <- nrow(object$iid)
  crit    <- qnorm(1 - (1 - level) / 2)

  exp_label <- switch(object$model[1],
    logit = "OR",
    exp   = "RR",
    lin   = NULL
  )

  strata_tables     <- vector("list", object$nstrata)
  strata_exp_tables <- vector("list", object$nstrata)

  for (s in seq_len(object$nstrata)) {
    cols  <- (s - 1L) * p + seq_len(p)
    b_s   <- object$coef[cols]
    iid_s <- object$iid[, cols, drop = FALSE]
    e_s   <- lava::estimate(NULL, coef = b_s, IC = iid_s * n)
    cm    <- e_s$coefmat
    rownames(cm) <- if (!is.null(nms)) nms else paste0("X", seq_len(p))
    strata_tables[[s]] <- cm

    if (!is.null(exp_label)) {
      e_exp <- lava::estimate(NULL, coef = b_s, IC = iid_s * n, f = exp)
      cm_exp <- e_exp$coefmat
      rownames(cm_exp) <- rownames(cm)
      strata_exp_tables[[s]] <- cm_exp
    }
  }
  names(strata_tables)     <- paste0("stratum", seq_len(object$nstrata))
  names(strata_exp_tables) <- paste0("stratum", seq_len(object$nstrata))

  structure(
    list(strata_tables     = strata_tables,
         strata_exp_tables = strata_exp_tables,
         exp_label         = exp_label,
         nstrata           = object$nstrata,
         p                 = p,
         model             = object$model[1],
         outcome           = object$outcome[1],
         n                 = object$n,
         nevent            = object$nevent,
         ncluster          = object$ncluster,
         cv.fold           = isTRUE(object$cv.fold),
         call              = object$call),
    class = "summary.binregStrata"
  )
} ## }}}

## --------------------------------------------------------------------------
## print.summary.binregStrata
## --------------------------------------------------------------------------
##' @export
print.summary.binregStrata <- function(x, digits = 3, ...) { ## {{{ 
  cat("\nCall:\n")
  if (!is.null(x$call)) print(x$call)
  cat("\nStratified binreg fit\n")
  cat("  Strata  :", x$nstrata,
      if (x$cv.fold) " [cv.fold -- leave-one-stratum-out CV]" else "", "\n")
  cat("  Model   :", x$model, "\n")
  cat("  Outcome :", x$outcome, "\n")
  cat("  n =", x$n,
      " | events =", x$nevent,
      " | clusters =", x$ncluster, "\n\n")

  for (s in seq_len(x$nstrata)) {
    cat(paste0("--- Stratum ", s,
               if (x$cv.fold)
                 paste0(" (fitted on strata != ", s, ")")
               else "",
               " ---\n"))
    printCoefmat(x$strata_tables[[s]], digits = digits, ...)

    if (!is.null(x$exp_label) && !is.null(x$strata_exp_tables[[s]])) {
      cat(x$exp_label, "= exp(coef):\n")
      em <- x$strata_exp_tables[[s]][, c("Estimate","2.5%","97.5%"), drop = FALSE]
      colnames(em)[1] <- x$exp_label
      print(round(em, digits))
    }
    cat("\n")
  }
  invisible(x)
} ## }}} 

## --------------------------------------------------------------------------
## print.binregStrata  -- brief one-liner per stratum
## --------------------------------------------------------------------------
##' @export
print.binregStrata <- function(x, digits = 3, ...) { ## {{{ 
  cat("\nStratified binreg fit\n")
  cat("  Model:", x$model[1], " | Outcome:", x$outcome[1],
      " | Strata:", x$nstrata,
      if (isTRUE(x$cv.fold)) " [cv.fold]" else "", "\n")
  cat("  n =", x$n,
      "| events =", x$nevent,
      "| clusters =", x$ncluster, "\n\n")

  p   <- length(x$coef) / x$nstrata
  nms <- colnames(x$design$x)
  n   <- nrow(x$iid)

  for (s in seq_len(x$nstrata)) {
    cols <- (s - 1L) * p + seq_len(p)
    b_s  <- x$coef[cols]
    se_s <- sqrt(diag(x$var)[cols])
    mat  <- cbind(Estimate = b_s, `Std.Err` = se_s,
                  `z`      = b_s / se_s,
                  `Pr(>|z|)` = 2 * pnorm(-abs(b_s / se_s)))
    if (!is.null(nms)) rownames(mat) <- nms
    cat(paste0("Stratum ", s,
               if (isTRUE(x$cv.fold))
                 paste0(" (train: strata != ", s, ")")
               else "", ":\n"))
    printCoefmat(mat, digits = digits, ...)
    cat("\n")
  }
  invisible(x)
} ## }}} 

##' @export
coef.binregStrata <- function(object, ...) object$coef

##' @export
vcov.binregStrata <- function(object, ...) object$var

##' @export
IC.binregStrata <- function(x, ...) x$iid * NROW(x$iid)


##' Cross-validated Brier score for competing risks, RMST and RMTL regression
##' using stratified leave-fold-out fits (\code{binreg} \code{binregStrata})
##'
##' Fits candidate models via \code{\link[mets]{binregStrata}} with
##' \code{cv.fold = TRUE} and evaluates each using apparent (in-sample) and
##' \emph{k}-fold cross-validated IPCW Brier scores relative to a null model.
##' Supports all outcome types available in \code{binreg}/\code{binregStrata}:
##' cause-specific cumulative incidence (competing risks), restricted mean
##' survival time (RMST), and restricted mean time lost (RMTL).
##'
##' Cross-validation is delegated directly to
##' \code{\link[mets]{binregStrata}} via its \code{strata} and
##' \code{cv.fold = TRUE} arguments: each stratum (fold) gets its own
##' coefficient vector \eqn{\beta_s}, estimated on the complement of that
##' fold. Out-of-fold predictions are obtained by routing each observation to
##' the \eqn{\beta_s} corresponding to its own fold, which was trained on the
##' other folds.
##'
##' The fold membership is stored in an internal \code{fold} column (1-based)
##' in \code{data} and passed to \code{binregStrata} as a 0-based
##' \code{strata} vector. Censoring weights for the underlying competing-risks
##' model are estimated according to \code{cens.model} as usual; since strata
##' already define the folds, \code{cens.model} need not include
##' \code{strata(fold)} explicitly.
##'
##' Brier scores are computed on the \strong{linear scale} via
##' \code{resmeanIPCW(..., model = "lin")} and then log-transformed using
##' \code{lava::estimate(..., f = log)} with the delta method. All stored
##' coefficients, standard errors, and influence functions are therefore on
##' the \strong{log scale}. Use \code{summary(fit, transform = exp)} or
##' \code{estimate(fit, f = exp)} to back-transform to the original Brier
##' scale.
##'
##' The delta-method correction for the cross-validated influence function
##' accounts for the block structure of \code{binregStrata}'s influence
##' function (\eqn{n_{id} \times (\code{nstrata} \times p)}, one
##' \eqn{p}-column block per stratum). For each stratum \eqn{s}, only
##' observations belonging to fold \eqn{s} have a non-zero derivative of
##' their out-of-fold prediction with respect to \eqn{\beta_s}, since
##' \eqn{\beta_s} is exactly the coefficient used to predict fold \eqn{s}.
##' The full gradient is therefore block-sparse, with one non-zero
##' \eqn{p}-length block per stratum, and the correction reduces to a single
##' matrix product against \code{iid(outff_fit)}.
##'
##' Left-truncation (delayed entry) is not supported and will raise an error.
##' Specifically, \code{Event(entry, time, status)} is rejected; use
##' \code{Event(time, status)} or \code{Event(time, status, cause)}.
##'
##' @param formula Two-sided formula of the form
##'   \code{Event(time, status) ~ covariates + cluster(id)} or
##'   \code{Event(time, status, cause) ~ covariates + cluster(id)}.
##'   The \code{cluster(id)} term is optional; if absent an internal
##'   identifier \code{id__} is added automatically assuming independent
##'   observations.
##' @param data A \code{data.frame} containing all variables in
##'   \code{formula} and \code{rhs}.
##' @param time Scalar or numeric vector of evaluation time points. When a
##'   vector is supplied the function loops over time points reusing the same
##'   fold structure and returns a named list of \code{binregCV} objects with
##'   class \code{"binregCV_list"}.
##' @param fold Integer. Number of cross-validation folds / strata
##'   (default 5).
##' @param rhs Named list of one-sided formulas specifying candidate models,
##'   e.g. \code{list(small = ~Z1, full = ~Z1 + Z2)}. If \code{NULL} (default),
##'   the right-hand side of \code{formula} (minus any \code{cluster()} term)
##'   is used as a single candidate named \code{"model1"}.
##' @param rhs0 One-sided formula for the null (intercept-only) model
##'   (default \code{~+1}).
##' @param cens.model Censoring model passed to \code{binregStrata} for the
##'   CV fold fits (default \code{~+1}). Since folds are handled internally
##'   by \code{binregStrata}'s \code{strata}/\code{cv.fold} mechanism, this
##'   does \emph{not} need to include \code{strata(fold)}.
##' @param brier.cens.model Censoring model passed to \code{resmeanIPCW} when
##'   computing Brier scores. Default \code{~+1} gives Kaplan-Meier weights.
##'   Supply a stratified model, e.g. \code{~strata(Z1, Z2)}, for
##'   covariate-adjusted IPCW weights.
##' @param brier.args Named list of additional arguments forwarded to every
##'   \code{resmeanIPCW} call (beyond \code{formula}, \code{data},
##'   \code{time}, \code{Ydirect}, \code{cens.model}, and \code{model} which
##'   are set internally). Typically empty.
##' @param cv.split.index List from \code{folds(n, fold)} or a vector of fold
##'   indices (1:fold) for each observation, to fix the cross-validation
##'   folds. \code{NULL} (default) generates folds internally via
##'   \code{folds(n, fold)}.
##' @param \dots Additional arguments passed to every \code{binreg} /
##'   \code{binregStrata} call, e.g. \code{cause = 2}, \code{outcome = "rmst"},
##'   or \code{offset}.
##'
##' @return An object of class \code{"binregCV"} (single time point) or
##'   \code{"binregCV_list"} (multiple time points), with the following
##'   structure:
##' \describe{
##'   \item{\code{null}}{List with \code{coef}, \code{se}, and \code{iid} for
##'     the null model. \code{coef} and \code{se} are named vectors with
##'     entries \code{bbrier} (apparent, log scale) and \code{bbrierCV}
##'     (cross-validated, log scale). \code{iid} contains the
##'     influence-function matrix for \code{bbrierCV} with delta-method
##'     correction for the stratified CV prediction step and the log
##'     transform.}
##'   \item{\code{models}}{Named list, one element per entry in \code{rhs}.
##'     Each element contains \code{coef} (\code{brier}, \code{brierCV}, log
##'     scale), \code{se}, \code{dcoef} (differences from null: log Brier
##'     ratios), \code{sed}, and \code{iid} (influence-function matrix for
##'     \code{brierCV}, log scale, with CV delta-method correction).}
##'   \item{\code{iid_mat}}{Matrix of dimension \eqn{n \times (1 +
##'     \code{length(rhs)})} with columns \code{null} and one per entry of
##'     \code{rhs}, each the log-scale CV Brier influence function for the
##'     corresponding model. Used by \code{\link{estimate.binregCV}} for joint
##'     inference.}
##'   \item{\code{cv.split.index}}{List of fold index vectors used.}
##'   \item{\code{rhs}, \code{rhs0}}{Processed candidate and null formulas
##'     (cluster terms removed).}
##'   \item{\code{model}}{Link function read from the \code{binreg} fit:
##'     \code{"logit"}, \code{"exp"}, or \code{"lin"}.}
##'   \item{\code{call}, \code{time}, \code{fold}, \code{cens.model},
##'     \code{brier.cens.model}}{Call and settings for printing and
##'     downstream use.}
##' }
##'
##' @section Outcome types:
##' The function supports all outcome types provided by \code{binreg}/
##' \code{binregStrata}:
##' \describe{
##'   \item{Competing risks CIF}{Default; use
##'     \code{Event(time, status, cause)} and pass \code{cause} via \code{...}.
##'     The Brier score measures squared error for the cause-specific CIF.}
##'   \item{RMST}{Pass \code{outcome = "rmst"} via \code{...}. The Brier score
##'     measures squared error for the restricted mean survival time up to
##'     \code{time}.}
##'   \item{RMTL}{Pass \code{outcome = "rmtl"} via \code{...}. The Brier score
##'     measures squared error for the restricted mean time lost.}
##' }
##'
##' @section Scale and inference:
##' All output is on the log scale. To obtain Brier scores on the original
##' scale with confidence intervals, use:
##' \preformatted{
##'   summary(fit, transform = exp)
##'   estimate(fit, f = exp)
##' }
##' Differences on the log scale correspond to log Brier ratios; exponentiate
##' to obtain Brier ratios. For absolute differences on the original scale,
##' compute contrasts after \code{estimate(fit, f = exp)}.
##'
##' @seealso
##'   \code{\link{estimate.binregCV}} for joint inference across models and
##'   time points;
##'   \code{\link{summary.binregCV}} and \code{\link{print.summary_binregCV_multi}}
##'   for formatted output;
##'   \code{\link[mets]{binregStrata}} for the underlying stratified
##'   regression model;
##'   \code{\link[mets]{resmeanIPCW}} for IPCW mean estimation;
##'   \code{\link[lava]{estimate}} for delta-method inference.
##'
##' @examples
##' \dontrun{
##' library(mets)
##' data(bmt)
##' bmt$id <- seq_len(nrow(bmt))
##'
##' ## --- competing risks CIF/Survival (default) ---
##' fit <- brier_binreg(
##'   Event(time, cause) ~ tcell + platelet + age + cluster(id),
##'   data = bmt, time = 50,
##'   rhs = list(small = ~age, full = ~tcell + platelet + age)
##' )
##' summary(fit)                        ## log scale
##' summary(fit, transform = exp)       ## Brier scale
##'
##' ## multiple time points
##' fitm <- brier_binreg(
##'   Event(time, cause) ~ tcell + platelet + age + cluster(id),
##'   data = bmt, time = c(50, 60),
##'   rhs = list(small = ~age, full = ~tcell + platelet + age)
##' )
##' summary(fitm, transform = exp)
##'
##' ## joint inference across models via lava::estimate
##' e <- estimate(fitm)
##' e
##' estimate(e, f = exp)               ## back-transform
##'
##' ## stratified IPCW weights for Brier score
##' fit_strat <- brier_binreg(
##'   Event(time, cause) ~ tcell + platelet + age + cluster(id),
##'   data = bmt, time = 50,
##'   rhs = list(small = ~age, full = ~tcell + platelet + age),
##'   brier.cens.model = ~strata(tcell)
##' )
##'
##' ## --- RMST --- survival case only
##' fit_rmst <- brier_binreg(
##'   Event(time, cause) ~ tcell + platelet + age + cluster(id),
##'   data = bmt, time = 50, cause = 1:2,
##'   rhs    = list(small = ~age, full = ~tcell + platelet + age),
##'   outcome = "rmst"
##' )
##' summary(fit_rmst, transform = exp)
##'
##' ## --- RMTL --- competing risks or survival (RMTL)
##' fit_rmtl <- brier_binreg(
##'   Event(time, cause) ~ tcell + platelet + age + cluster(id),
##'   data = bmt, time = 50,
##'   rhs    = list(small = ~age, full = ~tcell + platelet + age),
##'   outcome = "rmtl"
##' )
##' summary(fit_rmtl, transform = exp)
##'
##' ## --- no IPCW adjustment, regression ---------
##' fit_reg <- brier_binreg(
##'   time ~ tcell + platelet + age + cluster(id),
##'   data = bmt, 
##'   rhs    = list(small = ~age, full = ~tcell + platelet + age)
##' )
##' summary(fit_reg, transform = exp)
##'
##' ## --- no IPCW adjustment  binomial-regression  ---------------
##' bmt$binY  <- as.numeric(bmt$time>30)
##' fit_bin <- brier_binreg(
##'   binY ~ tcell + platelet + age + cluster(id),
##'   data = bmt, model="logit",
##'   rhs    = list(small = ~age, full = ~tcell + platelet + age)
##' )
##' summary(fit_bin, transform = exp)
##'
##'
##' }
##' @export
brier_binreg <- function(formula, data, time=NULL,
                                     fold             = 5,
                                     rhs              = NULL,
                                     rhs0             = ~+1,
                                     cens.model       = ~+1,
                                     brier.cens.model = ~+1,
                                     brier.args       = list(),
                                     cv.split.index   = NULL,
                                     ...) { ## {{{
  ## ------------------------------------------------------------------
  ## Small helpers (identical to brier_binregStrata)
  ## ------------------------------------------------------------------
  extract_cluster_id <- function(f) {
    tt     <- terms(f, specials = "cluster")
    cl_idx <- attr(tt, "specials")$cluster
    if (is.null(cl_idx)) return(NULL)
    cl_term <- attr(tt, "variables")[[cl_idx + 1L]]
    as.character(cl_term[[2]])
  }
  drop_cluster <- function(f) {
    tt     <- terms(f, specials = "cluster")
    cl_idx <- attr(tt, "specials")$cluster
    if (is.null(cl_idx)) return(f)
    cl_var <- attr(tt, "variables")[[cl_idx + 1L]]
    update(f, paste("~ . -", deparse(cl_var)))
  }

  ## ------------------------------------------------------------------
  ## Parse top-level formula
  ## ------------------------------------------------------------------
  lhs_str  <- deparse(formula[[2]])
  lhs_call <- formula[[2]]
  if (length(lhs_call) == 4L)
    stop("brier_binregStrata_loop does not support delayed entry (left-truncation).\n",
         "  Found: ", deparse(lhs_call))

  if ((length(lhs_call)>1) & is.null(time)) stop("must give time for binreg object with Event object responses\n"); 
  ## time not needed for continuous responses 
  noIPCW <- 0
  if (length(lhs_call)==1 &  is.null(time)) { time  <- 0; noIPCW <- 1}
###  time_var   <- as.character(lhs_call[[2]])
###  status_var <- as.character(lhs_call[[3]])

  id         <- extract_cluster_id(formula)
  if (is.null(id)) {
    id         <- "id__"
    data[[id]] <- seq_len(nrow(data))
    formula    <- update(formula, paste("~ . + cluster(id__)"))
  }
  base_resmean <- as.formula(paste(lhs_str, "~ +1 + cluster(", id, ")"))

  ## ------------------------------------------------------------------
  ## Parse rhs list
  ## ------------------------------------------------------------------
  formula_rhs_clean <- drop_cluster(
    reformulate(attr(terms(formula), "term.labels"), intercept = FALSE))

  if (is.null(rhs))        rhs <- list(model1 = formula_rhs_clean)
  if (!is.list(rhs))       rhs <- list(model1 = rhs)
  if (is.null(names(rhs)) || any(names(rhs) == ""))
    names(rhs) <- paste0("model", seq_along(rhs))
  rhs  <- lapply(rhs,  drop_cluster)
  rhs0 <- drop_cluster(rhs0)

  build_full <- function(f) {
    labs    <- attr(terms(f), "term.labels")
    rhs_str <- if (length(labs) == 0) "1" else paste(labs, collapse = " + ")
    as.formula(paste(lhs_str, "~", rhs_str, "+ cluster(", id, ")"))
  }
  if (!id %in% names(data)) stop("id variable '", id, "' not found in data")

  ## ------------------------------------------------------------------
  ## CV fold index -- computed once, reused for all time points
  ## ------------------------------------------------------------------
  as_fold_list <- function(x) {
    if (is.list(x))                      return(x)
    if (is.numeric(x) || is.integer(x)) return(split(seq_along(x), x))
    stop("Unsupported cv.split.index format")
  }
  n  <- nrow(data)
  dd <- if (is.null(cv.split.index)) folds(n, fold) else as_fold_list(cv.split.index)
  data$fold <- 0L
  for (k in seq_len(fold)) data$fold[dd[[k]]] <- k
  strata_vec <- data$fold - 1L
  build_strata_formula <- function(f) {
    labs    <- attr(terms(f), "term.labels")
    rhs_str <- if (length(labs) == 0) "1" else paste(labs, collapse = " + ")
    as.formula(paste(lhs_str, "~", rhs_str,
                     "+ strata(fold) + cluster(", id, ")"))
  }

  ## ------------------------------------------------------------------
  ## cv_iid_correctS -- block-sparse version (identical to brier_binregStrata)
  ## ------------------------------------------------------------------

cv_iid_correctS <- function(brier_fit, outff_fit, newdata,
                               pred_vec, ipcw_wt_vec, outcome_vec,
                               binreg_model) { ## {{{ 
    nstrata <- outff_fit$nstrata
    p       <- length(outff_fit$coef) / nstrata
    n       <- length(pred_vec)
    des        <- update_design(outff_fit$design, data = newdata,
                                       response = FALSE)
    X          <- des$x
    strata_idx <- if (!is.null(des$strata)) {
      as.integer(des$strata) - 1L
    } else if ("fold" %in% names(newdata)) {
      as.integer(newdata$fold) - 1L
    } else {
      stop("Cannot determine stratum for newdata: needs 'fold' column.")
    }
###    Xs <- matrix(0.0, n, nstrata * p)
###    for (s in seq_len(nstrata)) {
###      idx  <- which(strata_idx == (s - 1L))
###      cols <- (s - 1L) * p + seq_len(p)
###      Xs[idx, cols] <- X[idx, , drop = FALSE]
###    }
    mu_prime <- switch(binreg_model,
      logit = pred_vec * (1 - pred_vec),
      exp   = pred_vec,
      lin   = rep(1.0, n),
      stop("unrecognised binreg model '", binreg_model, "'")
    )
    resid       <- (outcome_vec - pred_vec) * ipcw_wt_vec

    Dbrier <- -2* Msumstrata(X * mu_prime * resid,strata_idx,nstrata)/n
    Dbrier_full  <-  c(t(Dbrier))

###   Dbrier_full <- -2 * colSums(Xs * mu_prime * resid) / n
###   print(Dbrier_full)

    iid(brier_fit) + iid(outff_fit) %*% Dbrier_full
  } ## }}} 

  eco <- function(bbrier) { ## {{{ 
    bbrier$design <- NULL
    bbrier$cens.weights <- bbrier$cens.strata <- bbrier$MGciid <- bbrier$call.id
    bbrier$name.id <- bbrier$iid.naive  <- bbrier$id <- NULL
    bbrier$design$y <- NULL
    bbrier$design$data <- NULL
    bbrier$design$x <- NULL
    bbrier$design$strata <- NULL
    return(bbrier)
  } ## }}} 

  ## ------------------------------------------------------------------
  ## resmeanIPCW wrapper / log-scale estimate helper
  ## ------------------------------------------------------------------
  n  <- nrow(data)
  do_resmean <- function(dat, Ydirect, tt) {
    do.call(binregStrata,
            c(list(base_resmean, dat, time = tt,
                   Ydirect    = Ydirect,
		   outcome="rmst",
                   cens.model = brier.cens.model,
                   model = "lin"),
              brier.args))
  }

  do_mean <- function(dat, brier_vec,tt) {
    dat$brier__ <- brier_vec
    do.call(lm, c(list(brier__ ~ 1, data = dat)))
  }

  log_est <- function(fit, iid_corrected = NULL) {
    if (!is.null(iid_corrected)) {
      nn <- nrow(iid_corrected)
      e <- lava::estimate(NULL,
                          coef = fit$coef,
                          IC   = iid_corrected * nn,
                          f    = log)
    } else {
      e <- estimate(fit, f = log)
    }
    list(coef = as.numeric(e$coef),
         se   = as.numeric(diag(vcov(e))^.5),
         iid  = iid(e))
  }

  ## ------------------------------------------------------------------
  ## Explicit nested for-loops: outer over time, inner over rhs models
  ## ------------------------------------------------------------------
  cl       <- match.call()
  ntime    <- length(time)
  nmodels  <- length(rhs)
  results  <- vector("list", ntime)
  names(results) <- as.character(time)

  for (ti in seq_len(ntime)) {
    tt <- time[ti]

    ## ---------------- null model: apparent ----------------
    outb         <- binregStrata(build_full(rhs0), data, time = tt,...)
    binreg_model <- outb$model
    outcome      <- as.numeric(outb$Y)
    ipcw_wt_b    <- as.numeric(outb$IPCW)
    predb        <- predict(outb, data, se = FALSE)
    data$brier_  <- (outcome - predb)^2
    if (noIPCW==0) bbrier       <- do_resmean(data, data$brier_, tt) else bbrier       <- do_mean(data, data$brier_, tt)
    bbrier <- eco(bbrier)
    e <- estimate(bbrier, f = log)
    bbrier_log   <- log_est(bbrier)

    ## ---------------- null model: CV ----------------
    f0_strata  <- build_strata_formula(rhs0)

    outfb      <- binregStrata(f0_strata, data = data, time = tt,
                                strata     = strata_vec,
                                cens.model = cens.model,
                                cv.fold    = TRUE, low.memory=TRUE,...)
    predbcv     <- predict(outfb, data,se=0)
    data$brier_ <- (outcome - predbcv)^2
    if (noIPCW==0) 
    bbrierCV    <- do_resmean(data, data$brier_, tt)
    else bbrierCV    <- do_mean(data, data$brier_, tt)
    bbrierCV    <- eco(bbrierCV)
    iid.bbrierCV <- cv_iid_correctS(bbrierCV, outfb, data,
                                     predbcv, ipcw_wt_b, outcome,
                                     binreg_model)
    outfb$iid <- NULL
    bbrierCV_log <- log_est(bbrierCV, iid_corrected = iid.bbrierCV)

    null_res <- list(
      coef = c(bbrier   = bbrier_log$coef,
               bbrierCV = bbrierCV_log$coef),
      se   = c(bbrier   = bbrier_log$se,
               bbrierCV = bbrierCV_log$se),
      iid  = list(bbrierCV = bbrierCV_log$iid)
    )

    ## ---------------- candidate models: explicit for-loop ----------------
    model_res        <- vector("list", nmodels)
    names(model_res) <- names(rhs)

    for (mi in seq_len(nmodels)) {
      f <- rhs[[mi]]

      ## apparent
      out         <- binregStrata(build_full(f), data, time = tt, ...)
      outcome_m   <- as.numeric(out$Y)
      ipcw_wt_m   <- as.numeric(out$IPCW)
      pred        <- predict(out, data, se = FALSE)
      data$brier_ <- (outcome_m - pred)^2
      out         <- NULL
      if (noIPCW==0) brier       <- do_resmean(data, data$brier_, tt) else brier   <- do_mean(data, data$brier_, tt)
      brier       <- eco(brier)
      brier_log   <- log_est(brier)
      dbrier      <- brier_log$coef - bbrier_log$coef
      se.dbrier   <- as.numeric(crossprod(brier_log$iid - bbrier_log$iid)^.5)

      ## CV
      f_strata <- build_strata_formula(f)
      outff    <- binregStrata(f_strata, data = data, time = tt,
                                strata     = strata_vec,
                                cens.model = cens.model,
                                cv.fold    = TRUE, low.memory=TRUE,...)
      pred        <- predict(outff, data,se=0)
      data$brier_ <- (outcome_m - pred)^2
      if (noIPCW==0) 
      brierCV     <- do_resmean(data, data$brier_, tt)
      else brierCV     <- do_mean(data, data$brier_, tt)
      brierCV     <- eco(brierCV)
      iid.brierCV <- cv_iid_correctS(brierCV, outff, data,
                                      pred, ipcw_wt_m, outcome_m, binreg_model)
      outff        <- NULL
      brierCV_log <- log_est(brierCV, iid_corrected = iid.brierCV)
      dbrierCV    <- brierCV_log$coef - bbrierCV_log$coef
      se.dbrierCV <- as.numeric(crossprod(brierCV_log$iid - bbrierCV_log$iid)^.5)

      model_res[[mi]] <- list(
        coef  = c(brier     = brier_log$coef,
                  brierCV   = brierCV_log$coef),
        se    = c(brier     = brier_log$se,
                  brierCV   = brierCV_log$se),
        dcoef = c(dbrier    = dbrier,
                  dbrierCV  = dbrierCV),
        sed   = c(se.dbrier   = se.dbrier,
                  se.dbrierCV = se.dbrierCV),
        iid   = list(brierCV = brierCV_log$iid)
      )
    } ## end inner for over models

    ## iid matrix: n x (1 + nmodels)
    iid_cols <- vector("list", nmodels + 1L)
    iid_cols[[1]] <- bbrierCV_log$iid
    bbrierCV_log$iid <- NULL
    for (mi in seq_len(nmodels)) {
	    iid_cols[[mi + 1L]] <- model_res[[mi]]$iid$brierCV
	    model_res[[mi]]$iid$brierCV <- NULL
    }
    iid_mat <- do.call(cbind, iid_cols)
    colnames(iid_mat) <- c("null", names(rhs))

    results[[ti]] <- structure(
      list(null             = null_res,
           models           = model_res,
           iid_mat          = iid_mat,
           cv.split.index   = dd,
           time             = tt,
           fold             = fold,
           rhs              = rhs,
           rhs0             = rhs0,
           model            = binreg_model,
           cens.model       = cens.model,
           brier.cens.model = brier.cens.model),
      class = "binregCV"
    )
  } ## end outer for over time

  ## ------------------------------------------------------------------
  ## Dispatch: single or multiple time points
  ## ------------------------------------------------------------------
  if (ntime == 1L) {
    out      <- results[[1]]
    out$call <- cl
    out
  } else {
    res        <- results
    res$call   <- cl
    res$rhs    <- rhs
    res$rhs0   <- rhs0
    class(res) <- "binregCV_list"
    res
  }
} ## }}}


## {{{  infrastructure for brier_binreg
## =============================================================================
## estimate.binregCV  -- lava-style joint estimate  (null + all CV models)
## =============================================================================



##' Joint estimate across CV models for a binregCV object
##'
##' @param x A `binregCV` object.
##' @param f Optional transformation function passed to `lava::estimate`.
##' @param ... Not currently used.
##' @export
estimate.binregCV <- function(x , f = NULL, ...) {
  coefs <- c(
    null = x$null$coef["bbrierCV"],
    setNames(
      sapply(x$models, function(m) m$coef["brierCV"]),
      names(x$models)
    )
  )
  iid_mat <- x$iid_mat   # n x (1 + nmodels)
  n       <- nrow(iid_mat)
  e       <- lava::estimate(NULL, coef = coefs, IC = iid_mat * n)
  if (!is.null(f)) e <- estimate(e, f = f)
  e
}


## =============================================================================
## summarize_brier_binregCV  (works for binregCV and binregCV_list)
## =============================================================================
summarize_brier_binregCV <- function(fit,
                                      transform  = NULL,
                                      conf.level = 0.95,
                                      type       = c("cv", "both")) {
  type <- match.arg(type)

  fit_call <- fit$call
  fit_rhs  <- fit$rhs
  fit_rhs0 <- fit$rhs0

  if (inherits(fit, "binregCV")) {
    tname    <- as.character(fit$time)
    fit_call <- fit$call
    fit_rhs  <- fit$rhs
    fit_rhs0 <- fit$rhs0
    fit      <- setNames(list(fit), tname)
  }

  times <- names(fit)[!(names(fit) %in% c("call", "rhs", "rhs0"))]
  if (length(times) == 0L) stop("No time points found in fit")

  pdiff <- lava::pairwise.diff

  col_se <- function(cm) grep("^Std", colnames(cm), value = TRUE)[1]
  col_lo <- function(cm) grep("%",    colnames(cm), value = TRUE)[1]
  col_hi <- function(cm) grep("%",    colnames(cm), value = TRUE)[2]
  col_p  <- function(cm) grep("^P",   colnames(cm), value = TRUE)[1]

  parse_lava_contrast <- function(rn) {
    m    <- gregexpr("(?<=\\[)[^\\]]+(?=\\])", rn, perl = TRUE)
    hits <- regmatches(rn, m)[[1]]
    c(reference = hits[1], model = hits[2])
  }
  strip_t <- function(x) sub("^t[^.]+\\.", "", x)

  coefmat_df <- function(e, time_val) {
    cm <- e$coefmat
    data.frame(
      time  = as.numeric(time_val),
      model = names(coef(e)),
      Brier = as.numeric(cm[, "Estimate"]),
      SE    = as.numeric(cm[, col_se(cm)]),
      lower = as.numeric(cm[, col_lo(cm)]),
      upper = as.numeric(cm[, col_hi(cm)]),
      stringsAsFactors = FALSE
    )
  }

  apparent_coefmat_df <- function(obj, time_val, tfun, cl) {
    z <- qnorm(1 - (1 - cl) / 2)
    coefs <- c(
      null = obj$null$coef["bbrier"],
      setNames(
        sapply(names(obj$models), function(nm) obj$models[[nm]]$coef["brier"]),
        names(obj$models)
      )
    )
    ses <- c(
      null = obj$null$se["bbrier"],
      setNames(
        sapply(names(obj$models), function(nm) obj$models[[nm]]$se["brier"]),
        names(obj$models)
      )
    )
    lower_log <- coefs - z * ses
    upper_log <- coefs + z * ses
    if (!is.null(tfun)) {
      eps    <- sqrt(.Machine$double.eps)
      fprime <- (tfun(coefs + eps) - tfun(coefs - eps)) / (2 * eps)
      data.frame(time  = as.numeric(time_val),
                 model = names(coefs),
                 Brier = as.numeric(tfun(coefs)),
                 SE    = as.numeric(abs(fprime) * ses),
                 lower = as.numeric(tfun(lower_log)),
                 upper = as.numeric(tfun(upper_log)),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(time  = as.numeric(time_val),
                 model = names(coefs),
                 Brier = as.numeric(coefs),
                 SE    = as.numeric(ses),
                 lower = as.numeric(lower_log),
                 upper = as.numeric(upper_log),
                 stringsAsFactors = FALSE)
    }
  }

  models_list          <- vector("list", length(times))
  comps_list           <- vector("list", length(times))
  apparent_models_list <- vector("list", length(times))
  names(models_list) <- names(comps_list) <-
    names(apparent_models_list) <- times

  for (tt in times) {
    obj <- fit[[tt]]

    e <- if (is.null(transform)) estimate(obj)
         else                    estimate(obj, f = transform)
    nmodels <- length(coef(e))
    models_list[[tt]] <- coefmat_df(e, tt)

    if (nmodels > 1L) {
      ec  <- estimate(e, -pdiff(nmodels))
      cmc <- ec$coefmat
      rn  <- rownames(cmc)
      parsed <- t(sapply(rn, parse_lava_contrast))
      comps_list[[tt]] <- data.frame(
        time        = as.numeric(tt),
        model       = strip_t(parsed[, "model"]),
        reference   = strip_t(parsed[, "reference"]),
        delta.Brier = as.numeric(cmc[, "Estimate"]),
        SE          = as.numeric(cmc[, col_se(cmc)]),
        lower       = as.numeric(cmc[, col_lo(cmc)]),
        upper       = as.numeric(cmc[, col_hi(cmc)]),
        P           = as.numeric(cmc[, col_p(cmc)]),
        stringsAsFactors = FALSE
      )
    }

    if (type == "both")
      apparent_models_list[[tt]] <- apparent_coefmat_df(obj, tt, transform, conf.level)
  }

  obj1       <- fit[[times[1]]]
  fit_fold   <- obj1$fold
  if (is.null(fit_rhs))  fit_rhs  <- obj1$rhs
  if (is.null(fit_rhs0)) fit_rhs0 <- obj1$rhs0
  model_list <- c(
    list(`null (rhs0)` = deparse(fit_rhs0)),
    setNames(lapply(fit_rhs, function(f) deparse(f[[length(f)]])),
             names(fit_rhs))
  )

  models_df <- do.call(rbind, models_list)
  comps_df  <- do.call(rbind, Filter(Negate(is.null), comps_list))
  if (is.null(comps_df)) comps_df <- data.frame()
  rownames(models_df) <- NULL
  if (nrow(comps_df) > 0) rownames(comps_df) <- NULL

  apparent_models_df <- if (type == "both")
    do.call(rbind, apparent_models_list) else NULL

  structure(
    list(models               = models_df,
         comparisons          = comps_df,
         apparent.models      = apparent_models_df,
         apparent.comparisons = NULL,
         type                 = type,
         transform            = transform,
         conf.level           = conf.level,
         times                = as.numeric(times),
         call                 = fit_call,
         fold                 = fit_fold,
         model_list           = model_list),
    class = "summary_binregCV_multi"
  )
}


## =============================================================================
## summary methods
## =============================================================================

##' Summary method for binregCV objects
##'
##' @param object A \code{binregCV} object, as returned by
##'   \code{\link{brier_binreg}} for a single evaluation time.
##' @param transform Optional back-transformation function applied to Brier
##'   scores and confidence limits, e.g. \code{exp} to undo the internal log
##'   scale.
##' @param conf.level Confidence level for the reported intervals (default
##'   0.95).
##' @param type Either \code{"cv"} (default) for cross-validated Brier scores
##'   only, or \code{"both"} to also include apparent (in-sample) scores.
##' @param ... Not currently used.
##' @aliases print.binregCV
##' @export
summary.binregCV <- function(object, transform = NULL,
                              conf.level = 0.95,
                              type = c("cv", "both"), ...) {
  summarize_brier_binregCV(object,
                            transform  = transform,
                            conf.level = conf.level,
                            type       = match.arg(type))
}

##' @export
print.binregCV <- function(x, transform = NULL, 
                              conf.level = 0.95, type = c("cv", "both"), ...) { ## {{{ 
  print(summarize_brier_binregCV(x,
                            transform  = transform,
                            conf.level = conf.level,
                            type       = match.arg(type)))
} ## }}} 


##' Summary method for binregCV_list objects
##'
##' @param object A \code{binregCV_list} object, as returned by
##'   \code{\link{brier_binreg}} when \code{time} has length greater than
##'   one.
##' @param transform Optional back-transformation function applied to Brier
##'   scores and confidence limits, e.g. \code{exp} to undo the internal log
##'   scale.
##' @param conf.level Confidence level for the reported intervals (default
##'   0.95).
##' @param type Either \code{"cv"} (default) for cross-validated Brier scores
##'   only, or \code{"both"} to also include apparent (in-sample) scores.
##' @param ... Not currently used.
##' @export
summary.binregCV_list <- function(object, transform = NULL,
                                   conf.level = 0.95,
                                   type = c("cv", "both"), ...) {
  summarize_brier_binregCV(object,
                            transform  = transform,
                            conf.level = conf.level,
                            type       = match.arg(type))
}


## =============================================================================
## print.summary_binregCV_multi
## =============================================================================

##' Print method for summary_binregCV_multi objects
##'
##' @param x A `summary_binregCV_multi` object.
##' @param digits Number of digits to round to.
##' @param ... Not currently used.
##' @export
print.summary_binregCV_multi <- function(x, digits = 4, ...) {
  sep_line <- paste(rep("-", 56), collapse = "")

  rd <- function(df) {
    num_cols <- sapply(df, is.numeric)
    df[num_cols] <- lapply(df[num_cols], round, digits)
    df
  }

  scale_note <- if (is.null(x$transform)) "[log scale]" else "[transformed scale]"
  cat(sprintf("\nCross-validated Brier score summary  %s\n", scale_note))
  if (!is.null(x$call)) { cat("Call: "); print(x$call) }
  cat(sprintf("Evaluation time(s): %s   Folds: %d   CI level: %g%%\n",
              paste(x$times, collapse = ", "), x$fold, 100 * x$conf.level))
  cat("Estimator: binregStrata  |  IID: block-diagonal delta-method correction\n")

  cat("\nModels:\n")
  for (nm in names(x$model_list)) {
    rhs_str <- trimws(gsub("- 1", "", gsub("\\s+", " ", x$model_list[[nm]])))
    cat(sprintf("  %-18s : %s\n", nm, rhs_str))
  }
  cat("\n")

  if (!is.null(x$type) && x$type == "both" && !is.null(x$apparent.models)) {
    cat("Apparent (no CV) Brier scores:\n\n")
    disp_app <- rd(x$apparent.models)[, c("model", "time", "Brier", "SE",
                                           "lower", "upper")]
    for (tt in unique(disp_app$time)) {
      block      <- disp_app[disp_app$time == tt, , drop = FALSE]
      block$time <- NULL
      cat(sprintf("  time = %g\n", tt))
      print(block, row.names = FALSE, quote = FALSE)
      cat(sep_line, "\n")
    }
    cat("Note: apparent comparisons not shown (IID not stored for apparent fits).\n\n")
  }

  cat("CV Brier scores by model:\n\n")
  disp_mod <- rd(x$models)[, c("model", "time", "Brier", "SE", "lower", "upper")]
  for (tt in unique(disp_mod$time)) {
    block      <- disp_mod[disp_mod$time == tt, , drop = FALSE]
    block$time <- NULL
    cat(sprintf("  time = %g\n", tt))
    print(block, row.names = FALSE, quote = FALSE)
    cat(sep_line, "\n")
  }

  if (!is.null(x$comparisons) && nrow(x$comparisons) > 0L) {
    cat("\nModel comparisons (delta CV Brier):\n\n")
    comps <- rd(x$comparisons)
    for (tt in unique(comps$time)) {
      block <- comps[comps$time == tt, , drop = FALSE]
      cat(sprintf("  time = %g\n", tt))
      print(block, row.names = FALSE, quote = FALSE)
      cat(sep_line, "\n")
    }
  }

  if (is.null(x$transform))
    cat("\nNote: use summary(fit, transform = exp) for Brier scores on the original scale.\n")

  invisible(x)
}

## }}} 

## =============================================================================
## estimate.binregCV_list -- joint estimate across ALL time points and models
## =============================================================================
##' @export
estimate.binregCV_list <- function(x, f = NULL, ...) { ## {{{ 
  times <- names(x)[!(names(x) %in% c("call", "rhs", "rhs0"))]
  if (length(times) == 0L) stop("No time points found in object")

  ## Collect per-time coefficient vectors and iid matrices, with
  ## time-qualified names so e.g. "null" at t=50 and t=60 don't collide.
  coef_list <- vector("list", length(times))
  iid_list  <- vector("list", length(times))
  names(coef_list) <- names(iid_list) <- times

  n_ref <- NULL
  for (tt in times) {
    obj <- x[[tt]]
    coefs <- c(
      null = obj$null$coef["bbrierCV"],
      setNames(
        sapply(obj$models, function(m) m$coef["brierCV"]),
        names(obj$models)
      )
    )
    nm <- paste0("t", tt, ".", names(coefs))
    names(coefs) <- nm

    im <- obj$iid_mat                 ## n x (1+nmodels), cols: null, model1, ...
    colnames(im) <- nm

    if (is.null(n_ref)) {
      n_ref <- nrow(im)
    } else if (nrow(im) != n_ref) {
      stop("iid_mat row counts differ across time points (", n_ref,
           " vs ", nrow(im), " at time ", tt, "); cannot combine into a ",
           "single joint estimate. This typically means rows correspond ",
           "to different clusters/ids across time points.")
    }

    coef_list[[tt]] <- coefs
    iid_list[[tt]]  <- im
  }

  coefs_all   <- do.call(c, coef_list)
  iid_mat_all <- do.call(cbind, iid_list)   ## n x (ntimes * (1+nmodels))

  n <- nrow(iid_mat_all)
  e <- lava::estimate(NULL, coef = coefs_all, IC = iid_mat_all * n)
  if (!is.null(f)) e <- estimate(e, f = f)
  e
} ## }}} 


brier_lm <- function(formula, data,
                      fold             = 5,
                      rhs              = NULL,
                      rhs0             = ~+1,
                      brier.args       = list(),
                      cv.split.index   = NULL,
                      ...) { ## {{{
  ## ------------------------------------------------------------------
  ## Small helpers (same as brier_binreg)
  ## ------------------------------------------------------------------
  extract_cluster_id <- function(f) {
    tt     <- terms(f, specials = "cluster")
    cl_idx <- attr(tt, "specials")$cluster
    if (is.null(cl_idx)) return(NULL)
    cl_term <- attr(tt, "variables")[[cl_idx + 1L]]
    as.character(cl_term[[2]])
  }
  drop_cluster <- function(f) {
    tt     <- terms(f, specials = "cluster")
    cl_idx <- attr(tt, "specials")$cluster
    if (is.null(cl_idx)) return(f)
    cl_var <- attr(tt, "variables")[[cl_idx + 1L]]
    update(f, paste("~ . -", deparse(cl_var)))
  }
  ## ------------------------------------------------------------------
  ## Parse top-level formula: plain Y ~ ..., no Surv/Event lhs
  ## ------------------------------------------------------------------
  lhs_str <- deparse(formula[[2]])
  id      <- extract_cluster_id(formula)
  if (is.null(id)) {
    id         <- "id__"
    data[[id]] <- seq_len(nrow(data))
    formula    <- update(formula, paste("~ . + cluster(", id, ")"))
  }
  if (!id %in% names(data)) stop("id variable '", id, "' not found in data")
  ## ------------------------------------------------------------------
  ## Parse rhs list
  ## ------------------------------------------------------------------
  formula_rhs_clean <- drop_cluster(
    reformulate(attr(terms(formula), "term.labels"), intercept = FALSE))
  if (is.null(rhs))        rhs <- list(model1 = formula_rhs_clean)
  if (!is.list(rhs))       rhs <- list(model1 = rhs)
  if (is.null(names(rhs)) || any(names(rhs) == ""))
    names(rhs) <- paste0("model", seq_along(rhs))
  rhs  <- lapply(rhs,  drop_cluster)
  rhs0 <- drop_cluster(rhs0)
###  build_lm_formula <- function(f) {
###    ## plain lm formula, no cluster/strata specials
###    labs    <- attr(terms(f), "term.labels")
###    rhs_str <- if (length(labs) == 0) "1" else paste(labs, collapse = " + ")
###    as.formula(paste(lhs_str, "~", rhs_str))
###  }
###  build_strata_formula <- function(f) {
###    ## for lmStrata: needs strata(fold) + cluster(id)
###    labs    <- attr(terms(f), "term.labels")
###    rhs_str <- if (length(labs) == 0) "1" else paste(labs, collapse = " + ")
###    as.formula(paste(lhs_str, "~", rhs_str,
###                     "+ strata(fold) + cluster(", id, ")"))
###  }

  build_lm_formula <- function(f) {
  ## plain lm formula, no cluster/strata specials
  tt        <- terms(f)
  labs      <- attr(tt, "term.labels")
  has_int   <- attr(tt, "intercept") == 1

  rhs_str <- if (length(labs) == 0) "1" else paste(labs, collapse = " + ")
  if (!has_int) rhs_str <- paste(rhs_str, "- 1")

  as.formula(paste(lhs_str, "~", rhs_str))
}

build_strata_formula <- function(f) {
  ## for lmStrata: needs strata(fold) + cluster(id)
  tt        <- terms(f)
  labs      <- attr(tt, "term.labels")
  has_int   <- attr(tt, "intercept") == 1

  rhs_str <- if (length(labs) == 0) "1" else paste(labs, collapse = " + ")
  if (!has_int) rhs_str <- paste(rhs_str, "- 1")

  as.formula(paste(lhs_str, "~", rhs_str,
                   "+ strata(fold) + cluster(", id, ")"))
}


  ## ------------------------------------------------------------------
  ## CV fold index -- computed once
  ## ------------------------------------------------------------------
  as_fold_list <- function(x) {
    if (is.list(x))                    return(x)
    if (is.numeric(x) || is.integer(x)) return(split(seq_along(x), x))
    stop("Unsupported cv.split.index format")
  }
  n  <- nrow(data)
  dd <- if (is.null(cv.split.index)) folds(n, fold) else as_fold_list(cv.split.index)
  data$fold <- 0L
  for (k in seq_len(fold)) data$fold[dd[[k]]] <- k
  strata_vec <- data$fold - 1L
  ## ------------------------------------------------------------------
  ## cv_iid_correctS -- lin-only: mu_prime == 1, no IPCW weighting
  ## ------------------------------------------------------------------
  cv_iid_correctS <- function(brier_fit, outff_fit, newdata,
                               pred_vec, outcome_vec,
                               binreg_model = "lin") { ## {{{
    nstrata <- outff_fit$nstrata
    p       <- length(outff_fit$coef) / nstrata
    n       <- length(pred_vec)
    des     <- update_design(outff_fit$design, data = newdata, response = FALSE)
    X       <- des$x
    strata_idx <- if (!is.null(des$strata)) {
      as.integer(des$strata) - 1L
    } else if ("fold" %in% names(newdata)) {
      as.integer(newdata$fold) - 1L
    } else {
      stop("Cannot determine stratum for newdata: needs 'fold' column.")
    }
    mu_prime <- switch(binreg_model,
      logit = pred_vec * (1 - pred_vec),
      exp   = pred_vec,
      lin   = rep(1.0, n),
      stop("unrecognised model '", binreg_model, "'")
    )
    resid  <- c(outcome_vec - pred_vec)          # no IPCW weight: no censoring
    Dbrier <- -2 * Msumstrata(X * mu_prime * resid, strata_idx, nstrata) / n
    Dbrier_full <- c(t(Dbrier))
    iid(brier_fit) + iid(outff_fit) %*% Dbrier_full
  } ## }}}
  eco <- function(bbrier) {
    bbrier$design        <- NULL
    bbrier$cens.weights   <- bbrier$cens.strata <- bbrier$MGciid <- NULL
    bbrier$name.id        <- bbrier$iid.naive   <- bbrier$id     <- NULL
    return(bbrier)
  }
  ## ------------------------------------------------------------------
  ## do_resmean: plain intercept-only lm on the brier residuals
  ## ------------------------------------------------------------------
  do_resmean <- function(dat, brier_vec) {
    dat$brier__ <- brier_vec
    do.call(lm, c(list(brier__ ~ 1, data = dat), brier.args))
  }

  log_est <- function(fit, iid_corrected = NULL) {
    if (!is.null(iid_corrected)) {
      nn <- nrow(iid_corrected)
      e  <- lava::estimate(NULL,
                           coef = coef(fit),
                           IC   = iid_corrected * nn,
                           f    = log)
    } else {
      e <- lava::estimate(fit, f = log)   ## uncorrected IC straight from lm
    }
    list(coef = as.numeric(e$coef),
         se   = as.numeric(diag(vcov(e))^.5),
         iid  = iid(e))
  }

  ## ------------------------------------------------------------------
  ## Single fit (no time index): null model
  ## ------------------------------------------------------------------
  cl <- match.call()
  y  <- data[[lhs_str]]
  ## ---------------- null model: apparent ----------------
  null_f <- print(build_lm_formula(rhs0))
  outb     <- lm(null_f, data = data)
  predb   <- as.numeric(predict(outb, data,se=0))
  data$brier_ <- (y - predb)^2
  bbrier      <- do_resmean(data, data$brier_)
  bbrier_log  <- log_est(bbrier)

  ## ---------------- null model: CV ----------------
  f0_strata <- build_strata_formula(rhs0)
  outfb     <- lmStrata(f0_strata, data = data,
                        strata     = strata_vec,
                     cv.fold = TRUE, low.memory = TRUE)
  predbcv     <- predict(outfb, data,se=0)
  data$brier_ <- (y - predbcv)^2
  bbrierCV    <- do_resmean(data, data$brier_)
  bbrierCV    <- eco(bbrierCV)
  iid.bbrierCV <- cv_iid_correctS(bbrierCV, outfb, data, predbcv, y, "lin")
  outfb$iid <- NULL
  bbrierCV_log <- log_est(bbrierCV, iid_corrected = iid.bbrierCV)
  null_res <- list(
    coef = c(bbrier = bbrier_log$coef, bbrierCV = bbrierCV_log$coef),
    se   = c(bbrier = bbrier_log$se,   bbrierCV = bbrierCV_log$se),
    iid  = list(bbrierCV = bbrierCV_log$iid)
  )

  ## ---------------- candidate models: explicit for-loop ----------------
  nmodels          <- length(rhs)
  model_res        <- vector("list", nmodels)
  names(model_res) <- names(rhs)
  for (mi in seq_len(nmodels)) {
    f <- rhs[[mi]]
    ## apparent
    f_strata <- build_strata_formula(f)
    out         <- lm(f_strata, data)
    pred        <- as.numeric(predict(out, data,se=0))
    data$brier_ <- (y - pred)^2
    out         <- NULL
    brier       <- do_resmean(data, data$brier_)
    brier_log   <- log_est(brier)
    dbrier      <- brier_log$coef - bbrier_log$coef
    se.dbrier   <- as.numeric(crossprod(brier_log$iid - bbrier_log$iid)^.5)

    ## CV
    f_strata <- build_strata_formula(f)
    outff    <- lmStrata(f_strata, data = data,
                            strata     = strata_vec,
                          cv.fold = TRUE, low.memory = TRUE)
    pred        <- predict(outff, data,se=0)
    data$brier_ <- (y - pred)^2
    brierCV     <- do_resmean(data, data$brier_)
    brierCV     <- eco(brierCV)
    iid.brierCV <- cv_iid_correctS(brierCV, outff, data, pred, y, "lin")
    outff       <- NULL
    brierCV_log <- log_est(brierCV, iid_corrected = iid.brierCV)
    dbrierCV    <- brierCV_log$coef - bbrierCV_log$coef
    se.dbrierCV <- as.numeric(crossprod(brierCV_log$iid - bbrierCV_log$iid)^.5)
    model_res[[mi]] <- list(
      coef  = c(brier    = brier_log$coef,   brierCV   = brierCV_log$coef),
      se    = c(brier    = brier_log$se,     brierCV   = brierCV_log$se),
      dcoef = c(dbrier   = dbrier,           dbrierCV  = dbrierCV),
      sed   = c(se.dbrier   = se.dbrier,     se.dbrierCV = se.dbrierCV),
      iid   = list(brierCV = brierCV_log$iid)
    )
  } ## end for over models
  ## iid matrix: n x (1 + nmodels)
  iid_cols      <- vector("list", nmodels + 1L)
  iid_cols[[1]] <- bbrierCV_log$iid
  bbrierCV_log$iid <- NULL
  for (mi in seq_len(nmodels)) {
    iid_cols[[mi + 1L]]         <- model_res[[mi]]$iid$brierCV
    model_res[[mi]]$iid$brierCV <- NULL
  }
  iid_mat <- do.call(cbind, iid_cols)
  colnames(iid_mat) <- c("null", names(rhs))

  out <- structure(
    list(null           = null_res,
         models         = model_res,
         iid_mat        = iid_mat,
         cv.split.index = dd,
         fold           = fold,
         rhs            = rhs,
         rhs0           = rhs0,
         model          = "lin",
         time           = 0,     ## numeric placeholder, not the string "NA"
         call           = cl),
    class = "binregCV"
  )
  out
} ## }}}

lmStrata <- function(formula, data, beta=NULL,
                      offset=NULL, weights=NULL, se=TRUE,
                      no.opt=FALSE, method="nr",
                      model=c("lin","logit","exp"),
                      strata=NULL, cv.fold=FALSE, low.memory=FALSE, ...)
{ ## {{{
  monotone <- TRUE
  cl <- match.call()
  ## ------------------------------------------------------------------
  ## Design
  ## ------------------------------------------------------------------
  des <- proc_design(
    formula, data = data,
    specials = c("offset","weights","cluster","strata"),
    intercept = TRUE
  )
  Y <- c(des$y)
  X <- des$x
  des.weights <- des$weights
  des.offset  <- des$offset
  id          <- des$cluster
  if (ncol(X) == 0) X <- matrix(nrow=0, ncol=0)
  call.id <- id
  conid   <- construct_id(id, nrow(X), namesX = rownames(X))
  name.id <- conid$name.id; id <- conid$id; nid <- conid$nid
  orig.id <- id

    ## ------------------------------------------------------------------
  ## Weights / offset  (explicit arg > formula-embedded special > default)
  ## ------------------------------------------------------------------
  if (!is.null(offset)) {
    offset <- offset
  } else if (!is.null(des.offset)) {
    offset <- des.offset
  } else {
    offset <- rep(0, nrow(X))
  }

  if (!is.null(weights)) {
    weights <- weights
  } else if (!is.null(des.weights)) {
    weights <- des.weights
  } else {
    weights <- rep(1, nrow(X))
  }

  ## ------------------------------------------------------------------
  ## Strata  (explicit arg > formula-embedded special > default)
  ## ------------------------------------------------------------------
  des.strata <- des$strata
  if (!is.null(strata)) {
    if (!is.numeric(strata)) stop("strata must be numeric\n")
    ustrata <- sort(unique(strata)); nstrata <- length(ustrata)
    strata  <- fast.approx(ustrata, strata) - 1
  } else if (!is.null(des.strata)) {
    strata  <- des.strata
    ustrata <- sort(unique(strata)); nstrata <- length(ustrata)
    strata  <- fast.approx(ustrata, strata) - 1
  } else {
    strata  <- rep(0, nrow(X)); nstrata <- 1
  }


  ## ------------------------------------------------------------------
  ## Response / design set-up
  ## ------------------------------------------------------------------
  expit <- function(z) 1 / (1 + exp(-z))
  p     <- ncol(X)
  if (is.null(beta)) beta <- rep(0, p)
  augmentation <- matrix(0, nstrata, p)
  X  <- as.matrix(X)
  X2 <- .Call("vecCPMat", X)$XX

  model <- model[1]
  ## ------------------------------------------------------------------
  ## obj(): full objective (used for final all=TRUE call and NR); this
  ## general path is only exercised for model %in% c("logit","exp") --
  ## "lin" is handled by the closed-form fast_lin_fit() below.
  ## ------------------------------------------------------------------
  obj <- function(pp, all = FALSE) { # {{{
    np <- length(pp[[1]])
    if (cv.fold) {
      ploglik_list  <- vector("list", nstrata)
      gradient_list <- vector("list", nstrata)
      D2log_list    <- vector("list", nstrata)
      for (s in 0:(nstrata - 1)) {
        idx_c <- which(strata != s)
        X_c   <- X[idx_c,  , drop=FALSE]
        X2_c  <- X2[idx_c, , drop=FALSE]
        Y_c   <- Y[idx_c]
        w_c   <- weights[idx_c]
        off_c <- offset[idx_c]
        aug_s <- augmentation[s + 1, ]
        lp_c  <- as.vector(X_c %*% pp[[s + 1]]) + off_c
        if (model == "exp") {
          pr_c  <- exp(lp_c)
          if (monotone) {
            Dl_c  <- w_c * X_c * c(Y_c - pr_c)
            D2l_c <- c(w_c * pr_c) * X2_c
          } else {
            Dl_c  <- w_c * X_c * c((Y_c - pr_c) * pr_c)
            D2l_c <- c(w_c * pr_c^2) * X2_c
          }
          pl_c <- 0
        } else if (model == "logit") {
          pr_c <- expit(lp_c)
          if (monotone) {
            Dl_c  <- w_c * X_c * c(Y_c - pr_c)
            D2l_c <- c(w_c * pr_c / (1 + exp(lp_c))) * X2_c
          } else {
            varp  <- pr_c / (1 + exp(lp_c))
            Dl_c  <- w_c * X_c * c((Y_c - pr_c) * varp)
            D2l_c <- c(w_c * varp^2) * X2_c
          }
          pl_c <- 0
        } else stop("obj() general path only supports logit/exp\n")
        ploglik_list[[s + 1]]  <- pl_c
        gradient_list[[s + 1]] <- -(colSums(Dl_c) + aug_s)
        D2log_list[[s + 1]]    <- colSums(D2l_c)
      }
      sym <- function(x)
        matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
                     PACKAGE="mets")$XXf, np, np)
      hessian <- lapply(D2log_list, sym)
      if (all) {
        ihess <- lapply(hessian, pinv)
        beta.iid <- do.call(cbind, lapply(seq_len(nstrata), function(s) {
          idx_c <- which(strata != (s - 1))
          X_c   <- X[idx_c,  , drop=FALSE]
          Y_c   <- Y[idx_c]
          w_c   <- weights[idx_c]
          off_c <- offset[idx_c]
          id_c  <- id[idx_c]
          lp_c  <- as.vector(X_c %*% pp[[s]]) + off_c
          Dl_c <- switch(model,
            exp = {
              pr_c <- exp(lp_c)
              if (monotone) w_c * X_c * c(Y_c - pr_c)
              else          w_c * X_c * c((Y_c - pr_c) * pr_c)
            },
            logit = {
              pr_c <- expit(lp_c)
              if (monotone) w_c * X_c * c(Y_c - pr_c)
              else {
                varp <- pr_c / (1 + exp(lp_c))
                w_c * X_c * c((Y_c - pr_c) * varp)
              }
            },
            stop("unknown model\n")
          )
          raw <- Dl_c %*% ihess[[s]]
          Msumstrata(raw, id_c, max(id) + 1)
        }))
        robvar <- crossprod(beta.iid)
        return(list(
          par      = pp,             coef      = unlist(pp),
          ploglik  = ploglik_list,   gradient  = gradient_list,
          hessian  = hessian,        ihessian  = ihess,
          id       = id,             Dlogl     = NULL,
          iid      = beta.iid,       robvar    = robvar,
          var      = robvar,         se.robust = diag(robvar)^.5
        ))
      }
      return(structure(ploglik_list, gradient = gradient_list, hessian = hessian))
    }
    ## ---------------- Standard non-CV ----------------
    loffset <- eta <- numeric(length(Y))
    for (s in 0:(nstrata - 1)) {
      idx          <- strata == s
      eta[idx]     <- as.vector(X[idx, , drop=FALSE] %*% pp[[s + 1]])
      loffset[idx] <- offset[idx]
    }
    lp <- eta + loffset
    if (model == "exp") {
      pr     <- exp(lp)
      if (monotone) {
        Dlogl  <- weights * X * c(Y - pr)
        D2logl <- c(weights * pr) * X2
      } else {
        Dlogl  <- weights * X * c((Y - pr) * pr)
        D2logl <- c(weights * pr^2) * X2
      }
    } else if (model == "logit") {
      pr <- expit(lp)
      if (monotone) {
        Dlogl  <- weights * X * c(Y - pr)
        D2logl <- c(weights * pr / (1 + exp(lp))) * X2
      } else {
        varp   <- pr / (1 + exp(lp))
        Dlogl  <- weights * X * c((Y - pr) * varp)
        D2logl <- c(weights * varp^2) * X2
      }
    } else stop("obj() general path only supports logit/exp\n")
    ploglik  <- matrix(0, nstrata, 1)
    D2log    <- Msumstrata(D2logl,  strata, nstrata)
    gradient <- Msumstrata(Dlogl,   strata, nstrata) + augmentation
    if (!is.matrix(D2log))    D2log    <- matrix(D2log,    1, length(D2log))
    if (!is.matrix(gradient)) gradient <- matrix(gradient, 1, np)
    gradient <- split(-gradient, seq_len(nrow(gradient)))
    ploglik  <- split(-ploglik,  seq_len(nrow(ploglik)))
    D2log    <- split(D2log,     seq_len(nrow(D2log)))
    sym <- function(x)
      matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
                   PACKAGE="mets")$XXf, np, np)
    hessian <- lapply(D2log, sym)
    if (all) {
      ihess    <- lapply(hessian, pinv)
      iid_list <- lapply(seq_len(nstrata), function(s) {
        idx  <- strata == (s - 1)
        id_s <- id[idx]
        Dl_s <- Dlogl[idx, , drop=FALSE]
        raw  <- Dl_s %*% ihess[[s]]
        Msumstrata(raw, id_s, max(id) + 1)
      })
      beta.iid <- do.call(cbind, iid_list)
      robvar   <- crossprod(beta.iid)
      return(list(
        par      = pp,           coef      = unlist(pp),
        ploglik  = ploglik,      gradient  = gradient,
        hessian  = hessian,      ihessian  = ihess,
        id       = id,           Dlogl     = Dlogl,
        iid      = beta.iid,     robvar    = robvar,
        var      = robvar,       se.robust = diag(robvar)^.5
      ))
    }
    structure(ploglik, gradient = gradient, hessian = hessian)
  } # }}}

  ## ------------------------------------------------------------------
  ## make_stratum_obj: per-stratum closure for NR (logit/exp only)
  ## ------------------------------------------------------------------
  make_stratum_obj <- function(s, pp) { # {{{
    idx       <- if (cv.fold) strata != s else strata == s
    X_s       <- X[idx,  , drop=FALSE]
    X_s2      <- X2[idx, , drop=FALSE]
    Y_s       <- Y[idx]
    weights_s <- weights[idx]
    offset_s  <- offset[idx]
    aug_s     <- augmentation[s + 1, ]
    sym <- function(x, np)
      matrix(.Call("XXMatFULL", matrix(x, nrow=1), np,
                   PACKAGE="mets")$XXf, np, np)
    function(theta) {
      lp <- as.vector(X_s %*% theta) + offset_s
      if (model == "exp") {
        pr <- exp(lp)
        if (monotone) {
          Dlogl  <- weights_s * X_s * c(Y_s - pr)
          D2logl <- c(weights_s * pr) * X_s2
        } else {
          Dlogl  <- weights_s * X_s * c((Y_s - pr) * pr)
          D2logl <- c(weights_s * pr^2) * X_s2
        }
      } else if (model == "logit") {
        pr <- expit(lp)
        if (monotone) {
          Dlogl  <- weights_s * X_s * c(Y_s - pr)
          D2logl <- c(weights_s * pr / (1 + exp(lp))) * X_s2
        } else {
          varp   <- pr / (1 + exp(lp))
          Dlogl  <- weights_s * X_s * c((Y_s - pr) * varp)
          D2logl <- c(weights_s * varp^2) * X_s2
        }
      } else stop("make_stratum_obj only supports logit/exp\n")
      np   <- length(theta)
      grad <- -(colSums(Dlogl) + aug_s)
      D2   <- colSums(D2logl)
      hess <- sym(D2, np)
      structure(0, gradient = grad, hessian = hess)
    }
  } # }}}

  ## ------------------------------------------------------------------
  ## fast_lin_fit(): closed-form path for model=="lin". No per-stratum
  ## subsetting of X/X2/Y/weights -- one grouped pass gives every
  ## stratum's own-sums; cv.fold complements are total - own; iid
  ## residuals for all strata computed via one n x nstrata matmul.
  ## ------------------------------------------------------------------
  fast_lin_fit <- function() { ## {{{ 
    np  <- ncol(X)
    wX2 <- c(weights) * X2
    wXY <- c(weights) * X * c(Y)

    own_XtWX <- Msumstrata(wX2, strata, nstrata)
    own_XtWY <- Msumstrata(wXY, strata, nstrata)
    tot_XtWX <- colSums(wX2)
    tot_XtWY <- colSums(wXY)

    if (!is.matrix(own_XtWX)) own_XtWX <- matrix(own_XtWX, 1, length(own_XtWX))
    if (!is.matrix(own_XtWY)) own_XtWY <- matrix(own_XtWY, 1, np)

    sym <- function(v) matrix(.Call("XXMatFULL", matrix(v, nrow=1), np,
                                     PACKAGE="mets")$XXf, np, np)

    hessian <- ihess <- beta_opt <- vector("list", nstrata)
    for (s in seq_len(nstrata)) {
      vecXtWX <- if (cv.fold) tot_XtWX - own_XtWX[s, ] else own_XtWX[s, ]
      XtWY    <- if (cv.fold) tot_XtWY - own_XtWY[s, ] else own_XtWY[s, ]
      XtWY    <- XtWY + augmentation[s, ]
      H       <- sym(vecXtWX)
      iH      <- pinv(H)
      hessian[[s]]  <- H
      ihess[[s]]    <- iH
      beta_opt[[s]] <- as.vector(iH %*% XtWY)
    }

    ## residuals for ALL strata betas computed at once: n x nstrata
    BetaMat <- do.call(cbind, beta_opt)
    LP      <- X %*% BetaMat
    R       <- c(Y) - LP

    mid <- max(id)
    iid_list <- vector("list", nstrata)
    for (s in seq_len(nstrata)) {
      mask <- if (cv.fold) (strata != (s - 1)) else (strata == (s - 1))
      Dl_s <- weights * X * (R[, s] * mask)
      raw  <- Dl_s %*% ihess[[s]]
      iid_list[[s]] <- Msumstrata(raw, id, mid + 1)
    }
    beta.iid <- do.call(cbind, iid_list)
    robvar   <- crossprod(beta.iid)

    list(par = beta_opt, coef = unlist(beta_opt),
         ploglik  = lapply(seq_len(nstrata), function(s) 0),
         gradient = lapply(beta_opt, function(b) rep(0, length(b))),
         hessian = hessian, ihessian = ihess,
         id = id, Dlogl = NULL,
         iid = beta.iid, robvar = robvar, var = robvar,
         se.robust = diag(robvar)^.5)
  } ## }}} 

  ## ------------------------------------------------------------------
  ## Optimise
  ## ------------------------------------------------------------------
  dots    <- list(...)
  control <- if (length(dots) == 0) {
    if (model == "exp") list(tol=1e-10, stepsize=0.5) else NULL
  } else dots[[1]]
  if (!is.matrix(beta)) beta <- matrix(beta, nstrata, p, byrow = TRUE)
  beta <- split(beta, seq_len(nrow(beta)))
  NR_list <- function() {
    Map(function(s) {
      lava::NR(start   = beta[[s + 1]],
                obj     = make_stratum_obj(s, beta),
                control = control)
    }, seq_len(nstrata) - 1L)
  }
  beta_opt <- beta   # fallback
  if (p > 0) {
    if (!no.opt) {
      if (model == "lin") {
        val      <- fast_lin_fit()
        beta_opt <- val$par
      } else if (tolower(method) == "nr") {
        tim      <- NR_list()
        beta_opt <- lapply(tim, function(r) r$par)
        cc       <- unlist(beta_opt)
        val      <- c(list(coef = cc), obj(beta_opt, all = TRUE))
        val$gradient <- lapply(tim, function(r) r$gradient)
      } else {
        stop("only NR of lava\n")
      }
    } else {
      if (model == "lin") {
        val <- fast_lin_fit()   # still uses supplied beta? no -- no.opt means
                                 # skip optimisation; keep supplied beta instead:
        val$par  <- beta
        val$coef <- unlist(beta)
        beta_opt <- beta
      } else {
        val <- c(list(coef = unlist(beta)), obj(beta, all = TRUE))
      }
    }
  } else {
    val <- if (model == "lin") fast_lin_fit() else obj(as.list(rep(0, nstrata)), all = TRUE)
  }

  ## ------------------------------------------------------------------
  ## Coefficient naming
  ## ------------------------------------------------------------------
  if (nstrata == 1 && length(val$coef) == length(colnames(X))) {
    names(val$coef) <- colnames(X)
  } else {
    cnames <- paste(
      rep(paste0("s", 1:nstrata), each = ncol(X)),
      rep(colnames(X), times = nstrata),
      sep = "."
    )
    names(val$coef)   <- cnames
    colnames(val$iid) <- cnames
  }
  ## ------------------------------------------------------------------
  ## Assemble output
  ## ------------------------------------------------------------------
  val <- c(val, list(
    n        = length(Y),
    ncluster = nid
  ))
  val$call      <- cl
  val$call.id   <- call.id
  val$id        <- orig.id
  val$name.id   <- name.id
  val$nid       <- nid
  val$iid.naive <- val$iid
  val$coef_list <- beta_opt
  if (!is.null(name.id)) val$iid <- nameme(val$iid, name.id)
  robvar        <- crossprod(val$iid)
  val$var       <- val$robvar <- robvar
  val$se.robust <- diag(robvar)^.5
  val$se.coef   <- diag(val$var)^.5

  val$augmentation <- augmentation
  val$model        <- model
  val$outcome      <- "identity"      # for print/summary compatibility
  val$nevent       <- NA              # not applicable outside survival outcomes
  val$cause        <- NULL
  val$Y            <- Y
  val$design       <- des
  val$nstrata      <- nstrata
  val$strata_orig  <- strata
  val$strata_call  <- des.strata
  val$cv.fold      <- cv.fold
  if (!se) {
    val$iid <- val$iid.naive <- val$robvar <- val$var <-
      val$se.robust <- val$se.coef <- NULL
  }
  if (low.memory) {
    val$design$y      <- NULL
    val$design$data   <- NULL
    val$design$x      <- NULL
    val$design$strata <- NULL
    val$strata_orig   <- NULL
    val$strata_call   <- NULL
    val$iid.naive     <- NULL
  }
  class(val) <- c("binregStrata","binreg")
  return(val)
} ## }}} 


