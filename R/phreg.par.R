###{{{ Weibull

info_phreg_weibull <- function(...) {
    list(npar = 2,
         start = c(-1, -1),
         name = "weibull",
         partrans = function(theta) {
             exp(theta)
         },
         cumhaz = function(t, theta) {
             (theta[1]*t)^theta[2]
         }
         )
}

logl_phreg_weibull <- function(theta, time, status,
                               X=NULL, theta.idx=NULL,
                              indiv=FALSE) {
  if (!is.null(theta.idx)) {
    offsets <- which(is.na(theta.idx))
    theta <- theta[theta.idx]
    theta[offsets] <- 1
  }
  lambda <- exp(theta[1])
  p <- exp(theta[2])
  if (is.null(X)) {
    eta <- 0
  } else {
    beta <- theta[-c(1:2)]
    eta <- X%*%beta
  }

  val <- status*log(lambda*p) +
    status*(p-1)*log(lambda*time) +
    status*eta - (lambda*time)^p*exp(eta)
  if (indiv)
    return(val)
  sum(val)
}

obj_phreg_weibull <- function(...)
  -logl_phreg_weibull(...)

score_phreg_weibull <- function(theta, time, status,
                                X=NULL, theta.idx=NULL, indiv=FALSE) {
  if (!is.null(theta.idx)) {
    offsets <- which(is.na(theta.idx))
    theta <- theta[theta.idx]
    theta[offsets] <- 1
  }
  lambda <- exp(theta[1])
  p <- exp(theta[2])
  lambdaT <- lambda*time
  loglambdaT <- log(lambdaT)
  lambdaTp <- exp(loglambdaT*p)

  if (is.null(X)) {
    eta <- 0
    expeta <- 1
    dbeta <- NULL
  } else {
    beta <- theta[-c(1:2)]
    eta <- X%*%beta
    expeta <- exp(eta)
    dbeta <- ((status-expeta*lambdaTp) %x% rbind(rep(1, NCOL(X))))*X
  }
  dp <- status*(1/p + loglambdaT) - loglambdaT*lambdaTp*expeta
  dlogp <- p*dp
  dlambda <- status*(p/lambda) - p*lambdaTp/lambda*expeta
  dloglambda <- lambda*dlambda
  S <- cbind(dloglambda, dlogp, dbeta)
  if (!is.null(theta.idx)) {
    u.idx <- na.omit(unique(theta.idx))
    newS <- matrix(0, ncol=length(u.idx), nrow=nrow(S))
    for (i in u.idx) {
      newS[, i] <- cbind(rowSums(S[, which(theta.idx==i), drop=FALSE]))
    }
    S <- newS
  }
  if (indiv)
    return(S)
  colSums(S)
}

hessian_phreg_weibull <- function(theta, time, status,
                                  X=NULL, theta.idx=NULL,
                                  all=FALSE) {
  if (!is.null(theta.idx)) {
    offsets <- which(is.na(theta.idx))
    theta <- theta[theta.idx]
    theta[offsets] <- 1
  }
  lambda <- exp(theta[1])
  p <- exp(theta[2])
  lambdaT <- lambda*time
  loglambdaT <- log(lambdaT)
  Tp <- time^p
  lambdaTp <- lambda^p*Tp
    if (is.null(X)) {
    eta <- 0
    expeta <- 1
    d2dlogpdbeta <- d2dlogpdbeta <- d2beta <- NULL
  } else {
    beta <- theta[-c(1:2)]
    eta <- X %*% beta
    expeta <- exp(eta)
    ## D(beta,beta)
    U <- ((expeta*Tp) %x% rbind(rep(1, NCOL(X))))*X
    d2beta <- -t(lambda^p*U) %*% X
    ## D(p,beta)
    d2dpdbeta <-  colSums(-((loglambdaT*lambdaTp*expeta) %x%
                            rbind(rep(1, NCOL(X))))*X)
    d2dlogpdbeta <- d2dpdbeta*p
    ## D(lambda,beta)
    d2dlambdadbeta <- -U*p*lambda^(p-1)
    d2dloglambdadbeta <- colSums(d2dlambdadbeta)*lambda
  }
  ## D(p,p)
  dp <- status*(1/p + loglambdaT) - loglambdaT*lambdaTp*expeta
  d2p <- -sum(status/(p^2) + loglambdaT^2*expeta*lambdaTp)
  dlogp <- p*dp
  d2logp <- sum(dlogp)+p^2*d2p
  ## D(lambda,lambda)
  dlambda <- status*(p/lambda) - p*lambdaTp/lambda*expeta
  d2lambda <- -sum(status*(p/lambda^2) + p*(p-1)*lambdaTp/(lambda^2)*expeta)
  dloglambda <- lambda*dlambda
  d2loglambda <- sum(dloglambda)+lambda^2*d2lambda
  ## D(p,lambda)
  d2dpdlambda <- -status/(p^2) - loglambdaT^2*expeta*lambdaTp
  d2dpdlambda <- status/lambda - lambdaTp/lambda*expeta -
    p*loglambdaT*lambdaTp/lambda*expeta
  d2dlogpdloglambda <- sum(d2dpdlambda)*p*lambda
  ## Hessian:
  H <- matrix(0, length(theta), length(theta))
  H[1, 1] <- d2loglambda
  H[2, 2] <- d2logp
  H[1, 2] <- H[2, 1] <- d2dlogpdloglambda
  if (!is.null(X)) {
    H[3:length(theta), 3:length(theta)] <- d2beta
    H[2, 3:length(theta)] <- H[3:length(theta), 2] <- d2dlogpdbeta
    H[1, 3:length(theta)] <- H[3:length(theta), 1] <- d2dloglambdadbeta
  }
  if (!is.null(theta.idx)) {
    u.idx <- na.omit(unique(theta.idx))
    newH <- matrix(0, length(u.idx), length(u.idx))
    for (i in u.idx) {
      for (j in u.idx) {
        newH[i, j] <- sum(H[which(theta.idx==i), which(theta.idx==j)])
      }
    }
    H <- newH
  }
  if (all) {
    ## Score:
    if (is.null(X)) dbeta <- NULL else dbeta <- status*X-lambda^p*U
    S <- cbind(dloglambda, dlogp, dbeta)
    if (!is.null(theta.idx)) {
      u.idx <- na.omit(unique(theta.idx))
      newS <- matrix(0, ncol=length(u.idx), nrow=nrow(S))
      for (i in u.idx) {
        newS[, i] <- cbind(rowSums(S[, which(theta.idx==i), drop=FALSE]))
      }
      S <- newS
    }
    attributes(H)$grad <- colSums(S)
    attributes(H)$score <- S
    ## LogLik
    attributes(H)$logL <- sum(status*log(lambda*p) +
                              status*(p-1)*loglambdaT +
                              status*eta - lambdaTp*expeta)
  }
  return(H)
}

###}}} Weibull

###{{{ Generalized-Gamma


## http://www.stanford.edu/~lutian/coursepdf/unit1.pdf

gengamma_phreg_f <- function(t, p, lambda, alpha, ...) {
    p*lambda*(lambda*t)^(alpha-1)*exp(-(lambda*t)^p)/gamma(alpha/p)
}

## incomplete gamma gamma(s,x) = pgamma(x,s)
gengamma_phreg_F <- function(t, p, lambda, alpha, ...) {
    pgamma((lambda*t)^p, alpha/p)/gamma(alpha/p)
}

## incomplete gamma gamma(s,x) = pgamma(x,s)
gengamma_phreg_h <- function(t, p, lambda, alpha, ...) {
    p*lambda*(lambda*t)^(alpha-1) *
      exp(-(lambda*t)^p)/pgamma((lambda*t)^p, alpha/p)
}

logl_phreg_gengamma <- function(theta, time, status, X=NULL, indiv=FALSE, ...) {
    ## suppressMessages(browser())
  p <- exp(theta[1])
  lambda <- exp(theta[2])
  alpha <- exp(theta[3])
  eta <- 0
  if (!is.null(X)) {
    beta <- theta[seq(length(theta)-3)+3]
    eta <- X %*% beta
  }
  res <- log(gengamma_phreg_f(time, p, lambda, alpha)) * status +
    status*eta +
    log(1-gengamma_phreg_F(time, p, lambda, alpha) * exp(eta))
  if (indiv) return(res)
  return(sum(res))
}

score_phreg_gengamma <- function(theta, ...) {
    numDeriv::jacobian(logl_phreg_gengamma, theta, ...)
}

hessian_phreg_gengamma <- function(theta, ...) {
    numDeriv::hessian(logl_phreg_gengamma, theta, ...)
}

info_phreg_gengamma <- function(...) {
    list(npar=3, start=c(-1, -1, -1), name="gengamma")
}

###}}} Generalized-Gamma

##' @export
predict.phreg.par <- function(object, p=coef(object),
                              X=object$X, time=object$time, ...) {
    info <- do.call(paste("info_phreg", object$model, sep="_"), list())
    cc <- coef(object)
    eta <- 0
    if (length(cc)>info$npar) {
        eta <- X %*% p[-seq(info$npar)]
    }
    exp(-(info$cumhaz(time, info$partrans(p))*exp(eta)))
}



###{{{ phreg.par + methods

##' @export
phreg.par <- function(formula, data=parent.frame(),
                      time, status, X=NULL, model="weibull",
                      theta.idx=NULL, theta0,
                      niter=100, tol1=1e-9, tol2=1e-9,
                      lambda1=0.5, lambda2=1, trace=0, ...) {

    if (!missing(formula)) {
        cl <- match.call()
        m <- match.call(expand.dots = TRUE)[1:3]
        special <- c("strata", "cluster")
        Terms <- terms(formula, special, data = data)
        m$formula <- Terms
        m[[1]] <- as.name("model.frame")
        m <- eval(m, parent.frame())
        Y <- model.extract(m, "response")
        if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
        if (ncol(Y)==2) {
            exit <- eval(Y[, 1], data)
            entry <- NULL ## rep(0,nrow(Y))
            status <- Y[, 2]
        } else {
            entry <- Y[, 1]
            exit <- Y[, 2]
            status <- Y[, 3]
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
        }

        X <- model.matrix(Terms, m)
        if (!is.null(intpos  <- attributes(Terms)$intercept))
            X <- X[, -intpos, drop=FALSE]
        if (ncol(X)==0) X <- NULL
        time <- exit
    }

    myinfo <- do.call(paste("info_phreg", model, sep="_"), list())
    if (is.null(X)) {
        theta0 <- myinfo$start
    } else {
        if (missing(theta0))
          theta0 <- rep(0,
                        ifelse(is.null(theta.idx),
                               myinfo$npar+NCOL(X),
                               length(unique(na.omit(theta.idx)))))
    }

  hess <- paste("hessian_phreg", model, sep="_")
  scor <- paste("score_phreg", model, sep="_")
  logl <- paste("logl_phreg", model, sep="_")
  thetas <- theta0
  logL <- c()
  for (i in 1:niter) {
      H <- do.call(hess, list(theta=theta0,
                              all=TRUE,
                              time=time,
                              status=status,
                              X=X,
                              theta.idx=theta.idx))
      if (!is.null(attributes(H)$score)) {
          S <- colSums(attributes(H)$score)
      } else {
          S <- do.call(scor, list(theta=theta0, time=time, status=status, X=X))
      }
      gamma <- 0
      if (lambda2>0) gamma <- lambda2*sqrt((t(S) %*% S)[1]) * diag(NROW(H))
      theta0 <- theta0 - lambda1*solve(H-gamma) %*% S
      thetas <- rbind(thetas, as.vector(theta0))
      if (!is.null(attributes(H)$logL))  {
          logL <- c(logL, attributes(H)$logL)
      } else {
          logL <- c(logL, do.call(logl,
                                  list(theta=theta0,
                                       time=time,
                                       status=status,
                                       X=X)))
      }
      if(trace>0)
          if (i%%trace==0) {
              cat("Iter=", i, ", logLik=", tail(logL, 1), "\n", sep="")
              cat("theta=(", paste(formatC(theta0),
                                   collapse=";"), ")\n", sep="")
          }
      if (i>1)
          if (sum(abs(S^2))<tol1 && abs(logL[i]-logL[i-1])<tol2) break
  }

  V <- solve(-H)
  coefs <- cbind(theta0, diag(V)^0.5)
  colnames(coefs) <- c("Estimate", "Std.Err")
  attributes(coefs)$score <- S
  attributes(coefs)$logLik <- tail(logL, 1)
  mynames <- c("log(scale)", "log(shape)")
  if (!is.null(X)) {
    xnames <- colnames(X)
    if (is.null(xnames))
      xnames <- paste("x", seq_len(NCOL(X)), sep="")
    mynames <- c(mynames, xnames)
  }
  if (!is.null(theta.idx)) {
    mynames <- mynames[na.omit(unique(theta.idx))]
  }
  rownames(coefs) <- mynames
  structure(list(call=match.call(), coef=coefs[, 1],
                 coefmat=coefs, vcov=V,
                 formula=eval(formula),
                 nevent=sum(status),
                 model=model,
                 time=time, status=status, X=X),
            class=c("phreg.par", "phreg"))
}

##' @export
print.phreg.par <- function(x, ...) {
    cat("Call:\n")
    dput(x$call)
    printCoefmat(x$coefmat, ...)
}

##' @export
vcov.phreg.par <- function(object, ...) {
    object$vcov
}

##' @export
coef.phreg.par <- function(object, ...) {
    object$coef
}

##' @export
IC.phreg.par <- function(x, p=coef(x), ...) {
    iI <- vcov(x)
    U <- do.call(paste("score_phreg", x$model, sep="_"),
                 list(theta=p, time=x$time, status=x$status,
                      X=x$X, indiv=TRUE))
    res <- (U)%*%iI*NROW(U)
    colnames(res) <- names(coef(x))
    structure(res, iI=iI)
}

##' @export
logLik.phreg.par <- function(object, p=coef(object), ...) {
  do.call(paste("logl_phreg", object$model, sep="_"),
          list(theta=p, object$time, object$status,
               object$X, ...))
}

##' @export
model.frame.phreg.par <- function(formula, ...) {
    formula$X
}

##' @export
predict.phreg.par <- function(object,
                              p = coef(object),
                              X = object$X,
                              time = object$time, ...) {
  info <- do.call(paste("info_phreg",
                        object$model, sep="_"), list())
    cc <- coef(object)
    eta <- 0
    if (length(cc)>info$npar) {
        eta <- X %*% p[-seq(info$npar)]
    }
    exp(-(info$cumhaz(time, info$partrans(p)) * exp(eta)))
}
