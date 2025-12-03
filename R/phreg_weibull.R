pred_weibull <- function(object, X, Z,
                         times,
                         individual.times = FALSE,
                         time.fun = NULL,
                         type = c("surv", "haz", "chaz", "lp"), ...) {
    p <- coef(object)
    V <- vcov(object)
    X <- rbind(X)
    Z <- rbind(Z)
    p1 <- X %*% p[seq_len(ncol(X))] # log-scale
    p2 <- Z %*% p[seq_len(ncol(Z)) + ncol(X)] # log-shape
    if (type[1] == "lp") {
      return(cbind(lograte = p1, logshape = p2))
    }
    if (nrow(X) != nrow(Z)) stop("incompatible X and Z matrix")
    n <- nrow(X)
    np <- ncol(X) + ncol(Z)
    xz <- matrix(0, ncol = np, nrow = 2)
    var_arr <- array(dim = c(2, 2, n))
    for (i in seq_len(n)) {
      xz[1, seq_len(ncol(X))] <- X[i, ]
      xz[2, seq_len(ncol(Z)) + ncol(X)] <- Z[i, ]
      var_arr[, , i] <- xz %*% V %*% t(xz)
    }
    par <- cbind(p1, p2)
    nt <- length(times)
    if (individual.times && nt != n) {
      stop("For individual time predictions 'times', 'X' and 'Z' should agree") # nolint
    }
    if (is.null(time.fun)) {
      time.fun <- structure(identity, grad = identity)
    }
    chaz <- function(par, ...) exp(par[1]) * (time.fun(newtime)**exp(par[2]))
    surv <- function(par, ...) exp(-chaz(par, ...))
    haz <- function(par) {
      res <- exp(par[1] + par[2] +
          (exp(par[2]) - 1) * log(time.fun(newtime)) +
          log(attr(time.fun, "grad")(newtime)))
      res[newtime == 0] <- 0
      return(res)
    }
    if (individual.times) {
      est_arr <- array(dim = c(nt, 4))
    } else {
      est_arr <- array(dim = c(nt, 4, n))
      newtime <- times
    }
    for (i in seq_len(n)) {
      if (individual.times) {
        newtime <- times[i]
      }
      pr <- estimate(
        coef = par[i, ],
        vcov = var_arr[, , i], f = get(type[1]), ...
        )$coefmat
      if (individual.times) {
        est_arr[i, ] <- pr[1:4]
      } else {
        est_arr[, , i] <- pr[, 1:4]
      }
    }
    if (individual.times) {
      dimnames(est_arr) <- list(seq_len(nt), colnames(pr)[1:4])
    } else {
      dimnames(est_arr) <- list(
        paste0("t", seq_len(nt)),
        colnames(pr)[1:4], seq_len(n)
      )
    }
  return(est_arr)
}

logl_weibull <- function(p, entry, exit, status,
                         X = NULL, Z = NULL,
                         time.fun = NULL) {
    if (is.null(X)) X <- cbind(rep(1, length(exit)))
    if (is.null(Z)) Z <- cbind(rep(1, length(exit)))
    dt <- 1
    if (!is.null(time.fun)) {
      if (is.null(attr(time.fun, "grad"))) {
          # complex-step derivative
          dt <- numDeriv::grad(time.fun, exit, method = "complex")
      } else {
          # numerical derivative
          dt <- attr(time.fun, "grad")(exit)
      }
      exit <- time.fun(exit)
      entry <- time.fun(entry)
    }
    if (length(p) != ncol(X) + ncol(Z)) stop("wrong parameter length")
    p1 <- X %*% p[seq_len(ncol(X))] # log-scale
    p2 <- Z %*% p[seq_len(ncol(Z)) + ncol(X)] # log-shape
    theta <- cbind(exp(p1), exp(p2))
    loghaz <- log(theta[, 1]) + log(theta[, 2]) +
      (theta[, 2] - 1) * log(exit) + log(dt)
    chaz <- function(t) theta[, 1] * t**theta[, 2]
    logl <- status * loghaz - chaz(exit) + chaz(entry)
    logl
}

score_weibull <- function(p, entry, exit, status,
                          X = NULL, Z = NULL,
                          time.fun = NULL) {
    if (is.null(X)) X <- cbind(rep(1, length(exit)))
    if (is.null(Z)) Z <- cbind(rep(1, length(exit)))
    if (!is.null(time.fun)) {
        entry = time.fun(entry)
        exit = time.fun(exit)
    }
    if (length(p) != ncol(X) + ncol(Z)) stop("wrong parameter length")
    p1 <- X %*% p[seq_len(ncol(X))]
    p2 <- Z %*% p[seq_len(ncol(Z)) + ncol(X)]
    theta <- exp(cbind(p1, p2))
    dloghaz.dtheta1 <- status / theta[, 1]
    dchaz.dtheta1 <- function(t) {
        idx0 <- which(t == 0)
        res <- t**theta[, 2]
        res[idx0] <- 0 # conv. log(x)*x =0 when x=0
        return(res)
    }
    dloghaz.dtheta2 <- status / theta[, 2] + status * log(exit)
    dchaz.dtheta2 <- function(t) {
        idx0 <- which(t == 0)
        res <- theta[, 1] * log(t) * t**theta[, 2]
        res[idx0] <- 0 # conv. log(x)*x =0 when x=0
        return(res)
    }
    res <- cbind(
        dloghaz.dtheta1 - dchaz.dtheta1(exit) + dchaz.dtheta1(entry),
        dloghaz.dtheta2 - dchaz.dtheta2(exit) + dchaz.dtheta2(entry)
    )
    for (i in seq_len(ncol(X))) X[, i] <- X[, i] * res[, 1] * theta[, 1]
    for (i in seq_len(ncol(Z))) Z[, i] <- Z[, i] * res[, 2] * theta[, 2]
    return(cbind(X, Z))
}

##' @description Fits a Cox-Weibull with cumulative hazard given by \deqn{
##'   \Lambda(t) = \lambda \cdot t^s } where \eqn{s} is the shape parameter, and
##'   \eqn{\lambda} the rate parameter. We here allow a regression model for
##'   both parameters \deqn{\lambda := \exp(\beta^\top X)} \deqn{s :=
##'   \exp(\gamma^\top Z)} as defined by `formula` and `shape.formula`
##'   respectively.
##' @details The parametrization
##' @title Weibull-Cox regression
##' @param formula Formula for proportional hazards. The right-handside must be
##'   an [Event] or [Surv] object (with right-censoring and possibly delayed
##'   entry).
##' @param shape.formula Formula for shape parameter
##' @param data data.frame
##' @param time.fun optional smooth function specifying transformation of
##'   time-variable. Here the cumulative hazard is assumed to be
##'   \eqn{H(t)=\lambda g(t)^s}
##' @param save.data if TRUE the data.frame is stored in the model object (for
##'   predictions and simulations)
##' @param control control arguments to optimization routine [stats::nlmbin]
##' @seealso [mets::phreg()]
##' @author Klaus Kähler Holst, Thomas Scheike
##' @examples
##' data(sTRACE, package="mets")
##' sTRACE$entry <- 0
##' fit1 <- phreg_weibull(Event(entry, time, status == 9) ~ age,
##'              shape.formula = ~age, data = sTRACE)
##' tt <- seq(0,10, length.out=100)
##' pr1 <- predict(fit1, newdata = sTRACE[1, ], times = tt)
##' fit2 <- phreg(Event(time, status == 9) ~ age, data = sTRACE)
##' pr2 <- predict(fit2, newdata = sTRACE[1, ], se = FALSE)
##' if (interactive()) {
##'    plot(pr2$times, pr2$surv, type="s")
##'    lines(tt, pr1[,1,1], col="red", lwd=2)
##' }
##' @export
##' @return `phreg.par` object
phreg_weibull <- function(formula,
                          shape.formula = ~1,
                          data,
                          time.fun = NULL,
                          save.data = TRUE,
                          control = list()) {
    cl <- match.call()
    des <- proc_design(
      formula,
      data = data,
      specials = c("offset", "weights", "cluster"),
      intercept = TRUE
    )
    des2 <- proc_design(shape.formula, data = data, intercept = TRUE)
    Y <- des$y
    if (!inherits(Y, c("Event", "Surv"))) {
        stop("Expected a 'Surv' or 'Event'-object")
    }
    if (ncol(Y) == 2) {
        exit <- Y[, 1]
        entry <- rep(0, nrow(Y))
        status <- Y[, 2]
    } else {
        entry <- Y[, 1]
        exit <- Y[, 2]
        status <- Y[, 3]
    }
    X <- des$x
    Z <- des2$x
    obj <- function(p) {
        -sum(logl_weibull(p, entry, exit, status,
            X = X, Z = Z, time.fun = time.fun
        ))
    }
    grad <- function(p, indiv = FALSE) {
        U <- -score_weibull(p, entry, exit, status,
            X = X, Z = Z, time.fun = time.fun
        )
        if (indiv) {
            return(U)
        }
        return(colSums(U))
    }
    if (!is.null(control$start)) {
        p0 <- control$start
        control$init <- NULL
    } else {
        p0 <- rep(0, ncol(X)+ncol(Z))
    }
    op <- stats::nlminb(p0, obj, grad, control = control)
    p <- op$par
    H <- numDeriv::jacobian(grad, p)
    U <- -grad(op$par, indiv=TRUE)
    n <- nrow(U)
    loglik <- structure(-obj(p),
        df = length(p), nobs = n, nall = n,
        class = "logLik"
        )
    ic <- t(solve(H) %*% t(U)) * n
    names(p) <- c(
      colnames(X),
      paste0("s:", colnames(Z))
    )
    est <- lava::estimate(coef = p, IC = ic)
    est$model.index <- list(ncol(X), ncol(X) + ncol(Z))
    res <- list(
        loglik = loglik,
        response = list(entry = entry, exit = exit, status = status),
        call = cl, opt = op, coef = p,
        hessian = H,
        rate.design = clean_design(des),
        shape.design = clean_design(des2),
        data = NULL,
        time.fun = time.fun,
        estimate = est
    )
    if (save.data) res$data <- data
    return(structure(res, class = "phreg.par"))
}

##' @export
print.phreg.par <- function(x, ...) {
    cat("\n- Weibull-Cox model -\n\n")
    cat("Call:\n")
    dput(x$call)
    cat("\n")
    cat(sprintf("log-Likelihood: %f", x$loglik), "\n\n")
    dur <- with(x$response, exit - entry)
    n <- length(x$response$status)
    events <- sum(x$response$status)
    nn <- rbind(c(n, events, sum(dur)))
    colnames(nn) <- c("n","events","obs.time")
    rownames(nn) <- ""
    print(nn, quote=FALSE)
    cat("\n")
    print(x$estimate, ...)
    cat("\n")
}

##' @export
vcov.phreg.par <- function(object, ...) {
    vcov(object$estimate)
}

##' @export
coef.phreg.par <- function(object, ...) {
    coef(object$estimate)
}

##' @export
simulate.phreg.par <- function(object,
                               nsim = nrow(data),
                               seed = NULL,
                               data = object$data,
                               cens.model,
                               var.names = c("time", "status"),
                               ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  # bootstrap covariates
  newd <- mets::dsample(size = nsim, data)
  # linear-predictors
  lp <- predict(object, newdata = newd, type = "lp")
  ## simulate event times
  time <- rweibullcox(nrow(lp),
                      rate = exp(lp[, 1]),
                      shape = exp(lp[, 2])
                      )
  if (missing(cens.model)) {
      y <- update_design(object$rate.design, data = data, response = TRUE)$y
      data$cens_ <- !y[, 2]
      data$time_ <- y[, 1]
      cens.model <- phreg_weibull(Surv(time_, cens_) ~ 1, data = data)
  }
  if (is.null(cens.model)) {
    cens <- rep(TRUE, nrow(newd))
  } else {
    cens.par <- exp(predict(cens.model, type = "lp", newdata = newd))
    cens <- rweibullcox(nsim, cens.par[, 1], cens.par[, 2])
  }
  newd[[var.names[1]]] <- pmin(time, cens)
  newd[[var.names[2]]] <- (time <= cens)
  return(newd)
}

##' @export
IC.phreg.par <- function(x, p=coef(x), ...) {
  IC(x$estimate)
}

##' @export
logLik.phreg.par <- function(object, ...) {
  object$loglik
}

##' @export
predict.phreg.par <- function(object,
                              newdata = object$data,
                              times,
                              individual.time = FALSE,
                              type = c("surv", "haz", "chaz", "lp"),
                              level = 0.05,
                              ...) {
  xy <- update_design(object$rate.design, data = newdata, response = TRUE)
  x <- xy$x
  z <- update_design(object$shape.design, data = newdata)$x
  if (type[1] != "lp" && missing(times)) {
    times <- xy$y[, 1]
  }
  pr <- pred_weibull(object,
                     X = x, Z = z,
                     times = times,
                     individual.time = individual.time,
                     level = level,
                     time.fun = object$time.fun,
                     type = type, ...
                     )
  return(pr)
}


##' @export
##' @description Simulate observations from the model with cumulative hazard
##'   given by \deqn{\Lambda(t) = \lambda\cdot t^s} where \eqn{\lambda} is the
##'   \emph{rate parameter} and \eqn{s} is the \emph{shape parameter}.
##' @details [stats::rweibull()] uses a different parametrization with
##'   cumulative hazard given by
##' \deqn{H(t) = (t/b)^a,}
##' i.e., the shape is the same \eqn{a:=s} but the scale paramter \eqn{b}
##' is related to rate paramter \eqn{r} by
##' \deqn{r := b^{-a}}
##' @seealso [stats::rweibull()]
##' @title Simulate observations from a Weibull distribution
##' @param n (integer) number of observations
##' @param rate (numeric) rate parameter (can be a vector of size n)
##' @param shape (numeric) shape parameter (can be a vector of size n)
rweibullcox <- function(n, rate, shape) {
  # F(t) = 1-exp(-H(t)) = 1-exp(-lambda t^s). Inversion=>
  (-log(1-runif(n))/rate)^(1/shape)
}
