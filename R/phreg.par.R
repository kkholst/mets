pred_weibull <- function(object, X, Z,
                         times,
                         individual.times = FALSE,
                         surv = TRUE, ...) {
    p <- coef(object)
    V <- vcov(object)
    X <- rbind(X)
    Z <- rbind(Z)
    p1 <- X %*% p[seq_len(ncol(X))] # log-scale
    p2 <- Z %*% p[seq_len(ncol(Z)) + ncol(X)] # log-shape
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
      stop("For individual time predictions 'times', 'X' and 'Z' should agree")
    }
    survf <- function(par, ...) exp(-exp(par[1]) * (newtime**exp(par[2])))
    haz <- function(par) {
        exp(par[1] + par[2] +
            (exp(par[2]) - 1) * log(times))
      }
    if (individual.times) {
      est_arr <- array(dim = c(nt, 4))
    } else {
      est_arr <- array(dim = c(nt, 4, n))
      newtime <- times
    }
    time <- times
    for (i in seq_len(n)) {
        if (individual.times) {
            newtime <- times[i]
        }
        if (surv) {
            pr <- estimate(
                coef = par[i, ],
                vcov = var_arr[, , i], f = survf, ...
            ) |>
             lava::parameter()
        } else {
            pr <- estimate(
                coef = par[i, ],
                vcov = var_arr[, , i], f = haz, ...
            ) |>
             lava::parameter()
        }
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

logl_weibull <- function(p, entry, exit, status, X = NULL, Z = NULL) {
    if (is.null(X)) X <- cbind(rep(1, length(exit)))
    if (is.null(Z)) Z <- cbind(rep(1, length(exit)))
    if (length(p) != ncol(X) + ncol(Z)) stop("wrong parameter length")
    p1 <- X %*% p[seq_len(ncol(X))] # log-scale
    p2 <- Z %*% p[seq_len(ncol(Z)) + ncol(X)] # log-shape
    theta <- cbind(exp(p1), exp(p2))
    loghaz <- log(theta[, 1]) + log(theta[, 2]) +
      (theta[, 2] - 1) * log(exit)
    chaz <- function(t) theta[, 1] * t**theta[, 2]
    logl <- status * loghaz - chaz(exit) + chaz(entry)
    logl
}

score_weibull <- function(p, entry, exit, status, X = NULL, Z = NULL) {
    if (is.null(X)) X <- cbind(rep(1, length(exit)))
    if (is.null(Z)) Z <- cbind(rep(1, length(exit)))
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

##' @description
##' Fits a Cox-Weibull with cumulative hazard given by
##' \deqn{
##'  \Lambda(t) = \alpha t^s
##' }
##' where \eqn{s} is the shape parameter, and \eqn{\alpha} the scale parameter.
##' We here allow a regression model for both parameters
##' \deqn{alpha := \beta^\top X}
##' \deqn{s := \gamma^\top Z}
##' as defined by `formula` and `shape.formula` respectively.
##' @title Weibull-Cox regression
##' @param formula Formula for proportional hazards. The right-handside must be
##'   an [Event] or [Surv] object.
##' @param shape.formula Formula for shape parameter
##' @param data data.frame
##' @param control control arguments to optimization routine [stats::nlmbin]
##' @seealso [mets::phreg()]
##' @author Klaus Kähler Holst, Thomas Scheike
##' @examples
##' data(sTRACE, package="mets")
##' sTRACE$entry <- 0
##' fit1 <- phreg_weibull(Event(entry, time, status == 9) ~ age,
##'              shape.formula = ~age, data = sTRACE)
##' tt <- seq(0,10, length.out=50)
##' pr1 <- predict(fit1, newdata = sTRACE[1, ], times = tt)
##' fit2 <- phreg(Event(time, status == 9) ~ age, data = sTRACE)
##' pr2 <- predict(fit2, newdata = sTRACE[1, ], se = FALSE)
##' if (interactive()) {
##'    plot(pr2$times, pr2$surv)
##'    points(tt, pr[,1,1], col="red", cex=0.5)
##' }
##' @export
##' @return `phreg.par` object
phreg_weibull <- function(formula,
                          shape.formula = ~1,
                          data,
                          control = list()) {
    cl <- match.call()
    des <- targeted::design(
        formula,
        data = data,
        specials = c("offset", "weights", "cluster"),
        intercept = TRUE
    )
    des2 <- targeted::design(shape.formula, data = data, intercept = TRUE)
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
        -sum(logl_weibull(p, entry, exit, status, X = X, Z = Z))
    }
    grad <- function(p, indiv=FALSE) {
      U <- -score_weibull(p, entry, exit, status, X = X, Z = Z)
      if (indiv) return(U)
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
      response = list(entry=entry, exit=exit, status=status),
      call = cl, opt = op, coef = p,
      hessian = H, design = des, shape.design = des2,
      estimate = est
    )
    return(structure(res, class = "phreg.par"))
}

###{{{ Weibull
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
    colnames(nn) <- c("n","events","exposure-time")
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
IC.phreg.par <- function(x, p=coef(x), ...) {
  IC(x$estimate)
}

##' @export
logLik.phreg.par <- function(object, ...) {
  object$loglik
}

##' @export
predict.phreg.par <- function(object,
                              newdata,
                              times,
                              individual.times = FALSE,
                              surv = TRUE,
                              level = 0.05,
                              ...) {
    x <- update(object$design, data = newdata)$x
    z <- update(object$shape.design, data = newdata)$x
    pr <- pred_weibull(object,
        X = x, Z = z,
        times = times,
        individual.times = individual.times,
        level = level,
        surv = surv, ...
    )
    return(pr)
}
