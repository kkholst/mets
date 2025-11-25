#' Event history object
#'
#' Constructur for Event History objects
#'
#' ... content for details
#'
#' @aliases Event as.character.Event as.matrix.Event [.Event format.Event
#' print.Event rbind.Event summary.Event
#' @param time Time
#' @param time2 Time 2
#' @param cause Cause
#' @param cens.code Censoring code (default 0)
#' @param ... Additional arguments
#' @return Object of class Event (a matrix)
#' @author Klaus K. Holst and Thomas Scheike
#' @examples
#'
#' 	t1 <- 1:10
#' 	t2 <- t1+runif(10)
#' 	ca <- rbinom(10,2,0.4)
#' 	(x <- Event(t1,t2,ca))
#' @export
Event <- function(time,time2=TRUE,cause=NULL,cens.code=0,...) {# {{{
  if (missing(cause)) {
    if (is.factor(time2) | is.character(time2)) {
      stop("cause should be numeric\n")
    }
  }
  else { if (is.factor(cause) | is.character(cause))
           stop("cause should be numeric\n");  }
  out <- cbind(time, time2, cause)
  if (any(is.na(out))) warning("missing values in Event object\n")
  if (!missing(cause)) {
    colnames(out) <- c("entry", "exit", "cause")
    tmp <- (out[,1]>out[,2])
    if (any(tmp)) warning("entry time later than exit time\n")
  } else {
    colnames(out) <- c("exit", "cause")
  }
  class(out) <- "Event"
  attr(out, "cens.code") <- cens.code
  return(out)
}
# }}}

# exported in zzz.R
as.matrix.Event <- function(x, ...) structure(x, class="matrix")

# exported in zzz.R
as.data.frame.Event <- function(x, ...) {
  as.data.frame.model.matrix(x, ...)
}

# exported in zzz.R
as.character.Event <- function(x, ...) {
  if (ncol(x) == 3) {
    res <- paste(
        "(", format(x[, 1], ...), ";",
        format(x[, 2], ...), ":",
        format(x[, 3], ...), "]",
        sep = ""
    )
  } else {
    res <- paste(format(x[, 1], ...), ":", format(x[, 2], ...), sep = "")
  }
  return(res)
}

# exported in zzz.R
format.Event <- function(x, ...) format(as.character.Event(x), ...)


# exported in zzz.R
print.Event <- function(x, ...) {
    print(as.character(x), ..., quote=FALSE)
}

# exported in zzz.R
summary.Event <- function(object,...) {
  cat(paste("cens.code=",attr(object,"cens.code"),"\n"))
  cat("causes:\n")
  print(table(object[,"cause"]))
  cat("exit:\n")
  print(summary(object[,"exit"]))
  if (ncol(object)==3) {
    cat("entry:\n")
    print(summary(object[,"entry"]))
    cat("exit-entry:\n")
    print(summary(object[,"exit"]- object[,"entry"]))
  }
}

# exported in zzz.R
"[.Event" <- function (x, i, j, drop = FALSE) {
  if (missing(j)) {
    atr <- attributes(x)
    class(x) <- "matrix"
    x <- x[i, , drop = FALSE]
    class(x) <- "Event"
    atr.keep <- c("cens.code", "entry")
    attributes(x)[atr.keep] <- atr[atr.keep]
    x
  }
  else {
    class(x) <- "matrix"
    NextMethod("[")
  }
}

# exported in zzz.R
rbind.Event <- function(...) {
    dots <- list(...)
    cens.code <- attributes(dots[[1]])$cens.code
    type <- attributes(dots[[1]])$type
    ncol <- dim(dots[[1]])[2]
    nrow <- unlist(lapply(dots, nrow))
    cnrow <- c(0, cumsum(nrow))
    M <- matrix(ncol = ncol, nrow = sum(nrow))
    for (i in seq_len(length(dots))) {
        M[(cnrow[i] + 1):cnrow[i + 1], ] <- dots[[i]]
    }
    x <- c()
    for (i in seq_len(ncol(M))) x <- c(x, list(M[, i]))
    x <- c(x, list(cens.code = cens.code))
    do.call("Event", x)
}

# exported in zzz.R
length.Event <- function(x) nrow(x)

# exported in zzz.R
c.Event <- function(...) {
    objects <- list(...)
    if (!all(unlist(lapply(objects, function(x) inherits(x, "Event")))) &&
        all.equal(unlist(lapply(objects, NCOL)))) {
        stop("All elements should be `Event` of the same type")
    }
    Reduce(rbind, objects)
}
