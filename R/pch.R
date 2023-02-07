##' Piecewise constant hazard distribution
##'
##' Piecewise constant hazard distribution
##' @aliases rpch ppch
##' @export
##' @param n sample size
##' @param lambda rate parameters
##' @param breaks time cut-points
##' @aliases rpch ppch
rpch <- function(n, lambda=1, breaks=c(0,Inf)) {
    res <- .Call("_mets_rpch",
                 n=n,
                 lambda=lambda,
                 time=breaks)
    return(res)
}

##' Piecewise constant hazard distribution
##'
##' Piecewise constant hazard distribution
##' @param base1 baseline 
##' @param rr relative risk terms 
##' @param entry entry times for left truncation 
##' @export
rchazC <- function(base1, rr, entry) {
    if (sum(abs(base1[1,]))>0) base1 <- rbind(0,base1)
    res <- .Call("_mets_rchazC",as.matrix(base1),rr,entry)
    return(res)
}

