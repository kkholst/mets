##' Piecewise constant hazard distribution
##'
##' Piecewise constant hazard distribution
##' @aliases rpch ppch
##' @export
##' @examples
##' 
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

rpch <- function(n, lambda=1, breaks=c(0,Inf)) {
    res <- .Call("_mets_rpch",
                 n=n,
                 lambda=lambda,
                 time=breaks)
    return(res)
}

