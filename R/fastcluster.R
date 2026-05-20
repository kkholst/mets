##' Fast Cluster Index Conversion
##'
##' Converts a cluster variable to consecutive numeric indices using compiled code.
##'
##' @param x integer vector of cluster identifiers.
##' @param ... additional arguments (not used).
##' @return Integer vector of 0-based cluster indices.
##' @export
fast.cluster <- function(x,...) {
    arglist <- list("FastCluster",
                    time=as.integer(x),
                    PACKAGE="mets")
    res <- do.call(".Call",arglist)
    return(as.vector(res))
}
