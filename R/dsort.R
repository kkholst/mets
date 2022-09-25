##' Sort data according to columns in data frame
##'
##' @title Sort data frame
##' @param data Data frame
##' @param x variable to order by
##' @param ... additional variables to order by
##' @param decreasing sort order (vector of length x)
##' @param return.order return order 
##' @return data.frame
##' @export
##' @examples
##' data(data="hubble",package="lava")
##' dsort(hubble, "sigma")
##' dsort(hubble, hubble$sigma,"v")
##' dsort(hubble,~sigma+v)
##' dsort(hubble,~sigma-v)
##'
##' ## with direct asignment
##' dsort(hubble) <- ~sigma-v
##' @aliases dsort2 dsort dsort<-
##' @export
dsort <- function(data,x,...,decreasing=FALSE,return.order=FALSE)
{# {{{
    if (missing(x)) return(data)
    if (inherits(x,"formula")) {
        xx <- lava::procformula(value=x)$res
        decr <- unlist(lapply(xx,function(x) substr(trim(x),1,1)=="-"))
        if (any(decr)) decreasing <- decr
        x <- all.vars(x)
    }
    if (is.character(x) && length(x)!=nrow(data)) x <- lapply(x,function(z) data[,z,drop=TRUE])
    dots <- list(...)
    args <- lapply(dots, function(x) {
        if (length(x)==1 && is.character(x)) x <- data[,x,drop=TRUE]
        x
    })
    if (!is.list(x)) x <- list(x)
    ord <- do.call("order",c(c(x,args),list(decreasing=decreasing,method="radix")))
    data[ord,]
}# }}}

##' @export
dsort2 <- function(data,x,...,decreasing=FALSE,return.order=FALSE,regex=FALSE)
{# {{{
    if (missing(x)) return(data)
    if (inherits(x,"formula")) {
        xx <- lava::procformula(value=x)$res
        yxzf <- procform(x,data=data,do.filter=FALSE,regex=regex)
        decr <- unlist(lapply(xx,function(x) substr(trim(x),1,1)=="-"))
        if (any(decr)) decreasing <- decr
        x <- all.vars(x)
	x <- yxzf$predictor
    }
    if (is.character(x) && length(x)!=nrow(data)) x <- lapply(x,function(z) data[,z,drop=TRUE])
    dots <- list(...)
    args <- lapply(dots, function(x) {
        if (length(x)==1 && is.character(x)) x <- data[,x,drop=TRUE]
        x
    })
    if (!is.list(x)) x <- list(x)
    ord <- do.call("order",c(c(x,args),list(decreasing=decreasing,method="radix")))
    data[ord,]
}# }}}

##' @export
"dsort<-" <- function(data,...,value)  dsort(data,value,...)
