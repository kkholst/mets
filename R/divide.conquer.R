##' Generate Random Fold Indices for Cross-Validation
##'
##' Splits \code{n} observations into random folds for cross-validation.
##'
##' @param n number of observations.
##' @param folds number of folds (default 10).
##' @return A list of integer vectors, each containing indices for one fold.
##' @export
folds<- function (n, folds = 10) 
{ ## {{{ 
   split(sample(1:n), rep(1:folds, length = n))
} ## }}} 

##' Split a data set and run function 
##'
##' @title Split a data set and run function 
##' @param func called function
##' @param data data-frame
##' @param size size of splits
##' @param splits number of splits (ignored if size is given)
##' @param id optional cluster variable
##' @param ... Additional arguments to lower level functions
##' @author Thomas Scheike, Klaus K. Holst
##' @export
##' @examples
##' ## avoid dependency on timereg
##' ## library(timereg)
##' ## data(TRACE)
##' ## res <- divide_conquer(prop.odds,TRACE,
##' ## 	     formula=Event(time,status==9)~chf+vf+age,n.sim=0,size=200)
divide_conquer <- function(func=NULL,data,size,splits,id=NULL,...)
{ ## {{{
    nn <- nrow(data)
    if (!is.null(id)) {
        if (is.character(id) && length(id)==1) id <- data[,id]
        if (length(id)!=nn) stop("Wrong length of id variable")
        cc <- cluster_index(id)
        if (!missing(size)) splits <- round(cc$uniqueclust/size)
        splits <- min(splits,cc$uniqueclust)
        all.folds <- folds(cc$uniqueclust,splits)
        res <- lapply(all.folds, function(i) {
            idx <- na.omit(as.vector(cc$idclustmat[i,]+1))
            do.call(func, c(list(data=data[idx,]),...))
        })        
        return(res)
    }    
    splits <- round(nn/size)
    all.folds <- folds(nn,splits)
    res <- lapply(all.folds, function(i) 
	do.call(func, c(list(data=data[i,]),...)))
    res
} ## }}} 



