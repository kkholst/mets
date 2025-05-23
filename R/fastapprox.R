##' Fast approximation
##' 
##' Fast approximation
##' @param time Original ordered time points
##' @param new.time New time points
##' @param equal If TRUE a list is returned with additional element
##' @param type Type of matching, nearest index, nearest greater than
##'     or equal (right), number of elements smaller than y otherwise
##'     the closest value above new.time is returned.
##' @param sorted Set to true if new.time is already sorted
##' @param ... Optional additional arguments
##' @author Klaus K. Holst
##' @examples
##' id <- c(1,1,2,2,7,7,10,10)
##' fast.approx(unique(id),id)
##' 
##' t <- 0:6
##' n <- c(-1,0,0.1,0.9,1,1.1,1.2,6,6.5)
##' fast.approx(t,n,type="left")
##' @aliases  indexstrata predictCumhaz cpred
##' @export
fast.approx <- function(time,new.time,equal=FALSE,type=c("nearest","right","left"),sorted=FALSE,...) {# {{{
    if (!sorted) {
        ord <- order(new.time,decreasing=FALSE)
        new.time <- new.time[ord]
    }
    A <- NULL
    if (NCOL(time)>1) {
        A <- time
        time <- A[,1,drop=TRUE]
    }
    if (is.unsorted(time)) warnings("'time' will be sorted")
    type <- agrep(type[1],c("nearest","right","left"))-1
    res <- .Call("FastApprox",
                 time=sort(time),
                 newtime=new.time,
                 equal=equal,
                 type=type,
                 PACKAGE="mets")
    if (!sorted) {
        oord <- order(ord)
        if (!equal) return(res[oord])
        res <- lapply(res,function(x) x[oord])
    }
    if (!is.null(A)) {
        A[res,,drop=FALSE]
    }
    return(res)
}# }}}

##' @export
predictCumhaz <- function (cum, new.time, type = c("left", "right", "nearest"),tminus=FALSE,tplus=FALSE,return.index=TRUE,...)
{# {{{
   if (NCOL(cum)>1) { cumh <- cum[,-1,drop=FALSE]; time <- cum[,1]} else time <- cum
   equal <- FALSE
   if (type[1]=="left" & tminus) equal <- TRUE
   if (type[1]=="right" & tplus) equal <- TRUE
   index <- fast.approx(time,new.time,type=type[1],equal=equal,...)
   if (equal) {
   equali <- which(index$eq!=0)
   index <- index$idx
   if (length(equali)>=1 & tminus & type[1]=="left") {
	   index[equali] <- index[equali] -1
   }
   if (length(equali)>=1 & tplus & type[1]=="right") {
	   index[equali] <- pmin(length(time),index[equali] +1)
   }
   if (any(index==0))  { 
      index[index==0] <- 1 
   }
   }
   if (NCOL(cum)>1) {
	   res <- cbind(new.time, cumh[index,])
   } else {
	   if (return.index) res <-  index else res <- time[index]
   }
   return(res)
}# }}}

##' @export
cpred <- function(...) predictCumhaz(...)

##' @export
indexstrata <- function(jump.times,eval.times,jump.strata=NULL,eval.strata=NULL,nstrata=NULL,start.time=NULL,...)
{   # {{{
	stopifnot(is.numeric(jump.times))
	stopifnot(is.numeric(eval.times))
	if (is.null(jump.strata)) jump.strata <- rep(0,length(jump.times))
	if (is.null(eval.strata)) eval.strata <- rep(0,length(eval.times))
###	if (any(eval.strata<0) | any(eval.strata>nstrata-1)) stop("eval.strata index not ok\n"); 
###	if (any(jump.strata<0) | any(jump.strata>nstrata-1)) stop("jump.strata index not ok\n"); 
	N <- length(jump.times)
	index <- rep(0,length(eval.times))
	
	for (i in unique(eval.strata)) {
	   wherej <- which(eval.strata==i)
	   whereJ <- which(jump.strata==i)
	   if (!is.null(start.time)) iindex <- predictCumhaz(c(start.time,jump.times[whereJ]),eval.times[wherej],...)
	   else iindex <- predictCumhaz(jump.times[whereJ],eval.times[wherej],...)
	iindexn0 <- which(iindex!=0)
        if (length(iindexn0)>0)
	   index[wherej[iindexn0]] <- whereJ[iindex[iindexn0]]
	}
	return(index)
}# }}}

