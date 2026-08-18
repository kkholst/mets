##' Split event-history records at one or more cut points (SurvSplit) 
##'
##' Splits each row's \code{(start, time]} interval at every cut point that
##' falls strictly inside it, replicating covariates onto the new rows and
##' setting the status of every intermediate (non-final) segment to
##' \code{cens.code}. Can also use a subject specific cut-point (event-based) 
##' when a variable name is given as the cuts.
##'
##' \code{cuts} can be:
##' \itemize{
##'   \item a single number: one global cut point applied to every row.
##'   \item a numeric vector: several global cut points, all applied in one
##'         call (equivalent to, but faster than, calling
##'         \code{event_splitC} once per cut point).
##'   \item a single column name (character): a per-row cut point, read from
##'         that column of \code{data} (so different rows can be cut at
##'         different times).
##'   \item a character vector of column names: several per-row cut points at
##'         once, one set of values per named column.
##' }
##' \code{NA} in a per-row cut column simply means "no cut" for that row/column
##' combination.
##'
##' @param data a data.frame with (at least) the time, status, and optionally
##'   start columns.
##' @param cuts cut point(s); see Details.
##' @param time name of the time (interval end) column.
##' @param status name of the status/event column.
##' @param name.id name of the subject-id column; created (as \code{1:n}) if
##'   not already present.
##' @param name.start name of the interval-start column. If numeric (e.g.
##'   \code{0}), that constant is used as the start time for every row and
##'   materialised into a new column called \code{start.<value>}.
##' @param cens.code status value assigned to every intermediate split
##'   segment (i.e. every segment except the last one for a given original
##'   row).
##' @param order.id if \code{TRUE}, sort the result by id then start time.
##' @param time.group if \code{TRUE}, add one \code{before.<cut>} 0/1
##'   indicator column per cut point/column, flagging whether the segment's
##'   start is before that particular cut.
##' @param group.var if not \code{NULL}, the name of a single integer column
##'   to add giving which time-axis bin the segment falls in: bin 1 is
##'   \code{[0, c1)}, bin 2 is \code{[c1, c2)}, ..., using every distinct cut
##'   value seen in \code{cuts} (sorted). Only meaningful when the data has
##'   actually been split at those same cut points.
##'
##' @return A data.frame with one row per split segment.
##'
##' @examples
##' d <- data.frame(event = round(5 * runif(5), 2),
##'                  start = 1:5, time = 2 * (1:5),
##'                  status = rbinom(5, 1, 0.5), x = 1:5)
##'
##' ## 1. single global cut point
##' event_split(d, cuts = 3.5)
##'
##' ## 2. several global cut points in one call
##' event_split(d, cuts = c(2, 4, 7))
##'
##' ## 3. per-row cut point taken from an existing column, with an
##' ##    explicit constant start time (start.0 is created automatically)
##' event_split(d, cuts = "event", name.start = 0)
##'
##' ## 4. several per-row cut points from multiple columns at once
##' d2 <- d
##' d2$cutA <- d2$event
##' d2$cutB <- d2$event + 1
##' event_split(d2, cuts = c("cutA", "cutB"), name.start = 0)
##' event_split(d2, cuts = c("cutA", "cutB"), name.start = "start")
##'
##' ## 5. group.var: a single time-axis bin label instead of several
##' ##    before/after indicators
##' event_split(d, cuts = c(2, 4, 7), group.var = "group")
##'
##' @export
event_split <- function(data, cuts,
                         time="time", status="status",
                         name.id="id", name.start="start",
                         cens.code=0, order.id=TRUE, time.group=FALSE,
                         group.var=NULL)
{ ## {{{
    n <- nrow(data)
    new.time   <- data[[time]]
    new.status <- data[[status]]

    ## --- resolve `cuts` into an n x m matrix of candidate cut values ---
    ## (built in-memory only; nothing is added to `data` here)
    if (is.character(cuts)) {
        ## cuts is (are) existing column name(s): value can differ per row
        missing.cols <- setdiff(cuts, names(data))
        if (length(missing.cols))
            stop("cuts refers to column(s) not in data: ",
                 paste(missing.cols, collapse=", "))
        cutnames <- cuts
        cutmat <- as.matrix(data[, cutnames, drop=FALSE])
    } else if (is.numeric(cuts)) {
        ## cuts is one or more fixed values, same for every row: no need to
        ## store these as columns, just broadcast them into the matrix
        cutnames <- paste("cut", cuts, sep=".")
        cutmat <- matrix(rep(cuts, each=n), nrow=n, ncol=length(cuts))
    } else {
        stop("cuts must be numeric, or a character vector of column names")
    }
    storage.mode(cutmat) <- "double"

    if (is.numeric(name.start)) {
        ## name.start is a constant value (e.g. 0), not a column name:
        ## materialize it as a real column, same as the original event_split
        start0 <- name.start
        name.start <- paste("start", name.start, sep=".")
        data[[name.start]] <- start0
    }

    if (name.start %in% names(data)) {
        new.start <- data[[name.start]]
    } else {
        new.start <- rep(0, n)
        data[[name.start]] <- new.start
    }
    if (any(new.start >= new.time)) cat("any(new.start>= new.time) is TRUE\n")

    if (name.id %in% names(data)) {
        idl <- data[[name.id]]
    } else {
        idl <- seq_len(n)
        data[[name.id]] <- idl
    }

    ## --- do the splitting in C++ ---
    res <- event_split_cpp(new.start, new.time, new.status, cutmat, cens.code)
    res <- as.data.frame(res)

    ## re-attach every other column (covariates, id, the cut column(s), ...)
    out <- data[res$row, , drop=FALSE]
    out[[time]]       <- res$time
    out[[status]]     <- res$status
    out[[name.start]] <- res$start
    rownames(out) <- NULL

    if (time.group) {
        if (is.character(cuts)) {
            for (cn in cutnames) {
                out[[paste("before", cn, sep=".")]] <- 1 * (out[[name.start]] < out[[cn]])
            }
        } else {
            for (j in seq_along(cuts)) {
                out[[paste("before", cutnames[j], sep=".")]] <- 1 * (out[[name.start]] < cuts[j])
            }
        }
    }

    if (!is.null(group.var)) {
        ## label each split segment with which time-axis bin it falls in:
        ## bin 1 = [0, c1), bin 2 = [c1, c2), ..., bin m+1 = [cm, Inf)
        ## uses every distinct cut value that was ever in play, whether the
        ## cuts came from fixed numbers or from (possibly per-row) columns
        sorted.cuts <- sort(unique(as.vector(cutmat)))
        sorted.cuts <- sorted.cuts[!is.na(sorted.cuts)]
        out[[group.var]] <- findInterval(out[[name.start]], sorted.cuts) + 1L
    }

    if (order.id) {
        out <- out[order(out[[name.id]], out[[name.start]]), ]
        rownames(out) <- NULL
    }
    out
} ## }}}


###event_split <- function(data,
###		time="time",status="status",cuts="cuts",name.id="id",
###		name.start="start", cens.code=0,order.id=TRUE, time.group=FALSE)
###{ ## {{{
###    n <- nrow(data)
###    new.time <- data[,time]
###    new.status <- data[,status]
###
###    if (is.numeric(cuts)) {
###        cutname <- paste("cut",cuts,sep=".")
###        data[,cutname] <- cuts
###    } else cutname <- cuts
###    new.cuts <- data[,cutname]
###
###    if (is.numeric(name.start)) {
###	    start0 <- name.start
###	    name.start <- paste("start",name.start,sep=".")
###            data[,name.start] <- start0
###    }  
###
###
###    if ((name.start %in% names(data))) {
###      new.start <- data[,name.start]
###    } else new.start <- rep(0,n)
###
###    if (any(new.start>= new.time)) cat("any(new.start>= new.time) is TRUE\n")
###
###    if ((name.id %in% names(data))) idl <- data[,name.id] else {
###	    idl <- 1:n
###	    data[,name.id] <- idl 
###    }
###
###    ## only split if cut not already among times
###    splits <- which(new.cuts<new.time & new.start<new.cuts)
###
###    if (length(splits)) {
###	    rows  <- c(1:n,splits)
###	    new.time <-   c(new.time,new.time[splits])
###	    new.start <-  c(new.start,new.cuts[splits])
###	    new.status <- c(new.status,new.status[splits])
###	    new.ccc <-    c(new.cuts,new.cuts[splits])
###	    idl <- c(idl,idl[splits])
###	    new.time[splits] <- new.cuts[splits]
###	    new.status[splits] <- cens.code
###	    data <- data[rows,]
###	    data[,time] <- new.time
###	    data[,status] <- new.status
###	    data[,name.start] <- new.start
###	    data[,name.id] <- idl
###    }
###
###    if (time.group) {
###      group.time <- paste("before",cutname,sep=".")
###      data[,group.time] <- 1*(data[,name.start]<data[,cutname]) ## sc(rep(1,n),rep(0,length(splits)))
###    } 
###
###    if (order.id) data <- data[order(idl,new.start),] 
###    rownames(data) <- NULL
###
###    return(data)
###}  ## }}}


##' Event split with two time-scales, time and gaptime 
##'
##' Cuts time for two time-scales, as event.split 
##'
##' @param data data to be split
##' @param time time variable.
##' @param status status variable.
##' @param entry name of entry variable.
##' @param cuts cuts variable or numeric cut (only one value)
##' @param name.id name of id variable.
##' @param gaptime gaptime variable.
##' @param gaptime.entry name of entry variable for gaptime.
##' @param cuttime to cut after time or gaptime
##' @param cens.code code for the censoring.
##' @param order.id order data after id and start.
##' @author Thomas Scheike
##' @keywords survival
##' @examples
##' rr  <- data.frame(time=c(500,1000),start=c(0,500),status=c(1,1),id=c(1,1))
##' rr$gaptime <-  rr$time-rr$start
##' rr$gapstart <- 0
##'
##' rr1 <- event_split2(rr,cuts=600,cuttime="time",   gaptime="gaptime",gaptime.entry="gapstart")
##' rr2 <- event_split2(rr1,cuts=100,cuttime="gaptime",gaptime="gaptime",gaptime.entry="gapstart")
##'
##' dlist(rr1,start-time+status+gapstart+gaptime~id)
##' dlist(rr2,start-time+status+gapstart+gaptime~id)
##'
##' @export 
event_split2 <- function(data,
		time="time",status="status",entry="start",cuts="cuts",name.id="id",
		gaptime=NULL,gaptime.entry=NULL,cuttime=c("time","gaptime"),
		cens.code=0,order.id=TRUE)
{
## {{{ 
    n <- nrow(data)
    new.time <- data[,time]
    new.status <- data[,status]

    if (!is.null(gaptime)) { new.gaptime <- data[,gaptime]; nngap <- new.gaptime;}
    if (!is.null(gaptime.entry)) new.gapstart <- data[,gaptime.entry]

    if (is.numeric(cuts)) {
	    new.cuts <- rep(cuts,nrow(data))
    } else { new.cuts  <-  data[,cuts] }

    if (is.numeric(entry)) {
	    start0 <- entry
	    name.entry <- paste("start",entry,sep=".")
            data[,name.entry] <- start0
    }  else name.entry <- entry

    if ((name.entry %in% names(data))) {
      new.start <- data[,name.entry]
    } else new.start <- rep(0,n)
    nnstart <- new.start

    if (any(new.start>= new.time)) cat("any(new.start>= new.time) is TRUE\n")

    if ((name.id %in% names(data))) idl <- data[,name.id] else {
	    idl <- 1:n
	    data[,name.id] <- idl 
    }

    if (cuttime[1]=="time") {# {{{
    splits <- which(new.cuts<new.time & new.start<new.cuts)

    if (length(splits)) {
	    rows  <- c(1:n,splits)
	    new.time <-   c(new.time,new.time[splits])
	    new.start <-  c(new.start,new.cuts[splits])
	    new.time[splits] <- new.cuts[splits]
	    new.status <- c(new.status,new.status[splits])
	    new.status[splits] <- cens.code
	    idl <- c(idl,idl[splits]) 

	    if (!is.null(gaptime)) {
	       new.gapstart <- c(new.gapstart,new.cuts[splits]-nnstart[splits]+new.gapstart[splits])
               new.gaptime <-   c(new.gaptime,new.gaptime[splits])
	       new.gaptime[splits] <- new.cuts[splits]- nnstart[splits]+new.gapstart[splits]
            }

	    data <- data[rows,]
	    data[,time] <- new.time
	    data[,status] <- new.status
	    data[,name.entry] <- new.start
	    data[,name.id] <- idl

	    if (!is.null(gaptime)) data[,gaptime]  <-  new.gaptime 
	    if (!is.null(gaptime)) data[,gaptime.entry] <- new.gapstart
    }

    if (order.id) data <- data[order(idl,new.start),] # }}}
    } else {# {{{
    splits <- which(new.cuts<new.gaptime & new.gapstart<new.cuts)

    if (length(splits)) {
	rows  <- c(1:n,splits)
	new.gaptime <-   c(new.gaptime,new.gaptime[splits])
        new.gapstart <-  c(new.gapstart,new.cuts[splits])
        new.gaptime[splits] <- new.cuts[splits]
        new.status <- c(new.status,new.status[splits])
        new.status[splits] <- cens.code
        idl <- c(idl,idl[splits]) 

        new.time  <- c(new.time,new.time[splits])
        new.time[splits] <- new.start[splits]+(new.gaptime[splits]-new.gapstart[splits])
        new.start  <-   c(new.start,new.start[splits]+(new.gaptime[splits]-new.gapstart[splits]))

        data <- data[rows,]
        data[,time] <- new.time
        data[,status] <- new.status
        data[,name.entry] <- new.start
        data[,name.id] <- idl
        ###
        data[,gaptime]  <-  new.gaptime 
        data[,gaptime.entry] <- new.gapstart
        if (order.id) data <- data[order(idl,new.start),] 
       }

    }# }}}

    return(data)
    ## }}} 
} 


