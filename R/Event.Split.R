##' Event split with two time-scales, time and gaptime 
##'
##' splits after cut times for the two time-scales. 
##'
##' @param data data to be split
##' @param time time variable.
##' @param status status variable.
##' @param entry name of entry variable.
##' @param cuts cuts variable or numeric cut (only one value)
##' @param name.id name of id variable.
##' @param gaptime gaptime variable.
##' @param gaptime.entry name of entry variable for gaptime.
##' @param cutttime to cut after time or gaptime
##' @param cens.code code for the censoring.
##' @param order.id order data after id and start.
##' @author Thomas Scheike
##' @keywords survival
##' @examples
##' rr  <- data.frame(time=c(500,1000),start=c(0,500),status=c(1,1),id=c(1,1))
##' rr$gaptime <-  rr$time-rr$start
##' rr$gapstart <- 0
##'
##' rr1 <- Event.Split(rr,cuts=600,cuttime="time",   gaptime="gaptime",gaptime.entry="gapstart")
##' rr2 <- Event.Split(rr1,cuts=100,cuttime="gaptime",gaptime="gaptime",gaptime.entry="gapstart")
##'
##' dlist(rr1,start-time+status+gapstart+gaptime~id)
##' dlist(rr2,start-time+status+gapstart+gaptime~id)
##'
##' @export 

Event.Split <- function(data,
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
    } else { new.cuts  <-  data[,cutname] }

    if (is.numeric(entry)) {
	    start0 <- entry
	    name.entry <- paste("start",entry,sep=".")
            data[,name.entry] <- start0
    }  else name.entry <- entry

    if ((name.entry %in% names(data))) {
      new.start <- data[,name.entry]
    } else new.start <- rep(0,n)
    nnstart <- new.start

    if (any(new.start>= new.time)) cat("any(new.start>= new.time) is TRUE\n"); 

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
        nord <- c(1:n,(1:n)[splits]+0.1)
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
