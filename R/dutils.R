##' Cutting, sorting, rm (removing), rename for data frames 
##' 
##' Cut variables, rm or rename 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param cuts vector of number of groups, 4 is default and gives quartiles.
##' @param breaks  possible breaks for cutting.
##' @param regex for regular expressions.
##' @param sep seperator for naming of cut names.
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike 
##' @examples
##' data("sTRACE",package="timereg")
##' sTRACE$age2 <- sTRACE$age^2
##' sTRACE$age3 <- sTRACE$age^3
##'
##' mm <- dcut(sTRACE,~age+wmi)
##' head(mm)
##'
##' mm <- dcut(sTRACE,~age+wmi,cuts=c(2,4))
##' head(mm)
##'
##' mm <- dcut(sTRACE,c("age","wmi"))
##' head(mm)
##'
##' mm <- dcut(sTRACE,~.)
##' head(mm)
#
##' mm <- dcut(sTRACE,c("age","wmi"),cuts=c(2,4))
##' head(mm)
##'
##' gx <- qcut(sTRACE$age)
##' head(gx)
##'
##'
##' ## Removes all cuts variables with these names wildcards
##' mm <- drm(mm,c("*.2","*.4"))
##' head(mm)
##'
##' ## wildcards, for age, age2, age4 and wmi 
##' head(dcut(mm,c("a*","?m*")))
##'
##' ## with direct asignment 
##' drm(mm) <- c("*.2","*.4")
##' head(mm)
##'
##' dcut(mm) <- c("age","*m*")
##' head(mm)
##'
##' ## renaming 
##' ############################
##'
##' head(mm)
##' drename(mm, ~age+wmi) <- c("Wmi","Age")
##' head(mm)
##'
##' ## A* to lower 
##' mm2 <-  drename(mm,c("A*","W*"))
##' head(mm2)
##' drename(mm,"A*") <- ~.
##' head(mm)
##' drename(mm) <- ~.
##' head(mm)
##' @aliases dcut dcut<-
##' @export
dcut <- function(data,x,cuts=4,breaks=NULL,regex=FALSE,sep=NULL,...)
{# {{{

 if (is.null(sep) & is.null(breaks))    sep <- "."
 if (is.null(sep) & (!is.null(breaks))) sep <- "b"

if (missing(x)) x<- ~.

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
     xnames <- x
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }


  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) {
			           if (length(x)==1 && is.character(x)) x <- data[,x]
			           x
	       })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)
  if (ll>1 & length(cuts)==1) cuts <- rep(cuts,ll)
  if (length(x)!=length(cuts)) stop("length of variables not consistent with cuts")


for (k in 1:ll) 
{
  xx <- x[[k]]
  if (is.numeric(xx)) {
  if (is.null(breaks)) {
     probs <- seq(0,1,length.out=cuts[k]+1)
     name<-paste(xnames[k],cuts[k],sep=sep)
     bb <- quantile(xx,probs)
     } else { bb <- breaks ; name<-paste(xnames[k],breaks[2],sep=sep) }

  if (sum(duplicated(bb))==0)
   data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE,...)
  }
}

return(data)
}# }}}

##' @export
"dcut<-" <- function(data,...,value)
{# {{{
 dcut(data,value,...)
}# }}}

##' @export
drm <- function(data,x,regex=FALSE)
{# {{{
 if (inherits(x,"formula")) {
     xnames <- all.vars(x)
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }

  data[,xnames] <- NULL

return(data)
}# }}}

##' @export
"drm<-" <- function(data,...,value)
{# {{{
 drm(data,value,...)
}# }}}

##' @export
drename <- function(data,var,value,regex=FALSE) 
{  # {{{

 if (missing(var)) var <- names(data)

 if (inherits(var,"formula")) {
      var <- all.vars(var)
   }

if (is.character(var)) {
        varxnames <- var 
        varnnames <- c()
        xxx<-c()
     for (xx in varxnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
        varnnames <- c(varnnames,n) 
     }
     varxnames <- xxx[!duplicated(xxx)]
     varnnames <- varnnames[!duplicated(varnnames)]
  }

 if (missing(value)) value <- tolower(varxnames)

 if (inherits(value,"formula")) {
     value <- all.vars(value)
     if (value[1]==".") value <- tolower(varxnames)
 }

   if (is.character(value)) {
        xnames <- value  
  }

### print(xnames); print(varxnames)

  if (length(xnames)!= length(varxnames)) stop("length of old and new variables must mach")

  names(data)[varnnames] <- xnames; 
  return(data) 
} # }}}

##' @export
"drename<-" <- function(x, var, value) 
{ # {{{
  drename(x,var,value)
}# }}}

##' @export
dkeep <- function(data,var,keep=TRUE,regex=FALSE) 
{  # {{{

 if (inherits(var,"formula")) { var <- all.vars(var) }

if (is.character(var)) {
        varxnames <- var 
        varnnames <- c()
        xxx<-c()
     for (xx in varxnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
        varnnames <- c(varnnames,n) 
     }
     varxnames <- xxx[!duplicated(xxx)]
     varnnames <- varnnames[!duplicated(varnnames)]
  }


if (keep) data <- data[,varnnames] else data <- data[,-1*varnnames]
  return(data) 
} # }}}

##' @export
"dkeep<-" <- function(x,...,value) 
{ # {{{
  dkeep(x,value,...)
}# }}}

##' @export
ddrop <- function(data,var,keep=FALSE) 
{  # {{{
	dkeep(data,var,keep=FALSE)
} # }}}

##' @export
"ddrop<-" <- function(x,...,value) 
{ # {{{
  dkeep(x,value,keep=FALSE,...)
}# }}}


##' Sort data according to columns in data frame
##'
##' @title Sort data frame
##' @param data Data frame
##' @param x variable to order by
##' @param ... additional variables to order by
##' @param decreasing sort order (vector of length x)
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
##' @export
dsort <- function(data,x,...,decreasing=FALSE) 
{# {{{
    if (missing(x)) return(data)
    if (inherits(x,"formula")) {
        xx <- lava::procformula(value=x)$res
        decreasing <- unlist(lapply(xx,function(x) substr(trim(x),1,1)=="-"))
        x <- all.vars(x)
    } 
    if (is.character(x) && length(x)<nrow(data)) x <- lapply(x,function(z) data[,z])
    dots <- list(...)
    args <- lapply(dots, function(x) {
        if (length(x)==1 && is.character(x)) x <- data[,x]
        x
    })
    if (!is.list(x)) x <- list(x)
    data[do.call("order",c(c(x,args),list(decreasing=decreasing,method="radix"))),]
}# }}}

##' @export
"dsort<-" <- function(data,...,value) 
{ # {{{
  dsort(data,value,...)
}# }}}


ddnames <- function(data,x,g,regex=FALSE)
{# {{{
if (missing(x)) x<- ~.

 if (inherits(x,"formula")) {
         avx <- all.vars(x)
         if (length(avx)==1 & avx[1]==".") x <- names(data)
	 else {
            rhs  <-  all.vars(update(x, 0~.))
	    lhs <-  all.vars(update(x, .~0))
	    if (lhs[1]!=".")              { x <- lhs;}
	    if (lhs[1]=="." & rhs[1]!=".") { x <- rhs;}
	    if (lhs[1]!="." & rhs[1]!=".") { x <- lhs;  g <- rhs;}
	 }
 } 


 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
     xnames <- x
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }


 group<-NULL
 if (!missing(g))
 if (inherits(g,"formula")) {
     g <- all.vars(g)
     if (g[1]==".") g <- names(data) 
     gnames <- g
     group <- data[,gnames]
  } else if  (is.character(g)) {
     gnames <- g
     xxx<-c()
     for (xx in gnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     gnames <- xxx[!duplicated(xxx)]
     group <- data[,gnames]
  }

 if (missing(g)) group <- gnames <- NULL 


 return(list(xnames=xnames,gnames=gnames,group=group))

}# }}}


##' summary, tables, and correlations for data frames 
##' 
##' summary, tables, and correlations for data frames 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param g possible group variable
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike 
##' @examples
##' data("sTRACE",package="timereg")
##' dt<- sTRACE
##' dt$time2 <- dt$time^2
##' dt$wmi2 <- dt$wmi^2
##' head(dt)
##' 
##' dcor(dt)
##' 
##' dcor(dt,~time+wmi)
##' dcor(dt,~time+wmi,~vf+chf)
##' dcor(dt,time+wmi~vf+chf)
##' 
##' dcor(dt,c("time*","wmi*"),~vf+chf)
##' 
##' @export
dcor <- function(data,x,g,regex=FALSE,...)
{# {{{

 ddname <- ddnames(data,x,g,regex=regex)
 xnames <- ddname$xnames 
 gnames <- ddname$gnames 
 group  <- ddname$group

 if (!is.null(group)) return(by(data[,xnames],group,cor,...))
 if (is.null(group))  return(cor(data[,xnames],...))

}# }}}

##' @export
dsummary <- function(data,x,g,regex=FALSE,...)
{# {{{
### x+y ~ group1+group2 to get correlations against groups defined by group1*group2

ddname <- ddnames(data,x,g,regex=regex)
xnames <- ddname$xnames 
gnames <- ddname$gnames 
group  <- ddname$group

 if (!is.null(group)) return(by(data[,xnames],group,summary,...))
 if (is.null(group))  return(summary(data[,xnames],...))

}# }}}


##' tables for data frames 
##' 
##' tables for data frames 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param g possible group variable
##' @param level  1 for all marginal tables, 2 for all 2 by 2 tables, and null for the full table, possible versus group variable
##' @param regex for regex wildcards otherwise unix style wildcards 
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike 
##' @examples
##' data("sTRACE",package="timereg")
##' dt<- sTRACE
##'
##' dtable(dt,~status)
##'
##' dtable(dt,~status+vf)
##' dtable(dt,~status+vf,level=1)
##'
##' dtable(dt,~status+vf,~chf+diabetes)
##'
##' dtable(dt,status+vf~chf+diabetes)
##' 
##' dtable(dt,c("*f*","status"),~diabetes)
##' dtable(dt,c("*f*","status"),~diabetes,level=NULL)
##' 
##' dtable(dt,c("*f*","status"),level=1)
##' @export
dtable<- function(data,x,g,level=2,regex=FALSE,...)
{# {{{

ddname <- ddnames(data,x,g,regex=regex)
xnames <- ddname$xnames 
gnames <- ddname$gnames 
group  <- ddname$group


 if (!is.null(level)) {
 ### all tables from xnames over g 
 ll<-list()
 k<-1
 nn <- length(xnames)

 if (level==1) {
	 for (i in seq(1,(nn)))  {
		x1<-xnames[i]
		if (!is.null(group)) llk<-by(data[,c(x1)],group,table,...) 
		else llk<-table(data[,x1],...)
		ll[[k]]<- list(name=x1,table=llk)
		k<-k+1
	 }
 } else {
 ### all 2 by 2 tables from xnames over g 
	if (nn>1) {
	 for (i in seq(1,(nn-1)))  {
	 for (j in seq((i+1),nn)) {
		 x1<-xnames[i]
		 x2<-xnames[j]
		if (!is.null(group)) llk<-by(data[,c(x1,x2)],group,table,...) 
		else llk<-table(data[,x1],data[,x2],...)
		ll[[k]]<- list(name=paste(x1,"x",x2,sep=""),table=llk)
		k<-k+1
	 }
	 }
	} else {
		if (!is.null(group)) llk<-by(data[,xnames],group,table,...) 
		else llk<-table(data[,xnames],...)
		ll[[k]]<- list(name=xnames,table=llk)
		k<-k+1
	}

 }
 } else { ## one big table 
     if (!is.null(group)) ll<-by(data[,xnames],group,table,...) 
     else ll<- table(data[,xnames],...)  
 }

 return(ll)

}# }}}

###dtable<- function(data,x,g,level=2,regex=FALSE,...)
###{# {{{
###
###ddname <- ddnames(data,x,g,regex=regex)
###xnames <- ddname$xnames 
###gnames <- ddname$gnames 
###group  <- ddname$group
###
###
### if (!is.null(level)) {
### ### all tables from xnames over g 
### ll<-list()
### k<-1
### nn <- length(xnames)
###
### if (level==1) {
###	 for (i in seq(1,(nn)))  {
###		x1<-xnames[i]
###		if (!is.null(group)) llk<-by(data[,c(x1,x2)],group,table,...) 
###		else llk<-table(data[,x1],data[,x2],...)
###		ll[[k]]<- list(name=x1,table=llk)
###		k<-k+1
###	 }
### } else {
### ### all 2 by 2 tables from xnames over g 
###	if (length(nn)>1) {
###	 for (i in seq(1,(nn-1)))  {
###	 for (j in seq((i+1),nn)) {
###		 x1<-xnames[i]
###		 x2<-xnames[j]
###		if (!is.null(group)) llk<-by(data[,c(x1,x2)],group,table,...) 
###		else llk<-table(data[,x1],data[,x2],...)
###		ll[[k]]<- list(name=paste(x1,"x",x2,sep=""),table=llk)
###		k<-k+1
###	 }
###	 }
###	} else {
###		if (!is.null(group)) llk<-by(data[,xnames],group,table,...) 
###		else llk<-table(data[,xnames],...)
###		ll[[k]]<- list(name=paste(x1),table=llk)
###		k<-k+1
###	}
###
######  if (!is.null(group)) llk<-by(data[,xnames],group,table,...) 
######	else llk<-table(data[,xnames],...)
######        ll[[k]]<- list(name=xnames,table=llk)
######	k<-k+1
###
### }
### } else { ## one big table 
###     if (!is.null(group)) ll<-by(data[,xnames],group,table,...) 
###     else ll<- table(data[,xnames],...)  
### }
###
### return(ll)
###
###}# }}}

##' @export
dhead <- function(data,x,regex=FALSE,...)
{# {{{

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
     xnames <- x
     formular <- 1
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }

return(head(data[,xnames],...))

}# }}}

##' @export
dreshape <- function(data,...)
{# {{{
    fast.reshape(data,...)
}# }}}



procform <- function(formula, sep, ...) {
    aa <- attributes(terms(formula,...))
    if (aa$response == 0) {
        res <- NULL
    } else {
        res <- paste(deparse(formula[[2]]), collapse = "")
    }
    if (!missing(sep) && length(aa$term.labels) > 0) {
        attributes(res)$x <- lapply(strsplit(aa$term.labels, 
           "\\|")[[1]], function(x) as.formula(paste0("~", x)))
        
    } else {
        attributes(res)$x <- aa$term.labels
    }
    return(res)
}

##' aggregating for for data frames 
##' 
##' aggregating for for data frames 
##' @examples
##' data("sTRACE",package="timereg")
##' daggregate(iris, "^.e.al", group="Species", fun=cor, regex=TRUE)
##' daggregate(iris, ~Sepal.Length+Petal.Length|Species, fun=summary)
##' daggregate(iris, ~log(Sepal.Length)+I(Petal.Length>1.5)|Species, fun=summary)
##' daggregate(iris, "*Length*", group="Species", fun=head)
##' daggregate(iris, "^.e.al", group="Species", fun=tail, regex=TRUE)
##' daggregate(sTRACE, status~ diabetes, fun=table)
##' daggregate(sTRACE, status~ diabetes+sex, fun=table)
##' daggregate(sTRACE, status~ diabetes+sex|vf+I(wmi>1.4), fun=table)
##' daggregate(iris, "^.e.al", group="Species",regex=TRUE)
##' daggregate(iris, I(Sepal.Length>7)~Species | I(Petal.Length>1.5))
##' daggregate(iris, I(Sepal.Length>7)~Species | I(Petal.Length>1.5), fun=table)
##' daggregate(iris, I(Sepal.Length>7)~Species | I(Petal.Length>1.5), fun=lava:::Print)
##' @export
daggregate <- function(data,x,y,...,group=NULL,fun="summary",regex=FALSE) 
{# {{{
    if (missing(x)) x <- colnames(data)
    if (inherits(x,"formula")) {
        xf <- procform(x,sep="|",data=data)
        xx <- attr(xf,"x")
        if (length(xx)>1) {
            if (length(xx)>0) {
                if (all.vars(xx[[1]])[1]==".") {
                    x <- as.formula(paste0("~.","-",paste0(all.vars(xx[[2]]),collapse="-")))
                } else {
                    x <- update(as.formula(paste0(xf,"~.")),xx[[1]])
                }
            }
            group <- model.frame(xx[[2]],data=data)
        } else {            
        }
                                        #
        if (length(xx)>1) {
            x <- subset(model.frame(x,data=data),select=-xx[[2]])
        } else {
            x <- model.frame(x,data=data)
        }
    } else {
        xx <- c()
        for (x0 in x) {
            if (!regex) x0 <- glob2rx(x0)
            n <- grep(x0,names(data),)
            xx <- union(xx,names(data)[n])
        }
        x <- data[,xx,drop=FALSE]
    }
    if (inherits(group,"character") && length(group)<NROW(data)) {
        group <- data[,group]
    }
    if (inherits(group,"formula")) {
        group <- model.frame(group,data=data)
    }
    ##browser()
    if (is.character(fun)) fun <- get(fun)
    if (!is.null(group)) {        
        return(by(x,group,fun,...))
    }
    
    do.call(fun, c(list(x),list(...)))
}# }}}


##' Lag operator
##'
##' Lag operator
##' @examples
##' d <- data.frame(y=1:10,x=c(10:1))
##' dlag(d,k=1:2)
##' dlag(d,~x,k=0:1)
##' dlag(d$x,k=1)
##' dlag(d$x,k=-1:2)
##' @export
##' @param data data.frame or vector
##' @param x optional column names or formula
##' @param k lag (vector of integers)
##' @param combine combine results with original data.frame
##' @param simplify Return vector if possible
##' @param names optional new column names
##' @param ... additional arguments to lower level functions
dlag <- function(data,x,k=1,combine=TRUE,simplify=TRUE,names,...) {
    isvec <- FALSE
    if (!is.data.frame(data)) {
        isvec <- is.vector(data)
        data <- as.data.frame(data)
        combine <- FALSE
    }
    if (missing(x)) x <- base::names(data)
    if (is.character(x)) x <- data[,x,drop=FALSE]
    if (inherits(x,"formula")) {
        x <- as.data.frame(model.matrix(update(x,~.-1),data=data,...))
    }
    kmin0 <- abs(min(k,0))
    kmax0 <- max(k,0)
    kmax <- kmin0+kmax0
    pad <- function(x) c(rep(NA,kmax0),x,rep(NA,kmin0))
    kidx <- match(k,seq(min(k),max(k)))
    val <- lapply(x,function(y) embed(pad(y),dimension=kmax+1)[,kidx,drop=FALSE])
    dval <- as.data.frame(val)
    if (!missing(names)) base::names(dval) <- names
    if (combine) {
        return(cbind(dval,data))
    }
    if (length(k)==1 && simplify && isvec) return(as.vector(val[[1]]))
    return(dval)
}