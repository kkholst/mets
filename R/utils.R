###{{{ RoundMat

RoundMat <- function(cc,digits = max(3, getOption("digits") - 2),na=TRUE,...) {
    res <- format(round(cc,max(1,digits)),digits=digits)
    if (na) return(res)
    res[grep("NA",res)] <- ""
    res
}

###}}} RoundMat

###{{{ multinomlogit

multinomlogit <- function(x,tr=exp,dtr=exp) {
  n <- length(x)
  ex <- tr(x)
  dex <- dtr(x)
  sx <- sum(ex)+1
  f <- c(ex,1)
  df <- c(dex,0)
  res <- f/sx
  dg <- -dex/sx^2   
  gradient <- matrix(ncol=n,nrow=n+1)
  I <- diag(n+1)
  for (i in seq_len(n)) {
    gradient[,i] <- df[i]*I[i,]/sx+dg[i]*f
  }
  attributes(res)$gradient <- gradient
  return(res)
}

###}}} multinomlogit

###{{{ grouptable

##' Create Group Contingency Table from Clustered Data
##'
##' Creates a contingency table by group from paired/clustered data,
##' optionally combining lower and upper triangles.
##'
##' @param data a data.frame.
##' @param id name of the cluster/pair identifier column.
##' @param group name of the grouping variable (e.g., zygosity).
##' @param var name of the outcome variable to tabulate.
##' @param lower logical; if TRUE, fold upper triangle into lower.
##' @param labels optional labels for levels of \code{var}.
##' @param order optional ordering of factor levels.
##' @param group.labels optional labels for the groups.
##' @param group.order optional ordering of groups.
##' @param combine separator for combining two groups (default \code{" & "}).
##' @param ... additional arguments.
##' @return A table or list of tables.
##' @export
grouptable <- function(data,id,group,var,lower=TRUE,
                       labels,order,
                       group.labels,group.order,
                       combine=" & ",...) {
    if (!missing(order) || !missing(labels)) {
        data[,var] <- as.factor(data[,var])
        if (missing(order)) order <- seq(length(labels))
        if (missing(labels)) labels <- levels(data[,var])
        data[,var] <- factor(data[,var],levels(data[,var])[order],labels=labels[order])
    }
    wide <- fast.reshape(data,id=id,varying=-group)    
    res <- lapply(split(wide,wide[,group]),
                  function(x) {
                      M <- with(x, table(get(paste(var,"1",sep="")),
                                         get(paste(var,"2",sep=""))))
                      if (lower) {
                          M[lower.tri(M)] <- M[lower.tri(M)]+M[upper.tri(M)]
                          M[upper.tri(M)] <- NA
                      }
                      return(M)
                  })
    if (!missing(group.order) && length(group.order)==length(res))
        res <- res[group.order]
    if (!missing(group.labels) && length(group.labels)==length(res))
        names(res) <- group.labels
    if (length(res)==2 && !is.null(combine)) {
        M <- res[[1]]
        M[upper.tri(M)] <- res[[2]][lower.tri(res[[2]])]
        diag(M) <- paste(diag(M),diag(res[[2]]),sep=combine)
        M <- cbind(rownames(M),M)
        M <- rbind(c("",rownames(M)),M)
        colnames(M) <- rownames(M) <- rep("",nrow(M))
        M[1,1] <- paste(names(res),collapse=combine)
        return(structure(M,class="table"))
        return(M);
    }    
    res
}

###}}} grouptable

## Pseudo-inverse
pinv <- function(x, ...) {
  lava::Inverse(x, ...)
}
