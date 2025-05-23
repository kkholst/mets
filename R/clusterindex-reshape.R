##' Finds subjects related to same cluster
##' 
##' @references
##' Cluster indeces 
##' @examples
##' i<-c(1,1,2,2,1,3)
##' d<- cluster.index(i)
##' print(d)
##'
##' type<-c("m","f","m","c","c","c")
##' d<- cluster.index(i,num=type,Rindex=1)
##' print(d)
##' @seealso familycluster.index familyclusterWithProbands.index
##' @author Klaus Holst, Thomas Scheike
##' @param clusters  list of indeces
##' @param index.type if TRUE then already list of integers of index.type
##' @param num to get numbering according to num-type in separate columns
##' @param Rindex index starts with 1, in C is it is 0
##' @param mat to return matrix of indeces
##' @param return.all return all arguments
##' @param code.na how to code missing values
##' @aliases countID pairRisk mystrata mystrata2index
##' @export
cluster.index <- function(clusters,index.type=FALSE,num=NULL,Rindex=0,mat=NULL,return.all=FALSE,code.na=NA)
{ ## {{{
  n <- length(clusters)

  if (index.type==FALSE)  {
    if (is.numeric(clusters)) clusters <-  fast.approx(unique(clusters),clusters)-1 else  {
      max.clust <- length(unique(clusters))
      clusters <- as.integer(factor(clusters, labels = seq(max.clust)))-1
    }
  }
  
  if ((!is.null(num))) { ### different types in different columns
    mednum <- 1
    if (is.numeric(num)) numnum <-  fast.approx(unique(num),num)-1
    else {
      numnum <- as.integer(factor(num, labels = seq(length(unique(clusters))))) -1
    }
  } else { numnum <- 0; mednum <- 0; }

  clustud <- .Call("clusterindexM",as.integer(clusters),as.integer(mednum), as.integer(numnum),mat,return.all,PACKAGE="mets")
  if (!is.null(mat) && !return.all) return(clustud)
  
  if (Rindex==1) clustud$idclust <- clustud$idclustmat+1
  if (Rindex==1) clustud$firstclustid <- clustud$firstclustid +1 
  ### avoid NA's for C call
  if (Rindex==0 & !is.na(code.na)) clustud$idclust[is.na(clustud$idclust)] <- code.na
  
  clustud
} ## }}}

##' @export
countID <- function(data,id="id",names.count="Count",total.count="Total",index.name="index",reverse=TRUE,sep="",addid=TRUE,sorted=FALSE)
{# {{{

clusters <- data[,id]
if (is.numeric(clusters)) {
   ## integeres from 0 to max.clust
   uc <- unique(clusters)
   max.clust <- length(uc)
   clusters <- fast.approx(uc, clusters) - 1
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

if (sorted) {
 clusters <- (0:(max.clust-1))[order(uc)][clusters+1]
}

###tabid <- table(clusters)
nclust <- table(clusters)

if (addid)  {
name1 <- paste(names.count,id,sep=sep)
name2 <- paste("index",id,sep=sep)
name3 <- paste(total.count,id,sep=sep)
name4 <- paste("reverse",names.count,id,sep=sep)
} else {
name1 <- names.count
name2 <- index.name
name3 <- total.count
name4 <- paste("reverse",names.count,sep=sep)
}

out <- data[,id,drop=FALSE]
out[,name1] <- c(cumsumstrata(rep(1,nrow(data)),clusters,max.clust))
out[,name2] <- clusters 
out[,name3] <- nclust[clusters+1]
out[,name4] <- 1+out[,name3]-out[,name1]

attr(out,"max.clust") <- max.clust

return(out)
}# }}}

##' @export
pairRisk <- function(start,stop,status,expo,clust,nsize=10,doublerisk=1)
{# {{{
    n <- length(start)
    id <- 1:n
    nclust <- length(unique(clust))
    sig <- c(rep(-1, each = n), rep(1, each = n))
    clust.seq <- c(as.numeric(as.factor(clust)), as.numeric(as.factor(clust)))
    clust <- c(clust, clust)
    expo <- c(expo, expo)
    id <- c(id, id)
    sstatus <- c(rep(0, length(start)), status)
    tts <- c(start, stop)
    ot <- order(clust.seq, tts, -rank(sstatus))
    tts <- tts[ot]
    sstatus <- sstatus[ot]
    sig <- sig[ot]
    expo <- expo[ot]
    clust <- clust[ot]
    clust.seq <- clust.seq[ot]
    id <- id[ot]
    cc <- c(mets::revcumsumstrata(sig * expo * nsize + sig * (expo == 0), clust.seq - 1, nclust))
    ### both under risk when cc>10
    if (doublerisk)
    pair.risk <- which(cc > nsize)
    else pair.risk <- which(cc>=1)
    clustpl <- clust[pair.risk]
    weightpl <- cc[pair.risk] - nsize
    caseweightpl <- rep(-1, length(weightpl))
    casepl <- expo[pair.risk]
    caseweightpl[casepl == 1] <- weightpl[casepl == 1]
    ttexit <- tts[pair.risk]
    ttentry <- tts[pair.risk - 1]
    ttid <- id[pair.risk]
    tstatus <- sstatus[pair.risk]
    timesout <- cbind(rep(ttentry, times = weightpl), rep(ttexit, times = weightpl))
    whichnotsame <- which(timesout[, 1] != timesout[, 2])
    caseweightrep <- rep(caseweightpl, times = weightpl) * (!duplicated(cbind(rep(ttexit,
        times = weightpl), rep(clustpl, time = weightpl)))) * 1
    weightstatusrep <- rep(tstatus, times = weightpl) * (caseweightrep != 0)
    idout <- rep(ttid, times = weightpl) * (caseweightrep != 0)
    clustplrep <- rep(clustpl, times = weightpl)
    out <- data.frame(timesout, weightstatusrep, caseweightrep,
        clustplrep, idout, stringsAsFactors = FALSE)[whichnotsame, ]
    out <- out[order(out[, 5]), ]
    return(out)
}# }}}

##' @export
mystrata <- function(ll,sort=TRUE) {# {{{
	for (j in seq(1,length(ll))) 
	if (!is.factor(ll[[j]])) ll[[j]] <- factor(ll[[j]])

	nll <- length(ll[[1]])
        ss <- rep(0,nll)
        nl <- unlist(lapply(ll,nlevels))
        poss <- c(exp(revcumsum(log(nl[-1]))))
	poss <- c(poss,1)
	for (j in seq(1,length(ll))) {
	     ss <- ss+as.numeric(ll[[j]])*poss[j]
	}
	uss <- unique(ss)
	nindex <- length(uss)
        sindex <- fast.approx(uss,ss)
        attr(sindex,"nlevel") <- nindex
        attr(sindex,"levels") <- nl 
	return(sindex)
} # }}}

##' @export
mystrata2index <- function(ll,sort=TRUE) {# {{{
	nll <- nrow(ll)
        ss <- rep(0,nll)
        nl <- apply(ll,2,max)
        poss <- c(exp(revcumsum(log(nl[-1]))))
	poss <- c(poss,1)
	for (j in seq(1,ncol(ll))) { ss <- ss+ll[,j]*poss[j] }
	uss <- unique(ss)
	nindex <- length(uss)
        sindex <- fast.approx(uss,ss)
        attr(sindex,"nlevel") <- nindex
        attr(sindex,"levels") <- nl 
	return(sindex)
} # }}}


##' Finds all pairs within a cluster (family)
##' 
##' @references
##' Cluster indeces 
##' @examples
##' i<-c(1,1,2,2,1,3)
##' d<- familycluster.index(i)
##' print(d)
##' @seealso cluster.index familyclusterWithProbands.index
##' @author Klaus Holst, Thomas Scheike
##' @param clusters  list of indeces 
##' @param index.type argument of cluster index 
##' @param num num 
##' @param Rindex index starts with 1 in R, and 0 in C
##' @export
familycluster.index <- function(clusters,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
  clusters <- cluster.index(clusters,Rindex=Rindex)
  totpairs <- sum(clusters$cluster.size*(clusters$cluster.size-1)/2)
  clustud <- .Call("familypairindex",clusters$idclust,clusters$cluster.size,as.integer(2*totpairs),PACKAGE="mets")
  clustud$pairs <- matrix(clustud$familypairindex,ncol=2,byrow=TRUE)
  clustud$clusters <- clusters$clusters[clustud$pairs[,2]]+1

  invisible(clustud)
} ## }}}

##' Finds all pairs within a cluster (famly)  with the proband (case/control) 
##' 
##' second column of pairs are the probands and the first column the related subjects
##' 
##' @references
##' Cluster indeces 
##' @examples
##' i<-c(1,1,2,2,1,3)
##' p<-c(1,0,0,1,0,1)
##' d<- familyclusterWithProbands.index(i,p)
##' print(d)
##' @author Klaus Holst, Thomas Scheike
##' @seealso familycluster.index cluster.index
##' @param clusters list of indeces giving the clusters (families)
##' @param probands list of 0,1 where 1 specifices which of the subjects that are probands 
##' @param index.type argument passed to other functions
##' @param num argument passed to other functions
##' @param Rindex index starts with 1, in C is it is 0
##' @export
familyclusterWithProbands.index <- function(clusters,probands,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
    famc <-familycluster.index(clusters,index.type=index.type,num=num,Rindex=Rindex)
    if (length(probands)!=length(clusters)) stop("clusters and probands not same length\n"); 
    index.probs <- (1:length(clusters))[probands==1]
    subfamsWprobands <-famc$subfamilyindex[ famc$familypairindex %in% index.probs ]
    indexWproband <- famc$subfamilyindex %in% subfamsWprobands 
    famc$subfamilyindex <- famc$subfamilyindex[indexWproband]
    famc$familypairindex <- famc$familypairindex[indexWproband]
    pairs <- matrix(famc$familypairindex,ncol=2,byrow=TRUE)
    ipi1 <- pairs[,1] %in% index.probs
    gem2 <- pairs[,2]
    pairs[ipi1,2] <- pairs[ipi1,1]
    pairs[ipi1,1] <- gem2[ipi1]
    famc$pairs <- pairs
    famc$clusters <-  famc$clusters[ipi1]

    famc$familypairindex <- c(t(pairs))

    invisible(famc)
} ## }}}

##' @export
coarse.clust <- function(clusters,max.clust=100)
{ ## {{{ 

if (is.numeric(clusters)) 
   clusters <-  fast.approx(unique(clusters),clusters)

  cluster.size <- length(unique(clusters))

qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
cclusters <-  as.integer(qqc)-1

return(cclusters)
} ## }}} 

##' @export
faster.reshape <- function(data,clusters,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
  if (NCOL(data)==1) data <- cbind(data)
 ### uses data.matrix 
  if (!is.matrix(data)) data <- data.matrix(data)
  if (is.character(clusters)) clusters <- data[,clusters]
  n <- length(clusters)

  if (nrow(data)!=n)  stop("nrow(data) and clusters of different lengths\n"); 

  if (index.type==FALSE)  {
    max.clust <- length(unique(clusters))
    if (is.numeric(clusters)) clusters <-  fast.approx(unique(clusters),clusters)-1 else 
    {
      max.clust <- length(unique(clusters))
      clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1
    }
  }

  if ((!is.null(num))) { ### different types in different columns
    if (length(num)!=n)  stop("clusters and num of different lengths\n"); 
    mednum <- 1
    if (is.numeric(num)) num <-  fast.approx(unique(num),num)-1
    else num <- as.integer(factor(num, labels = seq(length(unique(clusters))))) -1
  } else { num <- 0; mednum <- 0; }

  clustud <- .Call("clusterindexdata",as.integer(clusters),as.integer(mednum), as.integer(num),iddata=data,PACKAGE="mets")

  if (Rindex==1) clustud$idclust  <- clustud$idclust+1
###  if(Rindex==1) idclust[idclust==0] <- NA 
  maxclust <- clustud$maxclust

  xny <- clustud$iddata
  xnames <- colnames(data); 
  missingname <- (colnames(data)=="")
  xnames[missingname] <- paste(seq_len(maxclust))[missingname]
  xny <- data.frame(xny)
  mm <- as.vector(outer(xnames,seq_len(maxclust),function(...) paste(...,sep=".")))
  names(xny) <- mm

  return(xny); 
} ## }}}
