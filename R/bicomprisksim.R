
## n <- 1e4; ACE <- c(1,1,1)/3
## logscale <- -4.5; logshape <- .7
## p2 <- .065
## a2 <- -10; b2 <- 0.15 
## a1 <- 85; b1 <- 0.1

## pmvn(upper=c(q2,q2),sigma=diag(2)*(1-2/3)+2/3)
## pmvn(c(q2,q2),sigma=diag(2)*(1-2/3)+2/3)

## pnorm(q2,sd=1) ## Marginal / Perfect dependence
## pmvn(upper=c(q2,q2),sigma=diag(2)*(1-2/3)+2/3) ## Concordance
## pnorm(q2,sd=1)^2 ## Independence
## (lambdaR <- pmvn(upper=c(q2,q2),sigma=diag(2)*(1-2/3)+2/3)/pnorm(q2,sd=1)^2)

bicomprisksim <- function(n=1e4,
                          ACE=c(1/3,1/3,1/3),
                          logscale=-4.5,logshape=.7,
                          a1=85,b1=0.1,
                          a2=-10,b2=0.15,
                          p2=.065,
                          tt,
                          ...) {

    ACE <- ACE/sum(ACE)
    ialpha <- function(v,a,b) -(log(-v)+a)/b    
    q2 <- qnorm(p2) ## p2: Prostata cancer prevalence 
    p1 <- 1-p2; q1 <- qnorm(p1) ## Death without cancer
    alpha <- function(t,a,b) -exp(-(b*t+a))
    F2s <- function(t) pnorm(alpha(t,b=b2,a=a2)+q2) ## Marginal Cumulative Incidence (cancer)
### Random effects
    R <- diag(2)*.5+.5
    J <- matrix(1,ncol=2,nrow=2)
    I <- diag(2)
    zyg <- rep(c(0,1),each=n)
    A <- rbind(rmvn(n,sigma=R*ACE[1]),rmvn(n,sigma=J*ACE[1]))
    C <- rmvn(2*n,sigma=J*ACE[2])
### Random effects 'death'
    eta1 <- C
### Random effects 'cancer'
    eta2 <- A+C;
### Subject-specific probability of lifetime cancer 
    probcanc <- pnorm(q2+eta2,sd=ACE[3]^.5)   
### Cancer/Death without cancer realizations
    cancertrue <- (runif(length(probcanc))<probcanc)*1+1
### Event times given failure and random effect
    ## inverse P(T<t| eta,cause=1)
    iF1 <- function(u,eta,sd,a,b) (qnorm(u,sd=sd)-eta)/b+a 
    iF2 <- function(u,eta,pr,sd,q,a,b) ialpha(qnorm(u*pr,sd=sd)-eta-q,a=a,b=b)
    u <- runif(length(probcanc))
    t2 <- iF2(u,eta2,probcanc,sd=sqrt(1-sum(ACE[1:2])),q=q2,a=a2,b=b2)
    t1 <- iF1(u,eta1,sd=sqrt(1-sum(ACE[2])),a=a1,b=b1)   
    t0 <- t1; t0[which(cancertrue==2)] <- t2[which(cancertrue==2)] 
### Censoring distributon
    cens <- matrix(rep(rweibull(length(probcanc)/2,shape=exp(logshape),scale=1/exp(logscale)),2),
                   ncol=2)
    ##hist(cens,xlab=c(0,400),200)    
    cause <- cancertrue
    cause[t0>cens] <- 0
    t <- pmin(t0,cens)
    (censtab <- table(cause)/length(cause))    
### Data.frame
    d <- data.frame(t,cause,zyg,cancertrue-1,cens[,1],t0[,1],t0[,2]);
    names(d) <- c("time1","time2","cause1","cause2","zyg","cancertrue1","cancertrue2","cens.time","T01","T02")
### Long format
    dd <- fast.reshape(d)
    dd$cancer <- (dd$cause==2)*1

    if (missing(tt)) tt <- seq(0,max(dd$time))    
    Smz <- J*ACE[1]+J*ACE[2]+I*ACE[3]
    Sdz <- R*ACE[1]+J*ACE[2]+I*ACE[3]
    rr <- alpha(tt,b=b2,a=a2)+q2
    Cmz <- pmvn(upper=cbind(rr,rr),mu=matrix(0,ncol=2,nrow=length(rr)),sigma=Smz)
    Cdz <- pmvn(upper=cbind(rr,rr),mu=matrix(0,ncol=2,nrow=length(rr)),sigma=Sdz)
    
    true <- list(p2=p2,
                 p22mz=pmvn(lower=c(-Inf,-Inf),upper=c(q2,q2),sigma=Smz),
                 p12mz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Smz),
                 p12mz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Smz),
                 p11mz=pmvn(lower=c(q2,q2),upper=c(Inf,Inf),sigma=Smz),
                 p22dz=pmvn(lower=c(-Inf,-Inf),upper=c(q2,q2),sigma=Sdz),
                 p12dz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Sdz),
                 p12dz=pmvn(lower=c(q2,-Inf),upper=c(Inf,q2),sigma=Sdz),
                 p11dz=pmvn(lower=c(q2,q2),upper=c(Inf,Inf),sigma=Sdz),
                 time=tt,
                 F2=F2s(tt),
                 Cmz=Cmz,
                 Cdz=Cdz
                 )
    attributes(dd) <- c(attributes(dd),true)                  
    return(dd)
}


