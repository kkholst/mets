# Relative risk for additive gamma model

Computes the relative risk for additive gamma model at time 0

## Usage

``` r
EVaddGam(theta, x1, x2, thetades, ags)
```

## Arguments

- theta:

  theta

- x1:

  x1

- x2:

  x2

- thetades:

  thetades

- ags:

  ags

## References

Eriksson and Scheike (2015), Additive Gamma frailty models for competing
risks data, Biometrics (2015)

## Author

Thomas Scheike

## Examples

``` r
lam0 <- c(0.5,0.3)
pars <- c(1,1,1,1,0,1)
## genetic random effects, cause1, cause2 and overall
parg <- pars[c(1,3,5)]
## environmental random effects, cause1, cause2 and overall
parc <- pars[c(2,4,6)]

## simulate competing risks with two causes with hazards 0.5 and 0.3
## ace for each cause, and overall ace
out <- simCompete.twin.ace(10000,parg,parc,0,2,lam0=lam0,overall=1,all.sum=1)

## setting up design for running the model
mm <- familycluster.index(out$cluster)
head(mm$familypairindex,n=10)
#>  [1]  1  2  3  4  5  6  7  8  9 10
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
#>           [,1]  [,2]
#>  [9989,] 19977 19978
#>  [9990,] 19979 19980
#>  [9991,] 19981 19982
#>  [9992,] 19983 19984
#>  [9993,] 19985 19986
#>  [9994,] 19987 19988
#>  [9995,] 19989 19990
#>  [9996,] 19991 19992
#>  [9997,] 19993 19994
#>  [9998,] 19995 19996
#>  [9999,] 19997 19998
#> [10000,] 19999 20000
#
kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5

# dout <- make.pairwise.design.competing(pairs,kinship,
#          type="ace",compete=length(lam0),overall=1)
# head(dout$ant.rvs)
## MZ
# dim(dout$theta.des)
# dout$random.design[,,1]
## DZ
# dout$theta.des[,,nrow(pairs)]
# dout$random.design[,,nrow(pairs)]
#
# thetades <- dout$theta.des[,,1]
# x <- dout$random.design[,,1]
# x
##EVaddGam(rep(1,6),x[1,],x[3,],thetades,matrix(1,18,6))

# thetades <- dout$theta.des[,,nrow(out)/2]
# x <- dout$random.design[,,nrow(out)/2]
##EVaddGam(rep(1,6),x[1,],x[4,],thetades,matrix(1,18,6))
```
