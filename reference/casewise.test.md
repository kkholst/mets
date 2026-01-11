# Estimates the casewise concordance based on Concordance and marginal estimate using timereg and performs test for independence

Estimates the casewise concordance based on Concordance and marginal
estimate using timereg and performs test for independence

## Usage

``` r
casewise.test(conc, marg, test = "no-test", p = 0.01)
```

## Arguments

- conc:

  Concordance

- marg:

  Marginal estimate

- test:

  Type of test for independence assumption. "conc" makes test on
  concordance scale and "case" means a test on the casewise concordance

- p:

  check that marginal probability is greater at some point than p

## Details

Uses cluster based conservative standard errors for marginal and
sometimes only the uncertainty of the concordance estimates. This works
prettey well, alternatively one can use also the funcions Casewise for a
specific time point

## Author

Thomas Scheike

## Examples

``` r
 ## Reduce Ex.Timings
library("timereg")
data("prt",package="mets");
prt <- force.same.cens(prt,cause="status")

prt <- prt[which(prt$id %in% sample(unique(prt$id),7500)),]
### marginal cumulative incidence of prostate cancer
times <- seq(60,100,by=2)
outm <- timereg::comp.risk(Event(time,status)~+1,data=prt,cause=2,times=times)

cifmz <- predict(outm,X=1,uniform=0,resample.iid=1)
cifdz <- predict(outm,X=1,uniform=0,resample.iid=1)

### concordance for MZ and DZ twins
cc <- bicomprisk(Event(time,status)~strata(zyg)+id(id),
                 data=prt,cause=c(2,2))
#> Strata 'DZ'
#> Strata 'MZ'
cdz <- cc$model$"DZ"
cmz <- cc$model$"MZ"

### To compute casewise cluster argument must be passed on,
###  here with a max of 100 to limit comp-time
outm <- timereg::comp.risk(Event(time,status)~+1,data=prt,
                 cause=2,times=times,max.clust=100)
cifmz <- predict(outm,X=1,uniform=0,resample.iid=1)
cc <- bicomprisk(Event(time,status)~strata(zyg)+id(id),data=prt,
                cause=c(2,2),se.clusters=outm$clusters)
#> Strata 'DZ'
#> Strata 'MZ'
cdz <- cc$model$"DZ"
cmz <- cc$model$"MZ"

cdz <- casewise.test(cdz,cifmz,test="case") ## test based on casewise
cmz <- casewise.test(cmz,cifmz,test="conc") ## based on concordance

plot(cmz,ylim=c(0,0.7),xlim=c(60,100))
par(new=TRUE)
plot(cdz,ylim=c(0,0.7),xlim=c(60,100))


slope.process(cdz$casewise[,1],cdz$casewise[,2],iid=cdz$casewise.iid)
#> $intercept
#> (Intercept) 
#>   0.1520835 
#> 
#> $slope
#>      ctime 
#> 0.01849925 
#> 
#> $se.slope
#> (Intercept)       ctime 
#>  0.06273584  0.04470886 
#> 
#> $pval.slope
#>     ctime 
#> 0.6790415 
#> 

slope.process(cmz$casewise[,1],cmz$casewise[,2],iid=cmz$casewise.iid)
#> $intercept
#> (Intercept) 
#>   0.4456058 
#> 
#> $slope
#>        ctime 
#> 0.0004679675 
#> 
#> $se.slope
#> (Intercept)       ctime 
#>  0.08924422  0.07231706 
#> 
#> $pval.slope
#>     ctime 
#> 0.9948369 
#> 

```
