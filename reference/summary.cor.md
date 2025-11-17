# Summary for dependence models for competing risks

Computes concordance and casewise concordance for dependence models for
competing risks models of the type cor.cif, rr.cif or or.cif for the
given cumulative incidences and the different dependence measures in the
object.

## Usage

``` r
# S3 method for class 'cor'
summary(object, marg.cif = NULL, marg.cif2 = NULL, digits = 3, ...)
```

## Arguments

- object:

  object from cor.cif rr.cif or or.cif for dependence between competing
  risks data for two causes.

- marg.cif:

  a number that gives the cumulative incidence in one time point for
  which concordance and casewise concordance are computed.

- marg.cif2:

  the cumulative incidence for cause 2 for concordance and casewise
  concordance are computed. Default is that it is the same as marg.cif.

- digits:

  digits in output.

- ...:

  Additional arguments.

## Value

prints summary for dependence model.

- casewise:

  gives casewise concordance that is, probability of cause 2 (related to
  cif2) given that cause 1 (related to cif1) has occured.

- concordance:

  gives concordance that is, probability of cause 2 (related to cif2)
  and cause 1 (related to cif1).

- cif1:

  cumulative incidence for cause1.

- cif2:

  cumulative incidence for cause1.

## References

Cross odds ratio Modelling of dependence for Multivariate Competing
Risks Data, Scheike and Sun (2012), Biostatistics.

A Semiparametric Random Effects Model for Multivariate Competing Risks
Data, Scheike, Zhang, Sun, Jensen (2010), Biometrika.

## Author

Thomas Scheike

## Examples

``` r
## library("timereg")
## data("multcif",package="mets") # simulated data 
## multcif$cause[multcif$cause==0] <- 2
##  
## times=seq(0.1,3,by=0.1) # to speed up computations use only these time-points
## add <- timereg::comp.risk(Event(time,cause)~+1+cluster(id),
##                           data=multcif,n.sim=0,times=times,cause=1)
###
## out1<-cor.cif(add,data=multcif,cause1=1,cause2=1,theta=log(2+1))
## summary(out1)
## 
## pad <- predict(add,X=1,se=0,uniform=0)
## summary(out1,marg.cif=pad)
```
