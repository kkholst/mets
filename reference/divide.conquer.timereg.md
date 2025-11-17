# Split a data set and run function from timereg and aggregate

Split a data set and run function of cox-aalen type and aggregate
results

## Usage

``` r
divide.conquer.timereg(func = NULL, data, size, ...)
```

## Arguments

- func:

  called function

- data:

  data-frame

- size:

  size of splits

- ...:

  Additional arguments to lower level functions

## Author

Thomas Scheike, Klaus K. Holst

## Examples

``` r
## library(timereg)
## data(TRACE)
## a <- divide.conquer.timereg(prop.odds,TRACE,
##                            formula=Event(time,status==9)~chf+vf+age,n.sim=0,size=200)
## coef(a)
## a2 <- divide.conquer.timereg(prop.odds,TRACE,
##                              formula=Event(time,status==9)~chf+vf+age,n.sim=0,size=500)
## coef(a2)
## 
##if (interactive()) {
##par(mfrow=c(1,1))
##plot(a,xlim=c(0,8),ylim=c(0,0.01))
##par(new=TRUE)
##plot(a2,xlim=c(0,8),ylim=c(0,0.01))
##}
```
