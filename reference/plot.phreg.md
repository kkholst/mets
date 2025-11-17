# Plotting the baselines of stratified Cox

Plotting the baselines of stratified Cox

## Usage

``` r
# S3 method for class 'phreg'
plot(x, ...)
```

## Arguments

- x:

  phreg object

- ...:

  Additional arguments to baseplot funtion

## Author

Klaus K. Holst, Thomas Scheike

## Examples

``` r
data(TRACE)
dcut(TRACE) <- ~.
out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)

par(mfrow=c(2,2))
plot(out1)
plot(out1,stratas=c(0,3))
plot(out1,stratas=c(0,3),col=2:3,lty=1:2,se=TRUE)
plot(out1,stratas=c(0),col=2,lty=2,se=TRUE,polygon=FALSE)

plot(out1,stratas=c(0),col=matrix(c(2,1,3),1,3),lty=matrix(c(1,2,3),1,3),se=TRUE,polygon=FALSE)
```
