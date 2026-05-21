# Event history object

Constructur for Event History objects

## Usage

``` r
Event(time, time2 = TRUE, cause = NULL, cens.code = 0, ...)
```

## Arguments

- time:

  Time

- time2:

  Time 2

- cause:

  Cause

- cens.code:

  Censoring code (default 0)

- ...:

  Additional arguments

## Value

Object of class Event (a matrix)

## Details

... content for details

## Author

Klaus K. Holst and Thomas Scheike

## Examples

``` r

  t1 <- 1:10
  t2 <- t1+runif(10)
  ca <- rbinom(10,2,0.4)
  (x <- Event(t1,t2,ca))
#>  [1] ( 1; 1.324125:2] ( 2; 2.003554:1] ( 3; 3.168846:2] ( 4; 4.040927:0]
#>  [5] ( 5; 5.128933:1] ( 6; 6.928176:0] ( 7; 7.254323:0] ( 8; 8.685426:1]
#>  [9] ( 9; 9.649160:0] (10;10.267301:0]
```
