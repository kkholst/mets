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
#>  [1] ( 1; 1.171086:0] ( 2; 2.333062:0] ( 3; 3.183510:1] ( 4; 4.128726:1]
#>  [5] ( 5; 5.503349:1] ( 6; 6.767420:1] ( 7; 7.931043:1] ( 8; 8.471621:1]
#>  [9] ( 9; 9.953553:1] (10;10.700965:0]
```
