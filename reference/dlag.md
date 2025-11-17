# Lag operator

Lag operator

## Usage

``` r
dlag(data, x, k = 1, combine = TRUE, simplify = TRUE, names, ...)
```

## Arguments

- data:

  data.frame or vector

- x:

  optional column names or formula

- k:

  lag (vector of integers)

- combine:

  combine results with original data.frame

- simplify:

  Return vector if possible

- names:

  optional new column names

- ...:

  additional arguments to lower level functions

## Examples

``` r
d <- data.frame(y=1:10,x=c(10:1))
dlag(d,k=1:2)
#>     y  x y.1 y.2 x.1 x.2
#> 1   1 10  NA  NA  NA  NA
#> 2   2  9   1  NA  10  NA
#> 3   3  8   2   1   9  10
#> 4   4  7   3   2   8   9
#> 5   5  6   4   3   7   8
#> 6   6  5   5   4   6   7
#> 7   7  4   6   5   5   6
#> 8   8  3   7   6   4   5
#> 9   9  2   8   7   3   4
#> 10 10  1   9   8   2   3
dlag(d,~x,k=0:1)
#>     y  x x.0 x.1
#> 1   1 10  10  NA
#> 2   2  9   9  10
#> 3   3  8   8   9
#> 4   4  7   7   8
#> 5   5  6   6   7
#> 6   6  5   5   6
#> 7   7  4   4   5
#> 8   8  3   3   4
#> 9   9  2   2   3
#> 10 10  1   1   2
dlag(d$x,k=1)
#>  [1] NA 10  9  8  7  6  5  4  3  2
dlag(d$x,k=-1:2, names=letters[1:4])
#>        a  b  c  d
#>  [1,]  9 10 NA NA
#>  [2,]  8  9 10 NA
#>  [3,]  7  8  9 10
#>  [4,]  6  7  8  9
#>  [5,]  5  6  7  8
#>  [6,]  4  5  6  7
#>  [7,]  3  4  5  6
#>  [8,]  2  3  4  5
#>  [9,]  1  2  3  4
#> [10,] NA  1  2  3
```
