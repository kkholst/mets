# Split a data set and run function

Split a data set and run function

## Usage

``` r
divide.conquer(func = NULL, data, size, splits, id = NULL, ...)
```

## Arguments

- func:

  called function

- data:

  data-frame

- size:

  size of splits

- splits:

  number of splits (ignored if size is given)

- id:

  optional cluster variable

- ...:

  Additional arguments to lower level functions

## Author

Thomas Scheike, Klaus K. Holst

## Examples

``` r
## avoid dependency on timereg
## library(timereg)
## data(TRACE)
## res <- divide.conquer(prop.odds,TRACE,
##        formula=Event(time,status==9)~chf+vf+age,n.sim=0,size=200)
```
