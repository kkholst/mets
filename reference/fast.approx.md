# Fast approximation

Fast approximation

## Usage

``` r
fast.approx(
  time,
  new.time,
  equal = FALSE,
  type = c("nearest", "right", "left"),
  sorted = FALSE,
  ...
)
```

## Arguments

- time:

  Original ordered time points

- new.time:

  New time points

- equal:

  If TRUE a list is returned with additional element

- type:

  Type of matching, nearest index, nearest greater than or equal
  (right), number of elements smaller than y otherwise the closest value
  above new.time is returned.

- sorted:

  Set to true if new.time is already sorted

- ...:

  Optional additional arguments

## Author

Klaus K. Holst

## Examples

``` r
id <- c(1,1,2,2,7,7,10,10)
fast.approx(unique(id),id)
#> [1] 1 1 2 2 3 3 4 4

t <- 0:6
n <- c(-1,0,0.1,0.9,1,1.1,1.2,6,6.5)
fast.approx(t,n,type="left")
#> [1] 0 1 1 1 2 2 2 7 7
```
