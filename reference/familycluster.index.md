# Finds all pairs within a cluster (family)

Finds all pairs within a cluster (family)

## Usage

``` r
familycluster.index(clusters, index.type = FALSE, num = NULL, Rindex = 1)
```

## Arguments

- clusters:

  list of indeces

- index.type:

  argument of cluster index

- num:

  num

- Rindex:

  index starts with 1 in R, and 0 in C

## References

Cluster indeces

## See also

cluster.index familyclusterWithProbands.index

## Author

Klaus Holst, Thomas Scheike

## Examples

``` r
i<-c(1,1,2,2,1,3)
d<- familycluster.index(i)
print(d)
#> $familypairindex
#> [1] 1 2 1 5 2 5 3 4
#> 
#> $subfamilyindex
#> [1] 0 0 1 1 2 2 3 3
#> 
#> $pairs
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    5
#> [3,]    2    5
#> [4,]    3    4
#> 
#> $clusters
#> [1] 1 1 1 2
#> 
```
