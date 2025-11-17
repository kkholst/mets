# Finds all pairs within a cluster (famly) with the proband (case/control)

second column of pairs are the probands and the first column the related
subjects

## Usage

``` r
familyclusterWithProbands.index(
  clusters,
  probands,
  index.type = FALSE,
  num = NULL,
  Rindex = 1
)
```

## Arguments

- clusters:

  list of indeces giving the clusters (families)

- probands:

  list of 0,1 where 1 specifices which of the subjects that are probands

- index.type:

  argument passed to other functions

- num:

  argument passed to other functions

- Rindex:

  index starts with 1, in C is it is 0

## References

Cluster indeces

## See also

familycluster.index cluster.index

## Author

Klaus Holst, Thomas Scheike

## Examples

``` r
i<-c(1,1,2,2,1,3)
p<-c(1,0,0,1,0,1)
d<- familyclusterWithProbands.index(i,p)
print(d)
#> $familypairindex
#> [1] 2 1 5 1 3 4
#> 
#> $subfamilyindex
#> [1] 0 0 1 1 3 3
#> 
#> $pairs
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    5    1
#> [3,]    3    4
#> 
#> $clusters
#> [1] 1 1 2
#> 
```
