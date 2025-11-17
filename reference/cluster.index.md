# Finds subjects related to same cluster

Finds subjects related to same cluster

## Usage

``` r
cluster.index(
  clusters,
  index.type = FALSE,
  num = NULL,
  Rindex = 0,
  mat = NULL,
  return.all = FALSE,
  code.na = NA
)
```

## Arguments

- clusters:

  list of indeces

- index.type:

  if TRUE then already list of integers of index.type

- num:

  to get numbering according to num-type in separate columns

- Rindex:

  index starts with 1, in C is it is 0

- mat:

  to return matrix of indeces

- return.all:

  return all arguments

- code.na:

  how to code missing values

## References

Cluster indeces

## See also

familycluster.index familyclusterWithProbands.index

## Author

Klaus Holst, Thomas Scheike

## Examples

``` r
i<-c(1,1,2,2,1,3)
d<- cluster.index(i)
print(d)
#> $clusters
#> [1] 0 0 1 1 0 2
#> 
#> $maxclust
#> [1] 3
#> 
#> $idclustmat
#>      [,1] [,2] [,3]
#> [1,]    0    1    4
#> [2,]    2    3   NA
#> [3,]    5   NA   NA
#> 
#> $cluster.size
#> [1] 3 2 1
#> 
#> $uniqueclust
#> [1] 3
#> 
#> $firstclustid
#> [1] 0 2 5
#> 

type<-c("m","f","m","c","c","c")
d<- cluster.index(i,num=type,Rindex=1)
print(d)
#> $clusters
#> [1] 0 0 1 1 0 2
#> 
#> $maxclust
#> [1] 3
#> 
#> $idclustmat
#>      [,1] [,2] [,3]
#> [1,]    4    1    0
#> [2,]    3   NA    2
#> [3,]    5   NA   NA
#> 
#> $cluster.size
#> [1] 3 2 1
#> 
#> $uniqueclust
#> [1] 3
#> 
#> $firstclustid
#> [1] 1 3 6
#> 
#> $idclust
#>      [,1] [,2] [,3]
#> [1,]    5    2    1
#> [2,]    4   NA    3
#> [3,]    6   NA   NA
#> 
```
