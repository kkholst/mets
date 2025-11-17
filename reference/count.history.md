# Counts the number of previous events of two types for recurrent events processes

Counts the number of previous events of two types for recurrent events
processes

## Usage

``` r
count.history(
  data,
  status = "status",
  id = "id",
  types = 1,
  names.count = "Count",
  lag = TRUE,
  multitype = FALSE,
  marks = NULL
)
```

## Arguments

- data:

  data-frame

- status:

  name of status

- id:

  id

- types:

  types of the events (code) related to status (multiple values
  possible)

- names.count:

  name of Counts, for example Count1 Count2 when types=c(1,2)

- lag:

  if true counts previously observed, and if lag=FALSE counts up to know

- multitype, :

  if multitype is true then counts when status "in" types, otherwise
  counts for each value of type, types=c(1,2)

- marks:

  values related to status ("in" types), counts marks for types, only
  when multitype=TRUE

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
hf <- hfactioncpx12
dtable(hf,~status)
#> 
#> status
#>    0    1    2 
#>  617 1391  124 
#> 
rr <-  count.history(hf,types=1:2,id="id",status="status")
dtable(rr,~"Count*"+status,level=1)
#> 
#> Count1
#>   0   1   2   3   4   5   6   7 
#> 741 507 319 209 146  97  67  46 
#> 
#> Count2
#>    0 
#> 2132 
#> 
#> status
#>    0    1    2 
#>  617 1391  124 
#> 
```
