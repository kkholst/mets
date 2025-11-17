# Event split with two time-scales, time and gaptime

Cuts time for two time-scales, as event.split

## Usage

``` r
EventSplit2(
  data,
  time = "time",
  status = "status",
  entry = "start",
  cuts = "cuts",
  name.id = "id",
  gaptime = NULL,
  gaptime.entry = NULL,
  cuttime = c("time", "gaptime"),
  cens.code = 0,
  order.id = TRUE
)
```

## Arguments

- data:

  data to be split

- time:

  time variable.

- status:

  status variable.

- entry:

  name of entry variable.

- cuts:

  cuts variable or numeric cut (only one value)

- name.id:

  name of id variable.

- gaptime:

  gaptime variable.

- gaptime.entry:

  name of entry variable for gaptime.

- cuttime:

  to cut after time or gaptime

- cens.code:

  code for the censoring.

- order.id:

  order data after id and start.

## Author

Thomas Scheike

## Examples

``` r
rr  <- data.frame(time=c(500,1000),start=c(0,500),status=c(1,1),id=c(1,1))
rr$gaptime <-  rr$time-rr$start
rr$gapstart <- 0

rr1 <- EventSplit2(rr,cuts=600,cuttime="time",   gaptime="gaptime",gaptime.entry="gapstart")
rr2 <- EventSplit2(rr1,cuts=100,cuttime="gaptime",gaptime="gaptime",gaptime.entry="gapstart")

dlist(rr1,start-time+status+gapstart+gaptime~id)
#> id: 1
#>     start time status gapstart gaptime
#> 1     0    500 1        0      500    
#> 2   500    600 0        0      100    
#> 2.1 600   1000 1      100      500    
dlist(rr2,start-time+status+gapstart+gaptime~id)
#> id: 1
#>     start time status gapstart gaptime
#> 1     0    100 0        0      100    
#> 1.1 100    500 1      100      500    
#> 2   500    600 0        0      100    
#> 2.1 600   1000 1      100      500    
```
