# event_split (SurvSplit).

Constructs start stop formulation of event time data after a variable in
the data.set. Similar to SurvSplit of the survival package but can also
split after random time given in data frame.

## Usage

``` r
event_split(
  data,
  time = "time",
  status = "status",
  cuts = "cuts",
  name.id = "id",
  name.start = "start",
  cens.code = 0,
  order.id = TRUE,
  time.group = FALSE
)
```

## Arguments

- data:

  data to be split

- time:

  time variable.

- status:

  status variable.

- cuts:

  cuts variable or numeric cut (only one value)

- name.id:

  name of id variable.

- name.start:

  name of start variable in data, start can also be numeric "0"

- cens.code:

  code for the censoring.

- order.id:

  order data after id and start.

- time.group:

  make variable "before"."cut" that keeps track of wether start,stop is
  before (1) or after cut (0).

## Author

Thomas Scheike

## Examples

``` r

set.seed(1)
d <- data.frame(event=round(5*runif(5),2),start=1:5,time=2*1:5,
    status=rbinom(5,1,0.5),x=1:5)
d
#>   event start time status x
#> 1  1.33     1    2      1 1
#> 2  1.86     2    4      1 2
#> 3  2.86     3    6      1 3
#> 4  4.54     4    8      1 4
#> 5  1.01     5   10      0 5

d0 <- event_split(d,cuts="event",name.start=0)
d0
#>    event start  time status x start.0 id
#> 1   1.33     1  1.33      0 1    0.00  1
#> 2   1.33     1  2.00      1 1    1.33  1
#> 3   1.86     2  1.86      0 2    0.00  2
#> 4   1.86     2  4.00      1 2    1.86  2
#> 5   2.86     3  2.86      0 3    0.00  3
#> 6   2.86     3  6.00      1 3    2.86  3
#> 7   4.54     4  4.54      0 4    0.00  4
#> 8   4.54     4  8.00      1 4    4.54  4
#> 9   1.01     5  1.01      0 5    0.00  5
#> 10  1.01     5 10.00      0 5    1.01  5

dd <- event_split(d,cuts="event")
dd
#>   event start  time status x id
#> 1  1.33  1.00  1.33      0 1  1
#> 2  1.33  1.33  2.00      1 1  1
#> 3  1.86  2.00  4.00      1 2  2
#> 4  2.86  3.00  6.00      1 3  3
#> 5  4.54  4.00  4.54      0 4  4
#> 6  4.54  4.54  8.00      1 4  4
#> 7  1.01  5.00 10.00      0 5  5
ddd <- event_split(dd,cuts=3.5)
ddd
#>   event start  time status x id cut.3.5
#> 1  1.33  1.00  1.33      0 1  1     3.5
#> 2  1.33  1.33  2.00      1 1  1     3.5
#> 3  1.86  2.00  3.50      0 2  2     3.5
#> 4  1.86  3.50  4.00      1 2  2     3.5
#> 5  2.86  3.00  3.50      0 3  3     3.5
#> 6  2.86  3.50  6.00      1 3  3     3.5
#> 7  4.54  4.00  4.54      0 4  4     3.5
#> 8  4.54  4.54  8.00      1 4  4     3.5
#> 9  1.01  5.00 10.00      0 5  5     3.5
event_split(ddd,cuts=5.5)
#>    event start  time status x id cut.3.5 cut.5.5
#> 1   1.33  1.00  1.33      0 1  1     3.5     5.5
#> 2   1.33  1.33  2.00      1 1  1     3.5     5.5
#> 3   1.86  2.00  3.50      0 2  2     3.5     5.5
#> 4   1.86  3.50  4.00      1 2  2     3.5     5.5
#> 5   2.86  3.00  3.50      0 3  3     3.5     5.5
#> 6   2.86  3.50  5.50      0 3  3     3.5     5.5
#> 7   2.86  5.50  6.00      1 3  3     3.5     5.5
#> 8   4.54  4.00  4.54      0 4  4     3.5     5.5
#> 9   4.54  4.54  5.50      0 4  4     3.5     5.5
#> 10  4.54  5.50  8.00      1 4  4     3.5     5.5
#> 11  1.01  5.00  5.50      0 5  5     3.5     5.5
#> 12  1.01  5.50 10.00      0 5  5     3.5     5.5

### successive cutting for many values 
dd <- d
for  (cuts in seq(2,3,by=0.3)) dd <- event_split(dd,cuts=cuts)
dd
#>   event start time status x cut.2 id cut.2.3 cut.2.6 cut.2.9
#> 1  1.33   1.0  2.0      1 1     2  1     2.3     2.6     2.9
#> 2  1.86   2.0  2.3      0 2     2  2     2.3     2.6     2.9
#> 3  1.86   2.3  2.6      0 2     2  2     2.3     2.6     2.9
#> 4  1.86   2.6  2.9      0 2     2  2     2.3     2.6     2.9
#> 5  1.86   2.9  4.0      1 2     2  2     2.3     2.6     2.9
#> 6  2.86   3.0  6.0      1 3     2  3     2.3     2.6     2.9
#> 7  4.54   4.0  8.0      1 4     2  4     2.3     2.6     2.9
#> 8  1.01   5.0 10.0      0 5     2  5     2.3     2.6     2.9

###########################################################################
### same but for situation with multiple events along the time-axis
###########################################################################
d <- data.frame(event1=1:5+runif(5)*0.5,start=1:5,time=2*1:5,
    status=rbinom(5,1,0.5),x=1:5,start0=0)
d$event2 <- d$event1+0.2
d$event2[4:5] <- NA 
d
#>     event1 start time status x start0   event2
#> 1 1.102987     1    2      0 1      0 1.302987
#> 2 2.088278     2    4      1 2      0 2.288278
#> 3 3.343511     3    6      1 3      0 3.543511
#> 4 4.192052     4    8      0 4      0       NA
#> 5 5.384921     5   10      1 5      0       NA

d0 <- event_split(d,cuts="event1",name.start="start",time="time",status="status")
d0
#>      event1    start      time status x start0   event2 id
#> 1  1.102987 1.000000  1.102987      0 1      0 1.302987  1
#> 2  1.102987 1.102987  2.000000      0 1      0 1.302987  1
#> 3  2.088278 2.000000  2.088278      0 2      0 2.288278  2
#> 4  2.088278 2.088278  4.000000      1 2      0 2.288278  2
#> 5  3.343511 3.000000  3.343511      0 3      0 3.543511  3
#> 6  3.343511 3.343511  6.000000      1 3      0 3.543511  3
#> 7  4.192052 4.000000  4.192052      0 4      0       NA  4
#> 8  4.192052 4.192052  8.000000      0 4      0       NA  4
#> 9  5.384921 5.000000  5.384921      0 5      0       NA  5
#> 10 5.384921 5.384921 10.000000      1 5      0       NA  5
###
d00 <- event_split(d0,cuts="event2",name.start="start",time="time",status="status")
d00
#>      event1    start      time status x start0   event2 id
#> 1  1.102987 1.000000  1.102987      0 1      0 1.302987  1
#> 2  1.102987 1.102987  1.302987      0 1      0 1.302987  1
#> 3  1.102987 1.302987  2.000000      0 1      0 1.302987  1
#> 4  2.088278 2.000000  2.088278      0 2      0 2.288278  2
#> 5  2.088278 2.088278  2.288278      0 2      0 2.288278  2
#> 6  2.088278 2.288278  4.000000      1 2      0 2.288278  2
#> 7  3.343511 3.000000  3.343511      0 3      0 3.543511  3
#> 8  3.343511 3.343511  3.543511      0 3      0 3.543511  3
#> 9  3.343511 3.543511  6.000000      1 3      0 3.543511  3
#> 10 4.192052 4.000000  4.192052      0 4      0       NA  4
#> 11 4.192052 4.192052  8.000000      0 4      0       NA  4
#> 12 5.384921 5.000000  5.384921      0 5      0       NA  5
#> 13 5.384921 5.384921 10.000000      1 5      0       NA  5
```
