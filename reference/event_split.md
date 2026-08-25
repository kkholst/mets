# Split event-history records at one or more cut points (SurvSplit)

Splits each row's `(start, time]` interval at every cut point that falls
strictly inside it, replicating covariates onto the new rows and setting
the status of every intermediate (non-final) segment to `cens.code`. Can
also use a subject specific cut-point (event-based) when a variable name
is given as the cuts.

## Usage

``` r
event_split(
  data,
  cuts,
  time = "time",
  status = "status",
  name.id = "id",
  name.start = "start",
  cens.code = 0,
  order.id = TRUE,
  time.group = FALSE,
  group.var = NULL
)
```

## Arguments

- data:

  a data.frame with (at least) the time, status, and optionally start
  columns.

- cuts:

  cut point(s); see Details.

- time:

  name of the time (interval end) column.

- status:

  name of the status/event column.

- name.id:

  name of the subject-id column; created (as `1:n`) if not already
  present.

- name.start:

  name of the interval-start column. If numeric (e.g. `0`), that
  constant is used as the start time for every row and materialised into
  a new column called `start.<value>`.

- cens.code:

  status value assigned to every intermediate split segment (i.e. every
  segment except the last one for a given original row).

- order.id:

  if `TRUE`, sort the result by id then start time.

- time.group:

  if `TRUE`, add one `before.<cut>` 0/1 indicator column per cut
  point/column, flagging whether the segment's start is before that
  particular cut.

- group.var:

  if not `NULL`, the name of a single integer column to add giving which
  time-axis bin the segment falls in: bin 1 is `[0, c1)`, bin 2 is
  `[c1, c2)`, ..., using every distinct cut value seen in `cuts`
  (sorted). Only meaningful when the data has actually been split at
  those same cut points.

## Value

A data.frame with one row per split segment.

## Details

`cuts` can be:

- a single number: one global cut point applied to every row.

- a numeric vector: several global cut points, all applied in one call
  (equivalent to, but faster than, calling `event_splitC` once per cut
  point).

- a single column name (character): a per-row cut point, read from that
  column of `data` (so different rows can be cut at different times).

- a character vector of column names: several per-row cut points at
  once, one set of values per named column.

`NA` in a per-row cut column simply means "no cut" for that row/column
combination.

## Examples

``` r
d <- data.frame(event = round(5 * runif(5), 2),
                 start = 1:5, time = 2 * (1:5),
                 status = rbinom(5, 1, 0.5), x = 1:5)

## 1. single global cut point
event_split(d, cuts = 3.5)
#>   event start time status x id
#> 1  4.82   1.0  2.0      1 1  1
#> 2  2.21   2.0  3.5      0 2  2
#> 3  2.21   3.5  4.0      1 2  2
#> 4  1.85   3.0  3.5      0 3  3
#> 5  1.85   3.5  6.0      1 3  3
#> 6  0.85   4.0  8.0      1 4  4
#> 7  0.27   5.0 10.0      0 5  5

## 2. several global cut points in one call
event_split(d, cuts = c(2, 4, 7))
#>   event start time status x id
#> 1  4.82     1    2      1 1  1
#> 2  2.21     2    4      1 2  2
#> 3  1.85     3    4      0 3  3
#> 4  1.85     4    6      1 3  3
#> 5  0.85     4    7      0 4  4
#> 6  0.85     7    8      1 4  4
#> 7  0.27     5    7      0 5  5
#> 8  0.27     7   10      0 5  5

## 3. per-row cut point taken from an existing column, with an
##    explicit constant start time (start.0 is created automatically)
event_split(d, cuts = "event", name.start = 0)
#>   event start  time status x start.0 id
#> 1  4.82     1  2.00      1 1    0.00  1
#> 2  2.21     2  2.21      0 2    0.00  2
#> 3  2.21     2  4.00      1 2    2.21  2
#> 4  1.85     3  1.85      0 3    0.00  3
#> 5  1.85     3  6.00      1 3    1.85  3
#> 6  0.85     4  0.85      0 4    0.00  4
#> 7  0.85     4  8.00      1 4    0.85  4
#> 8  0.27     5  0.27      0 5    0.00  5
#> 9  0.27     5 10.00      0 5    0.27  5

## 4. several per-row cut points from multiple columns at once
d2 <- d
d2$cutA <- d2$event
d2$cutB <- d2$event + 1
event_split(d2, cuts = c("cutA", "cutB"), name.start = 0)
#>    event start  time status x cutA cutB start.0 id
#> 1   4.82     1  2.00      1 1 4.82 5.82    0.00  1
#> 2   2.21     2  2.21      0 2 2.21 3.21    0.00  2
#> 3   2.21     2  3.21      0 2 2.21 3.21    2.21  2
#> 4   2.21     2  4.00      1 2 2.21 3.21    3.21  2
#> 5   1.85     3  1.85      0 3 1.85 2.85    0.00  3
#> 6   1.85     3  2.85      0 3 1.85 2.85    1.85  3
#> 7   1.85     3  6.00      1 3 1.85 2.85    2.85  3
#> 8   0.85     4  0.85      0 4 0.85 1.85    0.00  4
#> 9   0.85     4  1.85      0 4 0.85 1.85    0.85  4
#> 10  0.85     4  8.00      1 4 0.85 1.85    1.85  4
#> 11  0.27     5  0.27      0 5 0.27 1.27    0.00  5
#> 12  0.27     5  1.27      0 5 0.27 1.27    0.27  5
#> 13  0.27     5 10.00      0 5 0.27 1.27    1.27  5
event_split(d2, cuts = c("cutA", "cutB"), name.start = "start")
#>   event start  time status x cutA cutB id
#> 1  4.82  1.00  2.00      1 1 4.82 5.82  1
#> 2  2.21  2.00  2.21      0 2 2.21 3.21  2
#> 3  2.21  2.21  3.21      0 2 2.21 3.21  2
#> 4  2.21  3.21  4.00      1 2 2.21 3.21  2
#> 5  1.85  3.00  6.00      1 3 1.85 2.85  3
#> 6  0.85  4.00  8.00      1 4 0.85 1.85  4
#> 7  0.27  5.00 10.00      0 5 0.27 1.27  5

## 5. group.var: a single time-axis bin label instead of several
##    before/after indicators
event_split(d, cuts = c(2, 4, 7), group.var = "group")
#>   event start time status x id group
#> 1  4.82     1    2      1 1  1     1
#> 2  2.21     2    4      1 2  2     2
#> 3  1.85     3    4      0 3  3     2
#> 4  1.85     4    6      1 3  3     3
#> 5  0.85     4    7      0 4  4     3
#> 6  0.85     7    8      1 4  4     4
#> 7  0.27     5    7      0 5  5     3
#> 8  0.27     7   10      0 5  5     4
```
