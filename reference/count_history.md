# Compute cumulative event counts as time-dependent covariates

For each subject and each row of `data`, counts the number of prior
events of specified types in the recurrent event history. The resulting
count columns can be used as time-dependent covariates in subsequent
models, e.g. to capture event-history dependence in the recurrent event
rate.

## Usage

``` r
count_history(
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

  A data frame in counting-process format, with one row per event
  interval per subject.

- status:

  Name of the column containing event status codes. Default is
  `"status"`.

- id:

  Name of the column containing subject identifiers. Default is `"id"`.

- types:

  Integer vector of status codes to count. Each value in `types`
  generates one new count column (when `multitype = FALSE`) or
  contributes to a single combined count (when `multitype = TRUE`).
  Default is `1`.

- names.count:

  Prefix for the names of the new count columns. The status code is
  appended, e.g. `"Count1"`, `"Count2"`. Default is `"Count"`.

- lag:

  Logical. If `TRUE` (default), the count at each row is the number of
  events strictly before the current time (\\N(t-)\\). If `FALSE`,
  events at the current time are included (\\N(t)\\).

- multitype:

  Logical. If `TRUE`, events with status in `types` are aggregated into
  a single count column, optionally weighted by `marks`. If `FALSE`
  (default), a separate count column is created for each value in
  `types`.

- marks:

  Optional numeric vector of weights applied to events when
  `multitype = TRUE`. If `NULL` (default), each event has weight 1.

## Value

The input data frame `data` with one or more new integer columns
appended. With `multitype = FALSE`, columns are named
`paste0(names.count, k)` for each `k` in `types`; with
`multitype = TRUE`, a single column named
`paste0(names.count, types[1])` is added. An internal bookkeeping column
`lbnr__id` is also added.

## Details

When `lag = TRUE` (default), the count at each row reflects events that
occurred strictly before the current time point (i.e. \\N(t-)\\), making
it suitable as a left-continuous covariate in counting-process models.

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
hf <- hfactioncpx12
dtable(hf, ~status)
#> 
#> status
#>    0    1    2 
#>  617 1391  124 
#> 

## Separate counts for event types 1 and 2
rr <- count_history(hf, types = 1:2, id = "id", status = "status")
dtable(rr, ~"Count*" + status, level = 1)
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
