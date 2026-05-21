# Break ties in event times for recurrent event data

Resolves tied event times in a counting-process dataset by adding a
small random perturbation to duplicated exit times of event rows. This
is a preprocessing step required by
[`recurrent_marginal`](http://kkholst.github.io/mets/reference/recurrent_marginal.md)
and related functions, which assume unique jump times within each
stratum.

## Usage

``` r
tie_breaker(
  data,
  stop = "time",
  start = "entry",
  status = "status",
  id = NULL,
  cause = NULL,
  cens.code = 0,
  exit.unique = TRUE,
  ddt = NULL,
  seed = NULL
)
```

## Arguments

- data:

  A data frame in counting-process format, sorted by subject and time.

- stop:

  Name of the column containing interval exit (stop) times. Default is
  `"time"`.

- start:

  Name of the column containing interval entry (start) times, used to
  update the following row when `id` is supplied. Default is `"entry"`.

- status:

  Name of the column containing event status codes. Default is
  `"status"`.

- id:

  Name of the column containing subject identifiers. If supplied, the
  start time of the next interval for the same subject is adjusted to
  match the perturbed stop time, preserving interval continuity. Default
  is `NULL` (no adjustment made).

- cause:

  Integer vector of status codes that identify events (non-censored
  rows). If `NULL` (default), all non-censoring status values are
  treated as events.

- cens.code:

  Integer code(s) for censoring. Rows with this status are never
  perturbed. Default is `0`.

- exit.unique:

  Logical. If `TRUE` (default), an event time is considered tied
  whenever it coincides with *any* exit time in the data (including
  censored rows). If `FALSE`, only ties among event rows are resolved.

- ddt:

  Maximum perturbation size. Tied event times are shifted by a uniform
  draw on \\\[0, \text{ddt}\]\\. If `NULL` (default), `ddt` is set to
  half the smallest positive gap between any two consecutive exit times
  in the data.

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing the
  perturbations, making results reproducible. Default is `NULL` (no seed
  set).

## Value

The input data frame `data` with tied event exit times perturbed. If
`id` is supplied, an additional logical column `tiebreaker` marks rows
whose start time was adjusted as a consequence of a perturbation in the
preceding row.

## Details

A tie is defined as an event exit time (rows where `status` is a cause
code) that coincides with another exit time in the dataset. When
`exit.unique = TRUE` (default), a tie is flagged whenever an event time
also appears among any other exit times (censored or event). When
`exit.unique = FALSE`, only exact ties between two event rows are
resolved.

Tied event times are perturbed by adding \\U \cdot \delta\\ where \\U
\sim \text{Uniform}(0, 1)\\ and \\\delta\\ is `ddt` (defaulting to half
the smallest observed positive gap between consecutive exit times). When
subject IDs are provided via `id`, the corresponding interval start time
of the immediately following row for the same subject is updated to
maintain a valid counting-process structure, and a logical `tiebreaker`
column is added to flag affected rows.

## See also

[`recurrent_marginal`](http://kkholst.github.io/mets/reference/recurrent_marginal.md),
[`test_logrankRecurrent`](http://kkholst.github.io/mets/reference/test_logrankRecurrent.md)

## Author

Thomas Scheike

## Examples

``` r
data(hfactioncpx12)
hf <- hfactioncpx12

## Check for ties in event exit times
ev <- hf[hf$status == 1, ]
any(duplicated(ev$time))
#> [1] FALSE

## Resolve ties before fitting the marginal mean model
hf_clean <- tie_breaker(hf, stop = "time", start = "entry",
                         status = "status", id = "id",
                         cause = 1, cens.code = 0)

out <- recurrent_marginal(Event(entry, time, status) ~ cluster(id),
                           data = hf_clean, cause = 1, death.code = 2)
summary(out, times = 1:5)
#> [[1]]
#>        new.time      mean         se   CI-2.5% CI-97.5% strata
#> 608           1 0.8282358 0.04844543 0.7385251 0.928844      0
#> 1053          2 1.5139493 0.07039884 1.3820710 1.658412      0
#> 1282          3 2.0244982 0.08351867 1.8672476 2.194992      0
#> 1392          4 2.5004732 0.10843166 2.2967320 2.722288      0
#> 1392.1        5 2.5004732 0.10843166 2.2967320 2.722288      0
#> 
```
