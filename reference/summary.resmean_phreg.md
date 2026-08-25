# Summary for resmean_phreg/cif_yearslost fits

If `times=NULL` was used when creating the object, the summary is the
full curve at all event times – exactly what `plot` would show – as a
data frame with the estimate, its standard error, confidence interval,
and (for `resmean_phreg`) a years-lost column (`time - rmst`); for
`cif_yearslost` one such table is given per cause.

## Usage

``` r
# S3 method for class 'resmean_phreg'
summary(object, contrast = NULL, cause = NULL, time = NULL, level = 0.95, ...)
```

## Arguments

- object:

  Object of class `"resmean_phreg"`.

- contrast:

  Optional contrast (matrix or function), applied via
  [`lava::estimate`](https://kkholst.github.io/lava/reference/estimate.default.html)
  before summarizing; triggers the usual `summary.estimate` test output.
  Only used when `times` was given.

- cause:

  Only relevant for `cif_yearslost` objects: restrict the summary to
  this cause (default: all causes).

- time:

  Optional: restrict the summary to a single time (only relevant if
  several `times` were given in the original call).

- level:

  Confidence level, only used for the no-`times` curve-summary.

- ...:

  Additional arguments passed to `summary.estimate` (only used together
  with `contrast`) or to `conftype` (when `times=NULL`).

## Details

If `times` was given, the object already stores the relevant `lava`
`"estimate"` object(s) (see
[`estimate.resmean_phreg`](http://kkholst.github.io/mets/reference/estimate.resmean_phreg.md)),
and this just returns those stored object(s) as-is – printing one shows
the usual Estimate/Std.Err/CI/P-value table. No joint null-hypothesis
test is added by default; pass `contrast` (a matrix or function, applied
via
[`lava::estimate`](https://kkholst.github.io/lava/reference/estimate.default.html))
if you want to test a specific linear combination, in which case the
usual `summary.estimate` output (including the test) is shown for that
contrast (see
[`summary.resmean_estimate`](http://kkholst.github.io/mets/reference/resmean_estimate.md)
for how this is applied when there is more than one time/cause). For
`cif_yearslost` objects, `cause=NULL` (the default) summarizes every
cause (returned as a named list); `cause=1` (say) restricts the summary
to that cause. `time` similarly restricts to a single time (only
relevant if `length(times)>1` in the original call).
