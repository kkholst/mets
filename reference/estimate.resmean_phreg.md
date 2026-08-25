# Extract the stored lava 'estimate' object(s) from a resmean_phreg/cif_yearslost fit

If `times` was given when the object was created, `resmean_phreg`/
`cif_yearslost` already built and stored the relevant `lava`
`"estimate"` object(s) (across strata, one per requested time). This
function simply returns those stored objects – it does not recompute
anything. For `cif_yearslost` objects (which have a `causes` component),
the causes are not jointly estimated (no covariance across causes is
available), so `cause` selects which cause's estimate object(s) to
return. If several `times` were requested, `time` can be used to pull
out the single (plain, class `"estimate"`) object for one time.

## Usage

``` r
# S3 method for class 'resmean_phreg'
estimate(x, contrast = NULL, cause = 1, time = NULL, ...)
```

## Arguments

- x:

  Object of class `"resmean_phreg"` (as returned by `resmean_phreg` or
  `cif_yearslost`).

- contrast:

  Optional contrast (matrix or function), applied via
  [`lava::estimate`](https://kkholst.github.io/lava/reference/estimate.default.html).
  Applied to every element if the result is a list of estimate objects.

- cause:

  Only relevant for `cif_yearslost` objects: which cause's estimate
  object to return.

- time:

  Optional: if several `times` were requested in the original call, pick
  out a single one (returned as a plain `"estimate"` object rather than
  a list).

- ...:

  Passed on to
  [`lava::estimate`](https://kkholst.github.io/lava/reference/estimate.default.html).

## Details

When more than one time (or, for `cif_yearslost`, more than one cause)
is present, the returned list of `"estimate"` objects carries an extra
class, `"resmean_estimate"`, so that
[`print()`](https://rdrr.io/r/base/print.html) shows each element the
same way an individual `"estimate"` object would, and
`summary(x,contrast=...)` applies the contrast to every element
separately (rather than to the list as a whole) – see
[`summary.resmean_estimate`](http://kkholst.github.io/mets/reference/resmean_estimate.md).
