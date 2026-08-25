# Print/summary for a list of estimate objects (resmean_phreg/cif_yearslost)

A `"resmean_estimate"` object is simply a named list of `lava`
`"estimate"` objects (one per time, or one per cause, or nested one per
cause per time), as returned by
`estimate.resmean_phreg`/`summary.resmean_phreg` when more than one such
object is present. `print` shows each element the same way a single
`"estimate"` object would (no extra "Length Class Mode" summary table).
`summary(...,contrast=)` applies the contrast to every element
separately (recursing into any nested `"resmean_estimate"` elements) and
shows the usual `summary.estimate` table (with its null-hypothesis test)
for each; with no `contrast`, it just returns the object unchanged.

## Usage

``` r
# S3 method for class 'resmean_estimate'
print(x, ...)

# S3 method for class 'resmean_estimate'
summary(object, contrast = NULL, ...)
```

## Arguments

- x:

  A `"resmean_estimate"` object (for `print`).

- ...:

  Passed to `print`/`summary` of each element.

- object:

  A `"resmean_estimate"` object (for `summary`).

- contrast:

  Optional contrast (matrix or function), applied via
  [`lava::estimate`](https://kkholst.github.io/lava/reference/estimate.default.html)
  to every element. Only used by `summary`.
