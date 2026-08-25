# Print method for resmean_phreg/cif_yearslost fits

If `times` was given at call time, printing shows the stored estimate
object(s) (same as `summary(x)`) since that output is already compact.

## Usage

``` r
# S3 method for class 'resmean_phreg'
print(x, ...)
```

## Arguments

- x:

  Object of class `"resmean_phreg"`.

- ...:

  Passed on to `summary.resmean_phreg` when `times` was given.

## Details

If `times=NULL`, the full curve (all event times) can be long, so
printing the object just gives a short description instead of dumping it
– use `summary(x)` (or `plot(x)`) to see the full rmst/years-lost curve.
