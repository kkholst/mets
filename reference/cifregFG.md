# Fine-Gray Cumulative Incidence Function Regression

Convenience wrapper for `cifreg` that specifically fits the Fine-Gray
model by setting `propodds=NULL`. This provides subdistribution hazard
ratio interpretations.

## Usage

``` r
cifregFG(formula, data, propodds = NULL, ...)
```

## Arguments

- formula:

  Formula with an 'Event' outcome.

- data:

  Data frame.

- propodds:

  Set to `NULL` to fit the Fine-Gray model (default behavior of this
  function).

- ...:

  Additional arguments passed to `cifreg`.

## Value

An object of class `"cifreg"` (extending `"phreg"`) with
`propodds=NULL`.

## See also

[`cifreg`](http://kkholst.github.io/mets/reference/cifreg.md)

## Author

Thomas Scheike
