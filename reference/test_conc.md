# Compare Two Concordance Estimates

Performs a Pepe-Mori type test to compare two concordance estimates
(e.g., MZ vs DZ twins). The test evaluates whether the concordance
functions differ significantly over time.

## Usage

``` r
test_conc(conc1, conc2, same.cluster = FALSE)
```

## Arguments

- conc1:

  Concordance estimate of group 1.

- conc2:

  Concordance estimate of group 2.

- same.cluster:

  Logical; if FALSE, groups are assumed independent. If TRUE, estimates
  are based on the same data (e.g., paired twins).

## Value

An object of class `"testconc"` containing:

- test:

  Matrix with cumulative difference, standard error, z-value, and
  p-value.

- mintime, maxtime:

  Time range used.

- same.cluster:

  Logical flag.

## See also

[`test_casewise`](http://kkholst.github.io/mets/reference/test_casewise.md)

## Author

Thomas Scheike
