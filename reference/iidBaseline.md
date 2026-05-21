# Influence Functions or IID Decomposition of Baseline

Computes the influence functions for the baseline cumulative hazard (and
optionally regression coefficients) for `phreg`, `recreg`, or `cifregFG`
objects.

## Usage

``` r
iidBaseline(
  object,
  time = NULL,
  ft = NULL,
  fixbeta = NULL,
  beta.iid = NULL,
  tminus = FALSE,
  ...
)
```

## Arguments

- object:

  Object of class `"phreg"`, `"recreg"`, or `"cifregFG"`.

- time:

  Time point for baseline IID (required).

- ft:

  Function to compute IID of baseline integrated against \\f(t)\\.

- fixbeta:

  Logical; if `TRUE`, fixes the coefficients (useful for specific
  tests).

- beta.iid:

  Optional matrix of beta influence functions to use.

- tminus:

  Logical; if `TRUE`, computes predictions at \\t-\\ (strictly before
  \\t\\), useful for IPCW techniques.

- ...:

  Additional arguments passed to lower-level functions.

## Value

An object of class `"iidBaseline"` containing:

- time:

  Time point.

- base.iid:

  Influence functions for the baseline.

- beta.iid:

  Influence functions for the regression coefficients.

- cumhaz.time:

  Cumulative hazard at the specified time.

- strata:

  Strata indices.

## Details

The decomposition is based on the formula: \$\$ \sum_i \int_0^t
\frac{f(s)}{S_0(s)} dM\_{ki}(s) - P(t) \beta_k \$\$ where \\k\\ denotes
the stratum and \\i\\ the cluster.

## See also

[`phreg`](http://kkholst.github.io/mets/reference/phreg.md),
[`recreg`](http://kkholst.github.io/mets/reference/recreg.md),
[`cifreg`](http://kkholst.github.io/mets/reference/cifreg.md)

## Author

Thomas Scheike
