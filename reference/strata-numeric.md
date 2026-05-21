# Stratified Cumulative and Summary Operations

Low-level helper functions for computing sums, cumulative sums, and
reverse cumulative sums within strata groups, as well as matrix
double-indexing.

## Usage

``` r
sumstrata(x, strata, nstrata)

cumsumstrata(x, strata, nstrata)

revcumsumstrata(x, strata, nstrata)

revcumsum(x)

matdoubleindex(x, rows, cols, xvec = NULL)

mdi(x, ...)
```

## Arguments

- x:

  numeric vector (or matrix for `matdoubleindex`).

- strata:

  integer vector of strata indices (0-based).

- nstrata:

  number of distinct strata.

- rows:

  row indices for `matdoubleindex`.

- cols:

  column indices for `matdoubleindex`.

- xvec:

  optional values to assign at indexed positions.

- ...:

  additional arguments (for `mdi`).

## Value

Numeric vector of the same length as `x` (or modified matrix).

## Author

Klaus K. Holst, Thomas Scheike
