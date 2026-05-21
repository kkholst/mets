# Generate Random Fold Indices for Cross-Validation

Splits `n` observations into random folds for cross-validation.

## Usage

``` r
folds(n, folds = 10)
```

## Arguments

- n:

  number of observations.

- folds:

  number of folds (default 10).

## Value

A list of integer vectors, each containing indices for one fold.
