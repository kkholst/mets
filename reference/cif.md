# Cumulative Incidence with Robust Standard Errors

Computes cumulative incidence functions with robust standard errors.

## Usage

``` r
cif(formula, data = data, cause = 1, cens.code = 0, death.code = NULL, ...)
```

## Arguments

- formula:

  Formula with 'Event' outcome and `strata` (only!).

- data:

  Data frame.

- cause:

  Cause of interest (default is `NULL`, which looks at all causes).

- cens.code:

  Censoring code (default is `"0"`).

- death.code:

  Alternative to `cens.code`; specifies codes of death.

- ...:

  Additional arguments passed to lower-level functions.

## Value

An object of class `"cif"` (extends `"phreg"`) containing:

- cumhaz:

  Cumulative incidence estimates.

- se.cumhaz:

  Standard errors.

- cause:

  Cause of interest.

## Author

Thomas Scheike

## Examples

``` r
data(bmt)
bmt$cluster <- sample(1:100, 408, replace = TRUE)
out1 <- cif(Event(time, cause) ~ +1, data = bmt, cause = 1)
out2 <- cif(Event(time, cause) ~ +1 + cluster(cluster), data = bmt, cause = 1)

par(mfrow = c(1, 2))
plot(out1, se = TRUE)
plot(out2, se = TRUE)
```
