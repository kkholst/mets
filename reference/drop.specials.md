# Remove Special Terms from a Formula

Removes terms such as
[`strata()`](https://rdrr.io/pkg/survival/man/strata.html) or
[`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html) from a
formula/terms object.

## Usage

``` r
drop.specials(x, components, ...)
```

## Arguments

- x:

  a formula object.

- components:

  character vector of special term names to remove.

- ...:

  additional arguments passed to
  [`lava::Specials`](https://kkholst.github.io/lava/reference/internal.html).

## Value

The updated formula with an attribute `"variables"` containing the
extracted variable names from removed specials.
