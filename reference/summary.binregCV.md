# Summary method for binregCV objects

Summary method for binregCV objects

## Usage

``` r
# S3 method for class 'binregCV'
summary(
  object,
  transform = NULL,
  conf.level = 0.95,
  type = c("cv", "both"),
  ...
)
```

## Arguments

- object:

  A `binregCV` object, as returned by
  [`brier_binreg`](http://kkholst.github.io/mets/reference/brier_binreg.md)
  for a single evaluation time.

- transform:

  Optional back-transformation function applied to Brier scores and
  confidence limits, e.g. `exp` to undo the internal log scale.

- conf.level:

  Confidence level for the reported intervals (default 0.95).

- type:

  Either `"cv"` (default) for cross-validated Brier scores only, or
  `"both"` to also include apparent (in-sample) scores.

- ...:

  Not currently used.
