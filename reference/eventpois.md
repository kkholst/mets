# Extract survival estimates from lifetable analysis

Summary for survival analyses via the 'lifetable' function

## Usage

``` r
eventpois(
  object,
  ...,
  timevar,
  time,
  int.len,
  confint = FALSE,
  level = 0.95,
  individual = FALSE,
  length.out = 25
)
```

## Arguments

- object:

  glm object (poisson regression)

- ...:

  Contrast arguments

- timevar:

  Name of time variable

- time:

  Time points (optional)

- int.len:

  Time interval length (optional)

- confint:

  If TRUE confidence limits are supplied

- level:

  Level of confidence limits

- individual:

  Individual predictions

- length.out:

  Length of time vector

## Details

Summary for survival analyses via the 'lifetable' function

## Author

Klaus K. Holst
